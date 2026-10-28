//! `Floor`, `Ceiling`, `Round`, `Sign`, `Max` and `Min` (CAS-DOC-15 §5,
//! CAS-241), Xcas's `floor`, `ceil`, `round`, `sign`, `max` and `min`.
//!
//! - `Floor`, `Ceiling` and `Round` of a rational exactly, `Round[x]` as
//!   `Floor[x + 1/2]` (giac's `round(-5/2)` is −2); of a real ball, or a
//!   numeric expression approximated at 64 to [`MAX_BITS`] bits, the
//!   integer where both ends of the ball give it, otherwise they stay (giac
//!   evaluates approximately, with a warning); `Floor[n + x] = n +
//!   Floor[x]` for an integer `n`, and `x` itself where the assumptions make
//!   it an integer; the infinities as they are.
//! - `Sign`: −1, 0 or 1 by the assumptions, a ball, or the approximation;
//!   `z/Abs[z]` of a non-real number; a rational coefficient's sign out,
//!   `Sign[-2 x]` is `-Sign[x]`.
//! - `Max` and `Min`: lists flattened, the arguments compared by the sign
//!   of their difference (the assumptions, a ball, the approximation), the
//!   dominated ones dropped; one left is the value, none gives `-Infinity`
//!   for `Max` and `Infinity` for `Min`, and others stay, as `max(x,1)`.

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::approximation;
use crate::assume::{Set, Sign, Truth, element, sign_of};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The bits a numeric expression is approximated to at most.
const MAX_BITS: u64 = 4096;

/// Which rounding.
#[derive(Clone, Copy)]
enum Rounding {
    Floor,
    Ceiling,
    Round,
}

impl Rounding {
    fn of(self, q: &Rational) -> Integer {
        match self {
            Rounding::Floor => q.floor(),
            Rounding::Ceiling => q.ceil(),
            Rounding::Round => (q + &Rational::new(Integer::ONE, Integer::from(2))).floor(),
        }
    }
}

/// The ends of a real ball, exactly; `None` for an infinite one.
fn ends(x: &Expr) -> Option<(Rational, Rational)> {
    let View::Ball(b) = x.view() else {
        return None;
    };
    let mid = b.mid().to_rational();
    let rad = b.rad().to_float()?.to_rational();
    Some((&mid - &rad, &mid + &rad))
}

/// `x`'s value as a real ball or a rational, at increasing precision until
/// `decided` tells from it; `None` past [`MAX_BITS`].
fn decide<T>(
    ctx: &Context,
    x: &Expr,
    mut decided: impl FnMut(&Rational, &Rational) -> Option<T>,
) -> Option<T> {
    if let Some(q) = x.to_rational() {
        return decided(&q, &q);
    }
    if let Some((lo, hi)) = ends(x) {
        return decided(&lo, &hi);
    }
    if !x.is_numeric() {
        return None;
    }
    let mut bits = 64;
    while bits <= MAX_BITS {
        let value = approximation(ctx, x, bits).ok()?;
        if let Some(q) = value.to_rational() {
            return decided(&q, &q);
        }
        let (lo, hi) = ends(&value)?;
        if let Some(t) = decided(&lo, &hi) {
            return Some(t);
        }
        bits *= 4;
    }
    None
}

/// `Floor[x]`, `Ceiling[x]` and `Round[x]`.
pub(super) fn round(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(f), [x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let rounding = match f {
        Symbol::FLOOR => Rounding::Floor,
        Symbol::CEILING => Rounding::Ceiling,
        _ => Rounding::Round,
    };
    let ctx = ev.context();
    if super::number::Infinity::of(x).is_some() {
        return Ok(Some(x.clone()));
    }
    if element(ctx, x, Set::Integers) == Truth::True {
        return Ok(Some(x.clone()));
    }
    if let Some(n) = decide(ctx, x, |lo, hi| {
        let (low, high) = (rounding.of(lo), rounding.of(hi));
        (low == high).then_some(low)
    }) {
        return Ok(Some(Expr::from(n)));
    }
    // An integer term out of a sum.
    if let (Some(Symbol::PLUS), args) = (x.head().and_then(Expr::as_symbol), x.args()) {
        let (integers, rest): (Vec<Expr>, Vec<Expr>) = args
            .iter()
            .cloned()
            .partition(|t| element(ctx, t, Set::Integers) == Truth::True);
        if !integers.is_empty() && !rest.is_empty() {
            let inner = Expr::call(f, [Expr::call(Symbol::PLUS, rest)]);
            return Ok(Some(Expr::call(Symbol::PLUS, integers) + inner));
        }
    }
    Ok(None)
}

/// The sign of a real `x`, from the assumptions, a ball or the
/// approximation; also the step's and the delta's (CAS-273).
pub(super) fn real_sign(ctx: &Context, x: &Expr) -> Option<Ordering> {
    match sign_of(ctx, x) {
        Sign::Positive => return Some(Ordering::Greater),
        Sign::Negative => return Some(Ordering::Less),
        Sign::Zero => return Some(Ordering::Equal),
        _ => {}
    }
    decide(ctx, x, |lo, hi| {
        let (a, b) = (lo.cmp(&Rational::ZERO), hi.cmp(&Rational::ZERO));
        (a == b).then_some(a)
    })
}

/// `Sign[x]`.
pub(super) fn sign(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let ctx = ev.context();
    if matches!(x.view(), View::Complex(_)) {
        return Ok(Some(x.clone() / Expr::call(Symbol::ABS, [x.clone()])));
    }
    if let Some(s) = real_sign(ctx, x) {
        return Ok(Some(Expr::from(match s {
            Ordering::Greater => 1,
            Ordering::Less => -1,
            Ordering::Equal => 0,
        })));
    }
    // A rational coefficient's sign out.
    if let (Some(Symbol::TIMES), [c, rest @ ..]) = (x.head().and_then(Expr::as_symbol), x.args())
        && let Some(c) = c.to_rational()
        && !rest.is_empty()
    {
        let inner = Expr::call(Symbol::SIGN, [Expr::call(Symbol::TIMES, rest.to_vec())]);
        return Ok(Some(if c.is_negative() { -inner } else { inner }));
    }
    Ok(None)
}

/// `Max[…]` and `Min[…]`.
pub(super) fn extremum(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(f) = e.head().and_then(Expr::as_symbol) else {
        return Ok(None);
    };
    let wanted = if f == Symbol::MAX {
        Ordering::Greater
    } else {
        Ordering::Less
    };
    // The arguments, lists flattened.
    let mut args = Vec::new();
    let mut stack: Vec<&Expr> = e.args().iter().rev().collect();
    while let Some(a) = stack.pop() {
        if a.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
            stack.extend(a.args().iter().rev());
        } else if !args.contains(a) {
            args.push(a.clone());
        }
    }
    if args.is_empty() {
        let infinity = Expr::from(Symbol::INFINITY);
        return Ok(Some(if f == Symbol::MAX {
            -infinity
        } else {
            infinity
        }));
    }
    // Drop each argument another one dominates.
    let mut kept: Vec<Expr> = Vec::new();
    for a in args {
        let mut dominated = false;
        let mut i = 0;
        while i < kept.len() {
            let difference = ev.eval(&(a.clone() - kept[i].clone()))?;
            match real_sign(ev.context(), &difference) {
                Some(s) if s == wanted => {
                    kept.remove(i);
                }
                Some(_) => {
                    dominated = true;
                    break;
                }
                None => i += 1,
            }
        }
        if !dominated {
            kept.push(a);
        }
    }
    if let [one] = kept.as_slice() {
        return Ok(Some(one.clone()));
    }
    let result = Expr::call(f, kept);
    Ok((result != *e).then_some(result))
}

#[cfg(test)]
mod tests;
