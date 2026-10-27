//! `Re`, `Im`, `Conjugate` and `Arg` (CAS-DOC-15 §4, CAS-240), Xcas's `re`,
//! `im`, `conj` and `arg`.
//!
//! An expression's parts `a + i b`, `a` and `b` real, are found
//! structurally ([`parts`]): a number by its parts; what the domain or the
//! assumptions make real (in the Real domain a symbol, `Pi`, `E`), `Abs`
//! and the parts themselves, as `a + 0 i`; sums by parts; products and
//! integer powers by complex arithmetic, a negative power through the
//! conjugate over `a² + b²` (`1/(x + i y)` is `(x − i y)/(x² + y²)`); `Exp[u]`
//! of `u = a + i b` as `e^a (cos b + i sin b)`; `Sin`, `Cos`, `Tan`,
//! `ArcTan`, `Sinh`, `Cosh`, `Tanh` and `ArcSinh` of a real argument, and
//! a power of a positive base by a real exponent, real. Where they are
//! found, `Re` is `a`, `Im` is `b`, `Conjugate` is `a − i b`, and `Arg` is
//! `0` or `Pi` for a real of known sign, `±Pi/2` for an imaginary one, and
//! `ArcTan[b/a]` for `a > 0`, plus or minus `Pi` for `a < 0` by the sign of
//! `b`; otherwise each stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of complex numbers: z = a + i b, to the power k"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use crate::assume::{Set, Sign, Truth, element, sign_of};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest integer power of a complex base expanded.
const MAX_POWER: i64 = 16;

/// `Re[z]`, `Im[z]` or `Conjugate[z]`.
pub(super) fn part(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(f), [z]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let Some((a, b)) = parts(ev.context(), z) else {
        return Ok(None);
    };
    let result = match f {
        Symbol::RE => a,
        Symbol::IM => b,
        _ => a - i() * b,
    };
    // Only when it is not the call again, as `Re[f[x]]` would be.
    Ok((result != *e).then_some(result))
}

/// `Arg[z]`.
pub(super) fn arg(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [z] = e.args() else {
        return Ok(None);
    };
    let Some((a, b)) = parts(ev.context(), z) else {
        return Ok(None);
    };
    let (a, b) = (ev.eval(&a)?, ev.eval(&b)?);
    let sign = |x: &Expr| sign(ev.context(), x);
    let pi = || Expr::from(Symbol::PI);
    let half_pi = || Expr::from(Rational::new(Integer::ONE, Integer::from(2))) * pi();
    Ok(match (sign(&a), sign(&b)) {
        (Sign::Positive | Sign::Zero | Sign::NonNegative, Sign::Zero) => Some(Expr::from(0)),
        (Sign::Negative, Sign::Zero) => Some(pi()),
        (Sign::Zero, Sign::Positive) => Some(half_pi()),
        (Sign::Zero, Sign::Negative) => Some(-half_pi()),
        (Sign::Positive, _) => Some(Expr::call(Symbol::ARC_TAN, [b / a])),
        (Sign::Negative, Sign::Positive) => Some(Expr::call(Symbol::ARC_TAN, [b / a]) + pi()),
        (Sign::Negative, Sign::Negative) => Some(Expr::call(Symbol::ARC_TAN, [b / a]) - pi()),
        _ => None,
    })
}

fn i() -> Expr {
    Expr::from(Symbol::I)
}

/// The sign of a real `x`: by the assumptions, or a ball's.
fn sign(ctx: &Context, x: &Expr) -> Sign {
    match x.view() {
        View::Ball(b) => match b.sign() {
            Some(Ordering::Greater) => Sign::Positive,
            Some(Ordering::Less) => Sign::Negative,
            Some(Ordering::Equal) => Sign::Zero,
            None => Sign::Unknown,
        },
        _ => sign_of(ctx, x),
    }
}

/// `(a, b)` with `z = a + i b`, both real, built; `None` when not found.
pub(super) fn parts(ctx: &Context, z: &Expr) -> Option<(Expr, Expr)> {
    let zero = || Expr::from(0);
    match z.view() {
        View::Complex(c) => {
            return Some((Expr::from(c.re().clone()), Expr::from(c.im().clone())));
        }
        View::ComplexBall(c) => {
            return Some((Expr::from(c.re().clone()), Expr::from(c.im().clone())));
        }
        View::Int(_) | View::Integer(_) | View::Rational(_) | View::Ball(_) => {
            return Some((z.clone(), zero()));
        }
        _ => {}
    }
    if element(ctx, z, Set::Reals) == Truth::True {
        return Some((z.clone(), zero()));
    }
    let head = z.head().and_then(Expr::as_symbol)?;
    let args = z.args();
    match (head, args) {
        (Symbol::ABS | Symbol::RE | Symbol::IM | Symbol::ARG, [_]) => Some((z.clone(), zero())),
        (Symbol::PLUS, _) => {
            let (mut re, mut im) = (Vec::new(), Vec::new());
            for t in args {
                let (a, b) = parts(ctx, t)?;
                re.push(a);
                im.push(b);
            }
            Some((Expr::call(Symbol::PLUS, re), Expr::call(Symbol::PLUS, im)))
        }
        (Symbol::TIMES, _) => {
            let mut acc = (Expr::from(1), zero());
            for f in args {
                acc = times(&acc, &parts(ctx, f)?);
            }
            Some(acc)
        }
        (Symbol::POWER, [b, n]) => power(ctx, z, b, n),
        (Symbol::EXP, [u]) => {
            let (a, b) = parts(ctx, u)?;
            let scale = Expr::call(Symbol::EXP, [a]);
            Some((
                scale.clone() * Expr::call(Symbol::COS, [b.clone()]),
                scale * Expr::call(Symbol::SIN, [b]),
            ))
        }
        (
            Symbol::SIN
            | Symbol::COS
            | Symbol::TAN
            | Symbol::ARC_TAN
            | Symbol::SINH
            | Symbol::COSH
            | Symbol::TANH
            | Symbol::ARC_SINH,
            [u],
        ) if element(ctx, u, Set::Reals) == Truth::True => Some((z.clone(), zero())),
        _ => None,
    }
}

/// `(a + i b)(c + i d)`, built.
fn times((a, b): &(Expr, Expr), (c, d): &(Expr, Expr)) -> (Expr, Expr) {
    let real = |x: &Expr| *x == Expr::from(0);
    if real(b) && real(d) {
        return (a.clone() * c.clone(), Expr::from(0));
    }
    (
        a.clone() * c.clone() - b.clone() * d.clone(),
        a.clone() * d.clone() + b.clone() * c.clone(),
    )
}

/// The parts of `z = base^n`.
fn power(ctx: &Context, z: &Expr, base: &Expr, n: &Expr) -> Option<(Expr, Expr)> {
    let zero = || Expr::from(0);
    // A positive base to a real power is real.
    if sign_of(ctx, base) == Sign::Positive && element(ctx, n, Set::Reals) == Truth::True {
        return Some((z.clone(), zero()));
    }
    let k = n.as_i64().filter(|k| k.abs() <= MAX_POWER)?;
    let (a, b) = parts(ctx, base)?;
    if b == zero() {
        return Some((z.clone(), zero()));
    }
    let mut acc = (Expr::from(1), zero());
    for _ in 0..k.unsigned_abs() {
        acc = times(&acc, &(a.clone(), b.clone()));
    }
    if k >= 0 {
        return Some(acc);
    }
    // 1/(a + i b)^m = (a − i b)^m/(a² + b²)^m: acc is (a + i b)^m.
    let norm = (a.pow(2) + b.pow(2)).pow(k);
    Some((acc.0 * norm.clone(), -acc.1 * norm))
}

#[cfg(test)]
mod tests;
