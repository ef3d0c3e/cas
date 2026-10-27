//! Radicals in the search (CAS-DOC-14 §6, CAS-232): square roots of
//! positive integers, `Power[n, 1/2]` in canonical form (CAS-DOC-4 §6.1).
//!
//! - Rationalized denominators: the denominator of `Together[e]`, while it
//!   has a square root `√n`, times its conjugate in `√n` (the root negated)
//!   and expanded, which squares the root away, at most [`MAX_ROOTS`]
//!   times; the numerator times the conjugates, expanded.
//!   `1/(1+sqrt(2))` is `sqrt(2)-1`.
//! - Denested square roots: `√(α + β √n)` (or its inverse), rationals `α`,
//!   `β` and `n` with `α + β √n > 0`, as `√((α + s)/2) ± √((α − s)/2)` when
//!   `α² − β² n = s²` (CAS-211's): `sqrt(3+2*sqrt(2))` is `1+sqrt(2)`.

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::super::solve::{surd_sign, surd_sqrt};
use super::evaluated;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The square roots a denominator is rid of at most.
const MAX_ROOTS: usize = 4;

/// Whether `e` is `±1/2`: `Some(true)` for `1/2`.
fn half(e: &Expr) -> Option<bool> {
    let q = e.to_rational()?;
    let h = Rational::new(Integer::ONE, Integer::from(2));
    if q == h {
        Some(true)
    } else if q == -&h {
        Some(false)
    } else {
        None
    }
}

/// A square root `Power[n, 1/2]` of a positive integer in `e`.
fn square_root(e: &Expr) -> Option<Expr> {
    let mut stack = vec![e];
    while let Some(x) = stack.pop() {
        if let View::Apply { head, args } = x.view() {
            if head.as_symbol() == Some(Symbol::POWER)
                && let [b, h] = args
                && half(h) == Some(true)
                && b.to_integer().is_some_and(|n| n.is_positive())
            {
                return Some(x.clone());
            }
            stack.extend(args);
        }
    }
    None
}

/// `e` with its denominator rid of square roots; `None` when it has none,
/// or keeps one.
pub(super) fn rationalized(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    if square_root(e).is_none() {
        return Ok(None);
    }
    let Some(together) = evaluated(ev, &Expr::call(Symbol::TOGETHER, [e.clone()]))? else {
        return Ok(None);
    };
    let part =
        |ev: &mut Eval<'_>, head: Symbol| evaluated(ev, &Expr::call(head, [together.clone()]));
    let (Some(mut num), Some(mut den)) =
        (part(ev, Symbol::NUMERATOR)?, part(ev, Symbol::DENOMINATOR)?)
    else {
        return Ok(None);
    };
    if square_root(&den).is_none() {
        return Ok(None);
    }
    for _ in 0..MAX_ROOTS {
        let Some(r) = square_root(&den) else {
            break;
        };
        let conjugate = den.substitute(&[(r.clone(), -r)]);
        num *= conjugate.clone();
        let Some(d) = evaluated(ev, &Expr::call(Symbol::EXPAND, [den * conjugate]))? else {
            return Ok(None);
        };
        den = d;
    }
    if square_root(&den).is_some() || den == Expr::from(0) {
        return Ok(None);
    }
    // Over a number, expanded through: evaluation leaves −(1 − √2).
    Ok(Some(if den.is_number() {
        Expr::call(Symbol::EXPAND, [num / den])
    } else {
        Expr::call(Symbol::EXPAND, [num]) / den
    }))
}

/// `α + β √n` as `(α, β, n)`, with `β ≠ 0`: rationals and rational
/// multiples of the square root of one positive integer.
fn surd(e: &Expr) -> Option<(Rational, Rational, Rational)> {
    let terms: &[Expr] = match e.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => e.args(),
        _ => std::slice::from_ref(e),
    };
    let mut alpha = Rational::ZERO;
    let mut beta = Rational::ZERO;
    let mut n: Option<Rational> = None;
    for t in terms {
        if let Some(q) = t.to_rational() {
            alpha = &alpha + &q;
            continue;
        }
        let (c, root) = match (t.head().and_then(Expr::as_symbol), t.args()) {
            (Some(Symbol::TIMES), [c, root]) => (c.to_rational()?, root),
            _ => (Rational::ONE, t),
        };
        let [base, exponent] = root.args() else {
            return None;
        };
        if root.head().and_then(Expr::as_symbol) != Some(Symbol::POWER)
            || half(exponent) != Some(true)
        {
            return None;
        }
        let base = base.to_rational().filter(Rational::is_positive)?;
        if n.as_ref().is_some_and(|n| *n != base) {
            return None;
        }
        n = Some(base);
        beta = &beta + &c;
    }
    let n = n?;
    (!beta.is_zero()).then_some((alpha, beta, n))
}

/// `e` with its square roots of positive surds denested where they can be;
/// `None` when nothing changes.
pub(super) fn denested(e: &Expr) -> Option<Expr> {
    let r = e.map_bottom_up(|x| {
        let (Some(Symbol::POWER), [b, h]) = (x.head().and_then(Expr::as_symbol), x.args()) else {
            return x.clone();
        };
        let (Some(positive), Some((alpha, beta, n))) = (half(h), surd(b)) else {
            return x.clone();
        };
        if surd_sign(&alpha, &beta, &n) != Ordering::Greater {
            return x.clone();
        }
        let root = surd_sqrt(&alpha, &beta, &n);
        if positive { root } else { root.pow(-1) }
    });
    (r != *e).then_some(r)
}
