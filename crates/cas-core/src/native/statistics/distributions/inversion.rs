//! Monotone inversion (CAS-DOC-23 §5, CAS-302), and `InverseErf` (§3).
//!
//! - [`invert`]: the `x` with `F(x) = q` for an increasing `F`, an
//!   expression in `x`, as a ball at the working precision. A bracket comes
//!   first, the ends where `F` is below and above `q`, doubled from
//!   `[−1, 1]`, or from `[low, low + 1]` on a support that starts at `low`.
//!   Then `FSolve`'s certified bisection (CAS-214) runs between them. It is
//!   for the inverses of continuous distributions without a closed form.
//! - `InverseErf[y]` (`erfinv`), the inverse of `Erf` on `(−1, 1)`, which
//!   giac lacks: 0 at 0, `±∞` at `±1`, odd (`erfinv(−y) = −erfinv(y)`). A
//!   real ball is the inversion of `Erf` at its precision, so that `N`
//!   reaches any target. An exact number stays, for `N`; one beyond `±1`
//!   stays with a warning.

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::super::super::number::{Num, Precision};
use super::super::super::solve::fsolve::{bisected, resolution, sign};
use super::super::super::special::negated;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The most doublings of a bracket.
const MAX_DOUBLINGS: usize = 64;

/// The `x` with `f = q` for `f` increasing in `x`, a ball at `bits` bits;
/// `low` the least value of `x` where `f` is defined, when there is one.
/// `None` when no bracket is found, or a value of `f` is not decided.
pub(in crate::native) fn invert(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    q: &Expr,
    low: Option<&Rational>,
    bits: u64,
) -> Option<Expr> {
    let g = f - q;
    let (mut lo, mut hi) = match low {
        Some(low) => (low.clone(), low + &Rational::ONE),
        None => (-Rational::ONE, Rational::ONE),
    };
    let two = Rational::from(Integer::from(2));
    for _ in 0..MAX_DOUBLINGS {
        let below = sign(ctx, &g, x, &lo, bits);
        let above = sign(ctx, &g, x, &hi, bits);
        match (below, above) {
            (Some(Ordering::Equal), _) => return Some(Expr::from(lo)),
            (_, Some(Ordering::Equal)) => return Some(Expr::from(hi)),
            (Some(Ordering::Less), Some(Ordering::Greater)) => {
                return bisected(ctx, &g, x, lo, hi, Ordering::Less, bits);
            }
            // Too far left: the bracket moves right, twice as wide.
            (Some(Ordering::Less), Some(Ordering::Less)) => {
                let width = &hi - &lo;
                lo = hi.clone();
                hi = &hi + &(&width * &two);
            }
            // Too far right, where there is room.
            (Some(Ordering::Greater), Some(Ordering::Greater)) if low.is_none() => {
                let width = &hi - &lo;
                hi = lo.clone();
                lo = &lo - &(&width * &two);
            }
            // An end within rounding of the root, an exact root there:
            // Cauchy's inverse at 3/4 is 1, the end of `[0, 1]` (CAS-306).
            (None, _) => return around(ctx, &g, x, &lo, bits),
            (_, None) => return around(ctx, &g, x, &hi, bits),
            _ => return None,
        }
    }
    None
}

/// The root of the increasing `g` within rounding of `point`, as the ball
/// of the resolution around it, when the signs there bracket it.
fn around(ctx: &Context, g: &Expr, x: &Expr, point: &Rational, bits: u64) -> Option<Expr> {
    let delta = resolution(point, bits);
    let (below, above) = (point - &delta, point + &delta);
    // At the bisection's own precision: the values there are of the
    // resolution's size.
    let bracketed = sign(ctx, g, x, &below, bits + 32) == Some(Ordering::Less)
        && sign(ctx, g, x, &above, bits + 32) == Some(Ordering::Greater);
    if bracketed {
        bisected(ctx, g, x, below, above, Ordering::Less, bits)
    } else {
        None
    }
}

/// `InverseErf[y]`.
pub(in crate::native) fn inverse_erf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [y] = e.args() else {
        return Ok(None);
    };
    if *y == Expr::from(0) {
        return Ok(Some(Expr::from(0)));
    }
    let infinity = Expr::from(Symbol::INFINITY);
    if *y == Expr::from(1) {
        return Ok(Some(infinity));
    }
    if *y == Expr::from(-1) {
        return Ok(Some(ev.eval(&-infinity)?));
    }
    if let Some(ref number @ Num::Ball(_)) = Num::of(y) {
        let bits = number.precision(Precision::of(ev.context()));
        let x = Expr::symbol("InverseErf$x");
        let erf = Expr::call(Symbol::ERF, [x.clone()]);
        return Ok(invert(ev.context(), &erf, &x, y, None, bits));
    }
    if let Some(z) = negated(y) {
        return Ok(Some(-Expr::call(Symbol::INVERSE_ERF, [z])));
    }
    if y.to_rational().is_some_and(|q| q.abs() > Rational::ONE) {
        ev.warn(Warning::Native {
            message: "erfinv is the inverse of erf on [-1, 1]".into(),
            expr: e.clone(),
        });
    }
    Ok(None)
}
