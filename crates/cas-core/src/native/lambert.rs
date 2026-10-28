//! `LambertW` (CAS-DOC-31 §3.3, CAS-356), Xcas's `LambertW(x)` and
//! `LambertW(x, k)`: the branches `k = 0` and `k = −1` of the inverse of
//! `w e^w`.
//!
//! - `LambertW[x, 0]` is `LambertW[x]`; another branch leaves the call.
//! - `W(0) = 0` on branch 0, and `W(z e^z) = z` for a rational `z ≥ −1` on
//!   branch 0 and `z ≤ −1` on branch −1: `W(e) = 1`, and `W(−1/e) = −1` on
//!   both.
//! - A real ball inside the branch's domain by `Ball::lambert_w`. Outside
//!   it, where giac's values are complex (`LambertW(-0.5)`), the call
//!   stays.
//! - `E^(k W(z)) = (z/W(z))^k` for an integer `k`, as `w e^w = z` on every
//!   branch ([`exponential`]).
//! - Anything else stays.

use cas_num::Rational;

use super::number::{Num, Precision};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `LambertW[x]` and `LambertW[x, k]`.
pub(super) fn lambert_w(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, lower) = match e.args() {
        [x] => (x, false),
        [x, k] => match k.as_i64() {
            Some(0) => return Ok(Some(Expr::call(Symbol::LAMBERT_W, [x.clone()]))),
            Some(-1) => (x, true),
            _ => return Ok(None),
        },
        _ => return Ok(None),
    };
    if *x == Expr::from(0) && !lower {
        return Ok(Some(Expr::from(0)));
    }
    if let Some(z) = inverse(x) {
        let minus_one = Rational::from(cas_num::Integer::from(-1));
        let on_branch = if lower {
            z <= minus_one
        } else {
            z >= minus_one
        };
        if on_branch {
            return Ok(Some(Expr::from(z)));
        }
    }
    if let Some(ref number @ Num::Ball(ref ball)) = Num::of(x) {
        let bits = number.precision(Precision::of(ev.context()));
        let value = ball.lambert_w(lower, bits + 16);
        if value.is_finite() {
            return Ok(Some(
                Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr(),
            ));
        }
    }
    Ok(None)
}

/// `z` for `x = z e^z` with a rational `z`: `E` itself, and `z E^z`.
fn inverse(x: &Expr) -> Option<Rational> {
    if x.as_symbol() == Some(Symbol::E) {
        return Some(Rational::ONE);
    }
    let [factor, power] = x.args() else {
        return None;
    };
    if x.head() != Some(&Expr::from(Symbol::TIMES)) {
        return None;
    }
    let z = factor.to_rational()?;
    let (Some(Symbol::POWER), [base, exponent]) =
        (power.head().and_then(Expr::as_symbol), power.args())
    else {
        return None;
    };
    (base.as_symbol() == Some(Symbol::E) && exponent.to_rational() == Some(z.clone())).then_some(z)
}

/// `E^n` for `n = k LambertW[z]` (or with a branch), `k` an integer: `(z/W)^k`.
pub(super) fn exponential(n: &Expr) -> Option<Expr> {
    let is_w = |w: &Expr| w.head() == Some(&Expr::from(Symbol::LAMBERT_W)) && !w.args().is_empty();
    let (k, w) = if is_w(n) {
        (1, n)
    } else {
        let [c, w] = n.args() else {
            return None;
        };
        if n.head() != Some(&Expr::from(Symbol::TIMES)) || !is_w(w) {
            return None;
        }
        (c.as_i64()?, w)
    };
    let z = w.args()[0].clone();
    Some((z / w.clone()).pow(k))
}

#[cfg(test)]
mod tests;
