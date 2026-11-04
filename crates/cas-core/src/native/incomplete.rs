//! The incomplete gamma and beta functions (CAS-DOC-46 §3, CAS-469): Xcas's
//! `igamma`, `ugamma` (`Gamma` with two arguments) and `Beta` with three.
//!
//! - `LowerGamma[a, x]` is `γ(a, x) = ∫₀ˣ t^(a−1) e^(−t) dt`, `Gamma[a, x]`
//!   is `Γ(a, x) = ∫ₓ^∞ t^(a−1) e^(−t) dt` and `Beta[a, b, x]` is `B_x(a, b)
//!   = ∫₀ˣ t^(a−1) (1 − t)^(b−1) dt`; a last argument 1 makes each
//!   regularized, `P(a, x)`, `Q(a, x)` and `I_x(a, b)`, as giac's
//!   `igamma(a, x, 1)`.
//! - **Balls**: real arguments with a ball among them, `a, b > 0`, `x ≥ 0`
//!   (`x ≤ 1` for `Beta`), by cas-num's `gamma_p`, `gamma_q` and `beta_i`,
//!   times `Γ(a)` or `B(a, b)` unless regularized.
//! - **Ends**: `γ(a, 0) = 0`, `Γ(a, 0) = Γ(a)`, `B_0 = 0`, `B_1(a, b) =
//!   B(a, b)`, the regularized ones 0 or 1.
//! - **Closed forms** at natural and half-integer `a` (and `b`) up to
//!   [`MAX_CLOSED`], those of the χ², Student and Fisher distributions
//!   (`statistics::distributions::incomplete`): giac's `igamma(2,x)` is `1 −
//!   e^(−x)(x + 1)`.
//!
//! Anything else stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the functions: γ(a, x), B_x(a, b), P, Q and I, w bits"
)]

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::number::{Num, Precision};
use super::statistics::distributions::incomplete::{
    beta as beta_closed, lower_gamma as lower_gamma_closed,
};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest `a` (and `b`) given in closed form.
const MAX_CLOSED: u64 = 20;

/// Whether the arguments end in the flag 1 of the regularized form, and
/// the others; `None` for a last argument that is another value.
fn split(args: &[Expr], n: usize) -> Option<(&[Expr], bool)> {
    match args.len() {
        k if k == n => Some((args, false)),
        k if k == n + 1 && args[n] == Expr::from(1) => Some((&args[..n], true)),
        _ => None,
    }
}

/// `2a` for a positive natural or half-integer `a` up to [`MAX_CLOSED`].
fn twice(a: &Expr) -> Option<u64> {
    let q = a.to_rational()?;
    let t = &q * &Rational::from(Integer::from(2));
    let t = t.is_integer().then(|| t.numer().to_i64())??;
    u64::try_from(t)
        .ok()
        .filter(|t| (1..=2 * MAX_CLOSED).contains(t))
}

/// The real balls of the arguments at the bits ball arithmetic keeps, when
/// one is a ball and all are real.
fn balls(ev: &Eval<'_>, args: &[Expr]) -> Option<(Vec<Ball>, u64)> {
    let numbers: Vec<Num> = args.iter().map(Num::of).collect::<Option<_>>()?;
    if !numbers.iter().any(|n| matches!(n, Num::Ball(_))) {
        return None;
    }
    let precision = Precision::of(ev.context());
    let bits = numbers
        .iter()
        .filter(|n| matches!(n, Num::Ball(_)))
        .map(|n| n.precision(precision))
        .min()?;
    let balls = numbers
        .iter()
        .map(|n| match n {
            Num::Real(q) => Some(Ball::from_rational(q, bits + 64)),
            Num::Ball(b) => Some(Ball::clone(b)),
            Num::Complex(_) | Num::ComplexBall(_) => None,
        })
        .collect::<Option<Vec<_>>>()?;
    Some((balls, bits))
}

/// A ball value, or the call left when it is indeterminate.
fn ball_value(value: &Ball, bits: u64) -> Option<Expr> {
    value
        .is_finite()
        .then(|| Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr())
}

fn gamma(a: &Expr) -> Expr {
    Expr::call(Symbol::GAMMA, [a.clone()])
}

/// `LowerGamma[a, x]` and `LowerGamma[a, x, 1]`.
pub(super) fn lower_gamma(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((args, regularized)) = split(e.args(), 2) else {
        return Ok(None);
    };
    let [a, x] = args else {
        return Ok(None);
    };
    if let Some((b, bits)) = balls(ev, &[a.clone(), x.clone()]) {
        let w = bits + 16;
        let p = Ball::gamma_p(&b[0], &b[1], w);
        let value = if regularized {
            p
        } else {
            p.mul(&b[0].gamma(w), w)
        };
        return Ok(ball_value(&value, bits));
    }
    if *x == Expr::from(0) {
        return Ok(Some(Expr::from(0)));
    }
    let Some(p) = twice(a).and_then(|t| lower_gamma_closed(t, x)) else {
        return Ok(None);
    };
    let result = if regularized { p } else { p * gamma(a) };
    ev.step("igamma", e, &result);
    Ok(Some(result))
}

/// `Gamma[a, x]` and `Gamma[a, x, 1]`, from `Gamma`'s native.
pub(super) fn upper_gamma(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((args, regularized)) = split(e.args(), 2) else {
        return Ok(None);
    };
    let [a, x] = args else {
        return Ok(None);
    };
    if let Some((b, bits)) = balls(ev, &[a.clone(), x.clone()]) {
        let w = bits + 16;
        let q = Ball::gamma_q(&b[0], &b[1], w);
        let value = if regularized {
            q
        } else {
            q.mul(&b[0].gamma(w), w)
        };
        return Ok(ball_value(&value, bits));
    }
    if *x == Expr::from(0) {
        return Ok(Some(if regularized { Expr::from(1) } else { gamma(a) }));
    }
    let Some(p) = twice(a).and_then(|t| lower_gamma_closed(t, x)) else {
        return Ok(None);
    };
    // Expanded: `Times[-1, Plus[…]]` would stay as it is.
    let q = Expr::call(Symbol::EXPAND, [1 - p]);
    let result = if regularized { q } else { q * gamma(a) };
    ev.step("ugamma", e, &result);
    Ok(Some(result))
}

/// `Beta[a, b, x]` and `Beta[a, b, x, 1]`, from `Beta`'s native.
pub(super) fn incomplete_beta(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((args, regularized)) = split(e.args(), 3) else {
        return Ok(None);
    };
    let [a, b, x] = args else {
        return Ok(None);
    };
    let complete = Expr::call(Symbol::BETA, [a.clone(), b.clone()]);
    if let Some((v, bits)) = balls(ev, &[a.clone(), b.clone(), x.clone()]) {
        let w = bits + 16;
        let i = Ball::beta_i(&v[0], &v[1], &v[2], w);
        let value = if regularized {
            i
        } else {
            let whole = v[0]
                .gamma(w)
                .mul(&v[1].gamma(w), w)
                .div(&v[0].add(&v[1], w).gamma(w), w);
            i.mul(&whole, w)
        };
        return Ok(ball_value(&value, bits));
    }
    if *x == Expr::from(0) {
        return Ok(Some(Expr::from(0)));
    }
    if *x == Expr::from(1) {
        return Ok(Some(if regularized { Expr::from(1) } else { complete }));
    }
    let Some(i) = twice(a)
        .zip(twice(b))
        .and_then(|(ta, tb)| beta_closed(ta, tb, x))
    else {
        return Ok(None);
    };
    let result = if regularized { i } else { i * complete };
    ev.step("Beta", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
