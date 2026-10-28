//! Trigonometric substitutions (CAS-DOC-36 §3, CAS-397): the rules of
//! `HalfTan`, `HalfTanHyperbolicToExp`, `TanToSinCos`, `SinToCosTan`,
//! `CosToSinTan`, `TanToSinCos2` and `TanToCosSin2`, giac's `halftan`,
//! `halftan_hyp2exp`, `tan2sincos`, `sin2costan`, `cos2sintan`,
//! `tan2sincos2` and `tan2cossin2`.
//!
//! - `HalfTan`: `sin u = 2t/(1 + t²)`, `cos u = (1 − t²)/(1 + t²)` and
//!   `tan u = 2t/(1 − t²)`, `t = tan(u/2)`.
//! - `HalfTanHyperbolicToExp`: `HalfTan`'s rule and `HyperbolicToExp`'s.
//! - `TanToSinCos`: `tan u = sin u/cos u`; `SinToCosTan`: `sin u = cos u
//!   tan u`; `CosToSinTan`: `cos u = sin u/tan u`.
//! - `TanToSinCos2`: `tan u = sin 2u/(1 + cos 2u)`; `TanToCosSin2`: `tan u
//!   = (1 − cos 2u)/sin 2u`.
//!
//! `Cot`, `Sec` and `Csc` are quotients of `Sin` and `Cos` once evaluated
//! (CAS-DOC-15 §3), so the rules meet them too.

use super::exponential::hyperbolic;
use super::head;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

fn sin(u: Expr) -> Expr {
    Expr::call(Symbol::SIN, [u])
}

fn cos(u: Expr) -> Expr {
    Expr::call(Symbol::COS, [u])
}

fn tan(u: Expr) -> Expr {
    Expr::call(Symbol::TAN, [u])
}

/// `part` as `(f, u)` for `Sin[u]`, `Cos[u]` or `Tan[u]`.
fn trigonometric(part: &Expr) -> Option<(Symbol, &Expr)> {
    match (head(part), part.args()) {
        (Some(f @ (Symbol::SIN | Symbol::COS | Symbol::TAN)), [u]) => Some((f, u)),
        _ => None,
    }
}

/// `HalfTan`'s rule.
pub(super) fn half_tan(_: &mut Eval<'_>, part: &Expr) -> Found {
    let Some((f, u)) = trigonometric(part) else {
        return Ok(None);
    };
    let t = tan(u.clone() / Expr::from(2));
    let one = || Expr::from(1);
    let square = t.pow(2);
    Ok(Some(match f {
        Symbol::SIN => Expr::from(2) * t / (one() + square),
        Symbol::COS => (one() - square.clone()) / (one() + square),
        _ => Expr::from(2) * t / (one() - square),
    }))
}

/// `HalfTanHyperbolicToExp`'s rule.
pub(super) fn half_tan_hyperbolic(ev: &mut Eval<'_>, part: &Expr) -> Found {
    match half_tan(ev, part)? {
        Some(form) => Ok(Some(form)),
        None => hyperbolic(ev, part),
    }
}

/// `TanToSinCos`'s rule.
pub(super) fn tan_to_sin_cos(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match trigonometric(part) {
        Some((Symbol::TAN, u)) => Some(sin(u.clone()) / cos(u.clone())),
        _ => None,
    })
}

/// `SinToCosTan`'s rule.
pub(super) fn sin_to_cos_tan(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match trigonometric(part) {
        Some((Symbol::SIN, u)) => Some(cos(u.clone()) * tan(u.clone())),
        _ => None,
    })
}

/// `CosToSinTan`'s rule.
pub(super) fn cos_to_sin_tan(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match trigonometric(part) {
        Some((Symbol::COS, u)) => Some(sin(u.clone()) / tan(u.clone())),
        _ => None,
    })
}

/// `TanToSinCos2`'s rule.
pub(super) fn tan_to_sin_cos_2(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match trigonometric(part) {
        Some((Symbol::TAN, u)) => {
            let double = Expr::from(2) * u.clone();
            Some(sin(double.clone()) / (Expr::from(1) + cos(double)))
        }
        _ => None,
    })
}

/// `TanToCosSin2`'s rule.
pub(super) fn tan_to_cos_sin_2(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match trigonometric(part) {
        Some((Symbol::TAN, u)) => {
            let double = Expr::from(2) * u.clone();
            Some((Expr::from(1) - cos(double.clone())) / sin(double))
        }
        _ => None,
    })
}

#[cfg(test)]
mod tests;
