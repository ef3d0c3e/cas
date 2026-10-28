//! Inverse trigonometric functions (CAS-DOC-36 §6, CAS-400): the rules of
//! `ArcCosToArcSin`, `ArcSinToArcCos`, `ArcCosToArcTan`, `ArcSinToArcTan`,
//! `ArcTanToArcCos`, `ArcTanToArcSin` and `ArcTrigToLog`, giac's
//! `acos2asin`, `asin2acos`, `acos2atan`, `asin2atan`, `atan2acos`,
//! `atan2asin` and `atrig2ln`.
//!
//! - `acos u = π/2 − asin u` and `asin u = π/2 − acos u`, everywhere.
//! - `acos u = π/2 − atan(u/√(1 − u²))` and `asin u = atan(u/√(1 −
//!   u²))`, formal: for real `u` in `(−1, 1)`.
//! - `atan u = π/2 − acos(u/√(1 + u²))` and `atan u = asin(u/√(1 + u²))`,
//!   formal: for real `u`.
//! - `ArcTrigToLog`: `asin u = −i ln(i u + √(1 − u²))`, `acos u = π/2 + i
//!   ln(i u + √(1 − u²))` and `atan u = (i/2) ln((i + u)/(i − u))`, the
//!   principal branches' logarithmic forms (DLMF 4.23(iv)). They agree off
//!   the branch cuts; on a cut, where this CAS takes the side above it, a
//!   composition of principal functions can take the other: `asin(7/3)` is
//!   `π/2 + i acosh(7/3)`, its logarithmic form `π/2 − i acosh(7/3)`. giac's
//!   `atrig2ln(asin(x))` is `i*ln(x+sqrt(x^2-1))+pi/2`, the same on `[−1,
//!   1]`.
//!
//! `ArcCot`, `ArcSec` and `ArcCsc` are `ArcTan`, `ArcCos` and `ArcSin` of
//! `1/u` once evaluated (CAS-DOC-15 §3).

use cas_num::{Integer, Rational};

use super::head;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

fn half() -> Expr {
    Expr::from(Rational::new(Integer::ONE, Integer::from(2)))
}

fn half_pi() -> Expr {
    half() * Expr::from(Symbol::PI)
}

fn i() -> Expr {
    Expr::from(Symbol::I)
}

/// `√(1 + s u²)`, built.
fn root(u: &Expr, s: i64) -> Expr {
    (Expr::from(1) + Expr::from(s) * u.pow(2)).pow(half())
}

/// `part` as `(f, u)` for `ArcSin[u]`, `ArcCos[u]` or `ArcTan[u]`.
fn inverse(part: &Expr) -> Option<(Symbol, &Expr)> {
    match (head(part), part.args()) {
        (Some(f @ (Symbol::ARC_SIN | Symbol::ARC_COS | Symbol::ARC_TAN)), [u]) => Some((f, u)),
        _ => None,
    }
}

/// The rule that rewrites `ArcSin`, `ArcCos` or `ArcTan`, the `from`
/// function, by `to`.
fn conversion(part: &Expr, from: Symbol, to: fn(&Expr) -> Expr) -> Found {
    Ok(match inverse(part) {
        Some((f, u)) if f == from => Some(to(u)),
        _ => None,
    })
}

/// `ArcCosToArcSin`'s rule.
pub(super) fn arc_cos_to_arc_sin(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_COS, |u| {
        half_pi() - Expr::call(Symbol::ARC_SIN, [u.clone()])
    })
}

/// `ArcSinToArcCos`'s rule.
pub(super) fn arc_sin_to_arc_cos(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_SIN, |u| {
        half_pi() - Expr::call(Symbol::ARC_COS, [u.clone()])
    })
}

/// `ArcCosToArcTan`'s rule.
pub(super) fn arc_cos_to_arc_tan(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_COS, |u| {
        half_pi() - Expr::call(Symbol::ARC_TAN, [u.clone() / root(u, -1)])
    })
}

/// `ArcSinToArcTan`'s rule.
pub(super) fn arc_sin_to_arc_tan(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_SIN, |u| {
        Expr::call(Symbol::ARC_TAN, [u.clone() / root(u, -1)])
    })
}

/// `ArcTanToArcCos`'s rule.
pub(super) fn arc_tan_to_arc_cos(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_TAN, |u| {
        half_pi() - Expr::call(Symbol::ARC_COS, [u.clone() / root(u, 1)])
    })
}

/// `ArcTanToArcSin`'s rule.
pub(super) fn arc_tan_to_arc_sin(_: &mut Eval<'_>, part: &Expr) -> Found {
    conversion(part, Symbol::ARC_TAN, |u| {
        Expr::call(Symbol::ARC_SIN, [u.clone() / root(u, 1)])
    })
}

/// `ArcTrigToLog`'s rule.
pub(super) fn arc_trig_to_log(_: &mut Eval<'_>, part: &Expr) -> Found {
    let Some((f, u)) = inverse(part) else {
        return Ok(None);
    };
    let log = |v: Expr| Expr::call(Symbol::LOG, [v]);
    // ln(i u + √(1 − u²)) = i asin u.
    let arcsine = || log(i() * u.clone() + root(u, -1));
    Ok(Some(match f {
        Symbol::ARC_SIN => -i() * arcsine(),
        Symbol::ARC_COS => half_pi() + i() * arcsine(),
        _ => half() * i() * log((i() + u.clone()) / (i() - u.clone())),
    }))
}

#[cfg(test)]
mod tests;
