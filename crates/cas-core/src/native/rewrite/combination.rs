//! Combinations and the Heaviside forms (CAS-DOC-36 §8, CAS-402): the rules
//! of `Combine`, giac's `combine`, with `exp`, and of `SignToHeaviside` and
//! `HeavisideToSign`, giac's `sign2Heaviside` and `Heaviside2sign`; and
//! `TranscendentalSimplify`, giac's `tsimplify`.
//!
//! - `Combine[e, exp]` merges the exponentials of a product, `e^a e^b =
//!   e^(a + b)`, and their integer powers, `(e^a)^k = e^(k a)`.
//! - `sign u = 2θ(u) − 1` and `θ(u) = (1 + sign u)/2`, for every real `u`
//!   as `θ(0) = 1/2` (CAS-DOC-18 §2).
//! - `TranscendentalSimplify`: `Together` of `ExpExpand` of
//!   `HyperbolicToExp` of `TrigToExp`, each evaluated, which writes the
//!   expression with as few exponentials as it can:
//!   `exp(2*x)+2*exp(x)+1` is `exp(x)^2+2*exp(x)+1`, and
//!   `sin(x)^2+cos(x)^2` is 1.

use super::exponential::{exp_expand, hyperbolic};
use super::{apply, head};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::trig::to_exp;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

/// `f` as the exponent of `Exp[u]` or `E^u`, `u`, or of `Exp[u]^k`, `k u`
/// for an integer `k`; with whether it was such a power.
fn exponent_of(f: &Expr) -> Option<(Expr, bool)> {
    match (head(f), f.args()) {
        (Some(Symbol::EXP), [u]) => Some((u.clone(), false)),
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => {
            Some((u.clone(), false))
        }
        (Some(Symbol::POWER), [b, k]) if k.to_integer().is_some() => match (head(b), b.args()) {
            (Some(Symbol::EXP), [u]) => Some((k.clone() * u.clone(), true)),
            _ => None,
        },
        _ => None,
    }
}

/// `Combine[e, exp]`'s rule.
pub(super) fn exp_combine(_: &mut Eval<'_>, part: &Expr) -> Found {
    let exp = |u: Expr| Expr::call(Symbol::EXP, [u]);
    Ok(match head(part) {
        Some(Symbol::POWER) => match exponent_of(part) {
            Some((u, true)) => Some(exp(u)),
            _ => None,
        },
        Some(Symbol::TIMES) => {
            let mut exponents = Vec::new();
            let mut rest = Vec::new();
            for f in part.args() {
                match exponent_of(f) {
                    Some((u, _)) => exponents.push(u),
                    None => rest.push(f.clone()),
                }
            }
            if exponents.len() < 2 {
                return Ok(None);
            }
            rest.push(exp(Expr::call(Symbol::PLUS, exponents)));
            Some(Expr::call(Symbol::TIMES, rest))
        }
        _ => None,
    })
}

/// `SignToHeaviside`'s rule.
pub(super) fn sign_to_heaviside(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match (head(part), part.args()) {
        (Some(Symbol::SIGN), [u]) => {
            Some(Expr::from(2) * Expr::call(Symbol::HEAVISIDE_THETA, [u.clone()]) - Expr::from(1))
        }
        _ => None,
    })
}

/// `HeavisideToSign`'s rule.
pub(super) fn heaviside_to_sign(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match (head(part), part.args()) {
        (Some(Symbol::HEAVISIDE_THETA), [u]) => {
            Some((Expr::from(1) + Expr::call(Symbol::SIGN, [u.clone()])) / Expr::from(2))
        }
        _ => None,
    })
}

/// `TranscendentalSimplify` of `x`, not a list or a relation, built.
pub(super) fn transcendental_simplify(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let exponential = ev.eval(&to_exp(x))?;
    let logarithmic = apply(ev, &exponential, hyperbolic)?;
    let logarithmic = ev.eval(&logarithmic)?;
    let expanded = apply(ev, &logarithmic, exp_expand)?;
    Ok(Expr::call(Symbol::TOGETHER, [expanded]))
}

#[cfg(test)]
mod tests;
