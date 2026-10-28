//! Pythagorean forms (CAS-DOC-36 §4, CAS-398): `TrigToCos`, `TrigToSin`
//! and `TrigToTan`, giac's `trigcos`, `trigsin` and `trigtan`, by `sin² +
//! cos² = 1` and `tan = sin/cos`, each favouring one function:
//!
//! - `TrigToCos`: `Tan[u]` as `Sin[u]/Cos[u]`, made `Together`, so that
//!   numerator and denominator are expanded; then `Sin[u]^n`, `|n| ≥ 2`, as
//!   `Sin[u]^(n mod 2) (1 − Cos[u]²)^(n div 2)`, the form `Simplify` tries
//!   (CAS-DOC-14 §8); then `Together` again. giac's
//!   `trigcos(sin(x)^4+sin(x)^2)` is `cos(x)^4-3*cos(x)^2+2`.
//! - `TrigToSin`: the same, `Sin` and `Cos` the other way.
//! - `TrigToTan`: `Sin[u]` as `Cos[u] Tan[u]`, made `Together`; then
//!   `Cos[u]^n`, `|n| ≥ 2`, as `Cos[u]^(n mod 2) (1 + Tan[u]²)^(−(n div
//!   2))`; then `Together` again. giac's `trigtan(sin(x)*cos(x))` is
//!   `tan(x)/(tan(x)^2+1)`.
//!
//! Lists and relations are rewritten entry by entry, as `Together` takes
//! one expression (`entrywise`).

use super::head;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::trig::{pythagorean, tan_out};
use crate::symbol::Symbol;

/// The function a form favours.
#[derive(Clone, Copy)]
pub(super) enum Favoured {
    Cos,
    Sin,
    Tan,
}

/// `e` with `Sin[u]` as `Cos[u] Tan[u]`, built.
fn sin_out(e: &Expr) -> Expr {
    e.map_bottom_up(|part| match (head(part), part.args()) {
        (Some(Symbol::SIN), [u]) => {
            Expr::call(Symbol::COS, [u.clone()]) * Expr::call(Symbol::TAN, [u.clone()])
        }
        _ => part.clone(),
    })
}

/// `e` with `Cos[u]^n`, `|n| ≥ 2`, as `Cos[u]^(n mod 2) (1 +
/// Tan[u]²)^(−(n div 2))`, built.
fn cos_powers(e: &Expr) -> Expr {
    e.map_bottom_up(|part| {
        let (Some(Symbol::POWER), [b, n]) = (head(part), part.args()) else {
            return part.clone();
        };
        let (Some(Symbol::COS), [u]) = (head(b), b.args()) else {
            return part.clone();
        };
        let Some(n) = n.as_i64().filter(|n| n.abs() >= 2) else {
            return part.clone();
        };
        let square = Expr::from(1) + Expr::call(Symbol::TAN, [u.clone()]).pow(2);
        b.pow(n.rem_euclid(2)) * square.pow(-n.div_euclid(2))
    })
}

/// `x`, not a list or a relation, in the form that favours `favoured`:
/// `Together` of it, built.
pub(super) fn form(ev: &mut Eval<'_>, x: &Expr, favoured: Favoured) -> Result<Expr, EvalError> {
    let first = match favoured {
        Favoured::Cos | Favoured::Sin => tan_out(x),
        Favoured::Tan => sin_out(x),
    };
    let first = ev.eval(&Expr::call(Symbol::TOGETHER, [first]))?;
    let second = match favoured {
        Favoured::Cos => pythagorean(&first, Symbol::SIN, Symbol::COS),
        Favoured::Sin => pythagorean(&first, Symbol::COS, Symbol::SIN),
        Favoured::Tan => cos_powers(&first),
    };
    Ok(Expr::call(Symbol::TOGETHER, [second]))
}

#[cfg(test)]
mod tests;
