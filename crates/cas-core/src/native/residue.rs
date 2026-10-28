//! `Residue` (CAS-DOC-25 §4, CAS-315): Xcas's `residue(f, z, a)` and
//! `residue(f, z = a)`, and Mathematica's `Residue[f, {z, a}]`, the
//! coefficient of `(z − a)⁻¹` in the Laurent series of `f` at `a`, from the
//! series engine (CAS-DOC-10) to the orders of [`ORDERS`] in turn: 0 where
//! `f` is analytic. At
//! infinity, `Infinity` with any direction or `ComplexInfinity`, it is
//! `−Res_{t=0} f(1/t)/t²`, Mathematica's sign: `Residue[1/z, z, Infinity]`
//! is −1, where giac answers 1.
//!
//! The call stays where there is no series (an essential singularity,
//! `e^(1/z)` at 0, which giac refuses too; an unknown head), at a branch
//! point (a ramification above 1, or `log t` in the coefficients: giac
//! answers `√z/z` at 0 with 0), or where a coefficient could not be told
//! from 0 (§3). Each call that answers takes one step, `"Residue"`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of complex analysis: f in z at a, the series s in t"
)]

use cas_num::{Integer, Rational};

use super::series::{expansion, fresh, log_t};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The orders tried, each until a series shows: a denominator whose
/// constant term cancels, `e^t − 1` or `1 + (i + t)²`, needs its terms past
/// the one that cancels, as `Limit`'s leading terms.
const ORDERS: [i64; 4] = [2, 4, 8, 16];

/// Whether `a` is a point at infinity.
fn infinite(a: &Expr) -> bool {
    a.contains(&Expr::from(Symbol::INFINITY))
        || a.contains(&Expr::from(Symbol::DIRECTED_INFINITY))
        || a.contains(&Expr::from(Symbol::COMPLEX_INFINITY))
}

/// `f`, `z` and the point of `Residue[f, z, a]`, `Residue[f, z == a]` or
/// `Residue[f, {z, a}]`, the variable a symbol; `None` for the point at
/// infinity.
fn parts(e: &Expr) -> Option<(&Expr, &Expr, Option<&Expr>)> {
    let (f, z, a) = match e.args() {
        [f, z, a] => (f, z, a),
        [f, s] => match (s.head().and_then(Expr::as_symbol), s.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL | Symbol::RULE | Symbol::LIST), [z, a]) => {
                (f, z, a)
            }
            _ => return None,
        },
        _ => return None,
    };
    z.as_symbol()?;
    Some((f, z, (!infinite(a)).then_some(a)))
}

pub(super) fn residue(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, z, a)) = parts(e) else {
        return Ok(None);
    };
    let t = fresh();
    if f.contains(&t) {
        return Ok(None);
    }
    let local = match a {
        Some(a) => f.substitute(&[(z.clone(), a.clone() + t.clone())]),
        None => -(f.substitute(&[(z.clone(), t.pow(-1))]) * t.pow(-2)),
    };
    let local = ev.eval(&local)?;
    let mut series = None;
    for order in ORDERS {
        series = expansion(ev, &local, &t, &Rational::from(Integer::from(order)))?;
        if series.is_some() {
            break;
        }
    }
    let Some(s) = series else {
        return Ok(None);
    };
    if s.ramification != 1 || s.coeffs.iter().any(|c| c.contains(&log_t())) {
        return Ok(None);
    }
    // The exponent of the k-th coefficient is the valuation plus k.
    let value = usize::try_from(-1 - s.valuation)
        .ok()
        .and_then(|k| s.coeffs.get(k).cloned())
        .unwrap_or_else(|| Expr::from(0));
    ev.step("Residue", e, &value);
    Ok(Some(value))
}

#[cfg(test)]
mod tests;
