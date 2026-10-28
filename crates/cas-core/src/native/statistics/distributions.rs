//! Probability distributions (CAS-DOC-23): Xcas's discrete and continuous
//! distributions, each with its mass or density, its cumulative
//! probability of a point or an interval `P(a ≤ X ≤ b)`, and its inverse.
//!
//! - The discrete ones ([`discrete`], §2, CAS-301): binomial, Poisson,
//!   geometric and negative binomial.
//! - The monotone inversion of continuous cumulative functions, and
//!   `InverseErf` ([`inversion`], §3, §5, CAS-302).
//! - The normal, exponential and uniform distributions ([`continuous`], §3,
//!   CAS-303).
//! - Student's t, χ² and Fisher's F ([`degrees`], §4, CAS-304), through the
//!   closed forms of the regularized incomplete gamma and beta functions
//!   at integer and half-integer parameters ([`incomplete`]).
//!
//! A parameter that is a number must be valid (a probability in `[0, 1]`,
//! a rate above 0, a natural count), as comparisons decide; otherwise the
//! call stays, with a warning where giac answers `undef`. A symbolic
//! parameter leaves the masses and the finite cumulative sums as
//! expressions, and the inverses as calls.

pub(in crate::native) mod continuous;
pub(in crate::native) mod degrees;
pub(in crate::native) mod discrete;
mod incomplete;
pub(in crate::native) mod inversion;

use std::cmp::Ordering;

use super::super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::trace::Warning;

/// Whether a parameter fits, as far as comparisons tell.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Fit {
    /// A number that fits.
    Number,
    /// Not a number: a symbol, or a number comparisons do not decide.
    Symbolic,
    /// A number that does not fit.
    Wrong,
}

/// How `x` compares with `y`, as a number; `None` when not decided.
fn compare(ev: &mut Eval<'_>, x: &Expr, y: &Expr) -> Result<Option<Ordering>, EvalError> {
    if let (Some(a), Some(b)) = (x.to_rational(), y.to_rational()) {
        return Ok(Some(a.cmp(&b)));
    }
    let difference = ev.eval(&(x - y))?;
    Ok(real_sign(ev.context(), &difference))
}

/// Whether `x` is a probability, in `[0, 1]`.
fn probability(ev: &mut Eval<'_>, x: &Expr) -> Result<Fit, EvalError> {
    if !x.is_number() && !numeric(x) {
        return Ok(Fit::Symbolic);
    }
    let low = compare(ev, x, &Expr::from(0))?;
    let high = compare(ev, x, &Expr::from(1))?;
    Ok(match (low, high) {
        (Some(Ordering::Greater | Ordering::Equal), Some(Ordering::Less | Ordering::Equal)) => {
            Fit::Number
        }
        (Some(_), Some(_)) => Fit::Wrong,
        _ => Fit::Symbolic,
    })
}

/// Whether `x` is a number above 0.
fn positive(ev: &mut Eval<'_>, x: &Expr) -> Result<Fit, EvalError> {
    if !x.is_number() && !numeric(x) {
        return Ok(Fit::Symbolic);
    }
    Ok(match compare(ev, x, &Expr::from(0))? {
        Some(Ordering::Greater) => Fit::Number,
        Some(_) => Fit::Wrong,
        None => Fit::Symbolic,
    })
}

/// Whether `x` is free of symbols that stand for values: an exact
/// numeric quantity such as `√2` or `π/4`, which comparisons may decide.
fn numeric(x: &Expr) -> bool {
    x.symbols().iter().all(|s| s.is_builtin())
}

/// The call stays, with a warning: a parameter giac answers `undef` for.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Result<Option<Expr>, EvalError> {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

#[cfg(test)]
mod tests;
