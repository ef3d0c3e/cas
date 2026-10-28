//! The Heaviside step and the Dirac delta (CAS-DOC-18 §2, CAS-273), Xcas's
//! `Heaviside` and `Dirac`.
//!
//! - `HeavisideTheta[x]`: 0 for `x < 0`, 1 for `x > 0`, `1/2` at 0 as
//!   giac's, the sign from the assumptions, a ball or the approximation, as
//!   `Sign` takes it; otherwise it stays.
//! - `DiracDelta[x]`: 0 for `x ≠ 0`; `ComplexInfinity` at 0, with the
//!   warning of `1/0` (giac's `infinity`); otherwise it stays.
//!
//! Their calculus is elsewhere: `D` of the step is the delta (`calculus`),
//! and `Integrate` has both at linear arguments and sifts the delta out of
//! definite integrals (`integrate`).

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `HeavisideTheta[x]`.
pub(super) fn heaviside(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(real_sign(ev.context(), x).map(|s| match s {
        Ordering::Greater => Expr::from(1),
        Ordering::Less => Expr::from(0),
        Ordering::Equal => Expr::from(Rational::new(Integer::ONE, Integer::from(2))),
    }))
}

/// `DiracDelta[x]`.
pub(super) fn dirac(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    match real_sign(ev.context(), x) {
        Some(Ordering::Equal) => {
            ev.warn(Warning::Infinite(e.clone()));
            Ok(Some(Symbol::COMPLEX_INFINITY.into()))
        }
        Some(_) => Ok(Some(Expr::from(0))),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests;
