//! `Zeta` (CAS-DOC-15 §7, CAS-245), Xcas's `Zeta`: the Riemann zeta
//! function.
//!
//! - At the integers: `ζ(2m) = (−1)^(m+1) B_{2m} 2^(2m−1) π^(2m)/(2m)!`
//!   (`Zeta[2] = Pi^2/6`), `ζ(−n) = −B_{n+1}/(n + 1)` (0 at the negative
//!   even integers), `ζ(0) = −1/2`, and `ComplexInfinity` at the pole 1
//!   with the warning of `1/0`, where giac answers `+infinity`; the odd
//!   positive integers stay, as giac's `Zeta(3)`, and so do the Bernoulli
//!   numbers past cas-num's table.
//! - A real ball by `Ball::zeta`, `Indeterminate` with a warning where it
//!   meets the pole.
//! - Anything else stays.

use cas_num::ball::bernoulli;
use cas_num::{Integer, Natural, Rational};

use super::number::{Num, Precision};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `Zeta[s]`.
pub(super) fn zeta(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [s] = e.args() else {
        return Ok(None);
    };
    if let Some(n) = s.as_i64() {
        return Ok(at_integer(ev, e, n));
    }
    let Some(number) = Num::of(s) else {
        return Ok(None);
    };
    let Num::Ball(ball) = &number else {
        return Ok(None);
    };
    let bits = number.precision(Precision::of(ev.context()));
    let value = ball.zeta(bits + 16).round_to_accuracy(bits);
    if !value.is_finite() {
        ev.warn(Warning::Indeterminate(e.clone()));
        return Ok(Some(Symbol::INDETERMINATE.into()));
    }
    Ok(Some(Num::Ball(Box::new(value)).to_expr()))
}

/// `ζ(n)` at an integer, where it is exact.
fn at_integer(ev: &mut Eval<'_>, e: &Expr, n: i64) -> Option<Expr> {
    match n {
        1 => {
            ev.warn(Warning::Infinite(e.clone()));
            Some(Symbol::COMPLEX_INFINITY.into())
        }
        0 => Some(Expr::from(Rational::new(
            Integer::from(-1),
            Integer::from(2),
        ))),
        _ if n < 0 => {
            // ζ(−m) = −B_{m+1}/(m + 1).
            let k = n.unsigned_abs() + 1;
            let b = bernoulli(k)?;
            Some(Expr::from(-(&b / &Rational::from(Integer::from(k)))))
        }
        _ if n % 2 == 0 => {
            let m = n.unsigned_abs();
            let b = bernoulli(m)?;
            // (−1)^(m/2 + 1) B_m 2^(m−1) / m!, times π^m.
            let scale = Rational::new(
                Integer::from(2).pow(m - 1),
                Integer::from(Natural::factorial(m)),
            );
            let mut coefficient = &b * &scale;
            if (m / 2) % 2 == 0 {
                coefficient = -coefficient;
            }
            Some(Expr::from(coefficient) * Expr::from(Symbol::PI).pow(n))
        }
        _ => None,
    }
}
