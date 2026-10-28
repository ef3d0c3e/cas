//! `ExpIntegralEi`, `SinIntegral`, `CosIntegral` and `LogIntegral`
//! (CAS-DOC-31 §3.1, CAS-354), Xcas's `Ei`, `Si`, `Ci` and `Li`.
//!
//! - `Ei[0]` and `Ci[0]` are `-Infinity`, with the warning of `Log[0]`;
//!   `Si[0] = 0`.
//! - At the infinities: `Ei[Infinity] = Infinity`, `Ei[-Infinity] = 0`,
//!   `Si[±Infinity] = ±Pi/2` and `Ci[Infinity] = 0`.
//! - `Si` is odd: the sign of a negative number or coefficient comes out,
//!   as `Erf`'s does (CAS-DOC-15 §8).
//! - A real ball by `Ball::ei`, `Ball::si` and `Ball::ci`. `Ci` of a
//!   negative real ball is `Ci[|x|] + I Pi`, giac's principal value
//!   (`Ci(-1.0)` is `0.3374…+3.1416…*i`). A ball meeting the pole 0 is
//!   `Indeterminate`, with a warning.
//! - `LogIntegral[x]` is `ExpIntegralEi[Log[x]]`, as giac writes `Li(x)`.
//! - Anything else stays, complex balls among them.

use cas_num::ball::{Ball, ComplexBall};

use super::number::{Infinity, Num, Precision};
use super::special::negated;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn half_pi() -> Expr {
    Expr::from(Symbol::PI) / 2
}

fn minus_infinity() -> Expr {
    Infinity::Directed(Num::one().negated()).to_expr()
}

/// `ExpIntegralEi[x]`, `SinIntegral[x]` and `CosIntegral[x]`.
pub(super) fn integral(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(f), [x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    if *x == Expr::from(0) {
        if f == Symbol::SIN_INTEGRAL {
            return Ok(Some(Expr::from(0)));
        }
        ev.warn(Warning::Infinite(e.clone()));
        return Ok(Some(minus_infinity()));
    }
    if let Some(Infinity::Directed(d)) = Infinity::of(x) {
        let up = d.is_one();
        if !up && !d.negated().is_one() {
            return Ok(None);
        }
        return Ok(match (f, up) {
            (Symbol::EXP_INTEGRAL_EI, true) => Some(Expr::from(Symbol::INFINITY)),
            (Symbol::EXP_INTEGRAL_EI, false) | (Symbol::COS_INTEGRAL, true) => Some(Expr::from(0)),
            (Symbol::SIN_INTEGRAL, true) => Some(half_pi()),
            (Symbol::SIN_INTEGRAL, false) => Some(-half_pi()),
            _ => None,
        });
    }
    if let Some(ref number @ Num::Ball(ref ball)) = Num::of(x) {
        let bits = number.precision(Precision::of(ev.context()));
        let p = bits + 16;
        let value = match f {
            Symbol::EXP_INTEGRAL_EI => ball.ei(p),
            Symbol::SIN_INTEGRAL => ball.si(p),
            _ if ball.mid().is_negative() => {
                // Ci(x) = Ci(|x|) + iπ for x < 0.
                let re = ball.abs().ci(p).round_to_accuracy(bits);
                if !re.is_finite() {
                    ev.warn(Warning::Indeterminate(e.clone()));
                    return Ok(Some(Symbol::INDETERMINATE.into()));
                }
                let im = Ball::const_pi(p).round_to_accuracy(bits);
                return Ok(Some(
                    Num::ComplexBall(Box::new(ComplexBall::new(re, im))).to_expr(),
                ));
            }
            _ => ball.ci(p),
        };
        if !value.is_finite() {
            ev.warn(Warning::Indeterminate(e.clone()));
            return Ok(Some(Symbol::INDETERMINATE.into()));
        }
        return Ok(Some(
            Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr(),
        ));
    }
    if f == Symbol::SIN_INTEGRAL
        && let Some(y) = negated(x)
    {
        return Ok(Some(-Expr::call(f, [y])));
    }
    Ok(None)
}

/// `LogIntegral[x]`: `ExpIntegralEi[Log[x]]`.
pub(super) fn log_integral(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(Expr::call(
        Symbol::EXP_INTEGRAL_EI,
        [Expr::call(Symbol::LOG, [x.clone()])],
    )))
}

#[cfg(test)]
mod tests;
