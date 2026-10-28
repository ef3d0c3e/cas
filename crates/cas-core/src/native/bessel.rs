//! The Bessel and Airy functions (CAS-DOC-34 §3, CAS-382): `BesselJ`,
//! `BesselY`, `BesselI` and `BesselK`, Xcas's `BesselJ(n, x)` and the like,
//! the order first; `AiryAi` and `AiryBi`, Xcas's `Airy_Ai` and `Airy_Bi`,
//! and their derivatives `AiryAiPrime` and `AiryBiPrime`.
//!
//! - **Orders**: an integer one. A negative order comes out: `J_(−n) =
//!   (−1)^n J_n`, `Y_(−n) = (−1)^n Y_n`, `I_(−n) = I_n` and `K_(−n) = K_n`.
//!   Another order stays.
//! - **At 0**: `J_0(0) = I_0(0) = 1` and `J_n(0) = I_n(0) = 0` for `n ≥
//!   1`; `Y_n(0)` and `K_n(0)` stay, where giac errs. Airy's functions are
//!   their closed forms: `Ai(0) = 3^(−2/3)/Γ(2/3)`, `Bi(0) =
//!   3^(−1/6)/Γ(2/3)`, `Ai′(0) = −3^(−1/3)/Γ(1/3)` and `Bi′(0) =
//!   3^(1/6)/Γ(1/3)`.
//! - **Negative numbers** come out of `J` and `I`: `J_n(−x) = (−1)^n
//!   J_n(x)` and `I_n(−x) = (−1)^n I_n(x)`, as `SinIntegral`'s sign does.
//! - **A real ball** by `Ball::bessel_j` and its kin, and `Ball::airy`.
//!   `Y` and `K` of a ball not shown above 0 stay, as their values are
//!   complex below it; an indeterminate value is `Indeterminate`, with a
//!   warning.
//! - Anything else stays: an exact argument is a number only under `N`.

use std::cmp::Ordering;

use cas_num::ball::{Airy, Ball};

use super::number::{Num, Precision};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `Power[3, e]` for a rational `e = p/q`.
fn power_of_three(p: i64, q: i64) -> Expr {
    Expr::call(
        Symbol::POWER,
        [Expr::from(3), Expr::from(p) / Expr::from(q)],
    )
}

/// `Gamma[1/3]` or `Gamma[2/3]`.
fn gamma_third(p: i64) -> Expr {
    Expr::call(Symbol::GAMMA, [Expr::from(p) / Expr::from(3)])
}

/// The kind of Bessel function.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    J,
    Y,
    I,
    K,
}

impl Kind {
    fn of(head: Symbol) -> Option<Kind> {
        Some(match head {
            Symbol::BESSEL_J => Kind::J,
            Symbol::BESSEL_Y => Kind::Y,
            Symbol::BESSEL_I => Kind::I,
            Symbol::BESSEL_K => Kind::K,
            _ => return None,
        })
    }

    /// Whether `C_(−n) = (−1)^n C_n`, else `C_(−n) = C_n`.
    fn odd_in_order(self) -> bool {
        matches!(self, Kind::J | Kind::Y)
    }

    fn of_ball(self, x: &Ball, n: i64, prec: u64) -> Ball {
        match self {
            Kind::J => x.bessel_j(n, prec),
            Kind::Y => x.bessel_y(n, prec),
            Kind::I => x.bessel_i(n, prec),
            Kind::K => x.bessel_k(n, prec),
        }
    }
}

/// `C` with a sign `(−1)^n` for an odd `n`.
fn signed(odd: bool, value: Expr) -> Expr {
    if odd { -value } else { value }
}

/// `BesselJ[n, x]`, `BesselY[n, x]`, `BesselI[n, x]` and `BesselK[n, x]`.
pub(super) fn bessel(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(head), [n, x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let Some(kind) = Kind::of(head) else {
        return Ok(None);
    };
    let Some(order) = n.to_integer() else {
        return Ok(None);
    };
    let odd = order.unsigned_abs().bit(0);
    if order.is_negative() {
        let positive = Expr::call(head, [Expr::from(-order), x.clone()]);
        return Ok(Some(signed(odd && kind.odd_in_order(), positive)));
    }
    if *x == Expr::from(0) {
        return Ok(match kind {
            Kind::J | Kind::I => Some(Expr::from(i64::from(order.is_zero()))),
            Kind::Y | Kind::K => None,
        });
    }
    let Some(number) = Num::of(x) else {
        return Ok(None);
    };
    match number {
        Num::Real(ref q) if q.is_negative() && matches!(kind, Kind::J | Kind::I) => {
            let mirrored = Expr::call(head, [n.clone(), Expr::from(-q.clone())]);
            Ok(Some(signed(odd, mirrored)))
        }
        Num::Ball(ref ball) => {
            let Some(m) = order.to_i64() else {
                return Ok(None);
            };
            // Y and K below 0 are complex.
            if matches!(kind, Kind::Y | Kind::K) && ball.sign() != Some(Ordering::Greater) {
                return Ok(None);
            }
            let bits = number.precision(Precision::of(ev.context()));
            let value = kind.of_ball(ball, m, bits + 16);
            if !value.is_finite() {
                ev.warn(Warning::Indeterminate(e.clone()));
                return Ok(Some(Symbol::INDETERMINATE.into()));
            }
            Ok(Some(
                Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr(),
            ))
        }
        _ => Ok(None),
    }
}

/// `AiryAi[x]`, `AiryBi[x]`, `AiryAiPrime[x]` and `AiryBiPrime[x]`.
pub(super) fn airy(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(head), [x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let which = match head {
        Symbol::AIRY_AI => Airy::Ai,
        Symbol::AIRY_AI_PRIME => Airy::AiPrime,
        Symbol::AIRY_BI => Airy::Bi,
        Symbol::AIRY_BI_PRIME => Airy::BiPrime,
        _ => return Ok(None),
    };
    if *x == Expr::from(0) {
        return Ok(Some(match which {
            Airy::Ai => power_of_three(-2, 3) / gamma_third(2),
            Airy::Bi => power_of_three(-1, 6) / gamma_third(2),
            Airy::AiPrime => -(power_of_three(-1, 3) / gamma_third(1)),
            Airy::BiPrime => power_of_three(1, 6) / gamma_third(1),
        }));
    }
    let Some(ref number @ Num::Ball(ref ball)) = Num::of(x) else {
        return Ok(None);
    };
    let bits = number.precision(Precision::of(ev.context()));
    let value = ball.airy(which, bits + 16);
    if !value.is_finite() {
        ev.warn(Warning::Indeterminate(e.clone()));
        return Ok(Some(Symbol::INDETERMINATE.into()));
    }
    Ok(Some(
        Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr(),
    ))
}

#[cfg(test)]
mod tests;
