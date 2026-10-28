//! The classical orthogonal polynomials (CAS-DOC-27 §5, CAS-331): Xcas's
//! `legendre`, `hermite`, `laguerre`, `tchebyshev1` and `tchebyshev2`, as
//! `LegendreP[n, x]`, `HermiteH[n, x]` (the physicists'), `LaguerreL[n, x]`
//! and `LaguerreL[n, x, a]`, `ChebyshevT[n, x]` and `ChebyshevU[n, x]`.
//!
//! `LaguerreL[n, x, a]` is `L_n^(a)(x)` in giac's order, the variable
//! before the parameter: Mathematica's `LaguerreL[n, a, x]` swaps them.
//! Without `a`, `a = 0`, as in Mathematica, where giac keeps a symbol `a`.
//!
//! For an integer degree `0 ≤ n ≤` [`MAX_DEGREE`], by the three-term
//! recurrences:
//!
//! - `P_{n+1} = ((2n + 1) x P_n − n P_{n−1})/(n + 1)`;
//! - `H_{n+1} = 2x H_n − 2n H_{n−1}`;
//! - `T_{n+1} = 2x T_n − T_{n−1}`, and `U` the same from `U_1 = 2x`;
//! - `L^(a)_{n+1} = ((2n + 1 + a − x) L_n − (n + a) L_{n−1})/(n + 1)`.
//!
//! The recurrence runs on the values when the argument is a number: the
//! recurrences are stable on `[−1, 1]`, where the expanded polynomials of
//! large degree cancel. Otherwise each step is expanded, so that `P_3(x)`
//! is giac's `5/2*x^3-3/2*x`.
//!
//! Ball arithmetic does not see that stability: a radius grows by `1 + √2`
//! at each step of `P`, `T` and `U`. So an inexact real ball `X` inside
//! `[−1, 1]` takes their value at its midpoint `m`, at `1.3n` bits more,
//! and the rest by Markov's inequality `|p′| ≤ n² max_{[−1,1]} |p|`, the
//! maximum 1 for `P_n` and `T_n` and `n + 1` for `U_n`:
//! `p(X) ⊂ p(m) ± n² max |p| rad X`. Other balls take the recurrence as it
//! is, which may give a wide ball. A negative or
//! symbolic degree, or one past the limit, leaves the call. Each call that
//! answers takes one step, its builtin's name.

use cas_num::Integer;
use cas_num::ball::{Ball, Mag};

use super::number::Num;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest degree computed.
const MAX_DEGREE: i64 = 4096;

/// A family of the module.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Family {
    Legendre,
    Hermite,
    Laguerre,
    ChebyshevT,
    ChebyshevU,
}

pub(super) fn legendre_p(ev: &mut Eval<'_>, e: &Expr) -> Native {
    orthogonal(ev, e, Family::Legendre)
}

pub(super) fn hermite_h(ev: &mut Eval<'_>, e: &Expr) -> Native {
    orthogonal(ev, e, Family::Hermite)
}

pub(super) fn laguerre_l(ev: &mut Eval<'_>, e: &Expr) -> Native {
    orthogonal(ev, e, Family::Laguerre)
}

pub(super) fn chebyshev_t(ev: &mut Eval<'_>, e: &Expr) -> Native {
    orthogonal(ev, e, Family::ChebyshevT)
}

pub(super) fn chebyshev_u(ev: &mut Eval<'_>, e: &Expr) -> Native {
    orthogonal(ev, e, Family::ChebyshevU)
}

fn orthogonal(ev: &mut Eval<'_>, e: &Expr, family: Family) -> Native {
    let (n, a, x) = match (family, e.args()) {
        (Family::Laguerre, [n, x, a]) => (n, a.clone(), x),
        (_, [n, x]) => (n, Expr::from(0), x),
        _ => return Ok(None),
    };
    let Some(n) = n.as_i64().filter(|n| (0..=MAX_DEGREE).contains(n)) else {
        return Ok(None);
    };
    if let View::Ball(b) = x.view()
        && !b.is_exact()
        && family != Family::Hermite
        && family != Family::Laguerre
        && let Some(value) = markov(
            family,
            n,
            b,
            ev.context()
                .working_precision()
                .max(u64::try_from(b.rel_accuracy_bits()).unwrap_or(0)),
        )
    {
        let value = Expr::from(&value);
        ev.step(family.name(), e, &value);
        return Ok(Some(value));
    }
    let numeric = Num::of(x).is_some() && Num::of(&a).is_some();
    let step = |ev: &mut Eval<'_>, value: Expr| -> Result<Expr, EvalError> {
        if numeric {
            ev.eval(&value)
        } else {
            ev.eval(&Expr::call(Symbol::EXPAND, [value]))
        }
    };
    let first = match family {
        Family::Legendre | Family::ChebyshevT => x.clone(),
        Family::Hermite | Family::ChebyshevU => Expr::from(2) * x.clone(),
        Family::Laguerre => Expr::from(1) + a.clone() - x.clone(),
    };
    let mut previous = Expr::from(1);
    let mut current = step(ev, first)?;
    if n == 0 {
        current = previous.clone();
    }
    for k in 1..n {
        let kk = Expr::from(k);
        let next = match family {
            Family::Legendre => {
                (Expr::from(2 * k + 1) * x.clone() * current.clone()
                    - kk.clone() * previous.clone())
                    / Expr::from(k + 1)
            }
            Family::Hermite => {
                Expr::from(2) * x.clone() * current.clone() - Expr::from(2 * k) * previous.clone()
            }
            Family::ChebyshevT | Family::ChebyshevU => {
                Expr::from(2) * x.clone() * current.clone() - previous.clone()
            }
            Family::Laguerre => {
                ((Expr::from(2 * k + 1) + a.clone() - x.clone()) * current.clone()
                    - (kk + a.clone()) * previous.clone())
                    / Expr::from(k + 1)
            }
        };
        let next = step(ev, next)?;
        previous = std::mem::replace(&mut current, next);
    }
    ev.step(family.name(), e, &current);
    Ok(Some(current))
}

impl Family {
    fn name(self) -> &'static str {
        match self {
            Family::Legendre => "LegendreP",
            Family::Hermite => "HermiteH",
            Family::Laguerre => "LaguerreL",
            Family::ChebyshevT => "ChebyshevT",
            Family::ChebyshevU => "ChebyshevU",
        }
    }
}

/// `P_n`, `T_n` or `U_n` over the inexact ball `b` within `[−1, 1]` (the
/// module's Markov bound); `None` when `b` reaches past `[−1, 1]`.
fn markov(family: Family, n: i64, b: &Ball, prec: u64) -> Option<Ball> {
    let one = Ball::from_integer(&Integer::from(1));
    let reach = Mag::from_float(&b.mid().abs()).add(b.rad());
    if reach > Mag::from_integer(&Integer::from(1)) {
        return None;
    }
    let count = u64::try_from(n).ok()?;
    let wp = prec + 32 + (13 * count).div_ceil(10);
    let m = Ball::exact(b.mid().clone());
    let two_m = m.mul_2exp(&Integer::from(1));
    let integer = |k: i64| Ball::from_integer(&Integer::from(k));
    let mut previous = one.clone();
    let mut current = match family {
        Family::ChebyshevU => two_m.clone(),
        _ => m.clone(),
    };
    if n == 0 {
        current = one;
    }
    for k in 1..n {
        let next = match family {
            Family::Legendre => integer(2 * k + 1)
                .mul(&m, wp)
                .mul(&current, wp)
                .sub(&integer(k).mul(&previous, wp), wp)
                .div(&integer(k + 1), wp),
            _ => two_m.mul(&current, wp).sub(&previous, wp),
        };
        previous = std::mem::replace(&mut current, next);
    }
    let most = if family == Family::ChebyshevU {
        n + 1
    } else {
        1
    };
    let slope = Mag::from_integer(&Integer::from(n.checked_mul(n)?.checked_mul(most)?));
    Some(current.add_error(&slope.mul(b.rad())))
}

#[cfg(test)]
mod tests;
