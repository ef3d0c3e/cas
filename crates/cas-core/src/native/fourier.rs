//! Fourier coefficients (CAS-DOC-25 §7, CAS-319): Xcas's `fourier_an`,
//! `fourier_bn` and `fourier_cn`, `FourierCosCoefficient`,
//! `FourierSinCoefficient` and `FourierCoefficient` here, of `f` in `x` of
//! period `T` over `[a, a + T]`, `a` 0 when left out:
//!
//! - `aₙ = (2/T) ∫ f(x) cos(2πnx/T) dx`, and `a₀ = (1/T) ∫ f`, the mean;
//! - `bₙ = (2/T) ∫ f(x) sin(2πnx/T) dx`, and `b₀ = 0`;
//! - `cₙ = (1/T) ∫ f(x) e^(−2iπnx/T) dx`, as `(aₙ − i bₙ)/2` from the same
//!   real kernels, so that products of sines and cosines linearise
//!   (CAS-199), and `c₀ = a₀`.
//!
//! The integrands `f` times a kernel are expanded, so that a sum integrates
//! term by term.
//!
//! A numeric `n` is used as it is, an integer or not, as giac does. A
//! symbolic `n` is taken as an integer: the integrals are evaluated in a
//! scope that assumes `n ∈ ℤ`, so that `sin(nπ)` is 0 and `cos(nπ)` is
//! `(−1)ⁿ` (CAS-313). Then the arguments of the sines, cosines and
//! exponentials of the result are expanded, so that `sin(2πn(a + T)/T)`
//! shows its `2πn`, and the result is evaluated there again and put
//! together. The answer holds for the `n` where the antiderivative does,
//! as giac's: `fourier_an(cos(x)^2, x, 2π, n, −π)` is 0, but `1/2` at
//! `n = 2`. A period that is not a number is assumed positive, so that the
//! breakpoints of `abs` and `sign` are ordered (CAS-316). A period 0, or an
//! integral that stays, leaves the call, with the integral's warnings.
//! Each call that answers takes one step, `"Fourier"`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Fourier series: f in x of period T from a, n"
)]

use cas_num::{Integer, Rational};

use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// Which coefficient.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    /// `aₙ`, against the cosine.
    Cosine,
    /// `bₙ`, against the sine.
    Sine,
    /// `cₙ`, against `e^(−iωnx)`.
    Exponential,
}

pub(super) fn fourier_cos_coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    coefficient(ev, e, Kind::Cosine)
}

pub(super) fn fourier_sin_coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    coefficient(ev, e, Kind::Sine)
}

pub(super) fn fourier_coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    coefficient(ev, e, Kind::Exponential)
}

/// `f`, `x`, `T`, `n` and `a` of a call, `x` a symbol and the others free
/// of it.
fn parts(e: &Expr) -> Option<[Expr; 5]> {
    let (f, x, period, n, a) = match e.args() {
        [f, x, period, n] => (f, x, period, n, Expr::from(0)),
        [f, x, period, n, a] => (f, x, period, n, a.clone()),
        _ => return None,
    };
    x.as_symbol()?;
    if [period, n, &a].into_iter().any(|g| g.contains(x)) {
        return None;
    }
    Some([f.clone(), x.clone(), period.clone(), n.clone(), a])
}

fn coefficient(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some([f, x, period, n, a]) = parts(e) else {
        return Ok(None);
    };
    if period == Expr::from(0) {
        return Ok(None);
    }
    let symbolic = !n.is_number();
    let mut assumptions = Vec::new();
    if symbolic {
        assumptions.push(Expr::call(
            Symbol::ELEMENT,
            [n.clone(), Expr::from(Symbol::INTEGERS)],
        ));
    }
    if !period.is_number() {
        assumptions.push(Expr::call(Symbol::GREATER, [period.clone(), Expr::from(0)]));
    }
    let mut scope: Option<Context> = None;
    for assumption in assumptions {
        scope = Some(scope.as_ref().unwrap_or(ev.context()).assuming(assumption));
    }
    // Expanded, so that a sum times the kernel integrates term by term.
    let integral = |kernel: Expr| {
        Expr::call(
            Symbol::INTEGRATE,
            [
                Expr::call(Symbol::EXPAND, [f.clone() * kernel]),
                x.clone(),
                a.clone(),
                a.clone() + period.clone(),
            ],
        )
    };
    let half = || Expr::from(Rational::new(Integer::from(1), Integer::from(2)));
    let value = if n == Expr::from(0) {
        match kind {
            Kind::Sine => Expr::from(0),
            _ => integral(Expr::from(1)) / period.clone(),
        }
    } else {
        let w = Expr::from(2) * Expr::from(Symbol::PI) * n.clone() * x.clone() / period.clone();
        let (cos, sin) = (
            integral(Expr::call(Symbol::COS, [w.clone()])),
            integral(Expr::call(Symbol::SIN, [w])),
        );
        let scale = Expr::from(2) / period.clone();
        match kind {
            Kind::Cosine => scale * cos,
            Kind::Sine => scale * sin,
            Kind::Exponential => half() * scale * (cos - Expr::from(Symbol::I) * sin),
        }
    };
    let Some(integrated) = evaluate(ev, scope.as_ref(), &value)? else {
        return Ok(None);
    };
    let together = Expr::call(Symbol::TOGETHER, [expanded_arguments(&integrated)]);
    let Some(result) = evaluate(ev, scope.as_ref(), &together)? else {
        return Ok(None);
    };
    ev.step("Fourier", e, &result);
    Ok(Some(result))
}

/// `value` evaluated in `scope`, or where the call is without one, the
/// warnings passed on; `None` when an integral stays.
fn evaluate(ev: &mut Eval<'_>, scope: Option<&Context>, value: &Expr) -> Native {
    let result = match scope {
        Some(scope) => {
            let (result, warnings) = ev.eval_in(scope, value);
            for w in warnings {
                ev.warn(w);
            }
            result?
        }
        None => ev.eval(value)?,
    };
    Ok((!result.contains(&Expr::from(Symbol::INTEGRATE))).then_some(result))
}

/// `e` with the arguments of its sines, cosines, tangents and exponentials
/// under `Expand`, built.
fn expanded_arguments(e: &Expr) -> Expr {
    let expand = |u: &Expr| Expr::call(Symbol::EXPAND, [u.clone()]);
    e.map_bottom_up(
        |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(h @ (Symbol::SIN | Symbol::COS | Symbol::TAN | Symbol::EXP)), [u]) => {
                Expr::call(h, [expand(u)])
            }
            (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => {
                Expr::call(Symbol::POWER, [b.clone(), expand(u)])
            }
            _ => part.clone(),
        },
    )
}

#[cfg(test)]
mod tests;
