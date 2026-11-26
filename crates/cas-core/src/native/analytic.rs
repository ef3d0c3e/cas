//! Analytic signals (CAS-635), Xcas's `hilbert` on a list, `instfreq`,
//! `instphase` and `emd`.
//!
//! - `HilbertMatrix[l]`, giac's `hilbert(l)` on a list of numbers (giac's
//!   `dht`, signalprocessing.cc:4957-4974): the analytic signal of the real
//!   parts, `X = fft(re l)`, `V_k = w_k X_k` for `k ≤ ⌊n/2⌋` with `w_k = 1`
//!   for `k = 0` and `k = ⌊(n+1)/2⌋`, else 2, `V_k = 0` above, the answer
//!   `ifft(V)`; the real parts are the samples, the imaginary parts their
//!   Hilbert transform. Through the FFT core of [`super::dft`] in complex
//!   balls at the working precision and its guard, so the lengths it takes;
//!   on exact samples an entry whose ball is exact is an exact number
//!   (`hilbert([1,0,0,0])` is `[1,i/2,0,−i/2]` where giac's doubles carry
//!   `3.06e-17`), else a ball. A nonzero imaginary part draws giac's warning
//!   "hilbert ignores the imaginary part of input". `HilbertMatrix[n]` stays
//!   the Hilbert matrix (permu.cc); giac's transform of an expression
//!   (`hilbert(sin(t),t)`, `hilbtrans`) is [`transform`]'s (CAS-648).
//! - `InstantaneousFrequency[z]`, giac's `instfreq` (:5583-5610), on a
//!   list: `ω_k = arg(z_(k+1) conj(z_k))`, `n − 1` values, built and
//!   evaluated, so exact on exact samples (`[π/2,π/2,π/2]` for
//!   `[1,i,−1,−i]`). `InstantaneousPhase[z]` (:5615-5646): `arg(z_0)`, then
//!   the sums with each `ω_k`, `n` values.
//! - `InstantaneousFrequency[u]` and `InstantaneousFrequency[u, x]` (`x` by
//!   default), `u` an expression: giac's `Im(conj(N) N')` with `N =
//!   u/|u|`, `x` real; with `u = a + i b` (`ComplexExpand`, every symbol
//!   real) that is `(a b' − b a')/(a² + b²)`, simplified.
//!   `InstantaneousPhase[u, x]`: `arg(u(0)) + ∫₀ˣ ω`, simplified (giac's
//!   `-2*pi*floor(1/2)+x^2` is our `x²` for `exp(i x²)`).
//! - `EmpiricalModeDecomposition[data, options…]`, giac's `emd` (:5131-5281,
//!   adapted by giac from Ivan Magrin-Chagnolleau's Matlab code): see
//!   [`emd`].
//!
//! Departures: giac's transform of an expression stays; `emd` of symbolic
//! data stays with giac's warning "Data must be numeric"; giac's errors
//! (a bad option, empty data) leave the call.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's and GSL's names: z, u, x, a, b, n, h, m"
)]

use std::borrow::Cow;

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Complex, Integer, Rational};

use super::dft;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

mod emd;
mod transform;

pub(super) use emd::empirical_mode_decomposition;
pub(super) use transform::hilbert_transform;

type Native = Result<Option<Expr>, EvalError>;

fn warn(ev: &mut Eval<'_>, e: &Expr, message: impl Into<Cow<'static, str>>) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

fn list(entries: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, entries)
}

/// The transform of `x` with the FFT core, or its inverse (scaled by
/// `1/n`); the length is checked by the caller.
fn transform(x: &[ComplexBall], inverse: bool, prec: u64) -> Vec<ComplexBall> {
    let n = x.len();
    let w = dft::roots(n, inverse, prec);
    let found = if n.is_power_of_two() {
        dft::radix_two(x, &w, prec)
    } else {
        dft::direct(x, &w, prec)
    };
    if !inverse {
        return found;
    }
    let scale = ComplexBall::from_real(Ball::from_rational(
        &Rational::new(
            Integer::ONE,
            Integer::from(i64::try_from(n).expect("a length")),
        ),
        prec,
    ));
    found.iter().map(|z| z.mul(&scale, prec)).collect()
}

/// An exact ball's value, the midpoint as a rational.
fn exact(b: &Ball) -> Option<Rational> {
    b.is_exact()
        .then(|| b.mid().checked_to_rational())
        .flatten()
}

/// `HilbertMatrix[l]` on a list of numbers, giac's `dht`.
pub(in crate::native) fn hilbert(ev: &mut Eval<'_>, e: &Expr, l: &Expr) -> Native {
    let n = l.args().len();
    if n == 0 {
        return Ok(None);
    }
    if !(n <= dft::MAX_DIRECT || (n.is_power_of_two() && n <= dft::MAX_RADIX_TWO)) {
        warn(
            ev,
            e,
            format!(
                "a length of {n}, not a power of 2, above {}",
                dft::MAX_DIRECT
            ),
        );
        return Ok(None);
    }
    let guard = 2 * u64::from(usize::BITS - n.leading_zeros()) + 16;
    let prec = ev.context().working_precision() + guard;
    let Some(x) = l
        .args()
        .iter()
        .map(|v| dft::ball(v, prec))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let exact_input = l.args().iter().all(|v| {
        matches!(
            v.view(),
            View::Int(_) | View::Integer(_) | View::Rational(_) | View::Complex(_)
        )
    });
    if x.iter()
        .any(|z| !(z.im().is_exact() && *z.im() == Ball::ZERO))
    {
        warn(ev, e, "hilbert ignores the imaginary part of input");
    }
    let real: Vec<ComplexBall> = x
        .iter()
        .map(|z| ComplexBall::from_real(z.re().clone()))
        .collect();
    let u = transform(&real, false, prec);
    let v: Vec<ComplexBall> = u
        .iter()
        .enumerate()
        .map(|(k, z)| {
            // giac's loop stops at n/2: for odd n, (n+1)/2 lies past it.
            if k > n / 2 {
                ComplexBall::ZERO
            } else if k == 0 || k == n.div_ceil(2) {
                z.clone()
            } else {
                z.add(z, prec)
            }
        })
        .collect();
    let found = transform(&v, true, prec);
    let entries = found
        .into_iter()
        .map(|z| {
            if exact_input && let (Some(re), Some(im)) = (exact(z.re()), exact(z.im())) {
                return Expr::from(Complex::new(re, im));
            }
            dft::number(z)
        })
        .collect();
    let result = list(entries);
    ev.step("hilbert", e, &result);
    Ok(Some(result))
}

/// `arg(b conj(a))`, unevaluated.
fn step_angle(a: &Expr, b: &Expr) -> Expr {
    Expr::call(
        Symbol::ARG,
        [Expr::call(
            Symbol::TIMES,
            [b.clone(), Expr::call(Symbol::CONJUGATE, [a.clone()])],
        )],
    )
}

/// The list forms: the frequencies, or the phases.
fn of_list(ev: &mut Eval<'_>, e: &Expr, z: &[Expr], phase: bool) -> Native {
    let n = z.len();
    let limit = ev.context().limits().terms;
    // A phase is a sum of up to n angles: n² terms unless every sample is
    // a ball, whose angles add up to one ball.
    let balls = z
        .iter()
        .all(|v| matches!(v.view(), View::Ball(_) | View::ComplexBall(_)));
    let size = if phase && !balls {
        n.checked_mul(n)
    } else {
        Some(n)
    };
    if size.is_none_or(|s| s > limit) {
        return Ok(None);
    }
    // Arg of a complex ball builds `atan(b/a) + Pi`: decimals answer balls.
    let inexact = z
        .iter()
        .any(|v| matches!(v.view(), View::Ball(_) | View::ComplexBall(_)));
    let angle = |ev: &mut Eval<'_>, a: Expr| {
        let a = ev.eval(&a)?;
        if inexact {
            ev.eval(&Expr::call(Symbol::N, [a]))
        } else {
            Ok(a)
        }
    };
    let mut out = Vec::with_capacity(n);
    if phase && n > 0 {
        out.push(angle(ev, Expr::call(Symbol::ARG, [z[0].clone()]))?);
    }
    for w in z.windows(2) {
        let omega = angle(ev, step_angle(&w[0], &w[1]))?;
        let value = match (phase, out.last()) {
            (true, Some(previous)) => {
                ev.eval(&Expr::call(Symbol::PLUS, [previous.clone(), omega]))?
            }
            _ => omega,
        };
        out.push(value);
    }
    let result = list(out);
    ev.step(if phase { "instphase" } else { "instfreq" }, e, &result);
    Ok(Some(result))
}

/// `ω = (a b' − b a')/(a² + b²)` for `u = a + i b`, `x` real.
fn frequency(ev: &mut Eval<'_>, u: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let part = |s: Symbol| Expr::call(Symbol::COMPLEX_EXPAND, [Expr::call(s, [u.clone()])]);
    let a = ev.eval(&part(Symbol::RE))?;
    let b = ev.eval(&part(Symbol::IM))?;
    let d = |f: &Expr| Expr::call(Symbol::D, [f.clone(), x.clone()]);
    let omega = (a.clone() * d(&b) - b.clone() * d(&a)) / (a.pow(2) + b.pow(2));
    ev.eval(&Expr::call(Symbol::SIMPLIFY, [omega]))
}

/// The expression forms: `ω`, or `arg(u(0)) + ∫₀ˣ ω`.
fn of_expression(ev: &mut Eval<'_>, e: &Expr, u: &Expr, x: &Expr, phase: bool) -> Native {
    let omega = frequency(ev, u, x)?;
    let result = if phase {
        let t = Expr::from(Symbol::new("InstantaneousPhase$t"));
        let at = |f: &Expr, v: &Expr| {
            Expr::call(
                Symbol::REPLACE_ALL,
                [f.clone(), Expr::call(Symbol::RULE, [x.clone(), v.clone()])],
            )
        };
        let start = Expr::call(Symbol::ARG, [at(u, &Expr::from(0))]);
        let integral = Expr::call(
            Symbol::INTEGRATE,
            [at(&omega, &t), t, Expr::from(0), x.clone()],
        );
        ev.eval(&Expr::call(Symbol::SIMPLIFY, [start + integral]))?
    } else {
        omega
    };
    ev.step(if phase { "instphase" } else { "instfreq" }, e, &result);
    Ok(Some(result))
}

fn instantaneous(ev: &mut Eval<'_>, e: &Expr, phase: bool) -> Native {
    match e.args() {
        [l] if is_list(l) => of_list(ev, e, l.args(), phase),
        [u] if !u.is_atom() && !is_list(u) => {
            let x = Expr::from(Symbol::new("x"));
            of_expression(ev, e, u, &x, phase)
        }
        [u, x] if !u.is_atom() && !is_list(u) && x.as_symbol().is_some() => {
            of_expression(ev, e, u, x, phase)
        }
        _ => Ok(None),
    }
}

/// `InstantaneousFrequency[z]`, `InstantaneousFrequency[u, x]`.
pub(super) fn instantaneous_frequency(ev: &mut Eval<'_>, e: &Expr) -> Native {
    instantaneous(ev, e, false)
}

/// `InstantaneousPhase[z]`, `InstantaneousPhase[u, x]`.
pub(super) fn instantaneous_phase(ev: &mut Eval<'_>, e: &Expr) -> Native {
    instantaneous(ev, e, true)
}

#[cfg(test)]
mod tests;
