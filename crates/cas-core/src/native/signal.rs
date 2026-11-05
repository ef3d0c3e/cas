//! Signal processing (CAS-521 …), Xcas's commands on lists of samples.
//!
//! Window functions (CAS-521), giac's `*_window`: `HannWindow[u]`,
//! `HannWindow[u, a]`, `HannWindow[u, n1 .. n2]` and `HannWindow[u, a, n1
//! .. n2]` multiply the samples `u_k` by the window `w_k`, `k = 0 … m`,
//! `m = N − 1` for `N` samples (at least 2), or those of indices `n1` to
//! `n2` (from 0), the answer the windowed samples alone. The windows,
//! giac's, value by value:
//!
//! - Bartlett–Hann: `31/50 − 12/25 |k/m − 1/2| − 19/50 cos(2πk/m)`;
//! - Blackman–Harris: `0.35875 − 0.48829 cos(2πk/m) + 0.14128 cos(4πk/m)
//!   − 0.01168 cos(6πk/m)`;
//! - Blackman, `a = 0.16`: `(1 − a)/2 − 1/2 cos(2πk/m) + a/2 cos(4πk/m)`;
//! - Bohman, `x = |2k/m − 1|`: `(1 − x) cos(πx) + sin(πx)/π`;
//! - cosine, `a = 1` (`a > 0`): `sin(πk/m)^a`;
//! - Gaussian, `a = 0.1` (`0 < a ≤ 1/2`): `exp(−1/2 ((2k − m)/(a m))²)`;
//! - Hamming: `27/50 − 23/50 cos(2πk/m)`;
//! - Hann–Poisson, `a = 1`: `sin(πk/m)² exp(−a |m − 2k|/m)`;
//! - Hann: `sin(πk/m)²`;
//! - Parzen, `n = k − m/2`, `h = m/2`: `1 − 6 (n/h)² (1 − |n|/h)` for `|n| ≤
//!   h/2`, else `2 (1 − |n|/h)³`;
//! - Poisson, `a = 1` (`a > 0`): `exp(−a |2k − m|/m)`;
//! - Riemann, `x = 2k/m − 1`: `sin(πx)/(πx)`, 1 at `x = 0`;
//! - triangle, `d = 0` (`d` in `{−1, 0, 1}`): `1 − |2k − m|/(N + d)`;
//! - Tukey, `a = 1/2` (`0 ≤ a ≤ 1`): `(1 + cos(π (2k/(a m) − 1)))/2` for `k <
//!   a m/2`, `(1 + cos(π (2k/(a m) − 2/a + 1)))/2` for `k > m (1 − a/2)`,
//!   else 1 (all 1 for `a = 0`);
//! - Welch: `1 − ((2k − m)/m)²`.
//!
//! RC filters (CAS-522), giac's `lowpass(s, c, [rate])` and `highpass(s, c,
//! [rate])`: `LowPass[s, c]`, `LowPass[s, c, rate]` and the same for
//! `HighPass`, first-order filters with the cutoff frequency `c > 0` at
//! the sample rate (default 44100), `α = 2πc/(rate + 2πc)`: `y_0 = x_0`,
//! then `y_k = y_{k−1} + α (x_k − y_{k−1})` for the low pass and `y_k = (1 −
//! α)(y_{k−1} + x_k − x_{k−1})` for the high pass. As α holds π, exact
//! recurrences would grow without end: the samples are real numbers and
//! the answer balls at the working precision, as giac's are doubles.
//!
//! The windows are exact (`sin(π/4)² = 1/2`, where giac's doubles leave
//! `1.5·10⁻³²` for 0), numbers with an inexact parameter (a decimal), and
//! evaluated with the samples, symbolic ones staying symbolic. A parameter is a real number; out of its range, as a
//! window of fewer than 2 samples, the call stays with a warning.

#![allow(
    clippy::many_single_char_names,
    reason = "the windows' names: samples k of m, parameter a, x, n, h"
)]

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

fn q(n: i64, d: i64) -> Expr {
    Expr::from(Rational::new(Integer::from(n), Integer::from(d)))
}

fn int(n: usize) -> Expr {
    Expr::from(Integer::from(n))
}

fn pi() -> Expr {
    Expr::from(Symbol::PI)
}

fn cos(e: Expr) -> Expr {
    Expr::call(Symbol::COS, [e])
}

fn sin(e: Expr) -> Expr {
    Expr::call(Symbol::SIN, [e])
}

fn exp(e: Expr) -> Expr {
    Expr::call(Symbol::EXP, [e])
}

fn pow(b: Expr, k: Expr) -> Expr {
    Expr::call(Symbol::POWER, [b, k])
}

/// The windows, and the default of the parameter of those that take one.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Window {
    BartlettHann,
    BlackmanHarris,
    Blackman,
    Bohman,
    Cosine,
    Gaussian,
    Hamming,
    HannPoisson,
    Hann,
    Parzen,
    Poisson,
    Riemann,
    Triangle,
    Tukey,
    Welch,
}

impl Window {
    fn of(s: Symbol) -> Option<Window> {
        Some(match s {
            Symbol::BARTLETT_HANN_WINDOW => Window::BartlettHann,
            Symbol::BLACKMAN_HARRIS_WINDOW => Window::BlackmanHarris,
            Symbol::BLACKMAN_WINDOW => Window::Blackman,
            Symbol::BOHMAN_WINDOW => Window::Bohman,
            Symbol::COSINE_WINDOW => Window::Cosine,
            Symbol::GAUSSIAN_WINDOW => Window::Gaussian,
            Symbol::HAMMING_WINDOW => Window::Hamming,
            Symbol::HANN_POISSON_WINDOW => Window::HannPoisson,
            Symbol::HANN_WINDOW => Window::Hann,
            Symbol::PARZEN_WINDOW => Window::Parzen,
            Symbol::POISSON_WINDOW => Window::Poisson,
            Symbol::RIEMANN_WINDOW => Window::Riemann,
            Symbol::TRIANGLE_WINDOW => Window::Triangle,
            Symbol::TUKEY_WINDOW => Window::Tukey,
            Symbol::WELCH_WINDOW => Window::Welch,
            _ => return None,
        })
    }

    /// The parameter's default, for a window that takes one.
    fn default(self) -> Option<Expr> {
        Some(match self {
            Window::Blackman => q(4, 25),
            Window::Cosine | Window::HannPoisson | Window::Poisson => q(1, 1),
            Window::Gaussian => q(1, 10),
            Window::Triangle => q(0, 1),
            Window::Tukey => q(1, 2),
            _ => return None,
        })
    }

    /// Whether `a` (its value `v`) is in the parameter's range; the reason
    /// when it is not.
    fn admits(self, a: &Expr, v: &Rational) -> Result<(), &'static str> {
        let zero = Rational::ZERO;
        let one = Rational::from(Integer::from(1));
        match self {
            Window::Cosine | Window::Poisson if *v <= zero => {
                Err("the window's parameter should be positive")
            }
            Window::Gaussian
                if *v <= zero || *v > Rational::new(Integer::from(1), Integer::from(2)) =>
            {
                Err("the Gaussian window's parameter should be in (0, 1/2]")
            }
            Window::Tukey if *v < zero || *v > one => {
                Err("the Tukey window's parameter should be in [0, 1]")
            }
            Window::Triangle if !matches!(a.to_rational(), Some(d) if d.is_integer() && d.abs() <= one) => {
                Err("the triangle window's parameter should be -1, 0 or 1")
            }
            _ => Ok(()),
        }
    }

    /// `w_k` for `N = m + 1` samples, with the parameter `a` (its value
    /// `v` for the branches).
    #[allow(clippy::too_many_lines, reason = "fifteen windows, one arm each")]
    fn weight(self, k: usize, m: usize, a: &Expr, v: &Rational) -> Expr {
        let (ke, me) = (int(k), int(m));
        let r = Rational::new(Integer::from(k), Integer::from(m));
        let ratio = Expr::from(r.clone());
        let angle = |c: i64| Expr::from(c) * pi() * ratio.clone();
        // 2k − m and its sign-free size, exact.
        let centered = Rational::from(Integer::from(2 * k)) - Rational::from(Integer::from(m));
        let half = Rational::new(Integer::from(1), Integer::from(2));
        match self {
            Window::BartlettHann => {
                q(31, 50) - q(12, 25) * Expr::from((&r - &half).abs()) - q(19, 50) * cos(angle(2))
            }
            Window::BlackmanHarris => {
                q(35_875, 100_000) - q(48_829, 100_000) * cos(angle(2))
                    + q(14_128, 100_000) * cos(angle(4))
                    - q(1_168, 100_000) * cos(angle(6))
            }
            Window::Blackman => {
                (q(1, 1) - a.clone()) / q(2, 1) - q(1, 2) * cos(angle(2))
                    + a.clone() / q(2, 1) * cos(angle(4))
            }
            Window::Bohman => {
                let x = Expr::from((&centered / &Rational::from(Integer::from(m))).abs());
                (q(1, 1) - x.clone()) * cos(pi() * x.clone()) + sin(pi() * x) / pi()
            }
            Window::Cosine => pow(sin(angle(1)), a.clone()),
            Window::Gaussian => {
                let x = Expr::from(centered) / (a.clone() * me);
                exp(q(-1, 2) * pow(x, int(2)))
            }
            Window::Hamming => q(27, 50) - q(23, 50) * cos(angle(2)),
            Window::HannPoisson => {
                let x = Expr::from((&centered / &Rational::from(Integer::from(m))).abs());
                pow(sin(angle(1)), int(2)) * exp(Expr::from(-1) * a.clone() * x)
            }
            Window::Hann => pow(sin(angle(1)), int(2)),
            Window::Parzen => {
                let h = Rational::new(Integer::from(m), Integer::from(2));
                let n = (&Rational::from(Integer::from(k)) - &h).abs();
                let x = &n / &h;
                let one = Rational::from(Integer::from(1));
                let w = if n <= &h * &half {
                    &one - &(&Rational::from(Integer::from(6)) * &(&(&x * &x) * &(&one - &x)))
                } else {
                    let y = &one - &x;
                    &Rational::from(Integer::from(2)) * &(&(&y * &y) * &y)
                };
                Expr::from(w)
            }
            Window::Poisson => {
                let x = Expr::from((&centered / &Rational::from(Integer::from(m))).abs());
                exp(Expr::from(-1) * a.clone() * x)
            }
            Window::Riemann => {
                let x = &centered / &Rational::from(Integer::from(m));
                if x.is_zero() {
                    q(1, 1)
                } else {
                    let x = Expr::from(x);
                    sin(pi() * x.clone()) / (pi() * x)
                }
            }
            Window::Triangle => {
                let d = a.to_rational().unwrap_or(Rational::ZERO);
                let width = Rational::from(Integer::from(m + 1)) + d;
                Expr::from(Rational::from(Integer::from(1)) - &centered.abs() / &width)
            }
            Window::Tukey => {
                if v.is_zero() {
                    return q(1, 1);
                }
                let kq = Rational::from(Integer::from(k));
                let mq = Rational::from(Integer::from(m));
                let rise = &(v * &mq) * &half;
                let fall = &mq * &(&Rational::from(Integer::from(1)) - &(v * &half));
                let x = q(2, 1) * ke / (a.clone() * me);
                let phase = if kq < rise {
                    x - q(1, 1)
                } else if kq > fall {
                    x - q(2, 1) / a.clone() + q(1, 1)
                } else {
                    return q(1, 1);
                };
                q(1, 2) * (q(1, 1) + cos(pi() * phase))
            }
            Window::Welch => {
                let x = &centered / &Rational::from(Integer::from(m));
                Expr::from(Rational::from(Integer::from(1)) - &x * &x)
            }
        }
    }
}

/// A real parameter's value as a rational: exact, or its ball's midpoint.
fn real(ev: &Eval<'_>, a: &Expr) -> Option<Rational> {
    if let Some(r) = a.to_rational() {
        return Some(r);
    }
    match approximation(ev.context(), a, 64).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(b.mid().to_rational()),
        _ => None,
    }
}

/// `n1 .. n2`'s ends as indices.
fn indices(r: &Expr) -> Option<(usize, usize)> {
    let (Some(Symbol::RANGE), [a, b]) = (head(r), r.args()) else {
        return None;
    };
    let index = |e: &Expr| {
        e.to_rational()
            .filter(Rational::is_integer)
            .and_then(|q| q.numer().to_i64())
            .and_then(|n| usize::try_from(n).ok())
    };
    Some((index(a)?, index(b)?))
}

/// `HannWindow[u, …]` and the other windows.
pub(super) fn window(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(kind) = head(e).and_then(Window::of) else {
        return Ok(None);
    };
    let [u, rest @ ..] = e.args() else {
        return Ok(None);
    };
    if head(u) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let (parameter, span) = match rest {
        [] => (None, None),
        [r] if head(r) == Some(Symbol::RANGE) => (None, Some(r)),
        [a] => (Some(a), None),
        [a, r] if head(r) == Some(Symbol::RANGE) => (Some(a), Some(r)),
        _ => return Ok(None),
    };
    let a = match (parameter, kind.default()) {
        (Some(a), Some(_)) => a.clone(),
        (None, Some(d)) => d,
        (None, None) => q(0, 1),
        (Some(_), None) => return Ok(None),
    };
    let Some(v) = real(ev, &a) else {
        return Ok(None);
    };
    if let Err(message) = kind.admits(&a, &v) {
        warn(ev, e, message);
        return Ok(None);
    }
    let samples = u.args();
    let (lo, hi) = match span {
        None => (0, samples.len().saturating_sub(1)),
        Some(r) => match indices(r) {
            Some((lo, hi)) if lo <= hi && hi < samples.len() => (lo, hi),
            _ => {
                warn(ev, e, "the indices should be within the list, in order");
                return Ok(None);
            }
        },
    };
    if samples.len() < 2 || hi == lo {
        warn(ev, e, "a window needs at least 2 samples");
        return Ok(None);
    }
    let m = hi - lo;
    // An inexact parameter, a decimal's ball, makes the weights numbers.
    let inexact = matches!(a.view(), View::Ball(_) | View::ComplexBall(_));
    let bits = ev.context().working_precision();
    let mut out = Vec::with_capacity(m + 1);
    for (k, x) in samples[lo..=hi].iter().enumerate() {
        let mut w = kind.weight(k, m, &a, &v);
        if inexact && let Ok(n) = approximation(ev.context(), &w, bits) {
            w = n;
        }
        out.push(ev.eval(&(w * x.clone()))?);
    }
    let result = list(out);
    ev.step("Window", e, &result);
    Ok(Some(result))
}

/// The default sample rate of giac's filters.
const SAMPLE_RATE: i64 = 44_100;

/// A real number's ball at `bits`.
fn ball(ev: &Eval<'_>, x: &Expr, bits: u64) -> Option<Ball> {
    match approximation(ev.context(), x, bits).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

/// `LowPass[s, c, rate]` and `HighPass[s, c, rate]`.
pub(super) fn filter(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let low = match head(e) {
        Some(Symbol::LOW_PASS) => true,
        Some(Symbol::HIGH_PASS) => false,
        _ => return Ok(None),
    };
    let (s, c, rate) = match e.args() {
        [s, c] => (s, c, Expr::from(SAMPLE_RATE)),
        [s, c, rate] => (s, c, rate.clone()),
        _ => return Ok(None),
    };
    if head(s) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let bits = ev.context().working_precision();
    let prec = bits + 32;
    let (Some(c), Some(rate)) = (ball(ev, c, prec), ball(ev, &rate, prec)) else {
        return Ok(None);
    };
    let positive = |b: &Ball| b.sign() == Some(std::cmp::Ordering::Greater);
    if !positive(&c) {
        warn(ev, e, "the cutoff frequency should be positive");
        return Ok(None);
    }
    if !positive(&rate) {
        warn(ev, e, "the sample rate should be positive");
        return Ok(None);
    }
    let mut xs = Vec::with_capacity(s.args().len());
    for x in s.args() {
        let Some(b) = ball(ev, x, prec) else {
            return Ok(None);
        };
        xs.push(b);
    }
    let Some(pi) = ball(ev, &Expr::from(Symbol::PI), prec) else {
        return Ok(None);
    };
    let two_pi_c = pi.mul_2exp(&Integer::from(1)).mul(&c, prec);
    let alpha = two_pi_c.div(&rate.add(&two_pi_c, prec), prec);
    let beta = Ball::from_rational(&Rational::from(Integer::from(1)), prec).sub(&alpha, prec);
    let mut ys: Vec<Ball> = Vec::with_capacity(xs.len());
    for (k, x) in xs.iter().enumerate() {
        let y = match ys.last() {
            None => x.clone(),
            Some(last) if low => last.add(&alpha.mul(&x.sub(last, prec), prec), prec),
            Some(last) => beta.mul(&last.add(&x.sub(&xs[k - 1], prec), prec), prec),
        };
        ys.push(y);
    }
    let result = list(
        ys.iter()
            .map(|y| Expr::from(y.round_to_accuracy(bits)))
            .collect(),
    );
    ev.step(if low { "LowPass" } else { "HighPass" }, e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
