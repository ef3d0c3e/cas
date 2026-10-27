//! The gamma function of real balls (CAS-DOC-15 §6, CAS-242).
//!
//! For an exact `t > 0`, `log Γ(t)` or `Γ(t)` by the cheaper of two series:
//!
//! - **Stirling's**, `log Γ(t) = (t − 1/2) log t − t + log(2π)/2 + Σ_{k<K}
//!   c_k / t^(2k−1)`, `c_k = B_{2k} / (2k (2k − 1))`, whose remainder for a
//!   real `t > 0` is at most the next term in absolute value (DLMF 5.11.ii);
//!   `K` the fewest terms that reach the precision, at most [`MAX_TERMS`]. The
//!   coefficients come from the tangent numbers `T_k` (Brent and Harvey's
//!   integer recurrence), `c_k = (−1)^(k−1) T_k / ((2k − 1) 4^k (4^k − 1))`,
//!   cached.
//! - **Brent's**, `Γ(t) = N^t e^(−N) Σ_{k≥0} N^k / (t (t + 1) ⋯ (t + k)) +
//!   Γ(t, N)`: the series of positive terms summed until `t + k ≥ 2N`, past
//!   which the terms halve, so the rest is at most the last term; `Γ(t, N) ≤
//!   2 N^(t−1) e^(−N)` for `N ≥ 2t` (`N^(t−1) e^(−N)/(1 − (t − 1)/N)`), and
//!   `N` large enough for that to be below the precision against `Γ(t) ≥
//!   0.885`. About `2N` steps of ball arithmetic.
//!
//! A ball's radius `ρ` adds `ρ sup |ψ|` to `log Γ` over it, with `|ψ(x)| ≤
//! |log x| + 1/x` for `x > 0` (DLMF 5.9.ii's `log x − 1/x < ψ(x) < log x`),
//! whose largest value is at an end. Then:
//!
//! - `x ≥ 1/2`: `Γ(x) = exp(log Γ(x))`, or Brent's `Γ` itself for an exact
//!   `x`;
//! - `0 < x < 1/2`: `Γ(x) = Γ(x + 1) / x`;
//! - `x < 0` and not a pole: `Γ(x) = Γ(x + n) / (x (x + 1) ⋯ (x + n − 1))`
//!   into `[1, 2)` for at most [`MAX_SHIFT`] steps, past that the reflection
//!   `Γ(x) = π e^(−log Γ(1 − x)) / sin(πx)`, finite however small;
//! - a ball meeting a pole, 0 or a negative integer, is indeterminate, and
//!   so is `Γ` past `exp`'s range (an exponent above `HUGE_EXP`).
//!
//! \[BH11\] Richard P. Brent and David Harvey, "Fast computation of Bernoulli,
//! tangent and secant numbers" (2011).

use std::cmp::Ordering;
use std::sync::{Arc, PoisonError, RwLock};

use super::exp::HUGE_EXP;
use super::{exponent, one, upper, working};
use crate::ball::{Ball, Float, Mag, Round};
use crate::{Integer, Rational};

/// The most terms of Stirling's series: Bernoulli numbers to `B_1200`.
const MAX_TERMS: usize = 600;

/// The most steps of the recurrence from a negative argument up to `[1, 2)`.
const MAX_SHIFT: u64 = 64;

/// The most steps of Brent's series.
const MAX_BRENT_STEPS: u64 = 1 << 24;

/// Stirling's coefficients `c_1, c_2, …`, cached.
static COEFFICIENTS: RwLock<Option<Arc<Vec<Rational>>>> = RwLock::new(None);

/// The tangent numbers `T_1 … T_n` (index `k − 1`), by Brent and Harvey's
/// recurrence (\[BH11\], `TangentNumbers`).
fn tangent_numbers(n: usize) -> Vec<Integer> {
    let mut t = vec![Integer::ZERO; n + 1];
    if n == 0 {
        return Vec::new();
    }
    t[1] = Integer::ONE;
    for k in 2..=n {
        t[k] = &t[k - 1] * &Integer::from(k - 1);
    }
    for k in 2..=n {
        for j in k..=n {
            t[j] = &(&t[j - 1] * &Integer::from(j - k)) + &(&t[j] * &Integer::from(j - k + 2));
        }
    }
    t.remove(0);
    t
}

/// The first `n` coefficients `c_k = B_{2k} / (2k (2k − 1))`, from the cache
/// or computed (at least twice what the cache had).
fn coefficients(n: usize) -> Arc<Vec<Rational>> {
    let have = {
        let c = COEFFICIENTS.read().unwrap_or_else(PoisonError::into_inner);
        if let Some(v) = c.as_ref()
            && v.len() >= n
        {
            return Arc::clone(v);
        }
        c.as_ref().map_or(0, |v| v.len())
    };
    let target = n.max(2 * have).clamp(1, MAX_TERMS);
    let tangent = tangent_numbers(target);
    let values: Vec<Rational> = tangent
        .iter()
        .enumerate()
        .map(|(i, t)| {
            let k = u64::try_from(i + 1).expect("few terms");
            let four = Integer::from(4).pow(k);
            let den = &(&Integer::from(2 * k - 1) * &four) * &(&four - &Integer::ONE);
            let c = Rational::new(t.clone(), den);
            if k % 2 == 0 { -c } else { c }
        })
        .collect();
    let values = Arc::new(values);
    let mut c = COEFFICIENTS.write().unwrap_or_else(PoisonError::into_inner);
    if c.as_ref().is_none_or(|v| v.len() < values.len()) {
        *c = Some(Arc::clone(&values));
    }
    values
}

/// An estimate of `log2 |c_k|`, from `|B_{2k}| ≈ 2 (2k)! / (2π)^(2k)`.
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn log2_coefficient(k: usize) -> f64 {
    let m = 2.0 * k as f64;
    // log2 m! by Stirling's formula, exact enough for a choice.
    let log2_factorial =
        (m * m.ln() - m + 0.5 * (2.0 * std::f64::consts::PI * m).ln()) / std::f64::consts::LN_2;
    1.0 + log2_factorial - m * (2.0 * std::f64::consts::PI).log2() - (m * (m - 1.0)).log2()
}

/// The fewest terms of Stirling's series at `t` (`log2 t ≈ log2_t`) whose
/// remainder is below `2^-target`; `None` past [`MAX_TERMS`].
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn stirling_terms(log2_t: f64, target: u64) -> Option<usize> {
    let target = target as f64 + 8.0;
    (1..=MAX_TERMS).find(|&k| log2_coefficient(k) - (2 * k - 1) as f64 * log2_t < -target)
}

/// `log Γ(t)` for an exact `t > 0` by Stirling's series with `k` terms, at
/// about `q` bits absolute.
fn stirling(t: &Float, q: u64, k: usize) -> Ball {
    let c = coefficients(k);
    let y = Ball::exact(t.clone());
    let half = Ball::exact(Float::from(1i64).mul_2exp(&-Integer::ONE));
    let two_pi = Ball::const_pi(q).mul_2exp(&Integer::ONE);
    let mut sum = y
        .sub(&half, q)
        .mul(&y.log(q), q)
        .sub(&y, q)
        .add(&two_pi.log(q).mul_2exp(&-Integer::ONE), q);
    let inverse = one().div(&y, q);
    let square = inverse.mul(&inverse, q);
    let mut power = inverse;
    for coefficient in &c[..k - 1] {
        sum = sum.add(&Ball::from_rational(coefficient, q).mul(&power, q), q);
        power = power.mul(&square, q);
    }
    // The remainder, at most |c_k| / t^(2k−1): `power` is 1/t^(2k−1).
    let last = &c[k - 1];
    let bound = Mag::from_integer(&last.numer().abs()).div(&Mag::from_integer_lower(last.denom()));
    sum.add_error(&bound.mul(&upper(&power)))
}

/// Brent's `N` for an exact `t > 0` at `wp` bits, and the steps of its
/// series; `None` past [`MAX_BRENT_STEPS`].
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "estimates of sizes"
)]
fn brent_size(t: f64, wp: u64) -> Option<(u64, u64)> {
    // N − (t − 1) log N ≥ (wp + 10) log 2, and N ≥ 2t + 2.
    let target = (wp as f64 + 10.0) * std::f64::consts::LN_2;
    let mut n = target;
    for _ in 0..8 {
        n = target + (t - 1.0).max(0.0) * n.max(2.0).ln();
    }
    let n = n.max(2.0 * t + 2.0).ceil();
    // The terms past t + k ≥ 2N halve.
    let steps = (2.0 * n - t).ceil().max(0.0) + 2.0;
    (steps.is_finite() && steps < MAX_BRENT_STEPS as f64).then_some((n as u64, steps as u64))
}

/// `Γ(t)` for an exact `t > 0` by Brent's series with `N = n`, `steps`
/// terms, at about `wp` bits.
fn brent(t: &Float, wp: u64, n: u64, steps: u64) -> Ball {
    let q = wp + 2 * u64::from(u64::BITS - steps.leading_zeros()) + 16;
    let big_n = Integer::from(n);
    let bn = Ball::from_integer(&big_n);
    let mut term = one().div(&Ball::exact(t.clone()), q);
    let mut sum = term.clone();
    for k in 1..steps {
        let k = i64::try_from(k).expect("fewer than 2^63 steps");
        let denominator = Ball::exact(t.add(&Float::from(k)));
        term = term.mul(&bn, q).div(&denominator, q);
        sum = sum.add(&term, q);
    }
    // N^t e^(−N).
    let exponent_ball = Ball::exact(t.clone()).mul(&bn.log(q), q).sub(&bn, q);
    let prefactor = exponent_ball.exp(q);
    // The series' rest, at most the last term, and Γ(t, N) ≤ 2 N^(t−1) e^(−N).
    let rest = upper(&term)
        .add(&Mag::from_integer(&Integer::from(2)).div(&Mag::from_integer_lower(&big_n)));
    prefactor
        .mul(&sum, q)
        .add_error(&upper(&prefactor).mul(&rest))
}

/// `log2 t`, saturated, for a positive float.
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn log2_of(t: &Float) -> f64 {
    let e = exponent(t);
    if e.abs() > 1000 {
        return e as f64;
    }
    t.to_f64().abs().log2()
}

/// `Γ(t)` or `log Γ(t)` for an exact `t > 0`: Stirling's `log Γ`, or
/// Brent's `Γ` where that is cheaper.
enum Value {
    Gamma(Ball),
    Log(Ball),
}

/// `t`'s value at `wp` bits, the cheaper series.
fn exact_value(t: &Float, wp: u64) -> Value {
    let log2_t = log2_of(t);
    // log Γ(t) ≈ t log t: its bits above the point, for its absolute
    // accuracy.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "an estimate"
    )]
    let above = if log2_t > 1.0 {
        (log2_t + log2_t.max(1.0).log2() + 2.0) as u64
    } else {
        2
    };
    let q = wp + above + 8;
    let stirling_k = stirling_terms(log2_t, q);
    let brent_steps = if log2_t < 64.0 {
        brent_size(t.to_f64(), wp)
    } else {
        None
    };
    match (stirling_k, brent_steps) {
        (Some(k), Some((_, steps))) if (k as u64) * 4 < steps => Value::Log(stirling(t, q, k)),
        (_, Some((n, steps))) => Value::Gamma(brent(t, wp, n, steps)),
        (Some(k), None) => Value::Log(stirling(t, q, k)),
        (None, None) => Value::Log(Ball::INDETERMINATE),
    }
}

/// `sup |ψ|` over `[a, b]`, `0 < a`: `|log x| + 1/x` at an end, from above.
fn digamma_bound(a: &Float, b: &Float) -> Mag {
    let at = |x: &Float| {
        let x = Ball::exact(x.clone());
        upper(&x.log(64)).add(&upper(&one().div(&x, 64)))
    };
    let (u, v) = (at(a), at(b));
    if u > v { u } else { v }
}

/// `log Γ(y)` for a ball `y` whose points are positive: the midpoint's,
/// then `ρ sup |ψ|`.
pub(super) fn log_gamma_positive(y: &Ball, wp: u64) -> Ball {
    let m = y.mid();
    let core = match exact_value(m, wp) {
        Value::Log(l) => l,
        Value::Gamma(g) => g.log(wp + 8),
    };
    if y.is_exact() {
        return core;
    }
    let Some(r) = y.rad().to_float() else {
        return Ball::INDETERMINATE;
    };
    let a = m.sub_round(&r, 64, Round::Floor).0;
    let b = m.add_round(&r, 64, Round::Ceil).0;
    core.add_error(&y.rad().mul(&digamma_bound(&a, &b)))
}

/// `Γ(y)` for a ball `y` whose points are positive.
fn gamma_positive(y: &Ball, wp: u64) -> Ball {
    if y.is_exact() {
        return match exact_value(y.mid(), wp) {
            Value::Gamma(g) => g,
            Value::Log(l) => l.exp(wp),
        };
    }
    log_gamma_positive(y, wp).exp(wp)
}

/// Whether the ball meets 0 or a negative integer.
fn meets_pole(x: &Ball) -> bool {
    let m = x.mid();
    let Some(r) = x.rad().to_float() else {
        return true;
    };
    let lo = m.sub_round(&r, 64, Round::Floor).0;
    if lo.sign() == Ordering::Greater {
        return false;
    }
    // The midpoint itself an integer, at most 0?
    let (mant, e) = m.to_integer_exp();
    let integral = !e.is_negative()
        || mant.is_zero()
        || mant
            .unsigned_abs()
            .trailing_zeros()
            .is_some_and(|z| Integer::from(z) >= -e.clone());
    if integral && m.sign() != Ordering::Greater {
        return true;
    }
    // A radius of at least 1/2 reaches an integer from lo ≤ 0.
    if *x.rad() >= Mag::power_of_two(&-Integer::ONE) {
        return true;
    }
    // Otherwise the midpoint has a fractional part, a short mantissa: the
    // integers between the ends, exactly.
    let hi = m.add_round(&r, 64, Round::Ceil).0;
    let (lo, hi) = (lo.to_rational(), hi.to_rational());
    let first = lo.ceil();
    let last = hi.floor().min(Integer::ZERO);
    first <= last
}

impl Ball {
    /// `Γ(x)` for every `x` in the ball, at `prec` bits; indeterminate where
    /// the ball meets a pole, 0 or a negative integer.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn gamma(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::gamma: precision 0");
        if !self.is_finite() || meets_pole(self) {
            return Ball::INDETERMINATE;
        }
        let wp = working(prec);
        let m = self.mid();
        let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
        let value = if *m >= half {
            // Past exp's range, Γ overflows it.
            if exponent(m) > HUGE_EXP {
                return Ball::INDETERMINATE;
            }
            gamma_positive(self, wp)
        } else if m.sign() == Ordering::Greater {
            // Γ(x) = Γ(x + 1)/x, x + 1 in (1, 3/2).
            gamma_positive(&self.add(&one(), wp), wp).div(self, wp)
        } else {
            negative(self, wp)
        };
        value.round_to_accuracy(prec)
    }
}

/// `Γ(x)` for a ball below 0 that meets no pole.
fn negative(x: &Ball, wp: u64) -> Ball {
    let m = x.mid();
    // Steps up to [1, 2): n = ⌈1 − m⌉ when that is small.
    let steps = if exponent(m) <= 7 {
        let need = (&Rational::ONE - &m.to_rational()).ceil();
        need.to_i64().and_then(|n| u64::try_from(n).ok())
    } else {
        None
    };
    if let Some(n) = steps.filter(|&n| n <= MAX_SHIFT) {
        let q = wp + 16;
        let mut product = one();
        for k in 0..n {
            product = product.mul(&x.add(&Ball::from_integer(&Integer::from(k)), q), q);
        }
        let shifted = x.add(&Ball::from_integer(&Integer::from(n)), q);
        return gamma_positive(&shifted, q).div(&product, q);
    }
    // π e^(−log Γ(1 − x)) / sin(πx); past exp's range e^(−L) is a bound
    // about 0, for which a low precision does.
    let q = if exponent(m) > HUGE_EXP { 64 } else { wp + 16 };
    let y = one().sub(x, q);
    let l = log_gamma_positive(&y, q);
    Ball::const_pi(q)
        .mul(&l.neg().exp(q), q)
        .div(&x.sin_pi(q), q)
}

#[cfg(test)]
mod tests;
