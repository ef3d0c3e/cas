//! The digamma and polygamma functions of real balls (CAS-DOC-31 §2.2,
//! CAS-351): `ψ = Γ′/Γ` and its derivatives `ψ^(m)`.
//!
//! - **Shift**: `ψ^(m)(x) = ψ^(m)(x + N) − (−1)^m m! Σ_{k<N} (x +
//!   k)^(−(m+1))`, the sum on the ball itself, `N` the fewest steps that
//!   bring `z = x + N` where the expansion below reaches the precision with
//!   at most [`MAX_TERMS`] terms.
//! - **Asymptotic expansion** at `z`: `ψ(z) = log z − 1/(2z) − Σ_{k<K}
//!   B_{2k}/(2k z^(2k)) + R` and, for `m ≥ 1`, `ψ^(m)(z) = (−1)^(m+1) (m −
//!   1)!/z^m (1 + m/(2z) + Σ_{k<K} B_{2k} C(2k + m − 1, m − 1)/z^(2k)) +
//!   R`. For a real `z > 0` the remainder is at most the first term left
//!   out: `ψ^(m)` is the Laplace transform of `t^(m−1) t/(1 − e^(−t))`
//!   (DLMF 5.15.1, 5.9.13 for `m = 0`), and `t/(1 − e^(−t)) − 1 − t/2 = Σ_j
//!   2t²/(t² + 4π² j²)`, whose terms expanded in `t²` leave remainders of
//!   the sign of the first term left out and at most its size. The ball
//!   `z` adds its radius times `sup |ψ^(m+1)| = |ψ^(m+1)(lo)| ≤ (m + 1)!
//!   (1/lo^(m+2) + 1/((m + 1) lo^(m+1)))`.
//! - **Reflection** below −[`MAX_SHIFT`]: `ψ^(m)(x) = (−1)^m ψ^(m)(1 − x)
//!   − π^(m+1) Q_m(cot πx)` (DLMF 5.15.6), with `Q_0(c) = c` and `Q_(m+1)(c)
//!   = −(1 + c²) Q_m′(c)`, so that `d^m/dx^m cot(πx) = π^m Q_m(cot πx)`.
//! - A ball meeting 0 or a negative integer is indeterminate, by
//!   `Ball::gamma`'s test, and so is an order above [`MAX_ORDER`].
//!
//! An exact input is computed again at more bits while its relative
//! accuracy falls short: `ψ` has a zero at 1.4616… and one between each
//! two poles.

use std::cmp::Ordering;

use super::gamma::{bernoulli, meets_pole};
use super::{exponent, one, refined, upper};
use crate::ball::{Ball, Float, Mag, Round};
use crate::{Integer, Natural, Rational};

/// The most terms of the asymptotic expansion: Bernoulli numbers to
/// `B_1200`, as `Ball::gamma`'s.
const MAX_TERMS: u64 = 600;

/// The highest order `m`.
const MAX_ORDER: u64 = 256;

/// How far below 0 the shift alone serves; the reflection past it.
const MAX_SHIFT: u64 = 64;

/// `log2 n!`, an estimate.
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn log2_factorial(n: u64) -> f64 {
    if n < 32 {
        return (2..=n).map(|k| (k as f64).log2()).sum();
    }
    let m = n as f64;
    (m * m.ln() - m + 0.5 * (2.0 * std::f64::consts::PI * m).ln()) / std::f64::consts::LN_2
}

/// The number `K` of terms at `z` (`log2 z ≈ log2_z`) whose first left out
/// is below `2^-target` of the leading term; `None` past [`MAX_TERMS`].
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn terms(log2_z: f64, m: u64, target: u64) -> Option<u64> {
    let target = target as f64;
    let log2_two_pi_z = (2.0 * std::f64::consts::PI).log2() + log2_z;
    (1..=MAX_TERMS).find(|&k| {
        // |B_{2k}| ≈ 2 (2k)!/(2π)^(2k); the term over (m − 1)!/z^m, or over
        // 1 for m = 0.
        let size = if m == 0 {
            1.0 + log2_factorial(2 * k - 1)
        } else {
            1.0 + log2_factorial(2 * k + m - 1) - log2_factorial(m - 1)
        };
        size - 2.0 * k as f64 * log2_two_pi_z < -target
    })
}

/// `log2 x` for a positive float, saturated.
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn log2_of(x: &Float) -> f64 {
    let e = exponent(x);
    if e.unsigned_abs() > 1000 {
        return e as f64;
    }
    x.to_f64().log2()
}

/// The steps `N` up from the midpoint `x` and the terms `K` at `x + N`;
/// `None` past `2^32` steps.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    reason = "estimates of sizes"
)]
fn plan(x: &Float, m: u64, wp: u64) -> Option<(u64, u64)> {
    let target = wp + 12;
    if x.sign() == Ordering::Greater
        && let Some(k) = terms(log2_of(x), m, target)
    {
        return Some((0, k));
    }
    // About wp/9 + m/2 from 0, more while the terms are not enough.
    let start = if exponent(x) > 60 { 0.0 } else { x.to_f64() };
    let mut z = (target as f64 / 9.0 + m as f64 / 2.0 + 2.0).max(start + 1.0);
    loop {
        if let Some(k) = terms(z.log2(), m, target) {
            let steps = (z - start).ceil().max(0.0);
            return (steps < 4_294_967_296.0).then_some((steps as u64, k));
        }
        z *= 1.5;
    }
}

/// `C(2k + m − 1, m − 1)`, the factor of `B_{2k}` in the expansion, from
/// the previous: `(2k + m − 1)(2k + m − 2)/((2k)(2k − 1))` times it.
fn next_binomial(previous: &Integer, k: u64, m: u64) -> Integer {
    let num = &(previous * &Integer::from(2 * k + m - 1)) * &Integer::from(2 * k + m - 2);
    let (q, r) = num.div_rem_euclid(&Integer::from(2 * k * (2 * k - 1)));
    debug_assert!(r.is_zero(), "a binomial coefficient");
    q
}

/// `ψ^(m)(z)` for an exact `z > 0` by the expansion with `k_terms` terms,
/// at about `wp` bits.
fn asymptotic(z: &Float, m: u64, wp: u64, k_terms: u64) -> Ball {
    let q = wp + 32;
    let zb = Ball::exact(z.clone());
    let w = one().div(&zb, q);
    let w2 = w.mul(&w, q);
    // Σ_{k<K} c_k w^(2k), and the first term left out, |c_K| w^(2K).
    let mut sum = Ball::ZERO;
    let mut power = w2.clone();
    let mut binomial = Integer::ONE;
    let mut rest = Mag::ZERO;
    for k in 1..=k_terms {
        let b = bernoulli(2 * k).expect("within the table");
        let c = if m == 0 {
            &b / &Rational::from(Integer::from(2 * k))
        } else {
            binomial = next_binomial(&binomial, k, m);
            &b * &Rational::from(binomial.clone())
        };
        if k == k_terms {
            rest = Mag::from_integer(&c.numer().abs())
                .div(&Mag::from_integer_lower(c.denom()))
                .mul(&upper(&power));
            break;
        }
        sum = sum.add(&Ball::from_rational(&c, q).mul(&power, q), q);
        power = power.mul(&w2, q);
    }
    let sum = sum.add_error(&rest);
    if m == 0 {
        // log z − w/2 − Σ.
        return zb.log(q).sub(&w.mul_2exp(&-Integer::ONE), q).sub(&sum, q);
    }
    // (−1)^(m+1) (m − 1)! w^m (1 + m w/2 + Σ).
    let half_m = Ball::from_integer(&Integer::from(m))
        .mul(&w, q)
        .mul_2exp(&-Integer::ONE);
    let bracket = one().add(&half_m, q).add(&sum, q);
    let factorial = Ball::from_integer(&Integer::from(Natural::factorial(m - 1)));
    let i = i64::try_from(m).expect("a small order");
    let value = factorial.mul(&w.pow(i, q), q).mul(&bracket, q);
    if m % 2 == 0 { value.neg() } else { value }
}

/// `(m + 1)! (1/lo^(m+2) + 1/((m + 1) lo^(m+1)))`, a bound of
/// `|ψ^(m+1)|` from `lo > 0` on.
fn slope(lo: &Mag, m: u64) -> Mag {
    let factorial = Mag::from_integer(&Integer::from(Natural::factorial(m + 1)));
    let inverse = Mag::ONE.div(lo);
    let first = inverse.pow(m + 2);
    let second = inverse
        .pow(m + 1)
        .div(&Mag::from_integer_lower(&Integer::from(m + 1)));
    factorial.mul(&first.add(&second))
}

/// `ψ^(m)` of a ball whose points are positive and large enough for `K`
/// terms at its midpoint: the midpoint's value and the radius's spread.
fn at_large(z: &Ball, m: u64, wp: u64, k_terms: u64) -> Ball {
    let core = asymptotic(z.mid(), m, wp, k_terms);
    if z.is_exact() {
        return core;
    }
    let lo = Mag::from_float_lower(z.mid()).sub_lower(z.rad());
    if lo.is_zero() {
        return Ball::INDETERMINATE;
    }
    core.add_error(&z.rad().mul(&slope(&lo, m)))
}

/// `ψ^(m)(x)` for a ball meeting no pole, at or above −[`MAX_SHIFT`], by
/// the shift.
fn shifted(x: &Ball, m: u64, wp: u64) -> Ball {
    let Some((steps, k_terms)) = plan(x.mid(), m, wp) else {
        return Ball::INDETERMINATE;
    };
    let q = wp + 2 * u64::from(u64::BITS - steps.leading_zeros()) + 16;
    let power = i64::try_from(m + 1).expect("a small order");
    let mut sum = Ball::ZERO;
    for k in 0..steps {
        let term = x
            .add(&Ball::from_integer(&Integer::from(k)), q)
            .pow(power, q);
        sum = sum.add(&one().div(&term, q), q);
    }
    let z = x.add(&Ball::from_integer(&Integer::from(steps)), q);
    let tail = at_large(&z, m, q, k_terms);
    // − (−1)^m m! Σ.
    let factorial = Ball::from_integer(&Integer::from(Natural::factorial(m)));
    let sum = factorial.mul(&sum, q);
    if m % 2 == 0 {
        tail.sub(&sum, q)
    } else {
        tail.add(&sum, q)
    }
}

/// The coefficients of `Q_m`, lowest first: `Q_0(c) = c`, `Q_(m+1)(c) =
/// −(1 + c²) Q_m′(c)`.
fn cot_polynomial(m: u64) -> Vec<Integer> {
    let mut p = vec![Integer::ZERO, Integer::ONE];
    for _ in 0..m {
        // The derivative, then −(1 + c²) times it.
        let derivative: Vec<Integer> = p
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, a)| a * &Integer::from(i))
            .collect();
        let mut next = vec![Integer::ZERO; derivative.len() + 2];
        for (i, a) in derivative.iter().enumerate() {
            next[i] = &next[i] - a;
            next[i + 2] = &next[i + 2] - a;
        }
        while next.last().is_some_and(Integer::is_zero) {
            next.pop();
        }
        p = next;
    }
    p
}

/// `ψ^(m)(x)` for a ball below −[`MAX_SHIFT`] meeting no pole, by the
/// reflection.
fn reflected(x: &Ball, m: u64, wp: u64) -> Ball {
    // 1 − x to the bits of x's integer part and more, as ψ(1 − x) ≈ log|x|.
    let above = u64::try_from(exponent(x.mid()).max(0)).unwrap_or(0);
    let q = wp + above + 16;
    let complement = one().sub(x, q + above);
    let mirrored = shifted(&complement, m, q);
    let cot = x.cos_pi(q).div(&x.sin_pi(q), q);
    let coefficients = cot_polynomial(m);
    let mut q_m = Ball::ZERO;
    for a in coefficients.iter().rev() {
        q_m = q_m.mul(&cot, q).add(&Ball::from_integer(a), q);
    }
    let i = i64::try_from(m + 1).expect("a small order");
    let scale = Ball::const_pi(q).pow(i, q);
    let mirrored = if m % 2 == 0 { mirrored } else { mirrored.neg() };
    mirrored.sub(&scale.mul(&q_m, q), q)
}

impl Ball {
    /// `ψ(x)`, the digamma function, for every `x` in the ball, at `prec`
    /// bits; indeterminate where the ball meets a pole.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn digamma(&self, prec: u64) -> Ball {
        self.polygamma(0, prec)
    }

    /// `ψ^(m)(x)`, the polygamma function of order `m`, for every `x` in
    /// the ball, at `prec` bits; indeterminate where the ball meets a pole,
    /// and for `m` above [`MAX_ORDER`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn polygamma(&self, m: u64, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::polygamma: precision 0");
        if !self.is_finite() || m > MAX_ORDER || meets_pole(self) {
            return Ball::INDETERMINATE;
        }
        let Some(r) = self.rad().to_float() else {
            return Ball::INDETERMINATE;
        };
        let lo = self.mid().sub_round(&r, 64, Round::Floor).0;
        let below = Float::from(-i64::try_from(MAX_SHIFT).expect("small"));
        let value = refined(self, prec, |wp| {
            if lo < below {
                reflected(self, m, wp)
            } else {
                shifted(self, m, wp)
            }
        });
        value.round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
