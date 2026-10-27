//! The error functions of real balls (CAS-DOC-15 §8, CAS-246): `erf` and
//! `erfc = 1 − erf`.
//!
//! For an exact `x ≥ 0`:
//!
//! - **The series** `erf x = 2x e^(−x²)/√π Σ_{k≥0} (2x²)^k/(1·3⋯(2k + 1))`
//!   (DLMF 7.6.2), of positive terms: summed until a term is below
//!   `2^-(wp + 12)` of the sum once they halve (`2x² ≤ (2k + 1)/2`), the
//!   rest at most that term.
//! - **The asymptotic expansion** `erfc x = e^(−x²)/(x√π) Σ_{n<N} (−1)^n
//!   (2n − 1)!!/(2x²)^n`, whose error for a real `x > 0` is at most the
//!   first term left out (DLMF 7.12(ii)): for `x² ≥ (wp + 16) log 2 + 4`,
//!   where its smallest term, about `√2 e^(−x²)`, is below the precision.
//!
//! `erf` is the series below that bound and `1 − erfc` above; `erfc` is the
//! expansion above it and `1 − erf` below, `erf` taken at `x² log2 e` more
//! bits for the cancellation. Negative arguments: `erf(−x) = −erf x`,
//! `erfc(−x) = 1 + erf x`. A ball's radius `ρ` adds `ρ (2/√π) e^(−d²)`, `d`
//! the ball's distance to 0, since `|erf′| = (2/√π) e^(−x²)`.

use super::{exponent, one, upper, working};
use crate::Integer;
use crate::ball::{Ball, Float, Mag};

/// `x²` at least where the expansion reaches `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "a threshold")]
fn large(x2: &Float, wp: u64) -> bool {
    if exponent(x2) > 62 {
        return true;
    }
    x2.to_f64() >= (wp as f64 + 16.0) * std::f64::consts::LN_2 + 4.0
}

/// `2/√π` at `q` bits.
fn two_over_root_pi(q: u64) -> Ball {
    Ball::from_integer(&Integer::from(2)).div(&Ball::const_pi(q).sqrt(q), q)
}

/// `erf x` for an exact `x ≥ 0` by the series, at about `wp` bits.
fn series(x: &Float, wp: u64) -> Ball {
    let x2 = x.mul(x);
    let q = wp + 32;
    let ratio = Ball::exact(x2.mul_2exp(&Integer::ONE));
    let small = -Integer::from(wp + 12);
    let mut term = one();
    let mut sum = one();
    let mut k: u64 = 0;
    loop {
        k += 1;
        let odd = Float::from(Integer::from(2 * k + 1));
        term = term.mul(&ratio, q).div(&Ball::exact(odd.clone()), q);
        sum = sum.add(&term, q);
        // Halving once 2·2x² ≤ 2k + 1.
        let halving = x2.mul_2exp(&Integer::from(2)) <= odd;
        if halving
            && upper(&term)
                <= Mag::from_float_lower(sum.mid()).mul_lower(&Mag::power_of_two(&small))
        {
            break;
        }
    }
    let rest = upper(&term);
    let scale = two_over_root_pi(q)
        .mul(&Ball::exact(x.clone()), q)
        .mul(&Ball::exact(x2.neg()).exp(q), q);
    scale.mul(&sum.add_error(&rest), q)
}

/// `erfc x` for an exact `x > 0`, `x²` past [`large`], by the asymptotic
/// expansion at about `wp` bits.
fn asymptotic(x: &Float, wp: u64) -> Ball {
    let x2 = x.mul(x);
    let q = wp + 32;
    let small = Mag::power_of_two(&-Integer::from(wp + 12));
    // −(2n − 1)/(2x²) from one term to the next.
    let twice = Ball::exact(x2.mul_2exp(&Integer::ONE));
    let mut term = one();
    let mut sum = one();
    let mut n: u64 = 0;
    let bound = loop {
        n += 1;
        let odd = Ball::exact(Float::from(Integer::from(2 * n - 1)));
        term = term.mul(&odd, q).div(&twice, q).neg();
        // The terms decrease while 2n − 1 < 2x²; stop at the precision, or
        // before they grow.
        let growing = x2.mul_2exp(&Integer::ONE) < Float::from(Integer::from(2 * n + 1));
        if upper(&term) <= small || growing {
            break upper(&term);
        }
        sum = sum.add(&term, q);
    };
    let scale = Ball::exact(x2.neg()).exp(q).div(
        &Ball::exact(x.clone()).mul(&Ball::const_pi(q).sqrt(q), q),
        q,
    );
    scale.mul(&sum.add_error(&bound), q)
}

/// `erf x` for an exact `x ≥ 0`.
fn erf_exact(x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    if large(&x.mul(x), wp) {
        return one().sub(&asymptotic(x, wp), wp);
    }
    series(x, wp)
}

/// `erfc x` for an exact `x ≥ 0`.
fn erfc_exact(x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return one();
    }
    let x2 = x.mul(x);
    if large(&x2, wp) {
        return asymptotic(x, wp);
    }
    // 1 − erf x cancels about x² log2 e bits.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a small x²"
    )]
    let extra = (x2.to_f64() * std::f64::consts::LOG2_E).ceil() as u64 + 8;
    one().sub(&series(x, wp + extra), wp + extra)
}

/// `ρ (2/√π) e^(−d²)` for the ball, `d` its distance to 0.
fn spread(b: &Ball) -> Mag {
    if b.is_exact() {
        return Mag::ZERO;
    }
    let distance = Mag::from_float_lower(&b.mid().abs()).sub_lower(b.rad());
    let d = distance.to_float().unwrap_or(Float::ZERO);
    let gauss = Ball::exact(d.mul(&d).neg()).exp(64);
    b.rad()
        .mul(&upper(&two_over_root_pi(64)))
        .mul(&upper(&gauss))
}

impl Ball {
    /// `erf x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn erf(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::erf: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let wp = working(prec);
        let m = self.mid();
        let value = erf_exact(&m.abs(), wp);
        let value = if m.is_negative() { value.neg() } else { value };
        value.add_error(&spread(self)).round_to_accuracy(prec)
    }

    /// `erfc x = 1 − erf x` for every `x` in the ball, at `prec` bits,
    /// relative to its value however small.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn erfc(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::erfc: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let wp = working(prec);
        let m = self.mid();
        let value = if m.is_negative() {
            one().add(&erf_exact(&m.abs(), wp), wp)
        } else {
            erfc_exact(m, wp)
        };
        value.add_error(&spread(self)).round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
