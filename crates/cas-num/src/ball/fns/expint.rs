//! The exponential, sine and cosine integrals of real balls (CAS-DOC-31
//! §2.1, CAS-350): `Ei(x) = PV ∫_{−∞}^x e^t/t dt`, `Si(x) = ∫_0^x sin t/t
//! dt` and `Ci(x) = γ + log x + ∫_0^x (cos t − 1)/t dt`.
//!
//! For an exact `x` with `|x|` below a threshold `T`, about `(wp + 16) log
//! 2`, the series:
//!
//! - `Ei(x) = γ + log |x| + Σ_{k≥1} x^k/(k k!)`;
//! - `Si(x) = Σ_{k≥0} (−1)^k x^(2k+1)/((2k + 1) (2k + 1)!)`;
//! - `Ci(x) = γ + log x + Σ_{k≥1} (−1)^k x^(2k)/(2k (2k)!)`;
//!
//! each summed until the terms halve and fall below the precision of the
//! largest, the rest at most the last term. The terms reach about `e^|x|`,
//! so `Si` and `Ci` take `|x| log2 e` more bits, and `Ei` of a negative `x`,
//! about `e^x/|x|`, twice that.
//!
//! From `T` on, the asymptotic expansions:
//!
//! - `Ei(x) = −E₁(z)` for `x = −z < 0`, `E₁(z) = e^(−z)/z (Σ_{k<n} (−1)^k
//!   k!/z^k + R)` with `|R| ≤ n!/z^n`: `E₁(z) = e^(−z) ∫_0^∞ e^(−zt)/(1 +
//!   t) dt`, and `1/(1 + t)` expanded leaves `(−t)^n/(1 + t)`, at most `t^n`;
//! - `Ei(x) = e^x/x (Σ_{k<n} k!/x^k + R)` for `x > 0` and `n ≤ x`, with
//!   `|R| ≤ n!/x^n + 4x² ((n − 1)/(e x))^(n−1)`: `e^(−x) Ei(x) = PV
//!   ∫_0^∞ e^(−xt)/(1 − t) dt`, and `1/(1 − t)` expanded leaves `PV ∫_0^∞
//!   g(t)/(1 − t) dt`, `g = e^(−xt) t^n`; past 2, `|1 − t| ≥ 1` bounds that
//!   part by `∫ g = n!/x^(n+1)`, and on `[0, 2]` the points `1 ± u` pair up,
//!   `|g(1 − u) − g(1 + u)| ≤ 2u sup |g′|`, `|g′| ≤ 2x ((n − 1)/(e x))^(n−1)`;
//! - `Si(x) = π/2 − f cos x − g sin x` and `Ci(x) = f sin x − g cos x` for
//!   `x > 0` (DLMF 6.5.3–4), `f = Σ_{k<n} (−1)^k (2k)!/x^(2k+1)` and `g =
//!   Σ_{k<n} (−1)^k (2k + 1)!/x^(2k+2)`, each within its first term left
//!   out: they are the Laplace transforms of `1/(1 + t²)` and `t/(1 + t²)`
//!   (DLMF 6.7.13–14), expanded with remainders at most `t^(2n)` and
//!   `t^(2n+1)`.
//!
//! `Si` is odd. `Ci` is for `x > 0`, and `Ei` and `Ci` are indeterminate on
//! a ball that meets 0, as `Ei` is past `exp`'s range. A radius `ρ` adds
//! `ρ` times `e^x/|x|`, `min(1, 1/|x|)` or `1/|x|`, bounded from the ball's
//! ends. An exact input is computed again at more bits while its relative
//! accuracy falls short, up to `4 (wp + bits) + 256` more: `Ei` and `Ci`
//! have zeros (`Ei(0.3725…) = 0`), next to which a fixed precision is only
//! absolute.

use std::f64::consts::{LN_2, LOG2_E};

use super::exp::HUGE_EXP;
use super::{exponent, one, refined, upper};
use crate::Integer;
use crate::ball::{Ball, Float, Mag, Round};

/// The arguments from which the asymptotic expansions serve `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "a threshold")]
fn threshold(wp: u64) -> f64 {
    let w = wp as f64 + 16.0;
    w * LN_2 + 2.0 * w.ln() + 8.0
}

/// `|x|`, saturated past `f64`'s range.
fn magnitude(x: &Float) -> f64 {
    if exponent(x) > 1000 {
        return f64::INFINITY;
    }
    x.to_f64().abs()
}

/// Extra bits for a series whose terms reach `2^bits` above its result.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a small positive estimate"
)]
fn lift(bits: f64) -> u64 {
    bits.max(0.0).ceil() as u64
}

fn small_integer(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

/// Whether the halving of the terms has started and `term` is below `2^-q`
/// of the largest term.
fn negligible(halving: bool, term: &Mag, largest: &Mag, q: u64) -> bool {
    halving && *term <= largest.mul(&Mag::power_of_two(&-Integer::from(q)))
}

/// `Ei(x)` for an exact `x ≠ 0` by the series, at about `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn ei_series(x: &Float, wp: u64) -> Ball {
    let ax = magnitude(x);
    let above = if x.is_negative() {
        2.0 * ax * LOG2_E + ax.max(1.0).log2()
    } else {
        0.0
    };
    let q = wp + lift(above) + 32;
    let xb = Ball::exact(x.clone());
    // x^k/k! and its term x^k/(k k!).
    let mut power = one();
    let mut sum = Ball::ZERO;
    let mut largest = Mag::ZERO;
    let mut k: u64 = 0;
    let rest = loop {
        k += 1;
        let kb = small_integer(k);
        power = power.mul(&xb, q).div(&kb, q);
        let term = power.div(&kb, q);
        sum = sum.add(&term, q);
        let size = upper(&term);
        if size > largest {
            largest = size.clone();
        }
        // |x| k/(k + 1)² ≤ 1/2 from k + 1 ≥ 2|x|: the rest is at most this
        // term.
        if negligible((k + 1) as f64 >= 2.0 * ax, &size, &largest, q) {
            break size;
        }
    };
    let log = Ball::exact(x.abs()).log(q);
    Ball::const_euler(q)
        .add(&log, q)
        .add(&sum.add_error(&rest), q)
}

/// `4x² ((n − 1)/(e x))^(n−1)`, from above, for `n ≥ 1` and `x > 0`.
fn principal_value_bound(x: &Float, n: u64) -> Mag {
    let x_up = Mag::from_float(x);
    let square = x_up.mul(&x_up).mul_2exp(&Integer::from(2));
    if n == 1 {
        return square;
    }
    // e > 2718/1000.
    let e_low = Mag::from_integer_lower(&Integer::from(2718))
        .div_lower(&Mag::from_integer(&Integer::from(1000)));
    let ratio =
        Mag::from_integer(&Integer::from(n - 1)).div(&Mag::from_float_lower(x).mul_lower(&e_low));
    square.mul(&ratio.pow(n - 1))
}

/// `Ei(x)` for an exact `x ≠ 0` by the asymptotic expansions, at about `wp`
/// bits; `None` where their terms grow before they reach it.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn ei_asymptotic(x: &Float, wp: u64) -> Option<Ball> {
    let q = wp + 32;
    let small = Mag::power_of_two(&-Integer::from(wp + 12));
    let negative = x.is_negative();
    let ax = magnitude(x);
    let z = Ball::exact(x.abs());
    // Σ_{k<n} (±1)^k k!/|x|^k, the next term and the bound on the rest.
    let mut sum = one();
    let mut term = one();
    let mut k: u64 = 0;
    let bound = loop {
        k += 1;
        // The terms grow past |x|, and the bound for x > 0 needs n ≤ x.
        if k as f64 > ax {
            return None;
        }
        let next = term.mul(&small_integer(k), q).div(&z, q);
        let next = if negative { next.neg() } else { next };
        // n = k terms are in the sum, and `next` is the first left out.
        let mut bound = upper(&next);
        if !negative {
            bound = bound.add(&principal_value_bound(x, k));
        }
        if bound <= small {
            break bound;
        }
        sum = sum.add(&next, q);
        term = next;
    };
    let scale = Ball::exact(x.clone()).exp(q).div(&z, q);
    let value = scale.mul(&sum.add_error(&bound), q);
    Some(if negative { value.neg() } else { value })
}

/// `Ei(x)` for an exact `x ≠ 0` with exponent at most `HUGE_EXP` if
/// positive.
fn ei_exact(x: &Float, wp: u64) -> Ball {
    if magnitude(x) >= threshold(wp)
        && let Some(value) = ei_asymptotic(x, wp)
    {
        return value;
    }
    ei_series(x, wp)
}

/// `Si(x)` or, for `cosine`, `Ci(x)` for an exact `x` (`x > 0` for `Ci`)
/// by the series, at about `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn sin_cos_series(x: &Float, wp: u64, cosine: bool) -> Ball {
    let ax = magnitude(x);
    let q = wp + lift(ax * LOG2_E + ax.max(1.0).log2()) + 32;
    let minus_square = Ball::exact(x.mul(x).neg());
    // x^j/j! for the j of the terms, odd for Si and even for Ci, and the
    // term x^j/(j j!) with its sign.
    let (mut power, mut sum, mut j) = if cosine {
        (one(), Ball::ZERO, 0u64)
    } else {
        let xb = Ball::exact(x.clone());
        (xb.clone(), xb, 1u64)
    };
    let mut largest = upper(&sum);
    let rest = loop {
        j += 2;
        power = power
            .mul(&minus_square, q)
            .div(&small_integer((j - 1) * j), q);
        let term = power.div(&small_integer(j), q);
        sum = sum.add(&term, q);
        let size = upper(&term);
        if size > largest {
            largest = size.clone();
        }
        // The next term is at most x²/((j + 1)(j + 2)) of this one.
        let halving = ((j + 1) * (j + 2)) as f64 >= 2.0 * ax * ax;
        if negligible(halving, &size, &largest, q) {
            break size;
        }
    };
    let sum = sum.add_error(&rest);
    if !cosine {
        return sum;
    }
    let log = Ball::exact(x.clone()).log(q);
    Ball::const_euler(q).add(&log, q).add(&sum, q)
}

/// The auxiliary `f(x)` and `g(x)` for an exact `x > 0` by their
/// asymptotic expansions, at about `wp` bits; `None` where the terms grow
/// before they reach it.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn auxiliary(x: &Float, wp: u64) -> Option<(Ball, Ball)> {
    let q = wp + 32;
    let small = Mag::power_of_two(&-Integer::from(wp + 12));
    let ax = magnitude(x);
    let xb = Ball::exact(x.clone());
    // b_j = j!/x^j: f takes the even j, g the odd ones, signs alternating.
    let mut b = one();
    let mut f = one();
    let mut g = Ball::ZERO;
    let mut j: u64 = 0;
    let bound = loop {
        j += 1;
        // The first terms left out are b_j and b_(j+1) ≤ b_j, for j < x.
        if (j + 1) as f64 > ax {
            return None;
        }
        b = b.mul(&small_integer(j), q).div(&xb, q);
        if upper(&b) <= small {
            break upper(&b);
        }
        let signed = if (j / 2) % 2 == 0 { b.clone() } else { b.neg() };
        if j % 2 == 0 {
            f = f.add(&signed, q);
        } else {
            g = g.add(&signed, q);
        }
    };
    let f = f.add_error(&bound).div(&xb, q);
    let g = g.add_error(&bound).div(&xb, q);
    Some((f, g))
}

/// `Si(x)` or `Ci(x)` for an exact `x > 0` by the asymptotic expansions.
fn sin_cos_asymptotic(x: &Float, wp: u64, cosine: bool) -> Option<Ball> {
    let q = wp + 32;
    let (f, g) = auxiliary(x, wp)?;
    let (s, c) = Ball::exact(x.clone()).sin_cos(q);
    Some(if cosine {
        f.mul(&s, q).sub(&g.mul(&c, q), q)
    } else {
        let half_pi = Ball::const_pi(q).mul_2exp(&-Integer::ONE);
        half_pi.sub(&f.mul(&c, q), q).sub(&g.mul(&s, q), q)
    })
}

/// `Si(x)` for an exact `x ≥ 0`, or `Ci(x)` for an exact `x > 0`.
fn sin_cos_exact(x: &Float, wp: u64, cosine: bool) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    if magnitude(x) >= threshold(wp)
        && let Some(value) = sin_cos_asymptotic(x, wp, cosine)
    {
        return value;
    }
    sin_cos_series(x, wp, cosine)
}

/// A lower bound of `|x|` over the ball, 0 where it meets 0.
fn distance_to_zero(b: &Ball) -> Mag {
    Mag::from_float_lower(&b.mid().abs()).sub_lower(b.rad())
}

/// The upper end of the ball, rounded up.
fn upper_end(b: &Ball) -> Option<Float> {
    let r = b.rad().to_float()?;
    Some(b.mid().add_round(&r, 64, Round::Ceil).0)
}

impl Ball {
    /// `Ei(x)` for every `x` in the ball, at `prec` bits; indeterminate on a
    /// ball meeting 0, and past `exp`'s range.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn ei(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::ei: precision 0");
        let distance = distance_to_zero(self);
        if !self.is_finite() || distance.is_zero() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        if !m.is_negative() && exponent(m) > HUGE_EXP {
            return Ball::INDETERMINATE;
        }
        let value = refined(self, prec, |wp| ei_exact(m, wp));
        if self.is_exact() {
            return value.round_to_accuracy(prec);
        }
        // ρ sup e^x/|x|: e^x at the upper end, 1/|x| at the end next to 0.
        let Some(hi) = upper_end(self) else {
            return Ball::INDETERMINATE;
        };
        let spread = self
            .rad()
            .mul(&upper(&Ball::exact(hi).exp(64)))
            .div(&distance);
        value.add_error(&spread).round_to_accuracy(prec)
    }

    /// `Si(x)` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn si(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::si: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        let value = refined(self, prec, |wp| {
            let v = sin_cos_exact(&m.abs(), wp, false);
            if m.is_negative() { v.neg() } else { v }
        });
        // |sin x/x| ≤ min(1, 1/|x|).
        let distance = distance_to_zero(self);
        let slope = if distance >= Mag::ONE {
            Mag::ONE.div(&distance)
        } else {
            Mag::ONE
        };
        value
            .add_error(&self.rad().mul(&slope))
            .round_to_accuracy(prec)
    }

    /// `Ci(x)` for every `x` in the ball, at `prec` bits; indeterminate on
    /// a ball that reaches 0 or below.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn ci(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::ci: precision 0");
        let distance = distance_to_zero(self);
        if !self.is_finite() || self.mid().is_negative() || distance.is_zero() {
            return Ball::INDETERMINATE;
        }
        let value = refined(self, prec, |wp| sin_cos_exact(self.mid(), wp, true));
        // |cos x/x| ≤ 1/|x|.
        value
            .add_error(&self.rad().div(&distance))
            .round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
