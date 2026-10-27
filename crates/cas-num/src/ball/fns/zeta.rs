//! The Riemann zeta function of real balls (CAS-DOC-15 §7, CAS-244).
//!
//! - `σ ≥ 1/2`: Borwein's algorithm \[Bor00, Algorithm 2\], `ζ(σ) = −1/(d_n
//!   (1 − 2^(1−σ))) Σ_{k<n} (−1)^k (d_k − d_n)/(k + 1)^σ + γ_n`, with `d_k =
//!   n Σ_{i≤k} (n + i − 1)! 4^i/((n − i)! (2i)!)` and `|γ_n| ≤ 3/((3 +
//!   √8)^n |1 − 2^(1−σ)|)` for a real `σ ≥ 1/2`: `n` about `0.39 wp` terms.
//!   The powers `m^(−σ)` are exponentials only at the primes, products of
//!   the smaller ones elsewhere; `1 − 2^(1−σ)` is `−expm1((1 − σ) log 2)`,
//!   whose relative accuracy holds next to the pole.
//! - `σ` above the precision: `1 + [0, 3·2^(−min(σ, wp + 64))]`, `Σ_{k≥2}
//!   k^(−σ) ≤ 2^(−σ) (1 + 2/(σ − 1))`.
//! - `s < 1/2`: the functional equation `ζ(s) = 2^s π^(s−1) sin(πs/2) Γ(1 −
//!   s) ζ(1 − s)`, the factors but the sine and `ζ(1 − s)` as one
//!   exponential of their logarithms (CAS-242's `log Γ`).
//! - Exactly: `ζ(0) = −1/2`, 0 at the negative even integers; a ball that
//!   meets the pole 1, or straddles `1/2`, is indeterminate, and so is one
//!   below `−2^1024`, where `|ζ|` grows as `|s|^|s|`.
//!
//! \[Bor00\] Peter Borwein, "An efficient algorithm for the Riemann zeta
//! function", *Canadian Mathematical Society Conference Proceedings* 27
//! (2000), 29–34.

use std::cmp::Ordering;

use super::gamma::log_gamma_positive;
use super::{exponent, one, working};
use crate::ball::{Ball, Float, Mag, Round};
use crate::{Integer, Rational};

/// `1 − x` for a ball: exact for an exact `x`, rounded at `q` bits.
fn one_minus(x: &Ball, q: u64) -> Ball {
    if x.is_exact() {
        Ball::exact(Float::from(1i64).sub(x.mid()))
    } else {
        one().sub(x, q)
    }
}

/// `m^(−σ)` for `m = 1 … n` (index `m − 1`): exponentials at the primes,
/// products elsewhere.
fn powers(n: u64, sigma: &Ball, q: u64) -> Vec<Ball> {
    let n = usize::try_from(n).expect("few terms");
    // The smallest prime factor of each m ≤ n.
    let mut least = vec![0usize; n + 1];
    for p in 2..=n {
        if least[p] == 0 {
            let mut k = p;
            while k <= n {
                if least[k] == 0 {
                    least[k] = p;
                }
                k += p;
            }
        }
    }
    let mut out: Vec<Ball> = Vec::with_capacity(n);
    out.push(one());
    for m in 2..=n {
        let p = least[m];
        let power = if p == m {
            Ball::from_integer(&Integer::from(m))
                .log(q)
                .mul(sigma, q)
                .neg()
                .exp(q)
        } else {
            out[p - 1].mul(&out[m / p - 1], q)
        };
        out.push(power);
    }
    out
}

/// `ζ(σ)` for a ball `σ ≥ 1/2`, `one_minus_sigma` its `1 − σ`, by Borwein's
/// algorithm at about `wp` bits.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    reason = "a term count from the precision"
)]
fn borwein(sigma: &Ball, one_minus_sigma: &Ball, wp: u64) -> Ball {
    // (3 + √8)^n ≥ 2^(wp + 10): log2(3 + √8) > 2.54.
    let n = ((wp as f64 + 10.0) / 2.54).ceil() as u64 + 1;
    let q = wp + 2 * u64::from(u64::BITS - n.leading_zeros()) + 16;
    // d_0 … d_n, from u_0 = 1, u_(i+1) = u_i 2 (n + i)(n − i)/((i + 1)(2i + 1)).
    let mut d = Vec::with_capacity(usize::try_from(n).expect("few terms") + 1);
    let mut u = Rational::ONE;
    d.push(Rational::ONE);
    for i in 0..n {
        let num = Integer::from(2 * (n + i)) * Integer::from(n - i);
        let den = Integer::from(i + 1) * Integer::from(2 * i + 1);
        u = &u * &Rational::new(num, den);
        let next = d.last().expect("d_0") + &u;
        d.push(next);
    }
    let dn = d.last().expect("d_n").clone();
    let power = powers(n, sigma, q);
    let mut sum = Ball::ZERO;
    for (k, p) in power.iter().enumerate() {
        let term = Ball::from_rational(&(&d[k] - &dn), q).mul(p, q);
        sum = if k % 2 == 0 {
            sum.add(&term, q)
        } else {
            sum.sub(&term, q)
        };
    }
    // 1 − 2^(1−σ) = −expm1((1 − σ) log 2).
    let factor = one_minus_sigma.mul(&Ball::const_log2(q), q).expm1(q).neg();
    let low = Mag::from_float_lower(&factor.mid().abs()).sub_lower(factor.rad());
    if low.is_zero() {
        return Ball::INDETERMINATE;
    }
    let value = sum
        .neg()
        .div(&Ball::from_rational(&dn, q).mul(&factor, q), q);
    // 3/((3 + √8)^n |1 − 2^(1−σ)|), with (3 + √8)^n ≥ 2^(2.54 n) as `n` was
    // chosen: 2^(5n/2) fell short of the precision by 0.04 n bits, past the
    // working bits from 3000 bits on (CAS-248).
    let tail = Mag::from_integer(&Integer::from(3))
        .mul(&Mag::power_of_two(&-Integer::from(254 * n / 100)))
        .div(&low);
    value.add_error(&tail)
}

/// `ζ(s)` for a ball below `1/2` by the functional equation.
fn functional(s: &Ball, wp: u64) -> Ball {
    // The logarithms are about |s| log |s|: their bits above the point.
    let above = u64::try_from(exponent(s.mid()).max(0)).unwrap_or(0);
    let q = wp + above + u64::from(u64::BITS - above.leading_zeros()) + 16;
    let complement = one_minus(s, q);
    let log_gamma = log_gamma_positive(&complement, q);
    // s log 2 + (s − 1) log π + log Γ(1 − s).
    let log_pi = Ball::const_pi(q).log(q);
    let exponent_ball = s
        .mul(&Ball::const_log2(q), q)
        .add(&s.sub(&one(), q).mul(&log_pi, q), q)
        .add(&log_gamma, q);
    let magnitude = exponent_ball.exp(q);
    let sine = s.mul_2exp(&-Integer::ONE).sin_pi(q);
    let rest = borwein(&complement, s, q);
    magnitude.mul(&sine, q).mul(&rest, q)
}

impl Ball {
    /// `ζ(s)` for every real `s` in the ball, at `prec` bits; indeterminate
    /// where the ball meets the pole 1 or straddles `1/2`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn zeta(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::zeta: precision 0");
        if !self.is_finite() || self.contains(&Float::from(1i64)) {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        if self.is_exact() {
            if m.is_zero() {
                return Ball::exact(Float::from(-1i64).mul_2exp(&-Integer::ONE));
            }
            // The trivial zeros: an even integer's lowest bit is 2 or above,
            // found without a rational, which a huge one has not (CAS-249).
            if m.is_negative() && *m.exponent() > Integer::from(m.bit_length()) {
                return Ball::ZERO;
            }
        }
        let wp = working(prec);
        let Some(r) = self.rad().to_float() else {
            return Ball::INDETERMINATE;
        };
        let lo = m.sub_round(&r, 64, Round::Floor).0;
        let hi = m.add_round(&r, 64, Round::Ceil).0;
        let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
        let value = if lo.cmp(&half) != Ordering::Less {
            // Above the precision: 1 + [0, 3·2^(−lo)].
            let large =
                exponent(&lo) > 64 || lo.to_rational() > Rational::from(Integer::from(wp + 4));
            if large {
                // 2^(−⌊lo⌋), the exponent at most wp + 64: still a bound, and
                // one a rational holds (3·2^(−2^40) was not, CAS-249).
                let cap = Integer::from(wp + 64);
                let e = if exponent(&lo) > 64 {
                    cap
                } else {
                    lo.to_rational().floor().min(cap)
                };
                let bound = Mag::from_integer(&Integer::from(3)).mul(&Mag::power_of_two(&-e));
                return one().add_error(&bound).round_to_accuracy(prec);
            }
            borwein(self, &one_minus(self, wp + 16), wp)
        } else if hi.cmp(&half) == Ordering::Less {
            // |ζ(s)| grows as |s|^|s|: past 2^1024 it is out of every range.
            if exponent(m) > 1024 {
                return Ball::INDETERMINATE;
            }
            functional(self, wp)
        } else {
            return Ball::INDETERMINATE;
        };
        value.round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
