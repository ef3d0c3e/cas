//! Bessel functions of real balls (CAS-DOC-34 §2, CAS-PLAN-28).
//!
//! - [`expansions`]: the asymptotic expansions for large arguments, for a
//!   real order `ν ≥ 0`, with the remainders DLMF proves (§2.1, §2.2,
//!   CAS-377).
//! - [`jy`]: `J_n` and `Y_n` of integer order (§2.1, CAS-378).
//! - [`ik`]: `I_n` and `K_n` of integer order (§2.2, CAS-379).
//!
//! Orders are integers up to [`MAX_ORDER`] in size: `J_(−n) = (−1)^n J_n`
//! and `Y_(−n) = (−1)^n Y_n`. An exact argument runs the Ziv loop; a ball
//! adds its radius times a bound of the derivative over it: `|J_n′| ≤ 1`,
//! as `|J_k| ≤ 1` for real arguments (DLMF 10.14.1) and `J_n′ = (J_(n−1)
//! − J_(n+1))/2`; `|Y_n′| ≤ (M_(n−1) + M_(n+1))/2` at the ball's lower
//! end, `M_k = √(J_k² + Y_k²)` falling with `x`; `|I_n′| ≤ I_(|n−1|)` at the
//! ball's largest `|x|`, as `I_n′ = (I_(n−1) + I_(n+1))/2` and `I_k` grows
//! with `|x|` and falls with `k`; `|K_n′| ≤ K_(n+1)` at its lower end, as
//! `K_n′ = −(K_(n−1) + K_(n+1))/2` and `K_k` falls with `x` and grows with
//! `k`. `Y_n` and `K_n` are indeterminate unless the ball lies above 0.

pub(super) mod expansions;
mod ik;
mod jy;

use super::{refined, upper, working};
use crate::ball::{Ball, Float, Round};

/// The largest order in size.
pub(super) const MAX_ORDER: u64 = 1024;

/// The lower end of a ball, rounded down; `None` unless it lies above 0.
fn positive_lower_end(b: &Ball) -> Option<Float> {
    let r = b.rad().to_float()?;
    let lo = b.mid().sub_round(&r, 64, Round::Floor).0;
    (!lo.is_negative() && !lo.is_zero()).then_some(lo)
}

/// `(−1)^n` applied to `b` for an odd negative `n`.
fn reflected(n: i64, b: Ball) -> Ball {
    if n < 0 && n % 2 != 0 { b.neg() } else { b }
}

impl Ball {
    /// `J_n(x)`, the Bessel function of the first kind, for every `x` in
    /// the ball, at `prec` bits; indeterminate for `|n|` above
    /// [`MAX_ORDER`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn bessel_j(&self, n: i64, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::bessel_j: precision 0");
        let m = n.unsigned_abs();
        if !self.is_finite() || m > MAX_ORDER {
            return Ball::INDETERMINATE;
        }
        let value = if self.is_exact() {
            refined(self, prec, |wp| jy::j_exact(m, self.mid(), wp))
        } else {
            jy::j_exact(m, self.mid(), working(prec)).add_error(self.rad())
        };
        reflected(n, value).round_to_accuracy(prec)
    }

    /// `Y_n(x)`, the Bessel function of the second kind, for every `x` in
    /// the ball, at `prec` bits; indeterminate unless the ball lies above
    /// 0, and for `|n|` above [`MAX_ORDER`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn bessel_y(&self, n: i64, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::bessel_y: precision 0");
        let m = n.unsigned_abs();
        if !self.is_finite() || m > MAX_ORDER {
            return Ball::INDETERMINATE;
        }
        let Some(lo) = positive_lower_end(self) else {
            return Ball::INDETERMINATE;
        };
        let value = if self.is_exact() {
            refined(self, prec, |wp| jy::y_exact(m, self.mid(), wp))
        } else {
            let slope = jy::modulus(m.abs_diff(1), &lo)
                .add(&jy::modulus(m + 1, &lo))
                .mul_2exp(&-crate::Integer::ONE);
            jy::y_exact(m, self.mid(), working(prec)).add_error(&self.rad().mul(&slope))
        };
        reflected(n, value).round_to_accuracy(prec)
    }

    /// `I_n(x)`, the modified Bessel function of the first kind, for every
    /// `x` in the ball, at `prec` bits; indeterminate for `|n|` above
    /// [`MAX_ORDER`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn bessel_i(&self, n: i64, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::bessel_i: precision 0");
        let m = n.unsigned_abs();
        if !self.is_finite() || m > MAX_ORDER {
            return Ball::INDETERMINATE;
        }
        if self.is_exact() {
            return refined(self, prec, |wp| ik::i_exact(m, self.mid(), wp))
                .round_to_accuracy(prec);
        }
        // I_(−n) = I_n; the derivative's bound at the largest |x|.
        let Some(r) = self.rad().to_float() else {
            return Ball::INDETERMINATE;
        };
        let top = self.mid().abs().add_round(&r, 64, Round::Ceil).0;
        let slope = upper(&ik::i_exact(m.abs_diff(1), &top, 64));
        ik::i_exact(m, self.mid(), working(prec))
            .add_error(&self.rad().mul(&slope))
            .round_to_accuracy(prec)
    }

    /// `K_n(x)`, the modified Bessel function of the second kind, for every
    /// `x` in the ball, at `prec` bits; indeterminate unless the ball lies
    /// above 0, and for `|n|` above [`MAX_ORDER`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn bessel_k(&self, n: i64, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::bessel_k: precision 0");
        let m = n.unsigned_abs();
        if !self.is_finite() || m > MAX_ORDER {
            return Ball::INDETERMINATE;
        }
        let Some(lo) = positive_lower_end(self) else {
            return Ball::INDETERMINATE;
        };
        // K_(−n) = K_n.
        if self.is_exact() {
            return refined(self, prec, |wp| ik::k_exact(m, self.mid(), wp))
                .round_to_accuracy(prec);
        }
        let slope = upper(&ik::k_exact(m + 1, &lo, 64));
        ik::k_exact(m, self.mid(), working(prec))
            .add_error(&self.rad().mul(&slope))
            .round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
