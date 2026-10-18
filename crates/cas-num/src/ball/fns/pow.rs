//! `pow` with a ball exponent and `root` (CAS-DOC-5 §7.3).
//!
//! `x^y` takes exact fast paths first: an exact integer `y` is
//! [`Ball::pow`] (any sign of `x`), an exact half-integer `n + 1/2` is
//! `sqrt(x)^(2n+1)`, an exact 0 base gives 0 or 1. Otherwise it is
//! `exp(y log x)` for `x > 0`, with `log x` at as many more bits as `y log x`
//! has above its point, so that `exp` sees an absolute error below its own
//! ulp; the radii of both balls propagate through the operations.
//!
//! `root(x, n)` is `exp(log(x) / n)` for `x > 0` and `-root(-x, n)` for odd
//! `n` and `x < 0` (the real `n`-th root, CAS-DOC-4 §6.1). An odd root of a
//! ball across 0 is `[0 ± R]` with `R` the larger root of the two ends; an
//! even root of a ball reaching below 0 is indeterminate. An exact `x` whose
//! root is exact (`c^n = x` for the midpoint `c` rounded to the root's
//! length) gives the exact ball.

use std::cmp::Ordering;

use super::working;
use crate::Integer;
use crate::ball::decimal::to_integer;
use crate::ball::{Ball, Float, Mag, Round};

/// `x` as an `i64` if it is an exact integer that fits.
fn as_i64(x: &Ball) -> Option<i64> {
    let m = x.mid();
    if !x.is_exact() {
        return None;
    }
    if m.is_zero() {
        return Some(0);
    }
    // |m| < 2^63 first, so that the rounded integer is small.
    if *m.exponent() > Integer::from(63i64) {
        return None;
    }
    let v = to_integer(m, false);
    (Float::from(&v) == *m).then(|| v.to_i64()).flatten()
}

/// `x^y` for a ball `x > 0` and any `y`, as `exp(y log x)`.
fn exp_log(x: &Ball, y: &Ball, prec: u64) -> Ball {
    let wp = working(prec);
    let mut extra = 0;
    loop {
        let l = x.log(wp + extra);
        let t = y.mul(&l, wp + extra);
        if !t.is_finite() {
            return Ball::INDETERMINATE;
        }
        // Bits of y log x above its point, which exp turns into relative
        // error of the result.
        let above = t
            .mid()
            .exponent()
            .to_i64()
            .map_or(0, |e| u64::try_from(e).unwrap_or(0));
        if t.mid().is_zero() || above <= extra {
            return t.exp(prec);
        }
        extra = above + 8;
    }
}

impl Ball {
    /// `x^y` for every `x` in `self` and `y` in `y`, at `prec` bits: exact
    /// fast paths for integer and half-integer `y`, `exp(y log x)` otherwise
    /// (which needs `x > 0`; otherwise indeterminate).
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn pow_ball(&self, y: &Ball, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::pow_ball: precision 0");
        if !self.is_finite() || !y.is_finite() {
            return Ball::INDETERMINATE;
        }
        // Squaring and multiplying rounds about 2 log2(n) times: that many
        // more bits, then one rounding.
        let wp =
            |n: i64| working(prec) + 2 * u64::from(u64::BITS - n.unsigned_abs().leading_zeros());
        if let Some(n) = as_i64(y) {
            return self.pow(n, wp(n)).round_to_accuracy(prec);
        }
        // Half-integers n + 1/2: sqrt(x)^(2n + 1), for x >= 0.
        if let Some(k) = as_i64(&y.mul_2exp(&Integer::ONE)) {
            let p = wp(k);
            return self.sqrt(p).pow(k, p).round_to_accuracy(prec);
        }
        if self.is_exact() && self.mid().is_zero() {
            // 0^y: 0 for y > 0, undefined otherwise.
            return match y.sign() {
                Some(Ordering::Greater) => Ball::ZERO,
                _ => Ball::INDETERMINATE,
            };
        }
        if self.sign() != Some(Ordering::Greater) {
            return Ball::INDETERMINATE;
        }
        exp_log(self, y, prec)
    }

    /// The real `n`-th root of every point of the ball, at `prec` bits: for
    /// odd `n` negative values have negative roots; an even root of a ball
    /// reaching below 0 is indeterminate. Exact for an exact perfect power.
    ///
    /// # Panics
    ///
    /// If `n` or `prec` is 0.
    #[must_use]
    pub fn root(&self, n: u64, prec: u64) -> Ball {
        assert!(n >= 1, "Ball::root: 0th root");
        assert!(prec >= 1, "Ball::root: precision 0");
        if n == 1 {
            return self.round_to_accuracy(prec);
        }
        if n == 2 {
            return self.sqrt(prec);
        }
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let odd = n % 2 == 1;
        match self.sign() {
            Some(Ordering::Equal) => Ball::ZERO,
            Some(Ordering::Greater) => root_positive(self, n, prec),
            Some(Ordering::Less) if odd => root_positive(&self.neg(), n, prec).neg(),
            _ if odd => {
                // Across 0: [-root(r - m), root(m + r)] for the ends m ± r.
                let r = self.rad().to_float().expect("finite");
                let hi = self.mid().add_round(&r, prec + 8, Round::Ceil).0;
                let lo = self.mid().sub_round(&r, prec + 8, Round::Floor).0;
                let up = |x: &Float| {
                    let b = root_positive(&Ball::exact(x.clone()), n, prec + 8);
                    Mag::from_float(b.mid()).add(b.rad())
                };
                // Both ends bound |root| on their side; the ball [0 ± max].
                Ball::new(Float::ZERO, up(&hi).max(up(&lo.neg())))
            }
            _ => Ball::INDETERMINATE,
        }
    }
}

/// The `n`-th root (`n >= 3`) of a ball `x > 0`.
fn root_positive(x: &Ball, n: u64, prec: u64) -> Ball {
    let wp = working(prec);
    let div = Ball::from_integer(&Integer::from(n));
    let r = x.log(wp).div(&div, wp).exp(prec);
    if x.is_exact() && r.is_finite() {
        // A perfect power: c^n = x exactly for c the midpoint rounded to the
        // root's length.
        let bits = x.mid().bit_length().div_ceil(n) + 2;
        let c = r.mid().round(bits, Round::Nearest).0;
        let cn = Ball::exact(c.clone()).pow(
            i64::try_from(n).unwrap_or(i64::MAX),
            x.mid().bit_length() + 64,
        );
        if cn.is_exact() && cn.mid() == x.mid() {
            return Ball::exact(c);
        }
    }
    r
}
