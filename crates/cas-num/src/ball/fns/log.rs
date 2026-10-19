//! `log` and `log1p` by "approximate, then certify" (CAS-DOC-5 §7.2).
//!
//! `y ≈ log(m)` comes from Newton's iteration on `exp` at the float level,
//! `y ← y + m e^-y - 1`, doubling the precision each step: no balls, no
//! error analysis. One ball evaluation then proves it:
//!
//! ```text
//! log(m) = y + log1p(u),   u = m e^-y - 1   (a ball),
//! |log1p(v) - v| <= v^2    for |v| <= 1/2,
//! ```
//!
//! since the series `-v^2/2 + v^3/3 - …` is at most `v^2 / (2 (1 - |v|))`.
//! So `log1p` of the ball `[c ± ρ]` lies in `c ± (ρ + (|c| + ρ)^2)`. The
//! certificate is itself a Newton step whose error `u^2` is proved, so
//! Newton only runs to half the precision: the cost is about one `exp` at
//! half precision and one at full.
//!
//! Arguments within `1/4` of 1 go through `log1p` of the exact `m - 1`,
//! whose Newton step `y ← y + e^-y (x - expm1(y))` keeps the precision
//! relative to the small result, with `e^-y = 1 / (1 + expm1(y))` from the
//! same `expm1`.

use super::exp::{exp_exact, expm1_exact};
use super::{exponent, one, working};
use crate::Integer;
use crate::ball::{Ball, Float, Mag, Round};

/// `log1p` of every point of `u` with `|u| <= 1/2`; `None` beyond that.
fn log1p_small(u: &Ball) -> Option<Ball> {
    let bound = Mag::from_float(u.mid()).add(u.rad());
    if bound > Mag::ONE.mul_2exp(&-Integer::ONE) {
        return None;
    }
    Some(u.add_error(&bound.mul(&bound)))
}

/// Newton precisions from a start accurate to `start` bits up to `target`:
/// each step about doubles the accuracy, so each precision is at most
/// twice the previous one less 4 (`p / 2 + 2` going down, which stops
/// falling at 4), the first at most `2 start`.
fn schedule(start: u64, target: u64) -> Vec<u64> {
    let mut v = vec![target];
    loop {
        let last = *v.last().expect("non-empty");
        if last <= 2 * start.max(4) {
            break;
        }
        v.push(last / 2 + 2);
    }
    v.reverse();
    v
}

/// The bits above the binary point of `y`: precision relative to `y` that
/// keeps an absolute accuracy of `2^-p`.
fn int_bits(y: &Float) -> u64 {
    y.exponent()
        .to_i64()
        .map_or(0, |e| u64::try_from(e).unwrap_or(0))
}

/// `y ≈ log(m)` for `m > 0`, to about `target` bits (absolute, beyond the
/// integer part).
fn log_approx(m: &Float, target: u64) -> Float {
    // y0 = e log 2 + log(f) for m = f 2^e, f in [1/2, 1), to ~50 bits.
    let e = m.exponent();
    let f = m.mul_2exp(&-e).to_f64();
    let l2 = Ball::const_log2(64 + e.bit_length());
    let mut y = Float::from(e)
        .mul_round(l2.mid(), 64 + e.bit_length(), Round::Nearest)
        .0
        .add(&Float::from_f64(f.ln()).expect("finite"));
    // y0 is within about 2^-50 (f64's log and the rounded e log 2).
    for p in schedule(48, target) {
        let p = p + 8 + int_bits(&y);
        let ey = exp_exact(&y.neg(), p);
        let u = m
            .mul_round(ey.mid(), p, Round::Nearest)
            .0
            .sub(&Float::from(1i64));
        y = y.add_round(&u, p, Round::Nearest).0;
    }
    y
}

/// `y ≈ log1p(x)` for exact `x` with `|x| < 1/4`, relative to `|x|`.
fn log1p_approx(x: &Float, target: u64) -> Float {
    // y0 = x - x^2/2 is within about |x|^3/3 of log1p(x), so relative to it
    // within |x|^2/2: 2(-e) - 1 bits for |x| < 2^e (at least 3 bits at e = -2).
    let mut y = x.sub(
        &x.mul_round(x, 64, Round::Nearest)
            .0
            .mul_2exp(&-Integer::ONE),
    );
    let e = x.exponent().to_i64().map_or(1, i64::unsigned_abs);
    for p in schedule((2 * e).saturating_sub(1).max(3), target) {
        let p = p + 8;
        let em1 = expm1_exact(&y, p);
        let d = x.sub_round(em1.mid(), p, Round::Nearest).0;
        // e^-y = 1 / (1 + expm1(y)).
        let den = em1.mid().add_round(&Float::from(1i64), p, Round::Nearest).0;
        y = y
            .add_round(&d.div_round(&den, p, Round::Nearest).0, p, Round::Nearest)
            .0;
    }
    y
}

/// A ball containing `log(m)` for an exact `m > 0`, at about `wp` bits.
pub(super) fn log_exact(m: &Float, wp: u64) -> Ball {
    // Within 1/4 of 1 only for exponents 0 and 1, where m - 1 is exact and
    // short (an exact m - 1 of a huge m would have as many bits as m's
    // exponent).
    if matches!(m.exponent().to_i64(), Some(0 | 1)) {
        let d = m.sub(&Float::from(1i64));
        if d.is_zero() {
            return Ball::ZERO;
        }
        if *d.exponent() <= Integer::from(-2i64) {
            return log1p_exact(&d, wp);
        }
    }
    let p = wp + 16;
    // To half the bits: the certificate's u^2 then stays below 2^-p.
    let y = log_approx(m, p / 2 + 16);
    let q = p + int_bits(&y);
    let ey = exp_exact(&y.neg(), q);
    let u = Ball::exact(m.clone()).mul(&ey, q).sub(&one(), q);
    if let Some(l) = log1p_small(&u) {
        Ball::exact(y).add(&l, q)
    } else {
        debug_assert!(false, "Newton missed log({m:?})");
        Ball::INDETERMINATE
    }
}

/// A ball containing `log1p(x)` for an exact `x > -1`, at about `wp` bits
/// relative to the result. From `|x| >= 1/4` it is `log` of the ball `1 + x`
/// (rounded: an exact sum of a huge `x` and 1 would be as long as its
/// exponent), whose midpoint stays at least `1/4` from 1; below `2^-wp`,
/// `log1p(x) ∈ [x - x^2, x]`.
pub(super) fn log1p_exact(x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    let e = exponent(x);
    if e >= -1 {
        return Ball::exact(x.clone()).add(&one(), wp + 16).log(wp);
    }
    if e < -i64::try_from(wp).unwrap_or(i64::MAX) {
        return Ball::new(x.clone(), Mag::from_float(x).mul(&Mag::from_float(x)));
    }
    let p = wp + 16;
    let y = log1p_approx(x, p / 2 + 16);
    // u = (1 + x) e^-y - 1 = (x - expm1(y)) / (1 + expm1(y)), relative to |x|.
    let em1 = expm1_exact(&y, p);
    let u = Ball::exact(x.clone())
        .sub(&em1, p)
        .div(&em1.add(&one(), p), p);
    if let Some(l) = log1p_small(&u) {
        Ball::exact(y).add(&l, p)
    } else {
        debug_assert!(false, "Newton missed log1p({x:?})");
        Ball::INDETERMINATE
    }
}

impl Ball {
    /// `log(x)` for every `x` in the ball, at `prec` bits; indeterminate
    /// unless the ball lies in `x > 0`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn log(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::log: precision 0");
        let m = self.mid();
        // m - r from below: 0 when the ball reaches 0.
        let gap = Mag::from_float_lower(m).sub_lower(self.rad());
        if !self.is_finite() || m.is_negative() || gap.is_zero() {
            return Ball::INDETERMINATE;
        }
        let core = log_exact(m, working(prec));
        // |log(m + δ) - log(m)| <= r / (m - r).
        core.add_error(&self.rad().div(&gap))
            .round_to_accuracy(prec)
    }

    /// `log(1 + x)` for every `x` in the ball, at `prec` bits, without
    /// cancellation near 0; indeterminate unless the ball lies in `x > -1`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn log1p(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::log1p: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        // 1 + m - r from below, through 1 + m rounded toward -∞.
        let m1 = m.add_round(&Float::from(1i64), 64, Round::Floor).0;
        let gap = Mag::from_float_lower(&m1).sub_lower(self.rad());
        if m1.is_negative() || gap.is_zero() {
            return Ball::INDETERMINATE;
        }
        let core = log1p_exact(m, working(prec));
        // |log1p(m + δ) - log1p(m)| <= r / (1 + m - r).
        core.add_error(&self.rad().div(&gap))
            .round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests {
    use super::schedule;

    #[test]
    fn schedules_end_and_double() {
        for start in [0, 1, 3, 5, 48, 1000] {
            for target in [1, 4, 9, 64, 100, 10_000, 1 << 30] {
                let v = schedule(start, target);
                assert_eq!(*v.last().unwrap(), target);
                assert!(
                    v[0] <= 2 * start.max(4) || v.len() == 1,
                    "{start} {target}: {v:?}"
                );
                for w in v.windows(2) {
                    assert!(w[1] + 4 <= 2 * w[0] + 4 && w[1] > w[0], "{v:?}");
                }
            }
        }
    }
}
