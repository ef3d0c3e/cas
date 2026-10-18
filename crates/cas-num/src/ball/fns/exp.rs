//! `exp` and `expm1` (CAS-DOC-5 §7.2).
//!
//! `exp(m)` for the exact midpoint: `m = k log 2 + t` with `|t| <= log 2 / 2`
//! (a ball, from `log 2` at enough bits for `k`), `t = x 2^h` with `|x| <=
//! 2^-r`, `r ≈ p^(1/3)`, the Taylor series of `exp(x)` summed in fixed point
//! by rectangular splitting ([`exp_sum`]), then squared back `h` times and
//! scaled by `2^k`, as balls. The input radius `ρ` adds `exp(m) (e^ρ - 1)`.
//!
//! Arguments of exponent above `HUGE_EXP` (`|m| >= 2^65536`) would need `log
//! 2` to that many bits: a positive one gives the indeterminate ball, a
//! negative one `[0 ± 2^(-2^(e - 2))]`, a bound on `exp` of every point.

use super::fixed::{exp_sum, shr};
use super::{one, upper, working};
use crate::ball::decimal::to_integer;
use crate::ball::{Ball, Float, Mag, Round};
use crate::{Integer, Natural};

/// Exponents of arguments that `exp` reduces; see the module comment.
pub(super) const HUGE_EXP: i64 = 1 << 16;

/// `x` truncated to a fixed-point integer `trunc(x 2^w)`.
pub(super) fn to_fixed(x: &Float, w: u64) -> Integer {
    let (m, e) = x.to_integer_exp();
    let shift = e + Integer::from(w);
    match shift.to_i64() {
        Some(s) if s >= 0 => {
            let (neg, mag) = m.into_parts();
            Integer::from_parts(neg, mag << s.unsigned_abs())
        }
        Some(s) => shr(m, s.unsigned_abs()),
        None if shift.is_negative() => Integer::ZERO,
        None => unreachable!("an argument below 1"),
    }
}

/// A lower bound of `n!`.
fn factorial_lower(n: u64) -> Mag {
    (2..=n).fold(Mag::ONE, |acc, k| {
        acc.mul_lower(&Mag::from_integer_lower(&Integer::from(k)))
    })
}

/// A ball containing `exp(t)` for an exact `|t| < 1`, at about `wp` bits.
///
/// With `|t| < 2^et`, `h = max(0, r + et)` halvings give `x = t / 2^h` with
/// `|x| < 2^ex`, `ex = et - h <= min(et, -r) <= -1`. In fixed point at `w`
/// bits, `X = trunc(x 2^w)` is within `2^-w` of `x`, which moves `exp` by at
/// most `exp(1/2) (e^(2^-w) - 1) < 4 · 2^-w`. The series to `n` terms leaves a
/// tail of at most `2 |x|^n / n! <= 2^(ex n + 1) / n!` (`|x| <= 1/2`), and
/// `n` is the first count with `ex n - log2(n!) <= -(w + 2)`. The `h`
/// squarings double the relative error each, so the sum has `w = wp + 2h +
/// 16` bits and the squarings `wp + h + 8`.
fn exp_small(t: &Float, wp: u64) -> Ball {
    if t.is_zero() {
        return one();
    }
    let et = t.exponent().to_i64().expect("|t| < 1");
    debug_assert!(et <= 0);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        clippy::cast_possible_wrap,
        reason = "a cube root of a precision, a small positive count"
    )]
    let r = ((wp as f64).cbrt().ceil() as i64).max(2);
    let h = (r + et).max(0);
    let ex = et - h;
    let x = t.mul_2exp(&Integer::from(-h));
    let hu = h.unsigned_abs();
    let w = wp + 2 * hu + 16;
    let xf = to_fixed(&x, w);
    // Terms: the first n with n·(-ex) + log2(n!) >= w + 2 (an estimate that
    // only sets the accuracy; the tail is bounded exactly below).
    #[allow(clippy::cast_precision_loss, reason = "an estimate")]
    let n = {
        let (mut n, mut bits) = (1u64, 0.0f64);
        let per = (-ex) as f64;
        while (n as f64) * per + bits < (w + 2) as f64 {
            n += 1;
            bits += (n as f64).log2();
        }
        n
    };
    let (s, err) = exp_sum(&xf, w, n);
    let unit = |k: Natural| Mag::from_integer(&Integer::from(k)).mul_2exp(&-Integer::from(w));
    let tail = Mag::power_of_two(&(Integer::from(ex) * Integer::from(n) + Integer::ONE))
        .div(&factorial_lower(n));
    let rad = unit(&err + &Natural::from(4u64)).add(&tail);
    let mut b = Ball::new(Float::from_integer_exp(&s, &-Integer::from(w)), rad);
    let sq = wp + hu + 8;
    for _ in 0..h {
        b = b.mul(&b, sq);
    }
    b
}

/// A ball containing `exp(m)` for an exact `m` with exponent at most
/// `HUGE_EXP`, at about `wp` bits.
pub(super) fn exp_exact(m: &Float, wp: u64) -> Ball {
    if m.is_zero() {
        return one();
    }
    let e = m.exponent().to_i64().expect("checked by the caller");
    if e <= -1 {
        return exp_small(m, wp);
    }
    // k = round(m / log 2) < 2^(e + 1), from log 2 at 64 bits beyond k's:
    // only |t| needs to be small, not k exact.
    let kbits = e.unsigned_abs() + 2;
    let approx = Ball::const_log2(kbits + 64);
    let k = to_integer(
        &m.div_round(approx.mid(), kbits + 8, Round::Nearest).0,
        false,
    );
    let p = wp + kbits + 8;
    let kl = Ball::const_log2(p).mul(&Ball::from_integer(&k), p);
    let t = Ball::exact(m.clone()).sub(&kl, p);
    // exp(t ± ρ) ⊂ exp(t) (1 ± (e^ρ - 1)).
    let s = exp_small(t.mid(), wp);
    let s = s.add_error(&upper(&s).mul(&t.rad().expm1()));
    s.mul_2exp(&k)
}

/// `exp` of every point of a ball whose midpoint has exponent above
/// `HUGE_EXP`: indeterminate unless the ball lies below `-2^(e-2)`.
fn exp_huge(b: &Ball, e: &Integer) -> Ball {
    let m = b.mid();
    // r < |m| / 2 keeps every point at most m / 2 <= -2^(e - 2).
    if !m.is_negative() || Mag::from_float(m).mul_2exp(&-Integer::ONE) <= *b.rad() {
        return Ball::INDETERMINATE;
    }
    // exp(x) <= 2^x for x <= 0. The bound's exponent is capped so that it
    // stays small; any bound below the true one would do.
    let k = (e - &Integer::from(2i64))
        .to_i64()
        .map_or(1 << 30, |k| k.min(1 << 30));
    let bound = Integer::from_parts(true, Natural::from(1u64) << k.unsigned_abs());
    Ball::new(Float::ZERO, Mag::power_of_two(&bound))
}

impl Ball {
    /// `exp(x)` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn exp(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::exp: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        if m.is_zero() && self.is_exact() {
            return one();
        }
        if !m.is_zero() && *m.exponent() > Integer::from(HUGE_EXP) {
            return exp_huge(self, m.exponent());
        }
        let core = exp_exact(m, working(prec));
        // |exp(m + δ) - exp(m)| = exp(m) |e^δ - 1| <= exp(m) (e^ρ - 1).
        let spread = upper(&core).mul(&self.rad().expm1());
        core.add_error(&spread).round_to_accuracy(prec)
    }

    /// `exp(x) - 1` for every `x` in the ball, at `prec` bits, without
    /// cancellation near 0: relative to the result, not to 1.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn expm1(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::expm1: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        if m.is_zero() {
            // |expm1(x)| <= expm1(|x|).
            return Ball::new(Float::ZERO, self.rad().expm1());
        }
        if *m.exponent() > Integer::from(HUGE_EXP) {
            return exp_huge(self, m.exponent()).sub(&one(), prec);
        }
        let wp = working(prec);
        let core = expm1_exact(m, wp);
        // exp(m) = 1 + expm1(m), and the spread as for exp.
        let spread = upper(&core.add(&one(), 64)).mul(&self.rad().expm1());
        core.add_error(&spread).round_to_accuracy(prec)
    }
}

/// A ball containing `exp(m) - 1` for an exact non-zero `m`, relative to the
/// result. From `|m| >= 1/4` the subtraction of 1 cancels at most 3 bits;
/// below, `exp` is taken at `-exp(m)` more bits to cover the cancellation
/// (`|expm1(m)| >= |m| / 2`); below `2^-wp`, `expm1(m)` lies in `[m, m + m^2]`
/// (`e^m - 1 - m = m^2/2 + m^3/6 + … ∈ [0, m^2]` for `|m| <= 1`).
pub(super) fn expm1_exact(m: &Float, wp: u64) -> Ball {
    let e = m.exponent().to_i64().expect("checked by the caller");
    let wpi = i64::try_from(wp).unwrap_or(i64::MAX);
    if e >= -1 {
        return exp_exact(m, wp + 4).sub(&one(), wp);
    }
    if e < -wpi {
        let sq = Mag::from_float(m).mul(&Mag::from_float(m));
        return Ball::new(m.clone(), sq);
    }
    exp_exact(m, wp + e.unsigned_abs() + 4).sub(&one(), wp)
}
