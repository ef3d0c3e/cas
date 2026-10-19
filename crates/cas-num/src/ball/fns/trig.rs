//! `sin`, `cos`, `tan`, `sin_pi`, `cos_pi` (CAS-DOC-5 §7.2).
//!
//! At the exact midpoint `m`: `m = k π/2 + t` with `|t| <= π/4` (a ball, from π
//! at `p + exp(m)` bits, more when `t` cancels: `sin(10^100)` needs about 333
//! more, retrying as in Ziv's strategy \[Ziv91\]), `t = x 2^h` with `|x| <=
//! 2^-r`, `sin x` and `cos x` from their Taylor series in `y = x^2` summed in
//! fixed point ([`hyp_sum`]), then `h` doublings `sin 2a = 2 sin a cos a`, `cos
//! 2a = 1 - 2 sin^2 a` as balls, and the quadrant `k mod 4`; tables or
//! bit-burst replace the halvings and the series at the precisions they serve
//! (see the [module comment](super)). Both are 1-Lipschitz, so the input radius
//! `ρ` adds `min(ρ, 2)`. Arguments of exponent above `HUGE_EXP`, an
//! indeterminate input, or a radius from 2 give `[0 ± 1]`.
//!
//! `sin_pi(x) = sin(π x)` and `cos_pi(x)` reduce `x` modulo 2 exactly (it
//! is dyadic) and fold it onto `[0, 1/2]` by exact differences, so integers
//! and half-integers give exact results, no π is needed for the reduction
//! and values next to 0 keep their relative accuracy; the rest is `sin` or
//! `cos` of `π b` with `0 < b < 1/4`.

use std::cmp::Ordering;

use super::exp::to_fixed;
use super::fixed::{hyp_sum, shr};
use super::{
    BIT_BURST_MIN_PREC, TABLE_MAX_PREC, bitburst, exponent, factorial_lower, one, tables, working,
};
use crate::ball::decimal::to_integer;
use crate::ball::{Ball, Float, Mag, Round};
use crate::{Integer, Natural};

/// Exponents of arguments that `sin` and `cos` reduce.
const HUGE_EXP: i64 = 1 << 16;

/// `[0 ± 1]`, which contains every sine and cosine.
fn unit_ball() -> Ball {
    Ball::new(Float::ZERO, Mag::ONE)
}

/// Balls containing `sin(t)` and `cos(t)` for an exact `|t| < 1`, at about
/// `wp` bits.
///
/// As for `exp`: `h = max(0, r + et)` halvings give `|x| < 2^ex <= 1/2`, so
/// `|y| = x^2 < 2^(2 ex) <= 1/4`. In fixed point at `w` bits, `X` is within
/// `2^-w` of `x` (which moves `sin` and `cos` by at most as much), and `Y =
/// trunc(X^2 / 2^w)` within one unit of `X'^2`; both series have
/// derivatives at most 1 in `y` for `|y| <= 1/4`, so that unit adds one. The
/// tails beyond `n` terms are at most `2 |y|^n / (2n)!` and `2 |y|^n /
/// (2n+1)!`. The `h` doublings double the absolute error each, hence `2h`
/// more bits in the sums and `h` in the doublings.
pub(super) fn sin_cos_small(t: &Float, wp: u64) -> (Ball, Ball) {
    if t.is_zero() {
        return (Ball::ZERO, one());
    }
    let et = exponent(t);
    debug_assert!(et <= 0);
    let wpi = i64::try_from(wp).unwrap_or(i64::MAX);
    if et < -wpi {
        // sin t ∈ t ± |t|^3/6, cos t ∈ [1 - t^2/2, 1].
        let t2 = Mag::from_float(t).mul(&Mag::from_float(t));
        return (
            Ball::new(t.clone(), t2.mul(&Mag::from_float(t))),
            one().add_error(&t2),
        );
    }
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
    let hu = h.unsigned_abs();
    // Fixed point relative to sin t ≈ t: |et| more bits (below 2^-r only,
    // where h = 0).
    let w = wp + 2 * hu + 24 + et.unsigned_abs();
    let xf = to_fixed(&t.mul_2exp(&Integer::from(-h)), w);
    let yf = shr(&xf * &xf, w);
    // Terms: the first n with n·(-2 ex) + log2((2n)!) >= w + 2.
    #[allow(clippy::cast_precision_loss, reason = "an estimate")]
    let n = {
        let (mut n, mut bits) = (1u64, 1.0f64);
        let per = (-2 * ex) as f64;
        while (n as f64) * per + bits < (w + 2) as f64 {
            n += 1;
            bits += ((2 * n - 1) as f64).log2() + ((2 * n) as f64).log2();
        }
        n
    };
    let small = |t: u64| i64::try_from(t).expect("small");
    let (c, ec) = hyp_sum(&yf, w, n, |t| -small((2 * t - 1) * 2 * t));
    let (s, es) = hyp_sum(&yf, w, n, |t| -small(2 * t * (2 * t + 1)));
    let unit = |k: Natural| Mag::from_integer(&Integer::from(k)).mul_2exp(&-Integer::from(w));
    let y_bound = Integer::from(2 * ex) * Integer::from(n) + Integer::ONE;
    let tail_c = Mag::power_of_two(&y_bound).div(&factorial_lower(2 * n));
    let tail_s = Mag::power_of_two(&y_bound).div(&factorial_lower(2 * n + 1));
    let at = |v: Integer| Float::from_integer_exp(&v, &-Integer::from(w));
    let x = Ball::exact(at(xf));
    // cos(X') and sin(X') / X', then sin(X'), each moved by up to 2^-w to x.
    let cos = Ball::new(at(c), unit(&ec + &Natural::from(2u64)).add(&tail_c));
    let sinc = Ball::new(at(s), unit(&es + &Natural::from(1u64)).add(&tail_s));
    let shift = unit(Natural::from(1u64));
    let q = wp + hu + 16;
    let mut sin = x.mul(&sinc, q).add_error(&shift);
    let mut cos = cos.add_error(&shift);
    for _ in 0..h {
        let s2 = sin.mul(&sin, q);
        sin = sin.mul(&cos, q).mul_2exp(&Integer::ONE);
        cos = one().sub(&s2.mul_2exp(&Integer::ONE), q);
    }
    (sin, cos)
}

/// Balls containing `sin(t)` and `cos(t)` for an exact `|t| < 0.789` (just
/// above `π/4`), at about `wp` bits relative to each, from the core for `wp`
/// (see the [module comment](super)).
fn sin_cos_reduced(t: &Float, wp: u64) -> (Ball, Ball) {
    if t.is_zero() {
        return (Ball::ZERO, one());
    }
    if wp >= BIT_BURST_MIN_PREC {
        return bitburst::sin_cos(t, wp);
    }
    if wp <= TABLE_MAX_PREC && exponent(t) >= -15 {
        if let Some(v) = tables::sin_cos(t, wp) {
            return v;
        }
    }
    sin_cos_small(t, wp)
}

/// Balls containing `sin(m)` and `cos(m)` for an exact `m` of exponent at
/// most `HUGE_EXP`, at about `wp` bits, relative to each result.
fn sin_cos_exact(m: &Float, wp: u64) -> (Ball, Ball) {
    if m.is_zero() {
        return (Ball::ZERO, one());
    }
    let e = exponent(m);
    if e <= -1 {
        return sin_cos_reduced(m, wp);
    }
    // k = round(m / (π/2)) < 2^(e + 1), from π at 64 bits beyond k's.
    let kbits = e.unsigned_abs() + 2;
    let half_pi = |p: u64| Ball::const_pi(p).mul_2exp(&-Integer::ONE);
    let k = to_integer(
        &m.div_round(half_pi(kbits + 64).mid(), kbits + 8, Round::Nearest)
            .0,
        false,
    );
    let wpi = i64::try_from(wp).unwrap_or(i64::MAX);
    // t = m - k π/2, with more bits of π while it cancels below wp bits.
    // A dyadic m stays away from multiples of π/2 by π's irrationality
    // measure [Mah53], so the extra stays within a few times m's length.
    let cap = 8 * (wp + m.bit_length() + kbits) + 256;
    let mut extra = 16;
    let t = loop {
        let p = wp + kbits + extra;
        let t = Ball::exact(m.clone()).sub(&half_pi(p).mul(&Ball::from_integer(&k), p), p);
        let known = t.rel_accuracy_bits();
        if known >= wpi + 8 || extra >= cap {
            break t;
        }
        extra += u64::try_from(wpi + 16 - known.max(0))
            .unwrap_or(cap)
            .max(extra);
    };
    let (s, c) = sin_cos_reduced(t.mid(), wp);
    let (s, c) = (s.add_error(t.rad()), c.add_error(t.rad()));
    // (sin, cos)(t + k π/2): rotate by the quadrant.
    let quadrant = k
        .div_rem_euclid(&Integer::from(4i64))
        .1
        .to_i64()
        .expect("0..4");
    match quadrant {
        0 => (s, c),
        1 => (c, s.neg()),
        2 => (s.neg(), c.neg()),
        _ => (c.neg(), s),
    }
}

/// `x mod 2` in `[-1, 1)` for a dyadic `x`, exactly. The work is bounded by
/// the mantissa: below 1 in magnitude `x` is already reduced (the shift
/// below would be as long as its exponent), and from 1 up its fraction has
/// fewer bits than the mantissa.
fn mod2(x: &Float) -> Float {
    if x.is_zero() || x.exponent().is_negative() || x.exponent().is_zero() {
        return x.clone();
    }
    let (m, e) = x.to_integer_exp();
    if !e.is_negative() {
        // An integer: m 2^e mod 2.
        if e.is_zero() && m.unsigned_abs().trailing_zeros() == Some(0) {
            return Float::from(1i64).neg();
        }
        return Float::ZERO;
    }
    // x = m / 2^s with s >= 1: reduce m modulo 2^(s+1) into [-2^s, 2^s).
    let s = (-e)
        .to_i64()
        .and_then(|s| u64::try_from(s).ok())
        .expect("bounded");
    let period = Integer::from_parts(false, Natural::from(1u64) << (s + 1));
    let half = Integer::from_parts(false, Natural::from(1u64) << s);
    let r = (&m + &half).div_rem_euclid(&period).1 - half;
    Float::from_integer_exp(&r, &-Integer::from(s))
}

/// `sin(π b)`, or `cos(π b)` when `cos`, for an exact `b` in `[0, 1/2]`,
/// at `prec` bits: exact at 0 and 1/2, else `sin` or `cos` of `π b'` with
/// `b' <= 1/4`, `b' = b` or `1/2 - b` (exact, and no longer than `b` since
/// `b > 1/4`), as `sin(π b) = cos(π (1/2 - b))`.
fn sin_pi_half(b: &Float, prec: u64, cos: bool) -> Ball {
    let quarter = Float::from(1i64).mul_2exp(&Integer::from(-2i64));
    let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
    if b.is_zero() {
        return if cos { one() } else { Ball::ZERO };
    }
    if *b == half {
        return if cos { Ball::ZERO } else { one() };
    }
    let wp = working(prec);
    let at = |c: &Float| Ball::const_pi(wp + 8).mul(&Ball::exact(c.clone()), wp + 8);
    let (t, take_cos) = if *b <= quarter {
        (at(b), cos)
    } else {
        (at(&half.sub(b)), !cos)
    };
    let (s, c) = sin_cos_reduced(t.mid(), wp);
    let v = if take_cos { c } else { s };
    v.add_error(t.rad()).round_to_accuracy(prec)
}

impl Ball {
    /// `(sin x, cos x)` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn sin_cos(&self, prec: u64) -> (Ball, Ball) {
        assert!(prec >= 1, "Ball::sin_cos: precision 0");
        let two = Mag::ONE.mul_2exp(&Integer::ONE);
        if !self.is_finite()
            || *self.rad() >= two
            || (!self.mid().is_zero() && *self.mid().exponent() > Integer::from(HUGE_EXP))
        {
            return (unit_ball(), unit_ball());
        }
        let (s, c) = sin_cos_exact(self.mid(), working(prec));
        (
            s.add_error(self.rad()).round_to_accuracy(prec),
            c.add_error(self.rad()).round_to_accuracy(prec),
        )
    }

    /// `sin x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn sin(&self, prec: u64) -> Ball {
        self.sin_cos(prec).0
    }

    /// `cos x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn cos(&self, prec: u64) -> Ball {
        self.sin_cos(prec).1
    }

    /// `tan x` for every `x` in the ball, at `prec` bits; indeterminate
    /// when the ball reaches a pole (its cosine contains 0).
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn tan(&self, prec: u64) -> Ball {
        let (s, c) = self.sin_cos(prec + 8);
        s.div(&c, prec)
    }

    /// `sin(π x)` for every `x` in the ball, at `prec` bits: exact at
    /// integers and half-integers.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn sin_pi(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::sin_pi: precision 0");
        if !self.is_finite() {
            return unit_ball();
        }
        // x mod 2 in [-1, 1); sin(π x) = sign · sin(π |x|), and sin(π a) =
        // sin(π (1 - a)) folds a in [1/2, 1] onto [0, 1/2].
        let r = mod2(self.mid());
        let neg = r.is_negative();
        let mut a = r.abs();
        let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
        if a > half {
            a = Float::from(1i64).sub(&a);
        }
        let v = sin_pi_half(&a, prec, false);
        let v = if neg { v.neg() } else { v };
        spread_pi(v, self.rad(), prec)
    }

    /// `cos(π x)` for every `x` in the ball, at `prec` bits: exact at
    /// integers and half-integers.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn cos_pi(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::cos_pi: precision 0");
        if !self.is_finite() {
            return unit_ball();
        }
        // x mod 2 in [-1, 1); cos(π x) = cos(π |x|), and cos(π a) = -cos(π (1
        // - a)) folds a in (1/2, 1] onto [0, 1/2), the difference exact.
        let mut a = mod2(self.mid()).abs();
        let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
        let neg = a > half;
        if neg {
            a = Float::from(1i64).sub(&a);
        }
        let v = sin_pi_half(&a, prec, true);
        let v = if neg { v.neg() } else { v };
        spread_pi(v, self.rad(), prec)
    }
}

/// `v` widened by `min(π ρ, 2)`; `[0 ± 1]` from `ρ >= 2/π`.
fn spread_pi(v: Ball, rho: &Mag, prec: u64) -> Ball {
    if rho.is_zero() {
        return v;
    }
    let spread = super::upper(&Ball::const_pi(64)).mul(rho);
    if spread.cmp(&Mag::ONE.mul_2exp(&Integer::ONE)) != Ordering::Less {
        return unit_ball();
    }
    v.add_error(&spread).round_to_accuracy(prec)
}
