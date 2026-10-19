//! `atan`, `asin`, `acos`, `atan2` (CAS-DOC-5 §7.2).
//!
//! `atan` at the exact midpoint: `|m| > 1` folds to `±π/2 - atan(1/m)`;
//! then `h` halvings `x ← x / (1 + √(1 + x^2))` (`atan x = 2 atan` of the
//! result; each step at least halves `|x|`) as balls, to `|x| <= 2^-r` with
//! `r ≈ p^(1/3) / 3` (a halving costs a root and a quotient, more than a
//! series term), and the series `x Σ (-x^2)^j / (2j+1)` in fixed point
//! ([`atan_sum`]). The input radius `ρ` adds `ρ / (1 + (|m| - ρ)^2)` (`ρ`
//! when the ball reaches 0).
//!
//! `asin x = atan(x / √((1 - x)(1 + x)))`, `acos x = 2 atan(√((1 - x) / (1 +
//! x)))` for `x >= 0` and `π - acos(-x)` below (no cancellation near ±1), as
//! ball compositions, whose radii propagate through every operation; exact
//! `±1` are special-cased. `atan2(y, x)` decides its quadrant from the signs
//! of the balls; a ball straddling the branch cut (`x < 0`, `y ∋ 0`), or
//! both balls containing 0, gives `[0 ± π]`, which holds every value.

use std::cmp::Ordering;

use super::exp::to_fixed;
use super::fixed::{atan_sum, shr};
use super::{BIT_BURST_MIN_PREC, TABLE_MAX_PREC, bitburst, exponent, one, tables, upper, working};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Natural};

/// `π/2` at `p` bits.
fn half_pi(p: u64) -> Ball {
    Ball::const_pi(p).mul_2exp(&-Integer::ONE)
}

/// `[0 ± π]`, the hull of every `atan2`.
fn around_pi() -> Ball {
    Ball::new(Float::ZERO, upper(&Ball::const_pi(64)))
}

/// A ball containing `atan(x)` for an exact `|x| <= 1/2`, at about `wp`
/// bits relative to the result.
///
/// With `|x| < 2^ex`, `y = x^2 < 2^(2 ex)`. In fixed point at `w` bits `X`
/// is within `2^-w` of `x` (moving `atan` by at most as much) and `Y =
/// trunc(X^2 / 2^w)` within a unit of `X'^2`, which moves the series (its
/// derivative in `y` is below 1) by at most a unit. The tail beyond `n` terms
/// is at most `|y|^n / (1 - |y|) <= 2^(2 ex n + 1)`. The sum's error count
/// grows like `n √n / 2`, so `w` has that many more bits, and `|ex|` more to
/// stay relative to `atan x ≈ x`. Below `2^-wp`, `atan x ∈ x ± |x|^3`.
pub(super) fn atan_series(x: &Float, wp: u64) -> Ball {
    odd_series(x, wp, false)
}

/// A ball containing `atanh(x)` for an exact `|x| <= 1/2`, at about `wp`
/// bits relative to the result: [`atan_series`]'s sum at `y = -x^2`, with the
/// same bounds (they only use `|y|`) except the shift from `x` to its fixed
/// point value, which moves `atanh` by up to `4/3` of a unit; below `2^-wp`, `atanh x ∈ x ± |x|^3`
/// (`atanh x - x = x^3/3 + x^5/5 + … <= |x|^3 / (3 (1 - x^2))`).
pub(super) fn atanh_series(x: &Float, wp: u64) -> Ball {
    odd_series(x, wp, true)
}

/// [`atan_series`], or [`atanh_series`] when `hyperbolic`.
fn odd_series(x: &Float, wp: u64, hyperbolic: bool) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    let ex = exponent(x);
    debug_assert!(ex <= 0);
    if ex < -i64::try_from(wp).unwrap_or(i64::MAX) {
        let x2 = Mag::from_float(x).mul(&Mag::from_float(x));
        return Ball::new(x.clone(), x2.mul(&Mag::from_float(x)));
    }
    let per = (-2 * ex).unsigned_abs().max(1);
    let w0 = wp + 24 + ex.unsigned_abs();
    let n = (w0 + 3).div_ceil(per).max(1);
    let w = w0 + u64::from(u64::BITS - (n * (n.isqrt() + 4)).leading_zeros());
    let xf = to_fixed(x, w);
    let yf = shr(&xf * &xf, w);
    let yf = if hyperbolic { -yf } else { yf };
    let (s, e) = atan_sum(&yf, w, n);
    let unit = |k: Natural| Mag::from_integer(&Integer::from(k)).mul_2exp(&-Integer::from(w));
    let tail = Mag::power_of_two(&(Integer::from(2 * ex) * Integer::from(n) + Integer::ONE));
    let at = |v: Integer| Float::from_integer_exp(&v, &-Integer::from(w));
    let ratio = Ball::new(at(s), unit(&e + &Natural::from(1u64)).add(&tail));
    // x' = X / 2^w is within a unit of x: atan' <= 1, atanh' <= 4/3 there.
    Ball::exact(at(xf))
        .mul(&ratio, wp + 8)
        .add_error(&unit(Natural::from(if hyperbolic { 2u64 } else { 1 })))
}

/// A ball containing `atan(x)` for every `x` in `b`, a ball with `|mid| <=
/// 1`: from the core for `wp` (see the [module comment](super)) at the
/// midpoint, then the radius (`atan` is 1-Lipschitz), or [`atan_halving`].
pub(super) fn atan_small(b: &Ball, wp: u64) -> Ball {
    let m = b.mid();
    if m.is_zero() {
        return Ball::new(Float::ZERO, b.rad().clone());
    }
    let e = exponent(m);
    let fast = |f: fn(&Float, u64) -> Ball| {
        let v = f(&m.abs(), wp);
        let v = if m.is_negative() { v.neg() } else { v };
        v.add_error(b.rad())
    };
    let within_one = m.cmp_abs(&Float::from(1i64)) != Ordering::Greater;
    if within_one && wp >= BIT_BURST_MIN_PREC {
        return fast(bitburst::atan);
    }
    if within_one && wp <= TABLE_MAX_PREC && e >= -15 {
        return fast(tables::atan);
    }
    atan_halving(b, wp)
}

/// [`atan_small`] by the baseline: halvings of the whole ball, then
/// [`atan_series`] (the tables are built by this).
pub(super) fn atan_halving(b: &Ball, wp: u64) -> Ball {
    let m = b.mid();
    if m.is_zero() {
        return Ball::new(Float::ZERO, b.rad().clone());
    }
    let e = exponent(m);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        clippy::cast_possible_wrap,
        reason = "a cube root of a precision, a small positive count"
    )]
    let r = (((wp as f64).cbrt() / 3.0).ceil() as i64).max(2);
    let h = (r + e + 1).max(0);
    let q = wp + 16 + u64::from(u64::BITS - h.unsigned_abs().leading_zeros());
    let mut x = b.clone();
    for _ in 0..h {
        let root = one().add(&x.mul(&x, q), q).sqrt(q);
        x = x.div(&one().add(&root, q), q);
    }
    if !x.is_finite() {
        return Ball::INDETERMINATE;
    }
    // atan is 1-Lipschitz: the halved ball's radius carries over.
    atan_series(x.mid(), wp)
        .add_error(x.rad())
        .mul_2exp(&Integer::from(h))
}

/// A ball containing `atan(m)` for an exact `m`, at about `wp` bits.
fn atan_exact(m: &Float, wp: u64) -> Ball {
    if m.is_zero() {
        return Ball::ZERO;
    }
    let one_f = Float::from(1i64);
    if m.cmp_abs(&one_f) != Ordering::Greater {
        return atan_small(&Ball::exact(m.clone()), wp);
    }
    // atan m = ±π/2 - atan(1/m): no cancellation, the difference is >= π/4.
    let inv = one().div(&Ball::exact(m.clone()), wp + 8);
    let r = atan_small(&inv, wp);
    let hp = half_pi(wp + 8);
    let hp = if m.is_negative() { hp.neg() } else { hp };
    hp.sub(&r, wp + 8)
}

impl Ball {
    /// `atan x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn atan(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::atan: precision 0");
        if !self.is_finite() {
            return Ball::new(Float::ZERO, upper(&half_pi(64)));
        }
        let core = atan_exact(self.mid(), working(prec));
        // |atan'| = 1 / (1 + x^2) <= 1 / (1 + (|m| - ρ)^2) on the ball.
        let rho = self.rad();
        let lo = Mag::from_float_lower(self.mid()).sub_lower(rho);
        let den = Mag::ONE.add_lower(&lo.mul_lower(&lo));
        core.add_error(&rho.div(&den)).round_to_accuracy(prec)
    }

    /// `asin x` for every `x` in the ball, at `prec` bits; indeterminate
    /// unless the ball lies in `[-1, 1]` and, if inexact, inside `(-1, 1)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn asin(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::asin: precision 0");
        if self.is_exact() && self.mid().cmp_abs(&Float::from(1i64)) == Ordering::Equal {
            let hp = half_pi(prec + 8).round_to_accuracy(prec);
            return if self.mid().is_negative() {
                hp.neg()
            } else {
                hp
            };
        }
        let wp = working(prec);
        let d = one().sub(self, wp).mul(&one().add(self, wp), wp);
        let q = self.div(&d.sqrt(wp), wp);
        // Outside the domain (or touching ±1 with a radius): indeterminate.
        if q.is_finite() {
            q.atan(prec)
        } else {
            Ball::INDETERMINATE
        }
    }

    /// `acos x` for every `x` in the ball, at `prec` bits; indeterminate
    /// unless the ball lies in `[-1, 1]` and, if inexact, inside `(-1, 1)`
    /// near `-1`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn acos(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::acos: precision 0");
        let wp = working(prec);
        if self.is_exact() && self.mid().cmp_abs(&Float::from(1i64)) == Ordering::Equal {
            return if self.mid().is_negative() {
                Ball::const_pi(prec)
            } else {
                Ball::ZERO
            };
        }
        // 2 atan(√((1 - x) / (1 + x))), from x >= 0; π - acos(-x) below.
        let half_angle = |x: &Ball| {
            let r = one().sub(x, wp).div(&one().add(x, wp), wp).sqrt(wp);
            // Outside the domain: indeterminate, not atan's hull.
            r.is_finite().then(|| r.atan(wp).mul_2exp(&Integer::ONE))
        };
        let v = if self.mid().is_negative() {
            half_angle(&self.neg()).map(|a| Ball::const_pi(wp).sub(&a, wp))
        } else {
            half_angle(self)
        };
        v.map_or(Ball::INDETERMINATE, |v| v.round_to_accuracy(prec))
    }

    /// `atan2(self, x)`, the angle of the point `(x, self)`, for every pair
    /// of points of the balls, at `prec` bits; `[0 ± π]` when the balls
    /// straddle the branch cut on the negative `x` axis or both contain 0.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn atan2(&self, x: &Ball, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::atan2: precision 0");
        let wp = working(prec);
        let y = self;
        let pi = || Ball::const_pi(wp);
        let v = match (x.sign(), y.sign()) {
            // Right half-plane: the principal atan.
            (Some(Ordering::Greater), _) => y.div(x, wp).atan(wp),
            // Left half-plane, off the cut (y = 0 exactly gives π).
            (Some(Ordering::Less), Some(Ordering::Greater | Ordering::Equal)) => {
                y.div(x, wp).atan(wp).add(&pi(), wp)
            }
            (Some(Ordering::Less), Some(Ordering::Less)) => y.div(x, wp).atan(wp).sub(&pi(), wp),
            // x may be 0: ±π/2 - atan(x / y) from the sign of y.
            (_, Some(Ordering::Greater)) => half_pi(wp).sub(&x.div(y, wp).atan(wp), wp),
            (_, Some(Ordering::Less)) => half_pi(wp).neg().sub(&x.div(y, wp).atan(wp), wp),
            _ => return around_pi(),
        };
        if v.is_finite() {
            v.round_to_accuracy(prec)
        } else {
            around_pi()
        }
    }
}
