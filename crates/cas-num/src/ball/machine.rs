//! `MachineBall`: an `f64` midpoint and radius, the fast path for
//! precisions up to 53 bits (CAS-DOC-5 §5.4).
//!
//! IEEE `+ − × ÷ √` round to nearest, so each midpoint is within half an ulp
//! of the exact result; that error, exactly when an error-free
//! transformation gives it (two-sum for sums, two-product residuals for the
//! others), goes into the radius. The radius itself is computed in
//! round-to-nearest with every step pushed up by one ulp ([`up`]): for a
//! correctly rounded `fl(x)`, `x <= next_up(fl(x))`, whatever the exponent.
//! Lower bounds (denominators) step down the same way.
//!
//! The exponent range is never trusted: an operation whose midpoint comes
//! out infinite or subnormal (or underflowed to zero), or whose radius
//! overflows, returns `None`, and the caller redoes it with [`Ball`] at 53
//! bits or more. A radius below the normal range is rounded up to
//! `f64::MIN_POSITIVE` instead, which is still an upper bound. So does an operation that a `Ball` would answer with the
//! indeterminate ball (division by a ball containing 0, the square root of
//! a ball reaching below 0), which a finite `f64` radius cannot hold.
//!
//! Every `f64` is dyadic, so [`MachineBall::to_ball`] is exact in the
//! midpoint (the radius becomes a [`Mag`], rounded up).
#![allow(
    clippy::many_single_char_names,
    clippy::float_cmp,
    reason = "midpoints a, b, radii r, s and the p, q, e of the error analysis, \
              whose comparisons are exact by design"
)]

use std::cmp::Ordering;

use super::{Ball, Float, Mag};

#[cfg(test)]
mod tests;

/// A real number enclosed as `mid ± rad` in `f64`s: both finite, normal or
/// zero, `rad >= 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MachineBall {
    mid: f64,
    rad: f64,
}

/// Results below this keep their error-free residuals exact: the residual of
/// a product or quotient is a multiple of about `|x| 2^-106`, representable
/// down to `2^-1074`.
const EFT_MIN: f64 = f64::from_bits((1023 - 900) << 52);

/// Factors from this up could overflow in Veltkamp's split (`2^27 + 1` times
/// the factor).
const SPLIT_MAX: f64 = f64::from_bits((1023 + 995) << 52);

/// `a b - p` for `p = fl(a b)`, exactly, by a fused multiply-add where the
/// target has one and Dekker's two-product otherwise (Rust's `mul_add` is a
/// software routine without the `fma` target feature); `None` where
/// exactness is not guaranteed (`|p| < EFT_MIN`, or a factor from
/// `SPLIT_MAX` up).
#[inline]
fn product_error(a: f64, b: f64, p: f64) -> Option<f64> {
    if !(p.abs() >= EFT_MIN && a.abs() < SPLIT_MAX && b.abs() < SPLIT_MAX) {
        return None;
    }
    if cfg!(target_feature = "fma") {
        return Some(a.mul_add(b, -p));
    }
    // Veltkamp's split: x = h + l, each half of at most 26 bits, so the four
    // partial products are exact.
    let split = |x: f64| {
        let c = 134_217_729.0 * x;
        let h = c - (c - x);
        (h, x - h)
    };
    let ((ah, al), (bh, bl)) = (split(a), split(b));
    Some(((ah * bh - p) + ah * bl + al * bh) + al * bl)
}

/// Whether `a = x y` exactly, for `p = fl(x y)` within a factor 2 of `a`
/// (so `a - p` is exact by Sterbenz's lemma); `false` when undecided.
#[inline]
fn is_exact_product(a: f64, x: f64, y: f64) -> bool {
    let p = x * y;
    a.abs() >= EFT_MIN && product_error(x, y, p).is_some_and(|e| a - p == e)
}

/// Normal or zero (and finite).
#[inline]
fn is_regular(x: f64) -> bool {
    x == 0.0 || x.is_normal()
}

/// The next `f64` above a finite `x` (`f64::next_up`, stable only from Rust
/// 1.86).
#[inline]
fn next_up(x: f64) -> f64 {
    let bits = x.to_bits();
    if x == 0.0 {
        f64::from_bits(1)
    } else if bits >> 63 == 0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

/// The next `f64` below a finite `x`.
#[inline]
fn next_down(x: f64) -> f64 {
    -next_up(-x)
}

/// An upper bound of the exact `x` from its correctly rounded `fl(x) >= 0`;
/// exact zeros (from zero operands) stay zero.
#[inline]
fn up(fl: f64, exact_zero: bool) -> f64 {
    if exact_zero { 0.0 } else { next_up(fl) }
}

/// `a + b` from above, for `a, b >= 0`.
#[inline]
fn add_up(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        next_up(a + b)
    }
}

/// `a × b` from above, for `a, b >= 0`.
#[inline]
fn mul_up(a: f64, b: f64) -> f64 {
    up(a * b, a == 0.0 || b == 0.0)
}

/// `a / b` from above, for `a >= 0`, `b > 0`.
#[inline]
fn div_up(a: f64, b: f64) -> f64 {
    up(a / b, a == 0.0)
}

/// A lower bound `>= 0` of the exact value whose correctly rounded value is
/// `fl`, for exact values `>= 0`; exact results (`exact`) are kept.
#[inline]
fn down(fl: f64, exact: bool) -> f64 {
    if exact || fl <= 0.0 {
        fl.max(0.0)
    } else {
        next_down(fl).max(0.0)
    }
}

/// Half an ulp of a normal `x`, `2^(e - 53)` for `|x|` in `[2^e, 2^(e+1))`,
/// as an `f64` (the smallest subnormal when that is below it).
#[inline]
fn half_ulp(x: f64) -> f64 {
    let biased = (x.to_bits() >> 52) & 0x7ff;
    if biased > 53 {
        f64::from_bits((biased - 53) << 52)
    } else {
        // 2^(biased - 1076) on the subnormal grid 2^-1074.
        f64::from_bits(1 << biased.saturating_sub(2))
    }
}

/// `x` exactly as a `Float`.
fn exact(x: f64) -> Float {
    Float::from_f64(x).expect("finite")
}

impl MachineBall {
    /// The exact zero ball.
    pub const ZERO: MachineBall = MachineBall { mid: 0.0, rad: 0.0 };

    /// `[mid ± rad]`; `None` unless both are finite and normal or zero and
    /// `rad >= 0`.
    #[must_use]
    #[inline]
    pub fn new(mid: f64, rad: f64) -> Option<MachineBall> {
        // `rad >= 0.0` also rules out NaN.
        (is_regular(mid) && is_regular(rad) && rad >= 0.0).then_some(MachineBall { mid, rad })
    }

    /// The exact ball `[x ± 0]`; `None` for non-finite or subnormal `x`.
    #[must_use]
    #[inline]
    pub fn exact(x: f64) -> Option<MachineBall> {
        MachineBall::new(x, 0.0)
    }

    #[must_use]
    #[inline]
    pub fn mid(&self) -> f64 {
        self.mid
    }

    #[must_use]
    #[inline]
    pub fn rad(&self) -> f64 {
        self.rad
    }

    /// Whether the radius is 0.
    #[must_use]
    #[inline]
    pub fn is_exact(&self) -> bool {
        self.rad == 0.0
    }

    /// The result ball from a midpoint, its error bound and the propagated
    /// radius, both upper bounds; `None` outside the regular range, except
    /// that a subnormal radius is rounded up to the smallest normal.
    #[inline]
    fn finish(mid: f64, err: f64, rad: f64) -> Option<MachineBall> {
        MachineBall::new(mid, normal_up(add_up(rad, err)))
    }

    /// `self + rhs`: radius `r + s + |e|`, `e` the exact rounding error of
    /// the midpoint (two-sum, exact for finite operands).
    #[must_use]
    #[inline]
    pub fn add(&self, rhs: &MachineBall) -> Option<MachineBall> {
        let (a, b) = (self.mid, rhs.mid);
        let (s, e) = two_sum(a, b);
        if !s.is_finite() {
            return None;
        }
        MachineBall::finish(s, e.abs(), add_up(self.rad, rhs.rad))
    }

    /// `self - rhs`, as [`MachineBall::add`].
    #[must_use]
    #[inline]
    pub fn sub(&self, rhs: &MachineBall) -> Option<MachineBall> {
        self.add(&rhs.neg())
    }

    /// `self × rhs`: radius `|a| s + |b| r + r s + |e|`, `e` the residual
    /// `a b - fl(a b)` (exact above `2^-900`, else bounded by half an
    /// ulp).
    #[must_use]
    #[inline]
    pub fn mul(&self, rhs: &MachineBall) -> Option<MachineBall> {
        let (a, b, r, s) = (self.mid, rhs.mid, self.rad, rhs.rad);
        let p = a * b;
        if p == 0.0 && a != 0.0 && b != 0.0 {
            return None;
        }
        let err = match product_error(a, b, p) {
            Some(e) => e.abs(),
            None if p == 0.0 => 0.0,
            None => half_ulp(p),
        };
        let rad = add_up(add_up(mul_up(a.abs(), s), mul_up(b.abs(), r)), mul_up(r, s));
        MachineBall::finish(p, err, rad)
    }

    /// `self / rhs`: radius `(|a| s + |b| r) / (|b| (|b| - s))` with lower
    /// bounds in the denominator, plus half an ulp unless the residual
    /// `a - q b` is 0; `None` if `rhs` may contain 0.
    #[must_use]
    #[inline]
    pub fn div(&self, rhs: &MachineBall) -> Option<MachineBall> {
        let (a, b, r, s) = (self.mid, rhs.mid, self.rad, rhs.rad);
        let bl = b.abs();
        // |b| - s from below; |b| - s is exact when s = 0.
        let gap = down(bl - s, s == 0.0);
        if gap == 0.0 {
            return None;
        }
        let q = a / b;
        if q == 0.0 && a != 0.0 {
            return None;
        }
        let err = if q == 0.0 || (q.is_normal() && is_exact_product(a, q, b)) {
            0.0
        } else {
            half_ulp(q)
        };
        let rad = if r == 0.0 && s == 0.0 {
            0.0
        } else {
            let num = add_up(mul_up(a.abs(), s), mul_up(bl, r));
            let den = down(bl * gap, false);
            if den == 0.0 {
                return None;
            }
            div_up(num, den)
        };
        MachineBall::finish(q, err, rad)
    }

    /// `1 / self`; `None` if `self` may contain 0.
    #[must_use]
    #[inline]
    pub fn inv(&self) -> Option<MachineBall> {
        MachineBall { mid: 1.0, rad: 0.0 }.div(self)
    }

    /// `sqrt(self)`: radius `r / (sqrt(a - r) + sqrt(a))` with lower bounds,
    /// plus half an ulp unless `a - s^2` is 0; `None` if the ball reaches
    /// below 0.
    #[must_use]
    #[inline]
    pub fn sqrt(&self) -> Option<MachineBall> {
        let (a, r) = (self.mid, self.rad);
        if a < r {
            return None;
        }
        let m = a.sqrt();
        let err = if a == 0.0 || is_exact_product(a, m, m) {
            0.0
        } else {
            half_ulp(m)
        };
        let rad = if r == 0.0 {
            0.0
        } else {
            let lo = down(a - r, false);
            let den = add_down(down(lo.sqrt(), lo == 0.0), down(m, false));
            div_up(r, den)
        };
        MachineBall::finish(m, err, rad)
    }

    /// `self^n` (`self^0 = 1`) by squaring and multiplying; a negative `n`
    /// inverts the power.
    #[must_use]
    pub fn pow(&self, n: i32) -> Option<MachineBall> {
        let k = n.unsigned_abs();
        let mut acc = MachineBall { mid: 1.0, rad: 0.0 };
        for i in (0..u32::BITS - k.leading_zeros()).rev() {
            acc = acc.mul(&acc)?;
            if (k >> i) & 1 == 1 {
                acc = acc.mul(self)?;
            }
        }
        if n < 0 { acc.inv() } else { Some(acc) }
    }

    /// `-self`, exactly.
    #[must_use]
    #[inline]
    pub fn neg(&self) -> MachineBall {
        MachineBall {
            mid: -self.mid,
            rad: self.rad,
        }
    }

    /// `[|m| ± r]`, which contains `|x|` for every `x` in the ball.
    #[must_use]
    #[inline]
    pub fn abs(&self) -> MachineBall {
        MachineBall {
            mid: self.mid.abs(),
            rad: self.rad,
        }
    }

    /// Whether `x` lies in the ball, exactly: two-sum gives `x - m = d + e`
    /// with `|e| <= ulp(d) / 2`, and as `r` is a double, `|d| < r` and
    /// `|d| > r` decide alone; `|d| = r` leaves the sign of `e`.
    #[must_use]
    #[inline]
    pub fn contains(&self, x: f64) -> bool {
        if !x.is_finite() {
            return false;
        }
        let (d, e) = two_sum(x, -self.mid);
        if !d.is_finite() {
            return false;
        }
        match d.abs().partial_cmp(&self.rad) {
            Some(Ordering::Less) => true,
            Some(Ordering::Equal) => e == 0.0 || (e < 0.0) != (d < 0.0),
            _ => false,
        }
    }

    /// Whether the balls share a point: `|a - b| <= r + s`, decided by
    /// rounded bounds and exactly in `Float`s when they cannot.
    #[must_use]
    pub fn overlaps(&self, other: &MachineBall) -> bool {
        let dist = (self.mid - other.mid).abs();
        let sum = self.rad + other.rad;
        if dist.is_finite() && sum.is_finite() {
            if next_up(dist) <= next_down(sum) {
                return true;
            }
            if next_down(dist) > next_up(sum) {
                return false;
            }
        }
        let dist = exact(self.mid).sub(&exact(other.mid)).abs();
        dist <= exact(self.rad).add(&exact(other.rad))
    }

    /// Whether `other` lies entirely in `self`: `|b - a| + s <= r`, as
    /// [`MachineBall::overlaps`].
    #[must_use]
    pub fn contains_ball(&self, other: &MachineBall) -> bool {
        if other.rad > self.rad {
            return false;
        }
        let dist = (self.mid - other.mid).abs();
        let room = self.rad - other.rad;
        if dist.is_finite() {
            if next_up(dist) <= down(room, other.rad == 0.0) {
                return true;
            }
            if next_down(dist) > next_up(room) {
                return false;
            }
        }
        let dist = exact(self.mid).sub(&exact(other.mid)).abs();
        dist <= exact(self.rad).sub(&exact(other.rad))
    }

    /// The sign of every point of the ball, if they share one: `Some(Equal)`
    /// only for the exact zero ball.
    #[must_use]
    #[inline]
    pub fn sign(&self) -> Option<Ordering> {
        if self.mid == 0.0 {
            return (self.rad == 0.0).then_some(Ordering::Equal);
        }
        let side = if self.mid > 0.0 {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        (self.mid.abs() > self.rad).then_some(side)
    }

    /// How every point of `self` compares with every point of `other`, if
    /// decided: the balls are disjoint, or both are the same exact value.
    #[must_use]
    #[allow(clippy::should_implement_trait, reason = "three-valued, as Ball::cmp")]
    pub fn cmp(&self, other: &MachineBall) -> Option<Ordering> {
        if self.is_exact() && other.is_exact() {
            return self.mid.partial_cmp(&other.mid);
        }
        (!self.overlaps(other)).then(|| self.mid.total_cmp(&other.mid))
    }

    /// The same ball as a [`Ball`]: the midpoint exactly, the radius rounded
    /// up to a [`Mag`].
    #[must_use]
    pub fn to_ball(&self) -> Ball {
        Ball::new(exact(self.mid), Mag::from_float(&exact(self.rad)))
    }

    /// A machine ball containing `b`: the midpoint rounded to the nearest
    /// `f64`, the radius grown by the rounding error; `None` when either is
    /// out of the regular range, or `b` is indeterminate.
    #[must_use]
    pub fn from_ball(b: &Ball) -> Option<MachineBall> {
        if !b.is_finite() {
            return None;
        }
        let mid = b.mid().to_f64();
        if !is_regular(mid) || (mid == 0.0 && !b.mid().is_zero()) {
            return None;
        }
        // rad + |b.mid - mid|, both exact Floats; the sum as a Mag has a
        // 30-bit mantissa, so its f64 is exact unless out of range.
        let err = Mag::from_float(&b.mid().sub(&exact(mid)));
        let total = b.rad().add(&err);
        let rad = total.to_float()?.to_f64();
        // Below the normal range the nearest f64 may be below the Mag.
        let rad = if total.is_zero() {
            0.0
        } else {
            rad.max(f64::MIN_POSITIVE)
        };
        MachineBall::new(mid, rad)
    }

    /// `f` evaluated on the ball through [`Ball`] at 64 bits, the fallback
    /// for functions without a machine kernel yet (CAS-88).
    #[must_use]
    pub fn via_ball(&self, f: impl FnOnce(&Ball, u64) -> Ball) -> Option<MachineBall> {
        MachineBall::from_ball(&f(&self.to_ball(), 64))
    }
}

/// `x >= 0` raised to `f64::MIN_POSITIVE` when it is positive and below it.
#[inline]
fn normal_up(x: f64) -> f64 {
    if x > 0.0 && x < f64::MIN_POSITIVE {
        f64::MIN_POSITIVE
    } else {
        x
    }
}

/// `a + b` from below, for `a, b >= 0`.
#[inline]
fn add_down(a: f64, b: f64) -> f64 {
    down(a + b, a == 0.0 || b == 0.0)
}

/// Knuth's two-sum: `(s, e)` with `s = fl(a + b)` and `s + e = a + b`
/// exactly when `s` is finite.
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}
