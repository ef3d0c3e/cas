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

mod dd;
mod kernels;
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
/// exactness is not guaranteed (`p` infinite, `|p| < EFT_MIN`, or a factor
/// from `SPLIT_MAX` up).
#[inline]
fn product_error(a: f64, b: f64, p: f64) -> Option<f64> {
    // A product rounded to ±∞ has no residual: `a - p == e` would compare
    // infinities (MAX / 5.6e17 passed as exact that way).
    if !(p.is_finite() && p.abs() >= EFT_MIN && a.abs() < SPLIT_MAX && b.abs() < SPLIT_MAX) {
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

    /// `exp`: the double-double kernel within `2^-90` and half an ulp, the
    /// spread `e^m (e^r - 1) <= |v| (r + r^2)` for `r <= 2^-20`; `Ball` at 64
    /// bits beyond the kernel's range or for wider balls.
    #[must_use]
    pub fn exp(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::exp(self.mid)?;
            let spread = self.small_rad()?;
            let spread = mul_up(spread, add_up(1.0, spread));
            MachineBall::from_kernel(v, |size| mul_up(size, spread))
        };
        fast().or_else(|| self.via_ball(Ball::exp))
    }

    /// `log`: the kernel, spread `r / (m - r)`; `None` unless the ball lies
    /// in `x > 0`.
    #[must_use]
    pub fn log(&self) -> Option<MachineBall> {
        let gap = down(self.mid - self.rad, self.rad == 0.0);
        if !(self.mid > 0.0 && gap > 0.0) {
            return None;
        }
        let fast = || {
            let v = kernels::log(self.mid)?;
            let spread = div_up(self.rad, gap);
            MachineBall::from_kernel(v, |_| spread)
        };
        fast().or_else(|| self.via_ball(Ball::log))
    }

    /// `(sin, cos)`: the kernel, spread `r` (both 1-Lipschitz) for `r <=
    /// 2^-20`; `Ball` for `|m| > 2^20`, next to multiples of π/2, or wider
    /// balls.
    #[must_use]
    pub fn sin_cos(&self) -> Option<(MachineBall, MachineBall)> {
        let fast = || {
            let (s, c) = kernels::sin_cos(self.mid)?;
            let spread = self.small_rad()?;
            Some((
                MachineBall::from_kernel(s, |_| spread)?,
                MachineBall::from_kernel(c, |_| spread)?,
            ))
        };
        fast().or_else(|| {
            let b = self.to_ball();
            let (s, c) = b.sin_cos(64);
            Some((MachineBall::from_ball(&s)?, MachineBall::from_ball(&c)?))
        })
    }

    /// `sin`, as [`MachineBall::sin_cos`].
    #[must_use]
    pub fn sin(&self) -> Option<MachineBall> {
        self.sin_cos().map(|(s, _)| s)
    }

    /// `cos`, as [`MachineBall::sin_cos`].
    #[must_use]
    pub fn cos(&self) -> Option<MachineBall> {
        self.sin_cos().map(|(_, c)| c)
    }

    /// `tan = sin / cos`; `None` across a pole.
    #[must_use]
    pub fn tan(&self) -> Option<MachineBall> {
        let (s, c) = self.sin_cos()?;
        s.div(&c)
    }

    /// `atan`: the kernel, spread `r` (`|atan'| <= 1`).
    #[must_use]
    pub fn atan(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::atan(self.mid)?;
            let spread = self.rad;
            MachineBall::from_kernel(v, |_| spread)
        };
        fast().or_else(|| self.via_ball(Ball::atan))
    }

    /// The radius if it is at most `2^-20` (the kernels' spreads assume a
    /// narrow ball).
    fn small_rad(&self) -> Option<f64> {
        (self.rad <= f64::from_bits((1023 - 20) << 52)).then_some(self.rad)
    }

    /// The ball of a kernel value `v` (relative error at most
    /// `KERNEL_ERR`): the nearest double, its distance to `v`, `2 KERNEL_ERR
    /// |v|` (`|f - v| <= KERNEL_ERR |f|` gives `|f - v| <= KERNEL_ERR |v| /
    /// (1 - KERNEL_ERR)`), and the spread from `spread(|v| bound)`.
    fn from_kernel(v: dd::Dd, spread: impl FnOnce(f64) -> f64) -> Option<MachineBall> {
        MachineBall::from_kernel_err(v, kernels::KERNEL_ERR, spread)
    }

    /// [`MachineBall::from_kernel`] for a kernel within `err <= 1/2`
    /// relative.
    fn from_kernel_err(
        v: dd::Dd,
        err: f64,
        spread: impl FnOnce(f64) -> f64,
    ) -> Option<MachineBall> {
        let mid = v.to_f64();
        if !mid.is_finite() {
            return None;
        }
        let (round, size) = kernel_size(v, mid);
        let kernel = mul_up(size, 2.0 * err);
        let rad = add_up(add_up(round, kernel), spread(size));
        MachineBall::finish(mid, 0.0, rad)
    }

    /// `exp(x) - 1`, through [`Ball`] at 64 bits.
    #[must_use]
    pub fn expm1(&self) -> Option<MachineBall> {
        self.via_ball(Ball::expm1)
    }

    /// `log(1 + x)`, through [`Ball`] at 64 bits; `None` unless the ball
    /// lies in `x > -1`.
    #[must_use]
    pub fn log1p(&self) -> Option<MachineBall> {
        self.via_ball(Ball::log1p)
    }

    /// `asin`, through [`Ball`] at 64 bits until CAS-88's kernels; `None` outside `[-1, 1]`.
    #[must_use]
    pub fn asin(&self) -> Option<MachineBall> {
        self.via_ball(Ball::asin)
    }

    /// `acos`, through [`Ball`] at 64 bits until CAS-88's kernels; `None` outside `[-1, 1]`.
    #[must_use]
    pub fn acos(&self) -> Option<MachineBall> {
        self.via_ball(Ball::acos)
    }

    /// `sinh`: the kernel, spread `r cosh(|m| + r) <= r C (1 + 2r)` for `r
    /// <= 2^-20` (`C` the kernel's `cosh m` from above, `cosh(a + r) <=
    /// cosh(a) e^r`); `Ball` beyond `|m| = 708` or for wider balls.
    #[must_use]
    pub fn sinh(&self) -> Option<MachineBall> {
        let fast = || {
            let (s, c) = kernels::sinh_cosh(self.mid)?;
            let r = self.small_rad()?;
            let slope = mul_up(kernel_bounds(c).1, add_up(1.0, 2.0 * r));
            MachineBall::from_kernel(s, |_| mul_up(r, slope))
        };
        fast().or_else(|| self.via_ball(Ball::sinh))
    }

    /// `cosh`: the kernel, spread `r (S + r C (1 + 2r))` for `r <= 2^-20`:
    /// on the ball `|sinh| <= sinh |m| + r cosh(|m| + r)`, with the kernel's
    /// `S`, `C` from above.
    #[must_use]
    pub fn cosh(&self) -> Option<MachineBall> {
        let fast = || {
            let (s, c) = kernels::sinh_cosh(self.mid)?;
            let r = self.small_rad()?;
            let far = mul_up(r, mul_up(kernel_bounds(c).1, add_up(1.0, 2.0 * r)));
            let slope = add_up(kernel_bounds(s).1, far);
            MachineBall::from_kernel(c, |_| mul_up(r, slope))
        };
        fast().or_else(|| self.via_ball(Ball::cosh))
    }

    /// `tanh`: the kernel, spread `r (1 - t^2)` for `t = T - r` (or 0), `T`
    /// the kernel's `|tanh m|` from below: `tanh' = 1 - tanh^2` falls with
    /// `|x|`, and `|tanh| >= T - r` on the ball.
    #[must_use]
    pub fn tanh(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::tanh(self.mid)?;
            let r = self.rad;
            let t = down(kernel_bounds(v).0 - r, false);
            let slope = up(1.0 - mul_down(t, t), false).min(1.0);
            MachineBall::from_kernel(v, |_| mul_up(r, slope))
        };
        fast().or_else(|| self.via_ball(Ball::tanh))
    }

    /// `asinh`: the kernel, spread `r / max(1, |m| - r)`, as `asinh' = 1 /
    /// √(1 + x^2)`.
    #[must_use]
    pub fn asinh(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::asinh(self.mid)?;
            let r = self.rad;
            let den = down(self.mid.abs() - r, r == 0.0).max(1.0);
            MachineBall::from_kernel(v, |_| div_up(r, den))
        };
        fast().or_else(|| self.via_ball(Ball::asinh))
    }

    /// `acosh`: the kernel, spread `r / (√(m - r - 1) √(m - r + 1))` from
    /// lower bounds; `Ball` when the ball reaches 1 (where `acosh'` is
    /// unbounded), `None` unless it lies in `x >= 1`.
    #[must_use]
    pub fn acosh(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::acosh(self.mid)?;
            let r = self.rad;
            if r == 0.0 {
                return MachineBall::from_kernel(v, |_| 0.0);
            }
            let lo = down(self.mid - r, false);
            let gap = down(lo - 1.0, false);
            let den = mul_down(
                down(gap.sqrt(), false),
                down(add_down(lo, 1.0).sqrt(), false),
            );
            if den == 0.0 {
                return None;
            }
            MachineBall::from_kernel(v, |_| div_up(r, den))
        };
        fast().or_else(|| self.via_ball(Ball::acosh))
    }

    /// `atanh`: the kernel, spread `r / (1 - h)` for `h = |m| + r` from
    /// above (`atanh' = 1 / ((1 - x)(1 + x))`); `Ball` when `h` is not
    /// below 1, `None` unless the ball lies in `(-1, 1)`.
    #[must_use]
    pub fn atanh(&self) -> Option<MachineBall> {
        let fast = || {
            let v = kernels::atanh(self.mid)?;
            let r = self.rad;
            let gap = down(1.0 - add_up(self.mid.abs(), r), false);
            if gap == 0.0 {
                return None;
            }
            MachineBall::from_kernel(v, |_| div_up(r, gap))
        };
        fast().or_else(|| self.via_ball(Ball::atanh))
    }

    /// `(sin(π x), cos(π x))`: the kernel, exact at integers and
    /// half-integers; spreads `π r min(1, C + π r)` for the sine and `π r
    /// min(1, S + π r)` for the cosine (`|sin(π (m + δ)) - sin(π m)| <= π |δ|
    /// max |cos(π x)|`, and `|cos(π x)| <= |cos(π m)| + π r` on the ball),
    /// with the kernel's `C`, `S` from above.
    #[must_use]
    pub fn sin_cos_pi(&self) -> Option<(MachineBall, MachineBall)> {
        let fast = || {
            let (s, c, exact) = kernels::sin_cos_pi(self.mid)?;
            let pr = mul_up(PI_UP, self.rad);
            let spread = |other: dd::Dd| mul_up(pr, add_up(kernel_bounds(other).1, pr).min(1.0));
            let (ss, sc) = (spread(c), spread(s));
            if exact {
                // 0 and ±1, exactly.
                return Some((
                    MachineBall::finish(s.hi, 0.0, ss)?,
                    MachineBall::finish(c.hi, 0.0, sc)?,
                ));
            }
            Some((
                MachineBall::from_kernel(s, |_| ss)?,
                MachineBall::from_kernel(c, |_| sc)?,
            ))
        };
        fast().or_else(|| {
            let b = self.to_ball();
            Some((
                MachineBall::from_ball(&b.sin_pi(64))?,
                MachineBall::from_ball(&b.cos_pi(64))?,
            ))
        })
    }

    /// `sin(π x)`, as [`MachineBall::sin_cos_pi`].
    #[must_use]
    pub fn sin_pi(&self) -> Option<MachineBall> {
        self.sin_cos_pi().map(|(s, _)| s)
    }

    /// `cos(π x)`, as [`MachineBall::sin_cos_pi`].
    #[must_use]
    pub fn cos_pi(&self) -> Option<MachineBall> {
        self.sin_cos_pi().map(|(_, c)| c)
    }

    /// `atan2(self, x)`, the angle of the point `(x, self)` for every pair
    /// of points of the balls, as [`Ball::atan2`]: `atan(y / x)` (`± π` on
    /// the left) for `|y| < |x|`, `±π/2 - atan(x / y)` otherwise, `[0 ± π]`
    /// across the cut on the negative `x` axis or when both balls contain 0.
    #[must_use]
    pub fn atan2(&self, x: &MachineBall) -> Option<MachineBall> {
        let y = self;
        let fast = || {
            // π and π/2 within half an ulp.
            let pi = MachineBall {
                mid: std::f64::consts::PI,
                rad: half_ulp(std::f64::consts::PI),
            };
            let hp = MachineBall {
                mid: std::f64::consts::FRAC_PI_2,
                rad: half_ulp(std::f64::consts::FRAC_PI_2),
            };
            let from_y = |up: bool| {
                let h = if up { hp } else { hp.neg() };
                h.sub(&x.div(y)?.atan()?)
            };
            match (x.sign(), y.sign()) {
                (_, Some(side @ (Ordering::Greater | Ordering::Less)))
                    if y.mid.abs() >= x.mid.abs() =>
                {
                    from_y(side == Ordering::Greater)
                }
                (Some(Ordering::Greater), _) => y.div(x)?.atan(),
                (Some(Ordering::Less), Some(Ordering::Greater | Ordering::Equal)) => {
                    y.div(x)?.atan()?.add(&pi)
                }
                (Some(Ordering::Less), Some(Ordering::Less)) => y.div(x)?.atan()?.sub(&pi),
                (_, Some(Ordering::Greater)) => from_y(true),
                (_, Some(Ordering::Less)) => from_y(false),
                _ => MachineBall::new(0.0, PI_UP),
            }
        };
        fast().or_else(|| MachineBall::from_ball(&y.to_ball().atan2(&x.to_ball(), 64)))
    }

    /// `x^y` for every `x` in `self` and `y` in `y`, as [`Ball::pow_ball`]:
    /// an exact integer `y` is [`MachineBall::pow`], an exact half-integer
    /// `n + 1/2` is `sqrt(x)^(2n + 1)`, an exact 0 base gives 0 for `y > 0`;
    /// otherwise, for `x > 0`, the kernel `exp(y log x)` (within `(1 + |y log
    /// x|) 2^-93`) and the spread `|v| δ (1 + 2δ)` for `δ = |n| λ + ρ_y (|log
    /// m| + λ) <= 2^-20`, `λ = ρ_x / (m - ρ_x)` (on the balls `|log x - log
    /// m| <= λ`, so `|y log x - n log m| <= δ`). `Ball` at 64 bits for the
    /// rest; `None` where that is indeterminate.
    #[must_use]
    pub fn pow_ball(&self, y: &MachineBall) -> Option<MachineBall> {
        let fast = || {
            if let Some(n) = exact_i32(y.mid).filter(|_| y.is_exact()) {
                return self.pow(n);
            }
            if let Some(k) = exact_i32(2.0 * y.mid).filter(|_| y.is_exact()) {
                return self.sqrt()?.pow(k);
            }
            if self.mid == 0.0 && self.rad == 0.0 {
                return (y.sign() == Some(Ordering::Greater)).then_some(MachineBall::ZERO);
            }
            let gap = down(self.mid - self.rad, self.rad == 0.0);
            if gap == 0.0 {
                return None;
            }
            let (v, err, l) = kernels::pow(self.mid, y.mid)?;
            let lambda = div_up(self.rad, gap);
            let log_m = kernel_bounds(l).1;
            let delta = add_up(
                mul_up(y.mid.abs(), lambda),
                mul_up(y.rad, add_up(log_m, lambda)),
            );
            if delta > f64::from_bits((1023 - 20) << 52) {
                return None;
            }
            let grow = mul_up(delta, add_up(1.0, 2.0 * delta));
            MachineBall::from_kernel_err(v, err, |size| mul_up(size, grow))
        };
        fast().or_else(|| MachineBall::from_ball(&self.to_ball().pow_ball(&y.to_ball(), 64)))
    }

    /// The real `n`-th root of every point of the ball, as [`Ball::root`]:
    /// for `x > 0` the kernel `exp(log(x) / n)` (exact for an exact perfect
    /// power), spread `|v| r / (n (m - r))` (`(x^(1/n))'` falls with `x`);
    /// odd roots of negative balls by symmetry; `Ball` across 0, `None` for
    /// an even root of a ball reaching below 0.
    ///
    /// # Panics
    ///
    /// If `n` is 0.
    #[must_use]
    pub fn root(&self, n: u64) -> Option<MachineBall> {
        assert!(n >= 1, "MachineBall::root: 0th root");
        match (n, self.sign()) {
            (1, _) => Some(*self),
            (2, _) => self.sqrt(),
            (_, Some(Ordering::Equal)) => Some(MachineBall::ZERO),
            (_, Some(Ordering::Greater)) => self.root_positive(n),
            (_, Some(Ordering::Less)) if n % 2 == 1 => self.neg().root_positive(n).map(|b| b.neg()),
            _ => self.via_ball(|b, p| b.root(n, p)),
        }
    }

    /// [`MachineBall::root`] of a ball in `x > 0`, `n >= 3`.
    fn root_positive(&self, n: u64) -> Option<MachineBall> {
        let fast = || {
            let (v, err) = kernels::root(self.mid, n)?;
            let r = self.rad;
            if r == 0.0 {
                // A perfect power: c^n = m exactly for the nearest double c.
                let c = MachineBall::exact(v.to_f64())?;
                let k = i32::try_from(n).ok();
                if k.and_then(|k| c.pow(k)) == Some(*self) {
                    return Some(c);
                }
                return MachineBall::from_kernel_err(v, err, |_| 0.0);
            }
            #[allow(clippy::cast_precision_loss, reason = "n < 2^53 (the kernel's range)")]
            let den = mul_down(n as f64, down(self.mid - r, false));
            if den == 0.0 {
                return None;
            }
            MachineBall::from_kernel_err(v, err, |size| div_up(mul_up(size, r), den))
        };
        fast().or_else(|| self.via_ball(|b, p| b.root(n, p)))
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

/// `a × b` from below, for `a, b >= 0`.
#[inline]
fn mul_down(a: f64, b: f64) -> f64 {
    down(a * b, a == 0.0 || b == 0.0)
}

/// For a kernel value `v` and its nearest double `mid`: `|v - mid|` and
/// `|v|` from above (`hi - mid` is exact by Sterbenz, the sum with `lo`
/// rounds once).
fn kernel_size(v: dd::Dd, mid: f64) -> (f64, f64) {
    let d = (v.hi - mid) + v.lo;
    let round = up(d.abs(), d == 0.0);
    (round, add_up(mid.abs(), round))
}

/// Bounds `lo <= |f| <= hi` on the true value `f` of a kernel value `v`
/// (relative error at most `KERNEL_ERR`, so `|f - v| <= 2 KERNEL_ERR |v|`).
fn kernel_bounds(v: dd::Dd) -> (f64, f64) {
    let mid = v.to_f64();
    let (round, size) = kernel_size(v, mid);
    let err = mul_up(size, 2.0 * kernels::KERNEL_ERR);
    let lo = down(down(mid.abs() - round, false) - err, false);
    (lo, add_up(size, err))
}

/// π from above.
const PI_UP: f64 = 3.141_592_653_589_793_6;

/// `x` as an `i32` if it is an integer in range.
#[allow(
    clippy::cast_possible_truncation,
    reason = "checked to be an integer in range"
)]
fn exact_i32(x: f64) -> Option<i32> {
    (x.trunc() == x && x.abs() <= f64::from(i32::MAX)).then_some(x as i32)
}

/// Knuth's two-sum: `(s, e)` with `s = fl(a + b)` and `s + e = a + b`
/// exactly when `s` is finite.
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}
