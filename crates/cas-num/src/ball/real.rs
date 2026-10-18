//! `Ball`: a real number enclosed as midpoint ± radius (CAS-DOC-5 §5).
//!
//! A ball `[m ± r]` stands for every real `x` with `|x - m| <= r`; the one
//! promise of every operation is that the result contains the exact result
//! for every choice of points in the inputs. Midpoints are [`Float`]s, rounded
//! toward zero at the requested precision (or computed by the approximate
//! entries, §5.2), and the rounding error goes into the [`Mag`] radius, so an
//! exact result stays exact (radius 0). An infinite radius is the
//! indeterminate ball (midpoint 0), the only "no information" state:
//! division by a ball containing 0, the square root of a ball reaching below
//! 0, or a radius that overflowed.
//!
//! Precision is passed to every operation (§5.1); [`RealBalls`] wraps it in a
//! ring object.

use std::cmp::Ordering;
use std::fmt;

use super::{Float, Mag, Round};
use crate::{Integer, Rational};

mod ring;
#[cfg(test)]
mod tests;

pub use ring::RealBalls;

/// A midpoint-radius enclosure of a real number.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Ball {
    mid: Float,
    rad: Mag,
}

/// `k` units of `2^(exp(x) - prec)` as a radius term (0 for `k = 0`).
fn ulps(x: &Float, prec: u64, k: u64) -> Mag {
    if k == 0 || x.is_zero() {
        return Mag::ZERO;
    }
    let ulp = Mag::power_of_two(&(x.exponent() - &Integer::from(prec)));
    match k {
        1 => ulp,
        2 => ulp.mul_2exp(&Integer::ONE),
        _ => ulp.mul(&Mag::from_integer(&Integer::from(k))),
    }
}

/// Compares `|x - a|` with the non-negative exact `r`: at 64 bits rounded
/// both ways first, exactly only when that cannot decide.
fn cmp_dist(x: &Float, a: &Float, r: &Float) -> Ordering {
    let hi = x.sub_round(a, 64, Round::Up).0.abs();
    if hi < *r {
        return Ordering::Less;
    }
    let lo = x.sub_round(a, 64, Round::Down).0.abs();
    if lo > *r {
        return Ordering::Greater;
    }
    x.sub(a).abs().cmp(r)
}

impl Ball {
    /// The exact zero ball.
    pub const ZERO: Ball = Ball {
        mid: Float::ZERO,
        rad: Mag::ZERO,
    };

    /// The indeterminate ball `[0 ± ∞]`: nothing known.
    pub const INDETERMINATE: Ball = Ball {
        mid: Float::ZERO,
        rad: Mag::INF,
    };

    /// `[mid ± rad]` (the indeterminate ball if `rad` is `+∞`).
    #[must_use]
    pub fn new(mid: Float, rad: Mag) -> Ball {
        if rad.is_inf() {
            Ball::INDETERMINATE
        } else {
            Ball { mid, rad }
        }
    }

    /// The exact ball `[x ± 0]`.
    #[must_use]
    pub fn exact(x: Float) -> Ball {
        Ball {
            mid: x,
            rad: Mag::ZERO,
        }
    }

    /// The exact ball of an integer.
    #[must_use]
    pub fn from_integer(n: &Integer) -> Ball {
        Ball::exact(Float::from(n))
    }

    /// A ball containing `q`: the midpoint rounded to nearest at `prec` bits,
    /// the radius half an ulp, or 0 when `q` is representable.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn from_rational(q: &Rational, prec: u64) -> Ball {
        let (m, o) = Float::from_rational(q, prec, Round::Nearest);
        Ball::from_nearest(m, o, prec)
    }

    /// The ball of a value rounded to nearest at `prec` bits as `m`, with
    /// ternary value `o`: radius half an ulp, 0 when `o` is `Equal`.
    pub(crate) fn from_nearest(m: Float, o: Ordering, prec: u64) -> Ball {
        let rad = if o == Ordering::Equal {
            Mag::ZERO
        } else {
            // |q - m| <= 2^(exp(q) - prec - 1) <= 2^(exp(m) - prec - 1).
            ulps(&m, prec + 1, 1)
        };
        Ball { mid: m, rad }
    }

    #[must_use]
    pub fn mid(&self) -> &Float {
        &self.mid
    }

    #[must_use]
    pub fn rad(&self) -> &Mag {
        &self.rad
    }

    #[must_use]
    pub fn into_parts(self) -> (Float, Mag) {
        (self.mid, self.rad)
    }

    /// Whether the radius is 0.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.rad.is_zero()
    }

    /// Whether the radius is finite (not the indeterminate ball).
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.rad.is_finite()
    }

    /// The ball with `err` added to its radius: it then also contains every
    /// point within `err` of a point of `self`.
    #[must_use]
    pub fn add_error(&self, err: &Mag) -> Ball {
        Ball::new(self.mid.clone(), self.rad.add(err))
    }

    /// `[m ± r]` from a rounded midpoint whose error is `k` ulps (`rad` the
    /// propagated radius).
    #[allow(
        clippy::needless_pass_by_value,
        reason = "callers build the radius as a temporary"
    )]
    fn with_error(mid: Float, rad: Mag, prec: u64, k: u64) -> Ball {
        let err = ulps(&mid, prec, k);
        Ball::new(mid, rad.add(&err))
    }

    /// `self + rhs` at `prec` bits: radius `r + s + ε`.
    #[must_use]
    pub fn add(&self, rhs: &Ball, prec: u64) -> Ball {
        if !self.is_finite() || !rhs.is_finite() {
            return Ball::INDETERMINATE;
        }
        let (m, o) = self.mid.add_round(&rhs.mid, prec, Round::Down);
        Ball::with_error(
            m,
            self.rad.add(&rhs.rad),
            prec,
            u64::from(o != Ordering::Equal),
        )
    }

    /// `self - rhs` at `prec` bits: radius `r + s + ε`.
    #[must_use]
    pub fn sub(&self, rhs: &Ball, prec: u64) -> Ball {
        if !self.is_finite() || !rhs.is_finite() {
            return Ball::INDETERMINATE;
        }
        let (m, o) = self.mid.sub_round(&rhs.mid, prec, Round::Down);
        Ball::with_error(
            m,
            self.rad.add(&rhs.rad),
            prec,
            u64::from(o != Ordering::Equal),
        )
    }

    /// `self × rhs` at `prec` bits: radius `|a| s + |b| r + r s + ε`.
    #[must_use]
    pub fn mul(&self, rhs: &Ball, prec: u64) -> Ball {
        if !self.is_finite() || !rhs.is_finite() {
            return Ball::INDETERMINATE;
        }
        let (m, o) = self.mid.mul_round(&rhs.mid, prec, Round::Down);
        let (r, s) = (&self.rad, &rhs.rad);
        let mut rad = Mag::ZERO;
        if !s.is_zero() {
            rad = rad.add(&Mag::from_float(&self.mid).mul(s));
        }
        if !r.is_zero() {
            rad = rad.add(&Mag::from_float(&rhs.mid).mul(r)).add(&r.mul(s));
        }
        Ball::with_error(m, rad, prec, u64::from(o != Ordering::Equal))
    }

    /// `self / rhs` at `prec` bits: radius `(|a| s + |b| r) / (|b| (|b| - s)) +
    /// ε`, with lower bounds in the denominator; indeterminate if `rhs`
    /// contains 0.
    #[must_use]
    pub fn div(&self, rhs: &Ball, prec: u64) -> Ball {
        if !self.is_finite() || !rhs.is_finite() {
            return Ball::INDETERMINATE;
        }
        // |b| - s, from below: 0 when rhs may contain 0.
        let gap = Mag::from_float_lower(&rhs.mid).sub_lower(&rhs.rad);
        if gap.is_zero() {
            return Ball::INDETERMINATE;
        }
        let Some((q, k)) = self.mid.div_mid(&rhs.mid, prec) else {
            return Ball::INDETERMINATE;
        };
        let (r, s) = (&self.rad, &rhs.rad);
        let rad = if r.is_zero() && s.is_zero() {
            Mag::ZERO
        } else {
            let num = Mag::from_float(&self.mid)
                .mul(s)
                .add(&Mag::from_float(&rhs.mid).mul(r));
            num.div(&Mag::from_float_lower(&rhs.mid).mul_lower(&gap))
        };
        Ball::with_error(q, rad, prec, k)
    }

    /// `1 / self` at `prec` bits.
    #[must_use]
    pub fn inv(&self, prec: u64) -> Ball {
        Ball::from_integer(&Integer::ONE).div(self, prec)
    }

    /// `sqrt(self)` at `prec` bits: radius `r / (sqrt(a - r) + sqrt(a)) + ε`,
    /// with lower bounds; indeterminate if the ball reaches below 0.
    #[must_use]
    pub fn sqrt(&self, prec: u64) -> Ball {
        if !self.is_finite() || self.mid.is_negative() {
            return Ball::INDETERMINATE;
        }
        let r = &self.rad;
        if !r.is_zero() {
            // Needs a - r >= 0, decided exactly.
            match r.to_float() {
                Some(rf) if self.mid.cmp_abs(&rf) != Ordering::Less => {}
                _ => return Ball::INDETERMINATE,
            }
        }
        let Some((s, k)) = self.mid.sqrt_mid(prec) else {
            return Ball::INDETERMINATE;
        };
        let rad = if r.is_zero() {
            Mag::ZERO
        } else {
            let a = Mag::from_float_lower(&self.mid);
            let den = a.sub_lower(r).sqrt_lower().add_lower(&a.sqrt_lower());
            r.div(&den)
        };
        Ball::with_error(s, rad, prec, k)
    }

    /// `self^n` at `prec` bits (`self^0 = 1`), by squaring and multiplying
    /// balls; a negative `n` inverts the power.
    #[must_use]
    pub fn pow(&self, n: i64, prec: u64) -> Ball {
        let k = n.unsigned_abs();
        let mut r = Ball::from_integer(&Integer::ONE);
        for i in (0..u64::BITS - k.leading_zeros()).rev() {
            r = r.mul(&r, prec);
            if (k >> i) & 1 == 1 {
                r = r.mul(self, prec);
            }
        }
        if n < 0 { r.inv(prec) } else { r }
    }

    /// `[|m| ± r]`, which contains `|x|` for every `x` in the ball.
    #[must_use]
    pub fn abs(&self) -> Ball {
        Ball {
            mid: self.mid.abs(),
            rad: self.rad.clone(),
        }
    }

    /// `-self`, exactly.
    #[must_use]
    pub fn neg(&self) -> Ball {
        Ball {
            mid: self.mid.neg(),
            rad: self.rad.clone(),
        }
    }

    /// `self × 2^k`, exactly.
    #[must_use]
    pub fn mul_2exp(&self, k: &Integer) -> Ball {
        Ball {
            mid: self.mid.mul_2exp(k),
            rad: self.rad.mul_2exp(k),
        }
    }

    /// Whether `x` lies in the ball.
    #[must_use]
    pub fn contains(&self, x: &Float) -> bool {
        match self.rad.to_float() {
            None => true,
            Some(r) => cmp_dist(x, &self.mid, &r) != Ordering::Greater,
        }
    }

    /// Whether `n` lies in the ball.
    #[must_use]
    pub fn contains_integer(&self, n: &Integer) -> bool {
        self.contains(&Float::from(n))
    }

    /// Whether `q` lies in the ball: `q` is bracketed by `p`-bit floats for
    /// `p = 64, 128, …` until both lie in the ball, or one lies outside on its
    /// side, or `q` is exact (a dyadic `q` is reached at its bit length; a
    /// non-dyadic one never equals an endpoint, which are dyadic).
    #[must_use]
    pub fn contains_rational(&self, q: &Rational) -> bool {
        if !self.is_finite() {
            return true;
        }
        let mut p = 64;
        loop {
            let (lo, o) = Float::from_rational(q, p, Round::Floor);
            if o == Ordering::Equal {
                return self.contains(&lo);
            }
            let (hi, _) = Float::from_rational(q, p, Round::Ceil);
            let (in_lo, in_hi) = (self.contains(&lo), self.contains(&hi));
            if in_lo && in_hi {
                return true;
            }
            if (!in_lo && lo > self.mid) || (!in_hi && hi < self.mid) {
                return false;
            }
            p *= 2;
        }
    }

    /// Whether `other` lies entirely in `self`: `|b - a| + s <= r`.
    #[must_use]
    pub fn contains_ball(&self, other: &Ball) -> bool {
        let Some(r) = self.rad.to_float() else {
            return true;
        };
        let Some(s) = other.rad.to_float() else {
            return false;
        };
        if s > r {
            return false;
        }
        let below = r.sub_round(&s, 64, Round::Down).0;
        if cmp_dist(&other.mid, &self.mid, &below) != Ordering::Greater {
            return true;
        }
        let above = r.sub_round(&s, 64, Round::Up).0;
        if cmp_dist(&other.mid, &self.mid, &above) == Ordering::Greater {
            return false;
        }
        cmp_dist(&other.mid, &self.mid, &r.sub(&s)) != Ordering::Greater
    }

    /// Whether the balls share a point: `|a - b| <= r + s`.
    #[must_use]
    pub fn overlaps(&self, other: &Ball) -> bool {
        let (Some(r), Some(s)) = (self.rad.to_float(), other.rad.to_float()) else {
            return true;
        };
        let below = r.add_round(&s, 64, Round::Down).0;
        if cmp_dist(&other.mid, &self.mid, &below) != Ordering::Greater {
            return true;
        }
        let above = r.add_round(&s, 64, Round::Up).0;
        if cmp_dist(&other.mid, &self.mid, &above) == Ordering::Greater {
            return false;
        }
        cmp_dist(&other.mid, &self.mid, &r.add(&s)) != Ordering::Greater
    }

    /// The sign of every point of the ball, if they share one: `Some(Equal)`
    /// only for the exact zero ball, `None` when the ball contains 0 and
    /// something else.
    #[must_use]
    pub fn sign(&self) -> Option<Ordering> {
        let r = self.rad.to_float()?;
        if self.mid.is_zero() {
            return r.is_zero().then_some(Ordering::Equal);
        }
        (self.mid.cmp_abs(&r) == Ordering::Greater).then(|| self.mid.sign())
    }

    /// How every point of `self` compares with every point of `other`, if
    /// that is decided: the balls are disjoint, or both are the same exact
    /// value.
    #[must_use]
    #[allow(
        clippy::should_implement_trait,
        reason = "three-valued: Ball is deliberately not Ord"
    )]
    pub fn cmp(&self, other: &Ball) -> Option<Ordering> {
        if self.is_exact() && other.is_exact() {
            return Some(self.mid.cmp(&other.mid));
        }
        (!self.overlaps(other)).then(|| self.mid.cmp(&other.mid))
    }

    /// The endpoints `m - r` and `m + r` rounded outward at `prec` bits.
    fn endpoints(&self, prec: u64) -> (Float, Float) {
        let r = self.rad.to_float().expect("finite");
        (
            self.mid.sub_round(&r, prec, Round::Floor).0,
            self.mid.add_round(&r, prec, Round::Ceil).0,
        )
    }

    /// The ball of the interval `[lo, hi]` (`lo <= hi`), at `prec` bits.
    fn from_endpoints(lo: &Float, hi: &Float, prec: u64) -> Ball {
        let m = lo
            .add_round(hi, prec, Round::Down)
            .0
            .mul_2exp(&-Integer::ONE);
        let up = |x: &Float, y: &Float| Mag::from_float(&x.sub_round(y, 64, Round::Up).0);
        let rad = up(hi, &m).max(up(&m, lo));
        Ball::new(m, rad)
    }

    /// A ball containing both balls, at `prec` bits.
    #[must_use]
    pub fn union(&self, other: &Ball, prec: u64) -> Ball {
        if !self.is_finite() || !other.is_finite() {
            return Ball::INDETERMINATE;
        }
        let ((a, b), (c, d)) = (self.endpoints(prec), other.endpoints(prec));
        Ball::from_endpoints(std::cmp::min(&a, &c), std::cmp::max(&b, &d), prec)
    }

    /// A ball containing every point the two balls share, at `prec` bits;
    /// `None` if they are disjoint.
    #[must_use]
    pub fn intersection(&self, other: &Ball, prec: u64) -> Option<Ball> {
        if !self.is_finite() {
            return Some(other.clone());
        }
        if !other.is_finite() {
            return Some(self.clone());
        }
        let ((a, b), (c, d)) = (self.endpoints(prec), other.endpoints(prec));
        // Outward rounding only grows the intervals: empty here is empty.
        let (lo, hi) = (std::cmp::max(&a, &c), std::cmp::min(&b, &d));
        (lo <= hi).then(|| Ball::from_endpoints(lo, hi, prec))
    }

    /// The number of bits known relative to the value, about
    /// `exp(m) - exp(r)`: `i64::MAX` for an exact ball, `i64::MIN` for the
    /// indeterminate ball or a ball centred at 0 with a non-zero radius.
    #[must_use]
    pub fn rel_accuracy_bits(&self) -> i64 {
        if self.rad.is_zero() {
            return i64::MAX;
        }
        if self.rad.is_inf() || self.mid.is_zero() {
            return i64::MIN;
        }
        let d = self.mid.exponent() - self.rad.exponent();
        d.to_i64()
            .unwrap_or(if d.is_negative() { i64::MIN } else { i64::MAX })
    }

    /// The ball with its midpoint rounded to nearest at `bits` bits and the
    /// radius grown to keep the enclosure (for `N`'s result, §9.2).
    ///
    /// # Panics
    ///
    /// If `bits` is 0.
    #[must_use]
    pub fn round_to_accuracy(&self, bits: u64) -> Ball {
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let (m, o) = self.mid.round(bits, Round::Nearest);
        let k = u64::from(o != Ordering::Equal);
        // |m - mid| <= 2^(exp(mid) - bits - 1) <= 2^(exp(m) - bits - 1).
        let err = ulps(&m, bits + 1, k);
        Ball::new(m, self.rad.add(&err))
    }

    /// `*self = *self + rhs` at `prec` bits.
    pub fn add_assign(&mut self, rhs: &Ball, prec: u64) {
        *self = self.add(rhs, prec);
    }

    /// `*self = *self - rhs` at `prec` bits.
    pub fn sub_assign(&mut self, rhs: &Ball, prec: u64) {
        *self = self.sub(rhs, prec);
    }

    /// `*self = *self × rhs` at `prec` bits.
    pub fn mul_assign(&mut self, rhs: &Ball, prec: u64) {
        *self = self.mul(rhs, prec);
    }

    /// `*self = *self / rhs` at `prec` bits.
    pub fn div_assign(&mut self, rhs: &Ball, prec: u64) {
        *self = self.div(rhs, prec);
    }
}

impl Default for Ball {
    fn default() -> Self {
        Ball::ZERO
    }
}

/// `[mid ± rad]` with the `Debug` forms of both.
impl fmt::Debug for Ball {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{:?} ± {:?}]", self.mid, self.rad)
    }
}
