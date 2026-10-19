//! `expm1`, `log1p` and `root` of complex balls (CAS-102):
//!
//! ```text
//! expm1(a + bi) = (expm1(a) cos b + (cos b - 1)) + i e^a sin b,
//!                 cos b - 1 = -sin^2 b / (1 + cos b) while cos b > 0
//! log1p(a + bi) = log1p(a (2 + a) + b^2) / 2 + i atan2(b, 1 + a)
//!                 (log((1 + a)^2 + b^2) / 2 away from |1 + z| = 1)
//! root(z, n)    = |z|^(1/n) (cos(θ/n) + i sin(θ/n)),  θ = arg z
//! ```
//!
//! Near 0 both parts of `expm1` and `log1p` are within a few ulps of `|z|`:
//! every term comes from `a` and `b` themselves, not from a rounded `1 + z`
//! or `e^z`, and what still cancels (`expm1(a) cos b` against `1 - cos b`,
//! `a (2 + a)` against `b^2`) are terms of size `|z|` where the part itself
//! vanishes. `e^a` is `expm1(a) + 1` from `a = -1/2` up and `exp(a)` below,
//! where `expm1(a) + 1` would cancel. `log1p`'s cut is the real axis below
//! `-1` (`arg` in `(-π, π]`), `root`'s the negative real axis: the root is
//! the principal one, `(-8)^(1/3) = 1 + i√3`.

#![allow(
    clippy::many_single_char_names,
    reason = "a, b, s, c, n as in the formulas of the module comment"
)]

use std::cmp::Ordering;

use super::{ComplexBall, half, one, sqrt_nonneg, square, upper};
use crate::Integer;
use crate::ball::{Ball, Float, Mag};

impl ComplexBall {
    /// `exp(self) - 1`, without cancellation near 0.
    #[must_use]
    pub fn expm1(&self, prec: u64) -> ComplexBall {
        let (a, b) = (&self.re, &self.im);
        if self.is_real() {
            return ComplexBall::from_real(a.expm1(prec));
        }
        let wp = prec + 8;
        let minus_half = Float::from(-1i64).mul_2exp(&-Integer::ONE);
        let (em, ea) = if *a.mid() >= minus_half {
            let em = a.expm1(wp);
            let ea = em.add(&one(), wp);
            (em, ea)
        } else {
            let ea = a.exp(wp);
            (ea.sub(&one(), wp), ea)
        };
        let (s, c) = b.sin_cos(wp);
        let cm1 = if c.sign() == Some(Ordering::Greater) {
            square(&s, wp).div(&one().add(&c, wp), wp).neg()
        } else {
            c.sub(&one(), wp)
        };
        ComplexBall::new(em.mul(&c, wp).add(&cm1, prec), ea.mul(&s, prec))
    }

    /// `log(1 + self)`, without cancellation near 0; cut on the real axis
    /// below `-1`, indeterminate if the rectangle may reach `-1`.
    #[must_use]
    pub fn log1p(&self, prec: u64) -> ComplexBall {
        let (a, b) = (&self.re, &self.im);
        let wp = prec + 8;
        let w = a.add(&one(), wp);
        if self.is_real() && w.sign() == Some(Ordering::Greater) {
            return ComplexBall::from_real(a.log1p(prec));
        }
        if w.contains(&Float::ZERO) && b.contains(&Float::ZERO) {
            return ComplexBall::INDETERMINATE;
        }
        // |1 + z|^2 - 1 = a (2 + a) + b^2.
        let two = Ball::from_integer(&Integer::from(2));
        let b2 = square(b, wp);
        let t = a.mul(&a.add(&two, wp), wp).add(&b2, wp);
        let log_abs = if upper(&t) <= Mag::ONE.mul_2exp(&-Integer::ONE) {
            t.log1p(prec + 1)
        } else {
            square(&w, wp).add(&b2, wp).log(prec + 1)
        };
        ComplexBall::new(half(&log_abs), b.atan2(&w, prec))
    }

    /// The principal `n`-th root `exp(log(self) / n)`, cut on the negative
    /// real axis; a real `x >= 0` gives the real root (exact for an exact
    /// perfect power).
    ///
    /// # Panics
    ///
    /// If `n` is 0.
    #[must_use]
    pub fn root(&self, n: u64, prec: u64) -> ComplexBall {
        assert!(n >= 1, "ComplexBall::root: 0th root");
        match n {
            1 => return self.round(prec),
            2 => return self.sqrt(prec),
            _ => {}
        }
        let (a, b) = (&self.re, &self.im);
        if self.is_real() && matches!(a.sign(), Some(Ordering::Greater | Ordering::Equal)) {
            return ComplexBall::from_real(a.root(n, prec));
        }
        let wp = prec + 16;
        if self.may_be_zero() {
            // Every root has modulus at most (|a| + |b|)^(1/n).
            let Some(u) = upper(a).add(&upper(b)).to_float() else {
                return ComplexBall::INDETERMINATE;
            };
            let bound = upper(&Ball::exact(u).root(n, 64));
            return ComplexBall::new(
                Ball::new(Float::ZERO, bound.clone()),
                Ball::new(Float::ZERO, bound),
            );
        }
        let modulus = sqrt_nonneg(&square(a, wp).add(&square(b, wp), wp), wp).root(n, wp);
        let theta = self.arg(wp).div(&Ball::from_integer(&Integer::from(n)), wp);
        let (s, c) = theta.sin_cos(wp);
        ComplexBall::new(modulus.mul(&c, prec), modulus.mul(&s, prec))
    }
}
