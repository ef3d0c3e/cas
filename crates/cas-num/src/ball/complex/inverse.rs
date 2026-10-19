//! `asin`, `acos`, `asinh`, `acosh` and `atanh` of complex balls (CAS-102).
//!
//! The first four follow Hull, Fairgrieve and Tang ("Implementing the
//! complex arcsine and arccosine functions using exception handling", ACM
//! TOMS 23, 1997): with `X = |x|`, `Y = |y|` for `z = x + yi`,
//!
//! ```text
//! A = (|z + 1| + |z - 1|) / 2 >= 1,   D = √(A^2 - X^2) = √((A - X)(A + X))
//! asin z  = atan2(x, D)          + i sgn(y) acosh A
//! acos z  = atan2(D, x)          - i sgn(y) acosh A
//! acosh z = acosh A              + i sgn(y) atan2(D, x)
//! asinh z = asinh's A, D at (|y|, |x|):  sgn(x) acosh A + i atan2(y, D)
//! acosh A = log1p((A - 1) + √((A - 1)(A + 1)))
//! ```
//!
//! `A` and `D` are continuous in `z`, so the only discontinuity is the
//! sign, across the cuts (`y = 0` with `X > 1` for `asin`, `acos`, where
//! `acosh A > 0`, and with `x < 1` for `acosh`, where `atan2(D, x) > 0`;
//! `x = 0`, `|y| > 1` for `asinh`). On a cut it is the `+0` side, and a
//! ball across the axis gets `[0 ± v]` for the part `sgn(y) v` (or `sgn(x)
//! v`), which holds both sides. `A - 1` and `A - X`
//! never cancel: with `r = |z + 1|`, `s = |z - 1|` and `r - (X + 1) = Y^2 /
//! (r + X + 1)`,
//!
//! ```text
//! X < 1:  A - 1 = (Y^2/(r + X + 1) + Y^2/(s + 1 - X)) / 2
//!         A - X = (Y^2/(r + X + 1) + s + 1 - X) / 2
//! X > 1:  A - 1 = (Y^2/(r + X + 1) + s + X - 1) / 2
//!         A - X = (Y^2/(r + X + 1) + Y^2/(s + X - 1)) / 2
//! ```
//!
//! (a ball across `X = 1` takes the sums `s ± (X - 1)` of both, whose
//! cancellation is bounded by the ball's width). Near 0, `acosh A ≈ |y|`
//! and `atan2(x, D) ≈ x` keep their relative accuracy; near `±1`, both are
//! as accurate as `A - 1` and `A - X`.
//!
//! `atanh` is direct, with `X = |x|`:
//!
//! ```text
//! Re atanh z = sgn(x) log1p(4X / ((1 - X)^2 + y^2)) / 4
//! Im atanh z = atan2(2y, (1 - X)(1 + X) - y^2) / 2
//! ```
//!
//! whose real part keeps its relative accuracy near 0 (and, by the
//! symmetry, near `-1`); the quotient, `>= 0`, is clamped to `[0, …]` when
//! its ball reaches below 0. The imaginary part's cancellation near `|z| =
//! 1` only moves an angle near `±π/2`. On the cut `atan2(+0, ·) = π` gives
//! the `+0` side, across it `atan2` gives `[0 ± π]`.

#![allow(
    clippy::many_single_char_names,
    reason = "x, y, r, s, a, d as in the formulas of the module comment"
)]

use std::cmp::Ordering;

use super::{ComplexBall, clamp0, fold, half, is_zero, one, sqrt_nonneg, square, unfold};
use crate::Integer;
use crate::ball::Ball;

/// Whether the ball lies in `[-1, 1]`, decided (the real `asin` and `acos`
/// domain).
fn within_one(x: &Ball) -> bool {
    matches!(
        one().sub(&x.abs(), 64).sign(),
        Some(Ordering::Greater | Ordering::Equal)
    )
}

/// `(D, acosh A)` of the module comment for `X = |x|`, `Y = |y|` (balls of
/// values `>= 0`), at `wp` bits.
fn hft(x: &Ball, y: &Ball, wp: u64) -> (Ball, Ball) {
    let y2 = square(y, wp);
    let xp1 = x.add(&one(), wp);
    let xm1 = x.sub(&one(), wp);
    let r = square(&xp1, wp).add(&y2, wp).sqrt(wp);
    let s = sqrt_nonneg(&square(&xm1, wp).add(&y2, wp), wp);
    let a = half(&r.add(&s, wp));
    let rx = y2.div(&r.add(&xp1, wp), wp);
    let (am1, amx) = match xm1.sign() {
        Some(Ordering::Less) => {
            let sp = s.sub(&xm1, wp);
            (rx.add(&y2.div(&sp, wp), wp), rx.add(&sp, wp))
        }
        Some(Ordering::Greater) => {
            let sp = s.add(&xm1, wp);
            (rx.add(&sp, wp), rx.add(&y2.div(&sp, wp), wp))
        }
        _ => (rx.add(&s.add(&xm1, wp), wp), rx.add(&s.sub(&xm1, wp), wp)),
    };
    let (am1, amx) = (half(&am1), half(&amx));
    let d = sqrt_nonneg(&amx.mul(&a.add(x, wp), wp), wp);
    let root = sqrt_nonneg(&am1.mul(&a.add(&one(), wp), wp), wp);
    (d, am1.add(&root, wp).log1p(wp))
}

impl ComplexBall {
    /// The principal `asin(self)`, cuts on the real axis beyond `±1`.
    #[must_use]
    pub fn asin(&self, prec: u64) -> ComplexBall {
        let (x, y) = (&self.re, &self.im);
        if self.is_real() && within_one(x) {
            return ComplexBall::from_real(x.asin(prec));
        }
        if is_zero(x) {
            // asin(iy) = i asinh(y).
            return ComplexBall::new(Ball::ZERO, y.asinh(prec));
        }
        let (d, acosh_a) = hft(&fold(x).0, &fold(y).0, prec + 16);
        ComplexBall::new(
            x.atan2(&d, prec),
            unfold(acosh_a, y.sign()).round_to_accuracy(prec),
        )
    }

    /// The principal `acos(self)`, cuts on the real axis beyond `±1`.
    #[must_use]
    pub fn acos(&self, prec: u64) -> ComplexBall {
        let (x, y) = (&self.re, &self.im);
        if self.is_real() && within_one(x) {
            return ComplexBall::from_real(x.acos(prec));
        }
        let (d, acosh_a) = hft(&fold(x).0, &fold(y).0, prec + 16);
        ComplexBall::new(
            d.atan2(x, prec),
            unfold(acosh_a, y.sign()).neg().round_to_accuracy(prec),
        )
    }

    /// The principal `asinh(self)`, cuts on the imaginary axis beyond `±i`.
    #[must_use]
    pub fn asinh(&self, prec: u64) -> ComplexBall {
        let (x, y) = (&self.re, &self.im);
        if self.is_real() {
            return ComplexBall::from_real(x.asinh(prec));
        }
        if is_zero(x) && within_one(y) {
            // asinh(iy) = i asin(y).
            return ComplexBall::new(Ball::ZERO, y.asin(prec));
        }
        let (d, acosh_a) = hft(&fold(y).0, &fold(x).0, prec + 16);
        ComplexBall::new(
            unfold(acosh_a, x.sign()).round_to_accuracy(prec),
            y.atan2(&d, prec),
        )
    }

    /// The principal `acosh(self)` (real part `>= 0`), cut on the real axis
    /// below 1.
    #[must_use]
    pub fn acosh(&self, prec: u64) -> ComplexBall {
        let (x, y) = (&self.re, &self.im);
        let one = one();
        if self.is_real() && matches!(x.cmp(&one), Some(Ordering::Greater | Ordering::Equal)) {
            return ComplexBall::from_real(x.acosh(prec));
        }
        let wp = prec + 16;
        let (d, acosh_a) = hft(&fold(x).0, &fold(y).0, wp);
        ComplexBall::new(
            acosh_a.round_to_accuracy(prec),
            unfold(d.atan2(x, wp), y.sign()).round_to_accuracy(prec),
        )
    }

    /// The principal `atanh(self)`, cuts on the real axis beyond `±1`;
    /// indeterminate if the rectangle reaches `±1`.
    #[must_use]
    pub fn atanh(&self, prec: u64) -> ComplexBall {
        let (x, y) = (&self.re, &self.im);
        let one = one();
        if self.is_real() && one.sub(&x.abs(), 64).sign() == Some(Ordering::Greater) {
            return ComplexBall::from_real(x.atanh(prec));
        }
        if is_zero(x) {
            // atanh(iy) = i atan(y).
            return ComplexBall::new(Ball::ZERO, y.atan(prec));
        }
        let wp = prec + 16;
        let (xa, sign) = fold(x);
        let omx = one.sub(&xa, wp);
        let y2 = square(y, wp);
        let ratio = clamp0(
            xa.mul_2exp(&Integer::from(2))
                .div(&square(&omx, wp).add(&y2, wp), wp),
        );
        let re = ratio.log1p(wp).mul_2exp(&Integer::from(-2));
        let den = omx.mul(&one.add(&xa, wp), wp).sub(&y2, wp);
        let im = half(&y.mul_2exp(&Integer::ONE).atan2(&den, wp));
        ComplexBall::new(
            unfold(re, sign).round_to_accuracy(prec),
            im.round_to_accuracy(prec),
        )
    }
}
