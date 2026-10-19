//! `tan`, `sinh`, `cosh`, `tanh`, `sin_pi` and `cos_pi` of complex balls
//! (CAS-102), from the real functions of the parts:
//!
//! ```text
//! tan(a + bi) = (sin a cos a + i sinh b cosh b) / (cos^2 a + sinh^2 b)
//!             = (4 e s c + i sgn(b) t (t + 2)) / (4 e c^2 + t^2),
//!               t = expm1(2|b|), e = t + 1, (s, c) = (sin a, cos a)
//! tanh z = -i tan(iz),  sinh z = -i sin(iz),  cosh z = cos(iz)
//! sin(π(a + bi)) = sin(πa) cosh(πb) + i cos(πa) sinh(πb)
//! cos(π(a + bi)) = cos(πa) cosh(πb) - i sin(πa) sinh(πb)
//! ```
//!
//! `tan` as `sin z / cos z` would cancel in its real part for large `|b|`
//! (`sin a cos a (cosh^2 b - sinh^2 b)` as a difference of two numbers near
//! `e^(2|b|)`); the form above (ours, CAS-102) has only sums of
//! non-negative terms, its denominator included (`cos 2a + cosh 2b` would
//! cancel next to the poles `a = π/2 + kπ`, `b = 0`). Beyond `exp`'s range,
//! `tan z` is `±i` within `2^(4 - 2|b|)` in each part. `sin_pi` and `cos_pi`
//! take the real functions' exact values at integers and half-integers.
//!
//! An inexact ball is evaluated at its exact midpoint and widened by a bound
//! on `|f'|` over the rectangle times its radii, from 64-bit balls of the
//! parts:
//!
//! ```text
//! |tan'| = 1 / (cos^2 a + sinh^2 b)   |sinh'| = √(sinh^2 a + cos^2 b)
//! |sin_pi'| = π √(cos^2(πa) + sinh^2(πb))   |cosh'| = √(sinh^2 a + sin^2 b)
//! ```
//!
//! (`|cos_pi'|` with `sin`). Composing the formulas on the balls instead
//! loses the correlations: `tan(-0.3 + 4i)`, which is flat next to `i`, came
//! out 4000 times wider than `|tan'| ρ`, and the real `cos`, `cosh` and
//! `sin_pi` bound their slope by 1, 1 and π where it vanishes.

#![allow(
    clippy::many_single_char_names,
    reason = "a, b, s, c, t, e as in the formulas of the module comment"
)]

use super::{ComplexBall, fold, is_zero, lower, one, square, unfold, upper};
use crate::Integer;
use crate::ball::{Ball, Float, Mag};

impl ComplexBall {
    /// `tan(self)`; indeterminate when the rectangle reaches a pole.
    #[must_use]
    pub fn tan(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.tan(prec));
        }
        if is_zero(&self.re) {
            return ComplexBall::new(Ball::ZERO, self.im.tanh(prec));
        }
        self.at_midpoint(prec, tan_at, || {
            // |tan'| = 1 / |cos w|^2 = 1 / (cos^2 a + sinh^2 b), with |cos a|
            // from below by the cosine ball and sinh |b| >= |b|.
            let c = lower(&self.re.cos(64));
            let b = lower(&self.im);
            let s = b
                .to_float()
                .map_or(Mag::ZERO, |f| lower(&Ball::exact(f).sinh(64)))
                .max(b);
            Mag::ONE.div(&c.mul_lower(&c).add_lower(&s.mul_lower(&s)))
        })
    }

    /// `tanh(self) = -i tan(i self)`; indeterminate when the rectangle
    /// reaches a pole.
    #[must_use]
    pub fn tanh(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.tanh(prec));
        }
        let t = self.mul_i().tan(prec);
        ComplexBall::new(t.im, t.re.neg())
    }

    /// `sinh(self) = -i sin(i self)`.
    #[must_use]
    pub fn sinh(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.sinh(prec));
        }
        let f = |m: &ComplexBall, p| {
            let s = m.mul_i().sin(p);
            ComplexBall::new(s.im, s.re.neg())
        };
        // |sinh'| = |cosh w| = √(sinh^2 a + cos^2 b).
        self.at_midpoint(prec, f, || hypot(&self.re, &self.im.cos(64)))
    }

    /// `cosh(self) = cos(i self)`.
    #[must_use]
    pub fn cosh(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.cosh(prec));
        }
        // |cosh'| = |sinh w| = √(sinh^2 a + sin^2 b).
        self.at_midpoint(
            prec,
            |m, p| m.mul_i().cos(p),
            || hypot(&self.re, &self.im.sin(64)),
        )
    }

    /// `sin(π self)`: exact at real integers and half-integers.
    #[must_use]
    pub fn sin_pi(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.sin_pi(prec));
        }
        let f = |m: &ComplexBall, p| {
            let wp = p + 8;
            let (ch, sh) = m.cosh_sinh_pi_im(wp);
            ComplexBall::new(m.re.sin_pi(wp).mul(&ch, p), m.re.cos_pi(wp).mul(&sh, p))
        };
        // |sin_pi'| = π |cos(π w)| = π √(cos^2(πa) + sinh^2(πb)).
        self.at_midpoint(prec, f, || pi_hypot(&self.im, &self.re.cos_pi(64)))
    }

    /// `cos(π self)`: exact at real integers and half-integers.
    #[must_use]
    pub fn cos_pi(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.cos_pi(prec));
        }
        let f = |m: &ComplexBall, p| {
            let wp = p + 8;
            let (ch, sh) = m.cosh_sinh_pi_im(wp);
            ComplexBall::new(
                m.re.cos_pi(wp).mul(&ch, p),
                m.re.sin_pi(wp).mul(&sh, p).neg(),
            )
        };
        // |cos_pi'| = π |sin(π w)| = π √(sin^2(πa) + sinh^2(πb)).
        self.at_midpoint(prec, f, || pi_hypot(&self.im, &self.re.sin_pi(64)))
    }

    /// `(cosh(π b), sinh(π b))` for the imaginary part `b`.
    fn cosh_sinh_pi_im(&self, wp: u64) -> (Ball, Ball) {
        let pb = self.im.mul(&Ball::const_pi(wp), wp);
        (pb.cosh(wp), pb.sinh(wp))
    }
}

/// `tan` at an exact `z` by the formula of the module comment.
fn tan_at(z: &ComplexBall, prec: u64) -> ComplexBall {
    let wp = prec + 16;
    let (s, c) = z.re.sin_cos(wp);
    let (b, sign) = fold(&z.im);
    let t = b.mul_2exp(&Integer::ONE).expm1(wp);
    if !t.is_finite() {
        return tan_far(&b, sign);
    }
    let e4 = t.add(&one(), wp).mul_2exp(&Integer::from(2));
    let den = e4.mul(&square(&c, wp), wp).add(&square(&t, wp), wp);
    let re = e4.mul(&s.mul(&c, wp), wp).div(&den, prec);
    let two = Ball::from_integer(&Integer::from(2));
    let im = t.mul(&t.add(&two, wp), wp).div(&den, prec);
    ComplexBall::new(re, unfold(im, sign))
}

/// An upper bound of `√(sinh^2 a + v^2)` for `a` in the ball `re` and `|v| <=
/// min(1, |w|)` for `w` in the ball `trig` (`sinh |a| <= expm1 |a|`).
fn hypot(re: &Ball, trig: &Ball) -> Mag {
    let sh = upper(re).expm1();
    let v = upper(trig).min(Mag::ONE);
    sh.mul(&sh).add(&v.mul(&v)).sqrt()
}

/// An upper bound of `π √(v^2 + sinh^2(π b))` for `b` in the ball `im` and
/// `|v| <= min(1, |w|)` for `w` in the ball `trig`.
fn pi_hypot(im: &Ball, trig: &Ball) -> Mag {
    let pi = upper(&Ball::const_pi(64));
    let sh = upper(im).mul(&pi).expm1();
    let v = upper(trig).min(Mag::ONE);
    sh.mul(&sh).add(&v.mul(&v)).sqrt().mul(&pi)
}

/// `tan` for `|Im z| = b` beyond `exp`'s range, `sign` the sign of `Im z`.
///
/// For `b >= 1`, `tan z ∓ i` is `(sin a cos a + i ((1 - e^(-2b)) / 2 - cos^2
/// a)) / (cos^2 a + sinh^2 b)`: both parts of the numerator are at most 1,
/// and `sinh^2 b >= 0.18 e^(2b)`, so each part is within `5.4 e^(-2b) <=
/// 2^(4 - 2b)` of `±i`'s. Any bound above that will do: from `b >= 2^30`
/// the radius is `2^(4 - 2^31)`.
fn tan_far(b: &Ball, sign: Option<std::cmp::Ordering>) -> ComplexBall {
    let low = Mag::from_float_lower(b.mid()).sub_lower(b.rad());
    if low < Mag::power_of_two(&Integer::from(30)) {
        return ComplexBall::INDETERMINATE;
    }
    let err = Mag::power_of_two(&Integer::from(4 - (1i64 << 31)));
    ComplexBall::new(
        Ball::new(Float::ZERO, err.clone()),
        unfold(one().add_error(&err), sign),
    )
}
