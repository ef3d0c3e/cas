//! `ComplexBall`: rectangular complex balls `re + im i` (CAS-DOC-5 §8), Arb's
//! `acb_t` representation \[Joh17\].
//!
//! Every operation is built from real ball operations on the parts, so it
//! inherits their enclosure property; the formulas are chosen so that
//! nothing cancels where it need not:
//!
//! ```text
//! (a + bi)(c + di) = (ac - bd) + (ad + bc) i        (4 products; from
//!                   = (ac - bd) + ((a+b)(c+d) - ac - bd) i  MUL3_PREC bits,
//!                                                    3 products)
//! (a + bi) / (c + di) = ((ac + bd) + (bc - ad) i) / (c^2 + d^2)
//! exp(a + bi) = e^a (cos b + i sin b)
//! log z = log|z| + i arg z,  log|z| = log1p((a - 1)(a + 1) + b^2) / 2
//! sqrt z = t + i b/(2t) (a >= 0),  |b|/(2t) + i sgn(b) t (a < 0),
//!          t = √((|z| + |a|) / 2)
//! sin(a + bi) = sin a cosh b + i cos a sinh b
//! cos(a + bi) = cos a cosh b - i sin a sinh b
//! atan z = (log(1 + iz) - log(1 - iz)) / (2i)
//! z^w = exp(w log z), integer w by squaring
//! ```
//!
//! `sqrt` is Kahan's \[Kah87\], which never subtracts.
//!
//! The other elementary functions (CAS-102) are in [`trig`] (`tan`, the
//! hyperbolic functions, `sin_pi`, `cos_pi`), [`inverse`] (`asin`, `acos`
//! and the inverse hyperbolic functions) and [`near`] (`expm1`, `log1p`,
//! `root`).
//!
//! Branch cuts are the principal ones (`log`, `sqrt`, `pow`, `root`, `log1p` on
//! the negative real axis, from `-1` for `log1p`; `atan`, `asinh` on the
//! imaginary axis beyond `±i`; `asin`, `acos`, `atanh` on the real axis beyond
//! `±1`, `acosh` below 1). A ball straddling a cut gets an enclosure of both
//! sides: `arg` is then `[0 ± π]`, and `sqrt`'s imaginary part `[0 ± t]`. A
//! point exactly on a cut takes the side a `+0` in the zero part selects in C99
//! (Annex G) and MPC \[C99\] \[MPC\], after Kahan \[Kah87\]: `arg` in `(-π, π]`
//! for `log` and `sqrt`, the side of `Re z > 0` for `atan` and `asinh`, of `Im
//! z > 0` for the others.

use std::cmp::Ordering;
use std::fmt;

use super::{Ball, Float, Mag};
use crate::{Complex, Integer};

mod inverse;
mod near;
mod trig;

#[cfg(test)]
mod fn_tests;
#[cfg(test)]
mod tests;

/// From this precision a complex product takes 3 real products instead of
/// 4 (the third costs additions and a wider radius from the cancellation in
/// `(a+b)(c+d) - ac - bd`, so it only pays once products dominate).
const MUL3_PREC: u64 = 4096;

/// A complex number enclosed as a rectangle `re + im i` of real balls.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ComplexBall {
    re: Ball,
    im: Ball,
}

fn one() -> Ball {
    Ball::from_integer(&Integer::ONE)
}

impl ComplexBall {
    /// The exact zero.
    pub const ZERO: ComplexBall = ComplexBall {
        re: Ball::ZERO,
        im: Ball::ZERO,
    };

    /// Indeterminate in both parts.
    pub const INDETERMINATE: ComplexBall = ComplexBall {
        re: Ball::INDETERMINATE,
        im: Ball::INDETERMINATE,
    };

    #[must_use]
    pub fn new(re: Ball, im: Ball) -> ComplexBall {
        ComplexBall { re, im }
    }

    /// The real ball `x + 0i`.
    #[must_use]
    pub fn from_real(x: Ball) -> ComplexBall {
        ComplexBall {
            re: x,
            im: Ball::ZERO,
        }
    }

    /// A ball containing the Gaussian rational `z`, each part as
    /// [`Ball::from_rational`] at `prec` bits.
    #[must_use]
    pub fn from_complex(z: &Complex, prec: u64) -> ComplexBall {
        ComplexBall {
            re: Ball::from_rational(z.re(), prec),
            im: Ball::from_rational(z.im(), prec),
        }
    }

    #[must_use]
    pub fn re(&self) -> &Ball {
        &self.re
    }

    #[must_use]
    pub fn im(&self) -> &Ball {
        &self.im
    }

    /// Whether both radii are 0.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.re.is_exact() && self.im.is_exact()
    }

    /// Whether both parts are finite balls.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }

    /// Whether the imaginary part is exactly 0.
    #[must_use]
    pub fn is_real(&self) -> bool {
        self.im.is_exact() && self.im.mid().is_zero()
    }

    /// Whether the rectangle contains the Gaussian rational `z`.
    #[must_use]
    pub fn contains_complex(&self, z: &Complex) -> bool {
        self.re.contains_rational(z.re()) && self.im.contains_rational(z.im())
    }

    /// Whether the rectangles share a point.
    #[must_use]
    pub fn overlaps(&self, other: &ComplexBall) -> bool {
        self.re.overlaps(&other.re) && self.im.overlaps(&other.im)
    }

    /// Whether the rectangle may contain 0.
    #[must_use]
    pub fn may_be_zero(&self) -> bool {
        self.re.contains(&Float::ZERO) && self.im.contains(&Float::ZERO)
    }

    #[must_use]
    pub fn neg(&self) -> ComplexBall {
        ComplexBall {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }

    #[must_use]
    pub fn conj(&self) -> ComplexBall {
        ComplexBall {
            re: self.re.clone(),
            im: self.im.neg(),
        }
    }

    /// `i · self`.
    fn mul_i(&self) -> ComplexBall {
        ComplexBall {
            re: self.im.neg(),
            im: self.re.clone(),
        }
    }

    #[must_use]
    pub fn add(&self, rhs: &ComplexBall, prec: u64) -> ComplexBall {
        ComplexBall {
            re: self.re.add(&rhs.re, prec),
            im: self.im.add(&rhs.im, prec),
        }
    }

    #[must_use]
    pub fn sub(&self, rhs: &ComplexBall, prec: u64) -> ComplexBall {
        ComplexBall {
            re: self.re.sub(&rhs.re, prec),
            im: self.im.sub(&rhs.im, prec),
        }
    }

    /// `self × rhs`: 4 real products, 3 from 4096 bits; a real
    /// factor multiplies the parts directly.
    #[must_use]
    pub fn mul(&self, rhs: &ComplexBall, prec: u64) -> ComplexBall {
        if rhs.is_real() {
            return ComplexBall {
                re: self.re.mul(&rhs.re, prec),
                im: self.im.mul(&rhs.re, prec),
            };
        }
        if self.is_real() {
            return rhs.mul(self, prec);
        }
        let (a, b, c, d) = (&self.re, &self.im, &rhs.re, &rhs.im);
        let wp = prec + 4;
        if prec < MUL3_PREC {
            let re = a.mul(c, wp).sub(&b.mul(d, wp), prec);
            let im = a.mul(d, wp).add(&b.mul(c, wp), prec);
            return ComplexBall { re, im };
        }
        let (ac, bd) = (a.mul(c, wp), b.mul(d, wp));
        let cross = a.add(b, wp).mul(&c.add(d, wp), wp);
        ComplexBall {
            re: ac.sub(&bd, prec),
            im: cross.sub(&ac, wp).sub(&bd, prec),
        }
    }

    /// `self / rhs`; indeterminate if `rhs` may be 0.
    #[must_use]
    pub fn div(&self, rhs: &ComplexBall, prec: u64) -> ComplexBall {
        if rhs.is_real() {
            return ComplexBall {
                re: self.re.div(&rhs.re, prec),
                im: self.im.div(&rhs.re, prec),
            };
        }
        let wp = prec + 8;
        let (a, b, c, d) = (&self.re, &self.im, &rhs.re, &rhs.im);
        let den = c.mul(c, wp).add(&d.mul(d, wp), wp);
        let re = a.mul(c, wp).add(&b.mul(d, wp), wp);
        let im = b.mul(c, wp).sub(&a.mul(d, wp), wp);
        ComplexBall {
            re: re.div(&den, prec),
            im: im.div(&den, prec),
        }
    }

    /// `1 / self`.
    #[must_use]
    pub fn inv(&self, prec: u64) -> ComplexBall {
        ComplexBall::from_real(one()).div(self, prec)
    }

    /// `|self|`, a real ball.
    #[must_use]
    pub fn abs(&self, prec: u64) -> Ball {
        if self.is_real() {
            return self.re.abs().round_to_accuracy(prec);
        }
        let wp = prec + 8;
        let n = self
            .re
            .mul(&self.re, wp)
            .add(&self.im.mul(&self.im, wp), wp);
        n.sqrt(prec)
    }

    /// The principal argument `atan2(im, re)` in `(-π, π]`, a real ball:
    /// `[0 ± π]` across the negative real axis or around 0.
    #[must_use]
    pub fn arg(&self, prec: u64) -> Ball {
        self.im.atan2(&self.re, prec)
    }

    /// `exp(self)`.
    #[must_use]
    pub fn exp(&self, prec: u64) -> ComplexBall {
        let wp = prec + 8;
        let ea = self.re.exp(wp);
        if self.is_real() {
            return ComplexBall::from_real(ea.round_to_accuracy(prec));
        }
        let (s, c) = self.im.sin_cos(wp);
        ComplexBall {
            re: ea.mul(&c, prec),
            im: ea.mul(&s, prec),
        }
    }

    /// The principal `log(self)`; indeterminate if `self` may be 0.
    #[must_use]
    pub fn log(&self, prec: u64) -> ComplexBall {
        if self.may_be_zero() {
            return ComplexBall::INDETERMINATE;
        }
        if self.is_real() && self.re.sign() == Some(Ordering::Greater) {
            return ComplexBall::from_real(self.re.log(prec));
        }
        let wp = prec + 8;
        let (a, b) = (&self.re, &self.im);
        let n = a.mul(a, wp).add(&b.mul(b, wp), wp);
        // log|z| = log(n) / 2; near |z| = 1, log1p(n - 1) / 2 with n - 1 as
        // (a - 1)(a + 1) + b^2, which does not cancel.
        let near_one = n.sub(&one(), 64).abs();
        let half = Mag::ONE.mul_2exp(&-Integer::ONE);
        let log_abs = if upper(&near_one) <= half {
            let t = a
                .sub(&one(), wp)
                .mul(&a.add(&one(), wp), wp)
                .add(&b.mul(b, wp), wp);
            t.log1p(prec + 1)
        } else {
            n.log(prec + 1)
        };
        ComplexBall {
            re: log_abs.mul_2exp(&-Integer::ONE),
            im: self.arg(prec),
        }
    }

    /// The principal square root (real part `>= 0`); across the negative real
    /// axis the imaginary part covers both sides.
    #[must_use]
    pub fn sqrt(&self, prec: u64) -> ComplexBall {
        let wp = prec + 8;
        let (a, b) = (&self.re, &self.im);
        if self.is_real() && a.sign() != Some(Ordering::Less) {
            return ComplexBall::from_real(a.sqrt(prec));
        }
        if b.is_exact() && b.mid().is_zero() && a.sign() == Some(Ordering::Less) {
            // sqrt of a negative real: i √|a| (the upper side of the cut).
            return ComplexBall {
                re: Ball::ZERO,
                im: a.neg().sqrt(prec),
            };
        }
        let t = self
            .abs(wp)
            .add(&a.abs(), wp)
            .mul_2exp(&-Integer::ONE)
            .sqrt(wp);
        let half_over_t = |x: &Ball| x.div(&t.mul_2exp(&Integer::ONE), prec);
        match a.sign() {
            Some(Ordering::Greater | Ordering::Equal) => ComplexBall {
                re: t.round_to_accuracy(prec),
                im: half_over_t(b),
            },
            Some(Ordering::Less) => match b.sign() {
                Some(Ordering::Greater) => ComplexBall {
                    re: half_over_t(&b.abs()),
                    im: t.round_to_accuracy(prec),
                },
                Some(Ordering::Less) => ComplexBall {
                    re: half_over_t(&b.abs()),
                    im: t.neg().round_to_accuracy(prec),
                },
                // Across the cut: both signs of the imaginary part.
                _ => ComplexBall {
                    re: half_over_t(&b.abs()),
                    im: Ball::new(Float::ZERO, upper(&t)),
                },
            },
            // a of either sign: with b away from 0, |z| >= |b| and the
            // direct √((|z| ± a) / 2) do not cancel; otherwise the hull, real
            // part in [0, t] and imaginary in [-t, t].
            None => match b.sign() {
                Some(side @ (Ordering::Greater | Ordering::Less)) => {
                    let m = self.abs(wp);
                    let re = m.add(a, wp).mul_2exp(&-Integer::ONE).sqrt(prec);
                    let im = m.sub(a, wp).mul_2exp(&-Integer::ONE).sqrt(prec);
                    ComplexBall {
                        re,
                        im: if side == Ordering::Less { im.neg() } else { im },
                    }
                }
                _ => ComplexBall {
                    re: Ball::new(Float::ZERO, upper(&t)),
                    im: Ball::new(Float::ZERO, upper(&t)),
                },
            },
        }
    }

    /// `self^n` for an integer `n` by squaring (`self^0 = 1`).
    #[must_use]
    pub fn pow_int(&self, n: i64, prec: u64) -> ComplexBall {
        let k = n.unsigned_abs();
        let wp = prec + 2 * u64::from(u64::BITS - k.leading_zeros()) + 8;
        let mut r = ComplexBall::from_real(one());
        for i in (0..u64::BITS - k.leading_zeros()).rev() {
            r = r.mul(&r, wp);
            if (k >> i) & 1 == 1 {
                r = r.mul(self, wp);
            }
        }
        let r = if n < 0 { r.inv(wp) } else { r };
        r.round(prec)
    }

    /// The principal `self^w = exp(w log self)`; integer `w` by squaring.
    #[must_use]
    pub fn pow(&self, w: &ComplexBall, prec: u64) -> ComplexBall {
        if w.is_real() && w.re.is_exact() {
            let m = w.re.mid();
            if *m.exponent() <= Integer::from(62i64) {
                let k = crate::ball::decimal::to_integer(m, false);
                if Float::from(&k) == *m {
                    if let Some(k) = k.to_i64() {
                        return self.pow_int(k, prec);
                    }
                }
            }
        }
        if self.is_exact() && self.re.mid().is_zero() && self.im.mid().is_zero() {
            return match w.re.sign() {
                Some(Ordering::Greater) => ComplexBall::ZERO,
                _ => ComplexBall::INDETERMINATE,
            };
        }
        let wp = prec + 16;
        w.mul(&self.log(wp), wp).exp(prec)
    }

    /// `sin(self)`.
    #[must_use]
    pub fn sin(&self, prec: u64) -> ComplexBall {
        let wp = prec + 8;
        if self.is_real() {
            return ComplexBall::from_real(self.re.sin(prec));
        }
        let (s, c) = self.re.sin_cos(wp);
        ComplexBall {
            re: s.mul(&self.im.cosh(wp), prec),
            im: c.mul(&self.im.sinh(wp), prec),
        }
    }

    /// `cos(self)`.
    #[must_use]
    pub fn cos(&self, prec: u64) -> ComplexBall {
        let wp = prec + 8;
        if self.is_real() {
            return ComplexBall::from_real(self.re.cos(prec));
        }
        let (s, c) = self.re.sin_cos(wp);
        ComplexBall {
            re: c.mul(&self.im.cosh(wp), prec),
            im: s.mul(&self.im.sinh(wp), prec).neg(),
        }
    }

    /// The principal `atan(self)`, cuts on the imaginary axis beyond `±i`.
    #[must_use]
    pub fn atan(&self, prec: u64) -> ComplexBall {
        if self.is_real() {
            return ComplexBall::from_real(self.re.atan(prec));
        }
        // Near 0 the two logarithms are ±iz + O(z^2) and their sum cancels
        // to about |z| relative: |log2 |z|| more bits.
        let size = Mag::from_float(self.re.mid()).max(Mag::from_float(self.im.mid()));
        let extra = if size.is_zero() {
            0
        } else {
            size.exponent()
                .to_i64()
                .map_or(0, |e| u64::try_from(-e).unwrap_or(0))
        };
        let wp = prec + 16 + extra;
        let iz = self.mul_i();
        let unit = ComplexBall::from_real(one());
        let d = unit
            .add(&iz, wp)
            .log(wp)
            .sub(&unit.sub(&iz, wp).log(wp), wp);
        // d / (2i) = -i d / 2.
        let h = d.mul_i().neg();
        let mut re = h.re.mul_2exp(&-Integer::ONE).round_to_accuracy(prec);
        // On a cut itself (real part exactly 0, |im| > 1) the value continuous
        // from the right half-plane, as C99's catan and MPC: real part +π/2.
        // The formula gives that above +i and -π/2 below -i.
        let beyond = |b: &Ball| b.abs().cmp(&one()) == Some(Ordering::Greater);
        if self.re.is_exact() && self.re.mid().is_zero() && beyond(&self.im) {
            re = Ball::const_pi(prec + 8)
                .mul_2exp(&-Integer::ONE)
                .round_to_accuracy(prec);
        }
        ComplexBall {
            re,
            im: h.im.mul_2exp(&-Integer::ONE).round_to_accuracy(prec),
        }
    }

    /// Both parts rounded to `prec` bits.
    fn round(&self, prec: u64) -> ComplexBall {
        ComplexBall {
            re: self.re.round_to_accuracy(prec),
            im: self.im.round_to_accuracy(prec),
        }
    }
}

/// An upper bound of `|x|` over the ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// Whether the ball is the exact 0.
fn is_zero(b: &Ball) -> bool {
    b.is_exact() && b.mid().is_zero()
}

/// The ball `[0, u]`.
fn hull0(u: &Mag) -> Ball {
    match u.to_float() {
        Some(f) => Ball::new(f.mul_2exp(&-Integer::ONE), u.mul_2exp(&-Integer::ONE)),
        None => Ball::INDETERMINATE,
    }
}

/// `x^2` as a ball of values `>= 0`: the product `x × x` when it stays
/// above 0, else `[0, (|m| + r)^2]` (the product of a ball near 0 with
/// itself reaches below 0).
fn square(x: &Ball, prec: u64) -> Ball {
    let p = x.mul(x, prec);
    if p.is_exact() || p.sign() == Some(Ordering::Greater) {
        p
    } else {
        let u = upper(x);
        hull0(&u.mul(&u))
    }
}

/// A ball whose values are known to be `>= 0`: itself, or `[0, |m| + r]`
/// if its radius reaches below 0 (a quotient of such terms, whose ball is
/// symmetric about the midpoint).
fn clamp0(x: Ball) -> Ball {
    if x.is_exact() || x.sign() == Some(Ordering::Greater) {
        x
    } else {
        hull0(&upper(&x))
    }
}

/// `√x` for a ball whose values are known to be `>= 0` but whose radius
/// may reach below 0 (a sum of such terms): `[0, √(|m| + r)]` then.
fn sqrt_nonneg(x: &Ball, prec: u64) -> Ball {
    if x.is_exact() || x.sign() == Some(Ordering::Greater) {
        x.sqrt(prec)
    } else {
        hull0(&upper(x).sqrt())
    }
}

/// `|x|` as a ball of values `>= 0` and the sign of `x`; for a ball across
/// 0 the hull `[0, |m| + r]` and `None`.
fn fold(x: &Ball) -> (Ball, Option<Ordering>) {
    match x.sign() {
        Some(Ordering::Less) => (x.neg(), Some(Ordering::Less)),
        None => (hull0(&upper(x)), None),
        sign => (x.clone(), sign),
    }
}

/// An odd function's value from its value `v` at `|x|`, for `x` of the
/// given sign ([`fold`]): `v`, `-v`, or `[0 ± |v|]` across 0. An exact 0
/// counts as `+0`.
fn unfold(v: Ball, sign: Option<Ordering>) -> Ball {
    match sign {
        Some(Ordering::Less) => v.neg(),
        Some(_) => v,
        None => Ball::new(Float::ZERO, upper(&v)),
    }
}

/// `x / 2`, exactly.
fn half(x: &Ball) -> Ball {
    x.mul_2exp(&-Integer::ONE)
}

/// A lower bound of `|x|` over the ball (0 if it contains 0).
fn lower(b: &Ball) -> Mag {
    Mag::from_float_lower(b.mid()).sub_lower(b.rad())
}

impl ComplexBall {
    /// `f` at the exact midpoint, widened in each part by `d √(ρ_re^2 +
    /// ρ_im^2)` for `d = deriv()` an upper bound of `|f'|` on the
    /// rectangle: `f` holomorphic there, the mean value inequality along the
    /// segment from the midpoint gives `|f(w) - f(m)| <= d |w - m|`. An exact
    /// ball is `f` itself.
    fn at_midpoint(
        &self,
        prec: u64,
        f: impl FnOnce(&ComplexBall, u64) -> ComplexBall,
        deriv: impl FnOnce() -> Mag,
    ) -> ComplexBall {
        if self.is_exact() {
            return f(self, prec);
        }
        if !self.is_finite() {
            return ComplexBall::INDETERMINATE;
        }
        let mid = ComplexBall {
            re: Ball::exact(self.re.mid().clone()),
            im: Ball::exact(self.im.mid().clone()),
        };
        let value = f(&mid, prec);
        let (rre, rim) = (self.re.rad(), self.im.rad());
        let err = deriv().mul(&rre.mul(rre).add(&rim.mul(rim)).sqrt());
        ComplexBall {
            re: value.re.add_error(&err),
            im: value.im.add_error(&err),
        }
    }
}

/// `[re] + [im] i` with the `Debug` forms of both parts.
impl fmt::Debug for ComplexBall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} + {:?} i", self.re, self.im)
    }
}

/// The parts in decimal as [`Ball`]'s `Display`, `re + im*I`.
impl fmt::Display for ComplexBall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.re, f)?;
        f.write_str(" + ")?;
        fmt::Display::fmt(&self.im, f)?;
        f.write_str("*I")
    }
}
