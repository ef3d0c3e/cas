//! Airy's functions of real balls (CAS-DOC-34 §2.3, CAS-380): `Ai`, `Ai′`,
//! `Bi` and `Bi′`, all four at once, with `ζ = (2/3)|x|^(3/2)`:
//!
//! - **Series** (DLMF 9.4.1–9.4.4) while `ζ` is below the expansions'
//!   [`threshold`]: `Ai = c₁ f − c₂ g`, `Ai′ = c₁ f′ − c₂ g′`, `Bi = √3 (c₁
//!   f + c₂ g)` and `Bi′ = √3 (c₁ f′ + c₂ g′)`, with `c₁ = Ai(0) =
//!   3^(−2/3)/Γ(2/3)` and `c₂ = −Ai′(0) = 3^(−1/3)/Γ(1/3)`. The terms:
//!   `f_0 = 1` and `f_k = f_(k−1) x³/(3k (3k − 1))`; `g_0 = x` and `g_k =
//!   g_(k−1) x³/((3k + 1) 3k)`; `f′_1 = x²/2` and `f′_k = f′_(k−1)
//!   x³/((3k − 1)(3k − 3))`; `g′_0 = 1` and `g′_k = g′_(k−1) x³/(3k (3k −
//!   2))`. Once the next least divisor, `3(k + 1)(3k + 1)`, is at least
//!   `2|x|³`, every next term is at most half the one before, and what
//!   follows a term at most it. `Ai` is about `e^(−ζ)` for
//!   `x > 0` where the terms reach `e^ζ`: `2ζ log₂ e` bits more, and `ζ
//!   log₂ e` for `x < 0`.
//! - **Through Bessel functions** of orders 1/3 and 2/3 beyond it (DLMF
//!   9.6.1–9.6.9), whose expansions take the ball `ζ`. For `x > 0`, `Ai(x) =
//!   (1/π) √(x/3) K_(1/3)(ζ)`, `Ai′(x) = −(x/(π√3)) K_(2/3)(ζ)`, `Bi(x) =
//!   √(x/3) (2 I_(1/3)(ζ) + (√3/π) K_(1/3)(ζ))` and `Bi′(x) = (x/√3) (2
//!   I_(2/3)(ζ) + (√3/π) K_(2/3)(ζ))`, by `I_(−ν) = I_ν + (2/π) sin(νπ)
//!   K_ν` (DLMF 10.27.2). For `x = −z < 0`, by `J_(−ν) = cos(νπ) J_ν −
//!   sin(νπ) Y_ν`: `Ai(−z) = (√z/2) (J_(1/3)(ζ) − Y_(1/3)(ζ)/√3)`,
//!   `Ai′(−z) = (z/2) (J_(2/3)(ζ) + Y_(2/3)(ζ)/√3)`, `Bi(−z) = −(√z/2)
//!   (J_(1/3)(ζ)/√3 + Y_(1/3)(ζ))` and `Bi′(−z) = (z/2) (J_(2/3)(ζ)/√3 −
//!   Y_(2/3)(ζ))`.
//!
//! An exact input runs the Ziv loop. A ball adds its radius times a bound
//! of the derivative over it: `Ai′` for `Ai`, `x Ai` for `Ai′`, and the same
//! for `Bi`, each evaluated on the ball itself at 64 bits, whose own radius
//! only widens the bound. Past an exponent of [`MAX_EXPONENT`] in `|x|` the
//! result is indeterminate: `ζ` would need `π` to that many bits.

use std::f64::consts::LOG2_E;

use super::bessel::expansions::{Order, hankel, i_expansion, k_expansion, threshold};
use super::{exponent, one, refined, upper, working};
use crate::ball::{Ball, Mag, Round};
use crate::{Integer, Rational};

/// The largest exponent of `|x|`.
const MAX_EXPONENT: i64 = 4096;

/// One of Airy's functions.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Airy {
    /// `Ai(x)`.
    Ai,
    /// `Ai′(x)`.
    AiPrime,
    /// `Bi(x)`.
    Bi,
    /// `Bi′(x)`.
    BiPrime,
}

impl Airy {
    fn index(self) -> usize {
        match self {
            Airy::Ai => 0,
            Airy::AiPrime => 1,
            Airy::Bi => 2,
            Airy::BiPrime => 3,
        }
    }
}

fn small_integer(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

fn rational(p: i64, q: i64) -> Rational {
    Rational::new(Integer::from(p), Integer::from(q))
}

/// `|x|` at the ball's midpoint in `f64`, saturated past its range.
fn magnitude(x: &Ball) -> f64 {
    if exponent(x.mid()) > 1000 {
        return f64::INFINITY;
    }
    x.mid().to_f64().abs()
}

/// `3^(−e)` for a rational `e`.
fn power_of_three(e: &Rational, q: u64) -> Ball {
    small_integer(3)
        .log(q)
        .mul(&Ball::from_rational(e, q), q)
        .neg()
        .exp(q)
}

/// The series' `[Ai, Ai′, Bi, Bi′]` at `x`, at about `wp` bits.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    reason = "term counts, a clamped estimate, and the series' names"
)]
fn series(x: &Ball, wp: u64) -> [Ball; 4] {
    let ax = magnitude(x);
    let zeta = 2.0 / 3.0 * ax * ax.sqrt();
    let above = if x.mid().is_negative() { 1.0 } else { 2.0 } * zeta * LOG2_E;
    let q = wp + above.clamp(0.0, 1e7).ceil() as u64 + 40;
    let cube = x.mul(x, q).mul(x, q);
    let (mut f, mut g) = (one(), x.clone());
    let (mut fp, mut gp) = (x.mul(x, q).mul_2exp(&-Integer::ONE), one());
    let (mut sf, mut sg, mut sfp, mut sgp) = (f.clone(), g.clone(), fp.clone(), gp.clone());
    let mut largest = upper(&f).add(&upper(&g)).add(&upper(&fp)).add(&upper(&gp));
    let twice_cube = 2.0 * ax * ax * ax;
    let small = Mag::power_of_two(&-Integer::from(q));
    let mut k = 0u64;
    let rest = loop {
        // f_k, g_k, g′_k and f′_(k+1).
        k += 1;
        f = f.mul(&cube, q).div(&small_integer(3 * k * (3 * k - 1)), q);
        g = g.mul(&cube, q).div(&small_integer((3 * k + 1) * 3 * k), q);
        gp = gp.mul(&cube, q).div(&small_integer(3 * k * (3 * k - 2)), q);
        fp = fp.mul(&cube, q).div(&small_integer((3 * k + 2) * 3 * k), q);
        sf = sf.add(&f, q);
        sg = sg.add(&g, q);
        sfp = sfp.add(&fp, q);
        sgp = sgp.add(&gp, q);
        let size = [&f, &g, &fp, &gp]
            .iter()
            .map(|t| upper(t))
            .fold(Mag::ZERO, |a, b| if b > a { b } else { a });
        if size > largest {
            largest = size.clone();
        }
        // The next step's least divisor is g′'s, 3(k + 1)(3k + 1).
        let halving = (3 * (k + 1) * (3 * k + 1)) as f64 >= twice_cube;
        if halving && size <= largest.mul(&small) {
            break size;
        }
    };
    let [sf, sg, sfp, sgp] = [sf, sg, sfp, sgp].map(|s| s.add_error(&rest));
    let c1 = power_of_three(&rational(2, 3), q)
        .div(&Ball::from_rational(&rational(2, 3), q).gamma(q), q);
    let c2 = power_of_three(&rational(1, 3), q)
        .div(&Ball::from_rational(&rational(1, 3), q).gamma(q), q);
    let root3 = small_integer(3).sqrt(q);
    let (a, b) = (c1.mul(&sf, q), c2.mul(&sg, q));
    let (ap, bp) = (c1.mul(&sfp, q), c2.mul(&sgp, q));
    [
        a.sub(&b, q),
        ap.sub(&bp, q),
        root3.mul(&a.add(&b, q), q),
        root3.mul(&ap.add(&bp, q), q),
    ]
}

/// `[Ai, Ai′, Bi, Bi′]` at `x` through Bessel functions of `ζ`; `None` when
/// the ball does not lie on one side of 0 or the expansions do not reach
/// `wp` bits.
fn through_bessel(x: &Ball, wp: u64) -> Option<[Ball; 4]> {
    let negative = x.mid().is_negative();
    let z = if negative { x.neg() } else { x.clone() };
    let r = z.rad().to_float()?;
    let lo = z.mid().sub_round(&r, 64, Round::Floor).0;
    if lo.is_negative() || lo.is_zero() {
        return None;
    }
    // ζ to the bits of its integer part more, as its sine and exponential
    // need it absolutely.
    let size = u64::try_from(exponent(z.mid()).max(0)).unwrap_or(0);
    let q = wp + 3 * size / 2 + 32;
    let root = z.sqrt(q);
    let zeta = z
        .mul(&root, q)
        .mul(&Ball::from_rational(&rational(2, 3), q), q);
    let third = Order::new(rational(1, 3));
    let two_thirds = Order::new(rational(2, 3));
    let pi = Ball::const_pi(q);
    let root3 = small_integer(3).sqrt(q);
    if negative {
        let (j1, y1) = hankel(&third, &zeta, q)?;
        let (j2, y2) = hankel(&two_thirds, &zeta, q)?;
        let half_root = root.mul_2exp(&-Integer::ONE);
        let half_z = z.mul_2exp(&-Integer::ONE);
        let ai = half_root.mul(&j1.sub(&y1.div(&root3, q), q), q);
        let ai_prime = half_z.mul(&j2.add(&y2.div(&root3, q), q), q);
        let bi = half_root.mul(&j1.div(&root3, q).add(&y1, q), q).neg();
        let bi_prime = half_z.mul(&j2.div(&root3, q).sub(&y2, q), q);
        return Some([ai, ai_prime, bi, bi_prime]);
    }
    let k1 = k_expansion(&third, &zeta, q)?;
    let k2 = k_expansion(&two_thirds, &zeta, q)?;
    let i1 = i_expansion(&third, &zeta, q)?;
    let i2 = i_expansion(&two_thirds, &zeta, q)?;
    let root_third = z.div(&small_integer(3), q).sqrt(q);
    let z_over_root3 = z.div(&root3, q);
    let ai = root_third.mul(&k1, q).div(&pi, q);
    let ai_prime = z_over_root3.mul(&k2, q).div(&pi, q).neg();
    let mixed = |i: &Ball, k: &Ball| {
        i.mul_2exp(&Integer::ONE)
            .add(&root3.div(&pi, q).mul(k, q), q)
    };
    let bi = root_third.mul(&mixed(&i1, &k1), q);
    let bi_prime = z_over_root3.mul(&mixed(&i2, &k2), q);
    Some([ai, ai_prime, bi, bi_prime])
}

/// An upper bound of `|f′|` over the ball, `f` the function `which`, from
/// the ball's ends and 0 alone. Below 0, by the energy `E = y′² − x y²` of
/// `y″ = x y`, which falls as `x` grows (`E′ = −y²`): on `[a, 0]`, `y′² ≤
/// E(a)` and `|x| y² ≤ E(a)`, so `|x y| ≤ √(|a| E(a))`. Above 0, by
/// monotony: `Ai` and `|Ai′|` fall there, `Bi` and `Bi′` rise. A wide ball
/// costs no more than its ends; the series over a ball that met 0 and
/// reached `|x| = 7·10⁵` took minutes (CAS-420). `None` past
/// [`MAX_EXPONENT`].
fn slope(x: &Ball, which: Airy) -> Option<Mag> {
    let r = x.rad().to_float()?;
    let lo = x.mid().sub_round(&r, 64, Round::Floor).0;
    let hi = x.mid().add_round(&r, 64, Round::Ceil).0;
    if exponent(&lo) > MAX_EXPONENT || exponent(&hi) > MAX_EXPONENT {
        return None;
    }
    let (y, dy) = match which {
        Airy::Ai | Airy::AiPrime => (0, 1),
        Airy::Bi | Airy::BiPrime => (2, 3),
    };
    let mut bound = Mag::ZERO;
    if lo.is_negative() {
        let a = Ball::exact(lo.clone());
        let v = airy_at(&a, 64);
        let square = |b: &Ball| b.mul(b, 64);
        let energy = upper(&square(&v[dy]).sub(&a.mul(&square(&v[y]), 64), 64));
        let part = match which {
            Airy::Ai | Airy::Bi => energy.sqrt(),
            Airy::AiPrime | Airy::BiPrime => Mag::from_float(&lo).mul(&energy).sqrt(),
        };
        bound = bound.max(part);
    }
    if !hi.is_negative() && !hi.is_zero() {
        let p = if lo.is_negative() {
            Ball::ZERO
        } else {
            Ball::exact(lo)
        };
        let (low, high) = (airy_at(&p, 64), airy_at(&Ball::exact(hi.clone()), 64));
        let far = Mag::from_float(&hi);
        let part = match which {
            Airy::Ai => upper(&low[1]),
            Airy::AiPrime => far.mul(&upper(&low[0])),
            Airy::Bi => upper(&high[3]),
            Airy::BiPrime => far.mul(&upper(&high[2])),
        };
        bound = bound.max(part);
    }
    Some(bound)
}

/// `[Ai, Ai′, Bi, Bi′]` at `x`, at about `wp` bits.
fn airy_at(x: &Ball, wp: u64) -> [Ball; 4] {
    let ax = magnitude(x);
    let zeta = 2.0 / 3.0 * ax * ax.sqrt();
    if zeta >= threshold(wp + 32)
        && let Some(values) = through_bessel(x, wp)
    {
        return values;
    }
    series(x, wp)
}

impl Ball {
    /// One of Airy's functions, `Ai`, `Ai′`, `Bi` or `Bi′`, for every `x`
    /// in the ball, at `prec` bits; indeterminate past an exponent of 4096
    /// in `|x|`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn airy(&self, which: Airy, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::airy: precision 0");
        if !self.is_finite() || exponent(self.mid()) > MAX_EXPONENT {
            return Ball::INDETERMINATE;
        }
        let i = which.index();
        if self.is_exact() {
            return refined(self, prec, |wp| airy_at(self, wp)[i].clone()).round_to_accuracy(prec);
        }
        // The derivative over the ball: Ai′, x Ai, Bi′, x Bi.
        let Some(slope) = slope(self, which) else {
            return Ball::INDETERMINATE;
        };
        let middle = Ball::exact(self.mid().clone());
        let value = airy_at(&middle, working(prec))[i].clone();
        value
            .add_error(&self.rad().mul(&slope))
            .round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
