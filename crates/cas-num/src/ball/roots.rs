//! The complex roots of a polynomial with exact coefficients (CAS-211), each
//! in a ball that holds it alone: the Weierstrass iteration, also Durand and
//! Kerner's \[Ker66\], on exact midpoints, then the inclusion disks of
//! \[BH73\] to validate them, as Arb's `acb_poly_find_roots` and
//! `acb_poly_validate_roots` \[Joh17\].
//!
//! - The polynomial `p = c₀ + c₁ x + … + cₙ xⁿ`, `cₙ ≠ 0`, has Gaussian
//!   rational coefficients and simple roots: it is square-free. A root at 0
//!   comes out first, exactly.
//! - A round replaces each approximation `zᵢ` in turn by the midpoint of
//!   `zᵢ − Wᵢ`, `Wᵢ = p(zᵢ) / (cₙ ∏_{j≠i} (zᵢ − zⱼ))`, from `zᵢ = ρ wⁱ`,
//!   `w = (13 + 29i)/32` near Arb's `0.4 + 0.9i` and `ρ` a power of 2 near
//!   the geometric mean of the roots' moduli.
//! - Validation: with `Wᵢ` in ball arithmetic at the midpoints, the disks
//!   `D(zᵢ, n |Wᵢ|)` hold every root, and a disk apart from all the others
//!   exactly one \[BH73\]. When all are apart, the root in `D(zᵢ, rᵢ)` is in
//!   the ball `zᵢ ± rᵢ`, in each part.
//! - Real coefficients: a disk apart from the real line holds a root that is
//!   not real. A disk whose mirror image meets no other disk holds a real
//!   root, its conjugate being in no other disk: its ball has the exact
//!   imaginary part 0. A disk that is neither waits for more precision.
//! - The roots are validated once each `rᵢ` is at most `|zᵢ| 2⁻ᵇⁱᵗˢ`. The
//!   iteration runs at `bits + 32` bits, then at twice as many while the
//!   roots are not validated, up to `max_bits`; `None` past it, as for a
//!   polynomial with a multiple root.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the method: p of degree n, roots z, corrections w, radii r"
)]

use super::{Ball, ComplexBall, Float, Mag};
use crate::{Complex, Integer, Rational};

/// The rounds at the first precision: this times the degree, plus
/// [`ROUNDS`].
const ROUNDS_PER_DEGREE: usize = 8;

/// The rounds at each precision, from approximations already close after
/// the first.
const ROUNDS: usize = 32;

/// Rounds between two validations.
const VALIDATE_EVERY: usize = 4;

/// The roots of `Σ coeffs[k] xᵏ`, each in a ball that holds it alone, to
/// `bits` bits relatively; a real root of a polynomial with real
/// coefficients has the exact imaginary part 0. In no particular order.
/// `None` for a constant, a leading coefficient 0, or roots not validated
/// within `max_bits` bits, as those of a multiple root never are.
///
/// ```
/// use cas_num::ball::polynomial_roots;
/// use cas_num::{Complex, Integer, Rational};
///
/// // x² − 2: ±√2, real.
/// let c = |n: i64| Complex::from(Rational::from(Integer::from(n)));
/// let roots = polynomial_roots(&[c(-2), c(0), c(1)], 53, 1024).unwrap();
/// assert_eq!(roots.len(), 2);
/// assert!(roots.iter().all(|z| z.is_real()));
/// ```
#[must_use]
pub fn polynomial_roots(coeffs: &[Complex], bits: u64, max_bits: u64) -> Option<Vec<ComplexBall>> {
    let n = coeffs.len().checked_sub(1)?;
    if n == 0 || coeffs[n].is_zero() {
        return None;
    }
    if coeffs[0].is_zero() {
        // A simple root at 0, at most.
        if coeffs.get(1).is_none_or(Complex::is_zero) {
            return None;
        }
        let mut roots = vec![ComplexBall::ZERO];
        if n > 1 {
            roots.extend(polynomial_roots(&coeffs[1..], bits, max_bits)?);
        }
        return Some(roots);
    }
    let real = coeffs.iter().all(Complex::is_real);
    let max_bits = max_bits.max(bits);
    let mut prec = bits.saturating_add(32).min(max_bits);
    let mut z = start(coeffs, prec);
    let mut rounds = ROUNDS_PER_DEGREE.saturating_mul(n).saturating_add(ROUNDS);
    loop {
        let p: Vec<ComplexBall> = coeffs
            .iter()
            .map(|c| ComplexBall::from_complex(c, prec))
            .collect();
        for round in 1..=rounds {
            step(&p, &mut z, prec);
            if (round % VALIDATE_EVERY == 0 || round == rounds)
                && let Some(roots) = validated(&p, &z, prec, bits, real)
            {
                return Some(roots);
            }
        }
        if prec >= max_bits {
            return None;
        }
        prec = prec.saturating_mul(2).min(max_bits);
        rounds = ROUNDS;
    }
}

/// The exact midpoint of `z`.
fn midpoint(z: &ComplexBall) -> ComplexBall {
    ComplexBall::new(
        Ball::exact(z.re().mid().clone()),
        Ball::exact(z.im().mid().clone()),
    )
}

/// The binary exponent of `|c|`, `2^(e − 1) ≤ |c| < 2^e` near enough.
fn exponent(c: &Complex) -> i64 {
    let a = ComplexBall::from_complex(c, 64).abs(64);
    a.mid().exponent().to_i64().unwrap_or(0)
}

/// The first approximations: `ρ wⁱ` for `i < n`, `w = (13 + 29i)/32` and
/// `ρ = 2^e`, `e` the mean of the exponents of the roots' moduli that the
/// first and last nonzero coefficients give.
fn start(coeffs: &[Complex], prec: u64) -> Vec<ComplexBall> {
    let n = coeffs.len() - 1;
    let (k, low) = coeffs
        .iter()
        .enumerate()
        .find(|(_, c)| !c.is_zero())
        .expect("a nonzero coefficient");
    let spread = i64::try_from(n - k).unwrap_or(i64::MAX).max(1);
    let e = exponent(low).saturating_sub(exponent(&coeffs[n])) / spread;
    let rho = Float::from_integer_exp(&Integer::ONE, &Integer::from(e));
    let over = |m: i64| Rational::new(Integer::from(m), Integer::from(32));
    let w = ComplexBall::from_complex(&Complex::new(over(13), over(29)), prec);
    let mut power = ComplexBall::from_real(Ball::exact(rho));
    let mut z = Vec::with_capacity(n);
    for _ in 0..n {
        let next = midpoint(&power.mul(&w, prec));
        z.push(power);
        power = next;
    }
    z
}

/// `p(z)`, by Horner's rule.
fn value(p: &[ComplexBall], z: &ComplexBall, prec: u64) -> ComplexBall {
    let (last, rest) = p.split_last().expect("a coefficient");
    rest.iter()
        .rev()
        .fold(last.clone(), |v, c| v.mul(z, prec).add(c, prec))
}

/// Weierstrass's correction of `z[i]`: `p(zᵢ) / (cₙ ∏_{j≠i} (zᵢ − zⱼ))`.
fn correction(p: &[ComplexBall], z: &[ComplexBall], i: usize, prec: u64) -> ComplexBall {
    let mut d = p[p.len() - 1].clone();
    for (j, zj) in z.iter().enumerate() {
        if j != i {
            d = d.mul(&z[i].sub(zj, prec), prec);
        }
    }
    value(p, &z[i], prec).div(&d, prec)
}

/// A round of the iteration, each approximation in turn.
fn step(p: &[ComplexBall], z: &mut [ComplexBall], prec: u64) {
    for i in 0..z.len() {
        let w = correction(p, z, i, prec);
        if w.is_finite() {
            z[i] = midpoint(&z[i].sub(&w, prec));
        }
    }
}

/// An upper bound of `|z|`.
fn upper(z: &ComplexBall, prec: u64) -> Mag {
    let a = z.abs(prec);
    Mag::from_float(a.mid()).add(a.rad())
}

/// A lower bound of `|z|`, 0 when the ball holds 0.
fn lower(z: &ComplexBall, prec: u64) -> Mag {
    let a = z.abs(prec);
    Mag::from_float_lower(a.mid()).sub_lower(a.rad())
}

/// The roots in balls, when the disks `D(zᵢ, n |Wᵢ|)` are apart and to
/// `bits` bits, and with real coefficients each root is known real or not.
fn validated(
    p: &[ComplexBall],
    z: &[ComplexBall],
    prec: u64,
    bits: u64,
    real: bool,
) -> Option<Vec<ComplexBall>> {
    let n = z.len();
    let degree = Mag::from_integer(&Integer::from(u64::try_from(n).ok()?));
    let scale = Mag::power_of_two(&Integer::from(bits));
    let mut radii = Vec::with_capacity(n);
    for (i, zi) in z.iter().enumerate() {
        let w = correction(p, z, i, prec);
        if !w.is_finite() {
            return None;
        }
        let r = upper(&w, prec).mul(&degree);
        if r.mul(&scale) > lower(zi, prec) {
            return None;
        }
        radii.push(r);
    }
    let apart =
        |a: &ComplexBall, j: usize, r: &Mag| lower(&a.sub(&z[j], prec), prec) > r.add(&radii[j]);
    for i in 0..n {
        if !(i + 1..n).all(|j| apart(&z[i], j, &radii[i])) {
            return None;
        }
    }
    let mut roots = Vec::with_capacity(n);
    for (i, (zi, r)) in z.iter().zip(&radii).enumerate() {
        let ball = |m: &Float| Ball::new(m.clone(), r.clone());
        let re = ball(zi.re().mid());
        if real && Mag::from_float_lower(zi.im().mid()) <= *r {
            // Touching the real line: real when its mirror image meets no
            // other disk.
            let mirror = zi.conj();
            if !(0..n).all(|j| j == i || apart(&mirror, j, r)) {
                return None;
            }
            roots.push(ComplexBall::from_real(re));
        } else {
            roots.push(ComplexBall::new(re, ball(zi.im().mid())));
        }
    }
    Some(roots)
}

#[cfg(test)]
mod tests;
