//! Ring objects: ℚ(i) as a field and ℤ\[i\] as a Euclidean domain.

use super::{Complex, quadrant, unit};
use crate::algebra::{EuclideanDomain, Field, Ring};
use crate::{Integer, Natural, Rational};

/// The field ℚ(i), with [`Complex`] elements.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct GaussianRationals;

/// The Euclidean domain ℤ\[i\] of the Gaussian integers, with [`Complex`]
/// elements whose parts are integers.
///
/// Division with remainder rounds each part of the exact quotient to the
/// nearest integer, so `norm(r) <= norm(b) / 2`: ℤ\[i\] is Euclidean
/// (\[HW08\], chapter XII). Normal forms are the
/// associates in the first quadrant: real part positive, imaginary part
/// non-negative.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct GaussianIntegers;

macro_rules! complex_ring {
    ($ring:ty) => {
        impl Ring for $ring {
            type Elem = Complex;

            fn zero(&self) -> Complex {
                Complex::ZERO
            }

            fn one(&self) -> Complex {
                Complex::ONE
            }

            fn add(&self, a: &Complex, b: &Complex) -> Complex {
                a + b
            }

            fn sub(&self, a: &Complex, b: &Complex) -> Complex {
                a - b
            }

            fn neg(&self, a: &Complex) -> Complex {
                -a
            }

            fn mul(&self, a: &Complex, b: &Complex) -> Complex {
                a * b
            }

            fn elem_integer(&self, n: &Integer) -> Complex {
                Complex::from(n.clone())
            }

            fn characteristic(&self) -> Natural {
                Natural::ZERO
            }

            fn is_zero(&self, a: &Complex) -> bool {
                a.is_zero()
            }

            fn is_one(&self, a: &Complex) -> bool {
                a.is_one()
            }
        }
    };
}

complex_ring!(GaussianRationals);
complex_ring!(GaussianIntegers);

impl EuclideanDomain for GaussianRationals {
    /// `(a / b, 0)`.
    fn div_rem(&self, a: &Complex, b: &Complex) -> (Complex, Complex) {
        (a / b, Complex::ZERO)
    }

    /// `(1, 1 / a)` for non-zero `a`.
    fn normal_form(&self, a: &Complex) -> (Complex, Complex) {
        match a.checked_recip() {
            None => (Complex::ZERO, Complex::ONE),
            Some(u) => (Complex::ONE, u),
        }
    }
}

impl Field for GaussianRationals {
    fn inv(&self, a: &Complex) -> Option<Complex> {
        a.checked_recip()
    }
}

impl EuclideanDomain for GaussianIntegers {
    /// `q` is the exact quotient with each part rounded to the nearest
    /// integer, and `r = a - q b`: each part of `a / b - q` is at most `1/2` in
    /// absolute value, so `norm(r) = norm(b) |a / b - q|^2 <= norm(b) / 2`.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    fn div_rem(&self, a: &Complex, b: &Complex) -> (Complex, Complex) {
        debug_assert!(a.is_gaussian_integer() && b.is_gaussian_integer());
        // a / b = a conj(b) / norm(b), rounded in integers: the exact
        // quotient as rationals reduced each part by a gcd as large as the
        // operands, at every step of Euclid's algorithm.
        let (ar, ai) = (a.re().numer(), a.im().numer());
        let (br, bi) = (b.re().numer(), b.im().numer());
        if let Some((re, im)) = leading_quotient([ar, ai, br, bi]) {
            let q = Complex::new(Rational::from(re), Rational::from(im));
            let r = a - &(&q * b);
            return (q, r);
        }
        let norm = &(br * br) + &(bi * bi);
        assert!(!norm.is_zero(), "division by zero");
        let re = &(ar * br) + &(ai * bi);
        let im = &(ai * br) - &(ar * bi);
        let q = Complex::new(
            Rational::from(round(&re, &norm)),
            Rational::from(round(&im, &norm)),
        );
        let r = a - &(&q * b);
        (q, r)
    }

    /// The associate `a u` in the first quadrant (`re > 0`, `im >= 0`), with
    /// the unit `u = i^-k` for `a` in quadrant `k`.
    fn normal_form(&self, a: &Complex) -> (Complex, Complex) {
        if a.is_zero() {
            return (Complex::ZERO, Complex::ONE);
        }
        let u = unit(4 - quadrant(a));
        (a * &u, u)
    }
}

/// Bits of `b` that [`leading_quotient`] keeps.
const KEPT: u64 = 64;

/// The largest excess of the bits of `a` over those of `b` for
/// [`leading_quotient`].
const GAP: u64 = 32;

/// Bits of the fixed-point quotient in [`leading_quotient`].
const FRACTION: u64 = 40;

/// The parts of `a / b` rounded to the nearest integers, for the parts
/// `[re a, im a, re b, im b]`, from the leading [`KEPT`] bits of `b`: or
/// `None` when `b` is small, `a` far larger, or a part near a half.
///
/// The parts truncated by `2^k` are off by less than 1 each and `|b'| >=
/// 2^(KEPT - 1)`, so `a' / b'` is within `(1 + |a' / b'|) 2 / |b'| <
/// 2^-28` of `a / b`, `|a' / b'| < 2^(GAP + 1)`; a fractional part farther
/// than `2^-20` from `1/2` rounds the same. Each step of Euclid's algorithm
/// then costs a product by the small quotient, not products of the operands.
fn leading_quotient([ar, ai, br, bi]: [&Integer; 4]) -> Option<(Integer, Integer)> {
    let lb = br.bit_length().max(bi.bit_length());
    if lb <= 2 * KEPT || ar.bit_length().max(ai.bit_length()) > lb + GAP {
        return None;
    }
    let k = lb - KEPT;
    let top = |x: &Integer| Integer::from_parts(x.is_negative(), &x.unsigned_abs() >> k);
    let (ar, ai, br, bi) = (top(ar), top(ai), top(br), top(bi));
    let norm = &(&br * &br) + &(&bi * &bi);
    let one = Integer::from(1u64 << FRACTION);
    let half = Integer::from(1u64 << (FRACTION - 1));
    let margin = Integer::from(1u64 << (FRACTION - 20));
    let part = |x: Integer| -> Option<Integer> {
        let (q, f) = (&x * &one).div_rem_floor(&norm).0.div_rem_floor(&one);
        if (&f - &half).abs() <= margin {
            return None;
        }
        Some(if f > half { q + Integer::ONE } else { q })
    };
    let re = part(&(&ar * &br) + &(&ai * &bi))?;
    let im = part(&(&ai * &br) - &(&ar * &bi))?;
    Some((re, im))
}

/// The nearest integer to `n / d`, `d > 0`, halves rounded away from zero
/// (as [`Rational::round`]).
fn round(n: &Integer, d: &Integer) -> Integer {
    // |n| / d + 1/2 = (2 |n| + d) / (2 d), floored.
    let m = n.abs();
    let m = (&(&m + &m) + d).div_rem_floor(&(d + d)).0;
    if n.is_negative() { -m } else { m }
}
