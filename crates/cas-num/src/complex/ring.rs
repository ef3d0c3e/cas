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
/// nearest integer, so `norm(r) <= norm(b) / 2`. Normal forms are the
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
        let exact = a / b;
        let q = Complex::new(
            Rational::from(exact.re().round()),
            Rational::from(exact.im().round()),
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
