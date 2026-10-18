//! Complex rationals ℚ(i), with the Gaussian integers ℤ\[i\] as the values whose
//! parts are integers (CAS-50).
//!
//! A [`Complex`] is a pair of canonical [`Rational`]s, so it is canonical too
//! and `Eq`/`Hash` are by value. The type allows a zero imaginary part
//! ([`Complex::is_real`]); the expression layer maps such values to the
//! rational itself, so each number has one form there. The ring objects
//! [`GaussianRationals`] (a field) and [`GaussianIntegers`] (a Euclidean
//! domain) live in [`crate::algebra`].

use std::cmp::Ordering;

use crate::{Integer, Rational};

mod fmt;
mod ops;
mod ring;

#[cfg(test)]
mod tests;

pub use ring::{GaussianIntegers, GaussianRationals};

/// A complex number with rational parts, `re + im·i`.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct Complex {
    re: Rational,
    im: Rational,
}

impl Complex {
    /// Zero.
    pub const ZERO: Complex = Complex {
        re: Rational::ZERO,
        im: Rational::ZERO,
    };

    /// One.
    pub const ONE: Complex = Complex {
        re: Rational::ONE,
        im: Rational::ZERO,
    };

    /// The imaginary unit.
    pub const I: Complex = Complex {
        re: Rational::ZERO,
        im: Rational::ONE,
    };

    #[must_use]
    pub fn new(re: Rational, im: Rational) -> Complex {
        Complex { re, im }
    }

    /// The real part.
    #[must_use]
    pub fn re(&self) -> &Rational {
        &self.re
    }

    /// The imaginary part.
    #[must_use]
    pub fn im(&self) -> &Rational {
        &self.im
    }

    /// `(re, im)`.
    #[must_use]
    pub fn into_parts(self) -> (Rational, Rational) {
        (self.re, self.im)
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }

    #[must_use]
    pub fn is_one(&self) -> bool {
        self.re.is_one() && self.im.is_zero()
    }

    /// Whether the imaginary part is zero.
    #[must_use]
    pub fn is_real(&self) -> bool {
        self.im.is_zero()
    }

    /// Whether both parts are integers: an element of ℤ\[i\].
    #[must_use]
    pub fn is_gaussian_integer(&self) -> bool {
        self.re.is_integer() && self.im.is_integer()
    }

    /// The complex conjugate `re - im·i`.
    #[must_use]
    pub fn conj(&self) -> Complex {
        Complex {
            re: self.re.clone(),
            im: -&self.im,
        }
    }

    /// The field norm `re² + im²` (the squared absolute value).
    #[must_use]
    pub fn norm(&self) -> Rational {
        &(&self.re * &self.re) + &(&self.im * &self.im)
    }

    /// `1 / self`, or `None` for zero: `conj(z) / norm(z)`.
    #[must_use]
    pub fn checked_recip(&self) -> Option<Complex> {
        let n = self.norm();
        let re = self.re.checked_div(&n)?;
        let im = -&self.im.checked_div(&n)?;
        Some(Complex { re, im })
    }

    /// `1 / self`.
    ///
    /// # Panics
    ///
    /// If `self` is zero.
    #[must_use]
    pub fn recip(&self) -> Complex {
        self.checked_recip().expect("attempt to divide by zero")
    }

    /// `self^e`, by binary exponentiation; a negative `e` raises the
    /// reciprocal, and `z^0 = 1` (also for zero).
    ///
    /// # Panics
    ///
    /// If `self` is zero and `e` is negative.
    #[must_use]
    pub fn pow(&self, e: i64) -> Complex {
        let base = if e < 0 { self.recip() } else { self.clone() };
        let mut r = Complex::ONE;
        let k = e.unsigned_abs();
        for i in (0..u64::BITS - k.leading_zeros()).rev() {
            r = &r * &r;
            if (k >> i) & 1 == 1 {
                r = &r * &base;
            }
        }
        r
    }

    /// Whether the canonical-form invariant holds (for tests).
    #[cfg(test)]
    pub(crate) fn is_canonical(&self) -> bool {
        self.re.is_canonical() && self.im.is_canonical()
    }
}

impl From<Rational> for Complex {
    fn from(re: Rational) -> Complex {
        Complex {
            re,
            im: Rational::ZERO,
        }
    }
}

impl From<Integer> for Complex {
    fn from(re: Integer) -> Complex {
        Complex::from(Rational::from(re))
    }
}

impl From<i64> for Complex {
    fn from(re: i64) -> Complex {
        Complex::from(Integer::from(re))
    }
}

/// The units of ℤ\[i\] and their order around the circle: `1, i, -1, -i`.
pub(crate) fn unit(k: u8) -> Complex {
    match k % 4 {
        0 => Complex::ONE,
        1 => Complex::I,
        2 => -Complex::ONE,
        _ => -Complex::I,
    }
}

/// The quadrant of a non-zero value, counting from the closed positive real
/// half-axis: 0 for `re > 0, im >= 0`, 1 for `re <= 0, im > 0`, and so on.
pub(crate) fn quadrant(z: &Complex) -> u8 {
    debug_assert!(!z.is_zero());
    match (z.re.sign(), z.im.sign()) {
        (Ordering::Greater, Ordering::Equal | Ordering::Greater) => 0,
        (Ordering::Less | Ordering::Equal, Ordering::Greater) => 1,
        (Ordering::Less, Ordering::Less | Ordering::Equal) => 2,
        _ => 3,
    }
}
