//! Arbitrary-precision rational numbers.

mod ops;
mod radix;

#[cfg(test)]
mod tests;

use std::cmp::Ordering;

use crate::{Integer, Natural, OutOfRangeError};

/// An arbitrary-precision rational number.
///
/// Always canonical: the denominator is positive, numerator and denominator
/// are coprime, and zero is `0/1`. Equality and hashing are therefore by
/// value. An integer value is still a `Rational` with denominator 1; the
/// expression layer decides whether to store it as an integer.
///
/// Arithmetic follows `std` conventions: operators are implemented for every
/// owned/borrowed combination and division by zero panics
/// ([`Rational::checked_div`] returns `None` instead). Sums and products
/// cancel common factors before multiplying (Henrici's method, Knuth §4.5.1),
/// so intermediate values stay small.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Rational {
    num: Integer,
    /// Positive and coprime to `num`.
    den: Integer,
}

impl Rational {
    /// Zero.
    pub const ZERO: Rational = Rational {
        num: Integer::ZERO,
        den: Integer::ONE,
    };

    /// One.
    pub const ONE: Rational = Rational {
        num: Integer::ONE,
        den: Integer::ONE,
    };

    /// `num / den`, reduced.
    ///
    /// # Panics
    ///
    /// If `den` is zero.
    #[must_use]
    pub fn new(num: Integer, den: Integer) -> Rational {
        Rational::checked_new(num, den).expect("attempt to divide by zero")
    }

    /// `num / den`, reduced, or `None` if `den` is zero.
    #[must_use]
    pub fn checked_new(num: Integer, den: Integer) -> Option<Rational> {
        if den.is_zero() {
            return None;
        }
        let g = num.gcd(&den);
        let (num, den) = if g.is_one() {
            (num, den)
        } else {
            (num.div_exact_unchecked(&g), den.div_exact_unchecked(&g))
        };
        Some(if den.is_negative() {
            Rational::from_canonical(-num, -den)
        } else {
            Rational::from_canonical(num, den)
        })
    }

    /// Wraps parts known to be canonical.
    fn from_canonical(num: Integer, den: Integer) -> Rational {
        let r = Rational { num, den };
        debug_assert!(r.is_canonical(), "non-canonical Rational {r:?}");
        r
    }

    /// The numerator; it carries the sign.
    #[must_use]
    pub fn numer(&self) -> &Integer {
        &self.num
    }

    /// The denominator, which is positive.
    #[must_use]
    pub fn denom(&self) -> &Integer {
        &self.den
    }

    /// Numerator and denominator.
    #[must_use]
    pub fn into_parts(self) -> (Integer, Integer) {
        (self.num, self.den)
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    #[must_use]
    pub fn is_one(&self) -> bool {
        self.num.is_one() && self.den.is_one()
    }

    /// Whether the denominator is 1.
    #[must_use]
    pub fn is_integer(&self) -> bool {
        self.den.is_one()
    }

    /// `Less`, `Equal` or `Greater` as the value is negative, zero or positive.
    #[must_use]
    pub fn sign(&self) -> Ordering {
        self.num.sign()
    }

    #[must_use]
    pub fn is_negative(&self) -> bool {
        self.num.is_negative()
    }

    #[must_use]
    pub fn is_positive(&self) -> bool {
        self.num.is_positive()
    }

    /// `|self|`.
    #[must_use]
    pub fn abs(&self) -> Rational {
        Rational::from_canonical(self.num.abs(), self.den.clone())
    }

    /// `1 / self`, or `None` for zero.
    #[must_use]
    pub fn checked_recip(&self) -> Option<Rational> {
        match self.num.sign() {
            Ordering::Equal => None,
            Ordering::Greater => Some(Rational::from_canonical(self.den.clone(), self.num.clone())),
            Ordering::Less => Some(Rational::from_canonical(-&self.den, -&self.num)),
        }
    }

    /// `1 / self`.
    ///
    /// # Panics
    ///
    /// If `self` is zero.
    #[must_use]
    pub fn recip(&self) -> Rational {
        self.checked_recip().expect("attempt to divide by zero")
    }

    /// `self^e`, with `x^0 = 1`. Powers of coprime parts stay coprime, so
    /// nothing is reduced.
    ///
    /// # Panics
    ///
    /// If `self` is zero and `e` negative, or if the bit length of the
    /// result overflows `u64`.
    #[must_use]
    pub fn pow(&self, e: i64) -> Rational {
        let base = if e < 0 { &self.recip() } else { self };
        let k = e.unsigned_abs();
        Rational::from_canonical(base.num.pow(k), base.den.pow(k))
    }

    /// The largest integer `<= self`.
    #[must_use]
    pub fn floor(&self) -> Integer {
        self.num.div_rem_floor(&self.den).0
    }

    /// The smallest integer `>= self`.
    #[must_use]
    pub fn ceil(&self) -> Integer {
        self.num.div_rem_ceil(&self.den).0
    }

    /// The integer part, rounding towards zero.
    #[must_use]
    pub fn trunc(&self) -> Integer {
        self.num.div_rem_trunc(&self.den).0
    }

    /// The nearest integer, halves rounded away from zero (as `f64::round`).
    #[must_use]
    pub fn round(&self) -> Integer {
        // |self| + 1/2 = (2|num| + den) / (2 den), floored.
        let twice_den = &self.den + &self.den;
        let num = self.num.abs();
        let m = (&num + &num + &self.den).div_rem_floor(&twice_den).0;
        if self.is_negative() { -m } else { m }
    }

    /// The nearest integer, halves rounded to even (as
    /// `f64::round_ties_even`).
    #[must_use]
    pub fn round_ties_even(&self) -> Integer {
        let (q, r) = self.num.div_rem_floor(&self.den);
        // 0 <= r < den; compare the fraction r / den with 1/2.
        let up = match (&r + &r).cmp(&self.den) {
            Ordering::Less => false,
            Ordering::Greater => true,
            Ordering::Equal => !q.is_divisible_by(&Integer::from(2)),
        };
        if up { q + Integer::ONE } else { q }
    }

    /// Whether the representation is canonical (checked by the tests on every
    /// result, and in debug builds on every construction).
    pub(crate) fn is_canonical(&self) -> bool {
        self.den.is_positive()
            && if self.num.is_zero() {
                self.den.is_one()
            } else {
                self.num.gcd(&self.den).is_one()
            }
    }
}

impl Default for Rational {
    fn default() -> Self {
        Self::ZERO
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        let sign = self.sign();
        if sign != other.sign() {
            return sign.cmp(&other.sign());
        }
        if self.den == other.den {
            return self.num.cmp(&other.num);
        }
        // With e = bits(num) - bits(den), 2^(e-1) < |x| < 2^(e+1); exponents
        // two apart decide without multiplying.
        let estimate =
            |x: &Rational| i128::from(x.num.bit_length()) - i128::from(x.den.bit_length());
        let by_magnitude = match estimate(self) - estimate(other) {
            d if d >= 2 => Some(Ordering::Greater),
            d if d <= -2 => Some(Ordering::Less),
            _ => None,
        };
        match by_magnitude {
            Some(o) if sign == Ordering::Less => o.reverse(),
            Some(o) => o,
            None => (&self.num * &other.den).cmp(&(&other.num * &self.den)),
        }
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl From<Integer> for Rational {
    fn from(n: Integer) -> Rational {
        Rational::from_canonical(n, Integer::ONE)
    }
}

impl From<&Integer> for Rational {
    fn from(n: &Integer) -> Rational {
        Rational::from(n.clone())
    }
}

impl From<Natural> for Rational {
    fn from(n: Natural) -> Rational {
        Rational::from(Integer::from(n))
    }
}

macro_rules! from_primitive {
    ($($t:ty),*) => {$(
        impl From<$t> for Rational {
            fn from(x: $t) -> Rational {
                Rational::from(Integer::from(x))
            }
        }
    )*};
}

from_primitive!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl TryFrom<&Rational> for Integer {
    type Error = OutOfRangeError;

    /// Fails unless the value is an integer.
    fn try_from(x: &Rational) -> Result<Integer, OutOfRangeError> {
        if x.is_integer() {
            Ok(x.num.clone())
        } else {
            Err(OutOfRangeError(()))
        }
    }
}

impl TryFrom<Rational> for Integer {
    type Error = OutOfRangeError;

    /// Fails unless the value is an integer.
    fn try_from(x: Rational) -> Result<Integer, OutOfRangeError> {
        if x.is_integer() {
            Ok(x.num)
        } else {
            Err(OutOfRangeError(()))
        }
    }
}
