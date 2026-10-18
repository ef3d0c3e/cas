//! Signed arbitrary-precision integers with a small-value fast path.

mod convert;
mod div;
mod gcd;
mod ops;
mod pow;
mod radix;

#[cfg(test)]
mod tests;

use std::borrow::Cow;
use std::cmp::Ordering;
use std::fmt;

use crate::Natural;

pub use div::Rounding;

/// An arbitrary-precision integer.
///
/// Values that fit an `i64` are stored inline and never allocate; larger values
/// are boxed sign and magnitude (one allocation up to 256 bits, layout E in the
/// numerics design). The representation is canonical (a value is boxed iff it
/// does not fit an `i64`), so equality and hashing are by value.
///
/// Arithmetic follows `std` conventions: operators are implemented for every
/// owned/borrowed combination; `/` and `%` truncate towards zero and panic on
/// a zero divisor. Other roundings are explicit: [`Integer::div_rem_floor`],
/// [`Integer::div_rem_ceil`], [`Integer::div_rem_euclid`].
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Integer(Repr);

#[derive(Clone, PartialEq, Eq, Hash)]
enum Repr {
    Small(i64),
    /// Never a value that fits an `i64`.
    Big(Box<Big>),
}

#[derive(Clone, PartialEq, Eq, Hash)]
struct Big {
    negative: bool,
    /// Non-zero (zero fits an `i64`).
    mag: Natural,
}

impl Integer {
    /// Zero.
    pub const ZERO: Integer = Integer(Repr::Small(0));

    /// One.
    pub const ONE: Integer = Integer(Repr::Small(1));

    /// Zero.
    #[must_use]
    pub const fn new() -> Self {
        Self::ZERO
    }

    /// The value with the given sign and magnitude; `negative` is ignored
    /// for zero.
    #[must_use]
    pub fn from_parts(negative: bool, mag: Natural) -> Self {
        match small_value(negative, &mag) {
            Some(x) => Integer(Repr::Small(x)),
            None => Integer(Repr::Big(Box::new(Big { negative, mag }))),
        }
    }

    /// A value the caller knows does not fit an `i64`, skipping the check
    /// [`Integer::from_parts`] makes.
    fn big(negative: bool, mag: Natural) -> Self {
        debug_assert!(small_value(negative, &mag).is_none());
        Integer(Repr::Big(Box::new(Big { negative, mag })))
    }

    /// Sign and magnitude (`false` for zero).
    #[must_use]
    pub fn into_parts(self) -> (bool, Natural) {
        match self.0 {
            Repr::Small(x) => (x < 0, Natural::from(x.unsigned_abs())),
            Repr::Big(b) => (b.negative, b.mag),
        }
    }

    /// `|self|` as a [`Natural`].
    #[must_use]
    pub fn unsigned_abs(&self) -> Natural {
        self.sign_mag().1.into_owned()
    }

    /// `|self|`.
    #[must_use]
    pub fn abs(&self) -> Integer {
        if self.is_negative() {
            -self
        } else {
            self.clone()
        }
    }

    /// `Less`, `Equal` or `Greater` as the value is negative, zero or positive.
    #[must_use]
    pub fn sign(&self) -> Ordering {
        match &self.0 {
            Repr::Small(x) => x.cmp(&0),
            Repr::Big(b) if b.negative => Ordering::Less,
            Repr::Big(_) => Ordering::Greater,
        }
    }

    #[must_use]
    pub fn is_negative(&self) -> bool {
        self.sign() == Ordering::Less
    }

    #[must_use]
    pub fn is_positive(&self) -> bool {
        self.sign() == Ordering::Greater
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        matches!(self.0, Repr::Small(0))
    }

    #[must_use]
    pub fn is_one(&self) -> bool {
        matches!(self.0, Repr::Small(1))
    }

    /// Bits of `|self|`, without leading zeros (0 for zero).
    #[must_use]
    pub fn bit_length(&self) -> u64 {
        match &self.0 {
            Repr::Small(x) => u64::from(u64::BITS - x.unsigned_abs().leading_zeros()),
            Repr::Big(b) => b.mag.bit_length(),
        }
    }

    /// The value as an `i64`, if it fits (the inline representation).
    #[must_use]
    pub fn to_i64(&self) -> Option<i64> {
        match self.0 {
            Repr::Small(x) => Some(x),
            Repr::Big(_) => None,
        }
    }

    /// Sign and magnitude, borrowing the magnitude of a boxed value.
    fn sign_mag(&self) -> (bool, Cow<'_, Natural>) {
        match &self.0 {
            Repr::Small(x) => (*x < 0, Cow::Owned(Natural::from(x.unsigned_abs()))),
            Repr::Big(b) => (b.negative, Cow::Borrowed(&b.mag)),
        }
    }

    /// Whether the representation is canonical (checked by the tests on every
    /// result).
    #[cfg(test)]
    pub(crate) fn is_canonical(&self) -> bool {
        match &self.0 {
            Repr::Small(_) => true,
            Repr::Big(b) => small_value(b.negative, &b.mag).is_none(),
        }
    }
}

/// `±mag` as an `i64`, if it fits.
fn small_value(negative: bool, mag: &Natural) -> Option<i64> {
    let m = u64::try_from(mag).ok()?;
    if negative {
        0i64.checked_sub_unsigned(m)
    } else {
        i64::try_from(m).ok()
    }
}

impl Default for Integer {
    fn default() -> Self {
        Self::ZERO
    }
}

impl Ord for Integer {
    fn cmp(&self, other: &Self) -> Ordering {
        match (&self.0, &other.0) {
            (Repr::Small(a), Repr::Small(b)) => a.cmp(b),
            // A boxed value lies outside the i64 range, on the side of its sign.
            (Repr::Small(_), Repr::Big(b)) => {
                if b.negative {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            }
            (Repr::Big(a), Repr::Small(_)) => {
                if a.negative {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
            (Repr::Big(a), Repr::Big(b)) => match (a.negative, b.negative) {
                (false, false) => a.mag.cmp(&b.mag),
                (true, true) => b.mag.cmp(&a.mag),
                (false, true) => Ordering::Greater,
                (true, false) => Ordering::Less,
            },
        }
    }
}

impl PartialOrd for Integer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::LowerHex for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (negative, mag) = self.sign_mag();
        f.pad_integral(!negative, "0x", &format!("{:x}", mag.as_ref()))
    }
}

impl fmt::UpperHex for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (negative, mag) = self.sign_mag();
        f.pad_integral(!negative, "0x", &format!("{:X}", mag.as_ref()))
    }
}
