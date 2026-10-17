//! Conversion to and from strings in radix 2 to 36.

use std::fmt;
use std::str::FromStr;

use super::Rational;
use crate::error::{ParseErrorKind, ParseNumberError};
use crate::natural::check_radix;
use crate::{Integer, Natural};

impl Rational {
    /// `num/den` in `radix`, or just `num` for an integer, lowercase, with a
    /// leading `-` for negative values.
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    #[must_use]
    pub fn to_string_radix(&self, radix: u32) -> String {
        let num = self.num.to_string_radix(radix);
        if self.is_integer() {
            num
        } else {
            format!("{num}/{}", self.den.to_string_radix(radix))
        }
    }

    /// Parses `num` or `num/den` in `radix` (either case) and reduces it.
    /// The numerator may have a leading `+` or `-`; the denominator has no
    /// sign. As for integers, nothing else is accepted: no whitespace,
    /// separators or prefix.
    ///
    /// # Errors
    ///
    /// [`ParseErrorKind::Empty`] for an empty string,
    /// [`ParseErrorKind::ZeroDenominator`] for a zero denominator, and
    /// [`ParseErrorKind::InvalidDigit`] for any other malformed input,
    /// including an empty numerator or denominator.
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    pub fn from_str_radix(s: &str, radix: u32) -> Result<Rational, ParseNumberError> {
        check_radix(radix);
        let Some((num, den)) = s.split_once('/') else {
            return Integer::from_str_radix(s, radix).map(Rational::from);
        };
        if num.is_empty() {
            return Err(ParseNumberError::new(ParseErrorKind::InvalidDigit));
        }
        let num = Integer::from_str_radix(num, radix)?;
        let den = Natural::from_digit_bytes(den.as_bytes(), radix)?;
        Rational::checked_new(num, Integer::from(den))
            .ok_or_else(|| ParseNumberError::new(ParseErrorKind::ZeroDenominator))
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = self.abs().to_string_radix(10);
        f.pad_integral(!self.is_negative(), "", &s)
    }
}

impl fmt::Debug for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl FromStr for Rational {
    type Err = ParseNumberError;

    /// Decimal, as [`Rational::from_str_radix`] with radix 10.
    fn from_str(s: &str) -> Result<Rational, ParseNumberError> {
        Rational::from_str_radix(s, 10)
    }
}
