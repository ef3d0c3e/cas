//! Conversion to and from strings in radix 2 to 36.

use std::fmt;
use std::str::FromStr;

use super::{Integer, Repr};
use crate::Natural;
use crate::error::{ParseErrorKind, ParseNumberError};
use crate::natural::check_radix;

#[cfg(test)]
mod tests;

impl Integer {
    /// The digits in `radix`, lowercase, with a leading `-` for negative
    /// values and no prefix or leading zeros.
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    #[must_use]
    pub fn to_string_radix(&self, radix: u32) -> String {
        check_radix(radix);
        let (negative, digits) = match self.0 {
            Repr::Small(x) if radix == 10 => (x < 0, x.unsigned_abs().to_string()),
            _ => {
                let (negative, mag) = self.sign_mag();
                (negative, mag.to_string_radix(radix))
            }
        };
        if negative {
            format!("-{digits}")
        } else {
            digits
        }
    }

    /// Parses digits in `radix` (either case), with an optional leading `+`
    /// or `-`. Like the primitive integers' `from_str_radix`, nothing else is
    /// accepted: no whitespace, separators or prefix. `-0` is zero.
    ///
    /// # Errors
    ///
    /// [`ParseErrorKind::Empty`] for an empty string,
    /// [`ParseErrorKind::InvalidDigit`] for any other character, or for a
    /// sign with no digits.
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    pub fn from_str_radix(s: &str, radix: u32) -> Result<Integer, ParseNumberError> {
        check_radix(radix);
        let (negative, digits) = match s.as_bytes() {
            [] => return Err(ParseNumberError::new(ParseErrorKind::Empty)),
            [b'+', rest @ ..] => (false, rest),
            [b'-', rest @ ..] => (true, rest),
            all => (false, all),
        };
        let mag = Natural::from_digit_bytes(digits, radix)?;
        Ok(Integer::from_parts(negative, mag))
    }
}

impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = self.to_string_radix(10);
        let (non_negative, digits) = match s.strip_prefix('-') {
            Some(d) => (false, d),
            None => (true, s.as_str()),
        };
        f.pad_integral(non_negative, "", digits)
    }
}

impl fmt::Debug for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Binary for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (negative, mag) = self.sign_mag();
        f.pad_integral(!negative, "0b", &mag.to_string_radix(2))
    }
}

impl fmt::Octal for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (negative, mag) = self.sign_mag();
        f.pad_integral(!negative, "0o", &mag.to_string_radix(8))
    }
}

impl FromStr for Integer {
    type Err = ParseNumberError;

    /// Decimal, as [`Integer::from_str_radix`] with radix 10.
    fn from_str(s: &str) -> Result<Integer, ParseNumberError> {
        Integer::from_str_radix(s, 10)
    }
}
