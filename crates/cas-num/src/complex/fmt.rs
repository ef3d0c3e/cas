//! Printing and parsing in the form `a + b*I`.
//!
//! `Display` writes the real part, then ` + ` or ` - ` and the magnitude of the
//! imaginary part as `b*I` (`I` when it is 1), leaving out a zero part:
//! `3`, `-1/2*I`, `I`, `1/2 - 3/4*I`. Parsing accepts exactly what `Display`
//! writes, with the parts in any form `Rational` parses (so `1 + 0*I` and
//! `2/4*I` are accepted and reduce).

use std::fmt;
use std::str::FromStr;

use super::Complex;
use crate::{ParseErrorKind, ParseNumberError, Rational};

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.im.is_zero() {
            return write!(f, "{}", self.re);
        }
        let mag = self.im.abs();
        let imag = |f: &mut fmt::Formatter<'_>| {
            if mag.is_one() {
                f.write_str("I")
            } else {
                write!(f, "{mag}*I")
            }
        };
        if self.re.is_zero() {
            if self.im.is_negative() {
                f.write_str("-")?;
            }
            return imag(f);
        }
        write!(f, "{}", self.re)?;
        f.write_str(if self.im.is_negative() { " - " } else { " + " })?;
        imag(f)
    }
}

impl fmt::Debug for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

/// The magnitude of an imaginary part written without a sign: `I` or `b*I`.
fn unsigned_imag(s: &str) -> Result<Rational, ParseNumberError> {
    if s == "I" {
        return Ok(Rational::ONE);
    }
    let invalid = || ParseNumberError::new(ParseErrorKind::InvalidDigit);
    let b = s.strip_suffix("*I").ok_or_else(invalid)?;
    if b.is_empty() || b.starts_with(['+', '-']) {
        return Err(invalid());
    }
    b.parse()
}

impl FromStr for Complex {
    type Err = ParseNumberError;

    fn from_str(s: &str) -> Result<Complex, ParseNumberError> {
        if s.is_empty() {
            return Err(ParseNumberError::new(ParseErrorKind::Empty));
        }
        if !s.ends_with('I') {
            return Ok(Complex::from(s.parse::<Rational>()?));
        }
        // The last separator, if any: rationals contain no spaces.
        let sep = [" + ", " - "].iter().filter_map(|p| s.rfind(p)).max();
        Ok(match sep {
            Some(i) => {
                let re = s[..i].parse()?;
                let mag = unsigned_imag(&s[i + 3..])?;
                let im = if &s[i..=i + 1] == " -" { -mag } else { mag };
                Complex::new(re, im)
            }
            None => match s.strip_prefix('-') {
                Some(rest) => Complex::new(Rational::ZERO, -unsigned_imag(rest)?),
                None => Complex::new(Rational::ZERO, unsigned_imag(s)?),
            },
        })
    }
}
