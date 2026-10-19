//! Exact hexadecimal float literals for [`Float`] and [`Mag`] (CAS-DOC-5
//! §6.2), C99's `%a` form \[C99\]: `0x1.921fb54442d18p+1`, `-0x1.8p-3`,
//! `0x0p+0`.
//!
//! The printed form is canonical: a leading `1`, the fraction without
//! trailing zeros, and a decimal exponent of any size. The parser also takes
//! other spellings of the same value (`0x18p-4`, `0X.8P1`, no exponent).
//! Both are exact; they are `FullForm`'s spelling of a ball's midpoint and
//! radius.

use super::{Float, Mag};
use crate::error::{ParseErrorKind, ParseNumberError};
use crate::{Integer, Natural};

#[cfg(test)]
mod tests;

fn invalid() -> ParseNumberError {
    ParseNumberError::new(ParseErrorKind::InvalidDigit)
}

impl Float {
    /// The exact value as a hexadecimal float literal.
    #[must_use]
    #[allow(
        clippy::missing_panics_doc,
        reason = "the expects hold for any float in memory"
    )]
    pub fn to_hex(&self) -> String {
        if self.is_zero() {
            return "0x0p+0".to_owned();
        }
        let (m, e) = self.to_integer_exp();
        let (neg, mut m) = m.into_parts();
        // Odd mantissa: the fraction then ends in a non-zero hex digit.
        let tz = m.trailing_zeros().expect("non-zero");
        m >>= tz;
        let len = m.bit_length();
        let frac_bits = len - 1;
        let pad = (4 - frac_bits % 4) % 4;
        let digits = usize::try_from((frac_bits + pad) / 4).expect("fits memory");
        let mut out = String::with_capacity(digits + 32);
        if neg {
            out.push('-');
        }
        out.push_str("0x1");
        if digits > 0 {
            let frac = (&m - &(&Natural::from(1u64) << frac_bits)) << pad;
            let hex = frac.to_string_radix(16);
            out.push('.');
            out.extend(std::iter::repeat_n('0', digits - hex.len()));
            out.push_str(&hex);
        }
        // value = 1.frac × 2^(e + tz + len - 1).
        let x = e + Integer::from(tz + frac_bits);
        out.push('p');
        if !x.is_negative() {
            out.push('+');
        }
        out.push_str(&x.to_string());
        out
    }

    /// Parses a hexadecimal float literal exactly: an optional sign, `0x`
    /// or `0X`, hex digits with an optional point (at least one digit), and
    /// an optional binary exponent `p` or `P` with a signed decimal integer.
    ///
    /// # Errors
    ///
    /// If `s` is empty or not of that form.
    pub fn from_hex(s: &str) -> Result<Float, ParseNumberError> {
        if s.is_empty() {
            return Err(ParseNumberError::new(ParseErrorKind::Empty));
        }
        let (neg, rest) = match s.as_bytes()[0] {
            b'-' => (true, &s[1..]),
            b'+' => (false, &s[1..]),
            _ => (false, s),
        };
        let rest = rest
            .strip_prefix("0x")
            .or_else(|| rest.strip_prefix("0X"))
            .ok_or_else(invalid)?;
        let (body, exp) = match rest.find(['p', 'P']) {
            Some(i) => (&rest[..i], Some(&rest[i + 1..])),
            None => (rest, None),
        };
        let (int, frac) = body.split_once('.').unwrap_or((body, ""));
        let hex = |t: &str| t.bytes().all(|c| c.is_ascii_hexdigit());
        if int.len() + frac.len() == 0 || !hex(int) || !hex(frac) {
            return Err(invalid());
        }
        let mut digits = String::with_capacity(int.len() + frac.len());
        digits.push_str(int);
        digits.push_str(frac);
        let m = Natural::from_str_radix(&digits, 16)?;
        let exp = match exp {
            None => Integer::ZERO,
            Some(t) => {
                // A decimal integer with an optional sign, nothing else.
                let unsigned = t.strip_prefix(['+', '-']).unwrap_or(t);
                if unsigned.is_empty() || !unsigned.bytes().all(|c| c.is_ascii_digit()) {
                    return Err(invalid());
                }
                Integer::from_str_radix(t.strip_prefix('+').unwrap_or(t), 10)?
            }
        };
        let e = exp - Integer::from(4 * frac.len());
        Ok(Float::from_integer_exp(&Integer::from_parts(neg, m), &e))
    }
}

impl Mag {
    /// The exact value as a hexadecimal float literal, `inf` for `+∞`.
    #[must_use]
    pub fn to_hex(&self) -> String {
        match self.to_float() {
            Some(x) => x.to_hex(),
            None => "inf".to_owned(),
        }
    }

    /// An upper bound of a non-negative hexadecimal float literal (exact up
    /// to 30 significant bits), or `inf`.
    ///
    /// # Errors
    ///
    /// If `s` is not of [`Float::from_hex`]'s form, or is negative.
    pub fn from_hex(s: &str) -> Result<Mag, ParseNumberError> {
        if s == "inf" {
            return Ok(Mag::INF);
        }
        let x = Float::from_hex(s)?;
        if x.is_negative() {
            return Err(invalid());
        }
        Ok(Mag::from_float(&x))
    }
}
