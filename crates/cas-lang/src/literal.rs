//! The values of literal tokens: numbers, strings, superscript exponents and
//! names.
//!
//! Each function accepts exactly the texts the lexer produces as one
//! error-free token of its kind, and nothing else.

use std::borrow::Cow;
use std::fmt;

use cas_num::{Integer, Natural, Rational};
use unicode_normalization::{UnicodeNormalization, is_nfc};

use crate::lexer::{lex, superscript_digit};
use crate::syntax::SyntaxKind;

/// The largest magnitude of a number's exponent: `1e1000000` is accepted.
///
/// The bound keeps a short literal from becoming a number of gigabytes;
/// larger powers are written `10^n`.
pub const MAX_DECIMAL_EXPONENT: u64 = 1_000_000;

/// Why [`number_value`] gives no value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum NumberError {
    /// The text is not exactly one error-free [`Number`](SyntaxKind::Number)
    /// token.
    Malformed,
    /// The exponent is larger in magnitude than [`MAX_DECIMAL_EXPONENT`].
    ExponentTooLarge,
}

impl fmt::Display for NumberError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NumberError::Malformed => f.write_str("not a number"),
            NumberError::ExponentTooLarge => write!(
                f,
                "exponent beyond ±{MAX_DECIMAL_EXPONENT}: write a larger power as 10^n"
            ),
        }
    }
}

impl std::error::Error for NumberError {}

/// The exact value of a number token. Decimals are exact (CAS-DOC-5 §6.4):
/// `0.1` is `1/10`.
///
/// ```
/// use cas_lang::number_value;
/// use cas_num::{Integer, Rational};
///
/// let q = |n: i64, d: i64| Rational::new(n.into(), d.into());
/// assert_eq!(number_value("0.1"), Ok(q(1, 10)));
/// assert_eq!(number_value("1.5e-3"), Ok(q(3, 2000)));
/// assert_eq!(number_value("1_000"), Ok(q(1000, 1)));
/// assert_eq!(number_value("0x1f"), Ok(Rational::from(Integer::from(31))));
/// ```
///
/// # Errors
///
/// [`NumberError::Malformed`] unless `text` lexes as exactly one error-free
/// number token, and [`NumberError::ExponentTooLarge`] for an exponent beyond
/// [`MAX_DECIMAL_EXPONENT`].
pub fn number_value(text: &str) -> Result<Rational, NumberError> {
    if !is_token(text, SyntaxKind::Number) {
        return Err(NumberError::Malformed);
    }
    let digits: String = text.chars().filter(|&c| c != '_').collect();
    let radix = match digits.get(..2) {
        Some("0x") => 16,
        Some("0o") => 8,
        Some("0b") => 2,
        _ => 10,
    };
    if radix != 10 {
        return natural(&digits[2..], radix).map(|n| Rational::from(Integer::from(n)));
    }
    let (mantissa, exponent) = match digits.split_once(['e', 'E']) {
        Some((mantissa, exponent)) => (mantissa, decimal_exponent(exponent)?),
        None => (digits.as_str(), 0),
    };
    let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let n = Integer::from(natural(&[int, frac].concat(), 10)?);
    let shift = exponent.saturating_sub_unsigned(u64::try_from(frac.len()).unwrap_or(u64::MAX));
    let scale = Integer::from(Natural::from(10u32).pow(shift.unsigned_abs()));
    Ok(if shift >= 0 {
        Rational::from(n * scale)
    } else {
        Rational::new(n, scale)
    })
}

/// The digits of a well-formed number in `radix`.
fn natural(digits: &str, radix: u32) -> Result<Natural, NumberError> {
    Natural::from_str_radix(digits, radix).map_err(|_| NumberError::Malformed)
}

/// The exponent of a decimal: an optional sign, then decimal digits.
fn decimal_exponent(text: &str) -> Result<i64, NumberError> {
    // Only overflow can fail: the lexer checked the digits.
    let e: i64 = text.parse().map_err(|_| NumberError::ExponentTooLarge)?;
    if e.unsigned_abs() > MAX_DECIMAL_EXPONENT {
        return Err(NumberError::ExponentTooLarge);
    }
    Ok(e)
}

/// The contents of a string token: the quotes removed, the escapes replaced.
///
/// ```
/// use cas_lang::string_value;
///
/// assert_eq!(string_value(r#""a\"b\\c""#).as_deref(), Some(r#"a"b\c"#));
/// assert_eq!(string_value(r#""\u{3c0}\n""#).as_deref(), Some("π\n"));
/// assert_eq!(string_value(r#""\q""#), None);
/// ```
///
/// `None` unless `text` lexes as exactly one error-free string token.
#[must_use]
pub fn string_value(text: &str) -> Option<String> {
    if !is_token(text, SyntaxKind::String) {
        return None;
    }
    let body = &text[1..text.len() - 1];
    let mut out = String::with_capacity(body.len());
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next()? {
            'n' => out.push('\n'),
            't' => out.push('\t'),
            'u' => {
                let rest = chars.as_str();
                let close = rest.find('}')?;
                let code = u32::from_str_radix(&rest[1..close], 16).ok()?;
                out.push(char::from_u32(code)?);
                chars = rest[close + 1..].chars();
            }
            // `\"` and `\\`.
            c => out.push(c),
        }
    }
    Some(out)
}

/// The exponent a superscript stands for: `²` is 2, `⁻¹` is −1.
///
/// ```
/// use cas_lang::superscript_value;
/// use cas_num::Integer;
///
/// assert_eq!(superscript_value("²"), Some(Integer::from(2)));
/// assert_eq!(superscript_value("⁻¹⁰"), Some(Integer::from(-10)));
/// assert_eq!(superscript_value("⁻"), None);
/// ```
///
/// `None` unless `text` lexes as exactly one error-free superscript token.
#[must_use]
pub fn superscript_value(text: &str) -> Option<Integer> {
    if !is_token(text, SyntaxKind::Superscript) {
        return None;
    }
    let negative = text.starts_with('⁻');
    let digits: String = text
        .chars()
        .filter_map(|c| superscript_digit(c).and_then(|d| char::from_digit(d, 10)))
        .collect();
    let n = Natural::from_str_radix(&digits, 10).ok()?;
    Some(Integer::from_parts(negative, n))
}

/// The name an identifier token stands for: its text in Unicode normal form
/// C, so that a name typed with precomposed or combining characters is one
/// name, and `pi`, `e` and `inf` for `π`, `ℯ` and `∞`.
///
/// ```
/// use cas_lang::name;
///
/// assert_eq!(name("θ"), "θ");
/// assert_eq!(name("e\u{301}"), "\u{e9}");
/// assert_eq!(name("π"), "pi");
/// ```
#[must_use]
pub fn name(text: &str) -> Cow<'_, str> {
    match text {
        "π" => Cow::Borrowed("pi"),
        "ℯ" => Cow::Borrowed("e"),
        "∞" => Cow::Borrowed("inf"),
        _ if is_nfc(text) => Cow::Borrowed(text),
        _ => Cow::Owned(text.nfc().collect()),
    }
}

/// Whether `text` lexes as one token of kind `kind`, without errors.
fn is_token(text: &str, kind: SyntaxKind) -> bool {
    let lexed = lex(text);
    lexed.len() == 1 && lexed.kind(0) == kind && lexed.errors().is_empty()
}

#[cfg(test)]
mod tests;
