//! The values of literal tokens: numbers, strings, superscript exponents and
//! names.
//!
//! Each function accepts exactly the texts the lexer produces as one
//! error-free token of its kind, and nothing else.

use std::borrow::Cow;
use std::fmt;

use cas_core::{Expr, Symbol, View};
use cas_num::ball::Ball;
use cas_num::{Integer, Natural};
use unicode_normalization::{UnicodeNormalization, is_nfc};

use crate::lexer::{lex, superscript_digit};
use crate::syntax::SyntaxKind;

/// The largest magnitude of a number's exponent: `1e10000` is accepted.
///
/// A decimal is a ball (CAS-148), whose exponent could be any size, but the
/// bound keeps printing it cheap: its exact digits are found with powers of
/// 10 of the size of the exponent, a few milliseconds at this bound and
/// seconds at `10^6`, CAS-117's bound for exact decimals. Larger powers are
/// written `10^n`.
pub const MAX_DECIMAL_EXPONENT: u64 = 10_000;

/// The value of a number token: an integer, exact, or a decimal, which is
/// approximate (CAS-DOC-5 §6.4, D1): a ball once its statement runs.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Number {
    /// `42`, `1_000_000`, `0x1f`.
    Integer(Integer),
    /// A decimal or scientific literal, `digits · 10^exponent`, with the
    /// digits as written: `2.50e-3` is `250 · 10^-5`, its trailing zero
    /// counted, since the digits typed set its precision.
    Decimal {
        /// The digits, without the point.
        digits: Integer,
        /// The power of 10 they are scaled by.
        exponent: i64,
    },
}

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

/// The value of a number token: an integer, or a decimal with its digits
/// as written (CAS-DOC-5 §6.4).
///
/// ```
/// use cas_lang::{Number, number_value};
/// use cas_num::Integer;
///
/// let decimal = |digits: i64, exponent| Ok(Number::Decimal { digits: digits.into(), exponent });
/// assert_eq!(number_value("0.1"), decimal(1, -1));
/// assert_eq!(number_value("1.50e-3"), decimal(150, -5));
/// assert_eq!(number_value("2E10"), decimal(2, 10));
/// assert_eq!(number_value("1_000"), Ok(Number::Integer(Integer::from(1000))));
/// assert_eq!(number_value("0x1f"), Ok(Number::Integer(Integer::from(31))));
/// ```
///
/// # Errors
///
/// [`NumberError::Malformed`] unless `text` lexes as exactly one error-free
/// number token, and [`NumberError::ExponentTooLarge`] for an exponent beyond
/// [`MAX_DECIMAL_EXPONENT`].
pub fn number_value(text: &str) -> Result<Number, NumberError> {
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
        return natural(&digits[2..], radix).map(|n| Number::Integer(Integer::from(n)));
    }
    let (mantissa, exponent) = match digits.split_once(['e', 'E']) {
        Some((mantissa, exponent)) => (mantissa, Some(decimal_exponent(exponent)?)),
        None => (digits.as_str(), None),
    };
    let point = mantissa.split_once('.');
    let n = Integer::from(natural(
        &point.map_or(Cow::Borrowed(mantissa), |(int, frac)| {
            Cow::Owned([int, frac].concat())
        }),
        10,
    )?);
    if point.is_none() && exponent.is_none() {
        return Ok(Number::Integer(n));
    }
    let frac = point.map_or(0, |(_, frac)| frac.len());
    let shift = exponent
        .unwrap_or(0)
        .saturating_sub_unsigned(u64::try_from(frac).unwrap_or(u64::MAX));
    Ok(Number::Decimal {
        digits: n,
        exponent: shift,
    })
}

/// `e` with each decimal literal of the language, `Decimal[digits,
/// exponent]`, made its ball at `working_precision` bits, or more for a
/// literal with more digits than that precision reads back
/// ([`Ball::from_literal`]). The session does this when it runs a
/// statement, held parts included, so a precision set by a statement
/// applies to the literals of the statements after it (CAS-148).
#[must_use]
pub fn decimal_balls(e: &Expr, working_precision: u64) -> Expr {
    e.replace(|part| {
        let (Some(Symbol::DECIMAL), [digits, exponent]) =
            (part.head().and_then(Expr::as_symbol), part.args())
        else {
            return None;
        };
        let View::Int(exponent) = exponent.view() else {
            return None;
        };
        let digits = digits.to_integer()?;
        Some(Expr::from(Ball::from_literal(
            &digits,
            exponent,
            working_precision,
        )))
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
