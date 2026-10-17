//! Error types.

use std::fmt;

/// A conversion failed because the value is out of range for the target type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct OutOfRangeError(pub(crate) ());

impl fmt::Display for OutOfRangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("value out of range for the target type")
    }
}

impl std::error::Error for OutOfRangeError {}

/// Parsing a number from a string failed.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ParseNumberError {
    kind: ParseErrorKind,
}

/// Why parsing a number failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ParseErrorKind {
    /// The string is empty.
    Empty,
    /// A character is not a digit of the radix, or a sign has no digits
    /// after it.
    InvalidDigit,
    /// A fraction has a zero denominator.
    ZeroDenominator,
}

impl ParseNumberError {
    pub(crate) fn new(kind: ParseErrorKind) -> Self {
        ParseNumberError { kind }
    }

    /// The reason parsing failed.
    #[must_use]
    pub fn kind(&self) -> ParseErrorKind {
        self.kind
    }
}

impl fmt::Display for ParseNumberError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self.kind {
            ParseErrorKind::Empty => "cannot parse number from empty string",
            ParseErrorKind::InvalidDigit => "invalid digit found in string",
            ParseErrorKind::ZeroDenominator => "denominator is zero",
        })
    }
}

impl std::error::Error for ParseNumberError {}
