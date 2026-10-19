//! The lexer: text to lossless tokens (CAS-DOC-7 §1).
//!
//! Every byte of the text belongs to exactly one token, trivia included, so the
//! token texts concatenate back to the text. Lexing never fails: a character
//! that starts no token becomes an [`Error`](SyntaxKind::Error) token, and each
//! problem, whether an error token or a fault inside a well-formed one (a bad
//! escape in a string), is reported as a [`LexError`] with its byte range.
//!
//! The rules that are not plain longest-match:
//! - A number takes `_` between digits; `0x`, `0o` and `0b` only when a digit
//!   of that base follows (`0x^2` is `0·x^2`); a fractional part only when a
//!   digit follows the point (`2.*x` is `2 .* x`, `1..2` a range); and an
//!   exponent only when a digit follows the `e` and its optional sign (`2e` is
//!   `2·e`, `2e-x` is `2·e - x`).
//! - `π` and `ℯ` never join a longer name, so `2πr` is three tokens, and `·`
//!   is a product inside names too.

use std::fmt;
use std::ops::Range;

use unicode_ident::{is_xid_continue, is_xid_start};

use crate::syntax::SyntaxKind;

/// The tokens of a text.
///
/// Token `i` covers [`range(i)`](Lexed::range). The ranges are non-empty,
/// contiguous and cover the whole text.
///
/// ```
/// use cas_lang::{SyntaxKind, lex};
///
/// let lexed = lex("2x² ≤ 1");
/// let tokens: Vec<_> = lexed.iter().filter(|(k, _)| *k != SyntaxKind::Whitespace).collect();
/// assert_eq!(
///     tokens,
///     [
///         (SyntaxKind::Number, "2"),
///         (SyntaxKind::Ident, "x"),
///         (SyntaxKind::Superscript, "²"),
///         (SyntaxKind::Le, "≤"),
///         (SyntaxKind::Number, "1"),
///     ]
/// );
/// assert!(lexed.errors().is_empty());
/// ```
#[derive(Clone, Debug)]
pub struct Lexed<'a> {
    text: &'a str,
    kinds: Vec<SyntaxKind>,
    /// Where each token starts, then the length of the text.
    starts: Vec<usize>,
    errors: Vec<LexError>,
}

impl<'a> Lexed<'a> {
    /// The text that was lexed.
    #[must_use]
    pub fn source(&self) -> &'a str {
        self.text
    }

    /// The number of tokens.
    #[must_use]
    pub fn len(&self) -> usize {
        self.kinds.len()
    }

    /// Whether there are no tokens, which is the case exactly for the empty
    /// text.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.kinds.is_empty()
    }

    /// The kind of token `i`.
    ///
    /// # Panics
    ///
    /// If `i` is not below [`len`](Lexed::len).
    #[must_use]
    pub fn kind(&self, i: usize) -> SyntaxKind {
        self.kinds[i]
    }

    /// The byte range of token `i`.
    ///
    /// # Panics
    ///
    /// If `i` is not below [`len`](Lexed::len).
    #[must_use]
    pub fn range(&self, i: usize) -> Range<usize> {
        self.starts[i]..self.starts[i + 1]
    }

    /// The text of token `i`.
    ///
    /// # Panics
    ///
    /// If `i` is not below [`len`](Lexed::len).
    #[must_use]
    pub fn text(&self, i: usize) -> &'a str {
        &self.text[self.range(i)]
    }

    /// The tokens as `(kind, text)`, in order.
    pub fn iter(&self) -> impl Iterator<Item = (SyntaxKind, &'a str)> {
        (0..self.len()).map(|i| (self.kind(i), self.text(i)))
    }

    /// The problems found, ordered by where they start.
    #[must_use]
    pub fn errors(&self) -> &[LexError] {
        &self.errors
    }
}

/// A problem in the text: where it is and what it is. `Display` gives a
/// message for the user.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LexError {
    /// The bytes concerned.
    pub range: Range<usize>,
    /// What is wrong.
    pub kind: LexErrorKind,
}

/// What is wrong in a [`LexError`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum LexErrorKind {
    /// A character that starts no token.
    UnexpectedChar(char),
    /// `'`: there is no prime notation; derivatives go through `diff`
    /// (CAS-DOC-7 L13).
    Apostrophe,
    /// `|`: absolute values are written `abs(x)`.
    Bar,
    /// `&`: the conjunction is `and`.
    Ampersand,
    /// `#`: comments start with `//`.
    Hash,
    /// A `.` that is part of no operator and no decimal.
    LoneDot,
    /// A `:` that is part of neither `::` nor `:=`.
    LoneColon,
    /// A decimal with no digit before the point: `.5`.
    LeadingPoint,
    /// A `_` in a number that is not between two digits.
    MisplacedUnderscore,
    /// A digit out of range for a `0x`, `0o` or `0b` number.
    InvalidDigit {
        /// 16, 8 or 2.
        radix: u32,
    },
    /// A string with no closing `"`.
    UnterminatedString,
    /// A block comment with no closing `*/`.
    UnterminatedComment,
    /// A backslash that starts none of `\"`, `\\`, `\n`, `\t` and `\u{…}` with
    /// 1 to 6 hexadecimal digits naming a character.
    InvalidEscape,
    /// A superscript sign with no superscript digit after it.
    MissingSuperscriptDigits,
}

impl fmt::Display for LexErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LexErrorKind::UnexpectedChar(c) => write!(f, "unexpected character {c:?}"),
            LexErrorKind::Apostrophe => f.write_str(
                "`'` is not part of the language: write derivatives with diff, as in diff(f(x), x)",
            ),
            LexErrorKind::Bar => f.write_str(
                "`|` is not part of the language: write absolute values with abs, as in abs(x)",
            ),
            LexErrorKind::Ampersand => {
                f.write_str("`&` is not part of the language: the conjunction is `and`")
            }
            LexErrorKind::Hash => {
                f.write_str("`#` is not part of the language: comments start with //")
            }
            LexErrorKind::LoneDot => f.write_str(
                "unexpected `.`: a decimal has digits on both sides of the point (0.5, 1.0), \
                 and the element-wise operators are .*, ./ and .^",
            ),
            LexErrorKind::LoneColon => {
                f.write_str("unexpected `:`: definitions use :=, and types use ::")
            }
            LexErrorKind::LeadingPoint => {
                f.write_str("a decimal needs a digit before the point: write 0.5, not .5")
            }
            LexErrorKind::MisplacedUnderscore => {
                f.write_str("`_` in a number must stand between two digits")
            }
            LexErrorKind::InvalidDigit { radix } => {
                write!(f, "digit out of range in a base-{radix} number")
            }
            LexErrorKind::UnterminatedString => f.write_str("string without its closing \""),
            LexErrorKind::UnterminatedComment => f.write_str("comment without its closing */"),
            LexErrorKind::InvalidEscape => f.write_str(
                "invalid escape: the escapes are \\\", \\\\, \\n, \\t, and \\u{…} with 1 to 6 \
                 hexadecimal digits",
            ),
            LexErrorKind::MissingSuperscriptDigits => {
                f.write_str("a superscript sign needs superscript digits after it")
            }
        }
    }
}

impl fmt::Display for LexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.kind.fmt(f)
    }
}

impl std::error::Error for LexError {}

/// Splits `text` into tokens.
#[must_use]
pub fn lex(text: &str) -> Lexed<'_> {
    let mut lexer = Lexer {
        text,
        pos: 0,
        errors: Vec::new(),
    };
    let mut kinds = Vec::new();
    let mut starts = Vec::new();
    while lexer.pos < text.len() {
        starts.push(lexer.pos);
        kinds.push(lexer.token());
    }
    starts.push(text.len());
    // A string's bad escapes are found before the missing quote that starts
    // with the string.
    lexer.errors.sort_by_key(|e| e.range.start);
    Lexed {
        text,
        kinds,
        starts,
        errors: lexer.errors,
    }
}

struct Lexer<'a> {
    text: &'a str,
    pos: usize,
    errors: Vec<LexError>,
}

impl Lexer<'_> {
    fn peek(&self) -> Option<char> {
        self.text[self.pos..].chars().next()
    }

    fn peek_nth(&self, n: usize) -> Option<char> {
        self.text[self.pos..].chars().nth(n)
    }

    fn bump(&mut self) -> Option<char> {
        let c = self.peek()?;
        self.pos += c.len_utf8();
        Some(c)
    }

    fn eat(&mut self, c: char) -> bool {
        let found = self.peek() == Some(c);
        if found {
            self.pos += c.len_utf8();
        }
        found
    }

    fn eat_while(&mut self, mut accept: impl FnMut(char) -> bool) {
        while let Some(c) = self.peek().filter(|&c| accept(c)) {
            self.pos += c.len_utf8();
        }
    }

    fn error(&mut self, range: Range<usize>, kind: LexErrorKind) {
        self.errors.push(LexError { range, kind });
    }

    fn error_token(&mut self, start: usize, kind: LexErrorKind) -> SyntaxKind {
        self.error(start..self.pos, kind);
        SyntaxKind::Error
    }

    /// Lexes the token at `self.pos`, which is before the end of the text.
    fn token(&mut self) -> SyntaxKind {
        let start = self.pos;
        let Some(c) = self.bump() else {
            unreachable!("token() is called before the end of the text")
        };
        match c {
            '\n' => SyntaxKind::Newline,
            '\r' => {
                self.eat('\n');
                SyntaxKind::Newline
            }
            c if is_space(c) => {
                self.eat_while(is_space);
                SyntaxKind::Whitespace
            }
            '0'..='9' => self.number(start, c),
            '.' => self.dot(start),
            '"' => self.string(start),
            '/' if self.eat('/') => {
                self.eat_while(|c| c != '\n' && c != '\r');
                SyntaxKind::LineComment
            }
            '/' if self.eat('*') => self.block_comment(start),
            'π' | 'ℯ' | '∞' => SyntaxKind::Ident,
            c if is_ident_start(c) => {
                self.eat_while(is_ident_continue);
                SyntaxKind::keyword(&self.text[start..self.pos]).unwrap_or(SyntaxKind::Ident)
            }
            c if is_superscript_sign(c) || superscript_digit(c).is_some() => {
                self.superscript(start, c)
            }
            c => self.operator(start, c),
        }
    }

    /// Lexes an operator or punctuation starting with `c`, already consumed,
    /// or an error token.
    fn operator(&mut self, start: usize, c: char) -> SyntaxKind {
        match c {
            '+' => SyntaxKind::Plus,
            '-' if self.eat('>') => SyntaxKind::Arrow,
            '-' | '\u{2212}' => SyntaxKind::Minus,
            '*' | '×' | '·' | '⋅' => SyntaxKind::Star,
            '/' | '÷' => SyntaxKind::Slash,
            '^' => SyntaxKind::Caret,
            '!' if self.eat('=') => SyntaxKind::Neq,
            '!' => SyntaxKind::Bang,
            '≠' => SyntaxKind::Neq,
            '=' if self.eat('=') => SyntaxKind::EqEq,
            '=' => SyntaxKind::Eq,
            '<' if self.eat('=') => SyntaxKind::Le,
            '<' => SyntaxKind::Lt,
            '≤' => SyntaxKind::Le,
            '>' if self.eat('=') => SyntaxKind::Ge,
            '>' => SyntaxKind::Gt,
            '≥' => SyntaxKind::Ge,
            ':' if self.eat(':') => SyntaxKind::ColonColon,
            ':' if self.eat('=') => SyntaxKind::ColonEq,
            '→' => SyntaxKind::Arrow,
            '∈' => SyntaxKind::InKw,
            '(' => SyntaxKind::LParen,
            ')' => SyntaxKind::RParen,
            '[' => SyntaxKind::LBracket,
            ']' => SyntaxKind::RBracket,
            '{' => SyntaxKind::LBrace,
            '}' => SyntaxKind::RBrace,
            ',' => SyntaxKind::Comma,
            ';' => SyntaxKind::Semicolon,
            '√' => SyntaxKind::Sqrt,
            '∛' => SyntaxKind::Cbrt,
            '∜' => SyntaxKind::FourthRoot,
            ':' => self.error_token(start, LexErrorKind::LoneColon),
            '\'' => self.error_token(start, LexErrorKind::Apostrophe),
            '|' => self.error_token(start, LexErrorKind::Bar),
            '&' => self.error_token(start, LexErrorKind::Ampersand),
            '#' => self.error_token(start, LexErrorKind::Hash),
            c => self.error_token(start, LexErrorKind::UnexpectedChar(c)),
        }
    }

    /// Lexes a number whose first digit, `first`, is consumed.
    fn number(&mut self, start: usize, first: char) -> SyntaxKind {
        if first == '0' {
            let radix = match self.peek() {
                Some('x') => 16,
                Some('o') => 8,
                Some('b') => 2,
                _ => 10,
            };
            if radix != 10 && self.peek_nth(1).is_some_and(|c| c.is_digit(radix)) {
                self.bump();
                self.digits(radix);
                let rest = self.pos;
                if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    self.digits(10);
                    self.error(rest..self.pos, LexErrorKind::InvalidDigit { radix });
                }
                return SyntaxKind::Number;
            }
        }
        self.pos = start;
        self.digits(10);
        if self.peek() == Some('.') && self.peek_nth(1).is_some_and(|c| c.is_ascii_digit()) {
            self.bump();
            self.digits(10);
        }
        self.exponent();
        SyntaxKind::Number
    }

    /// Consumes a run of digits in `radix` and `_`, which starts with a digit,
    /// and reports a `_` that is not between two digits.
    fn digits(&mut self, radix: u32) {
        let start = self.pos;
        self.eat_while(|c| c == '_' || c.is_digit(radix));
        let run = &self.text[start..self.pos];
        if run.ends_with('_') || run.contains("__") {
            self.error(start..self.pos, LexErrorKind::MisplacedUnderscore);
        }
    }

    /// Consumes an exponent, `e` or `E` with an optional sign, when a digit
    /// follows.
    fn exponent(&mut self) {
        if matches!(self.peek(), Some('e' | 'E')) {
            let sign = usize::from(matches!(self.peek_nth(1), Some('+' | '-')));
            if self.peek_nth(1 + sign).is_some_and(|c| c.is_ascii_digit()) {
                for _ in 0..=sign {
                    self.bump();
                }
                self.digits(10);
            }
        }
    }

    /// Lexes what follows a `.`: an element-wise operator, `..`, `...`, or a
    /// decimal missing its integer part.
    fn dot(&mut self, start: usize) -> SyntaxKind {
        match self.peek() {
            Some('*') => {
                self.bump();
                SyntaxKind::DotStar
            }
            Some('/') => {
                self.bump();
                SyntaxKind::DotSlash
            }
            Some('^') => {
                self.bump();
                SyntaxKind::DotCaret
            }
            Some('.') => {
                self.bump();
                if self.eat('.') {
                    SyntaxKind::Ellipsis
                } else {
                    SyntaxKind::DotDot
                }
            }
            Some(c) if c.is_ascii_digit() => {
                self.digits(10);
                self.exponent();
                self.error(start..self.pos, LexErrorKind::LeadingPoint);
                SyntaxKind::Number
            }
            _ => self.error_token(start, LexErrorKind::LoneDot),
        }
    }

    /// Lexes a string whose opening quote is consumed.
    fn string(&mut self, start: usize) -> SyntaxKind {
        loop {
            match self.bump() {
                None => {
                    self.error(start..self.pos, LexErrorKind::UnterminatedString);
                    break;
                }
                Some('"') => break,
                Some('\\') => self.escape(self.pos - 1),
                Some(_) => {}
            }
        }
        SyntaxKind::String
    }

    /// Consumes an escape whose backslash, at `start`, is consumed.
    fn escape(&mut self, start: usize) {
        let valid = match self.bump() {
            Some('"' | '\\' | 'n' | 't') => true,
            Some('u') => self.unicode_escape(),
            _ => false,
        };
        if !valid {
            self.error(start..self.pos, LexErrorKind::InvalidEscape);
        }
    }

    /// Consumes the `{…}` of a `\u` escape; whether it names a character.
    fn unicode_escape(&mut self) -> bool {
        if !self.eat('{') {
            return false;
        }
        let digits = self.pos;
        self.eat_while(|c| c.is_ascii_hexdigit());
        let hex = &self.text[digits..self.pos];
        self.eat('}')
            && (1..=6).contains(&hex.len())
            && u32::from_str_radix(hex, 16)
                .ok()
                .and_then(char::from_u32)
                .is_some()
    }

    /// Lexes a block comment whose `/*` is consumed.
    fn block_comment(&mut self, start: usize) -> SyntaxKind {
        let mut depth = 1usize;
        while depth > 0 {
            match self.bump() {
                None => {
                    self.error(start..self.pos, LexErrorKind::UnterminatedComment);
                    break;
                }
                Some('/') => {
                    if self.eat('*') {
                        depth += 1;
                    }
                }
                Some('*') => {
                    if self.eat('/') {
                        depth -= 1;
                    }
                }
                Some(_) => {}
            }
        }
        SyntaxKind::BlockComment
    }

    /// Lexes a superscript exponent whose first character, `first`, is
    /// consumed.
    fn superscript(&mut self, start: usize, first: char) -> SyntaxKind {
        self.eat_while(|c| superscript_digit(c).is_some());
        if is_superscript_sign(first) && self.pos == start + first.len_utf8() {
            self.error(start..self.pos, LexErrorKind::MissingSuperscriptDigits);
        }
        SyntaxKind::Superscript
    }
}

/// Whitespace other than line breaks. A byte order mark counts as whitespace.
fn is_space(c: char) -> bool {
    (c.is_whitespace() && c != '\n' && c != '\r') || c == '\u{feff}'
}

fn is_ident_start(c: char) -> bool {
    c == '_' || (is_xid_start(c) && !is_standalone(c))
}

fn is_ident_continue(c: char) -> bool {
    is_xid_continue(c) && !is_standalone(c) && c != '·'
}

/// Letters that are a token of their own, never part of a longer name.
fn is_standalone(c: char) -> bool {
    matches!(c, 'π' | 'ℯ')
}

fn is_superscript_sign(c: char) -> bool {
    matches!(c, '⁺' | '⁻')
}

/// The value of a superscript digit.
pub(crate) fn superscript_digit(c: char) -> Option<u32> {
    match c {
        '⁰' => Some(0),
        '¹' => Some(1),
        '²' => Some(2),
        '³' => Some(3),
        '\u{2074}'..='\u{2079}' => Some(u32::from(c) - 0x2070),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
