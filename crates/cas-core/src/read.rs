//! Reading full form, the inverse of [`Expr::full_form`] (CAS-52).
//!
//! The reader runs on an explicit stack, as the printers do, so text of any
//! depth reads in constant native stack.

use std::fmt;
use std::str::FromStr;

use cas_num::{Complex, Rational};

use crate::expr::Expr;
use crate::symbol::Symbol;

/// Why [`Expr::parse_full_form`] cannot read a text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FullFormError {
    /// The byte offset in the text where reading stopped.
    pub position: usize,
    /// Why.
    pub kind: FullFormErrorKind,
}

/// What a [`FullFormError`] is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FullFormErrorKind {
    /// The text ends where an expression, `,` or `]` should come.
    End,
    /// A character that cannot come here.
    Unexpected(char),
    /// A `-` without digits after it.
    MissingDigits,
    /// A rational with a zero denominator.
    ZeroDenominator,
    /// A string without its closing quote.
    UnterminatedString,
    /// An escape other than `\"`, `\\`, `\n`, `\t`, `\r` and `\u{…}`.
    InvalidEscape,
}

impl fmt::Display for FullFormError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "at byte {}: ", self.position)?;
        match self.kind {
            FullFormErrorKind::End => f.write_str("unexpected end of text"),
            FullFormErrorKind::Unexpected(c) => write!(f, "unexpected {c:?}"),
            FullFormErrorKind::MissingDigits => f.write_str("`-` without digits"),
            FullFormErrorKind::ZeroDenominator => f.write_str("zero denominator"),
            FullFormErrorKind::UnterminatedString => {
                f.write_str("string without its closing quote")
            }
            FullFormErrorKind::InvalidEscape => f.write_str("invalid escape"),
        }
    }
}

impl std::error::Error for FullFormError {}

impl Expr {
    /// Reads an expression in full form, as [`Expr::full_form`] and `Debug`
    /// print it: `head[args…]`, compound heads `f[a][b]`, and the atoms:
    ///
    /// - integers of any size and rationals, `-3/2`, reduced;
    /// - strings, with the escapes `\"`, `\\`, `\n`, `\t`, `\r` and `\u{…}`;
    /// - symbols: a letter, `_` or `$`, then letters, digits, `_` and `$`
    ///   (Unicode identifiers), so `$E`, the language's escaped `E`, reads
    ///   back;
    /// - complex numbers: `Complex[re, im]` with two numbers is the number,
    ///   never an application, so `Complex[2, 0]` is `2` (CAS-51).
    ///
    /// Whitespace between tokens is ignored. Nothing is evaluated or
    /// normalised: `Plus[1, 2]` stays a sum. Any text that [`Expr::full_form`]
    /// prints reads back as the same expression, if its symbols' names are
    /// identifiers.
    ///
    /// ```
    /// use cas_core::{Expr, build};
    ///
    /// let text = "Plus[x, Times[-3/2, Power[$E, 2]], Complex[0, 1]]";
    /// let e = Expr::parse_full_form(text).unwrap();
    /// assert_eq!(format!("{e:?}"), text);
    /// assert_eq!(e.args()[2], build::i());
    /// assert_eq!(Expr::parse_full_form("f[x").unwrap_err().position, 3);
    /// ```
    ///
    /// # Errors
    ///
    /// A [`FullFormError`] at the first byte that does not fit.
    pub fn parse_full_form(text: &str) -> Result<Expr, FullFormError> {
        Reader { text, pos: 0 }.read()
    }
}

/// Asserts that evaluating `input` gives `expected`, both texts in full
/// form ([`Expr::parse_full_form`]), in a new [`Context`](crate::Context) or
/// in the one given first. For readable evaluator tests (CAS-DOC-4 §10).
///
/// ```
/// use cas_core::{Context, assert_eval};
///
/// assert_eval!("f[x, Complex[0, 1]]", "f[x, Complex[0, 1]]");
/// assert_eval!(&Context::new(), "List[]", "List[]");
/// ```
///
/// # Panics
///
/// If a text does not read, evaluation fails or the result is not
/// `expected`.
#[macro_export]
macro_rules! assert_eval {
    ($input:expr, $expected:expr $(,)?) => {
        $crate::assert_eval!(&$crate::Context::new(), $input, $expected)
    };
    ($ctx:expr, $input:expr, $expected:expr $(,)?) => {{
        let read = |text: &str| {
            $crate::Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
        };
        let (input, expected) = (read($input), read($expected));
        let found = $ctx
            .eval(&input)
            .unwrap_or_else(|e| panic!("{input:?}: {e}"));
        assert_eq!(found, expected, "{input:?} evaluated to {found:?}");
    }};
}

struct Reader<'a> {
    text: &'a str,
    pos: usize,
}

impl Reader<'_> {
    fn read(mut self) -> Result<Expr, FullFormError> {
        // The applications whose arguments are being read: head, arguments so
        // far.
        let mut open: Vec<(Expr, Vec<Expr>)> = Vec::new();
        'expression: loop {
            let mut e = self.atom()?;
            loop {
                if self.eat('[') {
                    if self.eat(']') {
                        e = apply(&e, &[]);
                        continue;
                    }
                    open.push((e, Vec::new()));
                    continue 'expression;
                }
                // `e` is complete: an argument of the innermost open
                // application, or the whole text.
                let Some((_, args)) = open.last_mut() else {
                    self.skip_space();
                    return match self.peek() {
                        None => Ok(e),
                        Some(c) => Err(self.error(FullFormErrorKind::Unexpected(c))),
                    };
                };
                args.push(e);
                if self.eat(',') {
                    continue 'expression;
                }
                if !self.eat(']') {
                    return Err(self.unexpected());
                }
                let (head, args) = open.pop().expect("an open application");
                e = apply(&head, &args);
            }
        }
    }

    fn peek(&self) -> Option<char> {
        self.text[self.pos..].chars().next()
    }

    fn skip_space(&mut self) {
        let rest = &self.text[self.pos..];
        self.pos += rest.len() - rest.trim_start().len();
    }

    /// Consumes `c` after any whitespace, if it comes next.
    fn eat(&mut self, c: char) -> bool {
        self.skip_space();
        if self.peek() == Some(c) {
            self.pos += c.len_utf8();
            true
        } else {
            false
        }
    }

    fn error(&self, kind: FullFormErrorKind) -> FullFormError {
        FullFormError {
            position: self.pos,
            kind,
        }
    }

    /// The error for what comes next: the end, or a character out of place.
    fn unexpected(&self) -> FullFormError {
        match self.peek() {
            None => self.error(FullFormErrorKind::End),
            Some(c) => self.error(FullFormErrorKind::Unexpected(c)),
        }
    }

    /// A number, a string or a symbol.
    fn atom(&mut self) -> Result<Expr, FullFormError> {
        self.skip_space();
        match self.peek() {
            Some('"') => self.string(),
            Some(c) if c == '-' || c.is_ascii_digit() => self.number(),
            Some(c) if is_name_start(c) => {
                let start = self.pos;
                self.pos += self.text[start..]
                    .find(|c| !is_name_continue(c))
                    .unwrap_or(self.text.len() - start);
                Ok(Symbol::new(&self.text[start..self.pos]).into())
            }
            _ => Err(self.unexpected()),
        }
    }

    /// `-?digits(/digits)?`.
    fn number(&mut self) -> Result<Expr, FullFormError> {
        let start = self.pos;
        let digits = |r: &mut Self| {
            let from = r.pos;
            r.pos += r.text[from..]
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(r.text.len() - from);
            r.pos > from
        };
        if self.text[start..].starts_with('-') {
            self.pos += 1;
        }
        if !digits(self) {
            return Err(FullFormError {
                position: start,
                kind: FullFormErrorKind::MissingDigits,
            });
        }
        if self.text[self.pos..].starts_with('/') {
            self.pos += 1;
            if !digits(self) {
                return Err(self.unexpected());
            }
        }
        match Rational::from_str(&self.text[start..self.pos]) {
            Ok(q) => Ok(q.into()),
            Err(_) => Err(FullFormError {
                position: start,
                kind: FullFormErrorKind::ZeroDenominator,
            }),
        }
    }

    /// `"…"`, with the escapes the printer writes.
    fn string(&mut self) -> Result<Expr, FullFormError> {
        let start = self.pos;
        self.pos += 1;
        let mut s = String::new();
        loop {
            let Some(c) = self.peek() else {
                return Err(FullFormError {
                    position: start,
                    kind: FullFormErrorKind::UnterminatedString,
                });
            };
            self.pos += c.len_utf8();
            match c {
                '"' => return Ok(Expr::string(&s)),
                '\\' => s.push(self.escape(self.pos - 1)?),
                c => s.push(c),
            }
        }
    }

    /// The character an escape stands for; its backslash, at `start`, is
    /// consumed.
    fn escape(&mut self, start: usize) -> Result<char, FullFormError> {
        let invalid = FullFormError {
            position: start,
            kind: FullFormErrorKind::InvalidEscape,
        };
        let rest = &self.text[self.pos..];
        let (c, len) = match rest.chars().next() {
            Some('"') => ('"', 1),
            Some('\\') => ('\\', 1),
            Some('n') => ('\n', 1),
            Some('t') => ('\t', 1),
            Some('r') => ('\r', 1),
            Some('u') => {
                let close = rest.find('}').ok_or(invalid)?;
                let hex = rest[1..close].strip_prefix('{').ok_or(invalid)?;
                let valid =
                    (1..=6).contains(&hex.len()) && hex.chars().all(|c| c.is_ascii_hexdigit());
                let c = u32::from_str_radix(hex, 16)
                    .ok()
                    .filter(|_| valid)
                    .and_then(char::from_u32)
                    .ok_or(invalid)?;
                (c, close + 1)
            }
            _ => return Err(invalid),
        };
        self.pos += len;
        Ok(c)
    }
}

/// `head[args…]`; `Complex[re, im]` of two real numbers is the complex
/// number.
fn apply(head: &Expr, args: &[Expr]) -> Expr {
    if head.as_symbol() == Some(Symbol::COMPLEX) {
        if let [re, im] = args {
            if let (Some(re), Some(im)) = (re.to_rational(), im.to_rational()) {
                return Complex::new(re, im).into();
            }
        }
    }
    Expr::apply(head, args)
}

/// Whether `c` starts a symbol's name: a letter (Unicode `XID_Start`), `_`
/// or `$`.
fn is_name_start(c: char) -> bool {
    c == '_' || c == '$' || unicode_ident::is_xid_start(c)
}

/// Whether `c` continues a symbol's name: a letter, a digit, `_` (all
/// Unicode `XID_Continue`) or `$`.
fn is_name_continue(c: char) -> bool {
    c == '$' || unicode_ident::is_xid_continue(c)
}

/// Whether the reader reads `e`'s full form back as `e`: every symbol's name
/// is an identifier of the reader, and no application of `Complex` to two
/// real numbers looks like a complex number. For tests.
#[cfg(test)]
pub(crate) fn reads_back(e: &Expr) -> bool {
    use crate::expr::View;
    match e.view() {
        View::Symbol(s) => {
            let mut chars = s.name().chars();
            chars.next().is_some_and(is_name_start) && chars.all(is_name_continue)
        }
        View::Apply { head, args } => {
            let complex_like = head.as_symbol() == Some(Symbol::COMPLEX)
                && matches!(args, [re, im] if re.to_rational().is_some() && im.to_rational().is_some());
            !complex_like && reads_back(head) && args.iter().all(reads_back)
        }
        _ => true,
    }
}

#[cfg(test)]
mod tests;
