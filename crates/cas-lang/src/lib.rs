//! The surface language of the CAS: one math-like grammar for students and
//! programs. The design is in Lific page `CAS-DOC-7`.
//!
//! - [`lex`]: text to lossless [`Lexed`] tokens of a [`SyntaxKind`], with
//!   their [`LexError`]s. Lexing never fails, and the tokens' texts
//!   concatenate back to the input, as the language server needs.
//! - [`parse`]: text to a lossless syntax tree ([`SyntaxNode`], a `rowan`
//!   tree) with its [`SyntaxError`]s. Parsing never fails either: what has no
//!   place goes into error nodes, and recovery stays local to a line or a
//!   bracket.
//! - [`number_value`], [`string_value`], [`superscript_value`] and [`name`]:
//!   the values of literal tokens. Numbers are exact (`0.1` is `1/10`), and
//!   names are in Unicode normal form C.
//! - [`Names`]: the names the language reserves and what they stand for
//!   ([`Target`]), with the spellings printers use; English plus at most one
//!   locale, both data. User names that start with an ASCII capital are
//!   escaped ([`escape_name`]) so they never meet a builtin's.
//! - [`lower`]: the statements of a tree as expressions ([`Lowered`]):
//!   values to evaluate, or statement heads for the session (`Set`,
//!   `SetDelayed`, …), with [`Spans`] back into the text and
//!   [`Diagnostic`]s, some with a [`Fix`].
//! - [`print()`]: an expression as text that reads back as it (the input
//!   form), for output that can be pasted as input; a [`PrintError`] for
//!   what the language cannot write.
//! - [`Session`]: a program run one statement at a time against a
//!   `Context` (CAS-121): definitions, `assume` and `purge`, results in
//!   input form, and [`Message`]s that point into the text.

mod lexer;
mod literal;
mod lower;
mod names;
mod parser;
mod print;
mod session;
mod syntax;

pub use lexer::{LexError, LexErrorKind, Lexed, lex};
pub use literal::{
    MAX_DECIMAL_EXPONENT, Number, NumberError, decimal_balls, name, number_value, string_value,
    superscript_value,
};
pub use lower::{Diagnostic, Fix, LowerErrorKind, Lowered, Severity, Spans, Statement, lower};
pub use names::{
    Names, NamesError, NamesErrorKind, SyntaxForm, Target, TypeName, escape_name, unescape_name,
};
pub use parser::{MAX_DEPTH, Parse, ParseErrorKind, SyntaxError, SyntaxErrorKind, parse};
#[doc(hidden)]
pub use print::fold_arithmetic;
pub use print::{Balls, MAX_PRINT_DEPTH, PrintError, PrintErrorKind, print, print_with};
pub use session::{Level, Message, Output, Session, Shown};
pub use syntax::{CasLanguage, SyntaxElement, SyntaxKind, SyntaxNode, SyntaxToken};
