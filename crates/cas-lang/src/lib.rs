//! The surface language of the CAS: one math-like grammar for students and
//! programs. The design is in Lific page `CAS-DOC-7`.
//!
//! - [`lex`]: text to lossless [`Lexed`] tokens of a [`SyntaxKind`], with
//!   their [`LexError`]s. Lexing never fails, and the tokens' texts
//!   concatenate back to the input, as the language server needs.
//! - [`number_value`], [`string_value`], [`superscript_value`] and [`name`]:
//!   the values of literal tokens. Numbers are exact (`0.1` is `1/10`), and
//!   names are in Unicode normal form C.
//! - [`Names`]: the names the language reserves and what they stand for
//!   ([`Target`]), with the spellings printers use; English plus at most one
//!   locale, both data. User names that start with an ASCII capital are
//!   escaped ([`escape_name`]) so they never meet a builtin's.
//!
//! The parser, lowering to expressions and the session come with CAS-118,
//! CAS-120 and CAS-121.

mod lexer;
mod literal;
mod names;
mod syntax;

pub use lexer::{LexError, LexErrorKind, Lexed, lex};
pub use literal::{
    MAX_DECIMAL_EXPONENT, NumberError, name, number_value, string_value, superscript_value,
};
pub use names::{
    Names, NamesError, NamesErrorKind, SyntaxForm, Target, TypeName, escape_name, unescape_name,
};
pub use syntax::SyntaxKind;
