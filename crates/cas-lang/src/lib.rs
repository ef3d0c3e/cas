//! The surface language of the CAS: one math-like grammar for students and
//! programs. The design is in Lific page `CAS-DOC-7`.
//!
//! - [`lex`]: text to lossless [`Lexed`] tokens of a [`SyntaxKind`], with
//!   their [`LexError`]s. Lexing never fails, and the tokens' texts
//!   concatenate back to the input, as the language server needs.
//! - [`number_value`], [`string_value`], [`superscript_value`] and [`name`]:
//!   the values of literal tokens. Numbers are exact (`0.1` is `1/10`), and
//!   names are in Unicode normal form C.
//!
//! The parser, the name table, lowering to expressions and the session come
//! with CAS-118 to CAS-121.

mod lexer;
mod literal;
mod syntax;

pub use lexer::{LexError, LexErrorKind, Lexed, lex};
pub use literal::{
    MAX_DECIMAL_EXPONENT, NumberError, name, number_value, string_value, superscript_value,
};
pub use syntax::SyntaxKind;
