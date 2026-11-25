//! The hook into the surface language (CAS-614): what natives need to turn
//! expressions into text and back, as `string`, `expr`, `cat` and `join` do,
//! and later `latex` and `mathml`. The language lives in cas-lang, which
//! depends on this crate, so the context holds it behind a trait
//! ([`Context::set_language`](crate::Context::set_language)).
//!
//! It is part of the context's state: a value evaluated under one hook may
//! print differently under another, so setting one gives a new state id and
//! nothing cached under one is reused under the other.

use crate::context::Context;
use crate::expr::Expr;

/// The forms [`Language::print`] writes.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Form {
    /// The input form: text that [`Language::read`] reads back.
    Input,
    /// LaTeX, as giac's `latex()` writes it (CAS-615).
    Latex,
}

/// A language that prints expressions as text and reads text back.
pub trait Language: Send + Sync {
    /// `e` as text in `form`, its balls shown at the working precision of
    /// `context`.
    ///
    /// # Errors
    ///
    /// The reason, as text, when the form cannot write `e` or the text
    /// would be longer than `max_len` bytes.
    fn print(
        &self,
        e: &Expr,
        form: Form,
        context: &Context,
        max_len: usize,
    ) -> Result<String, String>;

    /// The expression `text` writes, unevaluated, its decimals read at the
    /// working precision of `context`.
    ///
    /// # Errors
    ///
    /// The reason, as text, when `text` is not one expression.
    fn read(&self, text: &str, context: &Context) -> Result<Expr, String>;
}
