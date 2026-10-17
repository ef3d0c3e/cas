//! Core symbolic expression layer for the CAS engine.
//!
//! Immutable, hash-consed expressions in a uniform `head(args…)` form. See the
//! "Expression design" page in Lific (`CAS-DOC-3`).
//!
//! - [`Expr`]: the expression handle, 16 bytes; symbols and `i64` integers are
//!   inline, everything else is interned.
//! - [`Symbol`]: permanent, interned names; the builtins (`Symbol::PLUS`, …)
//!   have fixed ids and fixed [`Attributes`].
//! - [`Context`]: values, user attributes and assumptions; scopes are child
//!   values. Evaluation goes through a context.

mod attributes;
mod builtin;
mod context;
mod expr;
mod intern;
mod node;
mod symbol;

pub use attributes::Attributes;
pub use context::{Context, ContextError, EvalError, Limits, StateId};
pub use expr::{Expr, View};
pub use intern::{InternerStats, interner_stats};
pub use symbol::Symbol;

#[cfg(test)]
mod tests;
