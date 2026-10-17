//! Core symbolic expression layer for the CAS engine.
//!
//! Immutable, hash-consed expressions in a uniform `head(args…)` form. See the
//! "Expression design" page in Lific (`CAS-DOC-3`).
//!
//! - [`Expr`]: the expression handle, 16 bytes; symbols and `i64` integers are
//!   inline, everything else is interned.
//! - [`Symbol`]: permanent, interned names.

mod expr;
mod intern;
mod node;
mod symbol;

pub use expr::{Expr, View};
pub use intern::{InternerStats, interner_stats};
pub use symbol::Symbol;

#[cfg(test)]
mod tests;
