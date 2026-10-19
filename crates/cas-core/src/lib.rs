//! Core symbolic expression layer for the CAS engine.
//!
//! Immutable, hash-consed expressions in a uniform `head(args…)` form. See the
//! "Expression design" page in Lific (`CAS-DOC-3`).
//!
//! - [`Expr`]: the expression handle, 16 bytes; symbols and `i64` integers are
//!   inline, everything else (larger [`Integer`]s, [`Rational`]s, strings,
//!   applications) is interned. `Ord` is the canonical order that
//!   evaluation sorts `Orderless` arguments into (see `order.rs`).
//!   It prints in infix with `Display` (`x + 2*y^3`) and in full form with
//!   `Debug` or [`Expr::full_form`] (`Plus[x, Times[2, Power[y, 3]]]`).
//!   Traversal (`part`, `positions`, `pre_order`, `distinct`, `contains`,
//!   `replace`, `substitute`, `map_bottom_up`) is structural and memoised by
//!   node identity, so it is linear in the distinct subexpressions.
//! - [`Symbol`]: permanent, interned names; the builtins (`Symbol::PLUS`, …)
//!   have fixed ids and fixed [`Attributes`].
//! - [`build`]: operators (`x + 2 * y.pow(3)`), [`Expr::call`] and the
//!   builtin functions, for building expressions from Rust without a context.
//! - [`Context`]: values, user attributes and assumptions; scopes are child
//!   values. Evaluation goes through a context.

mod attributes;
pub mod build;
mod builtin;
mod context;
mod expr;
mod intern;
mod node;
mod order;
mod print;
mod symbol;
mod traverse;

pub use attributes::Attributes;
pub use cas_num::{Complex, Integer, Rational};
pub use context::{Context, ContextError, EvalError, Limits, StateId};
pub use expr::{Expr, View};
pub use intern::{InternerStats, interner_stats};
pub use print::FullForm;
pub use symbol::Symbol;
pub use traverse::{Distinct, PostOrder, PreOrder};

#[cfg(test)]
mod tests;
