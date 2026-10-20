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
//!   `Debug` or [`Expr::full_form`] (`Plus[x, Times[2, Power[y, 3]]]`), which
//!   [`Expr::parse_full_form`] reads back; [`assert_eval!`] states evaluator
//!   tests in full form.
//!   Traversal (`part`, `positions`, `pre_order`, `distinct`, `contains`,
//!   `replace`, `substitute`, `map_bottom_up`) is structural and memoised by
//!   node identity, so it is linear in the distinct subexpressions.
//! - [`Symbol`]: permanent, interned names; the builtins (`Symbol::PLUS`, …)
//!   have fixed ids and fixed [`Attributes`].
//! - [`build`]: operators (`x + 2 * y.pow(3)`), [`Expr::call`] and the
//!   builtin functions, for building expressions from Rust without a context.
//! - [`Context`]: values, user attributes and assumptions; scopes are child
//!   values. [`Context::eval`] evaluates in a context (CAS-DOC-4 §2): own
//!   values, holds, normalisation of `Sequence`, `Flat`, `Listable` and
//!   `Orderless` heads, and native functions ([`NativeFn`], which receive
//!   the call's [`Eval`]: exact arithmetic with radicals and infinities,
//!   comparisons, logic, structure, `Replace`), to a fixpoint that is
//!   stamped, so evaluating it again is O(1). [`Pattern`] and [`Rule`] are
//!   compiled once and matched by backtracking. [`Limits`] bound recursion and rewrites;
//!   [`Limits::deep`] evaluates on a thread of its own with a large stack.
//!   [`Context::eval_with`] also reports [`Warning`]s and, as [`EvalOptions`]
//!   ask, a [`Trace`] of the steps taken, which [`Trace::verify`] checks at
//!   exact sample points.

mod assume;
mod attributes;
pub mod build;
mod builtin;
mod calculus;
mod context;
mod eval;
mod expr;
mod intern;
mod native;
mod node;
mod order;
mod pattern;
mod print;
mod read;
mod rules;
mod symbol;
mod symbol_map;
mod trace;
mod traverse;

pub use assume::{Domain, Sign, Truth};
pub use attributes::Attributes;
pub use cas_num::{Complex, Integer, Rational};
pub use context::{Context, ContextError, EvalError, Interrupt, Limits, StateId};
pub use eval::{Eval, NativeFn, STACK_PER_LEVEL};
pub use expr::{Expr, View};
pub use intern::{InternerStats, interner_stats};
pub use pattern::{Pattern, PatternError, PatternErrorKind, Rule};
pub use print::FullForm;
pub use read::{FullFormError, FullFormErrorKind};
pub use rules::{DefineError, UserNative};
pub use symbol::Symbol;
pub use trace::{
    Counterexample, EvalOptions, Origin, Report, RuleId, RuleKind, Step, StepKind, Trace,
    TraceLevel, Verification, Warning,
};
pub use traverse::{Distinct, PostOrder, PreOrder};

#[cfg(test)]
mod tests;
