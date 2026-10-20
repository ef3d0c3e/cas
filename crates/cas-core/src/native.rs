//! The native functions of the builtins (CAS-DOC-4 §6, CAS-54): arithmetic
//! and the infinities, comparisons and predicates, logic, structure, and
//! what the language lowers to (CAS-DOC-7 §9, CAS-123): indexing, sets,
//! factorials and `Unapply`.
//!
//! Each is a [`NativeFn`] in a table indexed by builtin, so dispatching to
//! one is an array read. Each receives the application after the
//! evaluator normalised it, and returns its rewrite or `None`.

#![allow(
    clippy::unnecessary_wraps,
    reason = "natives share the `NativeFn` signature, which can fail"
)]

mod arith;
mod factorial;
mod index;
mod linear;
mod logic;
mod number;
mod radical;
mod replace;
mod structure;

use crate::build;
use crate::builtin;
use crate::eval::NativeFn;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The native functions of the builtins, by id: dispatching to one is an
/// array read (CAS-DOC-4 §6).
static NATIVES: [Option<NativeFn>; builtin::COUNT] = {
    let entries: &[(Symbol, NativeFn)] = &[
        (Symbol::PLUS, arith::plus),
        (Symbol::TIMES, arith::times),
        (Symbol::POWER, arith::power),
        (Symbol::SQRT, arith::sqrt),
        (Symbol::DIRECTED_INFINITY, arith::directed_infinity),
        (Symbol::EQUAL, logic::equal_native),
        (Symbol::UNEQUAL, logic::unequal),
        (Symbol::LESS, logic::less),
        (Symbol::LESS_EQUAL, logic::less_equal),
        (Symbol::GREATER, logic::greater),
        (Symbol::GREATER_EQUAL, logic::greater_equal),
        (Symbol::SAME_Q, logic::same_q),
        (Symbol::FREE_Q, logic::free_q),
        (Symbol::INTEGER_Q, logic::integer_q),
        (Symbol::RATIONAL_Q, logic::rational_q),
        (Symbol::REAL_Q, logic::real_q),
        (Symbol::COMPLEX_Q, logic::complex_q),
        (Symbol::NUMBER_Q, logic::number_q),
        (Symbol::LIST_Q, logic::list_q),
        (Symbol::MATRIX_Q, logic::matrix_q),
        (Symbol::STRING_Q, logic::string_q),
        (Symbol::SYMBOL_Q, logic::symbol_q),
        (Symbol::POSITIVE, logic::sign_predicate),
        (Symbol::NEGATIVE, logic::sign_predicate),
        (Symbol::NON_NEGATIVE, logic::sign_predicate),
        (Symbol::NON_POSITIVE, logic::sign_predicate),
        (Symbol::ELEMENT, logic::element),
        (Symbol::ABS, arith::abs),
        (Symbol::ELEMENTWISE_TIMES, linear::elementwise_times),
        (Symbol::ELEMENTWISE_DIVIDE, linear::elementwise_divide),
        (Symbol::ELEMENTWISE_POWER, linear::elementwise_power),
        (Symbol::FACTORIAL, factorial::factorial),
        (Symbol::AND, logic::and),
        (Symbol::OR, logic::or),
        (Symbol::NOT, logic::not),
        (Symbol::IF, logic::if_native),
        (Symbol::LENGTH, structure::length),
        (Symbol::PART, structure::part),
        (Symbol::MAP, structure::map),
        (Symbol::APPLY, structure::apply),
        (Symbol::FINITE_SET, structure::finite_set),
        (Symbol::UNAPPLY, structure::unapply),
        (Symbol::AT, index::at),
        (Symbol::REPLACE_AT, index::replace_at),
        (Symbol::REPLACE_ALL, replace::replace_all),
        (Symbol::REPLACE, replace::replace),
        (Symbol::REPLACE_REPEATED, replace::replace_repeated),
        (Symbol::MATCH_Q, replace::match_q),
    ];
    let mut table: [Option<NativeFn>; builtin::COUNT] = [None; builtin::COUNT];
    let mut i = 0;
    while i < entries.len() {
        table[entries[i].0.id() as usize] = Some(entries[i].1);
        i += 1;
    }
    table
};

/// The native function of the builtin `s`, if it has one.
pub(crate) fn get(s: Symbol) -> Option<NativeFn> {
    NATIVES.get(s.id() as usize).copied().flatten()
}

/// The value of a builtin symbol: `I` is the complex number `i` (CAS-51).
pub(crate) fn value(s: Symbol) -> Option<Expr> {
    (s == Symbol::I).then(build::i)
}

#[cfg(test)]
mod tests;
