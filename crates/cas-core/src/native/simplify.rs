//! `Simplify` (CAS-DOC-14 §1–§3, CAS-229): Xcas's `simplify`, the least
//! complex form found of a value.
//!
//! - The complexity (§2), [`complexity`]: a symbol or a string 1; an integer
//!   its decimal digits (1 for 0), plus 1 when negative; a rational or a
//!   complex number its parts', plus 1; a ball 2; an application its head's
//!   and its arguments'.
//! - The search (§3): the arguments of functions, and the base and exponent
//!   of powers but integer ones, first, each by a `Simplify` of its own;
//!   then the candidates `Together`, `Factor` and `Expand` of the whole, and
//!   the least complex of it and them, the first on a tie, in at most
//!   [`ROUNDS`] rounds while the choice changes. A candidate is evaluated on
//!   its own: one that fails, short of an interrupt, is left out, and the
//!   warnings of all are.
//! - Lists, the sides of equations and comparisons, and the arguments of
//!   `And`, `Or` and `Not` are simplified each; a head that holds its
//!   arguments keeps them as they are.
//!
//! One step a call, `"simplify"`.

use crate::attributes::Attributes;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The rounds of candidates at most.
const ROUNDS: usize = 4;

/// The heads simplified argument by argument.
const THREADED: [Symbol; 11] = [
    Symbol::LIST,
    Symbol::EQUATION,
    Symbol::EQUAL,
    Symbol::UNEQUAL,
    Symbol::LESS,
    Symbol::LESS_EQUAL,
    Symbol::GREATER,
    Symbol::GREATER_EQUAL,
    Symbol::AND,
    Symbol::OR,
    Symbol::NOT,
];

/// The rewritings tried as candidates.
const CANDIDATES: [Symbol; 3] = [Symbol::TOGETHER, Symbol::FACTOR, Symbol::EXPAND];

/// `Simplify[e]`.
pub(super) fn simplify(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = simplified(ev, x)?;
    ev.step("simplify", e, &result);
    Ok(Some(result))
}

/// `Simplify[e]` evaluated, a level of its own.
fn nested(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    if e.head().is_none() {
        return Ok(e.clone());
    }
    ev.eval(&Expr::call(Symbol::SIMPLIFY, [e.clone()]))
}

/// `e`, evaluated, simplified.
fn simplified(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let Some(head) = e.head() else {
        return Ok(e.clone());
    };
    if head.as_symbol().is_some_and(|s| THREADED.contains(&s)) {
        let mut args = Vec::with_capacity(e.args().len());
        for a in e.args() {
            args.push(nested(ev, a)?);
        }
        return ev.eval(&Expr::call(head.clone(), args));
    }
    let mut best = inside(ev, e)?;
    for _ in 0..ROUNDS {
        let next = least(ev, &best)?;
        if next == best {
            break;
        }
        best = next;
    }
    Ok(best)
}

/// `e` with the arguments of its functions, and the base and exponent of
/// its powers but integer ones, simplified, and evaluated again; the
/// arguments of `Plus`, `Times` and integer powers looked into the same way.
fn inside(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let View::Apply { head, args } = e.view() else {
        return Ok(e.clone());
    };
    let mut new = Vec::with_capacity(args.len());
    match (head.as_symbol(), args) {
        (Some(Symbol::PLUS | Symbol::TIMES), _) => {
            for a in args {
                new.push(inside(ev, a)?);
            }
        }
        (Some(Symbol::POWER), [b, n]) if n.to_integer().is_some() => {
            new.push(inside(ev, b)?);
            new.push(n.clone());
        }
        (Some(s), _) if holds(ev.context().attributes(s)) => return Ok(e.clone()),
        _ => {
            for a in args {
                new.push(nested(ev, a)?);
            }
        }
    }
    if new.as_slice() == args {
        return Ok(e.clone());
    }
    ev.eval(&Expr::call(head.clone(), new))
}

/// Whether a head with `attributes` holds an argument.
fn holds(attributes: Attributes) -> bool {
    attributes.contains(Attributes::HOLD_FIRST) || attributes.contains(Attributes::HOLD_REST)
}

/// The least complex of `e` and its candidates, the first on a tie.
fn least(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let mut best = e.clone();
    let mut cost = complexity(e);
    for head in CANDIDATES {
        let Some(c) = candidate(ev, &Expr::call(head, [e.clone()]))? else {
            continue;
        };
        let k = complexity(&c);
        if k < cost {
            best = c;
            cost = k;
        }
    }
    Ok(best)
}

/// `call` evaluated on its own, without a trace or warnings; `None` when
/// that fails short of an interrupt, or leaves the call.
fn candidate(ev: &mut Eval<'_>, call: &Expr) -> Result<Option<Expr>, EvalError> {
    let context = ev.context();
    match ev.eval_in(context, call).0 {
        Ok(c) if c.head() == call.head() => Ok(None),
        Ok(c) => Ok(Some(c)),
        Err(EvalError::Interrupted) => Err(EvalError::Interrupted),
        Err(_) => Ok(None),
    }
}

/// The complexity of `e` (§2).
pub(super) fn complexity(e: &Expr) -> u64 {
    let mut total: u64 = 0;
    let mut stack = vec![e];
    while let Some(e) = stack.pop() {
        let own = match e.view() {
            View::Symbol(_) | View::String(_) => 1,
            View::Int(n) => small_digits(n),
            View::Integer(n) => digits(n.bit_length()) + u64::from(n.is_negative()),
            View::Rational(q) => rational(q),
            View::Complex(z) => rational(z.re()) + rational(z.im()) + 1,
            View::Ball(_) | View::ComplexBall(_) => 2,
            View::Apply { head, args } => {
                stack.push(head);
                stack.extend(args);
                0
            }
        };
        total = total.saturating_add(own);
    }
    total
}

/// The decimal digits of an `i64`, 1 for 0, plus 1 when negative.
fn small_digits(n: i64) -> u64 {
    let digits = n
        .unsigned_abs()
        .checked_ilog10()
        .map_or(1, |d| u64::from(d) + 1);
    digits + u64::from(n < 0)
}

/// The decimal digits of an integer of `bits` bits, about.
fn digits(bits: u64) -> u64 {
    bits.saturating_mul(30_103) / 100_000 + 1
}

/// The complexity of a rational: its parts', plus 1 if not an integer.
fn rational(q: &cas_num::Rational) -> u64 {
    let integer = |n: &cas_num::Integer| match n.to_i64() {
        Some(n) => small_digits(n),
        None => digits(n.bit_length()) + u64::from(n.is_negative()),
    };
    if q.is_integer() {
        return integer(q.numer());
    }
    integer(q.numer()) + integer(q.denom()) + 1
}

#[cfg(test)]
mod tests;
