//! Comparisons, predicates, `And`, `Or`, `Not` and `If` (CAS-DOC-4 §6.1).
//!
//! Comparisons decide on numbers, on the infinities `Infinity` and
//! `-Infinity`, and on expressions that are structurally equal; anything
//! else stays unevaluated until assumptions can decide it (CAS-58).
//! Predicates (`SameQ`, `FreeQ`, `IntegerQ`, `NumberQ`) always decide.

use std::cmp::Ordering;

use cas_num::Rational;

use super::number::{Infinity, Num};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn boolean(b: bool) -> Expr {
    Expr::from(if b { Symbol::TRUE } else { Symbol::FALSE })
}

fn truth(e: &Expr) -> Option<bool> {
    match e.as_symbol() {
        Some(Symbol::TRUE) => Some(true),
        Some(Symbol::FALSE) => Some(false),
        _ => None,
    }
}

/// Whether `a` and `b` are equal, when it can be told: structurally equal
/// expressions are; numbers, strings, booleans and infinities that differ are
/// not.
fn equal(a: &Expr, b: &Expr) -> Option<bool> {
    if a == b {
        return Some(true);
    }
    let known = |e: &Expr| {
        e.is_number()
            || e.as_str().is_some()
            || truth(e).is_some()
            || matches!(Infinity::of(e), Some(Infinity::Directed(_)))
    };
    (known(a) && known(b)).then_some(false)
}

/// A point of the extended real line, for ordering.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum Extended {
    NegativeInfinity,
    Real(Rational),
    Infinity,
}

fn extended(e: &Expr) -> Option<Extended> {
    if let Some(q) = e.to_rational() {
        return Some(Extended::Real(q));
    }
    match Infinity::of(e)? {
        Infinity::Directed(Num::Real(d)) if d.is_one() => Some(Extended::Infinity),
        Infinity::Directed(Num::Real(d)) if (-&d).is_one() => Some(Extended::NegativeInfinity),
        _ => None,
    }
}

/// A chain of tests between neighbours: `True` if all hold, `False` if one
/// fails, unevaluated otherwise.
fn chain(e: &Expr, test: impl Fn(&Expr, &Expr) -> Option<bool>) -> Native {
    let mut decided = true;
    for w in e.args().windows(2) {
        match test(&w[0], &w[1]) {
            Some(false) => return Ok(Some(boolean(false))),
            Some(true) => {}
            None => decided = false,
        }
    }
    Ok(decided.then(|| boolean(true)))
}

pub(super) fn equal_native(_: &mut Eval<'_>, e: &Expr) -> Native {
    chain(e, equal)
}

pub(super) fn unequal(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    Ok(equal(a, b).map(|eq| boolean(!eq)))
}

/// `Less`, `LessEqual`, `Greater` or `GreaterEqual`, by the orderings they
/// accept.
fn order(e: &Expr, accept: fn(Ordering) -> bool) -> Native {
    chain(e, |a, b| Some(accept(extended(a)?.cmp(&extended(b)?))))
}

pub(super) fn less(_: &mut Eval<'_>, e: &Expr) -> Native {
    order(e, Ordering::is_lt)
}

pub(super) fn less_equal(_: &mut Eval<'_>, e: &Expr) -> Native {
    order(e, Ordering::is_le)
}

pub(super) fn greater(_: &mut Eval<'_>, e: &Expr) -> Native {
    order(e, Ordering::is_gt)
}

pub(super) fn greater_equal(_: &mut Eval<'_>, e: &Expr) -> Native {
    order(e, Ordering::is_ge)
}

/// `SameQ[a, b, …]`: structural equality, O(1) per pair.
pub(super) fn same_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(Some(boolean(e.args().windows(2).all(|w| w[0] == w[1]))))
}

/// `FreeQ[e, x]`: whether `x` occurs nowhere in `e`, heads included.
pub(super) fn free_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [expr, x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(expr.free_of(x))))
}

pub(super) fn integer_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(matches!(
        x.view(),
        View::Int(_) | View::Integer(_)
    ))))
}

pub(super) fn number_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(x.is_number())))
}

/// `And` (`HoldAll`): its arguments are evaluated in order until one is
/// `False`; the `True` ones are dropped. `And[]` is `True` and `And[p]` is
/// `p`.
pub(super) fn and(ev: &mut Eval<'_>, e: &Expr) -> Native {
    connective(ev, e, false, Symbol::AND)
}

/// `Or` (`HoldAll`): the dual of [`and`].
pub(super) fn or(ev: &mut Eval<'_>, e: &Expr) -> Native {
    connective(ev, e, true, Symbol::OR)
}

/// `And` or `Or`: `stop` decides it, the other value is dropped.
fn connective(ev: &mut Eval<'_>, e: &Expr, stop: bool, head: Symbol) -> Native {
    let mut rest = Vec::new();
    for a in e.args() {
        let value = ev.eval(a)?;
        match truth(&value) {
            Some(b) if b == stop => return Ok(Some(boolean(stop))),
            Some(_) => {}
            None => rest.push(value),
        }
    }
    let result = match rest.len() {
        0 => boolean(!stop),
        1 => rest.pop().expect("one argument"),
        _ => Expr::apply(&head.into(), &rest),
    };
    Ok((result != *e).then_some(result))
}

/// `Not[True]` is `False`, `Not[False]` is `True`, `Not[Not[p]]` is `p`.
pub(super) fn not(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    if let Some(b) = truth(p) {
        return Ok(Some(boolean(!b)));
    }
    if p.head().and_then(Expr::as_symbol) == Some(Symbol::NOT) {
        if let [q] = p.args() {
            return Ok(Some(q.clone()));
        }
    }
    Ok(None)
}

/// `If[c, a]` and `If[c, a, b]` (`HoldRest`): the branch `c` chooses,
/// evaluated after, or `Null` without an `else`; unevaluated while `c` is
/// neither `True` nor `False`.
pub(super) fn if_native(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [c, a] => Ok(truth(c).map(|b| if b { a.clone() } else { Symbol::NULL.into() })),
        [c, a, b] => Ok(truth(c).map(|t| if t { a.clone() } else { b.clone() })),
        _ => Ok(None),
    }
}
