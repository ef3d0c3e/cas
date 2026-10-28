//! Comparisons, predicates, `And`, `Or`, `Not` and `If` (CAS-DOC-4 §6.1,
//! §7).
//!
//! Comparisons decide on numbers, on the infinities `Infinity` and
//! `-Infinity`, and on expressions that are structurally equal; on balls as
//! [`ball_relation`] says; otherwise
//! the assumptions and the number domain decide them when they can (`assume`,
//! CAS-58), then, for exact numeric quantities, approximations
//! ([`approximated`], CAS-285), and they stay unevaluated when not. So do `Positive`,
//! `Negative`, `NonNegative`, `NonPositive` and `Element`. The type tests
//! (`SameQ`, `FreeQ`, `IntegerQ` to `SymbolQ`) always decide; those of number
//! sets are true of a value the assumptions put in the set.

use std::cmp::Ordering;

use cas_num::Rational;
use cas_num::ball::{Ball, ComplexBall};

use super::number::{Infinity, Num};
use super::rounding::real_sign_within;
use crate::assume::{self, Set, Truth};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

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

/// How two numbers stand when one is a ball (CAS-DOC-5 §9.1, revised
/// 2026-10-20, D3).
enum Relation {
    /// They overlap, and both are accurate to within [`SLACK`] bits of the
    /// working precision: equal, as in Xcas and Mathematica.
    Equal,
    /// They are disjoint, ordered when both are real.
    Apart(Option<Ordering>),
    /// They overlap, and one has lost accuracy: the honest "don't know".
    Unknown,
}

/// The bits below the working precision that an accurate ball may lack:
/// Mathematica ignores the last 7 for `Equal`.
const SLACK: i64 = 8;

/// Whether a ball knows its value to within [`SLACK`] bits of `prec`; only
/// an exact one does when comparisons are `strict`, as inside `N`.
fn accurate(b: &Ball, prec: u64, strict: bool) -> bool {
    b.is_exact()
        || !strict && b.rel_accuracy_bits() >= i64::try_from(prec).unwrap_or(i64::MAX) - SLACK
}

/// The relation of overlapping or disjoint values: `Equal` when they
/// overlap and both are accurate.
fn relation(overlap: bool, accurate: bool, order: Option<Ordering>) -> Relation {
    if !overlap {
        Relation::Apart(order)
    } else if accurate {
        Relation::Equal
    } else {
        Relation::Unknown
    }
}

/// [`Relation`] of `a` and `b` if one is a ball and the other a number, at
/// the working precision of `ctx`, and strictly there if its comparisons
/// are (only exact balls may be equal).
fn ball_relation(a: &Expr, b: &Expr, ctx: &Context) -> Option<Relation> {
    // Most comparisons have no ball: no values to make for them.
    let ball = |e: &Expr| matches!(e.view(), View::Ball(_) | View::ComplexBall(_));
    if !ball(a) && !ball(b) {
        return None;
    }
    let (prec, strict) = (ctx.working_precision(), ctx.strict_comparisons());
    let (left, right) = (Num::of(a)?, Num::of(b)?);
    let close = |q: &Rational| Ball::from_rational(q, prec + 64);
    Some(match (&left, &right) {
        (Num::Ball(x), Num::Ball(y)) => relation(
            x.overlaps(y),
            accurate(x, prec, strict) && accurate(y, prec, strict),
            x.cmp(y),
        ),
        (Num::Ball(x), Num::Real(q)) => relation(
            x.contains_rational(q),
            accurate(x, prec, strict),
            x.cmp(&close(q)),
        ),
        (Num::Real(q), Num::Ball(y)) => relation(
            y.contains_rational(q),
            accurate(y, prec, strict),
            close(q).cmp(y),
        ),
        _ => complex_relation(&left, &right, prec, strict),
    })
}

/// Whether `a` and `b` are numbers, one of them a ball, that overlap
/// without being equal for sure, with the strict comparisons of `ctx`:
/// more precision may tell them apart, which `N` tries (CAS-DOC-5 §9.4).
fn undecided_balls(a: &Expr, b: &Expr, ctx: &Context) -> bool {
    matches!(ball_relation(a, b, ctx), Some(Relation::Unknown))
}

/// Whether a ball that is not exact occurs in `e`.
fn has_inexact_ball(e: &Expr) -> bool {
    match e.view() {
        View::Ball(b) => !b.is_exact(),
        View::ComplexBall(z) => !z.re().is_exact() || !z.im().is_exact(),
        View::Apply { head, args } => has_inexact_ball(head) || args.iter().any(has_inexact_ball),
        _ => false,
    }
}

/// [`ball_relation`] with a complex value: equal only part by part, never
/// ordered.
fn complex_relation(left: &Num, right: &Num, prec: u64, strict: bool) -> Relation {
    // Each side as a complex ball, and whether it is accurate.
    let side = |n: &Num| -> (ComplexBall, bool) {
        match n {
            Num::ComplexBall(z) => (
                ComplexBall::clone(z),
                accurate(z.re(), prec, strict) && accurate(z.im(), prec, strict),
            ),
            Num::Ball(x) => (
                ComplexBall::from_real(Ball::clone(x)),
                accurate(x, prec, strict),
            ),
            Num::Complex(z) => (ComplexBall::from_complex(z, prec + 64), true),
            Num::Real(q) => (
                ComplexBall::from_real(Ball::from_rational(q, prec + 64)),
                true,
            ),
        }
    };
    let ((l, l_accurate), (r, r_accurate)) = (side(left), side(right));
    let overlap = match (left, right) {
        (Num::Complex(z), _) => r.contains_complex(z),
        (_, Num::Complex(z)) => l.contains_complex(z),
        _ => l.overlaps(&r),
    };
    relation(overlap, l_accurate && r_accurate, None)
}

/// Whether `a` and `b` are equal in `ctx`, when it can be told:
/// structurally equal expressions are; numbers, strings, booleans and
/// infinities that differ are not; balls as [`ball_relation`] says. With
/// strict comparisons, equal expressions with an inexact ball are not equal
/// for sure: different values may have made the same ball.
fn equal(a: &Expr, b: &Expr, ctx: &Context) -> Option<bool> {
    if a == b && !(ctx.strict_comparisons() && has_inexact_ball(a)) {
        return Some(true);
    }
    if let Some(relation) = ball_relation(a, b, ctx) {
        return match relation {
            Relation::Equal => Some(true),
            Relation::Apart(_) => Some(false),
            Relation::Unknown => None,
        };
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

/// The bits [`approximated`] goes to at most, as Mathematica's
/// `$MaxExtraPrecision` (50 digits) past machine precision: closer values
/// are `N`'s to tell apart, raising the precision with its warnings
/// (CAS-DOC-5 §9.4).
const MAX_COMPARISON_BITS: u64 = 256;

/// Whether `a head b` holds for exact numeric quantities, `head` a
/// comparison, by the sign of `a − b` as `Sign` decides it, from
/// approximations at increasing precision up to [`MAX_COMPARISON_BITS`]
/// (CAS-285): `sqrt(2) < 3/2`, as giac's. `None` for other values, inexact
/// balls (which [`ball_relation`] rules), complex values, values closer
/// than the precision tells, and equal quantities in different forms,
/// which no approximation tells apart (`ln(1/2)` and `-ln(2)`).
fn approximated(
    ev: &mut Eval<'_>,
    head: Symbol,
    a: &Expr,
    b: &Expr,
) -> Result<Option<bool>, EvalError> {
    let exact = |e: &Expr| e.is_numeric() && !has_inexact_ball(e);
    if !exact(a) || !exact(b) {
        return Ok(None);
    }
    let difference = ev.eval(&(a.clone() - b.clone()))?;
    let Some(order) = real_sign_within(ev.context(), &difference, MAX_COMPARISON_BITS) else {
        return Ok(None);
    };
    Ok(match head {
        Symbol::LESS => Some(order.is_lt()),
        Symbol::LESS_EQUAL => Some(order.is_le()),
        Symbol::GREATER => Some(order.is_gt()),
        Symbol::GREATER_EQUAL => Some(order.is_ge()),
        Symbol::EQUAL => Some(order.is_eq()),
        Symbol::UNEQUAL => Some(order.is_ne()),
        _ => None,
    })
}

/// `true` or `false` for a decided truth.
fn decided(t: Truth) -> Option<bool> {
    match t {
        Truth::True => Some(true),
        Truth::False => Some(false),
        Truth::Unknown => None,
    }
}

/// A chain of tests between neighbours: `True` if all hold, `False` if one
/// fails, unevaluated otherwise. What `test` cannot tell, the assumptions
/// may. With strict comparisons, balls left undecided give
/// [`Warning::Undecided`], for `N` to raise the precision.
fn chain(ev: &mut Eval<'_>, e: &Expr, test: impl Fn(&Expr, &Expr) -> Option<bool>) -> Native {
    let head = e.head().and_then(Expr::as_symbol).expect("a comparison");
    let mut all = true;
    let mut undecided = false;
    for w in e.args().windows(2) {
        let result = match test(&w[0], &w[1]) {
            Some(b) => Some(b),
            None => match decided(assume::relation(ev, head, &w[0], &w[1])?) {
                Some(b) => Some(b),
                None => approximated(ev, head, &w[0], &w[1])?,
            },
        };
        match result {
            Some(false) => return Ok(Some(boolean(false))),
            Some(true) => {}
            None => {
                all = false;
                let ctx = ev.context();
                undecided |= ctx.strict_comparisons() && undecided_balls(&w[0], &w[1], ctx);
            }
        }
    }
    if undecided {
        warn_undecided(ev, e);
    }
    Ok(all.then(|| boolean(true)))
}

/// Says that the balls of `comparison` overlap too much to decide it; out
/// of line, as only `N` meets it.
#[cold]
#[inline(never)]
fn warn_undecided(ev: &mut Eval<'_>, comparison: &Expr) {
    let bits = ev.context().working_precision();
    ev.warn(Warning::Undecided {
        comparison: comparison.clone(),
        bits,
    });
}

pub(super) fn equal_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let ctx = ev.context();
    chain(ev, e, |a, b| equal(a, b, ctx))
}

pub(super) fn unequal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let ctx = ev.context();
    if let Some(eq) = equal(a, b, ctx) {
        return Ok(Some(boolean(!eq)));
    }
    let result = match decided(assume::relation(ev, Symbol::UNEQUAL, a, b)?) {
        Some(b) => Some(b),
        None => approximated(ev, Symbol::UNEQUAL, a, b)?,
    }
    .map(boolean);
    if result.is_none() && ctx.strict_comparisons() && undecided_balls(a, b, ctx) {
        warn_undecided(ev, e);
    }
    Ok(result)
}

/// `Less`, `LessEqual`, `Greater` or `GreaterEqual`, by the orderings they
/// accept; balls agree with `Equal` (D3): accurate overlapping balls are
/// equal, disjoint ones ordered, the others undecided.
fn order(ev: &mut Eval<'_>, e: &Expr, accept: fn(Ordering) -> bool) -> Native {
    let ctx = ev.context();
    chain(ev, e, |a, b| {
        if let Some(relation) = ball_relation(a, b, ctx) {
            return match relation {
                Relation::Equal => Some(accept(Ordering::Equal)),
                Relation::Apart(Some(o)) => Some(accept(o)),
                Relation::Apart(None) | Relation::Unknown => None,
            };
        }
        Some(accept(extended(a)?.cmp(&extended(b)?)))
    })
}

pub(super) fn less(ev: &mut Eval<'_>, e: &Expr) -> Native {
    order(ev, e, Ordering::is_lt)
}

pub(super) fn less_equal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    order(ev, e, Ordering::is_le)
}

pub(super) fn greater(ev: &mut Eval<'_>, e: &Expr) -> Native {
    order(ev, e, Ordering::is_gt)
}

pub(super) fn greater_equal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    order(ev, e, Ordering::is_ge)
}

/// `Positive`, `Negative`, `NonNegative` and `NonPositive`: by the sign of
/// their argument, which the assumptions may tell. A number that is not
/// real has none of these signs.
pub(super) fn sign_predicate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let head = e.head().and_then(Expr::as_symbol).expect("a predicate");
    let sign = assume::sign_of(ev.context(), x);
    Ok(decided(assume::sign_test(head, sign, x.is_exact_number())).map(boolean))
}

/// `Element[x, s]` for the number sets `Integers`, `Rationals`, `Reals` and
/// `Complexes`.
pub(super) fn element(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, s] = e.args() else {
        return Ok(None);
    };
    let Some(set) = Set::of(s) else {
        return Ok(None);
    };
    Ok(decided(assume::element(ev.context(), x, set)).map(boolean))
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

/// A type test of a number set: whether the value is known to be in it.
fn in_set(ev: &mut Eval<'_>, e: &Expr, set: Set) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(
        assume::element(ev.context(), x, set) == Truth::True,
    )))
}

pub(super) fn integer_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    in_set(ev, e, Set::Integers)
}

pub(super) fn rational_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    in_set(ev, e, Set::Rationals)
}

pub(super) fn real_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    in_set(ev, e, Set::Reals)
}

pub(super) fn complex_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    in_set(ev, e, Set::Complexes)
}

pub(super) fn number_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(x.is_number())))
}

/// A structural type test.
fn is(e: &Expr, test: impl Fn(&Expr) -> bool) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(boolean(test(x))))
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

pub(super) fn list_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    is(e, is_list)
}

/// `MatrixQ[m]`: a list of one or more rows, lists of one length, whose
/// entries are not lists.
pub(super) fn matrix_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    is(e, |m| {
        is_list(m)
            && !m.args().is_empty()
            && m.args().iter().all(|row| {
                is_list(row)
                    && row.args().len() == m.args()[0].args().len()
                    && !row.args().iter().any(is_list)
            })
    })
}

pub(super) fn string_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    is(e, |x| x.as_str().is_some())
}

pub(super) fn symbol_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    is(e, |x| x.as_symbol().is_some())
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
