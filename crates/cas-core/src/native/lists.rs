//! List operations (CAS-DOC-28, CAS-334 to CAS-336), Xcas's, next to
//! `Length`, `Part`, `Map` and `Apply` (CAS-DOC-4 §6.1):
//!
//! - **Building**: `Table[f, k, a, b]`, `Table[f, k, a, b, s]` and
//!   `Table[f, k == a .. b]` (`seq`), `f` evaluated at each `k = a + js`
//!   within `b`, `s` 1 by default and negative to count down;
//!   `MakeList[g, a, b, s]` (`makelist`) applies `g`. The bounds and the
//!   step are rationals, `s ≠ 0`, with at most [`MAX_ELEMENTS`] elements.
//!   `Table` holds its arguments, as giac's `seq` does (CAS-DOC-33 §2,
//!   CAS-370): `f` is evaluated anew at each value, so that `seq(rand(), k,
//!   1, 3)` draws three numbers and `seq(int(x^k, x), k, -1, 1)` integrates
//!   `x^-1` as it is, and an index with a value is bound anew. A value that
//!   still holds `k`, as `f`'s where `f` has the value `k^2`, has it put in
//!   again, giac's binding of the index.
//! - **Reordering**: `Sort[l]` (`sort`) puts the real numbers first, by
//!   value (their balls, overlapping ones in their order), and the rest in
//!   the canonical order; `Sort[l, p]` sorts by `p(u, v)`, true when `u`
//!   comes first, stably. `Reverse` (`revlist`, `reverse`), and `Rotate[l,
//!   n]` (`rotate`) `n` places to the left, `n = −1` by default, as giac's.
//! - **Joining and slicing**: `Join` (`concat`), a non-list counting as one
//!   element; `Append` and `Prepend`; `First` (`head`) and `Rest` (`tail`),
//!   which stay on an empty list.
//! - **Reshaping**: `Flatten` (`flatten`), and `Dimensions` (`dim`), the
//!   length of a vector and `[rows, columns]` of a matrix.
//! - **Filtering**: `Select[p, l]`, `Reject[p, l]` and `Count[p, l]`
//!   (`select`, `remove`, `count`) by a predicate, a `Function` or a
//!   builtin, true or false for each element; another `p` means equality
//!   with it. `MemberIndex[x, l]` (`member`) is the first index of `x`,
//!   from 1, or 0.
//! - **Combining**: `Zip[f, l, m]` (`zip`) up to the shorter list, and
//!   `Accumulate` (`cumsum`), the partial sums.
//!
//! The argument orders are giac's: Mathematica's `Select[l, p]` takes the
//! list first. A predicate neither true nor false leaves the call.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of lists: f in k from a to b by s, a list l, elements x, u and v"
)]

use std::cmp::Ordering;

use cas_num::ball::Ball;

use super::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most elements `Table` and `MakeList` build.
const MAX_ELEMENTS: i64 = 1 << 16;

/// The bits of the balls that order numbers in `Sort`.
const BITS: u64 = 128;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `a, a + s, …` within `b`, rationals; `None` for another range, or one
/// of more than [`MAX_ELEMENTS`] elements.
fn values(a: &Expr, b: &Expr, s: &Expr) -> Option<Vec<Expr>> {
    let (a, b, s) = (a.to_rational()?, b.to_rational()?, s.to_rational()?);
    if s.is_zero() {
        return None;
    }
    let steps = (&(&b - &a) / &s).floor();
    if steps.is_negative() {
        return Some(Vec::new());
    }
    let count = steps.to_i64()?.checked_add(1)?;
    if count > MAX_ELEMENTS {
        return None;
    }
    let mut out = Vec::with_capacity(usize::try_from(count).ok()?);
    let mut v = a;
    for _ in 0..count {
        out.push(Expr::from(v.clone()));
        v = &v + &s;
    }
    Some(out)
}

/// `k = a .. b` as its parts: the index and the ends.
fn index_range(spec: &Expr) -> Option<(&Expr, &Expr, &Expr)> {
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [k, range]) =
        (spec.head().and_then(Expr::as_symbol), spec.args())
    else {
        return None;
    };
    let (Some(Symbol::RANGE), [a, b]) = (range.head().and_then(Expr::as_symbol), range.args())
    else {
        return None;
    };
    Some((k, a, b))
}

/// The summand, the index and its values of `Table`'s forms, held: the
/// bounds and the step evaluated here, and the range of `Table[f, r]`
/// when `r` is no `k = a .. b` as written.
fn table_parts(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<(Expr, Expr, Vec<Expr>)>, EvalError> {
    let one = Expr::from(1);
    let (f, k, a, b, s) = match e.args() {
        [f, k, a, b] => (f.clone(), k.clone(), a.clone(), b.clone(), one),
        [f, k, a, b, s] => (f.clone(), k.clone(), a.clone(), b.clone(), s.clone()),
        [f, spec] | [f, spec, _] => {
            let evaluated;
            let spec = if index_range(spec).is_some() {
                spec
            } else {
                evaluated = ev.eval(spec)?;
                &evaluated
            };
            let Some((k, a, b)) = index_range(spec) else {
                return Ok(None);
            };
            let s = e.args().get(2).cloned().unwrap_or(one);
            (f.clone(), k.clone(), a.clone(), b.clone(), s)
        }
        _ => return Ok(None),
    };
    if k.as_symbol().is_none() {
        return Ok(None);
    }
    let (a, b, s) = (ev.eval(&a)?, ev.eval(&b)?, ev.eval(&s)?);
    Ok(values(&a, &b, &s).map(|values| (f, k, values)))
}

pub(super) fn table(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, k, values)) = table_parts(ev, e)? else {
        return Ok(None);
    };
    let mut items = Vec::with_capacity(values.len());
    for v in values {
        let at = [(k.clone(), v)];
        let mut item = ev.eval(&f.substitute(&at))?;
        if item.contains(&k) {
            item = ev.eval(&item.substitute(&at))?;
        }
        items.push(item);
    }
    Ok(Some(list(items)))
}

pub(super) fn make_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let one = Expr::from(1);
    let (g, a, b, s) = match e.args() {
        [g, a, b] => (g, a, b, &one),
        [g, a, b, s] => (g, a, b, s),
        _ => return Ok(None),
    };
    let Some(values) = values(a, b, s) else {
        return Ok(None);
    };
    let mut items = Vec::with_capacity(values.len());
    for v in values {
        items.push(ev.eval(&Expr::apply(g, &[v]))?);
    }
    Ok(Some(list(items)))
}

/// Whether `u` must come before `v`; `None` when undecided.
type Before<'a, T> = dyn FnMut(&T, &T) -> Result<Option<bool>, EvalError> + 'a;

/// A stable merge sort of `items` by `before(u, v)`; `None` when a
/// comparison is undecided. Unlike the standard sorts, it asks for no total
/// order: balls that overlap neither come first, and stay in their order.
fn merge_sort<T: Clone>(
    items: Vec<T>,
    before: &mut Before<'_, T>,
) -> Result<Option<Vec<T>>, EvalError> {
    if items.len() <= 1 {
        return Ok(Some(items));
    }
    let mut left = items;
    let right = left.split_off(left.len() / 2);
    let (Some(left), Some(right)) = (merge_sort(left, before)?, merge_sort(right, before)?) else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() && j < right.len() {
        // The right one first only when it comes strictly before.
        match before(&right[j], &left[i])? {
            Some(true) => {
                out.push(right[j].clone());
                j += 1;
            }
            Some(false) => {
                out.push(left[i].clone());
                i += 1;
            }
            None => return Ok(None),
        }
    }
    out.extend_from_slice(&left[i..]);
    out.extend_from_slice(&right[j..]);
    Ok(Some(out))
}

/// The ball of a real number, for `Sort`.
fn real_ball(ev: &Eval<'_>, e: &Expr) -> Option<Ball> {
    if !e.is_numeric() {
        return None;
    }
    let value = approximation(ev.context(), e, BITS).ok()?;
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

pub(super) fn sort(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, p) = match e.args() {
        [l] => (l, None),
        [l, p] => (l, Some(p)),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return Ok(None);
    }
    if let Some(p) = p {
        let mut before = |u: &Expr, v: &Expr| -> Result<Option<bool>, EvalError> {
            Ok(truth(&ev.eval(&Expr::apply(p, &[u.clone(), v.clone()]))?))
        };
        return Ok(merge_sort(l.args().to_vec(), &mut before)?.map(list));
    }
    let mut numbers: Vec<(Ball, Expr)> = Vec::new();
    let mut rest: Vec<Expr> = Vec::new();
    for x in l.args() {
        match real_ball(ev, x) {
            Some(b) => numbers.push((b, x.clone())),
            None => rest.push(x.clone()),
        }
    }
    // Strictly below comes first; an overlap keeps the order.
    let mut below =
        |(a, _): &(Ball, Expr), (b, _): &(Ball, Expr)| -> Result<Option<bool>, EvalError> {
            Ok(Some(a.sub(b, BITS).sign() == Some(Ordering::Less)))
        };
    let Some(numbers) = merge_sort(numbers, &mut below)? else {
        return Ok(None);
    };
    rest.sort();
    let mut all: Vec<Expr> = numbers.into_iter().map(|(_, x)| x).collect();
    all.extend(rest);
    Ok(Some(list(all)))
}

pub(super) fn reverse(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    Ok(Some(list(l.args().iter().rev().cloned().collect())))
}

pub(super) fn rotate(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, n) = match e.args() {
        [l] => (l, -1),
        [l, n] => match n.as_i64() {
            Some(n) => (l, n),
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    if !is_list(l) {
        return Ok(None);
    }
    let items = l.args();
    if items.is_empty() {
        return Ok(Some(l.clone()));
    }
    let len = i64::try_from(items.len()).expect("a length in memory");
    let start = usize::try_from(n.rem_euclid(len)).expect("in range");
    let mut out = items[start..].to_vec();
    out.extend_from_slice(&items[..start]);
    Ok(Some(list(out)))
}

pub(super) fn join(_: &mut Eval<'_>, e: &Expr) -> Native {
    if e.args().is_empty() || !e.args().iter().any(is_list) {
        return Ok(None);
    }
    let mut out = Vec::new();
    for x in e.args() {
        if is_list(x) {
            out.extend_from_slice(x.args());
        } else {
            out.push(x.clone());
        }
    }
    Ok(Some(list(out)))
}

pub(super) fn append(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, x] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let mut out = l.args().to_vec();
    out.push(x.clone());
    Ok(Some(list(out)))
}

pub(super) fn prepend(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, x] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let mut out = vec![x.clone()];
    out.extend_from_slice(l.args());
    Ok(Some(list(out)))
}

pub(super) fn first(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) => Ok(l.args().first().cloned()),
        _ => Ok(None),
    }
}

pub(super) fn rest(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) && !l.args().is_empty() => Ok(Some(list(l.args()[1..].to_vec()))),
        _ => Ok(None),
    }
}

fn flattened(x: &Expr, out: &mut Vec<Expr>) {
    if is_list(x) {
        for y in x.args() {
            flattened(y, out);
        }
    } else {
        out.push(x.clone());
    }
}

pub(super) fn flatten(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let mut out = Vec::new();
    flattened(l, &mut out);
    Ok(Some(list(out)))
}

pub(super) fn dimensions(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let rows = l.args();
    let count = |n: usize| Expr::from(i64::try_from(n).expect("a length in memory"));
    let columns = rows.first().filter(|r| is_list(r)).map(|r| r.args().len());
    Ok(Some(match columns {
        Some(m) if rows.iter().all(|r| is_list(r) && r.args().len() == m) => {
            list(vec![count(rows.len()), count(m)])
        }
        _ => count(rows.len()),
    }))
}

/// `True`, `False`, or neither.
fn truth(v: &Expr) -> Option<bool> {
    match v.as_symbol() {
        Some(Symbol::TRUE) => Some(true),
        Some(Symbol::FALSE) => Some(false),
        _ => None,
    }
}

/// Whether `p` holds for `x`: `p[x]` for a function or a builtin, equality
/// with `x` otherwise; `None` when undecided.
fn holds(ev: &mut Eval<'_>, p: &Expr, x: &Expr) -> Result<Option<bool>, EvalError> {
    let function = p.head() == Some(&Expr::from(Symbol::FUNCTION))
        || p.as_symbol().is_some_and(Symbol::is_builtin);
    if function {
        Ok(truth(&ev.eval(&Expr::apply(p, std::slice::from_ref(x)))?))
    } else {
        Ok(Some(p == x))
    }
}

/// The elements of `l` where `p` holds (`keep`) or not; `None` when one is
/// undecided.
fn filtered(
    ev: &mut Eval<'_>,
    p: &Expr,
    l: &Expr,
    keep: bool,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut out = Vec::new();
    for x in l.args() {
        match holds(ev, p, x)? {
            Some(h) if h == keep => out.push(x.clone()),
            Some(_) => {}
            None => return Ok(None),
        }
    }
    Ok(Some(out))
}

fn predicate_and_list(e: &Expr) -> Option<(&Expr, &Expr)> {
    match e.args() {
        [p, l] if is_list(l) => Some((p, l)),
        _ => None,
    }
}

pub(super) fn select(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, l)) = predicate_and_list(e) else {
        return Ok(None);
    };
    Ok(filtered(ev, p, l, true)?.map(list))
}

pub(super) fn reject(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, l)) = predicate_and_list(e) else {
        return Ok(None);
    };
    Ok(filtered(ev, p, l, false)?.map(list))
}

pub(super) fn count(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, l)) = predicate_and_list(e) else {
        return Ok(None);
    };
    Ok(filtered(ev, p, l, true)?
        .map(|kept| Expr::from(i64::try_from(kept.len()).expect("a length in memory"))))
}

pub(super) fn member_index(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let index = l.args().iter().position(|y| y == x).map_or(0, |i| i + 1);
    Ok(Some(Expr::from(
        i64::try_from(index).expect("a length in memory"),
    )))
}

pub(super) fn zip(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, l, m] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) || !is_list(m) {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(l.args().len().min(m.args().len()));
    for (u, v) in l.args().iter().zip(m.args()) {
        out.push(ev.eval(&Expr::apply(f, &[u.clone(), v.clone()]))?);
    }
    Ok(Some(list(out)))
}

pub(super) fn accumulate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let mut total = Expr::from(0);
    let mut out = Vec::with_capacity(l.args().len());
    for x in l.args() {
        total = ev.eval(&(total + x.clone()))?;
        out.push(total.clone());
    }
    Ok(Some(list(out)))
}

/// `Sum[l]` and `Product[l]` of one list: `Plus` and `Times` of its
/// elements, giac's `sum([1,2,3,a])` (CAS-336); `None` for another form.
pub(super) fn of_list(e: &Expr, head: Symbol) -> Option<Expr> {
    match e.args() {
        [l] if is_list(l) => Some(Expr::call(head, l.args().to_vec())),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
