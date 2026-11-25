//! Xcas's sorting commands (CAS-622): `SortA`/`sorta`, `SortD`/`sortd`
//! (giac ti89.cc `_SortA` to `_sortd`) and `sortperm` (optimization.cc
//! `_sortperm`), in `Sort`'s order (`lists::order`: numbers by value, then
//! the canonical order, strings in byte order; stable).
//!
//! - `SortAscending[l]` (`SortA`, `sorta`): `l` in `Sort`'s order; a
//!   sequence of values is a list (`SortA(3,4,2)` is `[2,3,4]`). A matrix
//!   (rows of equal length whose entries are not lists), or a sequence of
//!   such rows, has its columns reordered by its first row:
//!   `SortA([[3,4,2],[6,4,5]])` is `[[2,3,4],[5,6,4]]`. Anything else is
//!   sorted as a list, as giac's (`SortA([3,1],x)` sorts two elements).
//!   `SortDescending[l]` (`SortD`, `sortd`) is the same order reversed, as
//!   giac's: ties come out reversed too (`SortD([1.5,3/2])` is
//!   `[3/2,1.5]`). A lone non-list stays (giac's Bad Argument Value).
//!   giac's `SortA(v)` with `v` a variable sorts `v`'s value and stores it
//!   back into `v` (its arguments quoted); not ported: ours sorts the value
//!   and stores nothing.
//! - `SortPermutation[l]` (`sortperm`): the positions, from 0, that put `l`
//!   in order: `l[p[0]] ≤ l[p[1]] ≤ …`, stable (`sortperm([3,1,2,1])` is
//!   `[1,3,2,0]`). giac answers a one-element list `[array_start]` without
//!   comparing, `[0]` in Xcas's mode (1 in its Maple and TI modes); ours is
//!   `[0]`, the same as the order. `SortPermutation[l, p]`: `[l[p[0]],
//!   l[p[1]], …]` when `p` is a permutation of `0 … n − 1` for `l` of
//!   length `n`; another `p` stays (giac's Bad Argument Type).
//!
//! The answers are as long as the input, so no limit applies beyond the
//! input's own.

use super::{is_list, list, order};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// The rows of a matrix: at least one, lists of equal length, no entry a
/// list.
fn rows(items: &[Expr]) -> Option<Vec<&[Expr]>> {
    let first = items.first()?;
    if !is_list(first) {
        return None;
    }
    let n = first.args().len();
    items
        .iter()
        .map(|r| {
            (is_list(r) && r.args().len() == n && !r.args().iter().any(is_list)).then(|| r.args())
        })
        .collect()
}

/// `SortA` and `SortD`.
fn sorted(ev: &mut Eval<'_>, e: &Expr, ascending: bool) -> Native {
    let items = match e.args() {
        [l] if is_list(l) => l.args(),
        [] | [_] => return Ok(None),
        args => args,
    };
    let rows = rows(items);
    let key = match &rows {
        Some(rows) => rows[0],
        None => items,
    };
    let Some(mut order) = order(ev, key)? else {
        return Ok(None);
    };
    if !ascending {
        order.reverse();
    }
    let permuted = |r: &[Expr]| list(order.iter().map(|&i| r[i].clone()).collect());
    Ok(Some(match rows {
        Some(rows) => list(rows.into_iter().map(permuted).collect()),
        None => permuted(items),
    }))
}

pub(in crate::native) fn sort_ascending(ev: &mut Eval<'_>, e: &Expr) -> Native {
    sorted(ev, e, true)
}

pub(in crate::native) fn sort_descending(ev: &mut Eval<'_>, e: &Expr) -> Native {
    sorted(ev, e, false)
}

/// `p` as positions when it is a permutation of `0 … n − 1`.
fn permutation(p: &Expr, n: usize) -> Option<Vec<usize>> {
    if !is_list(p) || p.args().len() != n {
        return None;
    }
    let mut seen = vec![false; n];
    p.args()
        .iter()
        .map(|i| {
            let i = usize::try_from(i.as_i64()?).ok()?;
            (i < n && !std::mem::replace(&mut seen[i], true)).then_some(i)
        })
        .collect()
}

pub(in crate::native) fn sort_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) => {
            let Some(order) = order(ev, l.args())? else {
                return Ok(None);
            };
            Ok(Some(list(order.into_iter().map(integer).collect())))
        }
        [l, p] if is_list(l) => Ok(permutation(p, l.args().len())
            .map(|p| list(p.into_iter().map(|i| l.args()[i].clone()).collect()))),
        _ => Ok(None),
    }
}

fn integer(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("a list's position fits"))
}
