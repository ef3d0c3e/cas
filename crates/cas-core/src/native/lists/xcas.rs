//! Xcas's list commands (CAS-555 to CAS-557), giac's answers: positions
//! from 0, counts, differences and the small builders. giac's errors leave
//! the call, with the error as a warning.
//!
//! - `FindPositions[x, l]` (`find`, CAS-555): the positions of `x` in `l`,
//!   from 0 (misc.cc `_find`); the arguments exchanged when only the first
//!   is a list, as giac's. `FirstPosition[l, x]` (`index`, prog.cc
//!   `_index`) is the first of them, and leaves the call when there is none
//!   (giac's Bad Argument Value). Not `MemberIndex`, which counts from 1.
//!   An element is `x` when it is the same expression, or when both are
//!   numbers that `Equal` makes equal: giac finds `1` at `1.0`. The string
//!   forms stay (phase 1.6), as does `find`'s start position.
//! - `CountEqual[x, l]`, `CountLess[x, l]`, `CountGreater[x, l]`
//!   (`count_eq`, `count_inf`, `count_sup`, maple.cc `count_eq0` to
//!   `count_eq2`): the elements of `l` at every depth equal to, less than,
//!   greater than `x`, an `l` that is not a list counting itself; an order
//!   counts only where `Less` gives `True`, as giac's `is_strictly_greater`
//!   counts only what it decides (`count_inf(x,[1,2])` is 0). With `row`
//!   (`MatrixRow`, C0008) or `col` (`MatrixColumn`), the counts of each row
//!   or column of a matrix, in a list; another option is ignored.
//! - `Differences[l]` (`deltalist`, ti89.cc `_deltalist`): successive
//!   differences, row by row for a matrix (`[[1,2],[4,8]]` gives
//!   `[[3,6]]`); fewer than two elements is giac's Invalid dimension.
//! - `AccumulateHeadTail[l, p, q]` (`accumulate_head_tail`, maple.cc): the
//!   sum of the first `p` elements, the middle ones, the sum of the last
//!   `q`; `p, q ≥ 1` integers with `p + q ≤ length`.
//! - `Lengths[l]` (`sizes`, permu.cc `_sizes`): the lengths of the
//!   sublists.
//! - `MakeVector[a, …]` (`makevector`): the list of the arguments, one list
//!   itself; more than [`MAX_ELEMENTS`] arguments leave the call (C0019).
//!   `MakeSequence[l]` (`makesuite`): the sequence of `l`'s elements, a
//!   non-list `x` the list `[x]`. `NewList[n]` (`newList`): `|n|` zeros, at
//!   most [`MAX_ELEMENTS`] (giac's `LIST_SIZE_LIMIT`, larger). `Enumerate[l,
//!   s]` (`enumerate`, misc.cc): `[[s + k, l_k]]`, `s` 0 by default.
//!   `Copy[x]` (`copy`) is `x`.

use super::{MAX_ELEMENTS, is_list, list};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

const BAD_VALUE: &str = "bad argument value";
const BAD_TYPE: &str = "bad argument type";
const SIZE_LIMIT: &str = "calculation size limit exceeded";

/// giac's error for `e`: the call stays, with `why` as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

fn integer(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("a length in memory"))
}

fn is_string(e: &Expr) -> bool {
    e.as_str().is_some()
}

/// Whether `ev` decides `test` is `True`.
fn holds(ev: &mut Eval<'_>, test: Symbol, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    let value = ev.eval(&Expr::call(test, [a.clone(), b.clone()]))?;
    Ok(value.as_symbol() == Some(Symbol::TRUE))
}

/// giac's `==`: the same expression, or numbers `Equal` makes equal.
fn same(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    if a == b {
        return Ok(true);
    }
    if a.is_number() && b.is_number() {
        return holds(ev, Symbol::EQUAL, a, b);
    }
    Ok(false)
}

/// The positions of `x` in the list `l`, from 0.
fn positions(ev: &mut Eval<'_>, x: &Expr, l: &Expr) -> Result<Vec<usize>, EvalError> {
    let mut out = Vec::new();
    for (i, y) in l.args().iter().enumerate() {
        if same(ev, x, y)? {
            out.push(i);
        }
    }
    Ok(out)
}

/// `x` and the list to search, giac's exchange included; `None` when
/// neither is a list.
fn searched<'e>(x: &'e Expr, l: &'e Expr) -> Option<(&'e Expr, &'e Expr)> {
    if is_list(l) {
        Some((x, l))
    } else if is_list(x) {
        Some((l, x))
    } else {
        None
    }
}

/// `FindPositions[x, l]`.
pub(in crate::native) fn find_positions(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, l] = e.args() else {
        return Ok(None);
    };
    if is_string(x) || is_string(l) {
        return Ok(None);
    }
    let Some((x, l)) = searched(x, l) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let found = positions(ev, x, l)?;
    Ok(Some(list(found.into_iter().map(integer).collect())))
}

/// `FirstPosition[l, x]`.
pub(in crate::native) fn first_position(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, x] = e.args() else {
        return refuse(ev, e, BAD_VALUE);
    };
    if is_string(x) || is_string(l) {
        return Ok(None);
    }
    let Some((x, l)) = searched(x, l) else {
        return refuse(ev, e, BAD_VALUE);
    };
    match positions(ev, x, l)?.first() {
        Some(&i) => Ok(Some(integer(i))),
        None => refuse(ev, e, BAD_VALUE),
    }
}

/// What `count_eq`, `count_inf` and `count_sup` count.
#[derive(Clone, Copy)]
enum Relation {
    Equal,
    Less,
    Greater,
}

/// The elements of `v` at every depth in `relation` to `x`; `v` itself
/// when it is not a list.
fn counted(ev: &mut Eval<'_>, relation: Relation, x: &Expr, v: &Expr) -> Result<usize, EvalError> {
    if is_list(v) {
        let mut total = 0;
        for y in v.args() {
            total += counted(ev, relation, x, y)?;
        }
        return Ok(total);
    }
    let yes = match relation {
        Relation::Equal => same(ev, v, x)?,
        Relation::Less => holds(ev, Symbol::LESS, v, x)?,
        Relation::Greater => holds(ev, Symbol::LESS, x, v)?,
    };
    Ok(usize::from(yes))
}

/// giac's `ckmatrix`: a non-empty list of lists of one non-zero length.
fn is_matrix(m: &Expr) -> bool {
    let rows = m.args();
    is_list(m)
        && rows
            .first()
            .is_some_and(|r| is_list(r) && !r.args().is_empty())
        && rows
            .iter()
            .all(|r| is_list(r) && r.args().len() == rows[0].args().len())
}

fn count(ev: &mut Eval<'_>, e: &Expr, relation: Relation) -> Native {
    let (x, v, option) = match e.args() {
        [x, v] => (x, v, None),
        [x, v, o, ..] => (x, v, o.as_symbol()),
        _ => return refuse(ev, e, BAD_VALUE),
    };
    let rows: Vec<Expr> = match option {
        Some(Symbol::MATRIX_ROW | Symbol::MATRIX_COLUMN) if !is_matrix(v) => {
            return refuse(ev, e, BAD_TYPE);
        }
        Some(Symbol::MATRIX_ROW) => v.args().to_vec(),
        Some(Symbol::MATRIX_COLUMN) => (0..v.args()[0].args().len())
            .map(|j| list(v.args().iter().map(|r| r.args()[j].clone()).collect()))
            .collect(),
        _ => return Ok(Some(integer(counted(ev, relation, x, v)?))),
    };
    let mut out = Vec::with_capacity(rows.len());
    for r in &rows {
        out.push(integer(counted(ev, relation, x, r)?));
    }
    Ok(Some(list(out)))
}

/// `CountEqual[x, l]`.
pub(in crate::native) fn count_equal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Equal)
}

/// `CountLess[x, l]`.
pub(in crate::native) fn count_less(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Less)
}

/// `CountGreater[x, l]`.
pub(in crate::native) fn count_greater(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Greater)
}

/// `Differences[l]`.
pub(in crate::native) fn differences(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    if l.args().len() < 2 {
        return refuse(ev, e, "invalid dimension");
    }
    let mut out = Vec::with_capacity(l.args().len() - 1);
    for pair in l.args().windows(2) {
        out.push(ev.eval(&(pair[1].clone() - pair[0].clone()))?);
    }
    Ok(Some(list(out)))
}

/// `AccumulateHeadTail[l, p, q]`.
pub(in crate::native) fn accumulate_head_tail(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, p, q] = e.args() else {
        return Ok(None);
    };
    let length = |n: &Expr| n.to_integer().and_then(|n| n.to_i64());
    let (Some(p), Some(q)) = (length(p), length(q)) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let w = l.args();
    let s = i64::try_from(w.len()).expect("a length in memory");
    if !is_list(l) || p < 1 || q < 1 || p.saturating_add(q) > s {
        return refuse(ev, e, BAD_VALUE);
    }
    let (p, t) = (
        usize::try_from(p).expect("in range"),
        usize::try_from(s - q).expect("in range"),
    );
    let mut out = Vec::with_capacity(t - p + 2);
    out.push(ev.eval(&Expr::call(Symbol::PLUS, w[..p].to_vec()))?);
    out.extend_from_slice(&w[p..t]);
    out.push(ev.eval(&Expr::call(Symbol::PLUS, w[t..].to_vec()))?);
    Ok(Some(list(out)))
}

/// `Lengths[l]`.
pub(in crate::native) fn lengths(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_TYPE);
    }
    if !l.args().iter().all(is_list) {
        return refuse(ev, e, BAD_VALUE);
    }
    Ok(Some(list(
        l.args().iter().map(|r| integer(r.args().len())).collect(),
    )))
}

/// `MakeVector[a, …]`.
pub(in crate::native) fn make_vector(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) => Ok(Some(l.clone())),
        args if i64::try_from(args.len()).is_ok_and(|n| n <= MAX_ELEMENTS) => {
            Ok(Some(list(args.to_vec())))
        }
        _ => refuse(ev, e, SIZE_LIMIT),
    }
}

/// `MakeSequence[l]`.
pub(in crate::native) fn make_sequence(_: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(Some(match e.args() {
        [l] if is_list(l) => Expr::call(Symbol::SEQUENCE, l.args().to_vec()),
        [x] => list(vec![x.clone()]),
        args => Expr::call(Symbol::SEQUENCE, args.to_vec()),
    }))
}

/// `NewList[n]`.
pub(in crate::native) fn new_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return refuse(ev, e, BAD_VALUE);
    };
    let Some(n) = n.to_integer() else {
        return refuse(ev, e, BAD_VALUE);
    };
    match n.to_i64().and_then(i64::checked_abs) {
        Some(n) if n <= MAX_ELEMENTS => Ok(Some(list(vec![
            Expr::from(0);
            usize::try_from(n).expect("in range")
        ]))),
        _ => refuse(ev, e, SIZE_LIMIT),
    }
}

/// `Enumerate[l, s]`.
pub(in crate::native) fn enumerate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, start) = match e.args() {
        [] => return Ok(Some(list(Vec::new()))),
        [l] => (l, Expr::from(0)),
        [l, s] => (l, s.clone()),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_TYPE);
    }
    let mut out = Vec::with_capacity(l.args().len());
    for (k, x) in l.args().iter().enumerate() {
        let index = ev.eval(&(start.clone() + integer(k)))?;
        out.push(list(vec![index, x.clone()]));
    }
    Ok(Some(list(out)))
}

/// `Copy[x]`.
pub(in crate::native) fn copy(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(x.clone()))
}
