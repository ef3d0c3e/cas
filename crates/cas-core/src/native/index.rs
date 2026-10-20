//! Indexing from 0 (CAS-DOC-7 §7): `At[v, k…]`, the language's `v[k…]`,
//! and `ReplaceAt[v, k…, x]`, a copy of `v` with what `v[k…]` selects
//! replaced by `x`, which `v[k…] := x` assigns.
//!
//! An index is an integer, from 0 for the first element, or negative,
//! counting from the end, -1 for the last; or a slice `Range[a, b]` of such
//! integers, from `a` to `b` included, as in Xcas. A slice whose end comes
//! before its start is empty, so that `v[1..-1]`, all but the first, is
//! empty for one element. Otherwise its elements are in the list.
//!
//! Each index goes down a level: `m[i, j]` is an entry of a matrix, `m[i]`
//! a row and `m[0..1, j]` part of a column. An integer takes its level away
//! and a slice keeps it, as a list. `m[i, j]` is `m[i][j]`.
//!
//! Only lists are indexed. Anything else, as `x[0]` for an unassigned `x`
//! (an indexed name), and an index that is neither an integer nor a slice,
//! as `v[k]`, stay unevaluated. An index out of range is
//! [`EvalError::IndexOutOfRange`]: unlike Xcas, a slice past the end is not
//! cut short and an assignment past the end does not extend the list.

use std::ops::Range;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// What an index selects in a list.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Index {
    /// One element.
    One(usize),
    /// The elements of a slice, in a list.
    Slice(Range<usize>),
}

/// The value of an integer index: one beyond `i64` is out of range of any
/// list, so its sign is enough.
fn integer(e: &Expr) -> Option<i128> {
    let n = e.to_integer()?;
    Some(n.to_i64().map_or(
        if n.is_negative() {
            i128::MIN
        } else {
            i128::MAX
        },
        i128::from,
    ))
}

impl Index {
    /// What `index` selects in a list of `len` elements; `None` when it is
    /// neither an integer nor a slice of integers.
    fn of(index: &Expr, len: usize) -> Result<Option<Index>, EvalError> {
        let len = i128::try_from(len).expect("fits");
        let from_end = |k: i128| if k < 0 { k + len } else { k };
        if let Some(k) = integer(index) {
            let k = from_end(k);
            if !(0..len).contains(&k) {
                return Err(EvalError::IndexOutOfRange);
            }
            return Ok(Some(Index::One(usize::try_from(k).expect("in range"))));
        }
        if index.head().and_then(Expr::as_symbol) != Some(Symbol::RANGE) {
            return Ok(None);
        }
        let [a, b] = index.args() else {
            return Ok(None);
        };
        let (Some(a), Some(b)) = (integer(a), integer(b)) else {
            return Ok(None);
        };
        let (a, b) = (from_end(a), from_end(b));
        // An empty slice may start just past the end, or end just before
        // the start.
        if !(0..=len).contains(&a) || !(-1..len).contains(&b) {
            return Err(EvalError::IndexOutOfRange);
        }
        let start = usize::try_from(a).expect("in range");
        let end = usize::try_from(b + 1).expect("in range").max(start);
        Ok(Some(Index::Slice(start..end)))
    }
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

fn list(items: &[Expr]) -> Expr {
    Expr::apply(&Symbol::LIST.into(), items)
}

/// `At[v, k…]`: what the indices `k…` select in the list `v`.
pub(super) fn at(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v, indices @ ..] = e.args() else {
        return Ok(None);
    };
    if indices.is_empty() {
        return Ok(None);
    }
    select(ev, v, indices)
}

/// What `indices`, at least one, select in `v`, if it is a list.
fn select(ev: &mut Eval<'_>, v: &Expr, indices: &[Expr]) -> Native {
    let Some((index, rest)) = indices.split_first() else {
        return Ok(Some(v.clone()));
    };
    if !is_list(v) {
        return Ok(None);
    }
    let elements = v.args();
    let range = match Index::of(index, elements.len())? {
        None => return Ok(None),
        Some(Index::One(k)) => return below(ev, &elements[k], rest).map(Some),
        Some(Index::Slice(range)) => range,
    };
    let mut items = Vec::with_capacity(range.len());
    for x in &elements[range] {
        items.push(below(ev, x, rest)?);
    }
    Ok(Some(list(&items)))
}

/// What `rest` selects in `x`, a level down, which takes a level of
/// `Limits::recursion`; when `x` is not a list, `At` of what is left,
/// unevaluated, so `v[i, j]` is `v[i][j]`.
fn below(ev: &mut Eval<'_>, x: &Expr, rest: &[Expr]) -> Result<Expr, EvalError> {
    if rest.is_empty() {
        return Ok(x.clone());
    }
    ev.enter()?;
    let selected = select(ev, x, rest);
    ev.leave();
    Ok(selected?.unwrap_or_else(|| {
        let mut args = vec![x.clone()];
        args.extend_from_slice(rest);
        Expr::apply(&Symbol::AT.into(), &args)
    }))
}

/// `ReplaceAt[v, k…, x]`: `v` with what the indices `k…` select replaced
/// by `x`; a slice by a list of its length, element by element. A value of
/// another shape for a slice is [`EvalError::Shape`]: a scalar is not
/// spread over it (CAS-DOC-4 O8), unlike in Xcas.
pub(super) fn replace_at(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v, indices @ .., x] = e.args() else {
        return Ok(None);
    };
    if indices.is_empty() {
        return Ok(None);
    }
    replace(ev, v, indices, x)
}

/// `v` with what `indices` select replaced by `x`; `None` past a level that
/// is not a list.
fn replace(ev: &mut Eval<'_>, v: &Expr, indices: &[Expr], x: &Expr) -> Native {
    let Some((index, rest)) = indices.split_first() else {
        return Ok(Some(x.clone()));
    };
    if !is_list(v) {
        return Ok(None);
    }
    let mut elements = v.args().to_vec();
    match Index::of(index, elements.len())? {
        None => return Ok(None),
        Some(Index::One(k)) => match deeper(ev, &elements[k], rest, x)? {
            Some(new) => elements[k] = new,
            None => return Ok(None),
        },
        Some(Index::Slice(range)) => {
            if !is_list(x) || x.args().len() != range.len() {
                return Err(EvalError::Shape);
            }
            for (k, new) in range.zip(x.args()) {
                match deeper(ev, &elements[k], rest, new)? {
                    Some(new) => elements[k] = new,
                    None => return Ok(None),
                }
            }
        }
    }
    Ok(Some(list(&elements)))
}

/// [`replace`] a level down, which takes a level of `Limits::recursion`.
fn deeper(ev: &mut Eval<'_>, v: &Expr, indices: &[Expr], x: &Expr) -> Native {
    if indices.is_empty() {
        return Ok(Some(x.clone()));
    }
    ev.enter()?;
    let replaced = replace(ev, v, indices, x);
    ev.leave();
    replaced
}

#[cfg(test)]
mod tests;
