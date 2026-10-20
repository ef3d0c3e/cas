//! `Length`, `Part`, `Map` and `Apply` (CAS-DOC-4 §6.1).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;

type Native = Result<Option<Expr>, EvalError>;

/// `Length[e]`: the number of arguments; 0 for an atom.
pub(super) fn length(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(Expr::from(
        i64::try_from(x.args().len()).expect("fits"),
    )))
}

/// `Part[e, i, j, …]`, in the `Part` convention (CAS-25): 0 is the head, 1
/// the first argument, and a negative index counts from the end, -1 the
/// last. An index out of range leaves it unevaluated.
pub(super) fn part(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [whole, indexes @ ..] = e.args() else {
        return Ok(None);
    };
    let mut current = whole;
    for index in indexes {
        let Some(index) = index.as_i64() else {
            return Ok(None);
        };
        let len = i64::try_from(current.args().len()).expect("fits");
        // A negative index counts arguments only: `-len` is the first.
        let (position, first) = if index < 0 {
            (len + 1 + index, 1)
        } else {
            (index, 0)
        };
        if !(first..=len).contains(&position) || (position == 0 && current.is_atom()) {
            return Ok(None);
        }
        let position = usize::try_from(position).expect("in range");
        let Some(p) = current.part(&[position]) else {
            return Ok(None);
        };
        current = p;
    }
    Ok(Some(current.clone()))
}

/// `Map[f, e]`: `f` applied to each argument of `e`; an atom stays.
pub(super) fn map(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, x] = e.args() else {
        return Ok(None);
    };
    if x.is_atom() {
        return Ok(Some(x.clone()));
    }
    Ok(Some(
        x.map_args(|a| Expr::apply(f, std::slice::from_ref(a))),
    ))
}

/// `Apply[f, e]`: `e` with `f` for its head; an atom stays.
pub(super) fn apply(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, x] = e.args() else {
        return Ok(None);
    };
    if x.is_atom() {
        return Ok(Some(x.clone()));
    }
    Ok(Some(Expr::apply(f, x.args())))
}
