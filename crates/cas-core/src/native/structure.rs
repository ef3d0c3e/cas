//! `Length`, `Part`, `Map` and `Apply` (CAS-DOC-4 §6.1); `FiniteSet` and
//! `Unapply` (CAS-DOC-7 §4, §9).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

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

/// `FiniteSet[…]`, the language's `{…}`: its elements without repeats.
/// `Orderless` sorted them, so a repeat comes right after its element.
pub(super) fn finite_set(_: &mut Eval<'_>, e: &Expr) -> Native {
    let elements = e.args();
    if elements.windows(2).all(|w| w[0] != w[1]) {
        return Ok(None);
    }
    let mut unique = elements.to_vec();
    unique.dedup();
    Ok(Some(Expr::apply(&Symbol::FINITE_SET.into(), &unique)))
}

/// `Unapply[body, x…]`: the function `Function[{x…}, body]` of at least one
/// variable, user symbols without repeats. Its body is the evaluated
/// argument, so it is computed once, when the function is made:
/// `unapply(diff(x^3, x), x)` is `x -> 3 x^2`, where `x -> diff(x^3, x)`
/// would substitute first, and at 2 differentiate 8 by 2. Other variables
/// leave it unevaluated.
pub(super) fn unapply(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [body, variables @ ..] = e.args() else {
        return Ok(None);
    };
    let user = |v: &Expr| {
        v.as_symbol()
            .is_some_and(|s| s.builtin_attributes().is_none())
    };
    if variables.is_empty()
        || !variables.iter().all(user)
        || (1..variables.len()).any(|i| variables[..i].contains(&variables[i]))
    {
        return Ok(None);
    }
    let parameters = Expr::apply(&Symbol::LIST.into(), variables);
    Ok(Some(Expr::apply(
        &Symbol::FUNCTION.into(),
        &[parameters, body.clone()],
    )))
}
