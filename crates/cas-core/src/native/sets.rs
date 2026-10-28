//! Set operations (CAS-DOC-32 §2, CAS-365): `Union`, `Intersection` and
//! `Complement` of finite sets and of lists read as sets, giac's infix
//! `union`, `intersect` and `minus`; `IncludedQ`, giac's `is_included`; and
//! `Element` of a finite set or a list.
//!
//! - The operands are `FiniteSet`s and `List`s; anything else among them,
//!   as a symbol, leaves the call.
//! - The results are `FiniteSet`s, whose elements are the same when they
//!   are the same expression after evaluation, as `FiniteSet` removes
//!   repeats: `{1, 1.0}` keeps both.
//! - `IncludedQ[A, B]` is `True` when each element of `A` is in `B`, and
//!   `False` when one is decided not to be.
//! - `Element[x, A]` is `True` when an element of `A` is `x`, `False` when
//!   each is decided different from it (numbers, or by `Unequal`), and
//!   stays otherwise: `x in {y}` for symbols.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The elements of a finite set or a list.
pub(super) fn elements(e: &Expr) -> Option<&[Expr]> {
    let head = e.head()?.as_symbol()?;
    matches!(head, Symbol::FINITE_SET | Symbol::LIST).then(|| e.args())
}

/// The finite set of `items`, its repeats removed by evaluation.
fn set(ev: &mut Eval<'_>, items: Vec<Expr>) -> Native {
    ev.eval(&Expr::call(Symbol::FINITE_SET, items)).map(Some)
}

/// The operands' elements, or `None` when one is no set or list.
fn operands(e: &Expr) -> Option<Vec<&[Expr]>> {
    e.args().iter().map(elements).collect()
}

/// `Union[A, B, …]`.
pub(super) fn union(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(sets) = operands(e) else {
        return Ok(None);
    };
    set(ev, sets.concat())
}

/// `Intersection[A, B, …]`: the elements of the first that are in every
/// other.
pub(super) fn intersection(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(sets) = operands(e) else {
        return Ok(None);
    };
    let Some((first, rest)) = sets.split_first() else {
        return Ok(None);
    };
    let kept = first
        .iter()
        .filter(|x| rest.iter().all(|s| s.contains(x)))
        .cloned()
        .collect();
    set(ev, kept)
}

/// `Complement[A, B, …]`: the elements of the first in no other.
pub(super) fn complement(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(sets) = operands(e) else {
        return Ok(None);
    };
    let Some((first, rest)) = sets.split_first() else {
        return Ok(None);
    };
    let kept = first
        .iter()
        .filter(|x| !rest.iter().any(|s| s.contains(x)))
        .cloned()
        .collect();
    set(ev, kept)
}

/// Whether `x` is one of `items`: `Some(true)` when one is `x`,
/// `Some(false)` when each is decided different from it, `None` otherwise.
pub(super) fn member(
    ev: &mut Eval<'_>,
    x: &Expr,
    items: &[Expr],
) -> Result<Option<bool>, EvalError> {
    if items.contains(x) {
        return Ok(Some(true));
    }
    for item in items {
        if x.is_number() && item.is_number() {
            // Two numbers that are not the same: equal only as balls that
            // overlap, which `Unequal` tells.
            if x.to_rational().is_some() && item.to_rational().is_some() {
                continue;
            }
        }
        let unequal = ev.eval(&Expr::call(Symbol::UNEQUAL, [x.clone(), item.clone()]))?;
        if unequal.as_symbol() != Some(Symbol::TRUE) {
            return Ok(None);
        }
    }
    Ok(Some(false))
}

fn boolean(b: bool) -> Expr {
    Expr::from(if b { Symbol::TRUE } else { Symbol::FALSE })
}

/// `IncludedQ[A, B]`: whether `A ⊆ B`.
pub(super) fn included_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let (Some(a), Some(b)) = (elements(a), elements(b)) else {
        return Ok(None);
    };
    let mut known = true;
    for x in a {
        match member(ev, x, b)? {
            Some(true) => {}
            Some(false) => return Ok(Some(boolean(false))),
            None => known = false,
        }
    }
    Ok(known.then(|| boolean(true)))
}

/// `Element[x, A]` for a finite set or a list `A`; `None` for another.
pub(super) fn element(ev: &mut Eval<'_>, x: &Expr, s: &Expr) -> Native {
    let Some(items) = elements(s) else {
        return Ok(None);
    };
    Ok(member(ev, x, items)?.map(boolean))
}

#[cfg(test)]
mod tests;
