//! Linear algebra (CAS-DOC-12): matrices as lists of lists, a vector a list
//! whose entries are not lists ([`linear::shape`]), and the products of
//! CAS-DOC-4 §6.2 in [`linear`].
//!
//! - `Transpose[m]`, Xcas's `tran`: the transpose, a vector as a column
//!   matrix.
//! - `IdentityMatrix[n]`, Xcas's `idn`: the `n × n` identity, for
//!   `1 ≤ n` with `n²` within `Limits::terms`.
//! - `Tr[m]`, Xcas's `trace`: the sum of the diagonal of a square matrix;
//!   another stays with a warning, as giac's stays.
//! - `Dot[u, v]`, Xcas's `dot`: `Σ conj(uₖ) vₖ` of two vectors of one
//!   length, and entrywise of two matrices of one shape, the conjugate of
//!   the first as giac's (`dot([i,0],[1,0])` is `-i`). An entry is its own
//!   conjugate when it is known real; another leaves the call.
//! - `Cross[u, v]`, Xcas's `cross`: the cross product of 3-vectors, a
//!   2-vector taken with a third entry 0 (`cross([1,2],[3,4])` is
//!   `[0,0,-2]`).
//!
//! Dimensions that do not fit are [`EvalError::Shape`], as in products.
//! Results are built unevaluated, for the loop to evaluate; each call takes
//! one step, named after the giac function.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear algebra: vectors u and v, an m × n matrix"
)]

use super::linear::{Shape, is_list, shape};
use crate::assume::{Set, Truth, element};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn entry(m: &Expr, i: usize, j: usize) -> &Expr {
    &m.args()[i].args()[j]
}

/// A step and its result.
fn stepped(ev: &mut Eval<'_>, name: &'static str, e: &Expr, result: Expr) -> Native {
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// Says that a square matrix is needed.
fn not_square(ev: &mut Eval<'_>, e: &Expr) {
    ev.warn(Warning::Native {
        message: "a square matrix is needed".into(),
        expr: e.clone(),
    });
}

/// `Transpose[m]`.
pub(super) fn transpose(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let result = match shape(m) {
        Some(Shape::Matrix(rows, columns)) => {
            list((0..columns).map(|j| list((0..rows).map(|i| entry(m, i, j).clone()))))
        }
        Some(Shape::Vector(n)) if n > 0 => list(m.args().iter().map(|x| list([x.clone()]))),
        _ => return Ok(None),
    };
    stepped(ev, "tran", e, result)
}

/// `IdentityMatrix[n]`.
pub(super) fn identity_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.as_i64().and_then(|n| usize::try_from(n).ok()) else {
        return Ok(None);
    };
    if n == 0 || n.saturating_mul(n) > ev.context().limits().terms {
        return Ok(None);
    }
    let result = list((0..n).map(|i| list((0..n).map(|j| Expr::from(i64::from(i == j))))));
    stepped(ev, "idn", e, result)
}

/// `Tr[m]`.
pub(super) fn tr(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    if rows != columns {
        not_square(ev, e);
        return Ok(None);
    }
    let result = Expr::call(
        Symbol::PLUS,
        (0..rows)
            .map(|i| entry(m, i, i).clone())
            .collect::<Vec<_>>(),
    );
    stepped(ev, "trace", e, result)
}

/// The conjugate of the entry `x`: of a number, exactly or of its ball; `x`
/// itself when it is known real; `None` otherwise.
fn conjugate(ctx: &Context, x: &Expr) -> Option<Expr> {
    match x.view() {
        View::ComplexBall(z) => return Some(Expr::from(z.conj())),
        View::Ball(_) => return Some(x.clone()),
        _ => {}
    }
    if let Some(z) = x.to_complex() {
        return Some(Expr::from(&z.conj()));
    }
    (element(ctx, x, Set::Reals) == Truth::True).then(|| x.clone())
}

/// `Dot[u, v]`.
///
/// # Errors
///
/// [`EvalError::Shape`] for two vectors or matrices of other dimensions.
pub(super) fn dot(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [u, v] = e.args() else {
        return Ok(None);
    };
    let pairs: Vec<(&Expr, &Expr)> = match (shape(u), shape(v)) {
        (Some(Shape::Vector(n)), Some(Shape::Vector(k))) if n == k => {
            u.args().iter().zip(v.args()).collect()
        }
        (Some(Shape::Matrix(m, n)), Some(Shape::Matrix(k, p))) if (m, n) == (k, p) => u
            .args()
            .iter()
            .zip(v.args())
            .flat_map(|(a, b)| a.args().iter().zip(b.args()))
            .collect(),
        (Some(_), Some(_)) => return Err(EvalError::Shape),
        _ => return Ok(None),
    };
    let ctx = ev.context();
    let mut terms = Vec::with_capacity(pairs.len());
    for (a, b) in pairs {
        let Some(a) = conjugate(ctx, a) else {
            return Ok(None);
        };
        terms.push(Expr::call(Symbol::TIMES, [a, b.clone()]));
    }
    stepped(ev, "dot", e, Expr::call(Symbol::PLUS, terms))
}

/// `Cross[u, v]`.
///
/// # Errors
///
/// [`EvalError::Shape`] for vectors of other lengths.
pub(super) fn cross(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [u, v] = e.args() else {
        return Ok(None);
    };
    if !is_list(u) || !is_list(v) {
        return Ok(None);
    }
    let three = |w: &Expr| -> Option<[Expr; 3]> {
        match (shape(w)?, w.args()) {
            (Shape::Vector(3), [a, b, c]) => Some([a.clone(), b.clone(), c.clone()]),
            (Shape::Vector(2), [a, b]) => Some([a.clone(), b.clone(), Expr::from(0)]),
            _ => None,
        }
    };
    let (Some([u1, u2, u3]), Some([v1, v2, v3])) = (three(u), three(v)) else {
        return if shape(u).is_some() && shape(v).is_some() {
            Err(EvalError::Shape)
        } else {
            Ok(None)
        };
    };
    let minus = |a: Expr, b: Expr, c: Expr, d: Expr| a * b - c * d;
    let result = list([
        minus(u2.clone(), v3.clone(), u3.clone(), v2.clone()),
        minus(u3, v1.clone(), u1.clone(), v3),
        minus(u1, v2, u2, v1),
    ]);
    stepped(ev, "cross", e, result)
}

#[cfg(test)]
mod tests;
