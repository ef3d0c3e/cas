//! Initial and boundary conditions (CAS-DOC-16 §8, CAS-264):
//! `DSolve[{eq, c₁, …, cₖ}, y[x]]`, each condition `y[x₀] = a` or
//! `Derivative[k][y][x₀] = b` (which `subst(diff(y(x), x), x = x₀)` gives),
//! `x₀` free of `x`, at one point or several.
//!
//! The general solution and its derivatives at the conditions' points make
//! equations in the constants:
//!
//! - a linear equation's are a linear system (`LinSolve`): its one solution
//!   gives the expression, free constants staying where the conditions are
//!   fewer than the order; no solution gives `{}`, as giac's `[]`;
//! - each branch of a nonlinear one has its constant from `Solve`, one
//!   branch per root, and the branches without one are dropped.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the conditions: y(x₀) = a, l = r"
)]

use super::super::solve::linear_system;
use super::derivative_head;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// `general` with its constants fixed by `conditions`; `None` when a
/// condition is not one, or its equations are not solved here.
pub(super) fn apply(
    ev: &mut Eval<'_>,
    e: &Expr,
    general: &Expr,
    conditions: &[Expr],
    y: &Expr,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let list = Expr::from(Symbol::LIST);
    let used: Vec<Expr> = constants
        .iter()
        .filter(|c| general.contains(c))
        .cloned()
        .collect();
    if general.head() != Some(&list) {
        let Some(equations) = equations(ev, general, conditions, y, x)? else {
            return Ok(None);
        };
        let Some(found) = linear_system(ev, e, &equations, &used)? else {
            return Ok(None);
        };
        return match found.args() {
            [] => Ok(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new()))),
            [values] if values.args().len() == used.len() => {
                let rules: Vec<(Expr, Expr)> = used
                    .iter()
                    .cloned()
                    .zip(values.args().iter().cloned())
                    .collect();
                ev.eval(&general.substitute(&rules)).map(Some)
            }
            _ => Ok(None),
        };
    }
    // The branches of a first-order equation, each with its constant.
    let [c] = used.as_slice() else {
        return Ok(None);
    };
    let mut kept = Vec::new();
    for branch in general.args() {
        let Some(equations) = equations(ev, branch, conditions, y, x)? else {
            return Ok(None);
        };
        let [equation] = equations.as_slice() else {
            return Ok(None);
        };
        let roots = ev.eval(&Expr::call(Symbol::SOLVE, [equation.clone(), c.clone()]))?;
        if roots.head() != Some(&list) || roots.contains(&Expr::from(Symbol::SOLVE)) {
            return Ok(None);
        }
        for root in roots.args() {
            kept.push(ev.eval(&branch.substitute(&[(c.clone(), root.clone())]))?);
        }
    }
    Ok(Some(Expr::call(Symbol::LIST, kept)))
}

/// Each condition `l = r` as `l − r = 0`, `y` and its derivatives at their
/// points taken from `solution`; `None` when one is not a condition at a
/// point.
fn equations(
    ev: &mut Eval<'_>,
    solution: &Expr,
    conditions: &[Expr],
    y: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let order = conditions
        .iter()
        .flat_map(|c| c.distinct())
        .filter_map(|part| match part.view() {
            View::Apply { head, .. } => derivative_head(head, y),
            _ => None,
        })
        .max()
        .unwrap_or(0);
    // The solution's derivatives, 0 to the highest a condition asks for.
    let mut derivatives = vec![solution.clone()];
    for k in 1..=order {
        let previous = derivatives[k - 1].clone();
        derivatives.push(ev.eval(&Expr::call(Symbol::D, [previous, x.clone()]))?);
    }
    let mut out = Vec::with_capacity(conditions.len());
    for condition in conditions {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) =
            (condition.head().and_then(Expr::as_symbol), condition.args())
        else {
            return Ok(None);
        };
        let at = (l.clone() - r.clone()).replace(|part| {
            let View::Apply { head, args: [a] } = part.view() else {
                return None;
            };
            if a.contains(x) {
                return None;
            }
            let k = if head == y {
                0
            } else {
                derivative_head(head, y)?
            };
            Some(derivatives[k].substitute(&[(x.clone(), a.clone())]))
        });
        if at.contains(y) || at.contains(x) {
            return Ok(None);
        }
        let value = ev.eval(&at)?;
        out.push(Expr::call(Symbol::EQUATION, [value, Expr::from(0)]));
    }
    Ok(Some(out))
}
