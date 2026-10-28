//! Polynomial systems by lexicographic Gröbner bases (CAS-DOC-21 §4,
//! CAS-289): `Solve[{eqs}, {x₁, …, xₙ}]`, `n ≥ 2`, for polynomial equations
//! with rational coefficients in the unknowns that the resultant of two
//! equations (CAS-215) does not serve.
//!
//! - The `plex` basis of the equations, `x₁ > … > xₙ` (cas-poly's
//!   `groebner_basis`): `{1}` means no solution, `{}`.
//! - Zero-dimensional when for each unknown some leading monomial is a
//!   power of it alone; otherwise the solutions are not finitely many and
//!   the call stays.
//! - From `xₙ` up, the basis is triangular: at each partial solution, the
//!   elements in `xᵢ, …, xₙ` that hold `xᵢ`, the later unknowns' values
//!   substituted, are polynomials in `xᵢ`; the roots of one of the least
//!   degree (by `Solve`'s `roots`, real ones for `Solve`) that are roots of
//!   the others too, by the check of CAS-213, extend it.
//! - Each tuple once, checked in the equations, by the unknowns' values
//!   ascending when all are real numbers: `{{a₁, …, aₙ}, …}`.

use cas_num::ball::Float;
use cas_poly::{Cancel, Distributed, GroebnerError, MonomialOrder, groebner_basis};

use super::super::approximation;
use super::super::groebner::{distributed, expression};
use super::transcendental::holds;
use super::{BITS, roots, sides};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The solutions of the polynomial system `equations` in `unknowns`, as
/// `{{a₁, …, aₙ}, …}`; `None` when it is not one, its solutions are not
/// finitely many, or they are not found here.
pub(super) fn system(
    ev: &mut Eval<'_>,
    e: &Expr,
    equations: &[Expr],
    unknowns: &[Expr],
    real: bool,
) -> Result<Option<Expr>, EvalError> {
    let n = unknowns.len();
    let mut residuals = Vec::with_capacity(equations.len());
    let mut polys = Vec::with_capacity(equations.len());
    for eq in equations {
        let f = match sides(eq) {
            Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
            None => eq.clone(),
        };
        let Some(p) = distributed(ev, &f, unknowns, MonomialOrder::Lex)? else {
            return Ok(None);
        };
        residuals.push(f);
        polys.push(p);
    }
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let basis = match groebner_basis(&polys, MonomialOrder::Lex, cancel) {
        Ok(basis) => basis,
        Err(GroebnerError::Cancelled) => return Err(EvalError::Interrupted),
        Err(GroebnerError::TooLarge) => return Ok(None),
    };
    let leads: Vec<&[u64]> = basis
        .iter()
        .filter_map(Distributed::leading_monomial)
        .collect();
    if leads.iter().any(|m| m.iter().all(|&d| d == 0)) {
        return Ok(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new())));
    }
    // Zero-dimensional: a pure power of each unknown leads an element.
    let pure = |i: usize| {
        leads.iter().any(|m| {
            m.iter()
                .enumerate()
                .all(|(j, &d)| if j == i { d > 0 } else { d == 0 })
        })
    };
    if basis.is_empty() || !(0..n).all(pure) {
        return Ok(None);
    }
    // Each element as an expression, with its first unknown and its degree
    // in it.
    let mut elements: Vec<(usize, u64, Expr)> = Vec::with_capacity(basis.len());
    for g in &basis {
        let first = (0..n)
            .find(|&i| g.terms().iter().any(|(m, _)| m[i] > 0))
            .unwrap_or(n);
        let degree = g
            .terms()
            .iter()
            .map(|(m, _)| m.get(first).copied().unwrap_or(0))
            .max()
            .unwrap_or(0);
        elements.push((first, degree, ev.eval(&expression(g, unknowns))?));
    }
    // Partial solutions, the values of xᵢ … xₙ, from the last.
    let mut partial: Vec<Vec<Expr>> = vec![Vec::new()];
    for i in (0..n).rev() {
        let Some(next) = extend(ev, e, &elements, unknowns, i, &partial, real)? else {
            return Ok(None);
        };
        partial = next;
    }
    let mut tuples: Vec<Vec<Expr>> = Vec::with_capacity(partial.len());
    for t in partial {
        let rules: Vec<(Expr, Expr)> = unknowns.iter().cloned().zip(t.iter().cloned()).collect();
        let ctx = ev.context();
        if residuals.iter().all(|f| holds(ctx, f, &rules)) && !tuples.contains(&t) {
            tuples.push(t);
        }
    }
    let tuples = ordered(ev.context(), tuples);
    Ok(Some(Expr::call(
        Symbol::LIST,
        tuples.into_iter().map(|t| Expr::call(Symbol::LIST, t)),
    )))
}

/// The partial solutions, values of `xᵢ₊₁ … xₙ`, extended by `xᵢ`: the
/// roots of the least-degree element whose first unknown is `xᵢ`, the
/// values substituted, that are roots of the others too. An element left
/// free of `xᵢ` that is not 0 ends its branch; `None` when roots are not
/// found.
fn extend(
    ev: &mut Eval<'_>,
    e: &Expr,
    elements: &[(usize, u64, Expr)],
    unknowns: &[Expr],
    i: usize,
    partial: &[Vec<Expr>],
    real: bool,
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let x = &unknowns[i];
    let mut here: Vec<&(usize, u64, Expr)> =
        elements.iter().filter(|(first, ..)| *first == i).collect();
    here.sort_by_key(|(_, degree, _)| *degree);
    let mut next = Vec::new();
    'branches: for values in partial {
        let rules: Vec<(Expr, Expr)> = unknowns[i + 1..]
            .iter()
            .cloned()
            .zip(values.iter().cloned())
            .collect();
        let mut hs = Vec::with_capacity(here.len());
        for (_, _, g) in &here {
            let h = ev.eval(&g.substitute(&rules))?;
            if h.contains(x) {
                hs.push(h);
            } else if !holds(ev.context(), &h, &[]) {
                continue 'branches;
            }
        }
        let Some((h, rest)) = hs.split_first() else {
            return Ok(None);
        };
        let Some(found) = roots(ev, e, h, x, real)? else {
            return Ok(None);
        };
        for r in found {
            if r.contains(x) {
                return Ok(None);
            }
            let rule = [(x.clone(), r.clone())];
            let ctx = ev.context();
            if rest.iter().all(|h| holds(ctx, h, &rule)) {
                let mut extended = Vec::with_capacity(values.len() + 1);
                extended.push(r);
                extended.extend(values.iter().cloned());
                next.push(extended);
            }
        }
    }
    Ok(Some(next))
}

/// The tuples by their values ascending, the first unknown's first, when
/// all are real numbers; otherwise as they are.
fn ordered(ctx: &Context, tuples: Vec<Vec<Expr>>) -> Vec<Vec<Expr>> {
    let mid = |e: &Expr| -> Option<Float> {
        match approximation(ctx, e, BITS).ok()?.view() {
            View::Ball(b) if b.is_finite() => Some(b.mid().clone()),
            _ => None,
        }
    };
    let keys: Option<Vec<Vec<Float>>> =
        tuples.iter().map(|t| t.iter().map(mid).collect()).collect();
    let Some(keys) = keys else {
        return tuples;
    };
    let mut keyed: Vec<(Vec<Float>, Vec<Expr>)> = keys.into_iter().zip(tuples).collect();
    keyed.sort_by(|(p, _), (q, _)| p.cmp(q));
    keyed.into_iter().map(|(_, t)| t).collect()
}

#[cfg(test)]
mod tests;
