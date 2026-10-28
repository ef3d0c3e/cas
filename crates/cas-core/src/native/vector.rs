//! Vector calculus (CAS-DOC-37): the differential operators (§2, CAS-408),
//! Xcas's `grad`, `divergence`, `curl`, `laplacian` and `hessian`.
//!
//! - `Grad[f, {x₁, …}]` is `{∂f/∂x₁, …}`, `Divergence[F, {x₁, …}]` is `Σ
//!   ∂F_k/∂x_k`, `Laplacian[f, {x₁, …}]` is `Σ ∂²f/∂x_k²` and
//!   `Hessian[f, {x₁, …}]` is `(∂²f/∂x_j ∂x_k)`, in any number of
//!   variables; `Curl[F, {x, y, z}]` is `{∂F₃/∂y − ∂F₂/∂z, ∂F₁/∂z −
//!   ∂F₃/∂x, ∂F₂/∂x − ∂F₁/∂y}`, in three.
//! - `Laplacian[n]`, for an integer `1 ≤ n ≤` [`MAX_SIZE`], is the `n × n`
//!   matrix of the discrete Laplacian, 2 on the diagonal and −1 beside it,
//!   as giac's `laplacian(3)`.
//!
//! The variables are distinct symbols, as many as the field's entries.
//! Otherwise the call stays, as `D`'s does for a variable that is not a
//! symbol. The derivatives are `D`'s, evaluated with the answer.
//!
//! - [`potential`]: `Potential` and `VectorPotential` (§3).
//! - [`curves`]: `Curvature`, `Frenet`, `Evolute` and `ArcLength` (§4).
//! - [`implicit`]: `ImplicitD` (§5).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of vector calculus: F = (p, q, r) in x, y and z"
)]

pub(super) mod curves;
pub(super) mod implicit;
pub(super) mod potential;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest discrete Laplacian built.
const MAX_SIZE: i64 = 1024;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn d(f: &Expr, x: &Expr) -> Expr {
    Expr::call(Symbol::D, [f.clone(), x.clone()])
}

/// The entries of a list, or `None`.
pub(super) fn entries(e: &Expr) -> Option<&[Expr]> {
    (e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)).then(|| e.args())
}

/// The variables of a call: a list of distinct symbols that are not
/// builtins.
pub(super) fn variables(e: &Expr) -> Option<&[Expr]> {
    let vars = entries(e)?;
    let symbol = |v: &Expr| v.as_symbol().is_some_and(|s| !s.is_builtin());
    let distinct = vars.iter().enumerate().all(|(i, v)| !vars[..i].contains(v));
    (!vars.is_empty() && vars.iter().all(symbol) && distinct).then_some(vars)
}

/// `Grad[f, vars]`.
pub(super) fn grad(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, vars] = e.args() else {
        return Ok(None);
    };
    let Some(vars) = variables(vars) else {
        return Ok(None);
    };
    if entries(f).is_some() {
        return Ok(None);
    }
    let result = list(vars.iter().map(|x| d(f, x)).collect());
    ev.step("grad", e, &result);
    Ok(Some(result))
}

/// `Divergence[F, vars]`.
pub(super) fn divergence(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [field, vars] = e.args() else {
        return Ok(None);
    };
    let (Some(field), Some(vars)) = (entries(field), variables(vars)) else {
        return Ok(None);
    };
    if field.len() != vars.len() {
        return Ok(None);
    }
    let terms: Vec<Expr> = field.iter().zip(vars).map(|(f, x)| d(f, x)).collect();
    let result = Expr::call(Symbol::PLUS, terms);
    ev.step("divergence", e, &result);
    Ok(Some(result))
}

/// `Curl[F, {x, y, z}]`.
pub(super) fn curl(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [field, vars] = e.args() else {
        return Ok(None);
    };
    let (Some([p, q, r]), Some([x, y, z])) = (entries(field), variables(vars)) else {
        return Ok(None);
    };
    let result = list(vec![
        d(r, y) - d(q, z),
        d(p, z) - d(r, x),
        d(q, x) - d(p, y),
    ]);
    ev.step("curl", e, &result);
    Ok(Some(result))
}

/// `Laplacian[f, vars]`, or `Laplacian[n]` for the discrete one.
pub(super) fn laplacian(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let result = match e.args() {
        [n] => {
            let Some(n) = n.as_i64().filter(|n| (1..=MAX_SIZE).contains(n)) else {
                return Ok(None);
            };
            let row = |j: i64| -> Expr {
                list(
                    (0..n)
                        .map(|k| {
                            Expr::from(match (j - k).abs() {
                                0 => 2,
                                1 => -1,
                                _ => 0,
                            })
                        })
                        .collect(),
                )
            };
            list((0..n).map(row).collect())
        }
        [f, vars] => {
            let Some(vars) = variables(vars) else {
                return Ok(None);
            };
            if entries(f).is_some() {
                return Ok(None);
            }
            Expr::call(Symbol::PLUS, vars.iter().map(|x| d(&d(f, x), x)))
        }
        _ => return Ok(None),
    };
    ev.step("laplacian", e, &result);
    Ok(Some(result))
}

/// `Hessian[f, vars]`.
pub(super) fn hessian(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, vars] = e.args() else {
        return Ok(None);
    };
    let Some(vars) = variables(vars) else {
        return Ok(None);
    };
    if entries(f).is_some() {
        return Ok(None);
    }
    let result = list(
        vars.iter()
            .map(|x| list(vars.iter().map(|y| d(&d(f, x), y)).collect()))
            .collect(),
    );
    ev.step("hessian", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
