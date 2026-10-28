//! Potentials (CAS-DOC-37 §3, CAS-409): `Potential`, giac's `potential`,
//! and `VectorPotential`, giac's `vpotential`.
//!
//! - `Potential[F, {x₁, …, x_n}]` is `U` with `grad U = F`, by successive
//!   integration: `U₁ = ∫ F₁ dx₁`, then `U_k = U_(k−1) + ∫ (F_k −
//!   ∂U_(k−1)/∂x_k) dx_k`, the integrand made free of `x₁ … x_(k−1)` by
//!   `Simplify`. The field is checked first: `∂F_j/∂x_k − ∂F_k/∂x_j` is 0
//!   by `Simplify` for every pair.
//! - `VectorPotential[F, {x, y, z}]` is `A` with `curl A = F` and `A₃ = 0`:
//!   `A₁ = ∫ F₂ dz`, and `A₂ = −∫ F₁ dz + ∫ φ dx`, where `φ = F₃ −
//!   (∂/∂x(−∫ F₁ dz) − ∂A₁/∂y)`, which is free of `z` when `div F = 0`. The
//!   field is checked first: `div F` is 0 by `Simplify`.
//!
//! A field that fails its check leaves the call with a warning, as giac's
//! "Not a potential"; so does a check `Simplify` cannot decide, or an
//! integral that stays.

use super::{entries, variables};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn d(f: &Expr, x: &Expr) -> Expr {
    Expr::call(Symbol::D, [f.clone(), x.clone()])
}

/// `Simplify[e]`, evaluated.
fn simplified(ev: &mut Eval<'_>, e: Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::SIMPLIFY, [e]))
}

/// `∫ f dx`, evaluated; `None` when it stays.
fn integral(ev: &mut Eval<'_>, f: Expr, x: &Expr) -> Result<Option<Expr>, EvalError> {
    let result = ev.eval(&Expr::call(Symbol::INTEGRATE, [f, x.clone()]))?;
    Ok((!result.contains(&Expr::from(Symbol::INTEGRATE))).then_some(result))
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// `Potential[F, vars]`.
pub(in crate::native) fn potential(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [field, vars] = e.args() else {
        return Ok(None);
    };
    let (Some(field), Some(vars)) = (entries(field), variables(vars)) else {
        return Ok(None);
    };
    if field.len() != vars.len() {
        return Ok(None);
    }
    for j in 0..vars.len() {
        for k in j + 1..vars.len() {
            let curl = d(&field[j], &vars[k]) - d(&field[k], &vars[j]);
            if simplified(ev, curl)? != Expr::from(0) {
                warn(ev, e, "not a gradient: its curl is not 0");
                return Ok(None);
            }
        }
    }
    let Some(mut u) = integral(ev, field[0].clone(), &vars[0])? else {
        return Ok(None);
    };
    for k in 1..vars.len() {
        let rest = simplified(ev, field[k].clone() - d(&u, &vars[k]))?;
        if vars[..k].iter().any(|x| rest.contains(x)) {
            warn(ev, e, "not a gradient, or not shown to be one");
            return Ok(None);
        }
        let Some(more) = integral(ev, rest, &vars[k])? else {
            return Ok(None);
        };
        u = ev.eval(&(u + more))?;
    }
    ev.step("potential", e, &u);
    Ok(Some(u))
}

/// `VectorPotential[F, {x, y, z}]`.
pub(in crate::native) fn vector_potential(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [field, vars] = e.args() else {
        return Ok(None);
    };
    let (Some([p, q, r]), Some([x, y, z])) = (entries(field), variables(vars)) else {
        return Ok(None);
    };
    let divergence = d(p, x) + d(q, y) + d(r, z);
    if simplified(ev, divergence)? != Expr::from(0) {
        warn(ev, e, "not a curl: its divergence is not 0");
        return Ok(None);
    }
    let (Some(first), Some(minus)) = (integral(ev, q.clone(), z)?, integral(ev, p.clone(), z)?)
    else {
        return Ok(None);
    };
    let second = -minus;
    let rest = simplified(ev, r.clone() - (d(&second, x) - d(&first, y)))?;
    if rest.contains(z) {
        warn(ev, e, "not a curl, or not shown to be one");
        return Ok(None);
    }
    let Some(more) = integral(ev, rest, x)? else {
        return Ok(None);
    };
    let result = Expr::call(Symbol::LIST, [first, second + more, Expr::from(0)]);
    ev.step("vpotential", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
