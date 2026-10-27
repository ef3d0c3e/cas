//! Higher orders without `y` (CAS-DOC-16 §7, CAS-265): an equation in `x`
//! and `y⁽ᵏ⁾ … y⁽ⁿ⁾` only, `1 ≤ k < n`, is one of order `n − k` in
//! `p = y⁽ᵏ⁾`, solved by the other methods, its constants first; then `p`,
//! or each of its branches, is integrated `k` times, a constant after each:
//! `y″ = y′²` is `−ln|x − c_0| + c_1` (giac's `ln(x − c_0)`, without the
//! absolute value).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equation: p = y⁽ᵏ⁾, order n"
)]

use super::{Ode, solve};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The general solution of an equation without `y`, from that of `p`;
/// `None` when `y` occurs, or `p` is not solved.
pub(super) fn without_y(
    ev: &mut Eval<'_>,
    e: &Expr,
    ode: &Ode,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let n = ode.order;
    let Some(k) = (0..=n).find(|&j| ode.residual.contains(&ode.ys[j])) else {
        return Ok(None);
    };
    if k == 0 || k == n {
        return Ok(None);
    }
    // y⁽ᵏ⁺ʲ⁾ as p⁽ʲ⁾.
    let rename: Vec<(Expr, Expr)> = (k..=n)
        .map(|j| (ode.ys[j].clone(), ode.ys[j - k].clone()))
        .collect();
    let reduced = Ode {
        residual: ode.residual.substitute(&rename),
        ys: ode.ys[..=n - k].to_vec(),
        order: n - k,
    };
    let Some(p) = solve(ev, e, &reduced, x)? else {
        return Ok(None);
    };
    let integrated = |ev: &mut Eval<'_>, branch: &Expr| -> Result<Expr, EvalError> {
        let mut y = branch.clone();
        for c in &constants[n - k..n] {
            let integral = Expr::call(Symbol::INTEGRATE, [y, x.clone()]);
            y = ev.eval(&(integral + c.clone()))?;
        }
        Ok(y)
    };
    if p.head() == Some(&Expr::from(Symbol::LIST)) {
        let mut branches = Vec::with_capacity(p.args().len());
        for branch in p.args() {
            branches.push(integrated(ev, branch)?);
        }
        return Ok(Some(Expr::call(Symbol::LIST, branches)));
    }
    integrated(ev, &p).map(Some)
}
