//! Implicit differentiation (CAS-DOC-37 §5, CAS-412): `ImplicitD`, giac's
//! `implicitdiff`.
//!
//! `ImplicitD[E, y, x]` reads `E` as `F(x, y) = 0`, an `Equation` or
//! `Equal` as its sides' difference and another expression as itself,
//! with `y` a function of `x` and every other symbol a constant:
//!
//! - `dy/dx = −F_x/F_y`, made `Together`: `implicitdiff(x^2+y^2=1,y,x)` is
//!   `-x/y`;
//! - the `n`-th derivative, `ImplicitD[E, y, x, …, x]` with `x` given `n`
//!   times or `ImplicitD[E, y, {x, n}]` as `D` takes it, for `1 ≤ n ≤`
//!   [`MAX_ORDER`]: the `(n−1)`-th differentiated in `x`, `y′` replaced by
//!   the first, made `Together`. giac writes the order as `x$2`, which this
//!   CAS's language has not.
//!
//! `x` and `y` are distinct symbols; otherwise the call stays.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The highest order.
const MAX_ORDER: i64 = 16;

fn d(f: &Expr, x: &Expr) -> Expr {
    Expr::call(Symbol::D, [f.clone(), x.clone()])
}

fn together(ev: &mut Eval<'_>, e: Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::TOGETHER, [e]))
}

/// The variable and the order of the arguments after `y`: `x` given `n`
/// times, or `{x, n}`.
fn order(rest: &[Expr]) -> Option<(&Expr, i64)> {
    match rest {
        [spec] if spec.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) => {
            let [x, n] = spec.args() else {
                return None;
            };
            let n = n.as_i64().filter(|n| (1..=MAX_ORDER).contains(n))?;
            Some((x, n))
        }
        [x, more @ ..] if more.iter().all(|v| v == x) => {
            let n = i64::try_from(rest.len()).ok()?;
            (n <= MAX_ORDER).then_some((x, n))
        }
        _ => None,
    }
}

/// `ImplicitD[E, y, x…]`.
pub(in crate::native) fn implicit_d(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [equation, y, rest @ ..] = e.args() else {
        return Ok(None);
    };
    let Some((x, n)) = order(rest) else {
        return Ok(None);
    };
    let symbol = |v: &Expr| v.as_symbol().is_some_and(|s| !s.is_builtin());
    if !symbol(x) || !symbol(y) || x == y {
        return Ok(None);
    }
    let f = match (equation.head().and_then(Expr::as_symbol), equation.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
        (Some(Symbol::LIST), _) => return Ok(None),
        _ => equation.clone(),
    };
    let first = together(ev, -(d(&f, x) / d(&f, y)))?;
    let mut derivative = first.clone();
    for _ in 1..n {
        let next = d(&derivative, x) + d(&derivative, y) * first.clone();
        derivative = together(ev, next)?;
    }
    ev.step("implicitdiff", e, &derivative);
    Ok(Some(derivative))
}

#[cfg(test)]
mod tests;
