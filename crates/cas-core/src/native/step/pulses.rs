//! Xcas's pulse shapes (CAS-648), ports of giac's signalprocessing.cc
//! (:2938-2975), built of `Sin` and the step, `HeavisideTheta[0] = 1/2` as
//! libgiac's `Heaviside(0)`:
//!
//! - `Sinc[g]`, giac's `sinc`: 1 when `g` is zero (`0` or `0.0`, giac's
//!   `is_zero`), else `Sin[g]/g`, so `Sinc[x]` is `Sin[x]/x` at once as
//!   giac's `sin(x)/x`, and exact input stays exact (`sinc(π/2)` = `2/π`).
//! - `UnitTriangle[g]`, giac's `tri`: `(1 − g)(θ(g) − θ(g − 1)) + (1 + g)(θ(−g)
//!   − θ(−g − 1))`, 1 at 0, `1/2` at `±1/2`, 0 at `±1` and beyond.
//! - `Boxcar[a, b, x]`, giac's `boxcar`: `θ(x − a) − θ(x − b)`, `1/2` at the
//!   ends, `−1` inside for `a > b` as giac's. Exactly three arguments, else
//!   the call stays with giac's warning "Bad Argument Type".
//!
//! Each builds giac's formula and evaluates it; on lists through their
//! `Listable` attribute (giac's `sinc([0,π,x])` is `[undef]`, its
//! `tri([0,1/2])` a `poly1`: departures).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn theta(u: Expr) -> Expr {
    Expr::call(Symbol::HEAVISIDE_THETA, [u])
}

/// giac's `is_zero` on a number: `0` or a decimal `0.0`.
fn is_zero(g: &Expr) -> bool {
    match g.view() {
        View::Int(n) => n == 0,
        View::Ball(b) => b.is_exact() && b.mid().is_zero(),
        _ => false,
    }
}

/// `Sinc[g]`.
pub(in crate::native) fn sinc(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [g] = e.args() else {
        return Ok(None);
    };
    if is_zero(g) {
        return Ok(Some(Expr::from(1)));
    }
    let value = Expr::call(Symbol::SIN, [g.clone()]) / g.clone();
    ev.eval(&value).map(Some)
}

/// `UnitTriangle[g]`.
pub(in crate::native) fn unit_triangle(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [g] = e.args() else {
        return Ok(None);
    };
    let one = || Expr::from(1);
    let value = (one() - g.clone()) * (theta(g.clone()) - theta(g.clone() - one()))
        + (one() + g.clone()) * (theta(-g.clone()) - theta(-g.clone() - one()));
    ev.eval(&value).map(Some)
}

/// `Boxcar[a, b, x]`.
pub(in crate::native) fn boxcar(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b, x] = e.args() else {
        ev.warn(Warning::Native {
            message: "Bad Argument Type".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    let value = theta(x.clone() - a.clone()) - theta(x.clone() - b.clone());
    ev.eval(&value).map(Some)
}

#[cfg(test)]
mod tests;
