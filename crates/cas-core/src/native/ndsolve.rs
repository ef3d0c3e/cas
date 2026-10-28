//! `NDSolve` (CAS-DOC-35, CAS-PLAN-29), Xcas's `odesolve`: initial value
//! problems `y′ = f(t, y)`, `y(t0) = y0`, integrated by Taylor's method on
//! balls, the answer an enclosure of `y(t1)` (§1, CAS-393):
//!
//! - `NDSolve[f, {t, y}, {t0, y0}, t1]`, giac's `odesolve(f, [t, y], [t0,
//!   y0], t1)`, is `{y(t1)}`, a list as giac's;
//! - `NDSolve[f, t == t0 .. t1, y, y0]`, giac's other order;
//! - a system: `f`, `y` and `y0` lists, `{{y1(t1), …}}`;
//! - a last argument `curve`: `{{t, {y}}, …}` at the steps taken, `t0`
//!   first and `t1` last.
//!
//! `t0`, `t1` and `y0` are real numbers, by `N`'s values; the right side is
//! compiled ([`rhs`]) and integrated ([`taylor`]) at the working precision
//! and 32 bits more, the answer kept to the working precision. A right side
//! that does not compile, or an integration that fails, leaves the call,
//! the latter with a warning.

mod rhs;
mod taylor;

use cas_num::ball::Ball;

use self::rhs::Node;
use super::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// What a call asks: the right sides, `t`, the unknowns, `t0`, `y0`, `t1`,
/// whether it is a system, and whether it asks for the curve.
struct Problem {
    f: Vec<Expr>,
    t: Expr,
    ys: Vec<Expr>,
    t0: Expr,
    y0: Vec<Expr>,
    t1: Expr,
    system: bool,
    curve: bool,
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The parts of a list, or the one expression.
fn items(e: &Expr) -> Vec<Expr> {
    if is_list(e) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    }
}

/// The call's problem; `None` for another form.
fn problem(e: &Expr) -> Option<Problem> {
    let args = e.args();
    let curve = args
        .last()
        .and_then(Expr::as_symbol)
        .is_some_and(|s| s.name() == "curve");
    let args = if curve { &args[..args.len() - 1] } else { args };
    let (f, t, ys, t0, y0, t1) = match args {
        // NDSolve[f, {t, y}, {t0, y0}, t1].
        [f, vars, start, t1] if is_list(vars) && is_list(start) => {
            let ([t, y], [t0, y0]) = (vars.args(), start.args()) else {
                return None;
            };
            (f, t, y, t0, y0, t1.clone())
        }
        // NDSolve[f, t == t0 .. t1, y, y0].
        [f, span, y, y0] => {
            let (Some(Symbol::EQUATION | Symbol::EQUAL), [t, range]) =
                (span.head().and_then(Expr::as_symbol), span.args())
            else {
                return None;
            };
            let (Some(Symbol::RANGE), [t0, t1]) =
                (range.head().and_then(Expr::as_symbol), range.args())
            else {
                return None;
            };
            (f, t, y, t0, y0, t1.clone())
        }
        _ => return None,
    };
    t.as_symbol()?;
    let system = is_list(ys);
    let (fs, unknowns, starts) = (items(f), items(ys), items(y0));
    if fs.len() != unknowns.len() || unknowns.len() != starts.len() || fs.is_empty() {
        return None;
    }
    if unknowns.iter().any(|u| u.as_symbol().is_none()) {
        return None;
    }
    Some(Problem {
        f: fs,
        t: t.clone(),
        ys: unknowns,
        t0: t0.clone(),
        y0: starts,
        t1,
        system,
        curve,
    })
}

/// A real ball of `e` at `prec` bits, by `N`'s value.
fn real(ev: &Eval<'_>, e: &Expr, prec: u64) -> Option<Ball> {
    if let Some(q) = e.to_rational() {
        return Some(Ball::from_rational(&q, prec));
    }
    let value = approximation(ev.context(), e, prec).ok()?;
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => value.to_rational().map(|q| Ball::from_rational(&q, prec)),
    }
}

/// `NDSolve[…]`.
pub(super) fn ndsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = problem(e) else {
        return Ok(None);
    };
    let prec = ev.context().working_precision();
    let q = prec + 32;
    let (Some(t0), Some(t1)) = (real(ev, &p.t0, q), real(ev, &p.t1, q)) else {
        return Ok(None);
    };
    let Some(y0) =
        p.y0.iter()
            .map(|y| real(ev, y, q))
            .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let mut variables = vec![p.t.clone()];
    variables.extend(p.ys.iter().cloned());
    let Some(f) =
        p.f.iter()
            .map(|fi| Node::compile(ev.context(), fi, &variables, q))
            .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let Some((y1, points)) = taylor::integrate(&f, &t0, &y0, &t1, prec) else {
        ev.warn(Warning::Native {
            message: "the integration failed: a solution that blows up, or too many steps".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    let kept = |b: &Ball| Expr::from(b.round_to_accuracy(prec));
    let value = |ys: &[Ball]| -> Expr {
        if p.system {
            Expr::call(Symbol::LIST, ys.iter().map(kept))
        } else {
            kept(&ys[0])
        }
    };
    if p.curve {
        let rows = points.iter().map(|(t, ys)| {
            Expr::call(
                Symbol::LIST,
                [kept(t), Expr::call(Symbol::LIST, ys.iter().map(kept))],
            )
        });
        return Ok(Some(Expr::call(Symbol::LIST, rows)));
    }
    Ok(Some(Expr::call(Symbol::LIST, [value(&y1)])))
}

#[cfg(test)]
mod tests;
