//! `DSolve` (CAS-DOC-16, CAS-259), Xcas's `desolve`: ordinary differential
//! equations in one unknown function of one variable.
//!
//! - `DSolve[eq, y[x]]`, giac's `DSolve[eq, x, y]` with the variable first,
//!   and `DSolve[eq, y]` with `x` the one argument `y` has in `eq`. The
//!   equation is `l = r` (`Equation` or `Equal`) or an expression meaning
//!   `= 0`; its derivatives are `Derivative[n][y][x]`, which `D` makes of
//!   `diff(y(x), x)` (§1).
//! - The normal form (§3): `y[x]` and its derivatives become the symbols
//!   `DSolve$y0`, `DSolve$y1`, …, `Exp[u]` becomes `E^u`, and the residual
//!   `l − r` is cleared of a denominator free of them; another occurrence of
//!   `y` leaves the call.
//! - A linear equation `Σ aₖ y⁽ᵏ⁾ = f` is one whose residual has
//!   derivatives in the `yₖ` free of all of them.
//! - First order, linear (§4.1): `y′ + p y = q` is `y = (∫ μ q dx + c_0)/μ`
//!   with `μ = e^(∫p dx)`, a logarithm of an absolute value in `∫p` taken
//!   without it, since any constant multiple of `μ` serves (`y′ = 2y/x` is
//!   `c_0 x²`, as giac's); a term `c ln u` of `∫p` makes `u^c` of `μ`.
//! - The answer is expanded; the constants `c_0`, `c_1`, … skip the names of
//!   the call. What no method solves stays, with a warning.
//! - `DSolve[{eq, conditions…}, y[x]]` fixes the constants by the
//!   conditions ([`conditions`], §8).
//! - `DSolve[{eq₁, …}, {y₁[x], …}]`, a list of unknowns, is a system
//!   ([`system`], CAS-DOC-20).
//!
//! Each call that solves takes one step, `"DSolve"`.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `DSolve[eq, y[x]]`, `DSolve[eq, x, y]` or `DSolve[eq, y]`; with a list
/// of unknowns, a system ([`system`], CAS-DOC-20).
pub(super) fn dsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let solution = if let Some(problem) = system::problem(e) {
        system::solve(ev, e, &problem)?
    } else if let Some((eq, y, x)) = problem(e) {
        one(ev, e, eq, &y, &x)?
    } else {
        return Ok(None);
    };
    let Some(solution) = solution else {
        ev.warn(Warning::Native {
            message: "unable to solve the differential equation".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    ev.step("DSolve", e, &solution);
    Ok(Some(solution))
}

/// The solution of one equation `eq`, or of `{eq, conditions…}`, for `y`
/// of `x`; `None` when it is not solved here.
fn one(
    ev: &mut Eval<'_>,
    e: &Expr,
    eq: Expr,
    y: &Expr,
    x: &Expr,
) -> Result<Option<Expr>, EvalError> {
    // Conditions come in a list after the equation (§8).
    let (eq, conditions) = match (eq.head() == Some(&Expr::from(Symbol::LIST)), eq.args()) {
        (true, [first, rest @ ..]) => (first.clone(), rest.to_vec()),
        (true, []) => return Ok(None),
        _ => (eq, Vec::new()),
    };
    let Some(ode) = Ode::of(ev, &eq, y, x)? else {
        return Ok(None);
    };
    match solve(ev, e, &ode, x)? {
        Some(general) if !conditions.is_empty() => {
            let constants = constants(e, ode.order);
            conditions::apply(ev, e, &general, &conditions, y, x, &constants)
        }
        general => Ok(general),
    }
}

/// The general solution of `ode`, by the first method of §4–§7 that applies.
fn solve(ev: &mut Eval<'_>, e: &Expr, ode: &Ode, x: &Expr) -> Result<Option<Expr>, EvalError> {
    let constants = constants(e, ode.order);
    if let Some((a, f)) = ode.linear(ev)? {
        if let [a0, a1] = a.as_slice() {
            return first_order_linear(ev, x, a0, a1, &f, &constants[0]).map(Some);
        }
        if let Some(y) = linear::constant_coefficients(ev, e, &a, &f, x, &constants)? {
            return Ok(Some(y));
        }
        if let Some(y) = euler::euler(ev, e, &a, &f, x, &constants)? {
            return Ok(Some(y));
        }
    } else if ode.order == 1 {
        return first::nonlinear(ev, ode, x, &constants[0]);
    }
    reduced::without_y(ev, e, ode, x, &constants)
}

/// The equation, the unknown `y` and the variable `x` of a call.
fn problem(e: &Expr) -> Option<(Expr, Expr, Expr)> {
    let unknown = |y: &Expr| y.as_symbol().is_some_and(|s| !s.is_builtin());
    let variable = |x: &Expr| x.as_symbol().is_some_and(|s| !s.is_builtin());
    match e.args() {
        [eq, y] if unknown(y) => {
            // The equation's argument, not the conditions' points.
            let ode = match eq.args() {
                [first, ..] if eq.head() == Some(&Expr::from(Symbol::LIST)) => first,
                _ => eq,
            };
            let x = argument_of(ode, y)?;
            Some((eq.clone(), y.clone(), x))
        }
        [eq, f] => {
            let (y, [x]) = (f.head()?, f.args()) else {
                return None;
            };
            (unknown(y) && variable(x)).then(|| (eq.clone(), y.clone(), x.clone()))
        }
        [eq, x, y] if variable(x) && unknown(y) => Some((eq.clone(), y.clone(), x.clone())),
        _ => None,
    }
}

/// The one argument `y` is applied to in `eq`, a symbol, as `y[x]` or
/// `Derivative[n][y][x]`.
fn argument_of(eq: &Expr, y: &Expr) -> Option<Expr> {
    let mut found: Option<&Expr> = None;
    for part in eq.distinct() {
        let View::Apply { head, args: [a] } = part.view() else {
            continue;
        };
        if head != y && derivative_head(head, y).is_none() {
            continue;
        }
        if found.is_some_and(|f| f != a) {
            return None;
        }
        found = Some(a);
    }
    found
        .filter(|x| x.as_symbol().is_some_and(|s| !s.is_builtin()))
        .cloned()
}

/// `n` for `Derivative[n][y]`, `n ≥ 1`.
fn derivative_head(head: &Expr, y: &Expr) -> Option<usize> {
    let View::Apply { head: d, args: [f] } = head.view() else {
        return None;
    };
    let View::Apply { head: d, args: [n] } = d.view() else {
        return None;
    };
    (d.as_symbol() == Some(Symbol::DERIVATIVE) && f == y)
        .then(|| usize::try_from(n.as_i64()?).ok())
        .flatten()
        .filter(|&n| n >= 1)
}

/// `y⁽ᵏ⁾(x)` in the normal form: `DSolve$yk`.
fn derivative_symbol(k: usize) -> Expr {
    Expr::from(Symbol::new(&format!("DSolve$y{k}")))
}

/// The constants `c_0`, `c_1`, …, `count` of them, skipping names of `e`.
/// Also `RSolve`'s (CAS-310).
pub(super) fn constants(e: &Expr, count: usize) -> Vec<Expr> {
    (0u64..)
        .map(|k| Expr::from(Symbol::new(&format!("c_{k}"))))
        .filter(|c| !e.contains(c))
        .take(count.max(1))
        .collect()
}

/// An equation in the normal form of §3: `residual = 0` in `x` and the
/// symbols `ys[k]` for `y⁽ᵏ⁾(x)`, `k ≤ order`.
struct Ode {
    residual: Expr,
    ys: Vec<Expr>,
    order: usize,
}

impl Ode {
    /// `eq` in the normal form, for the unknown `y` of `x`; `None` when `y`
    /// occurs otherwise, or not differentiated.
    fn of(ev: &mut Eval<'_>, eq: &Expr, y: &Expr, x: &Expr) -> Result<Option<Ode>, EvalError> {
        let residual = match (eq.head().and_then(Expr::as_symbol), eq.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
            _ => eq.clone(),
        };
        // Exp[u] as E^u, which products merge: x Exp[x] E^(−x) would stay.
        let residual = residual.map_bottom_up(|part| {
            match (part.head().and_then(Expr::as_symbol), part.args()) {
                (Some(Symbol::EXP), [u]) => Expr::from(Symbol::E).pow(u.clone()),
                _ => part.clone(),
            }
        });
        let mut order = 0;
        let replaced = residual.replace(|part| {
            let View::Apply { head, args: [a] } = part.view() else {
                return None;
            };
            if a != x {
                return None;
            }
            let k = if head == y {
                0
            } else {
                derivative_head(head, y)?
            };
            order = order.max(k);
            Some(derivative_symbol(k))
        });
        if order == 0 || replaced.contains(y) {
            return Ok(None);
        }
        let ys: Vec<Expr> = (0..=order).map(derivative_symbol).collect();
        let together = ev.eval(&Expr::call(Symbol::TOGETHER, [replaced]))?;
        let denominator = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together.clone()]))?;
        let residual = if ys.iter().any(|y| denominator.contains(y)) {
            together
        } else {
            ev.eval(&Expr::call(Symbol::NUMERATOR, [together]))?
        };
        Ok(Some(Ode {
            residual,
            ys,
            order,
        }))
    }

    /// The coefficients `a₀ … aₙ` and the right side `f` of a linear
    /// equation `Σ aₖ y⁽ᵏ⁾ = f`, all free of the `yₖ`; `None` when it is not
    /// linear, or its leading coefficient is 0.
    fn linear(&self, ev: &mut Eval<'_>) -> Result<Option<(Vec<Expr>, Expr)>, EvalError> {
        let mut a = Vec::with_capacity(self.ys.len());
        for yk in &self.ys {
            let ak = ev.eval(&Expr::call(Symbol::D, [self.residual.clone(), yk.clone()]))?;
            if self.ys.iter().any(|y| ak.contains(y)) {
                return Ok(None);
            }
            a.push(ak);
        }
        if a.last() == Some(&Expr::from(0)) {
            return Ok(None);
        }
        let zero: Vec<(Expr, Expr)> = self.ys.iter().map(|y| (y.clone(), Expr::from(0))).collect();
        let f = ev.eval(&-self.residual.substitute(&zero))?;
        Ok(Some((a, f)))
    }
}

/// `a₁ y′ + a₀ y = f`: `y = (∫ μ q dx + c)/μ`, `μ = e^(∫ a₀/a₁ dx)` and
/// `q = f/a₁` (§4.1).
#[allow(
    clippy::many_single_char_names,
    reason = "the names of the equation and its solution"
)]
fn first_order_linear(
    ev: &mut Eval<'_>,
    x: &Expr,
    a0: &Expr,
    a1: &Expr,
    f: &Expr,
    c: &Expr,
) -> Result<Expr, EvalError> {
    let p = ev.eval(&(a0.clone() / a1.clone()))?;
    let q = ev.eval(&(f.clone() / a1.clone()))?;
    let big_p = ev.eval(&Expr::call(Symbol::INTEGRATE, [p, x.clone()]))?;
    let mu = exponential(ev, &without_absolute_values(&big_p))?;
    let integral = ev.eval(&Expr::call(Symbol::INTEGRATE, [mu.clone() * q, x.clone()]))?;
    ev.eval(&Expr::call(Symbol::EXPAND, [(integral + c.clone()) / mu]))
}

/// `e` with `Log[Abs[u]]` as `Log[u]`: integrating factors up to a constant
/// factor, which their sign is on each interval.
fn without_absolute_values(e: &Expr) -> Expr {
    e.map_bottom_up(
        |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(Symbol::LOG), [u]) if u.head() == Some(&Expr::from(Symbol::ABS)) => {
                Expr::call(Symbol::LOG, [u.args()[0].clone()])
            }
            _ => part.clone(),
        },
    )
}

/// `e^p`, a term `c Log[u]` of `p` (`c` a rational) as `u^c`: `x²` for
/// `2 Log[x]`.
fn exponential(ev: &mut Eval<'_>, p: &Expr) -> Result<Expr, EvalError> {
    let terms = if p.head() == Some(&Expr::from(Symbol::PLUS)) {
        p.args().to_vec()
    } else {
        vec![p.clone()]
    };
    let factors: Vec<Expr> = terms
        .iter()
        .map(|t| match logarithm_multiple(t) {
            Some((c, u)) => u.pow(c),
            None => Expr::from(Symbol::E).pow(t.clone()),
        })
        .collect();
    ev.eval(&Expr::call(Symbol::TIMES, factors))
}

/// `(c, u)` for a term `c Log[u]` with `c` a rational, or `Log[u]`.
fn logarithm_multiple(t: &Expr) -> Option<(Expr, Expr)> {
    let log = |e: &Expr| match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::LOG), [u]) => Some(u.clone()),
        _ => None,
    };
    if let Some(u) = log(t) {
        return Some((Expr::from(1), u));
    }
    match (t.head().and_then(Expr::as_symbol), t.args()) {
        (Some(Symbol::TIMES), [c, l]) if c.to_rational().is_some() => {
            log(l).map(|u| (c.clone(), u))
        }
        _ => None,
    }
}

mod conditions;
mod euler;
mod first;
mod linear;
mod particular;
mod reduced;
mod system;

#[cfg(test)]
mod tests;
