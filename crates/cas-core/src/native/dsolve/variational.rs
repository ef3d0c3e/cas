//! The calculus of variations (CAS-575), ported from giac's
//! optimization.cc: `_euler_lagrange`, `_jacobi_equation` with
//! `parse_functional`, and `_conjugate_equation` with `strip_constants`.
//!
//! - `EulerLagrange[L]`, `[L, t]`, `[L, y[t]]`, `[L, {x[t], y[t]}]`,
//!   `[L, t, y]`, `[L, t, {u, v}]`: the unknowns `y` (default) of `t`
//!   (default `x`, or `t` when `L` has no `x`). In `L`, `y[t]` becomes the
//!   symbol `y` (a bare `y` also stands for it, as in giac) and `y′(t)` the
//!   fresh `EulerLagrange$du0`, as `DSolve`'s normal form does
//!   ([`super::derivative_symbol`]). The equations, giac's order: with one
//!   unknown and `L` free of `t`, the first integral `y′ L_y′ − L = K`;
//!   then for each unknown `L_y′ = K` when `L` is free of it, else
//!   `L_y − d/dt L_y′ = 0`. Each side is `Simplify`d (giac's `simp`, a
//!   recursive `normal`, then `radsimp`: `Together` alone leaves
//!   `(u² − (1 + u²))/√(1 + u²)`); the constants are `K_0`, `K_1`, … per
//!   call, skipping names of the call (giac counts across calls). One
//!   equation is not in a list.
//! - `d/dt L_y′` is the total derivative in every unknown; giac's
//!   differentiates only the unknown of the equation, so with coupled
//!   unknowns its equations differ. giac's solving for `y″` never happens
//!   (its test `solve_quiet(…).empty()` is inverted), so ours does not
//!   solve either.
//! - `JacobiEquation[L, t, y, Y, h, a]` or `[L, y[t], Y, h, a]`: with the
//!   second partials of `L` at the stationary `Y`, `P = L_y′y′`,
//!   `d = L_yy − d/dt L_yy′`, the equation `P′ h′ + P h″ − d h = 0` and
//!   `DSolve`'s solution with `h(a) = 0`, as a `Sequence` (the partials by
//!   `Together`, giac's `simp`); the equation alone when `DSolve` leaves
//!   it. With `P = 0`, whether `d = 0`.
//! - `ConjugateEquation[y, {α, β}, {α₀, β₀}, t, a]`: with `y₁ = ∂y/∂α`,
//!   `y₂ = ∂y/∂β`, `y₁ y₂(a) − y₂ y₁(a)` at the parameter values, `Together`
//!   then `FactorRational` (giac's `collect`), its numeric factors and sign
//!   dropped.
//!
//! What does not fit these forms stays.

#![allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    reason = "giac's names: L, t, y, u, h, a and the partials L_u, L_du, L_yy, L_ydy"
)]

use super::derivative_head;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::pattern::substitute;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// A symbol a user may name: not a builtin.
fn name(e: &Expr) -> Option<Symbol> {
    e.as_symbol().filter(|s| !s.is_builtin())
}

/// `f` and `s` of `f[s]`, both names.
fn applied(e: &Expr) -> Option<(Symbol, Symbol)> {
    let View::Apply { head, args: [s] } = e.view() else {
        return None;
    };
    Some((name(head)?, name(s)?))
}

/// `Derivative[k][u][t]`.
fn derivative(k: i64, u: Symbol, t: Symbol) -> Expr {
    let d = Expr::call(Symbol::DERIVATIVE, [Expr::from(k)]);
    Expr::apply(&Expr::apply(&d, &[u.into()]), &[t.into()])
}

/// `l` with `u[t]` as `u` and `u′[t]` as `du`, for each pair `(u, du)`;
/// `None` when another application of a `u` remains.
fn functional(l: &Expr, t: Symbol, unknowns: &[(Symbol, Expr)]) -> Option<Expr> {
    let t = Expr::from(t);
    let replaced = l.replace(|part| {
        let View::Apply { head, args: [a] } = part.view() else {
            return None;
        };
        if *a != t {
            return None;
        }
        unknowns.iter().find_map(|(u, du)| {
            let u = Expr::from(*u);
            if *head == u {
                Some(u)
            } else {
                (derivative_head(head, &u) == Some(1)).then(|| du.clone())
            }
        })
    });
    let left = replaced.distinct().any(|part| {
        let View::Apply { head, .. } = part.view() else {
            return false;
        };
        unknowns.iter().any(|(u, _)| {
            let u = Expr::from(*u);
            *head == u || derivative_head(head, &u).is_some()
        })
    });
    (!left).then_some(replaced)
}

/// `∂f/∂x`.
fn d(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))
}

/// giac's `simp` then `radsimp`: `Simplify`.
fn simp(ev: &mut Eval<'_>, f: Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::SIMPLIFY, [f]))
}

/// giac's `simp`, a recursive `normal`: `Together`.
fn normal(ev: &mut Eval<'_>, f: Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::TOGETHER, [f]))
}

/// The variable and unknowns of an `EulerLagrange` call, giac's defaults.
fn unknowns_of(e: &Expr) -> Option<(Expr, Symbol, Vec<Symbol>)> {
    let y = Symbol::new("y");
    let is_list = |a: &Expr| a.head() == Some(&Expr::from(Symbol::LIST));
    match e.args() {
        [l] => {
            let x = Symbol::new("x");
            let t = if l.contains(&x.into()) {
                x
            } else {
                Symbol::new("t")
            };
            Some((l.clone(), t, vec![y]))
        }
        [l, a] if is_list(a) => {
            let pairs: Option<Vec<_>> = a.args().iter().map(applied).collect();
            let pairs = pairs.filter(|p| !p.is_empty())?;
            let t = pairs[0].1;
            pairs
                .iter()
                .all(|p| p.1 == t)
                .then(|| (l.clone(), t, pairs.iter().map(|p| p.0).collect()))
        }
        [l, a] => match applied(a) {
            Some((u, t)) => Some((l.clone(), t, vec![u])),
            None => Some((l.clone(), name(a)?, vec![y])),
        },
        [l, t, u] => {
            let us = if is_list(u) {
                u.args().iter().map(name).collect::<Option<Vec<_>>>()?
            } else {
                vec![name(u)?]
            };
            let t = name(t)?;
            (!us.is_empty()).then(|| (l.clone(), t, us))
        }
        _ => None,
    }
}

/// `EulerLagrange[L, …]` (giac's `_euler_lagrange`).
pub(in crate::native) fn euler_lagrange(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((l, t, us)) = unknowns_of(e) else {
        return Ok(None);
    };
    let fresh = |what: &str, i: usize| Expr::symbol(&format!("EulerLagrange${what}{i}"));
    let du: Vec<Expr> = (0..us.len()).map(|i| fresh("du", i)).collect();
    let d2u: Vec<Expr> = (0..us.len()).map(|i| fresh("d2u", i)).collect();
    let pairs: Vec<(Symbol, Expr)> = us.iter().copied().zip(du.iter().cloned()).collect();
    let Some(l) = functional(&l, t, &pairs) else {
        return Ok(None);
    };
    let tx = Expr::from(t);
    let mut constants = (0u64..)
        // `$`: the user's capitalised names are escaped (cas-lang `escape_name`).
        .map(|k| Expr::symbol(&format!("$K_{k}")))
        .filter(|k| !e.contains(k));
    let mut next = || constants.next().expect("names without end");
    let mut equations = Vec::new();
    if us.len() == 1 && !l.contains(&tx) {
        let l_du = d(ev, &l, &du[0])?;
        let first = simp(ev, du[0].clone() * l_du - l.clone())?;
        equations.push(Expr::call(Symbol::EQUATION, [first, next()]));
    }
    for (i, u) in us.iter().enumerate() {
        let u = Expr::from(*u);
        let l_du = d(ev, &l, &du[i])?;
        if !l.contains(&u) {
            let side = simp(ev, l_du)?;
            equations.push(Expr::call(Symbol::EQUATION, [side, next()]));
            continue;
        }
        // d/dt L_u′: through t and every unknown and its derivative.
        let mut total = d(ev, &l_du, &tx)?;
        for (j, v) in us.iter().enumerate() {
            total = total
                + d(ev, &l_du, &Expr::from(*v))? * du[j].clone()
                + d(ev, &l_du, &du[j])? * d2u[j].clone();
        }
        let l_u = d(ev, &l, &u)?;
        let side = simp(ev, l_u - total)?;
        equations.push(Expr::call(Symbol::EQUATION, [side, Expr::from(0)]));
    }
    let mut back = Vec::new();
    for (i, u) in us.iter().enumerate() {
        back.push((*u, Expr::call(*u, [tx.clone()])));
        for (k, s) in [(1, &du[i]), (2, &d2u[i])] {
            back.push((s.as_symbol().expect("a symbol"), derivative(k, *u, t)));
        }
    }
    let result = match <[Expr; 1]>::try_from(equations) {
        Ok([one]) => one,
        Err(all) => Expr::call(Symbol::LIST, all),
    };
    ev.eval(&substitute(&result, &back)).map(Some)
}

/// `JacobiEquation[L, t, y, Y, h, a]` or `[L, y[t], Y, h, a]` (giac's
/// `_jacobi_equation`).
pub(in crate::native) fn jacobi_equation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, t, y, big_y, h, a) = match e.args() {
        [l, f, big_y, h, a] => {
            let Some((y, t)) = applied(f) else {
                return Ok(None);
            };
            (l, t, y, big_y, h, a)
        }
        [l, t, y, big_y, h, a] => {
            let (Some(t), Some(y)) = (name(t), name(y)) else {
                return Ok(None);
            };
            (l, t, y, big_y, h, a)
        }
        _ => return Ok(None),
    };
    let Some(h) = name(h) else {
        return Ok(None);
    };
    if !a.is_numeric() {
        return Ok(None);
    }
    let dy = Expr::symbol("JacobiEquation$dy");
    let Some(l) = functional(l, t, &[(y, dy.clone())]) else {
        return Ok(None);
    };
    let (tx, yx) = (Expr::from(t), Expr::from(y));
    let dy_t = d(ev, big_y, &tx)?;
    let at = [
        (y, big_y.clone()),
        (dy.as_symbol().expect("a symbol"), dy_t),
    ];
    // A second partial of L, at y = Y.
    let second = |ev: &mut Eval<'_>, p: &Expr, q: &Expr| -> Result<Expr, EvalError> {
        let first = d(ev, &l, p)?;
        let second = d(ev, &first, q)?;
        let second = normal(ev, second)?;
        ev.eval(&substitute(&second, &at))
    };
    let p = second(ev, &dy, &dy)?;
    let l_yy = second(ev, &yx, &yx)?;
    let l_ydy = second(ev, &yx, &dy)?;
    let l_ydy_t = d(ev, &l_ydy, &tx)?;
    let d_coefficient = normal(ev, l_yy - l_ydy_t)?;
    if p == Expr::from(0) {
        let zero = d_coefficient == Expr::from(0);
        return Ok(Some(Expr::from(if zero {
            Symbol::TRUE
        } else {
            Symbol::FALSE
        })));
    }
    let p_t = d(ev, &p, &tx)?;
    let p_t = normal(ev, p_t)?;
    let lhs = p_t * derivative(1, h, t) + p * derivative(2, h, t)
        - d_coefficient * Expr::call(h, [tx.clone()]);
    let lhs = ev.eval(&lhs)?;
    let equation = Expr::call(Symbol::EQUATION, [lhs, Expr::from(0)]);
    let condition = Expr::call(
        Symbol::EQUATION,
        [Expr::call(h, [a.clone()]), Expr::from(0)],
    );
    let problem = Expr::call(
        Symbol::DSOLVE,
        [
            Expr::call(Symbol::LIST, [equation.clone(), condition]),
            Expr::call(h, [tx]),
        ],
    );
    let solution = ev.eval(&problem)?;
    if solution.head() == Some(&Expr::from(Symbol::DSOLVE)) {
        return Ok(Some(equation));
    }
    Ok(Some(Expr::call(Symbol::SEQUENCE, [equation, solution])))
}

/// `g` without its numeric factors and sign (giac's `strip_constants`).
fn strip_constants(g: &Expr) -> Expr {
    match (g.head().and_then(Expr::as_symbol), g.args()) {
        (Some(Symbol::TIMES), factors) => Expr::call(
            Symbol::TIMES,
            factors
                .iter()
                .filter(|f| !f.is_numeric())
                .map(strip_constants)
                .collect::<Vec<_>>(),
        ),
        (Some(Symbol::POWER), [b, n]) if *n == Expr::from(-1) => strip_constants(b).pow(-1),
        _ => g.clone(),
    }
}

/// `ConjugateEquation[y, {α, β}, {α₀, β₀}, t, a]` (giac's
/// `_conjugate_equation`).
pub(in crate::native) fn conjugate_equation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [y, parameters, values, t, a] = e.args() else {
        return Ok(None);
    };
    let pair = |l: &Expr| {
        (l.head() == Some(&Expr::from(Symbol::LIST)))
            .then(|| <[Expr; 2]>::try_from(l.args().to_vec()).ok())
            .flatten()
    };
    let (Some([alpha, beta]), Some([alpha0, beta0]), Some(t)) =
        (pair(parameters), pair(values), name(t))
    else {
        return Ok(None);
    };
    let (Some(alpha_s), Some(beta_s)) = (name(&alpha), name(&beta)) else {
        return Ok(None);
    };
    if y.head().is_none() || !a.is_numeric() {
        return Ok(None);
    }
    let y1 = d(ev, y, &alpha)?;
    let y2 = d(ev, y, &beta)?;
    let at_a = |f: &Expr| substitute(f, &[(t, a.clone())]);
    let w = y1.clone() * at_a(&y2) - y2 * at_a(&y1);
    let w = substitute(&w, &[(alpha_s, alpha0), (beta_s, beta0)]);
    let w = normal(ev, w)?;
    let w = ev.eval(&Expr::call(Symbol::FACTOR_RATIONAL, [w]))?;
    ev.eval(&strip_constants(&w)).map(Some)
}
