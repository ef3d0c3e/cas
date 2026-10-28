//! Systems of differential equations (CAS-DOC-20, CAS-283): linear systems
//! with constant coefficients in several unknown functions of one
//! variable, by Laplace transforms.
//!
//! - `DSolve[{eq₁, …, eqₙ}, {y₁[x], …, yₙ[x]}]`, `DSolve[{…}, x, {y₁, …}]`
//!   and `DSolve[{…}, {y₁, …}]`, `x` the one argument the unknowns have in
//!   the first item (§1). The items with an unknown at `x` are the
//!   equations, as many as the unknowns; the others are conditions.
//! - The normal form (§2): each residual `l − r`, the unknowns and their
//!   derivatives at `x` taken as symbols, is affine in them with
//!   coefficients free of `x`; `kⱼ` is the highest order of `yⱼ`, and the
//!   coefficients of the `yⱼ^(kⱼ)` make a nonsingular matrix.
//! - Laplace transforms (§3): `L[yⱼ^(k)] = s^k Yⱼ − Σ_(m<k) s^(k−1−m)
//!   yⱼ^(m)(0)`, the initial values the constants `c_0, c_1, …` (`y₁(0)`,
//!   `y₁′(0)`, …, `y₂(0)`, …), the forcing's by `LaplaceTransform`; the
//!   system in the `Yⱼ` by `LinSolve`, and each `Yⱼ`, split by the
//!   constants, inverted part by part by `InverseLaplaceTransform`.
//! - The answer is giac's for the matrix form `desolve(y'=A*y+b,y)`: the
//!   list of the one solution, `{{v₁, …, vₙ}}`, each component expanded.
//! - Conditions (§4, CAS-284): the unknowns and their derivatives at the
//!   conditions' points, taken from the components, make linear equations
//!   in the constants (`LinSolve`); no solution gives `{}`, as giac's `[]`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the system: a, f, s, x, yⱼ^(k)"
)]

use super::super::solve::linear_system;
use super::{constants, derivative_head};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// A system: its equations and conditions, its unknowns and its variable.
pub(super) struct Problem {
    equations: Vec<Expr>,
    conditions: Vec<Expr>,
    ys: Vec<Expr>,
    x: Expr,
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// A symbol that may name an unknown or the variable.
fn name(e: &Expr) -> bool {
    e.as_symbol().is_some_and(|s| !s.is_builtin())
}

/// The system of a call whose unknowns are a list; `None` for another call.
pub(super) fn problem(e: &Expr) -> Option<Problem> {
    let (items, x, ys) = match e.args() {
        [items, x, ys] if name(x) && is_list(ys) && ys.args().iter().all(name) => {
            (items, Some(x.clone()), ys.args().to_vec())
        }
        [items, ys] if is_list(ys) && ys.args().iter().all(name) => {
            (items, None, ys.args().to_vec())
        }
        [items, fs] if is_list(fs) => {
            // y₁[x], …, yₙ[x], one variable.
            let mut x: Option<&Expr> = None;
            let mut ys = Vec::with_capacity(fs.args().len());
            for f in fs.args() {
                let View::Apply { head, args: [a] } = f.view() else {
                    return None;
                };
                if !name(head) || !name(a) || x.is_some_and(|x| x != a) {
                    return None;
                }
                x = Some(a);
                ys.push(head.clone());
            }
            (items, x.cloned(), ys)
        }
        _ => return None,
    };
    let distinct = ys.iter().enumerate().all(|(i, y)| !ys[..i].contains(y));
    if ys.is_empty() || !distinct || !is_list(items) {
        return None;
    }
    let x = match x {
        Some(x) => x,
        None => argument(items.args().first()?, &ys)?,
    };
    if ys.contains(&x) {
        return None;
    }
    let (equations, conditions) = items
        .args()
        .iter()
        .cloned()
        .partition(|item| at_variable(item, &ys, &x));
    Some(Problem {
        equations,
        conditions,
        ys,
        x,
    })
}

/// Whether `head` is one of `ys` or a derivative of one.
fn unknown_head(head: &Expr, ys: &[Expr]) -> bool {
    ys.iter()
        .any(|y| head == y || derivative_head(head, y).is_some())
}

/// The one argument the unknowns are applied to in `item`, a symbol.
fn argument(item: &Expr, ys: &[Expr]) -> Option<Expr> {
    let mut found: Option<&Expr> = None;
    for part in item.distinct() {
        let View::Apply { head, args: [a] } = part.view() else {
            continue;
        };
        if !unknown_head(head, ys) {
            continue;
        }
        if found.is_some_and(|f| f != a) {
            return None;
        }
        found = Some(a);
    }
    found.filter(|x| name(x)).cloned()
}

/// Whether an unknown, or a derivative of one, is applied to `x` in `item`.
fn at_variable(item: &Expr, ys: &[Expr], x: &Expr) -> bool {
    item.distinct().any(|part| match part.view() {
        View::Apply { head, args: [a] } => a == x && unknown_head(head, ys),
        _ => false,
    })
}

/// The solution of the system, `{{v₁, …, vₙ}}`; `None` when it is not
/// solved here.
pub(super) fn solve(
    ev: &mut Eval<'_>,
    e: &Expr,
    problem: &Problem,
) -> Result<Option<Expr>, EvalError> {
    if problem.equations.len() != problem.ys.len() {
        return Ok(None);
    }
    let Some(system) = Linear::of(ev, problem)? else {
        return Ok(None);
    };
    let count: usize = system.orders.iter().sum();
    let constants = constants(e, count);
    let constants = &constants[..count];
    let Some(components) = laplace(ev, e, &system, &problem.x, constants)? else {
        return Ok(None);
    };
    if problem.conditions.is_empty() {
        return Ok(Some(solution(components)));
    }
    fix(ev, e, &components, problem, constants)
}

/// `{{v₁, …, vₙ}}`.
fn solution(components: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, [Expr::call(Symbol::LIST, components)])
}

/// The general solution's `components` with their constants fixed by the
/// conditions (§4); `{}` when no solution meets them; `None` when a
/// condition is not one, at a point free of `x`.
fn fix(
    ev: &mut Eval<'_>,
    e: &Expr,
    components: &[Expr],
    problem: &Problem,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let (ys, x) = (&problem.ys, &problem.x);
    // The components' derivatives, as far as the conditions ask.
    let mut derivatives: Vec<Vec<Expr>> = components.iter().map(|v| vec![v.clone()]).collect();
    let mut equations = Vec::with_capacity(problem.conditions.len());
    for condition in &problem.conditions {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) =
            (condition.head().and_then(Expr::as_symbol), condition.args())
        else {
            return Ok(None);
        };
        for part in condition.distinct() {
            let View::Apply { head, args: [_] } = part.view() else {
                continue;
            };
            for (j, y) in ys.iter().enumerate() {
                let Some(k) = (if head == y {
                    Some(0)
                } else {
                    derivative_head(head, y)
                }) else {
                    continue;
                };
                while derivatives[j].len() <= k {
                    let last = derivatives[j].last().expect("the component").clone();
                    derivatives[j].push(ev.eval(&Expr::call(Symbol::D, [last, x.clone()]))?);
                }
            }
        }
        let at = (l.clone() - r.clone()).replace(|part| {
            let View::Apply { head, args: [a] } = part.view() else {
                return None;
            };
            ys.iter().enumerate().find_map(|(j, y)| {
                let k = if head == y {
                    0
                } else {
                    derivative_head(head, y)?
                };
                Some(derivatives[j][k].substitute(&[(x.clone(), a.clone())]))
            })
        });
        if at.contains(x) || ys.iter().any(|y| at.contains(y)) {
            return Ok(None);
        }
        let value = ev.eval(&at)?;
        equations.push(Expr::call(Symbol::EQUATION, [value, Expr::from(0)]));
    }
    let used: Vec<Expr> = constants
        .iter()
        .filter(|c| components.iter().any(|v| v.contains(c)))
        .cloned()
        .collect();
    if used.is_empty() {
        return Ok(None);
    }
    let Some(found) = linear_system(ev, e, &equations, &used)? else {
        return Ok(None);
    };
    match found.args() {
        [] => Ok(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new()))),
        [values] if values.args().len() == used.len() => {
            let rules: Vec<(Expr, Expr)> = used
                .iter()
                .cloned()
                .zip(values.args().iter().cloned())
                .collect();
            let mut fixed = Vec::with_capacity(components.len());
            for v in components {
                fixed.push(ev.eval(&Expr::call(Symbol::EXPAND, [v.substitute(&rules)]))?);
            }
            Ok(Some(solution(fixed)))
        }
        _ => Ok(None),
    }
}

/// A linear system with constant coefficients: equation `i` is
/// `Σⱼₖ a[i][j][k] yⱼ^(k) = f[i]`.
struct Linear {
    a: Vec<Vec<Vec<Expr>>>,
    f: Vec<Expr>,
    /// `kⱼ`, the highest order of `yⱼ`.
    orders: Vec<usize>,
}

/// `yⱼ^(k)(x)` in the normal form: `DSolve$y{j}d{k}`.
fn derivative_symbol(j: usize, k: usize) -> Expr {
    Expr::from(Symbol::new(&format!("DSolve$y{j}d{k}")))
}

impl Linear {
    /// The system of `problem` in the normal form of §2; `None` when it is
    /// not linear with constant coefficients and nonsingular.
    fn of(ev: &mut Eval<'_>, problem: &Problem) -> Result<Option<Linear>, EvalError> {
        let (ys, x) = (&problem.ys, &problem.x);
        let n = ys.len();
        let mut orders = vec![0usize; n];
        let mut residuals = Vec::with_capacity(n);
        for eq in &problem.equations {
            let residual = match (eq.head().and_then(Expr::as_symbol), eq.args()) {
                (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
                _ => eq.clone(),
            };
            let replaced = residual.replace(|part| {
                let View::Apply { head, args: [a] } = part.view() else {
                    return None;
                };
                if a != x {
                    return None;
                }
                let (j, k) = ys.iter().enumerate().find_map(|(j, y)| {
                    let k = if head == y {
                        0
                    } else {
                        derivative_head(head, y)?
                    };
                    Some((j, k))
                })?;
                orders[j] = orders[j].max(k);
                Some(derivative_symbol(j, k))
            });
            if ys.iter().any(|y| replaced.contains(y)) {
                return Ok(None);
            }
            residuals.push(ev.eval(&replaced)?);
        }
        if orders.iter().all(|&k| k == 0) {
            return Ok(None);
        }
        let symbols: Vec<Expr> = (0..n)
            .flat_map(|j| (0..=orders[j]).map(move |k| derivative_symbol(j, k)))
            .collect();
        let mut a = Vec::with_capacity(n);
        let mut f = Vec::with_capacity(n);
        for residual in &residuals {
            let mut row = Vec::with_capacity(n);
            for (j, &order) in orders.iter().enumerate() {
                let mut coefficients = Vec::with_capacity(order + 1);
                for k in 0..=order {
                    let v = derivative_symbol(j, k);
                    let c = ev.eval(&Expr::call(Symbol::D, [residual.clone(), v]))?;
                    if c.contains(x) || symbols.iter().any(|s| c.contains(s)) {
                        return Ok(None);
                    }
                    coefficients.push(c);
                }
                row.push(coefficients);
            }
            a.push(row);
            let zero: Vec<(Expr, Expr)> =
                symbols.iter().map(|s| (s.clone(), Expr::from(0))).collect();
            f.push(ev.eval(&-residual.substitute(&zero))?);
        }
        // The highest derivatives' coefficients, nonsingular.
        let leading: Vec<Expr> = a
            .iter()
            .map(|row: &Vec<Vec<Expr>>| {
                Expr::call(
                    Symbol::LIST,
                    row.iter().zip(&orders).map(|(c, &k)| c[k].clone()),
                )
            })
            .collect();
        let det = ev.eval(&Expr::call(
            Symbol::DET,
            [Expr::call(Symbol::LIST, leading)],
        ))?;
        if det == Expr::from(0) || det.head() == Some(&Expr::from(Symbol::DET)) {
            return Ok(None);
        }
        Ok(Some(Linear { a, f, orders }))
    }
}

/// The components of the general solution of `system` by Laplace
/// transforms (§3), the initial values `constants`, in order.
fn laplace(
    ev: &mut Eval<'_>,
    e: &Expr,
    system: &Linear,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Vec<Expr>>, EvalError> {
    let n = system.orders.len();
    let s = Expr::from(Symbol::new("DSolve$s"));
    let transforms: Vec<Expr> = (0..n)
        .map(|j| Expr::from(Symbol::new(&format!("DSolve$Y{j}"))))
        .collect();
    // The constant of yⱼ^(m)(0) is initial[j][m].
    let mut rest = constants.iter();
    let initial: Vec<Vec<Expr>> = system
        .orders
        .iter()
        .map(|&k| rest.by_ref().take(k).cloned().collect())
        .collect();
    let mut equations = Vec::with_capacity(n);
    for (row, f) in system.a.iter().zip(&system.f) {
        let mut left = Vec::new();
        let mut right = Vec::new();
        if *f != Expr::from(0) {
            let transform = ev.eval(&Expr::call(
                Symbol::LAPLACE_TRANSFORM,
                [f.clone(), x.clone(), s.clone()],
            ))?;
            if transform.contains(&Expr::from(Symbol::LAPLACE_TRANSFORM)) {
                return Ok(None);
            }
            right.push(transform);
        }
        for (j, coefficients) in row.iter().enumerate() {
            for (k, c) in coefficients.iter().enumerate() {
                if *c == Expr::from(0) {
                    continue;
                }
                let power = |p: usize| s.pow(i64::try_from(p).expect("a small order"));
                left.push(c.clone() * power(k) * transforms[j].clone());
                for (m, value) in initial[j].iter().enumerate().take(k) {
                    right.push(c.clone() * power(k - 1 - m) * value.clone());
                }
            }
        }
        equations.push(Expr::call(
            Symbol::EQUATION,
            [
                Expr::call(Symbol::PLUS, left),
                Expr::call(Symbol::PLUS, right),
            ],
        ));
    }
    let Some(found) = linear_system(ev, e, &equations, &transforms)? else {
        return Ok(None);
    };
    let [values] = found.args() else {
        return Ok(None);
    };
    if values.args().len() != n {
        return Ok(None);
    }
    let mut components = Vec::with_capacity(n);
    for v in values.args() {
        // Yⱼ = Hⱼ + Σ c Gⱼc, linear in the constants: Gⱼc its derivative in
        // c, which no difference leaves uncancelled; each part inverted.
        let zero: Vec<(Expr, Expr)> = constants
            .iter()
            .map(|c| (c.clone(), Expr::from(0)))
            .collect();
        let base = ev.eval(&v.substitute(&zero))?;
        let Some(h) = inverse(ev, &base, &s, x)? else {
            return Ok(None);
        };
        let mut terms = vec![h];
        for c in constants {
            let part = ev.eval(&Expr::call(Symbol::D, [v.clone(), c.clone()]))?;
            let Some(g) = inverse(ev, &part, &s, x)? else {
                return Ok(None);
            };
            terms.push(c.clone() * g);
        }
        components.push(ev.eval(&Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::PLUS, terms)],
        ))?);
    }
    Ok(Some(components))
}

/// The inverse transform of `big_f` from `s` into `x`; `None` when it
/// stays.
fn inverse(ev: &mut Eval<'_>, big_f: &Expr, s: &Expr, x: &Expr) -> Result<Option<Expr>, EvalError> {
    if *big_f == Expr::from(0) {
        return Ok(Some(Expr::from(0)));
    }
    let f = ev.eval(&Expr::call(
        Symbol::INVERSE_LAPLACE_TRANSFORM,
        [big_f.clone(), s.clone(), x.clone()],
    ))?;
    Ok((!f.contains(&Expr::from(Symbol::INVERSE_LAPLACE_TRANSFORM))).then_some(f))
}

#[cfg(test)]
mod tests;
