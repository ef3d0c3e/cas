//! Optimization (CAS-DOC-38): local extrema (§2, CAS-415), Xcas's
//! `extrema`; [`global`] optima (§3); linear programs ([`lpsolve`], by
//! [`simplex`] and, with integers, [`integer`], §4).
//!
//! `Extrema[f, vars]` and `Extrema[f, constraints, vars]`:
//!
//! - The variables are `x`, `x = a..b` or a list of them, distinct
//!   symbols; a range keeps the points within it, ends included.
//! - The critical points are the real solutions of `grad f = 0` by `Solve`
//!   (CAS-DOC-11), or with equality constraints `g_k = c_k` those of the
//!   Lagrangian `f − Σ λ_k (g_k − c_k)` in the variables and the `λ_k`.
//!   Points that are not tuples of real numbers are dropped; a system that
//!   `Solve` leaves leaves the call.
//! - In one variable a point is classified by the first derivative of order
//!   `k ≥ 2` that is not 0 there, up to [`MAX_ORDER`]: a minimum for `k`
//!   even and a positive value, a maximum for a negative one, neither for
//!   `k` odd. In several, by the Hessian's leading principal minors
//!   (Sylvester): all positive, a minimum; alternating from negative, a
//!   maximum; otherwise neither, a saddle or a singular Hessian.
//! - With `m` constraints in `n` variables, by the bordered Hessian `(0, J;
//!   Jᵀ, H_L)`, `J` the constraints' Jacobian and `H_L` the Lagrangian's
//!   Hessian: its leading principal minors of orders `2m + 1 … m + n` all
//!   of the sign of `(−1)^m` give a minimum, of the sign of `(−1)^(k − m)`
//!   for the order `k` a maximum.
//! - The answer is giac's `{{minima}, {maxima}}`: points as numbers in one
//!   variable, as lists in several. Signs are the assumptions' or the
//!   balls'; one not decided leaves the point out, as giac reports a
//!   saddle or an indeterminate point and leaves it out.

pub(super) mod global;
mod integer;
pub(super) mod lpsolve;
mod simplex;

use std::cmp::Ordering;

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The highest derivative tried in one variable.
const MAX_ORDER: i64 = 12;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn d(f: &Expr, x: &Expr) -> Expr {
    Expr::call(Symbol::D, [f.clone(), x.clone()])
}

/// A variable and its range, ends included; an infinite end is `None`.
#[derive(Clone, Debug)]
pub(super) struct Variable {
    pub(super) x: Expr,
    pub(super) low: Option<Expr>,
    pub(super) high: Option<Expr>,
}

/// Whether `e` is infinite: `Infinity`, `-Infinity`, `ComplexInfinity`.
fn infinite(e: &Expr) -> bool {
    e.may_contain(Symbol::INFINITY)
        || e.may_contain(Symbol::DIRECTED_INFINITY)
        || e.may_contain(Symbol::COMPLEX_INFINITY)
}

/// `x`, `x = a..b`, or a list of them: distinct symbols with their
/// ranges.
pub(super) fn variables(e: &Expr) -> Option<Vec<Variable>> {
    let symbol = |v: &Expr| v.as_symbol().is_some_and(|s| !s.is_builtin());
    let one = |v: &Expr| -> Option<Variable> {
        if symbol(v) {
            return Some(Variable {
                x: v.clone(),
                low: None,
                high: None,
            });
        }
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, range]) = (head(v), v.args()) else {
            return None;
        };
        let (Some(Symbol::RANGE), [a, b]) = (head(range), range.args()) else {
            return None;
        };
        let end = |c: &Expr| (!infinite(c)).then(|| c.clone());
        symbol(x).then(|| Variable {
            x: x.clone(),
            low: end(a),
            high: end(b),
        })
    };
    let items = match head(e) {
        Some(Symbol::LIST) => e.args().to_vec(),
        _ => vec![e.clone()],
    };
    let vars: Vec<Variable> = items.iter().map(one).collect::<Option<_>>()?;
    let distinct = vars
        .iter()
        .enumerate()
        .all(|(i, v)| vars[..i].iter().all(|w| w.x != v.x));
    (!vars.is_empty() && distinct).then_some(vars)
}

/// The equations of a constraint or a list of them, as `g − c`; `None`
/// for anything else.
pub(super) fn equalities(e: &Expr) -> Option<Vec<Expr>> {
    let one = |c: &Expr| match (head(c), c.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => Some(l.clone() - r.clone()),
        _ => None,
    };
    match head(e) {
        Some(Symbol::LIST) => e.args().iter().map(one).collect(),
        _ => one(e).map(|g| vec![g]),
    }
}

/// The sign of a real number, `None` when not decided.
pub(super) fn sign(ev: &Eval<'_>, e: &Expr) -> Option<Ordering> {
    if *e == Expr::from(0) {
        return Some(Ordering::Equal);
    }
    real_sign(ev.context(), e)
}

/// Whether `c`, a real number, lies within the variable's range.
pub(super) fn within(ev: &mut Eval<'_>, v: &Variable, c: &Expr) -> Result<bool, EvalError> {
    for (end, above) in [(&v.low, true), (&v.high, false)] {
        if let Some(end) = end {
            let gap = ev.eval(&if above {
                c.clone() - end.clone()
            } else {
                end.clone() - c.clone()
            })?;
            if !matches!(sign(ev, &gap), Some(Ordering::Greater | Ordering::Equal)) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Whether `c` is a real number: free of symbols but constants, and real.
fn real_number(ev: &Eval<'_>, c: &Expr) -> bool {
    c.symbols().iter().all(|s| s.is_builtin())
        && crate::assume::element(ev.context(), c, crate::assume::Set::Reals)
            == crate::assume::Truth::True
}

/// The real solutions of `equations = 0` in `unknowns`, as tuples; `None`
/// when `Solve` leaves the system.
pub(super) fn solutions(
    ev: &mut Eval<'_>,
    equations: Vec<Expr>,
    unknowns: &[Expr],
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let equations: Vec<Expr> = equations
        .into_iter()
        .map(|g| Expr::call(Symbol::EQUATION, [g, Expr::from(0)]))
        .collect();
    let solved = ev.eval(&Expr::call(
        Symbol::SOLVE,
        [list(equations), list(unknowns.to_vec())],
    ))?;
    if head(&solved) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let mut tuples = Vec::new();
    for t in solved.args() {
        let point = match head(t) {
            Some(Symbol::LIST) if t.args().len() == unknowns.len() => t.args().to_vec(),
            _ if unknowns.len() == 1 => vec![t.clone()],
            _ => return Ok(None),
        };
        if point.iter().all(|c| real_number(ev, c)) {
            tuples.push(point);
        }
    }
    Ok(Some(tuples))
}

/// What a critical point is.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Minimum,
    Maximum,
    Neither,
}

/// `e` evaluated, then with the variables at the point: a derivative is
/// taken before its variable is replaced.
fn at(ev: &mut Eval<'_>, e: &Expr, xs: &[Expr], point: &[Expr]) -> Result<Expr, EvalError> {
    let rules: Vec<(Expr, Expr)> = xs.iter().cloned().zip(point.iter().cloned()).collect();
    let e = ev.eval(e)?;
    ev.eval(&e.substitute(&rules))
}

/// A critical point of `f` in one variable, by its derivatives there.
fn in_one(ev: &mut Eval<'_>, f: &Expr, x: &Expr, c: &Expr) -> Result<Kind, EvalError> {
    for k in 2..=MAX_ORDER {
        let derivative = Expr::call(Symbol::D, [f.clone(), list(vec![x.clone(), Expr::from(k)])]);
        let value = at(
            ev,
            &derivative,
            std::slice::from_ref(x),
            std::slice::from_ref(c),
        )?;
        let even = k % 2 == 0;
        match sign(ev, &value) {
            Some(Ordering::Equal) => {}
            Some(Ordering::Greater) if even => return Ok(Kind::Minimum),
            Some(Ordering::Less) if even => return Ok(Kind::Maximum),
            _ => return Ok(Kind::Neither),
        }
    }
    Ok(Kind::Neither)
}

/// The signs of the leading principal minors of `m` from the order `from`
/// on; `None` when one is not decided.
fn minor_signs(
    ev: &mut Eval<'_>,
    m: &[Vec<Expr>],
    from: usize,
) -> Result<Option<Vec<Ordering>>, EvalError> {
    let mut signs = Vec::new();
    for k in from..=m.len() {
        let block = list(m[..k].iter().map(|row| list(row[..k].to_vec())).collect());
        let det = ev.eval(&Expr::call(Symbol::DET, [block]))?;
        let Some(s) = sign(ev, &det) else {
            return Ok(None);
        };
        signs.push(s);
    }
    Ok(Some(signs))
}

/// The matrix of second derivatives of `f` in `xs` at the point, with the
/// unknowns `all` at `point`, evaluated.
fn hessian(
    ev: &mut Eval<'_>,
    f: &Expr,
    xs: &[Expr],
    all: &[Expr],
    point: &[Expr],
) -> Result<Vec<Vec<Expr>>, EvalError> {
    let mut rows = Vec::new();
    for x in xs {
        let mut row = Vec::new();
        for y in xs {
            row.push(at(ev, &d(&d(f, x), y), all, point)?);
        }
        rows.push(row);
    }
    Ok(rows)
}

/// A critical point by the Hessian's minors, or with `m` constraints the
/// bordered Hessian's.
fn in_several(
    ev: &mut Eval<'_>,
    lagrangian: &Expr,
    constraints: &[Expr],
    xs: &[Expr],
    all: &[Expr],
    point: &[Expr],
) -> Result<Kind, EvalError> {
    let (n, m) = (xs.len(), constraints.len());
    // As many constraints as variables leave isolated points.
    if m >= n {
        return Ok(Kind::Neither);
    }
    let h = hessian(ev, lagrangian, xs, all, point)?;
    if m == 0 {
        let Some(signs) = minor_signs(ev, &h, 1)? else {
            return Ok(Kind::Neither);
        };
        let positive = signs.iter().all(|&s| s == Ordering::Greater);
        let alternating = signs.iter().enumerate().all(|(i, &s)| {
            s == if i % 2 == 0 {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        });
        return Ok(if positive {
            Kind::Minimum
        } else if alternating {
            Kind::Maximum
        } else {
            Kind::Neither
        });
    }
    // (0, J; Jᵀ, H_L).
    let mut bordered = vec![vec![Expr::from(0); m + n]; m + n];
    for (j, g) in constraints.iter().enumerate() {
        for (k, x) in xs.iter().enumerate() {
            let entry = at(ev, &d(g, x), all, point)?;
            bordered[j][m + k] = entry.clone();
            bordered[m + k][j] = entry;
        }
    }
    for (i, row) in h.into_iter().enumerate() {
        for (k, entry) in row.into_iter().enumerate() {
            bordered[m + i][m + k] = entry;
        }
    }
    let Some(signs) = minor_signs(ev, &bordered, 2 * m + 1)? else {
        return Ok(Kind::Neither);
    };
    let parity = |p: usize| {
        if p % 2 == 0 {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    };
    let minimum = signs.iter().all(|&s| s == parity(m));
    let maximum = signs
        .iter()
        .enumerate()
        .all(|(i, &s)| s == parity(2 * m + 1 + i - m));
    Ok(if minimum {
        Kind::Minimum
    } else if maximum {
        Kind::Maximum
    } else {
        Kind::Neither
    })
}

/// `Extrema[f, vars]` and `Extrema[f, constraints, vars]`.
pub(super) fn extrema(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, constraints, vars) = match e.args() {
        [f, vars] => (f, Vec::new(), vars),
        [f, constraints, vars] => {
            let Some(gs) = equalities(constraints) else {
                return Ok(None);
            };
            (f, gs, vars)
        }
        _ => return Ok(None),
    };
    let Some(vars) = variables(vars) else {
        return Ok(None);
    };
    let xs: Vec<Expr> = vars.iter().map(|v| v.x.clone()).collect();
    let lambdas: Vec<Expr> = (1..=constraints.len())
        .map(|k| Expr::symbol(&format!("λ#{k}")))
        .collect();
    let mut lagrangian = f.clone();
    for (g, l) in constraints.iter().zip(&lambdas) {
        lagrangian -= l.clone() * g.clone();
    }
    let all: Vec<Expr> = xs.iter().chain(&lambdas).cloned().collect();
    let gradient: Vec<Expr> = all.iter().map(|u| d(&lagrangian, u)).collect();
    let Some(points) = solutions(ev, gradient, &all)? else {
        return Ok(None);
    };
    let (mut minima, mut maxima) = (Vec::new(), Vec::new());
    'points: for point in points {
        for (v, c) in vars.iter().zip(&point) {
            if !within(ev, v, c)? {
                continue 'points;
            }
        }
        let kind = if xs.len() == 1 && constraints.is_empty() {
            in_one(ev, f, &xs[0], &point[0])?
        } else {
            in_several(ev, &lagrangian, &constraints, &xs, &all, &point)?
        };
        let shown = if xs.len() == 1 {
            point[0].clone()
        } else {
            list(point[..xs.len()].to_vec())
        };
        match kind {
            Kind::Minimum => minima.push(shown),
            Kind::Maximum => maxima.push(shown),
            Kind::Neither => {}
        }
    }
    let result = list(vec![list(minima), list(maxima)]);
    ev.step("extrema", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
