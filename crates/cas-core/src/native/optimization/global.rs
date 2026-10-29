//! Global optima (CAS-DOC-38 §3, CAS-416): `Minimize` and `Maximize`,
//! Xcas's `minimize` and `maximize`, and `FMin` and `FMax`, its `fMin` and
//! `fMax`.
//!
//! `Minimize[f, vars]`, `Minimize[f, constraints, vars]`, and the same with a
//! last argument `point` or `coordinates`:
//!
//! - The constraints are equations `g = c`, and inequalities `g ≤ c`, `g ≥
//!   c`, `<` and `>` taken as their closures; the variables' ranges are
//!   inequalities too. More than [`MAX_INEQUALITIES`] leave the call.
//! - The candidates, by Karush, Kuhn and Tucker: for each set of
//!   inequalities taken as equations, with the equations, at most as many
//!   as the variables, the critical points of the Lagrangian (`Extrema`'s
//!   `solutions`) that satisfy every inequality. A bound `x = c` among them
//!   fixes `x` instead of taking a multiplier, so an interval's ends need
//!   no system. That takes in the interior's critical points, those of
//!   each face, and the corners.
//! - The region must be bounded: every variable in a finite range, or a
//!   constraint whose function is quadratic with a definite Hessian (a
//!   disc, a ball, an ellipse); or the function quadratic with a definite
//!   Hessian of the optimum's side, whose level sets are bounded, as
//!   `x² + y²` on `x + y = 1`. In one variable without constraints, an
//!   infinite end is compared by its limit (CAS-DOC-10): `-Infinity` there
//!   is the answer, and a finite limit below the candidates' values is an
//!   infimum not attained, with a warning. Otherwise the call stays, with
//!   a warning: giac's `minimize(x^3,x)` is 0, its critical value.
//! - The least value, or with the option `{value, {points}}`, each point a
//!   list, giac's `[-1/4,[[sqrt(2)/2],[-sqrt(2)/2]]]`. Values are compared
//!   exactly, `N`'s balls telling them apart, `Simplify` a difference that
//!   they do not; one that neither tells leaves the call.
//! - `FMin` and `FMax` are the first point's abscissa, in one variable.

use std::cmp::Ordering;

use super::{Variable, sign, solutions, variables};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The most inequalities, the ranges' ends among them.
const MAX_INEQUALITIES: usize = 8;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn d(f: &Expr, x: &Expr) -> Expr {
    Expr::call(Symbol::D, [f.clone(), x.clone()])
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// A problem: the function, the variables, `h = 0` and `h ≤ 0`, whether
/// the points are asked for.
struct Problem {
    f: Expr,
    vars: Vec<Variable>,
    equalities: Vec<Expr>,
    inequalities: Vec<Expr>,
    point: bool,
}

/// The constraints of a relation or a list of them.
fn relations(e: &Expr, equalities: &mut Vec<Expr>, inequalities: &mut Vec<Expr>) -> Option<()> {
    let items = match head(e) {
        Some(Symbol::LIST) => e.args().to_vec(),
        _ => vec![e.clone()],
    };
    for c in items {
        let (Some(h), [l, r]) = (head(&c), c.args()) else {
            return None;
        };
        match h {
            Symbol::EQUATION | Symbol::EQUAL => equalities.push(l.clone() - r.clone()),
            Symbol::LESS_EQUAL | Symbol::LESS => inequalities.push(l.clone() - r.clone()),
            Symbol::GREATER_EQUAL | Symbol::GREATER => {
                inequalities.push(r.clone() - l.clone());
            }
            _ => return None,
        }
    }
    Some(())
}

/// The problem of a call.
fn problem(e: &Expr) -> Option<Problem> {
    let mut args = e.args();
    let mut point = false;
    if let [rest @ .., last] = args
        && (*last == Expr::symbol("point") || *last == Expr::symbol("coordinates"))
    {
        point = true;
        args = rest;
    }
    let (f, constraints, vars) = match args {
        [f, vars] => (f, None, vars),
        [f, constraints, vars] => (f, Some(constraints), vars),
        _ => return None,
    };
    let vars = variables(vars)?;
    let (mut equalities, mut inequalities) = (Vec::new(), Vec::new());
    if let Some(c) = constraints {
        relations(c, &mut equalities, &mut inequalities)?;
    }
    for v in &vars {
        if let Some(low) = &v.low {
            inequalities.push(low.clone() - v.x.clone());
        }
        if let Some(high) = &v.high {
            inequalities.push(v.x.clone() - high.clone());
        }
    }
    (inequalities.len() <= MAX_INEQUALITIES).then(|| Problem {
        f: f.clone(),
        vars,
        equalities,
        inequalities,
        point,
    })
}

/// Whether the constraint `g` (`g = 0`, or `g ≤ 0` when `inequality`)
/// bounds every variable: `g` quadratic with a definite Hessian, positive
/// for an inequality.
fn bounds_all(
    ev: &mut Eval<'_>,
    g: &Expr,
    xs: &[Expr],
    inequality: bool,
) -> Result<bool, EvalError> {
    let mut h = Vec::new();
    for x in xs {
        let mut row = Vec::new();
        for y in xs {
            let entry = ev.eval(&d(&d(g, x), y))?;
            if xs.iter().any(|v| entry.contains(v)) {
                return Ok(false);
            }
            row.push(entry);
        }
        h.push(row);
    }
    let Some(signs) = super::minor_signs(ev, &h, 1)? else {
        return Ok(false);
    };
    let positive = signs.iter().all(|&s| s == Ordering::Greater);
    let negative = signs.iter().enumerate().all(|(i, &s)| {
        s == if i % 2 == 0 {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    });
    Ok(positive || negative && !inequality)
}

/// Whether the optimum is shown to lie among the candidates: the region
/// bounded, or the function quadratic with a definite Hessian of the
/// optimum's side, which its level sets bound.
fn bounded(ev: &mut Eval<'_>, p: &Problem, maximum: bool) -> Result<bool, EvalError> {
    if p.vars.iter().all(|v| v.low.is_some() && v.high.is_some()) {
        return Ok(true);
    }
    let xs: Vec<Expr> = p.vars.iter().map(|v| v.x.clone()).collect();
    let side = if maximum { -p.f.clone() } else { p.f.clone() };
    if bounds_all(ev, &side, &xs, true)? {
        return Ok(true);
    }
    for g in &p.equalities {
        if bounds_all(ev, g, &xs, false)? {
            return Ok(true);
        }
    }
    for g in &p.inequalities {
        if bounds_all(ev, g, &xs, true)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `h`, evaluated, as a bound `s (x_i − c)`, `s = ±1` and `c` free of the
/// variables: `(i, c)`.
fn bound(ev: &mut Eval<'_>, h: &Expr, xs: &[Expr]) -> Result<Option<(usize, Expr)>, EvalError> {
    let mut slopes = Vec::with_capacity(xs.len());
    for x in xs {
        slopes.push(ev.eval(&d(h, x))?);
    }
    let one = Expr::from(1);
    let minus_one = Expr::from(-1);
    let zero = Expr::from(0);
    let Some(i) = slopes.iter().position(|s| *s == one || *s == minus_one) else {
        return Ok(None);
    };
    if slopes.iter().enumerate().any(|(j, s)| j != i && *s != zero) {
        return Ok(None);
    }
    // c = x_i − h/s.
    let c = ev.eval(&(xs[i].clone() - h.clone() / slopes[i].clone()))?;
    Ok(xs.iter().all(|x| !c.contains(x)).then_some((i, c)))
}

/// The candidates, as tuples of the variables' values; `None` when a
/// system stays. The active bounds fix their variables, the other active
/// constraints get multipliers: giac's `maximize(sin(x),x=0..4)` needs no
/// system in `x` and a multiplier.
fn candidates(ev: &mut Eval<'_>, p: &Problem) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let xs: Vec<Expr> = p.vars.iter().map(|v| v.x.clone()).collect();
    let n = xs.len();
    let mut found: Vec<Vec<Expr>> = Vec::new();
    'masks: for mask in 0u32..(1 << p.inequalities.len()) {
        let active = p
            .inequalities
            .iter()
            .enumerate()
            .filter(|(k, _)| mask >> k & 1 == 1)
            .map(|(_, h)| h);
        let equations: Vec<&Expr> = p.equalities.iter().chain(active).collect();
        if equations.len() > n {
            continue;
        }
        let mut fixed: Vec<Option<Expr>> = vec![None; n];
        let mut rest = Vec::new();
        for g in equations {
            match bound(ev, g, &xs)? {
                Some((i, c)) => match &fixed[i] {
                    None => fixed[i] = Some(c),
                    Some(other) if *other == c => {}
                    Some(_) => continue 'masks,
                },
                None => rest.push(g.clone()),
            }
        }
        let rules: Vec<(Expr, Expr)> = xs
            .iter()
            .zip(&fixed)
            .filter_map(|(x, c)| c.as_ref().map(|c| (x.clone(), c.clone())))
            .collect();
        let free: Vec<Expr> = xs
            .iter()
            .zip(&fixed)
            .filter(|(_, c)| c.is_none())
            .map(|(x, _)| x.clone())
            .collect();
        let rest: Vec<Expr> = rest.iter().map(|g| g.substitute(&rules)).collect();
        let points = if free.is_empty() {
            vec![Vec::new()]
        } else {
            if rest.len() > free.len() {
                continue;
            }
            let lambdas: Vec<Expr> = (1..=rest.len())
                .map(|k| Expr::symbol(&format!("λ#{k}")))
                .collect();
            let mut lagrangian = p.f.substitute(&rules);
            for (g, l) in rest.iter().zip(&lambdas) {
                lagrangian -= l.clone() * g.clone();
            }
            let unknowns: Vec<Expr> = free.iter().chain(&lambdas).cloned().collect();
            let gradient: Vec<Expr> = unknowns.iter().map(|u| d(&lagrangian, u)).collect();
            let Some(points) = solutions(ev, gradient, &unknowns)? else {
                return Ok(None);
            };
            points
        };
        'points: for point in points {
            let mut values = point.into_iter();
            let point: Vec<Expr> = fixed
                .iter()
                .map(|c| c.clone().or_else(|| values.next()).expect("a value"))
                .collect();
            let rules: Vec<(Expr, Expr)> = xs.iter().cloned().zip(point.iter().cloned()).collect();
            // Bounds alone fix the point: the equations are not solved.
            let unsolved = if free.is_empty() {
                &p.equalities[..]
            } else {
                &[]
            };
            for g in unsolved {
                let value = ev.eval(&g.substitute(&rules))?;
                if sign(ev, &value) != Some(Ordering::Equal) {
                    continue 'points;
                }
            }
            for h in &p.inequalities {
                let value = ev.eval(&h.substitute(&rules))?;
                if !matches!(sign(ev, &value), Some(Ordering::Less | Ordering::Equal)) {
                    continue 'points;
                }
            }
            if !found.contains(&point) {
                found.push(point);
            }
        }
    }
    Ok(Some(found))
}

/// Whether `a` is better than `b` (`Less`), as good (`Equal`) or worse,
/// for a minimum or a maximum; `None` when not told.
fn compare(
    ev: &mut Eval<'_>,
    a: &Expr,
    b: &Expr,
    maximum: bool,
) -> Result<Option<Ordering>, EvalError> {
    let difference = ev.eval(&(a.clone() - b.clone()))?;
    let order = if let Some(o) = sign(ev, &difference) {
        Some(o)
    } else {
        let simplified = ev.eval(&Expr::call(Symbol::SIMPLIFY, [difference]))?;
        sign(ev, &simplified)
    };
    Ok(order.map(|o| if maximum { o.reverse() } else { o }))
}

/// An optimum: its value and the points where it is reached, none for an
/// infimum not attained.
type Optimum = (Expr, Vec<Vec<Expr>>);

/// The optimum of a call, or `None` to leave it.
fn optimum(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Problem,
    maximum: bool,
) -> Result<Option<Optimum>, EvalError> {
    let one = p.vars.len() == 1
        && p.equalities.is_empty()
        && p.inequalities.len()
            == usize::from(p.vars[0].low.is_some()) + usize::from(p.vars[0].high.is_some());
    let finite = bounded(ev, p, maximum)?;
    if !finite && !one {
        warn(ev, e, "the region is not shown to be bounded");
        return Ok(None);
    }
    let Some(points) = candidates(ev, p)? else {
        return Ok(None);
    };
    let xs: Vec<Expr> = p.vars.iter().map(|v| v.x.clone()).collect();
    let mut best: Option<Optimum> = None;
    for point in points {
        let rules: Vec<(Expr, Expr)> = xs.iter().cloned().zip(point.iter().cloned()).collect();
        let value = ev.eval(&p.f.substitute(&rules))?;
        best = match best {
            None => Some((value, vec![point])),
            Some((b, mut at)) => match compare(ev, &value, &b, maximum)? {
                Some(Ordering::Less) => Some((value, vec![point])),
                Some(Ordering::Equal) => {
                    at.push(point);
                    Some((b, at))
                }
                Some(Ordering::Greater) => Some((b, at)),
                None => {
                    warn(ev, e, "values not told apart");
                    return Ok(None);
                }
            },
        };
    }
    if finite {
        return Ok(best);
    }
    // One variable with an infinite end: its limit.
    let v = &p.vars[0];
    let infinity = Expr::from(Symbol::INFINITY);
    for (end, infinite) in [(&v.low, -infinity.clone()), (&v.high, infinity.clone())] {
        if end.is_some() {
            continue;
        }
        let limit = ev.eval(&Expr::call(
            Symbol::LIMIT,
            [p.f.clone(), v.x.clone(), infinite],
        ))?;
        let plus = ev.eval(&infinity)?;
        let minus = ev.eval(&-infinity.clone())?;
        let (worst, unbounded) = if maximum {
            (minus, plus)
        } else {
            (plus, minus)
        };
        if limit == worst {
            continue;
        }
        if limit == unbounded {
            return Ok(Some((limit, Vec::new())));
        }
        if !super::real_number(ev, &limit) {
            warn(ev, e, "no limit at an infinite end");
            return Ok(None);
        }
        let below = match &best {
            None => true,
            Some((b, _)) => match compare(ev, &limit, b, maximum)? {
                Some(Ordering::Less) => true,
                Some(_) => false,
                None => {
                    warn(ev, e, "values not told apart");
                    return Ok(None);
                }
            },
        };
        if below {
            warn(ev, e, "an infimum or a supremum not attained");
            best = Some((limit, Vec::new()));
        }
    }
    Ok(best)
}

/// `Minimize` and `Maximize`.
fn extremum(ev: &mut Eval<'_>, e: &Expr, maximum: bool) -> Native {
    let Some(p) = problem(e) else {
        return Ok(None);
    };
    let Some((value, points)) = optimum(ev, e, &p, maximum)? else {
        return Ok(None);
    };
    let result = if p.point {
        list(vec![value, list(points.into_iter().map(list).collect())])
    } else {
        value
    };
    ev.step(if maximum { "maximize" } else { "minimize" }, e, &result);
    Ok(Some(result))
}

/// `Minimize[…]`.
pub(in crate::native) fn minimize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    extremum(ev, e, false)
}

/// `Maximize[…]`.
pub(in crate::native) fn maximize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    extremum(ev, e, true)
}

/// `FMin` and `FMax`: the first point's abscissa, in one variable.
fn abscissa(ev: &mut Eval<'_>, e: &Expr, maximum: bool) -> Native {
    let Some(p) = problem(e) else {
        return Ok(None);
    };
    if p.vars.len() != 1 {
        return Ok(None);
    }
    let Some((_, points)) = optimum(ev, e, &p, maximum)? else {
        return Ok(None);
    };
    let Some(first) = points.first() else {
        return Ok(None);
    };
    let result = first[0].clone();
    ev.step(if maximum { "fMax" } else { "fMin" }, e, &result);
    Ok(Some(result))
}

/// `FMin[…]`.
pub(in crate::native) fn f_min(ev: &mut Eval<'_>, e: &Expr) -> Native {
    abscissa(ev, e, false)
}

/// `FMax[…]`.
pub(in crate::native) fn f_max(ev: &mut Eval<'_>, e: &Expr) -> Native {
    abscissa(ev, e, true)
}

#[cfg(test)]
mod tests;
