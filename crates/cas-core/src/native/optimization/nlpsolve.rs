//! Nonlinear programs (CAS-519): `NLPSolve`, Xcas's `nlpsolve`, a local
//! optimum of a smooth function under equations, inequalities and bounds.
//!
//! `NLPSolve[f, constraints, bounds…, options…]`:
//!
//! - The constraints are `g = c`, `g ≤ c`, `g ≥ c`, and `<`, `>` taken as
//!   their closures, in a list that may hold bounds too; a bound is `x =
//!   a..b`, ends real or infinite. The options are `maximize`, `maximize =
//!   true` or `false`, `nlp_nonnegative = true` (every variable at least 0
//!   unless bounded) and `nlp_initialpoint = {x = v, …}`. The variables are
//!   the symbols of `f`, the constraints and the bounds, sorted, giac's
//!   order.
//! - The functions must compile to ball operations, as `NDSolve`'s right
//!   sides (CAS-DOC-35 §3): no parameters, no `Abs` or steps.
//! - The search, on doubles: an augmented Lagrangian (Powell, Hestenes,
//!   Rockafellar) whose subproblems a spectral projected gradient (Birgin,
//!   Martínez, Raydan) solves within the bounds, from the initial point, or
//!   each variable's range midpoint, an end ± 1, or 0. Like giac's, it is
//!   local.
//! - Its point is then certified: the bounds and inequalities it meets are
//!   taken as equations, the bounds fixing their variables at the exact
//!   ends, and the Karush–Kuhn–Tucker system (the Lagrangian's gradient in
//!   the free variables, and the active constraints) is solved by Newton's
//!   iteration and Krawczyk's test, as `FSolve`'s systems (CAS-518). An
//!   inequality's multiplier that is certainly negative, a fixed
//!   variable's derivative of the wrong sign, or an inactive inequality
//!   certainly violated, means that the point is no optimum: the call
//!   stays, with a warning, as it does for a system that is not certified.
//! - The answer is giac's `{value, {x = …, …}}`, balls or exact ends, the
//!   value `f`'s at the point. A search that runs off to infinity or finds
//!   no feasible point leaves the call with a warning.

#![allow(
    clippy::many_single_char_names,
    reason = "the method's names: x, f, g, h, multipliers, steps s and y"
)]

use cas_num::Rational;
use cas_num::ball::{Ball, Float};

use super::super::approximation;
use super::super::ndsolve::rhs::Node;
use super::super::solve::fsolve::systems::certified_root;
use super::lpsolve::symbols_in_order;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The most variables.
const MAX_VARIABLES: usize = 16;

/// The most constraints.
const MAX_CONSTRAINTS: usize = 32;

/// The precision of the search's ball evaluations.
const SEARCH_BITS: u64 = 64;

/// The most outer (multiplier) iterations.
const MAX_OUTER: usize = 60;

/// The most projected gradient steps per subproblem.
const MAX_INNER: usize = 3000;

/// The largest penalty.
const MAX_PENALTY: f64 = 1e12;

/// A coordinate this large means the search ran off.
const RUNAWAY: f64 = 1e12;

/// Within this (relative) distance a bound or inequality is active.
const ACTIVE: f64 = 1e-6;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// A bound's end: exact, with its double.
#[derive(Clone)]
struct End {
    exact: Expr,
    value: f64,
}

/// A program as read: the objective (negated for a maximum), the
/// equations `h = 0` and inequalities `g ≤ 0`, the bounds and the start.
struct Program {
    xs: Vec<Expr>,
    objective: Expr,
    minimized: Expr,
    equalities: Vec<Expr>,
    inequalities: Vec<Expr>,
    lower: Vec<Option<End>>,
    upper: Vec<Option<End>>,
    start: Vec<Option<Expr>>,
}

/// `NLPSolve[f, constraints, bounds…, options…]`.
pub(in crate::native) fn nlp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(program) = program(ev, e)? else {
        return Ok(None);
    };
    let Some(compiled) = compile(ev, &program)? else {
        return Ok(None);
    };
    let Some(x0) = start(ev.context(), &program) else {
        return Ok(None);
    };
    let lo: Vec<f64> = program
        .lower
        .iter()
        .map(|b| b.as_ref().map_or(f64::NEG_INFINITY, |b| b.value))
        .collect();
    let hi: Vec<f64> = program
        .upper
        .iter()
        .map(|b| b.as_ref().map_or(f64::INFINITY, |b| b.value))
        .collect();
    if lo.iter().zip(&hi).any(|(a, b)| a > b) {
        warn(ev, e, "a variable's range is empty");
        return Ok(None);
    }
    let search = match lagrangian(&compiled, x0, &lo, &hi) {
        Ok(found) => found,
        Err(Failure::Runaway) => {
            warn(
                ev,
                e,
                "the search ran off to infinity: the program may be unbounded",
            );
            return Ok(None);
        }
        Err(Failure::Infeasible) => {
            warn(ev, e, "the search found no feasible point");
            return Ok(None);
        }
        Err(Failure::Domain) => {
            warn(ev, e, "the search left the functions' domain");
            return Ok(None);
        }
    };
    let Some(answer) = certify(ev, &program, &compiled, &search, &lo, &hi)? else {
        warn(
            ev,
            e,
            "the search's point could not be certified as an optimum",
        );
        return Ok(None);
    };
    ev.step("NLPSolve", e, &answer);
    Ok(Some(answer))
}

/// The arguments after the objective, sorted out.
#[derive(Default)]
struct Items<'a> {
    maximum: bool,
    nonnegative: bool,
    initial: &'a [Expr],
    relations: Vec<(Symbol, &'a Expr, &'a Expr)>,
    ranges: Vec<(&'a Expr, &'a Expr, &'a Expr)>,
}

/// A `true` or `false` option's value.
fn flag(r: &Expr) -> Option<bool> {
    match r.as_symbol() {
        Some(Symbol::TRUE) => Some(true),
        Some(Symbol::FALSE) => Some(false),
        _ => None,
    }
}

/// The constraints, bounds and options; `None` for anything else.
fn items(rest: &[Expr]) -> Option<Items<'_>> {
    let mut all: Vec<&Expr> = Vec::new();
    for a in rest {
        if head(a) == Some(Symbol::LIST) {
            all.extend(a.args());
        } else {
            all.push(a);
        }
    }
    let mut found = Items::default();
    for item in all {
        if item.as_symbol() == Some(Symbol::MAXIMIZE) {
            found.maximum = true;
            continue;
        }
        let (Some(h), [l, r]) = (head(item), item.args()) else {
            return None;
        };
        match h {
            Symbol::EQUATION | Symbol::EQUAL => {
                if l.as_symbol() == Some(Symbol::MAXIMIZE) {
                    found.maximum = flag(r)?;
                } else if *l == Expr::symbol("nlp_nonnegative") {
                    found.nonnegative = flag(r)?;
                } else if *l == Expr::symbol("nlp_initialpoint") {
                    if head(r) != Some(Symbol::LIST) {
                        return None;
                    }
                    found.initial = r.args();
                } else if let (Some(Symbol::RANGE), [a, b]) = (head(r), r.args()) {
                    if l.as_symbol().is_none_or(Symbol::is_builtin) {
                        return None;
                    }
                    found.ranges.push((l, a, b));
                } else {
                    found.relations.push((Symbol::EQUATION, l, r));
                }
            }
            Symbol::LESS | Symbol::LESS_EQUAL | Symbol::GREATER | Symbol::GREATER_EQUAL => {
                found.relations.push((h, l, r));
            }
            _ => return None,
        }
    }
    Some(found)
}

/// The program of a call; `None` when it is not one.
fn program(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Program>, EvalError> {
    let [objective, rest @ ..] = e.args() else {
        return Ok(None);
    };
    let Some(Items {
        maximum,
        nonnegative,
        initial,
        relations,
        ranges,
    }) = items(rest)
    else {
        return Ok(None);
    };
    let mut xs = Vec::new();
    symbols_in_order(objective, &mut xs);
    for (_, l, r) in &relations {
        symbols_in_order(l, &mut xs);
        symbols_in_order(r, &mut xs);
    }
    for (x, _, _) in &ranges {
        symbols_in_order(x, &mut xs);
    }
    xs.sort();
    if xs.is_empty() || xs.len() > MAX_VARIABLES || relations.len() > MAX_CONSTRAINTS {
        return Ok(None);
    }
    let mut equalities = Vec::new();
    let mut inequalities = Vec::new();
    for (h, l, r) in relations {
        let (a, b) = match h {
            Symbol::GREATER | Symbol::GREATER_EQUAL => (r, l),
            _ => (l, r),
        };
        let g = ev.eval(&(a.clone() - b.clone()))?;
        if h == Symbol::EQUATION {
            equalities.push(g);
        } else {
            inequalities.push(g);
        }
    }
    let n = xs.len();
    let zero = End {
        exact: Expr::from(0),
        value: 0.0,
    };
    let mut lower: Vec<Option<End>> = vec![nonnegative.then(|| zero.clone()); n];
    let mut upper: Vec<Option<End>> = vec![None; n];
    for (x, a, b) in ranges {
        let j = xs.iter().position(|v| v == x).expect("collected");
        let (Ok(low), Ok(high)) = (end(ev.context(), a), end(ev.context(), b)) else {
            return Ok(None);
        };
        lower[j] = low;
        upper[j] = high;
    }
    let mut start = vec![None; n];
    for p in initial {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, v]) = (head(p), p.args()) else {
            return Ok(None);
        };
        let Some(j) = xs.iter().position(|y| y == x) else {
            return Ok(None);
        };
        start[j] = Some(v.clone());
    }
    let objective = ev.eval(objective)?;
    let minimized = if maximum {
        ev.eval(&(Expr::from(0) - objective.clone()))?
    } else {
        objective.clone()
    };
    Ok(Some(Program {
        xs,
        objective,
        minimized,
        equalities,
        inequalities,
        lower,
        upper,
        start,
    }))
}

/// A bound's end: `Ok(None)` for an infinite one, `Err` for one that is no
/// real number.
fn end(ctx: &Context, a: &Expr) -> Result<Option<End>, ()> {
    if head(a) == Some(Symbol::DIRECTED_INFINITY) {
        return Ok(None);
    }
    let value = double(ctx, a).ok_or(())?;
    Ok(Some(End {
        exact: a.clone(),
        value,
    }))
}

/// A real number's double, by its ball.
fn double(ctx: &Context, a: &Expr) -> Option<f64> {
    match approximation(ctx, a, SEARCH_BITS).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(b.mid().to_f64()).filter(|v| v.is_finite()),
        _ => None,
    }
}

/// A function and its gradient, compiled.
struct Smooth {
    value: Node,
    gradient: Vec<Node>,
}

/// The program's functions compiled: the minimized objective, the
/// equations and the inequalities.
struct Compiled {
    objective: Smooth,
    equalities: Vec<Smooth>,
    inequalities: Vec<Smooth>,
}

fn compile(ev: &mut Eval<'_>, p: &Program) -> Result<Option<Compiled>, EvalError> {
    let smooth = |ev: &mut Eval<'_>, f: &Expr| -> Result<Option<Smooth>, EvalError> {
        let Some(value) = Node::compile(ev.context(), f, &p.xs, SEARCH_BITS) else {
            return Ok(None);
        };
        let mut gradient = Vec::with_capacity(p.xs.len());
        for x in &p.xs {
            let d = ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?;
            let Some(node) = Node::compile(ev.context(), &d, &p.xs, SEARCH_BITS) else {
                return Ok(None);
            };
            gradient.push(node);
        }
        Ok(Some(Smooth { value, gradient }))
    };
    let Some(objective) = smooth(ev, &p.minimized)? else {
        return Ok(None);
    };
    let mut equalities = Vec::with_capacity(p.equalities.len());
    for h in &p.equalities {
        let Some(s) = smooth(ev, h)? else {
            return Ok(None);
        };
        equalities.push(s);
    }
    let mut inequalities = Vec::with_capacity(p.inequalities.len());
    for g in &p.inequalities {
        let Some(s) = smooth(ev, g)? else {
            return Ok(None);
        };
        inequalities.push(s);
    }
    Ok(Some(Compiled {
        objective,
        equalities,
        inequalities,
    }))
}

/// The start: the initial point's values, else each range's midpoint, an
/// end moved 1 inside, or 0; projected into the bounds.
fn start(ctx: &Context, p: &Program) -> Option<Vec<f64>> {
    let mut x = Vec::with_capacity(p.xs.len());
    for j in 0..p.xs.len() {
        let (lo, hi) = (p.lower[j].as_ref(), p.upper[j].as_ref());
        let v = match (&p.start[j], lo, hi) {
            (Some(v), _, _) => double(ctx, v)?,
            (None, Some(a), Some(b)) => f64::midpoint(a.value, b.value),
            (None, Some(a), None) => a.value + 1.0,
            (None, None, Some(b)) => b.value - 1.0,
            (None, None, None) => 0.0,
        };
        let v = lo.map_or(v, |a| v.max(a.value));
        x.push(hi.map_or(v, |b| v.min(b.value)));
    }
    Some(x)
}

/// The balls of a point's coordinates.
fn balls(x: &[f64]) -> Option<Vec<Ball>> {
    x.iter()
        .map(|&v| Float::from_f64(v).map(Ball::exact))
        .collect()
}

/// A node's value at `at`, as a finite double.
fn at(node: &Node, at: &[Ball]) -> Option<f64> {
    let b = node.ball(at, SEARCH_BITS)?;
    Some(b.mid().to_f64()).filter(|v| b.is_finite() && v.is_finite())
}

/// A function's value and gradient at `at`.
fn evaluate(s: &Smooth, point: &[Ball]) -> Option<(f64, Vec<f64>)> {
    let v = at(&s.value, point)?;
    let g = s
        .gradient
        .iter()
        .map(|d| at(d, point))
        .collect::<Option<Vec<_>>>()?;
    Some((v, g))
}

/// The multipliers and penalty of the augmented Lagrangian.
struct Multipliers {
    lambda: Vec<f64>,
    mu: Vec<f64>,
    rho: f64,
}

/// The augmented Lagrangian at a point: its value and gradient, and the
/// constraints' values.
struct Merit {
    value: f64,
    gradient: Vec<f64>,
    equalities: Vec<f64>,
    inequalities: Vec<f64>,
}

/// The augmented Lagrangian at `x`.
fn merit(c: &Compiled, m: &Multipliers, x: &[f64]) -> Option<Merit> {
    let point = balls(x)?;
    let (mut v, mut grad) = evaluate(&c.objective, &point)?;
    let mut hs = Vec::with_capacity(c.equalities.len());
    for (s, lambda) in c.equalities.iter().zip(&m.lambda) {
        let (h, dh) = evaluate(s, &point)?;
        let w = lambda + m.rho * h;
        v += lambda * h + m.rho / 2.0 * h * h;
        for (g, d) in grad.iter_mut().zip(&dh) {
            *g += w * d;
        }
        hs.push(h);
    }
    let mut gs = Vec::with_capacity(c.inequalities.len());
    for (s, mu) in c.inequalities.iter().zip(&m.mu) {
        let (g, dg) = evaluate(s, &point)?;
        let t = mu + m.rho * g;
        if t > 0.0 {
            v += (t * t - mu * mu) / (2.0 * m.rho);
            for (gr, d) in grad.iter_mut().zip(&dg) {
                *gr += t * d;
            }
        } else {
            v -= mu * mu / (2.0 * m.rho);
        }
        gs.push(g);
    }
    v.is_finite().then_some(Merit {
        value: v,
        gradient: grad,
        equalities: hs,
        inequalities: gs,
    })
}

fn project(x: &mut [f64], lo: &[f64], hi: &[f64]) {
    for ((v, a), b) in x.iter_mut().zip(lo).zip(hi) {
        *v = v.clamp(*a, *b);
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Why a search failed.
enum Failure {
    Runaway,
    Infeasible,
    Domain,
}

/// What the search found: the point and the multipliers.
struct Search {
    x: Vec<f64>,
    lambda: Vec<f64>,
    mu: Vec<f64>,
}

/// The spectral projected gradient on the augmented Lagrangian from `x`:
/// the subproblem's minimum within the bounds, to `tol`.
fn subproblem(
    c: &Compiled,
    m: &Multipliers,
    mut x: Vec<f64>,
    lo: &[f64],
    hi: &[f64],
    tol: f64,
) -> Result<Vec<f64>, Failure> {
    project(&mut x, lo, hi);
    let first = merit(c, m, &x).ok_or(Failure::Domain)?;
    let (mut v, mut g) = (first.value, first.gradient);
    let mut alpha = 1.0 / g.iter().fold(1.0_f64, |a, d| a.max(d.abs()));
    for _ in 0..MAX_INNER {
        let mut p: Vec<f64> = x.iter().zip(&g).map(|(a, d)| a - d).collect();
        project(&mut p, lo, hi);
        let stationary = p
            .iter()
            .zip(&x)
            .fold(0.0_f64, |a, (p, x)| a.max((p - x).abs()));
        if stationary <= tol {
            return Ok(x);
        }
        let mut d: Vec<f64> = x.iter().zip(&g).map(|(a, gr)| a - alpha * gr).collect();
        project(&mut d, lo, hi);
        for (di, xi) in d.iter_mut().zip(&x) {
            *di -= xi;
        }
        let slope = dot(&g, &d);
        let mut t = 1.0;
        let mut next = None;
        for _ in 0..60 {
            let xn: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + t * b).collect();
            if let Some(found) = merit(c, m, &xn)
                && found.value <= v + 1e-4 * t * slope
            {
                next = Some((xn, found.value, found.gradient));
                break;
            }
            t /= 2.0;
        }
        let Some((xn, vn, gn)) = next else {
            return Ok(x);
        };
        let s: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        alpha = if sy > 0.0 {
            (dot(&s, &s) / sy).clamp(1e-12, 1e12)
        } else {
            1e12_f64.min(alpha * 10.0)
        };
        if xn.iter().any(|v| v.abs() > RUNAWAY) {
            return Err(Failure::Runaway);
        }
        (x, v, g) = (xn, vn, gn);
    }
    Ok(x)
}

/// The augmented Lagrangian's outer loop from `x0`.
fn lagrangian(c: &Compiled, x0: Vec<f64>, lo: &[f64], hi: &[f64]) -> Result<Search, Failure> {
    let mut m = Multipliers {
        lambda: vec![0.0; c.equalities.len()],
        mu: vec![0.0; c.inequalities.len()],
        rho: 10.0,
    };
    let mut x = x0;
    let mut previous = f64::INFINITY;
    let mut tol = 1e-6;
    for _ in 0..MAX_OUTER {
        x = subproblem(c, &m, x, lo, hi, tol)?;
        let Merit {
            equalities: hs,
            inequalities: gs,
            ..
        } = merit(c, &m, &x).ok_or(Failure::Domain)?;
        let mut violation = hs.iter().fold(0.0_f64, |a, h| a.max(h.abs()));
        for (g, mu) in gs.iter().zip(&m.mu) {
            violation = violation.max(g.max(-mu / m.rho).abs());
        }
        for (lambda, h) in m.lambda.iter_mut().zip(&hs) {
            *lambda += m.rho * h;
        }
        for (mu, g) in m.mu.iter_mut().zip(&gs) {
            *mu = (*mu + m.rho * g).max(0.0);
        }
        if violation < 1e-10 && tol <= 1e-11 {
            break;
        }
        if violation > previous / 4.0 && m.rho < MAX_PENALTY {
            m.rho *= 10.0;
        }
        previous = violation;
        tol = (tol / 10.0).max(1e-11);
    }
    let point = balls(&x).ok_or(Failure::Domain)?;
    for s in &c.equalities {
        if at(&s.value, &point).ok_or(Failure::Domain)?.abs() > 1e-6 {
            return Err(Failure::Infeasible);
        }
    }
    for s in &c.inequalities {
        if at(&s.value, &point).ok_or(Failure::Domain)? > 1e-6 {
            return Err(Failure::Infeasible);
        }
    }
    Ok(Search {
        x,
        lambda: m.lambda,
        mu: m.mu,
    })
}

fn rational(v: f64) -> Option<Rational> {
    Float::from_f64(v).map(|f| f.to_rational())
}

/// Whether `v` is within [`ACTIVE`] of `end`.
fn near(v: f64, end: f64) -> bool {
    (v - end).abs() <= ACTIVE * (1.0 + end.abs())
}

/// The sign of `e`'s ball: `Some(true)` for certainly positive,
/// `Some(false)` for certainly negative, `None` for undecided.
fn certain_sign(ctx: &Context, e: &Expr) -> Option<bool> {
    let bits = ctx.working_precision();
    let v = approximation(ctx, e, bits).ok()?;
    match v.view() {
        View::Ball(b) => match b.sign()? {
            std::cmp::Ordering::Greater => Some(true),
            std::cmp::Ordering::Less => Some(false),
            std::cmp::Ordering::Equal => None,
        },
        _ => None,
    }
}

/// The bounds and inequalities the search's point meets: each fixed
/// variable's end and side (`true` for the lower), and the active
/// inequalities.
struct Active {
    fixed: Vec<Option<(Expr, bool)>>,
    inequalities: Vec<usize>,
}

fn active(p: &Program, c: &Compiled, s: &Search, lo: &[f64], hi: &[f64]) -> Option<Active> {
    let n = p.xs.len();
    let mut fixed: Vec<Option<(Expr, bool)>> = vec![None; n];
    for j in 0..n {
        if let Some(a) = &p.lower[j]
            && near(s.x[j], lo[j])
        {
            fixed[j] = Some((a.exact.clone(), true));
        } else if let Some(b) = &p.upper[j]
            && near(s.x[j], hi[j])
        {
            fixed[j] = Some((b.exact.clone(), false));
        }
    }
    let point = balls(&s.x)?;
    let mut inequalities = Vec::new();
    for (k, g) in c.inequalities.iter().enumerate() {
        let v = at(&g.value, &point)?;
        if v >= -ACTIVE * (1.0 + s.mu[k].abs()) || s.mu[k] > 0.0 {
            inequalities.push(k);
        }
    }
    Some(Active {
        fixed,
        inequalities,
    })
}

/// The search's point certified as a Karush–Kuhn–Tucker point: the
/// answer, or `None`.
fn certify(
    ev: &mut Eval<'_>,
    p: &Program,
    c: &Compiled,
    s: &Search,
    lo: &[f64],
    hi: &[f64],
) -> Result<Option<Expr>, EvalError> {
    let Some(on) = active(p, c, s, lo, hi) else {
        return Ok(None);
    };
    let (fixed, active) = (&on.fixed, &on.inequalities);
    let n = p.xs.len();
    let rules: Vec<(Expr, Expr)> =
        p.xs.iter()
            .zip(fixed)
            .filter_map(|(x, f)| f.as_ref().map(|(v, _)| (x.clone(), v.clone())))
            .collect();
    let free: Vec<usize> = (0..n).filter(|&j| fixed[j].is_none()).collect();
    let constraints: Vec<Expr> = p
        .equalities
        .iter()
        .chain(active.iter().map(|&k| &p.inequalities[k]))
        .map(|g| g.substitute(&rules))
        .collect();
    if constraints.len() > free.len() {
        return Ok(None);
    }
    let multipliers: Vec<Expr> = (0..constraints.len())
        .map(|i| Expr::symbol(&format!("NLPSolve$lambda{i}")))
        .collect();
    // The Lagrangian in every variable, and with the fixed ones put in.
    let mut full = p.minimized.clone();
    for (lambda, g) in multipliers.iter().zip(
        p.equalities
            .iter()
            .chain(active.iter().map(|&k| &p.inequalities[k])),
    ) {
        full += lambda.clone() * g.clone();
    }
    let full = ev.eval(&full)?;
    let lagrangian = ev.eval(&full.substitute(&rules))?;
    let mut equations = Vec::with_capacity(free.len() + constraints.len());
    let mut unknowns = Vec::with_capacity(free.len() + constraints.len());
    let mut guess = Vec::with_capacity(free.len() + constraints.len());
    for &j in &free {
        equations.push(ev.eval(&Expr::call(
            Symbol::D,
            [lagrangian.clone(), p.xs[j].clone()],
        ))?);
        unknowns.push(p.xs[j].clone());
        let Some(q) = rational(s.x[j]) else {
            return Ok(None);
        };
        guess.push(q);
    }
    let estimates = s.lambda.iter().chain(active.iter().map(|&k| &s.mu[k]));
    for ((lambda, g), &estimate) in multipliers.iter().zip(&constraints).zip(estimates) {
        equations.push(g.clone());
        unknowns.push(lambda.clone());
        let Some(q) = rational(estimate) else {
            return Ok(None);
        };
        guess.push(q);
    }
    let values: Vec<Expr> = if equations.is_empty() {
        Vec::new()
    } else {
        match certified_root(ev, equations, &unknowns, guess)? {
            Some(v) => v,
            None => return Ok(None),
        }
    };
    let solution: Vec<(Expr, Expr)> = unknowns.iter().cloned().zip(values).collect();
    let everything: Vec<(Expr, Expr)> = rules.iter().cloned().chain(solution.clone()).collect();
    let multiplied = &solution[free.len() + p.equalities.len()..];
    if !optimal(ev, p, &on, &full, multiplied, &everything)? {
        return Ok(None);
    }
    let point: Vec<(Expr, Expr)> = everything
        .into_iter()
        .filter(|(x, _)| p.xs.contains(x))
        .collect();
    let bits = ev.context().working_precision();
    let Ok(value) = approximation(ev.context(), &p.objective.substitute(&point), bits) else {
        return Ok(None);
    };
    let value = if matches!(value.view(), View::Ball(_)) {
        value
    } else {
        return Ok(None);
    };
    let coordinates =
        p.xs.iter()
            .map(|x| {
                let v = point
                    .iter()
                    .find(|(y, _)| y == x)
                    .map(|(_, v)| v.clone())
                    .expect("every variable has a value");
                Expr::call(Symbol::EQUATION, [x.clone(), v])
            })
            .collect();
    Ok(Some(list(vec![value, list(coordinates)])))
}

/// Whether the certified Karush–Kuhn–Tucker point is no evident
/// non-optimum: the inequalities' multipliers not certainly negative, no
/// fixed variable's derivative of the Lagrangian `full` certainly of the
/// side that would move it inside, and no inactive inequality certainly
/// violated, at `everything`'s values.
fn optimal(
    ev: &mut Eval<'_>,
    p: &Program,
    on: &Active,
    full: &Expr,
    multiplied: &[(Expr, Expr)],
    everything: &[(Expr, Expr)],
) -> Result<bool, EvalError> {
    if multiplied
        .iter()
        .any(|(_, v)| certain_sign(ev.context(), v) == Some(false))
    {
        return Ok(false);
    }
    for (j, f) in on.fixed.iter().enumerate() {
        let Some((_, lower)) = f else {
            continue;
        };
        let d = ev.eval(&Expr::call(Symbol::D, [full.clone(), p.xs[j].clone()]))?;
        if certain_sign(ev.context(), &d.substitute(everything)) == Some(!*lower) {
            return Ok(false);
        }
    }
    Ok(!p.inequalities.iter().enumerate().any(|(k, g)| {
        !on.inequalities.contains(&k)
            && certain_sign(ev.context(), &g.substitute(everything)) == Some(true)
    }))
}
