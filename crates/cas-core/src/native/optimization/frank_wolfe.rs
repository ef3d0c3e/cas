//! giac's `frank_wolfe` (CAS-609), after optimization.cc `_frank_wolfe`,
//! `FW_backtracking` and `FW_step_size`: the conditional gradient method
//! (M. Frank and P. Wolfe, *An algorithm for quadratic programming*, 1956)
//! with the adaptive step of F. Pedregosa, G. Négiar, A. Askari and M. Jaggi,
//! *Linearly convergent Frank–Wolfe with backtracking line-search*, 2020.
//!
//! - **`FrankWolfe[f, A, …]`**, `f` an expression in at least two symbols
//!   that are no builtins (in canonical order, giac's `_lname` order for
//!   `x1, x2`), `A` a list of at least two points of that dimension with
//!   rational or decimal coordinates: `f` minimised on the convex hull of
//!   `A`. Then at most three options in any order: a list is the start `x0`
//!   (the points' mean), a positive integer `maxiter` (100000), another real
//!   `tol > 0` (`10^−3`), bare `exact` the line search of `find_minimum`
//!   (global, 100 evaluations) and `exact = n` the local one in `n`
//!   evaluations, bare `rand` 1000 random starts and `rand = n` (`n ≥ 2`)
//!   `n` of them. Anything else leaves the call.
//! - Each iteration, giac's: the gradient `∇f(x)` from `D[f, x_i]`; the
//!   vertex `a` of `A` least by `⟨∇f(x), a⟩` gives the direction
//!   `d = a − x`; the step `γ ∈ [0, 1]` by Lipschitz backtracking (giac's
//!   first estimate `L = ‖∇f(x) − ∇f(x + tol·d)‖/(tol·‖d‖)`, then the last
//!   accepted `M`; `τ = 2`, `η = 0.9`, at most 100 halvings), or with
//!   `exact` by [`search`] on `γ ↦ f(x + γd)`; then `x + γd`. In doubles.
//! - The stop, ours: the gap `max_a ⟨∇f(x), x − a⟩` in balls at the working
//!   precision and 32 guard bits, the gradient at the iterate's doubles
//!   taken exactly; the answer once its upper bound is below `tol`. For a
//!   convex `f`, `f(y) ≥ f(x) + ⟨∇f(x), y − x⟩` on the hull, so
//!   `f(x) − min f ≤ gap < tol` is certified (even when the iterate's
//!   doubles lie a rounding outside the hull).
//! - The answer: the iterate's coordinates as exact decimals.
//! - Limits (C0019): `|A|·n` at most `terms`; `maxiter` and `rand`'s count
//!   clamped to `terms/(|A|·n)`; all runs share `terms` evaluations (an
//!   iteration's gap counts `|A|·n`, a gradient `n`, a value of `f` one);
//!   past them a warning and the call stays, or with `rand` the best
//!   certified answer so far, with the warning. The host's interrupt is
//!   polled at each iteration.
//!
//! Departures from giac:
//!
//! - The stop is the certified gap: giac's gap is in doubles.
//! - No certified answer within `maxiter`, or a step that does not move the
//!   iterate: the warning `frank_wolfe: unable to minimize at the given
//!   precision` and the call stays (giac answers the list
//!   `["Unable to minimize at given precision, last value", x]`).
//! - `rand`'s starts are random convex combinations of `A`'s points drawn
//!   from the session's generator, the best certified answer by `f`'s value
//!   (giac draws points with a positive definite Hessian and starts at
//!   their k-means centres).
//! - `f` as giac's list of programs `[f, df(, hf)]` is not accepted, and
//!   `f` must compile to ball operations without `Abs` or `Sign`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the method: f, g, d, x, a, m, n"
)]

use std::cell::Cell;
use std::cmp::Ordering;
use std::sync::atomic::{AtomicBool, Ordering as Atomic};

use cas_num::Rational;
use cas_num::ball::{Ball, Float, Round};
use cas_num::random::Generator;

use super::find_minimum::{Stop, interrupted, search, upper_end, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::ndsolve::rhs::Node;
use crate::symbol::Symbol;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// Bits beyond the working precision for the gap.
const GUARD: u64 = 32;
/// giac's defaults: `maxiter`, `tol`, bare `exact`'s evaluations, bare
/// `rand`'s starts.
const MAXITER: usize = 100_000;
const TOL: f64 = 1e-3;
const EXACT_STEPS: usize = 100;
const STARTS: usize = 1000;
/// giac's backtracking: the factor of `M`, `η`, the most tries.
const TAU: f64 = 2.0;
const ETA: f64 = 0.9;
const BACKTRACKS: usize = 100;

const UNABLE: &str = "frank_wolfe: unable to minimize at the given precision";
const LIMITS: &str = "frank_wolfe: past the limits";

/// The step size rule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Step {
    /// giac's Lipschitz backtracking.
    Backtracking,
    /// [`search`] on `[0, 1]` in `iterations` evaluations.
    Search { iterations: usize, global: bool },
}

impl Step {
    /// The search's evaluations clamped to `terms`.
    fn clamped(self, terms: usize) -> Step {
        match self {
            Step::Search { iterations, global } => Step::Search {
                iterations: iterations.min(terms),
                global,
            },
            Step::Backtracking => Step::Backtracking,
        }
    }
}

/// How a run ended.
enum Outcome {
    /// The certified gap below `tol` at this point.
    Converged(Vec<f64>),
    /// `maxiter` reached, or a step that does not move the iterate.
    Stuck,
    /// The shared evaluations spent.
    Limits,
    /// `f` or its gradient not finite.
    NotReal,
    /// The host's interrupt.
    Interrupted,
}

/// The options after `A`.
struct Options {
    x0: Option<Expr>,
    maxiter: usize,
    tol: f64,
    step: Step,
    starts: usize,
}

/// A positive integer, saturated to `usize::MAX`; `None` for anything else.
fn count(e: &Expr) -> Option<usize> {
    let n = e.to_integer()?;
    if !n.is_positive() {
        return None;
    }
    Some(
        n.to_i64()
            .and_then(|k| usize::try_from(k).ok())
            .unwrap_or(usize::MAX),
    )
}

fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// A finite real number as a double: a rational or a finite decimal.
fn real(e: &Expr) -> Option<f64> {
    let v = match e.view() {
        View::Ball(b) if b.is_finite() => b.mid().to_f64(),
        View::Ball(_) => return None,
        _ => to_f64(&e.to_rational()?),
    };
    v.is_finite().then_some(v)
}

/// Whether `e` is the option `name`, as its builtin (C0008) or bare.
fn is_option(e: &Expr, builtin: Symbol, name: &str) -> bool {
    e.as_symbol()
        .is_some_and(|s| s == builtin || s.name() == name)
}

fn options(given: &[Expr]) -> Option<Options> {
    let mut o = Options {
        x0: None,
        maxiter: MAXITER,
        tol: TOL,
        step: Step::Backtracking,
        starts: 1,
    };
    for g in given {
        let head = g.head().and_then(Expr::as_symbol);
        if let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) = (head, g.args()) {
            if is_option(l, Symbol::EXACT, "exact") {
                o.step = Step::Search {
                    iterations: count(r)?,
                    global: false,
                };
            } else if is_option(l, Symbol::RANDOM, "rand") {
                o.starts = count(r).filter(|&n| n >= 2)?;
            } else {
                return None;
            }
        } else if is_option(g, Symbol::EXACT, "exact") {
            o.step = Step::Search {
                iterations: EXACT_STEPS,
                global: true,
            };
        } else if is_option(g, Symbol::RANDOM, "rand") {
            o.starts = STARTS;
        } else if head == Some(Symbol::LIST) {
            o.x0 = Some(g.clone());
        } else if g.to_integer().is_some() {
            o.maxiter = count(g)?;
        } else {
            let tol = real(g)?;
            if tol.total_cmp(&0.0) != Ordering::Greater {
                return None;
            }
            o.tol = tol;
        }
    }
    Some(o)
}

/// A coordinate: its exact value (a decimal's midpoint as a double), its
/// ball at `prec` and its double; `None` for anything but a rational or a
/// finite decimal (C0028).
fn coordinate(e: &Expr, prec: u64) -> Option<(Rational, Ball, f64)> {
    if let View::Ball(b) = e.view() {
        let d = b.mid().to_f64();
        if !(b.is_finite() && d.is_finite()) {
            return None;
        }
        return Some((Float::from_f64(d)?.to_rational(), b.clone(), d));
    }
    let q = e.to_rational()?;
    let d = to_f64(&q);
    d.is_finite()
        .then(|| (q.clone(), Ball::from_rational(&q, prec), d))
}

/// A list of `n` coordinates.
fn point(e: &Expr, n: usize, prec: u64) -> Option<Vec<(Rational, Ball, f64)>> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) || e.args().len() != n {
        return None;
    }
    e.args().iter().map(|c| coordinate(c, prec)).collect()
}

/// The exact balls of a point's doubles.
fn balls(x: &[f64]) -> Option<Vec<Ball>> {
    x.iter()
        .map(|&v| Float::from_f64(v).map(Ball::exact))
        .collect()
}

/// The problem compiled, and the evaluations shared by its runs.
struct Problem<'a> {
    f: &'a Node,
    gradient: &'a [Node],
    /// `A`'s points as balls (the gap) and as doubles (the direction).
    exact: &'a [Vec<Ball>],
    vertices: &'a [Vec<f64>],
    tol: f64,
    tol_float: Float,
    maxiter: usize,
    step: Step,
    prec: u64,
    flag: &'a AtomicBool,
    spent: Cell<usize>,
    budget: usize,
}

impl Problem<'_> {
    /// Counts `k` evaluations; whether they are within the budget.
    fn charge(&self, k: usize) -> bool {
        let spent = self.spent.get().saturating_add(k);
        self.spent.set(spent);
        spent <= self.budget
    }

    /// `f(x)` in doubles; NaN outside its domain, past the budget or once
    /// interrupted.
    fn value(&self, x: &[f64]) -> f64 {
        if !self.charge(1) || self.flag.load(Atomic::Relaxed) {
            return f64::NAN;
        }
        balls(x)
            .and_then(|b| self.f.ball(&b, 64))
            .map_or(f64::NAN, |v| v.mid().to_f64())
    }

    /// `∇f` at the exact balls `x`, as balls.
    fn gradient_balls(&self, x: &[Ball]) -> Option<Vec<Ball>> {
        self.gradient
            .iter()
            .map(|d| d.ball(x, self.prec).filter(Ball::is_finite))
            .collect()
    }

    /// `∇f(x)` in doubles.
    fn gradient_at(&self, x: &[f64]) -> Option<Vec<f64>> {
        self.charge(x.len());
        let g = self.gradient_balls(&balls(x)?)?;
        g.iter()
            .map(|b| Some(b.mid().to_f64()).filter(|v| v.is_finite()))
            .collect()
    }

    /// Why a value was not finite.
    fn why(&self) -> Outcome {
        if self.flag.load(Atomic::Relaxed) {
            Outcome::Interrupted
        } else if self.spent.get() > self.budget {
            Outcome::Limits
        } else {
            Outcome::NotReal
        }
    }

    /// Whether the certified gap at `x` is below `tol`.
    fn certified(&self, x: &[Ball], g: &[Ball]) -> bool {
        let prec = self.prec;
        let mut most: Option<Float> = None;
        for a in self.exact {
            let mut s = Ball::ZERO;
            for ((gi, xi), ai) in g.iter().zip(x).zip(a) {
                s = s.add(&gi.mul(&xi.sub(ai, prec), prec), prec);
            }
            let Some(up) = upper_end(&s, prec) else {
                return false;
            };
            if most.as_ref().is_none_or(|m| up > *m) {
                most = Some(up);
            }
        }
        most.is_some_and(|m| m < self.tol_float)
    }

    /// giac's backtracking step on `line` (`FW_step_size`), updating the
    /// last accepted `M` and value in `last`.
    fn backtrack(&self, line: &Line<'_>, last: &mut Last) -> Result<f64, Outcome> {
        let (gap, dn2, fk) = (line.gap, line.dn2, line.fk);
        let estimate = if let Some(l) = last.estimate {
            l
        } else {
            let Some(shifted) = self.gradient_at(&line.at(self.tol)) else {
                return Err(self.why());
            };
            let diff: f64 = line
                .g
                .iter()
                .zip(&shifted)
                .map(|(u, v)| (u - v) * (u - v))
                .sum();
            diff.sqrt() / (self.tol * dn2.sqrt())
        };
        let mut m = match last.value {
            None => (ETA + 1.0) * estimate / 2.0,
            Some(fp) => (ETA * estimate).max(estimate.min(gap * gap / (2.0 * (fp - fk) * dn2))),
        };
        if !m.is_finite() {
            return Err(Outcome::NotReal);
        }
        let mut gamma = (gap / (m * dn2)).min(1.0);
        for _ in 0..BACKTRACKS {
            let q = fk - gamma * gap + gamma * gamma * m * dn2 / 2.0;
            if q >= self.value(&line.at(gamma)) {
                break;
            }
            m *= TAU;
            gamma = (gap / (m * dn2)).min(1.0);
        }
        if self.spent.get() > self.budget || self.flag.load(Atomic::Relaxed) {
            return Err(self.why());
        }
        *last = Last {
            estimate: Some(m),
            value: Some(fk),
        };
        Ok(gamma)
    }

    /// The step on `line` by [`search`] on `[0, 1]`.
    fn exact_step(&self, line: &Line<'_>, iterations: usize, global: bool) -> Result<f64, Outcome> {
        let along = |gamma: f64| self.value(&line.at(gamma));
        let pops = (self.budget.saturating_sub(self.spent.get()) / iterations).max(1);
        match search(&along, 0.0, 1.0, self.tol, iterations, global, pops) {
            Ok((gamma, _)) => Ok(gamma),
            Err(Stop::Limits) => Err(Outcome::Limits),
            Err(Stop::NotReal) => Err(self.why()),
        }
    }

    /// The certified answer at `x`, or the line from `x` toward the vertex
    /// least by `⟨∇f(x), a⟩` (the first of ties).
    fn direction(&self, x: &[f64], fk: f64) -> Result<(Vec<f64>, Vec<f64>, f64), Outcome> {
        let n = x.len();
        if !self.charge(self.exact.len() * n + n) {
            return Err(Outcome::Limits);
        }
        let Some(xb) = balls(x) else {
            return Err(Outcome::NotReal);
        };
        let Some(gb) = self.gradient_balls(&xb) else {
            return Err(Outcome::NotReal);
        };
        if self.certified(&xb, &gb) {
            return Err(Outcome::Converged(x.to_vec()));
        }
        let g: Vec<f64> = gb.iter().map(|b| b.mid().to_f64()).collect();
        if !(fk.is_finite() && g.iter().all(|v| v.is_finite())) {
            return Err(Outcome::NotReal);
        }
        let mut best: Option<(f64, Vec<f64>)> = None;
        for a in self.vertices {
            let d: Vec<f64> = a.iter().zip(x).map(|(ai, xi)| ai - xi).collect();
            let t: f64 = g.iter().zip(&d).map(|(gi, di)| gi * di).sum();
            if best.as_ref().is_none_or(|(b, _)| t < *b) {
                best = Some((t, d));
            }
        }
        let Some((t, d)) = best else {
            return Err(Outcome::NotReal);
        };
        Ok((g, d, -t))
    }

    /// The runs from `origins`: the certified answer least by `f`'s value,
    /// whether the limits stopped them, and whether one was stuck.
    fn best(&self, origins: Vec<Vec<f64>>) -> (Option<Vec<f64>>, bool, bool) {
        let mut best: Option<(f64, Vec<f64>)> = None;
        let (mut limited, mut stuck) = (false, false);
        for origin in origins {
            match self.run(origin) {
                Outcome::Converged(x) => {
                    let v = self.value(&x);
                    if best
                        .as_ref()
                        .is_none_or(|(b, _)| v.total_cmp(b) == Ordering::Less)
                    {
                        best = Some((v, x));
                    }
                }
                Outcome::Stuck => stuck = true,
                Outcome::NotReal => {}
                Outcome::Limits => {
                    limited = true;
                    break;
                }
                Outcome::Interrupted => break,
            }
        }
        (best.map(|(_, x)| x), limited, stuck)
    }

    /// One run of giac's `FW_backtracking` from `x`.
    fn run(&self, mut x: Vec<f64>) -> Outcome {
        let mut last = Last {
            estimate: None,
            value: None,
        };
        for _ in 0..self.maxiter {
            if self.flag.load(Atomic::Relaxed) {
                return Outcome::Interrupted;
            }
            let fk = self.value(&x);
            if !fk.is_finite() {
                return self.why();
            }
            let (g, d, gap) = match self.direction(&x, fk) {
                Ok(found) => found,
                Err(end) => return end,
            };
            let dn2: f64 = d.iter().map(|v| v * v).sum();
            if !(gap > 0.0 && dn2 > 0.0 && dn2.is_finite()) {
                return Outcome::Stuck;
            }
            let line = Line {
                x: &x,
                d: &d,
                g: &g,
                gap,
                dn2,
                fk,
            };
            let gamma = match self.step {
                Step::Search { iterations, global } => self.exact_step(&line, iterations, global),
                Step::Backtracking => self.backtrack(&line, &mut last),
            };
            let gamma = match gamma {
                Ok(gamma) => gamma,
                Err(end) => return end,
            };
            let next = line.at(gamma);
            let moved = next
                .iter()
                .zip(&x)
                .any(|(a, b)| a.total_cmp(b) != Ordering::Equal);
            if !(gamma > 0.0 && moved && next.iter().all(|v| v.is_finite())) {
                return Outcome::Stuck;
            }
            x = next;
        }
        Outcome::Stuck
    }
}

/// The line `x + γd` of an iteration: the gradient `g` at `x`, the gap
/// `⟨g, x − a⟩` in doubles, `‖d‖²` and `f(x)`.
struct Line<'a> {
    x: &'a [f64],
    d: &'a [f64],
    g: &'a [f64],
    gap: f64,
    dn2: f64,
    fk: f64,
}

impl Line<'_> {
    fn at(&self, gamma: f64) -> Vec<f64> {
        self.x
            .iter()
            .zip(self.d)
            .map(|(xi, di)| xi + gamma * di)
            .collect()
    }
}

/// The backtracking's last accepted `M` and `f`'s last value.
struct Last {
    estimate: Option<f64>,
    value: Option<f64>,
}

/// A uniform double in `[0, 1)` from 64 random bits.
#[expect(clippy::cast_precision_loss, reason = "53 bits, exact in a double")]
fn uniform(bits: u64) -> f64 {
    (bits >> 11) as f64 * f64::powi(2.0, -53)
}

/// The start: `x0` given, else the mean of the points' exact values.
fn start(
    x0: Option<&Expr>,
    parsed: &[Vec<(Rational, Ball, f64)>],
    n: usize,
    prec: u64,
) -> Option<Vec<f64>> {
    if let Some(x0) = x0 {
        return Some(point(x0, n, prec)?.into_iter().map(|c| c.2).collect());
    }
    let m = Rational::from(i64::try_from(parsed.len()).ok()?);
    (0..n)
        .map(|j| {
            let mut sum = Rational::from(0);
            for p in parsed {
                sum = &sum + &p[j].0;
            }
            Some(to_f64(&(&sum / &m))).filter(|v| v.is_finite())
        })
        .collect()
}

/// `f` and its gradient `D[f, x_i]` compiled in `vars`.
fn compile(
    ev: &mut Eval<'_>,
    f: &Expr,
    vars: &[Expr],
    prec: u64,
) -> Result<Option<(Node, Vec<Node>)>, EvalError> {
    let Some(tree) = Node::compile(ev.context(), f, vars, prec) else {
        return Ok(None);
    };
    let mut gradient = Vec::with_capacity(vars.len());
    for x in vars {
        let d = ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?;
        let Some(node) = Node::compile(ev.context(), &d, vars, prec) else {
            return Ok(None);
        };
        gradient.push(node);
    }
    Ok(Some((tree, gradient)))
}

/// `count` random convex combinations of the vertices, their weights
/// uniform in `[0, 1)` from the session's generator (`x0` for all zero).
fn draw(ev: &mut Eval<'_>, vertices: &[Vec<f64>], x0: &[f64], count: usize) -> Vec<Vec<f64>> {
    (0..count)
        .map(|_| {
            let weights: Vec<f64> = (0..vertices.len())
                .map(|_| uniform(ev.draw(Generator::next_u64)))
                .collect();
            let total: f64 = weights.iter().sum();
            (0..x0.len())
                .map(|j| {
                    let sum: f64 = weights.iter().zip(vertices).map(|(w, a)| w * a[j]).sum();
                    if total > 0.0 { sum / total } else { x0[j] }
                })
                .collect()
        })
        .collect()
}

/// giac's `frank_wolfe`: the certified minimiser's coordinates; `None`
/// where the call stays.
pub(in crate::native) fn frank_wolfe(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, points, given @ ..] = e.args() else {
        return Ok(None);
    };
    if given.len() > 3 || f.head() == Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let Some(o) = options(given) else {
        return Ok(None);
    };
    let vars: Vec<Expr> = f
        .symbols()
        .into_iter()
        .filter(|s| !s.is_builtin())
        .map(Expr::from)
        .collect();
    let n = vars.len();
    if n < 2 || points.head() != Some(&Expr::from(Symbol::LIST)) || points.args().len() < 2 {
        return Ok(None);
    }
    let wp = ev.context().working_precision() + GUARD;
    let Some(parsed) = points
        .args()
        .iter()
        .map(|p| point(p, n, wp))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let Some(x0) = start(o.x0.as_ref(), &parsed, n, wp) else {
        return Ok(None);
    };
    // The work bounded before any, as a product (C0019).
    let terms = ev.context().limits().terms;
    let per = match parsed.len().checked_mul(n) {
        Some(k) if k <= terms => k,
        _ => {
            warn(ev, e, LIMITS);
            return Ok(None);
        }
    };
    let maxiter = o.maxiter.min(terms / per);
    let starts = o.starts.min(terms / per);
    let Some((tree, gradient)) = compile(ev, f, &vars, wp)? else {
        return Ok(None);
    };
    let exact: Vec<Vec<Ball>> = parsed
        .iter()
        .map(|p| p.iter().map(|c| c.1.clone()).collect())
        .collect();
    let vertices: Vec<Vec<f64>> = parsed
        .iter()
        .map(|p| p.iter().map(|c| c.2).collect())
        .collect();
    let Some(tol_float) = Float::from_f64(o.tol) else {
        return Ok(None);
    };
    // The random starts, drawn before any run.
    let origins = if starts > 1 {
        draw(ev, &vertices, &x0, starts)
    } else {
        vec![x0]
    };
    let flag = ev.context().interrupt_flag();
    let problem = Problem {
        f: &tree,
        gradient: &gradient,
        exact: &exact,
        vertices: &vertices,
        tol: o.tol,
        tol_float,
        maxiter,
        step: o.step.clamped(terms),
        prec: wp,
        flag,
        spent: Cell::new(0),
        budget: terms,
    };
    let (best, limited, stuck) = problem.best(origins);
    interrupted(ev)?;
    if limited {
        warn(ev, e, LIMITS);
    }
    let Some(x) = best else {
        if !limited && stuck {
            warn(ev, e, UNABLE);
        }
        return Ok(None);
    };
    let Some(coordinates) = x
        .into_iter()
        .map(|v| Float::from_f64(v).map(|z| Expr::from(Ball::exact(z))))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    Ok(Some(Expr::call(Symbol::LIST, coordinates)))
}
