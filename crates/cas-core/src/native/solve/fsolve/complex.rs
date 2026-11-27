//! `CFSolve` (CAS-536), Xcas's `cfsolve`: numeric complex roots, as balls
//! at the working precision that each hold one root.
//!
//! - `CFSolve[eqn, x]`, of a rational function of `x` alone with exact
//!   coefficients: each root of its numerator's square-free parts once, by
//!   `polynomial_roots` as `FSolve`'s, ordered by real part then imaginary
//!   part, the real ones real (CAS-211): giac's `cfsolve(x^2+1=0,x)` is
//!   `[-i,i]`. A list, `[]` without a root.
//! - Otherwise, and for `CFSolve[{eqns}, {vars}]`, candidates: the start
//!   point `x = x0` or `{x = x0, …}`, and when it is real and gives no
//!   root, the start moved by `i/2`, off the real axis Newton's iteration
//!   of a real function keeps to; else `CSolve`'s solutions,
//!   approximated; else [`START`] for each unknown. Each is refined by
//!   Newton's iteration and certified by Krawczyk's test on the `2n` real
//!   system of the real and imaginary parts ([`super::systems`]); a root
//!   whose imaginary parts' balls hold 0 is real when the test certifies it
//!   on the real system. Candidates that are not certified go, with a
//!   warning; roots found twice count once.
//! - The roots in the candidates' order: a list of numbers for one unknown,
//!   of tuples for a list of unknowns; from a start point, the root alone,
//!   as giac's `cfsolve(x^2+1=0,x=1+i)` is `i`. No root: `[]` when `CSolve`
//!   found none, else the call stays.
//! - Parameters, a system that is not square, or unknowns that are not
//!   distinct symbols leave the call.
//!
//! CAS-650: `point` refuses what is not numeric before `approximation`.
//! CAS-650: a function holding a nested own-head call (depth >= 2) leaves the call.

use cas_num::ball::ComplexBall;
use cas_num::{Integer, Rational};

use super::super::super::{approximation, holds_nested};
use super::super::sides;
use super::systems::{MAX_UNKNOWNS, certified, complex, distinct_symbols, listed};
use super::{list, polynomial};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The start of Newton's iteration for each unknown without one, off the
/// real axis: `1/3 + i/2`.
const START: (i64, i64, i64, i64) = (1, 3, 1, 2);

/// `CFSolve[eqn, x]`, `CFSolve[eqn, x == x0]`, `CFSolve[{eqns}, {vars}]` or
/// `CFSolve[{eqns}, {x == x0, …}]`.
pub(in crate::native) fn cfsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [eqns, spec] = e.args() else {
        return Ok(None);
    };
    let one = listed(spec).is_none();
    let specs = listed(spec).unwrap_or(std::slice::from_ref(spec));
    let eqns = match listed(eqns) {
        Some(eqns) if !one => eqns,
        None if one => std::slice::from_ref(eqns),
        _ => return Ok(None),
    };
    let starts: Vec<Option<&Expr>> = specs.iter().map(|s| sides(s).map(|(_, x0)| x0)).collect();
    let vars: Vec<Expr> = specs
        .iter()
        .map(|s| sides(s).map_or(s, |(x, _)| x).clone())
        .collect();
    let n = vars.len();
    let with_start = starts.iter().all(Option::is_some);
    if n == 0
        || n > MAX_UNKNOWNS
        || eqns.len() != n
        || !distinct_symbols(&vars)
        || !(with_start || starts.iter().all(Option::is_none))
    {
        return Ok(None);
    }
    let mut fs = Vec::with_capacity(n);
    for eqn in eqns {
        fs.push(match sides(eqn) {
            Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
            None => eqn.clone(),
        });
    }
    if fs.iter().any(|f| holds_nested(f, Symbol::CFSOLVE)) || !numeric(ev.context(), &fs, &vars) {
        return Ok(None);
    }
    let bits = ev.context().working_precision();
    if one
        && !with_start
        && let Some(roots) = polynomial(ev, &fs[0], &vars[0], bits, true)?
    {
        let result = list(roots);
        ev.step("CFSolve", e, &result);
        return Ok(Some(result));
    }
    let (candidates, none) = if with_start {
        let Some(candidates) = started(ev.context(), starts) else {
            return Ok(None);
        };
        (candidates, false)
    } else {
        solved(ev, eqns, spec, n)?
    };
    let roots = found(ev, e, &fs, &vars, candidates, with_start)?;
    let number = |z: ComplexBall| {
        if z.is_real() {
            Expr::from(z.re().clone())
        } else {
            Expr::from(z)
        }
    };
    let mut values: Vec<Expr> = roots
        .into_iter()
        .map(|r| {
            let mut r = r.into_iter().map(number);
            if one {
                r.next().unwrap_or_else(|| Expr::from(0))
            } else {
                list(r.collect())
            }
        })
        .collect();
    let result = match values.len() {
        0 if none => list(Vec::new()),
        0 => return Ok(None),
        1 if with_start => values.pop().expect("one root"),
        _ => list(values),
    };
    ev.step("CFSolve", e, &result);
    Ok(Some(result))
}

/// The candidates from the start points `x0`, `[Re x0, Im x0, …]`: the
/// point, and when it is real, the point moved by `i/2` in each unknown.
fn started(ctx: &Context, starts: Vec<Option<&Expr>>) -> Option<Vec<Vec<Rational>>> {
    let mut q = Vec::with_capacity(2 * starts.len());
    for x0 in starts.into_iter().flatten() {
        let (re, im) = point(ctx, x0)?;
        q.extend([re, im]);
    }
    let mut candidates = vec![q.clone()];
    if q.iter().skip(1).step_by(2).all(Rational::is_zero) {
        let (_, _, c, d) = START;
        let im = Rational::new(Integer::from(c), Integer::from(d));
        for b in q.iter_mut().skip(1).step_by(2) {
            *b = im.clone();
        }
        candidates.push(q);
    }
    Some(candidates)
}

/// The certified roots from `candidates`, each once; from a start point
/// the first. A warning for the candidates that give none, unless a start
/// point's second gives one.
fn found(
    ev: &mut Eval<'_>,
    e: &Expr,
    fs: &[Expr],
    vars: &[Expr],
    candidates: Vec<Vec<Rational>>,
    with_start: bool,
) -> Result<Vec<Vec<ComplexBall>>, EvalError> {
    let mut roots: Vec<Vec<ComplexBall>> = Vec::new();
    let mut dropped = 0_usize;
    for q in candidates {
        if with_start && !roots.is_empty() {
            break;
        }
        match root(ev, fs, vars, q)? {
            Some(r) => {
                if !roots
                    .iter()
                    .any(|s| s.iter().zip(&r).all(|(a, b)| a.overlaps(b)))
                {
                    roots.push(r);
                }
            }
            None => dropped += 1,
        }
    }
    if dropped > 0 && (roots.is_empty() || !with_start) {
        ev.warn(Warning::Native {
            message: format!("Newton's iteration certified no root from {dropped} candidate(s)")
                .into(),
            expr: e.clone(),
        });
    }
    Ok(roots)
}

/// Whether every `f` is a finite number at the point `1/3 + i/2, …`: no
/// parameters.
fn numeric(ctx: &Context, fs: &[Expr], vars: &[Expr]) -> bool {
    let rules: Vec<(Expr, Expr)> = vars.iter().map(|x| (x.clone(), start())).collect();
    fs.iter()
        .all(|f| complex(ctx, &f.substitute(&rules), 64).is_some())
}

/// [`START`] as a number.
fn start() -> Expr {
    let (a, b, c, d) = START;
    Expr::from(Rational::new(Integer::from(a), Integer::from(b)))
        + Expr::call(
            Symbol::TIMES,
            [
                Expr::from(Rational::new(Integer::from(c), Integer::from(d))),
                Expr::from(cas_num::Complex::I),
            ],
        )
}

/// A number's real and imaginary parts as rationals: its ball's midpoints.
fn point(ctx: &Context, e: &Expr) -> Option<(Rational, Rational)> {
    if !e.is_numeric() {
        return None;
    }
    match approximation(ctx, e, 128).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some((b.mid().to_rational(), Rational::ZERO)),
        View::ComplexBall(z) if z.is_finite() => {
            Some((z.re().mid().to_rational(), z.im().mid().to_rational()))
        }
        _ => None,
    }
}

/// The candidates from `CSolve`'s solutions, `[Re x_0, Im x_0, …]` each,
/// or [`START`] when it gives none that is a number; and whether `CSolve`
/// found no solution.
#[allow(clippy::type_complexity, reason = "the candidates and a flag")]
fn solved(
    ev: &mut Eval<'_>,
    eqns: &[Expr],
    spec: &Expr,
    n: usize,
) -> Result<(Vec<Vec<Rational>>, bool), EvalError> {
    let eqn = if n == 1 && listed(spec).is_none() {
        eqns[0].clone()
    } else {
        list(eqns.to_vec())
    };
    let exact = ev.eval(&Expr::call(Symbol::CSOLVE, [eqn, spec.clone()]))?;
    let mut candidates = Vec::new();
    let mut none = false;
    if let Some(solutions) = listed(&exact) {
        none = solutions.is_empty();
        for s in solutions {
            let parts = if listed(spec).is_some() {
                listed(s).filter(|p| p.len() == n)
            } else {
                Some(std::slice::from_ref(s))
            };
            let Some(parts) = parts else {
                continue;
            };
            let mut q = Vec::with_capacity(2 * n);
            for part in parts {
                let Some((re, im)) = point(ev.context(), part) else {
                    break;
                };
                q.extend([re, im]);
            }
            if q.len() == 2 * n {
                candidates.push(q);
            }
        }
    }
    if candidates.is_empty() && !none {
        let (a, b, c, d) = START;
        let z = [
            Rational::new(Integer::from(a), Integer::from(b)),
            Rational::new(Integer::from(c), Integer::from(d)),
        ];
        candidates.push(z.iter().cycle().take(2 * n).cloned().collect());
    }
    Ok((candidates, none))
}

/// The root Newton's iteration reaches from `q`, certified: each unknown's
/// ball, real when the real system's test certifies it.
fn root(
    ev: &mut Eval<'_>,
    fs: &[Expr],
    vars: &[Expr],
    q: Vec<Rational>,
) -> Result<Option<Vec<ComplexBall>>, EvalError> {
    let Some(balls) = certified(ev, fs.to_vec(), vars, q, true)? else {
        return Ok(None);
    };
    let zs: Vec<ComplexBall> = balls
        .chunks(2)
        .map(|c| ComplexBall::new(c[0].clone(), c[1].clone()))
        .collect();
    if zs.iter().all(|z| z.im().contains_rational(&Rational::ZERO)) {
        let re = zs.iter().map(|z| z.re().mid().to_rational()).collect();
        if let Some(real) = certified(ev, fs.to_vec(), vars, re, false)? {
            return Ok(Some(
                real.into_iter()
                    .map(ComplexBall::from_real)
                    .collect::<Vec<_>>(),
            ));
        }
    }
    Ok(Some(zs))
}
