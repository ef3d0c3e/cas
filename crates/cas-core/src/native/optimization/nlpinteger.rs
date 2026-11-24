//! Integer and binary variables in nonlinear programs (CAS-605): `NLPSolve`
//! with `assume = integer` (`nlp_integer`), `nlp_nonnegint`, `nlp_binary`,
//! `nlp_integervariables` or `nlp_binaryvariables` on a program that is
//! not linear (a linear one is `LPSolve`'s, [`super::nlpsolve`]).
//!
//! - The integer variables' ends are integral already: a binary variable's
//!   within `0..1`, the lower end's ceiling, the upper end's floor (giac's,
//!   optimization.cc:6861-6876).
//! - When every integer variable's range is finite and their box, the
//!   product of the ranges' sizes, is at most [`MAX_BOX`] points
//!   ([`MAX_MIXED_BOX`] when continuous variables remain), every point is
//!   tried in lexicographic order, the variables sorted and the first the
//!   most significant. At a point the constraints in the integers alone
//!   are decided exactly (a value that is no rational by its certified
//!   ball's sign), the objective evaluated exactly; with continuous
//!   variables the rest of the program, the integers put in, is
//!   `NLPSolve`'s certified search's, its value a ball. The first strictly
//!   best point is the answer, values compared exactly or by their balls,
//!   an undecided comparison keeping the earlier.
//! - Otherwise a depth-first branch and bound on the relaxations, each
//!   solved by the certified search, at most [`MAX_NODES`] of them (and
//!   the evaluator's iterations): a relaxation whose integer variables are
//!   all within `10⁻⁶·(1 + |m|)` of integers is rounded there and the
//!   point tried as above; else the first that is not is branched on, its
//!   floor's side first; a relaxation certainly no better than the best
//!   point is pruned. As giac's (optimization.cc:5972-6102), it is local.
//! - The answer is giac's `{value, {x = …, …}}`, the integer coordinates
//!   exact, the value exact where the objective is at an integer point. No
//!   feasible integer point leaves the call with the warning "no integer
//!   point is feasible"; a point a constraint or the objective leaves
//!   undecided is skipped with a warning; a continuous search that is not
//!   certified is dropped with a warning that the answer is local; the node
//!   cap answers the best point with a warning, or leaves the call with one.
//!   A root relaxation that fails, or a search that runs off, leaves the
//!   call as `NLPSolve` without integers does.

use std::cmp::Ordering;

use cas_num::Rational;

use super::lpsolve::symbols_in_order;
use super::nlpsolve::{
    End, Native, Program, Solved, certain_sign, continuous, double, list, rational, warn,
};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most points enumerated when every variable is an integer.
const MAX_BOX: u64 = 10_000;

/// The most points enumerated when continuous variables remain, each a
/// certified search.
const MAX_MIXED_BOX: u64 = 64;

/// The most relaxations the branch and bound solves.
const MAX_NODES: u32 = 256;

/// Within this (relative) distance a relaxation's coordinate is an
/// integer.
const NEAR: f64 = 1e-6;

/// What an integer point gives.
enum Point {
    /// Its value and coordinates, `x = v` in the program's order.
    Value(Expr, Vec<Expr>),
    /// Infeasible: certainly, or its search found no feasible point.
    Infeasible,
    /// A constraint or the objective undecided there.
    Undecided,
    /// Its continuous search not certified.
    Dropped,
    /// Its continuous search ran off to infinity.
    Runaway(&'static str),
}

/// How a search ended.
enum Ended {
    Searched,
    /// At the node cap.
    Capped,
    /// The call stays, with this warning or none.
    Stays(Option<&'static str>),
}

/// The best point so far, and what the answer warns of.
#[derive(Default)]
struct Best {
    point: Option<(Expr, Vec<Expr>)>,
    undecided: bool,
    dropped: bool,
}

/// `NLPSolve` on a program with an integer variable.
pub(super) fn nlp_integer(ev: &mut Eval<'_>, e: &Expr, p: &Program) -> Native {
    let ints: Vec<usize> = (0..p.xs.len()).filter(|&j| p.integer[j]).collect();
    let cap = if ints.len() == p.xs.len() {
        MAX_BOX
    } else {
        MAX_MIXED_BOX
    };
    let mut best = Best::default();
    let ended = match points(p, &ints, cap) {
        Some((lows, sizes, total)) => enumerate(ev, p, &ints, &lows, &sizes, total, &mut best)?,
        None => branch_and_bound(ev, p, &ints, &mut best)?,
    };
    if let Ended::Stays(why) = ended {
        if let Some(why) = why {
            warn(ev, e, why);
        }
        return Ok(None);
    }
    if best.undecided {
        warn(
            ev,
            e,
            "a constraint or the objective is undecided at an integer point: the point is skipped",
        );
    }
    if best.dropped {
        warn(
            ev,
            e,
            "a search could not be certified and was dropped: the answer is local",
        );
    }
    let capped = matches!(ended, Ended::Capped);
    let Some((value, coordinates)) = best.point else {
        warn(
            ev,
            e,
            if capped {
                "the node cap was reached before an integer point was found"
            } else {
                "no integer point is feasible"
            },
        );
        return Ok(None);
    };
    if capped {
        warn(
            ev,
            e,
            "the node cap was reached; the best integer point found is answered",
        );
    }
    let answer = list(vec![value, list(coordinates)]);
    ev.step("NLPSolve", e, &answer);
    Ok(Some(answer))
}

/// The integer variables' box when it is finite and at most `cap` points:
/// each lower end, each range's size and their product.
fn points(p: &Program, ints: &[usize], cap: u64) -> Option<(Vec<Rational>, Vec<u64>, u64)> {
    let mut lows = Vec::with_capacity(ints.len());
    let mut sizes = Vec::with_capacity(ints.len());
    let mut total: u64 = 1;
    for &j in ints {
        let low = p.lower[j].as_ref()?.exact.to_rational()?;
        let high = p.upper[j].as_ref()?.exact.to_rational()?;
        let size = &(&high - &low) + &Rational::ONE;
        let size = if size.is_positive() {
            u64::try_from(size.numer().to_i64()?).ok()?
        } else {
            0
        };
        total = total.checked_mul(size)?;
        lows.push(low);
        sizes.push(size);
    }
    (total <= cap).then_some((lows, sizes, total))
}

/// `Err(Interrupted)` when the host asked to stop.
fn interrupted(ev: &Eval<'_>) -> Result<(), EvalError> {
    if ev
        .context()
        .interrupt_flag()
        .load(std::sync::atomic::Ordering::Relaxed)
    {
        Err(EvalError::Interrupted)
    } else {
        Ok(())
    }
}

/// Every point of the box, in lexicographic order.
fn enumerate(
    ev: &mut Eval<'_>,
    p: &Program,
    ints: &[usize],
    lows: &[Rational],
    sizes: &[u64],
    total: u64,
    best: &mut Best,
) -> Result<Ended, EvalError> {
    for k in 0..total {
        interrupted(ev)?;
        let mut rest = k;
        let mut values = vec![Rational::ZERO; ints.len()];
        for i in (0..ints.len()).rev() {
            values[i] = &lows[i] + &Rational::from(rest % sizes[i]);
            rest /= sizes[i];
        }
        if let Some(why) = consider(ev, p, ints, &values, best)? {
            return Ok(Ended::Stays(Some(why)));
        }
    }
    Ok(Ended::Searched)
}

/// The point `values` tried, the best updated; the warning when its
/// search ran off.
fn consider(
    ev: &mut Eval<'_>,
    p: &Program,
    ints: &[usize],
    values: &[Rational],
    best: &mut Best,
) -> Result<Option<&'static str>, EvalError> {
    match point(ev, p, ints, values)? {
        Point::Value(value, coordinates) => {
            let better = match &best.point {
                None => true,
                Some((incumbent, _)) => {
                    let d = ev.eval(&(value.clone() - incumbent.clone()))?;
                    sign(ev, &d)
                        == Some(if p.maximum {
                            Ordering::Greater
                        } else {
                            Ordering::Less
                        })
                }
            };
            if better {
                best.point = Some((value, coordinates));
            }
        }
        Point::Undecided => best.undecided = true,
        Point::Dropped => best.dropped = true,
        Point::Runaway(why) => return Ok(Some(why)),
        Point::Infeasible => {}
    }
    Ok(None)
}

/// `v`'s sign: exact for a rational, else its certified ball's; `None`
/// when undecided.
fn sign(ev: &Eval<'_>, v: &Expr) -> Option<Ordering> {
    if let Some(q) = v.to_rational() {
        return Some(if q.is_zero() {
            Ordering::Equal
        } else if q.is_positive() {
            Ordering::Greater
        } else {
            Ordering::Less
        });
    }
    certain_sign(ev.context(), v).map(|positive| {
        if positive {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    })
}

/// The program at the integer point `values` (of the variables `ints`).
fn point(
    ev: &mut Eval<'_>,
    p: &Program,
    ints: &[usize],
    values: &[Rational],
) -> Result<Point, EvalError> {
    let rules: Vec<(Expr, Expr)> = ints
        .iter()
        .zip(values)
        .map(|(&j, v)| (p.xs[j].clone(), Expr::from(v.clone())))
        .collect();
    let mut equalities = Vec::new();
    let mut inequalities = Vec::new();
    let all = p
        .equalities
        .iter()
        .map(|g| (g, true))
        .chain(p.inequalities.iter().map(|g| (g, false)));
    for (g, equality) in all {
        let g = ev.eval(&g.substitute(&rules))?;
        let mut symbols = Vec::new();
        symbols_in_order(&g, &mut symbols);
        if !symbols.is_empty() {
            if equality {
                equalities.push(g);
            } else {
                inequalities.push(g);
            }
            continue;
        }
        match (sign(ev, &g), equality) {
            (None, _) => return Ok(Point::Undecided),
            (Some(Ordering::Equal), _) | (Some(Ordering::Less), false) => {}
            _ => return Ok(Point::Infeasible),
        }
    }
    let free: Vec<usize> = (0..p.xs.len()).filter(|&j| !p.integer[j]).collect();
    let objective = ev.eval(&p.objective.substitute(&rules))?;
    if free.is_empty() {
        if double(ev.context(), &objective).is_none() {
            return Ok(Point::Undecided);
        }
        let coordinates = rules
            .into_iter()
            .map(|(x, v)| Expr::call(Symbol::EQUATION, [x, v]))
            .collect();
        return Ok(Point::Value(objective, coordinates));
    }
    let pick = |ends: &[Option<End>]| free.iter().map(|&j| ends[j].clone()).collect();
    let rest = Program {
        xs: free.iter().map(|&j| p.xs[j].clone()).collect(),
        objective,
        minimized: ev.eval(&p.minimized.substitute(&rules))?,
        equalities,
        inequalities,
        lower: pick(&p.lower),
        upper: pick(&p.upper),
        start: free.iter().map(|&j| p.start[j].clone()).collect(),
        limit: p.limit,
        precision: p.precision,
        maximum: p.maximum,
        integer: vec![false; free.len()],
    };
    let answer = match continuous(ev, &rest)? {
        Solved::Optimum(answer) => answer,
        Solved::Infeasible(_) => return Ok(Point::Infeasible),
        Solved::Runaway(why) => return Ok(Point::Runaway(why)),
        Solved::Failed(_) => return Ok(Point::Dropped),
    };
    let [value, found] = answer.args() else {
        return Ok(Point::Dropped);
    };
    let mut fixed = rules.into_iter();
    let mut solved = found.args().iter();
    let mut coordinates = Vec::with_capacity(p.xs.len());
    for &integer in &p.integer {
        let c = if integer {
            fixed
                .next()
                .map(|(x, v)| Expr::call(Symbol::EQUATION, [x, v]))
        } else {
            solved.next().cloned()
        };
        let Some(c) = c else {
            return Ok(Point::Dropped);
        };
        coordinates.push(c);
    }
    Ok(Point::Value(value.clone(), coordinates))
}

/// An integral double's end.
fn end(v: f64) -> Option<End> {
    Some(End {
        exact: Expr::from(rational(v)?),
        value: v,
    })
}

/// The depth-first branch and bound on the certified relaxations.
fn branch_and_bound(
    ev: &mut Eval<'_>,
    p: &Program,
    ints: &[usize],
    best: &mut Best,
) -> Result<Ended, EvalError> {
    let cap = MAX_NODES.min(ev.context().limits().iterations);
    let mut stack = vec![(p.lower.clone(), p.upper.clone())];
    let mut nodes: u32 = 0;
    while let Some((lower, upper)) = stack.pop() {
        if nodes == cap {
            return Ok(Ended::Capped);
        }
        interrupted(ev)?;
        let root = nodes == 0;
        nodes += 1;
        let relaxation = Program {
            lower: lower.clone(),
            upper: upper.clone(),
            ..p.clone()
        };
        let answer = match continuous(ev, &relaxation)? {
            Solved::Optimum(answer) => answer,
            Solved::Infeasible(why) | Solved::Runaway(why) | Solved::Failed(Some(why)) if root => {
                return Ok(Ended::Stays(Some(why)));
            }
            Solved::Failed(None) if root => return Ok(Ended::Stays(None)),
            Solved::Infeasible(_) => continue,
            Solved::Runaway(_) | Solved::Failed(_) => {
                best.dropped = true;
                continue;
            }
        };
        let [value, found] = answer.args() else {
            best.dropped = true;
            continue;
        };
        if let Some((incumbent, _)) = &best.point {
            let d = ev.eval(&(value.clone() - incumbent.clone()))?;
            let worse = if p.maximum {
                Ordering::Less
            } else {
                Ordering::Greater
            };
            if matches!(sign(ev, &d), Some(s) if s == worse || s == Ordering::Equal) {
                continue;
            }
        }
        let mids: Option<Vec<f64>> = ints
            .iter()
            .map(|&j| {
                let v = found.args().get(j)?.args().get(1)?;
                double(ev.context(), v)
            })
            .collect();
        let Some(mids) = mids else {
            best.dropped = true;
            continue;
        };
        let Some(i) = mids
            .iter()
            .position(|m| (m - m.round()).abs() > NEAR * (1.0 + m.abs()))
        else {
            let values: Option<Vec<Rational>> = mids.iter().map(|m| rational(m.round())).collect();
            let Some(values) = values else {
                best.dropped = true;
                continue;
            };
            if let Some(why) = consider(ev, p, ints, &values, best)? {
                return Ok(Ended::Stays(Some(why)));
            }
            continue;
        };
        let j = ints[i];
        let f = mids[i].floor();
        let (Some(down), Some(up)) = (end(f), end(f + 1.0)) else {
            best.dropped = true;
            continue;
        };
        let mut ceiling = (lower.clone(), upper.clone());
        ceiling.0[j] = Some(up);
        let mut floor = (lower, upper);
        floor.1[j] = Some(down);
        stack.push(ceiling);
        stack.push(floor);
    }
    Ok(Ended::Searched)
}
