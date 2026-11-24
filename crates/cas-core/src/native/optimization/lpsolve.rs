//! Linear programs (CAS-DOC-38 §4, CAS-417): `LPSolve`, Xcas's `lpsolve`.
//!
//! `LPSolve[c, constraints, bounds…, options…]`:
//!
//! - `c` and the constraints' sides are linear in the variables with
//!   rational coefficients; decimals in `c`, the constraints and the ends
//!   of the ranges are made exact first by `Exact`, as giac's
//!   `make_problem_exact` (lpsolve.cc:1306; CAS-601), and a ball that
//!   stays leaves the call. The constraints are `≤`, `≥`, `=`, and `<`, `>`
//!   taken as their closures, in a list; the variables are the symbols of
//!   `c` and the constraints in the order they first appear. Anything else
//!   leaves the call.
//! - The bounds are `x = a..b` and `[x, y, …] = a..b`, ends rational or
//!   infinite, of the program's variables; each intersects the variable's
//!   bounds, and a variable is free otherwise. The options, giac's
//!   (CAS-600, lpsolve.cc:2364-2504), in any order, a later one deciding:
//!   - `maximize` or `lp_maximize` alone, or `maximize = b`,
//!     `lp_maximize = b` with `b` `true` or `false`, maximize; `minimize`
//!     alone minimizes.
//!   - `assume = k` or `lp_assume = k`, the last one, applies to every
//!     variable after the other options: `integer` and `lp_integer`
//!     integers, `lp_nonnegint` (`nonnegint`) nonnegative integers,
//!     `lp_binary` integers in `0..1`, `lp_nonnegative` at least 0
//!     ([`super::integer::kind`]). Ours also takes `lp_nonnegative = true`
//!     (`false` changes nothing), which giac ignores.
//!   - `lp_integervariables = [x, …]` makes those variables integers,
//!     `lp_binaryvariables = [x, …]` integers in `0..1`; giac's indices
//!     and index ranges there leave the call.
//!   - `lp_method = m` with `m` `exact`, `float`, `lp_simplex` or
//!     `lp_interiorpoint` (CAS-601): every one is our exact simplex, and
//!     the answer is exact, as giac's answers are with each of them.
//!   - `lp_verbose` and `lp_presolve`, alone or `= true` or `= false`, do
//!     nothing: ours prints nothing more and has no presolve step.
//!   - Another of giac's words alone (`lp_integer`, `lp_simplex`, …) is
//!     ignored, as giac's, with a warning that names it.
//!   - The keys of branch and bound (`lp_nodeselect`, `lp_varselect`,
//!     `lp_maxcuts`, `lp_heuristic`, `lp_depthlimit`, `lp_nodelimit`,
//!     `lp_iterationlimit`, `lp_timelimit`, `lp_gaptolerance`) leave the
//!     call for now (T3), and so does any other option, a value of another
//!     kind, or a range of a symbol that is not a variable.
//! - Each variable becomes nonnegative: `x = a + x′` for a lower end `a`,
//!   with `x′ ≤ b − a` for an upper one, `x = b − x′` for an upper end
//!   alone, `x = x⁺ − x⁻` for none. The program goes to [`simplex`], which
//!   minimizes, `−c` for a maximum.
//! - The answer is giac's `{optimum, {x = …, …}}`. An unbounded program is
//!   `{-Infinity, {}}`, or `{Infinity, {}}` for a maximum, and an
//!   infeasible one `{}`, each with a warning, as giac's.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear programs: c·x, rows a·x ≤ b, sides l and r"
)]

use cas_num::Rational;

use super::simplex::{Outcome, Relation, Row, minimize};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

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

/// The symbols of `e` that are not builtins, in the order they first
/// appear, depth first and left to right, after those in `found`.
pub(super) fn symbols_in_order(e: &Expr, found: &mut Vec<Expr>) {
    if let Some(s) = e.as_symbol() {
        if !s.is_builtin() && !found.contains(e) {
            found.push(e.clone());
        }
        return;
    }
    for a in e.args() {
        symbols_in_order(a, found);
    }
}

/// `e` as `Σ a_j x_j + k`, rational coefficients, evaluated; `None` when
/// it is not linear so.
fn linear(
    ev: &mut Eval<'_>,
    e: &Expr,
    xs: &[Expr],
) -> Result<Option<(Vec<Rational>, Rational)>, EvalError> {
    let mut coefficients = Vec::with_capacity(xs.len());
    for x in xs {
        let a = ev.eval(&Expr::call(Symbol::D, [e.clone(), x.clone()]))?;
        let Some(a) = a.to_rational() else {
            return Ok(None);
        };
        coefficients.push(a);
    }
    let zero: Vec<(Expr, Expr)> = xs.iter().map(|x| (x.clone(), Expr::from(0))).collect();
    let Some(k) = ev.eval(&e.substitute(&zero))?.to_rational() else {
        return Ok(None);
    };
    let mut rebuilt = Expr::from(k.clone());
    for (a, x) in coefficients.iter().zip(xs) {
        rebuilt += Expr::from(a.clone()) * x.clone();
    }
    let gap = ev.eval(&Expr::call(Symbol::EXPAND, [e.clone() - rebuilt]))?;
    Ok((gap == Expr::from(0)).then_some((coefficients, k)))
}

/// A rational end of a range, `None` for an infinite one; `Err(())` for
/// anything else.
fn end(e: &Expr) -> Result<Option<Rational>, ()> {
    if let Some(q) = e.to_rational() {
        return Ok(Some(q));
    }
    let infinite = e.may_contain(Symbol::INFINITY) || e.may_contain(Symbol::DIRECTED_INFINITY);
    if infinite { Ok(None) } else { Err(()) }
}

/// Whether `e` holds a ball, real or complex.
fn holds_ball(e: &Expr) -> bool {
    matches!(e.view(), View::Ball(_) | View::ComplexBall(_)) || e.args().iter().any(holds_ball)
}

/// `e` with its decimals exact, `Exact[e]` as giac's `make_problem_exact`
/// (lpsolve.cc:1306); `None` when a ball stays.
fn exact(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    if !holds_ball(e) {
        return Ok(Some(e.clone()));
    }
    let e = ev.eval(&Expr::call(Symbol::EXACT, [e.clone()]))?;
    Ok((!holds_ball(&e)).then_some(e))
}

/// The row of the constraint `c`, its decimals exact first; `None` when it
/// is not a linear one.
fn row(ev: &mut Eval<'_>, c: &Expr, xs: &[Expr]) -> Result<Option<Row>, EvalError> {
    let Some(c) = exact(ev, c)? else {
        return Ok(None);
    };
    let (Some(h), [l, r]) = (head(&c), c.args()) else {
        return Ok(None);
    };
    let relation = match h {
        Symbol::LESS_EQUAL | Symbol::LESS => Relation::AtMost,
        Symbol::GREATER_EQUAL | Symbol::GREATER => Relation::AtLeast,
        Symbol::EQUATION | Symbol::EQUAL => Relation::Equal,
        _ => return Ok(None),
    };
    let Some((a, k)) = linear(ev, &(l.clone() - r.clone()), xs)? else {
        return Ok(None);
    };
    Ok(Some((a, relation, -k)))
}

/// Narrows the ends `b` of a variable to their intersection with
/// `low..high`, an absent end infinite.
pub(super) fn tighten(
    b: &mut (Option<Rational>, Option<Rational>),
    low: Option<&Rational>,
    high: Option<&Rational>,
) {
    if let Some(low) = low
        && b.0.as_ref().is_none_or(|l| l < low)
    {
        b.0 = Some(low.clone());
    }
    if let Some(high) = high
        && b.1.as_ref().is_none_or(|h| h > high)
    {
        b.1 = Some(high.clone());
    }
}

/// `True` or `False`.
fn boolean(e: &Expr) -> Option<bool> {
    match e.as_symbol() {
        Some(Symbol::TRUE) => Some(true),
        Some(Symbol::FALSE) => Some(false),
        _ => None,
    }
}

/// giac's spelling of one of `lpsolve`'s option words; `None` for another
/// symbol.
fn lp_word(s: Symbol) -> Option<&'static str> {
    Some(match s {
        Symbol::LP_ASSUME => "lp_assume",
        Symbol::LP_MAXIMIZE => "lp_maximize",
        Symbol::LP_NON_NEGATIVE => "lp_nonnegative",
        Symbol::LP_NON_NEG_INT => "lp_nonnegint",
        Symbol::LP_INTEGER => "lp_integer",
        Symbol::LP_BINARY => "lp_binary",
        Symbol::LP_INTEGER_VARIABLES => "lp_integervariables",
        Symbol::LP_BINARY_VARIABLES => "lp_binaryvariables",
        Symbol::LP_METHOD => "lp_method",
        Symbol::LP_SIMPLEX => "lp_simplex",
        Symbol::LP_INTERIOR_POINT => "lp_interiorpoint",
        Symbol::LP_VERBOSE => "lp_verbose",
        Symbol::LP_PRESOLVE => "lp_presolve",
        Symbol::LP_NODE_SELECT => "lp_nodeselect",
        Symbol::LP_DEPTH_FIRST => "lp_depthfirst",
        Symbol::LP_BREADTH_FIRST => "lp_breadthfirst",
        Symbol::LP_BEST_LOCAL_BOUND => "lp_bestlocalbound",
        Symbol::LP_BEST_PROJECTION => "lp_bestprojection",
        Symbol::LP_HYBRID => "lp_hybrid",
        Symbol::LP_VAR_SELECT => "lp_varselect",
        Symbol::LP_FIRST_FRACTIONAL => "lp_firstfractional",
        Symbol::LP_LAST_FRACTIONAL => "lp_lastfractional",
        Symbol::LP_MOST_FRACTIONAL => "lp_mostfractional",
        Symbol::LP_PSEUDO_COST => "lp_pseudocost",
        Symbol::LP_MAX_CUTS => "lp_maxcuts",
        Symbol::LP_HEURISTIC => "lp_heuristic",
        Symbol::LP_DEPTH_LIMIT => "lp_depthlimit",
        Symbol::LP_NODE_LIMIT => "lp_nodelimit",
        Symbol::LP_ITERATION_LIMIT => "lp_iterationlimit",
        Symbol::LP_TIME_LIMIT => "lp_timelimit",
        Symbol::LP_GAP_TOLERANCE => "lp_gaptolerance",
        _ => return None,
    })
}

/// A program: the objective's coefficients and constant, the rows, the
/// variables' bounds, and whether it maximizes.
pub(super) struct Program {
    pub(super) xs: Vec<Expr>,
    pub(super) objective: Vec<Rational>,
    pub(super) constant: Rational,
    pub(super) rows: Vec<Row>,
    pub(super) bounds: Vec<(Option<Rational>, Option<Rational>)>,
    pub(super) maximum: bool,
    pub(super) integer: Vec<bool>,
}

/// The options of a call, as giac reads them; ranges and lists of
/// variables are checked against the program's later.
#[derive(Default)]
struct Options<'a> {
    maximum: bool,
    nonnegative: bool,
    assumption: Option<&'a Expr>,
    /// The lists of `lp_integervariables` (false) and `lp_binaryvariables`
    /// (true).
    typed: Vec<(&'a Expr, bool)>,
    ranges: Vec<(&'a Expr, &'a Expr, &'a Expr)>,
    /// The words alone that do nothing, giac's spelling.
    ignored: Vec<&'static str>,
}

/// The options `rest`; `None` when one is not an option.
fn options(rest: &[Expr]) -> Option<Options<'_>> {
    let mut o = Options::default();
    for option in rest {
        if let Some(s) = option.as_symbol() {
            match s {
                Symbol::MAXIMIZE | Symbol::LP_MAXIMIZE => o.maximum = true,
                Symbol::MINIMIZE => o.maximum = false,
                // Ours prints nothing more and has no presolve step.
                Symbol::LP_VERBOSE | Symbol::LP_PRESOLVE => {}
                _ => o.ignored.push(lp_word(s)?),
            }
            continue;
        }
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) = (head(option), option.args()) else {
            return None;
        };
        match (l.as_symbol(), head(r), r.args()) {
            (Some(Symbol::MAXIMIZE | Symbol::LP_MAXIMIZE), _, _) => o.maximum = boolean(r)?,
            (Some(Symbol::ASSUME | Symbol::LP_ASSUME), _, _) => o.assumption = Some(r),
            (Some(Symbol::LP_NON_NEGATIVE), _, _) => o.nonnegative |= boolean(r)?,
            (Some(Symbol::LP_INTEGER_VARIABLES), _, _) => o.typed.push((r, false)),
            (Some(Symbol::LP_BINARY_VARIABLES), _, _) => o.typed.push((r, true)),
            // Every method is our exact simplex.
            (Some(Symbol::LP_METHOD), _, _) => {
                let method = matches!(
                    r.as_symbol(),
                    Some(Symbol::EXACT | Symbol::LP_SIMPLEX | Symbol::LP_INTERIOR_POINT)
                ) || *r == Expr::symbol("float");
                if !method {
                    return None;
                }
            }
            (Some(Symbol::LP_VERBOSE | Symbol::LP_PRESOLVE), _, _) => {
                boolean(r)?;
            }
            (_, Some(Symbol::RANGE), [a, b]) => o.ranges.push((l, a, b)),
            _ => return None,
        }
    }
    Some(o)
}

/// Warns of each word in `ignored`.
fn warn_ignored(ev: &mut Eval<'_>, e: &Expr, ignored: &[&str]) {
    for word in ignored {
        let message = if matches!(
            *word,
            "lp_integer" | "lp_binary" | "lp_nonnegative" | "lp_nonnegint"
        ) {
            format!("the option {word} alone is ignored; assume={word} applies it")
        } else {
            format!("the option {word} alone is ignored")
        };
        warn(ev, e, &message);
    }
}

/// The program of a call; `None` when it is not one.
pub(super) fn program(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Program>, EvalError> {
    let [objective, rest @ ..] = e.args() else {
        return Ok(None);
    };
    let (constraints, rest) = match rest {
        [c, more @ ..] if head(c) == Some(Symbol::LIST) => (c.args(), more),
        _ => (&[][..], rest),
    };
    let Some(o) = options(rest) else {
        return Ok(None);
    };
    let mut xs = Vec::new();
    symbols_in_order(objective, &mut xs);
    for c in constraints {
        symbols_in_order(c, &mut xs);
    }
    // The position of a program variable, `None` for anything else.
    let index = |x: &Expr| xs.iter().position(|v| v == x);
    let mut spans: Vec<(usize, &Expr, &Expr)> = Vec::new();
    for (l, a, b) in o.ranges {
        let vars = if head(l) == Some(Symbol::LIST) {
            l.args()
        } else {
            std::slice::from_ref(l)
        };
        for x in vars {
            let Some(j) = index(x) else {
                return Ok(None);
            };
            spans.push((j, a, b));
        }
    }
    let mut integer = vec![false; xs.len()];
    let mut binary = vec![false; xs.len()];
    for (list, is_binary) in o.typed {
        if head(list) != Some(Symbol::LIST) {
            return Ok(None);
        }
        for x in list.args() {
            let Some(j) = index(x) else {
                return Ok(None);
            };
            integer[j] = true;
            binary[j] |= is_binary;
        }
    }
    let Some(objective) = exact(ev, objective)? else {
        return Ok(None);
    };
    let Some((objective, constant)) = linear(ev, &objective, &xs)? else {
        return Ok(None);
    };
    let mut rows = Vec::new();
    for c in constraints {
        let Some(row) = row(ev, c, &xs)? else {
            return Ok(None);
        };
        rows.push(row);
    }
    // Every option narrows the bounds; the assumption comes last, as
    // giac's (lpsolve.cc:2672-2689).
    let mut bounds = vec![(None, None); xs.len()];
    for (j, a, b) in spans {
        let (Some(a), Some(b)) = (exact(ev, a)?, exact(ev, b)?) else {
            return Ok(None);
        };
        let (Ok(low), Ok(high)) = (end(&a), end(&b)) else {
            return Ok(None);
        };
        tighten(&mut bounds[j], low.as_ref(), high.as_ref());
    }
    for (b, _) in bounds.iter_mut().zip(&binary).filter(|(_, is)| **is) {
        tighten(b, Some(&Rational::ZERO), Some(&Rational::ONE));
    }
    if o.nonnegative {
        for b in &mut bounds {
            tighten(b, Some(&Rational::ZERO), None);
        }
    }
    if let Some(k) = o.assumption
        && !super::integer::kind(k, &mut bounds, &mut integer)
    {
        return Ok(None);
    }
    warn_ignored(ev, e, &o.ignored);
    Ok(Some(Program {
        xs,
        objective,
        constant,
        rows,
        bounds,
        maximum: o.maximum,
        integer,
    }))
}

/// How a variable goes to nonnegative ones: `x = a + x′`, `x = b − x′`, or
/// `x = x⁺ − x⁻`, by the column of `x′` or `x⁺`.
#[derive(Clone, Copy)]
enum Shift {
    From(usize),
    Below(usize),
    Free(usize),
}

/// The continuous relaxation's best point and value, in the program's
/// variables.
pub(super) fn relaxation(p: &Program, bounds: &[(Option<Rational>, Option<Rational>)]) -> Outcome {
    let mut shifts = Vec::with_capacity(p.xs.len());
    let mut columns = 0;
    for (low, high) in bounds {
        shifts.push(match (low, high) {
            (Some(_), _) => Shift::From(columns),
            (None, Some(_)) => Shift::Below(columns),
            (None, None) => Shift::Free(columns),
        });
        columns += if matches!((low, high), (None, None)) {
            2
        } else {
            1
        };
    }
    // a·x in the new columns, and the constant it adds.
    let transform = |a: &[Rational]| -> (Vec<Rational>, Rational) {
        let mut row = vec![Rational::ZERO; columns];
        let mut constant = Rational::ZERO;
        for ((coefficient, shift), (low, high)) in a.iter().zip(&shifts).zip(bounds) {
            match *shift {
                Shift::From(c) => {
                    row[c] = coefficient.clone();
                    constant = &constant + &(coefficient * low.as_ref().expect("a lower end"));
                }
                Shift::Below(c) => {
                    row[c] = -coefficient.clone();
                    constant = &constant + &(coefficient * high.as_ref().expect("an upper end"));
                }
                Shift::Free(c) => {
                    row[c] = coefficient.clone();
                    row[c + 1] = -coefficient.clone();
                }
            }
        }
        (row, constant)
    };
    let sign = if p.maximum {
        -Rational::ONE
    } else {
        Rational::ONE
    };
    let objective: Vec<Rational> = p.objective.iter().map(|a| a * &sign).collect();
    let (cost, _) = transform(&objective);
    let mut rows: Vec<Row> = Vec::new();
    for (a, relation, b) in &p.rows {
        let (row, constant) = transform(a);
        rows.push((row, *relation, b - &constant));
    }
    for (shift, (low, high)) in shifts.iter().zip(bounds) {
        if let (Shift::From(c), Some(low), Some(high)) = (*shift, low, high) {
            let mut row = vec![Rational::ZERO; columns];
            row[c] = Rational::ONE;
            rows.push((row, Relation::AtMost, high - low));
        }
    }
    match minimize(&cost, &rows) {
        Outcome::Optimal(y, _) => {
            let x: Vec<Rational> = shifts
                .iter()
                .zip(bounds)
                .map(|(shift, (low, high))| match *shift {
                    Shift::From(c) => low.as_ref().expect("a lower end") + &y[c],
                    Shift::Below(c) => high.as_ref().expect("an upper end") - &y[c],
                    Shift::Free(c) => &y[c] - &y[c + 1],
                })
                .collect();
            let value = p
                .objective
                .iter()
                .zip(&x)
                .fold(p.constant.clone(), |s, (a, v)| &s + &(a * v));
            Outcome::Optimal(x, value)
        }
        other => other,
    }
}

/// `LPSolve[…]`.
pub(in crate::native) fn lp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = program(ev, e)? else {
        return Ok(None);
    };
    let outcome = if p.integer.iter().any(|&b| b) {
        super::integer::branch_and_bound(&p)
    } else {
        relaxation(&p, &p.bounds)
    };
    let result = match outcome {
        Outcome::Optimal(x, value) => list(vec![
            Expr::from(value),
            list(
                p.xs.iter()
                    .zip(x)
                    .map(|(v, q)| Expr::call(Symbol::EQUATION, [v.clone(), Expr::from(q)]))
                    .collect(),
            ),
        ]),
        Outcome::Unbounded => {
            warn(ev, e, "the program is unbounded");
            let infinity = Expr::from(Symbol::INFINITY);
            let value = if p.maximum { infinity } else { -infinity };
            list(vec![value, list(Vec::new())])
        }
        Outcome::Infeasible => {
            warn(ev, e, "the program has no feasible point");
            list(Vec::new())
        }
        Outcome::GaveUp => {
            warn(ev, e, "the program takes too many steps");
            return Ok(None);
        }
    };
    ev.step("lpsolve", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
