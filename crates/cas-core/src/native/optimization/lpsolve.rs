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
//!   - The keys of branch and bound (CAS-602, CAS-603; lpsolve.cc:2405-2443,
//!     :2459-2463, :2489-2493; [`super::integer`]): `lp_nodeselect = s`
//!     with `s` `lp_depthfirst`, `lp_breadthfirst`, `lp_bestlocalbound`,
//!     `lp_bestprojection` or `lp_hybrid`; `lp_varselect = s` with `s`
//!     `lp_firstfractional`, `lp_lastfractional`, `lp_mostfractional` or
//!     `lp_pseudocost`; `lp_heuristic = b`; `lp_gaptolerance = g`, `g` a
//!     nonnegative rational or a decimal made exact; and the limits
//!     `lp_maxcuts`, `lp_depthlimit`, `lp_nodelimit`, `lp_iterationlimit`
//!     and `lp_timelimit`, each a nonnegative integer or `Infinity`, which
//!     is 0, no limit (giac's `parse_limit`, :2350-2362; for `lp_maxcuts`, 0
//!     is no cut, giac's too). Ours has no clock: `lp_timelimit` is read
//!     and changes nothing; the evaluator's limits and its interrupt bound
//!     the work. `lp_iterationlimit` counts simplex iterations over the
//!     whole call, a program with no integer variable too
//!     ([`super::simplex::Work`]).
//!   - Any other option, a value of another kind, or a range of a symbol
//!     that is not a variable leaves the call.
//! - Each variable becomes nonnegative: `x = a + x′` for a lower end `a`,
//!   with `x′ ≤ b − a` for an upper one, `x = b − x′` for an upper end
//!   alone, `x = x⁺ − x⁻` for none. The program goes to [`simplex`], which
//!   minimizes, `−c` for a maximum.
//! - The relaxation of an integer program takes giac's Gomory
//!   mixed-integer cuts (lpsolve.cc:715-790), at most `lp_maxcuts` of them:
//!   from each row of the optimal tableau whose basic column is an integer
//!   one at a fraction, in row order, kept by giac's tests (a fraction at
//!   least 0.08 from an integer, parallel enough to the objective, 0.86,
//!   integralized coefficients at most 10^6, in doubles as giac's but the
//!   last), written over the program's columns, slacks substituted, as a
//!   `≥` row, and the program solved again from scratch. A column is an
//!   integer one when its variable is an integer and its shift end is
//!   (both of `x⁺ − x⁻`); slacks are continuous. giac adds the cut with its
//!   coefficients integralized and its right side kept 1, another
//!   inequality whenever that scale is not 1, which may cut off integer
//!   points; ours adds the cut itself. A cut with a part past
//!   [`MAX_BITS`] bits is dropped, and a relaxation's point so gives up.
//! - The answer is giac's `{optimum, {x = …, …}}`. An unbounded program is
//!   `{-Infinity, {}}`, or `{Infinity, {}}` for a maximum, and an
//!   infeasible one `{}`, each with a warning, as giac's. A search a user's
//!   limit or gap stops answers its best integer point, or `{}` (giac's
//!   `[]`) when it has none, with a warning that names the option.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear programs: c·x, rows a·x ≤ b, sides l and r"
)]

use cas_num::ball::{Float, Round};
use cas_num::{Integer, Rational};

use super::integer::{NodeSelect, Settings, VarSelect, branch_and_bound, fraction};
use super::simplex::{Final, Outcome, Relation, Row, Work, minimize};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The most bits of a numerator or denominator of a cut or a point.
const MAX_BITS: u64 = 1024;

/// giac's `LP_MIN_AWAY`, `LP_MIN_PARALLELISM` and `LP_MAX_MAGNITUDE`
/// (lpsolve.h:36-38).
const MIN_AWAY: f64 = 0.08;
const MIN_PARALLELISM: f64 = 0.86;
const MAX_MAGNITUDE: i64 = 1_000_000;

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
    pub(super) settings: Settings,
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
    /// Branch and bound's, `lp_gaptolerance` made exact later.
    settings: Settings,
    gap_tolerance: Option<&'a Expr>,
}

/// giac's limit (`parse_limit`, lpsolve.cc:2350-2362): a nonnegative
/// integer, `Infinity` 0, no limit; a huge one the most a `usize` holds.
fn limit(e: &Expr) -> Option<usize> {
    if e.as_symbol() == Some(Symbol::INFINITY) {
        return Some(0);
    }
    let q = e.to_rational()?;
    if !q.is_integer() || q.is_negative() {
        return None;
    }
    Some(usize::try_from(q.numer()).unwrap_or(usize::MAX))
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
            (Some(Symbol::LP_NODE_SELECT), _, _) => {
                o.settings.node_select = match r.as_symbol()? {
                    Symbol::LP_DEPTH_FIRST => NodeSelect::DepthFirst,
                    Symbol::LP_BREADTH_FIRST => NodeSelect::BreadthFirst,
                    Symbol::LP_BEST_LOCAL_BOUND => NodeSelect::BestLocalBound,
                    Symbol::LP_BEST_PROJECTION => NodeSelect::BestProjection,
                    Symbol::LP_HYBRID => NodeSelect::Hybrid,
                    _ => return None,
                };
            }
            (Some(Symbol::LP_VAR_SELECT), _, _) => {
                o.settings.var_select = match r.as_symbol()? {
                    Symbol::LP_FIRST_FRACTIONAL => VarSelect::FirstFractional,
                    Symbol::LP_LAST_FRACTIONAL => VarSelect::LastFractional,
                    Symbol::LP_MOST_FRACTIONAL => VarSelect::MostFractional,
                    Symbol::LP_PSEUDO_COST => VarSelect::PseudoCost,
                    _ => return None,
                };
            }
            (Some(Symbol::LP_HEURISTIC), _, _) => o.settings.heuristic = boolean(r)?,
            (Some(Symbol::LP_MAX_CUTS), _, _) => o.settings.max_cuts = limit(r)?,
            (Some(Symbol::LP_DEPTH_LIMIT), _, _) => o.settings.depth_limit = limit(r)?,
            (Some(Symbol::LP_NODE_LIMIT), _, _) => o.settings.node_limit = limit(r)?,
            (Some(Symbol::LP_ITERATION_LIMIT), _, _) => o.settings.iteration_limit = limit(r)?,
            // Ours has no clock: the evaluator's limits and its interrupt
            // bound the work.
            (Some(Symbol::LP_TIME_LIMIT), _, _) => {
                limit(r)?;
            }
            (Some(Symbol::LP_GAP_TOLERANCE), _, _) => o.gap_tolerance = Some(r),
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
    let mut settings = o.settings;
    if let Some(g) = o.gap_tolerance {
        match exact(ev, g)?.as_ref().and_then(Expr::to_rational) {
            Some(q) if !q.is_negative() => settings.gap_tolerance = q,
            _ => return Ok(None),
        }
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
        settings,
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

/// The nearest `f64`.
fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// Whether the numerator and denominator of `q` have at most [`MAX_BITS`]
/// bits.
fn small(q: &Rational) -> bool {
    q.numer().bit_length() <= MAX_BITS && q.denom().bit_length() <= MAX_BITS
}

/// giac's Gomory mixed-integer cuts (lpsolve.cc:715-773) of the optimal
/// tableau `last` of `rows`, its first `n` columns the program's, those
/// `integer` integer ones, `cost` the objective minimized: each over the
/// program's columns, the slacks substituted, in row order.
fn cuts(last: &Final, rows: &[Row], integer: &[bool], cost: &[Rational]) -> Vec<Row> {
    let n = integer.len();
    let real = last.real;
    let norm = cost.iter().map(|c| to_f64(c).powi(2)).sum::<f64>().sqrt();
    let mut found = Vec::new();
    for (row, &basic) in last.rows.iter().zip(&last.basis) {
        if basic >= n || !integer[basic] {
            continue;
        }
        let f0 = fraction(&row[row.len() - 1]);
        let g0 = &Rational::ONE - &f0;
        if f0.is_zero() || to_f64(&f0).min(to_f64(&g0)) < MIN_AWAY {
            continue;
        }
        // The cut Σ g_k x_k ≥ 1 over the columns before the artificials,
        // which are 0; a basic column's entry is 0 or 1, its g 0.
        let g: Vec<Rational> = row[..real]
            .iter()
            .enumerate()
            .map(|(k, a)| {
                if k < n && integer[k] {
                    let fj = fraction(a);
                    if fj > f0 {
                        &(&Rational::ONE - &fj) / &g0
                    } else {
                        &fj / &f0
                    }
                } else if a.is_negative() {
                    &(-a.clone()) / &g0
                } else {
                    a / &f0
                }
            })
            .collect();
        if g.iter().all(Rational::is_zero) {
            continue;
        }
        // Parallel enough to the objective, over the program's columns, in
        // doubles as giac's: a cut with no such part passes (NaN or ∞).
        let along = g[..n]
            .iter()
            .zip(cost)
            .fold(Rational::ZERO, |s, (a, c)| &s + &(a * c));
        let squares = g[..n].iter().fold(Rational::ZERO, |s, a| &s + &(a * a));
        if (to_f64(&along) / (norm * to_f64(&squares).sqrt())).abs() < MIN_PARALLELISM {
            continue;
        }
        // giac's `integralize`: times the denominators' lcm, over the
        // numerators' gcd, then at most LP_MAX_MAGNITUDE, exactly.
        let lcm = g
            .iter()
            .filter(|a| !a.is_zero())
            .fold(Integer::from(1), |l, a| l.lcm(a.denom()));
        let lcm = Rational::from(lcm);
        let scaled: Vec<Integer> = g
            .iter()
            .map(|a| Integer::try_from(&(a * &lcm)).expect("an integer"))
            .collect();
        let gcd = scaled.iter().fold(Integer::from(0), |d, a| d.gcd(a));
        let most = scaled
            .iter()
            .map(Integer::abs)
            .max()
            .expect("a coefficient");
        if &most / &gcd > Integer::from(MAX_MAGNITUDE) {
            continue;
        }
        // A slack of row (a, ≤, b) is b − a·x, a surplus of (a, ≥, b)
        // a·x − b.
        let mut a: Vec<Rational> = g[..n].to_vec();
        let mut b = Rational::ONE;
        for (k, gk) in g.iter().enumerate().take(real).skip(n) {
            if gk.is_zero() {
                continue;
            }
            let (ar, relation, br) = &rows[last.slack_rows[k - n]];
            let sigma = if *relation == Relation::AtMost {
                gk.clone()
            } else {
                -gk.clone()
            };
            for (aj, arj) in a.iter_mut().zip(ar) {
                *aj = &*aj - &(&sigma * arj);
            }
            b = &b - &(&sigma * br);
        }
        if a.iter().chain([&b]).all(small) {
            found.push((a, Relation::AtLeast, b));
        }
    }
    found
}

/// The integer columns of `p` in `bounds` shifted by `shifts`: an integer
/// variable's whose shift end is an integer, both of `x⁺ − x⁻`.
fn integer_columns(
    p: &Program,
    shifts: &[Shift],
    bounds: &[(Option<Rational>, Option<Rational>)],
    columns: usize,
) -> Vec<bool> {
    let mut integer = vec![false; columns];
    for ((shift, (low, high)), &is) in shifts.iter().zip(bounds).zip(&p.integer) {
        if !is {
            continue;
        }
        match *shift {
            Shift::From(c) => integer[c] = low.as_ref().is_some_and(Rational::is_integer),
            Shift::Below(c) => integer[c] = high.as_ref().is_some_and(Rational::is_integer),
            Shift::Free(c) => {
                integer[c] = true;
                integer[c + 1] = true;
            }
        }
    }
    integer
}

/// The continuous relaxation's best point and value, in the program's
/// variables, of `p` in `bounds`, with at most `max_cuts` of giac's cuts
/// when it has integer variables (lpsolve.cc:592-838): solved, cut, and
/// solved again from scratch until no cut is found or `max_cuts` were
/// added (`cuts_applied - cuts0 >= max_cuts`, :714).
pub(super) fn relaxation(
    p: &Program,
    bounds: &[(Option<Rational>, Option<Rational>)],
    work: &mut Work<'_>,
    max_cuts: usize,
) -> Outcome {
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
    let integer = integer_columns(p, &shifts, bounds, columns);
    let cut = p.integer.iter().any(|&b| b);
    let mut applied = 0;
    let y = loop {
        let (y, last) = match minimize(&cost, &rows, work) {
            (Outcome::Optimal(y, _), Some(last)) => (y, last),
            (outcome, _) => return outcome,
        };
        if !cut || applied >= max_cuts {
            break y;
        }
        let found = cuts(&last, &rows, &integer, &cost);
        if found.is_empty() {
            break y;
        }
        for row in found {
            rows.push(row);
            applied += 1;
            if applied == max_cuts {
                break;
            }
        }
    };
    let x: Vec<Rational> = shifts
        .iter()
        .zip(bounds)
        .map(|(shift, (low, high))| match *shift {
            Shift::From(c) => low.as_ref().expect("a lower end") + &y[c],
            Shift::Below(c) => high.as_ref().expect("an upper end") - &y[c],
            Shift::Free(c) => &y[c] - &y[c + 1],
        })
        .collect();
    if !x.iter().all(small) {
        return Outcome::GaveUp;
    }
    let value = p
        .objective
        .iter()
        .zip(&x)
        .fold(p.constant.clone(), |s, (a, v)| &s + &(a * v));
    Outcome::Optimal(x, value)
}

/// `LPSolve[…]`.
pub(in crate::native) fn lp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = program(ev, e)? else {
        return Ok(None);
    };
    let context = ev.context();
    let mut work = Work::new(
        p.settings.iteration_limit,
        context.limits().terms,
        context.interrupt_flag(),
    );
    let (outcome, stop) = if p.integer.iter().any(|&b| b) {
        branch_and_bound(&p, &mut work)
    } else {
        match relaxation(&p, &p.bounds, &mut work, 0) {
            Outcome::Limit => (Outcome::Infeasible, Some("lp_iterationlimit")),
            outcome => (outcome, None),
        }
    };
    let result = match outcome {
        Outcome::Optimal(x, value) => {
            if let Some(option) = stop {
                let message = format!("{option} reached; the best integer point found is answered");
                warn(ev, e, &message);
            }
            list(vec![
                Expr::from(value),
                list(
                    p.xs.iter()
                        .zip(x)
                        .map(|(v, q)| Expr::call(Symbol::EQUATION, [v.clone(), Expr::from(q)]))
                        .collect(),
                ),
            ])
        }
        Outcome::Unbounded => {
            warn(ev, e, "the program is unbounded");
            let infinity = Expr::from(Symbol::INFINITY);
            let value = if p.maximum { infinity } else { -infinity };
            list(vec![value, list(Vec::new())])
        }
        Outcome::Infeasible => {
            let message = match stop {
                Some(option) => format!("{option} reached before any point was found"),
                None => "the program has no feasible point".to_owned(),
            };
            warn(ev, e, &message);
            list(Vec::new())
        }
        Outcome::Interrupted => return Err(EvalError::Interrupted),
        Outcome::GaveUp | Outcome::Limit => {
            warn(ev, e, "the program takes too many steps");
            return Ok(None);
        }
    };
    ev.step("lpsolve", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
