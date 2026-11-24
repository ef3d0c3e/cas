//! Nonlinear programs (CAS-519, CAS-604): `NLPSolve`, Xcas's `nlpsolve`, a
//! local optimum of a smooth function under equations, inequalities and
//! bounds.
//!
//! `NLPSolve[f, constraints, bounds…, options…]`:
//!
//! - The constraints are `g = c`, `g ≤ c`, `g ≥ c`, and `<`, `>` taken as
//!   their closures, in a list that may hold bounds too; a bound is `x =
//!   a..b` or `[x, y, …] = a..b`, ends real or infinite, several of one
//!   variable intersected. The variables are the symbols of `f`, the
//!   constraints and the bounds, sorted, giac's order.
//! - The options, giac's (optimization.cc:6650-6832), matched by their
//!   builtins, a later one deciding:
//!   - `maximize` or `nlp_maximize` alone, or `maximize = b`,
//!     `nlp_maximize = b` with `b` `true` or `false`, maximize.
//!   - `assume = nlp_nonnegative` makes every variable at least 0 unless
//!     bounded higher; ours also takes `nlp_nonnegative = true`, which giac
//!     rejects. `assume = integer` or `nlp_integer` makes every variable an
//!     integer, `assume = nlp_nonnegint` (`nonnegint`) a nonnegative one,
//!     `assume = nlp_binary` one in `0..1`; `nlp_integervariables = [x, …]`
//!     and `nlp_binaryvariables = [x, …]` those variables. A nonlinear
//!     program with an integer variable leaves the call.
//!   - `nlp_initialpoint = [x = v, …]` starts the search there, `[x = a..b,
//!     …]` at the region's midpoint; giac's several points (`x = [v, …]`)
//!     leave the call.
//!   - `nlp_iterationlimit = n`, a positive integer, caps the outer
//!     iterations at `n` (at most [`MAX_OUTER`]); a search it stops before
//!     it converges leaves the call with a warning.
//!   - `nlp_precision = r`, a positive real: the certified balls' radii
//!     are at most `r` times their midpoints' magnitudes, or `r` below
//!     magnitude 1 (a zero optimum is then reachable), the certification
//!     raising its precision for it up to [`MAX_BITS`] bits, else the call
//!     stays with a warning.
//!   - `nlp_tolerance = r`, a positive real, is accepted and changes
//!     nothing: ours checks the constraints at the certified point.
//!   - `nlp_method = m`, `m` one of giac's strings `"automatic"`
//!     (`"auto"`), `"ipt"` (`"interior-point"`), `"nm"`
//!     (`"nelder-mead"`), `"de"` (`"differential-evolution"`), `"cobyla"`
//!     (`"cbl"`), or a list `["name", "parameter" = v, …]`, the
//!     parameters unchecked beyond their form: every method is our one
//!     certified method below. Another string, or a symbol, leaves the
//!     call with a warning that names it.
//!   - `nlp_presolve` and `nlp_verbose`, alone (giac rejects these) or `=
//!     b`, `nlp_verbose = n` with `n` in `0..4`, change nothing.
//!   - Anything else, or an option with a value of another kind, leaves
//!     the call; an equation whose left side is another builtin is never a
//!     constraint.
//! - A linear objective under linear constraints, rational coefficients
//!   once decimals are exact, is giac's linear program
//!   (optimization.cc:6080-6085): the call becomes `LPSolve[…]` with the
//!   same constraints, bounds, `maximize`, `assume = lp_nonnegative`,
//!   `lp_integervariables` and `lp_binaryvariables`, and its exact answer
//!   is `NLPSolve`'s, the point in this order.
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
use super::super::solve::fsolve::systems::certified_root_at;
use super::lpsolve::{exact, linear, lp_solve, symbols_in_order};
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

/// The most bits the certification raises its precision to for
/// `nlp_precision`.
const MAX_BITS: u64 = 4096;

/// giac's names of its methods (optimization.cc:6785-6797, optimization.h:
/// 109-113).
const METHODS: [&str; 10] = [
    "automatic",
    "auto",
    "ipt",
    "interior-point",
    "nm",
    "nelder-mead",
    "de",
    "differential-evolution",
    "cobyla",
    "cbl",
];

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

/// Each variable's lower or upper end, `None` for an infinite one.
type Ends = Vec<Option<End>>;

/// A bound's end: exact, with its double.
#[derive(Clone)]
struct End {
    exact: Expr,
    value: f64,
}

/// A program as read: the objective (negated for a maximum), the
/// equations `h = 0` and inequalities `g ≤ 0`, the bounds, the start, the
/// outer iterations' cap and the precision asked.
struct Program {
    xs: Vec<Expr>,
    objective: Expr,
    minimized: Expr,
    equalities: Vec<Expr>,
    inequalities: Vec<Expr>,
    lower: Vec<Option<End>>,
    upper: Vec<Option<End>>,
    start: Vec<Option<Expr>>,
    limit: usize,
    precision: Option<f64>,
}

/// A call as read.
enum Read {
    /// A nonlinear program.
    Nonlinear(Program),
    /// A linear one, as `LPSolve[…]`, and the variables in order.
    Linear(Expr, Vec<Expr>),
    /// None: the call stays, with this warning or none.
    Not(Option<String>),
}

/// `NLPSolve[f, constraints, bounds…, options…]`.
pub(in crate::native) fn nlp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let program = match program(ev, e)? {
        Read::Nonlinear(p) => p,
        Read::Linear(call, xs) => return linear_program(ev, &call, &xs),
        Read::Not(why) => {
            if let Some(message) = why {
                warn(ev, e, &message);
            }
            return Ok(None);
        }
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
    let search = match lagrangian(&compiled, x0, &lo, &hi, program.limit) {
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
        Err(Failure::Unconverged) => {
            warn(
                ev,
                e,
                "nlp_iterationlimit reached before the search converged",
            );
            return Ok(None);
        }
    };
    let cap = MAX_BITS.min(ev.context().limits().number_bits);
    let mut bits = ev.context().working_precision();
    let answer = loop {
        let Some(answer) = certify(ev, &program, &compiled, &search, &lo, &hi, bits)? else {
            warn(
                ev,
                e,
                "the search's point could not be certified as an optimum",
            );
            return Ok(None);
        };
        match program.precision {
            Some(r) if !within(&answer, r) => {}
            _ => break answer,
        }
        bits = bits.saturating_mul(2);
        if bits > cap {
            warn(ev, e, "nlp_precision could not be met");
            return Ok(None);
        }
    };
    ev.step("NLPSolve", e, &answer);
    Ok(Some(answer))
}

/// Whether every ball of `answer` has a radius at most `r` times its
/// midpoint's magnitude, or `r` below magnitude 1.
fn within(answer: &Expr, r: f64) -> bool {
    match answer.view() {
        View::Ball(b) => {
            let radius = b.rad().to_float().map_or(f64::INFINITY, |f| f.to_f64());
            radius <= r * b.mid().to_f64().abs().max(1.0)
        }
        _ => answer.args().iter().all(|a| within(a, r)),
    }
}

/// A linear program's answer, `LPSolve`'s, its point in the order `xs`.
fn linear_program(ev: &mut Eval<'_>, call: &Expr, xs: &[Expr]) -> Native {
    let Some(answer) = lp_solve(ev, call)? else {
        return Ok(None);
    };
    let (Some(Symbol::LIST), [value, point]) = (head(&answer), answer.args()) else {
        return Ok(Some(answer));
    };
    let mut point = point.args().to_vec();
    point.sort_by_key(|p| {
        p.args()
            .first()
            .and_then(|x| xs.iter().position(|y| y == x))
    });
    Ok(Some(list(vec![value.clone(), list(point)])))
}

/// The arguments after the objective, sorted out.
#[derive(Default)]
struct Items<'a> {
    maximum: bool,
    nonnegative: bool,
    /// Every variable an integer (`assume = integer`, `nlp_integer` or
    /// `nlp_nonnegint`): `Some(false)`; in `0..1` (`nlp_binary`):
    /// `Some(true)`.
    every: Option<bool>,
    /// The lists of `nlp_integervariables` (false) and
    /// `nlp_binaryvariables` (true).
    typed: Vec<(&'a Expr, bool)>,
    initial: &'a [Expr],
    iterations: Option<usize>,
    precision: Option<f64>,
    relations: Vec<(Symbol, &'a Expr, &'a Expr)>,
    /// One range per variable.
    ranges: Vec<(&'a Expr, &'a Expr, &'a Expr)>,
}

/// Why a call is no program: `None` silently, else a warning.
type Reject = Option<String>;

/// `o`'s value, or a silent rejection.
fn silent<T>(o: Option<T>) -> Result<T, Reject> {
    o.ok_or(None)
}

/// A `true` or `false` option's value.
fn flag(r: &Expr) -> Option<bool> {
    match r.as_symbol() {
        Some(Symbol::TRUE) => Some(true),
        Some(Symbol::FALSE) => Some(false),
        _ => None,
    }
}

/// A positive real option's value, as a double.
fn positive(ctx: &Context, r: &Expr) -> Option<f64> {
    double(ctx, r).filter(|v| *v > 0.0)
}

/// The constraints, bounds and options.
fn items<'a>(ctx: &Context, rest: &'a [Expr]) -> Result<Items<'a>, Reject> {
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
        if let Some(s) = item.as_symbol() {
            match s {
                Symbol::MAXIMIZE | Symbol::NLP_MAXIMIZE => found.maximum = true,
                // Ours prints nothing more and has no presolve step.
                Symbol::NLP_PRESOLVE | Symbol::NLP_VERBOSE => {}
                _ => return Err(None),
            }
            continue;
        }
        let (Some(h), [l, r]) = (head(item), item.args()) else {
            return Err(None);
        };
        match h {
            Symbol::EQUATION | Symbol::EQUAL => equation(ctx, &mut found, l, r)?,
            Symbol::LESS | Symbol::LESS_EQUAL | Symbol::GREATER | Symbol::GREATER_EQUAL => {
                found.relations.push((h, l, r));
            }
            _ => return Err(None),
        }
    }
    Ok(found)
}

/// One `l = r`: an option when `l` is a builtin, else a range of variables
/// or an equation.
fn equation<'a>(
    ctx: &Context,
    found: &mut Items<'a>,
    l: &'a Expr,
    r: &'a Expr,
) -> Result<(), Reject> {
    let Some(key) = l.as_symbol().filter(|s| s.is_builtin()) else {
        if let (Some(Symbol::RANGE), [a, b]) = (head(r), r.args()) {
            let xs = if head(l) == Some(Symbol::LIST) {
                l.args()
            } else {
                std::slice::from_ref(l)
            };
            for x in xs {
                if x.as_symbol().is_none_or(Symbol::is_builtin) {
                    return Err(None);
                }
                found.ranges.push((x, a, b));
            }
        } else if head(l) == Some(Symbol::LIST) {
            return Err(None);
        } else {
            found.relations.push((Symbol::EQUATION, l, r));
        }
        return Ok(());
    };
    match key {
        Symbol::MAXIMIZE | Symbol::NLP_MAXIMIZE => found.maximum = silent(flag(r))?,
        Symbol::ASSUME => match r.as_symbol() {
            Some(Symbol::NLP_NON_NEGATIVE) => found.nonnegative = true,
            Some(Symbol::NLP_INTEGER) => _ = found.every.get_or_insert(false),
            Some(Symbol::NLP_NON_NEG_INT | Symbol::LP_NON_NEG_INT) => {
                found.every.get_or_insert(false);
                found.nonnegative = true;
            }
            Some(Symbol::NLP_BINARY) => found.every = Some(true),
            // `integer` is a name of [types], not a builtin.
            _ if *r == Expr::symbol("integer") => _ = found.every.get_or_insert(false),
            _ => return Err(None),
        },
        Symbol::NLP_NON_NEGATIVE => found.nonnegative |= silent(flag(r))?,
        Symbol::NLP_INITIAL_POINT => {
            if head(r) != Some(Symbol::LIST) {
                return Err(None);
            }
            found.initial = r.args();
        }
        Symbol::NLP_ITERATION_LIMIT => {
            let q = silent(r.to_rational())?;
            if !q.is_integer() || !q.is_positive() {
                return Err(None);
            }
            found.iterations = Some(usize::try_from(q.numer()).unwrap_or(usize::MAX));
        }
        Symbol::NLP_PRECISION => found.precision = Some(silent(positive(ctx, r))?),
        // Ours checks the constraints at the certified point.
        Symbol::NLP_TOLERANCE => {
            silent(positive(ctx, r))?;
        }
        Symbol::NLP_METHOD => method(r)?,
        // Ours prints nothing more and has no presolve step.
        Symbol::NLP_PRESOLVE => {
            silent(flag(r))?;
        }
        Symbol::NLP_VERBOSE => {
            if flag(r).is_none() && !(0..=4).any(|n| *r == Expr::from(n)) {
                return Err(None);
            }
        }
        Symbol::NLP_INTEGER_VARIABLES => found.typed.push((r, false)),
        Symbol::NLP_BINARY_VARIABLES => found.typed.push((r, true)),
        _ => return Err(None),
    }
    Ok(())
}

/// `nlp_method`'s value: one of giac's names, alone or first in a list of
/// `"parameter" = v`; every one runs ours.
fn method(r: &Expr) -> Result<(), Reject> {
    let name = match (head(r), r.args()) {
        (Some(Symbol::LIST), [name, parameters @ ..]) => {
            let shaped = parameters.iter().all(|p| {
                matches!(
                    (head(p), p.args()),
                    (Some(Symbol::EQUATION | Symbol::EQUAL), [k, _]) if k.as_str().is_some()
                )
            });
            if !shaped {
                return Err(None);
            }
            name
        }
        _ => r,
    };
    match name.as_str() {
        Some(s) if METHODS.contains(&s) => Ok(()),
        Some(s) => Err(Some(format!(
            "nlp_method = \"{s}\" is not one of giac's methods"
        ))),
        None if name.as_symbol().is_some() => Err(Some(format!(
            "nlp_method = {name} is no method: giac's are strings, as \"{name}\""
        ))),
        None => Err(None),
    }
}

/// Whether the objective and every relation are linear in `xs`, rational
/// coefficients once decimals are exact, as `LPSolve` reads them.
fn is_linear(
    ev: &mut Eval<'_>,
    objective: &Expr,
    relations: &[(Symbol, &Expr, &Expr)],
    xs: &[Expr],
) -> Result<bool, EvalError> {
    let sides = relations
        .iter()
        .map(|(_, l, r)| (*l).clone() - (*r).clone());
    for f in std::iter::once(objective.clone()).chain(sides) {
        let Some(f) = exact(ev, &f)? else {
            return Ok(false);
        };
        if linear(ev, &f, xs)?.is_none() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// A linear program as `LPSolve[…]` (giac's, optimization.cc:6080-6085):
/// the constraints, bounds, maximum, nonnegativity, integer and binary
/// variables.
fn lp_call(
    objective: &Expr,
    relations: &[(Symbol, &Expr, &Expr)],
    ranges: &[(&Expr, &Expr, &Expr)],
    maximum: bool,
    nonnegative: bool,
    integers: Vec<Expr>,
    binaries: Vec<Expr>,
) -> Expr {
    let constraints = relations
        .iter()
        .map(|&(h, l, r)| Expr::call(h, [l.clone(), r.clone()]))
        .collect();
    let mut args = vec![objective.clone(), list(constraints)];
    for &(x, a, b) in ranges {
        let range = Expr::call(Symbol::RANGE, [a.clone(), b.clone()]);
        args.push(Expr::call(Symbol::EQUATION, [x.clone(), range]));
    }
    if maximum {
        args.push(Expr::from(Symbol::MAXIMIZE));
    }
    if nonnegative {
        let kind = Expr::from(Symbol::LP_NON_NEGATIVE);
        args.push(Expr::call(
            Symbol::EQUATION,
            [Expr::from(Symbol::ASSUME), kind],
        ));
    }
    for (xs, key) in [
        (integers, Symbol::LP_INTEGER_VARIABLES),
        (binaries, Symbol::LP_BINARY_VARIABLES),
    ] {
        if !xs.is_empty() {
            args.push(Expr::call(Symbol::EQUATION, [Expr::from(key), list(xs)]));
        }
    }
    Expr::call(Symbol::LP_SOLVE, args)
}

/// Which variables are integers, and which of them in `0..1`: all by
/// `every`, and those `typed` lists name; `None` for a list naming no
/// variable.
fn integers(
    xs: &[Expr],
    every: Option<bool>,
    typed: Vec<(&Expr, bool)>,
) -> Option<(Vec<bool>, Vec<bool>)> {
    let mut whole = vec![every.is_some(); xs.len()];
    let mut zero_one = vec![every == Some(true); xs.len()];
    for (named, is_binary) in typed {
        if head(named) != Some(Symbol::LIST) {
            return None;
        }
        for x in named.args() {
            let j = xs.iter().position(|v| v == x)?;
            whole[j] = true;
            zero_one[j] |= is_binary;
        }
    }
    Some((whole, zero_one))
}

/// The program of a call.
fn program(ev: &mut Eval<'_>, e: &Expr) -> Result<Read, EvalError> {
    let [objective, rest @ ..] = e.args() else {
        return Ok(Read::Not(None));
    };
    let Items {
        maximum,
        nonnegative,
        every,
        typed,
        initial,
        iterations,
        precision,
        relations,
        ranges,
    } = match items(ev.context(), rest) {
        Ok(found) => found,
        Err(why) => return Ok(Read::Not(why)),
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
        return Ok(Read::Not(None));
    }
    let Some((whole, zero_one)) = integers(&xs, every, typed) else {
        return Ok(Read::Not(None));
    };
    if is_linear(ev, objective, &relations, &xs)? {
        let pick = |integers: bool| -> Vec<Expr> {
            xs.iter()
                .zip(whole.iter().zip(&zero_one))
                .filter(|(_, (w, b))| **w && **b != integers)
                .map(|(x, _)| x.clone())
                .collect()
        };
        let call = lp_call(
            objective,
            &relations,
            &ranges,
            maximum,
            nonnegative,
            pick(true),
            pick(false),
        );
        return Ok(Read::Linear(call, xs));
    }
    // A nonlinear program with integer variables: not yet.
    if whole.contains(&true) {
        return Ok(Read::Not(None));
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
    let Some((lower, upper)) = bounds(ev.context(), &xs, ranges, nonnegative) else {
        return Ok(Read::Not(None));
    };
    let Some(start) = initial_point(ev, &xs, initial)? else {
        return Ok(Read::Not(None));
    };
    let objective = ev.eval(objective)?;
    let minimized = if maximum {
        ev.eval(&(Expr::from(0) - objective.clone()))?
    } else {
        objective.clone()
    };
    Ok(Read::Nonlinear(Program {
        xs,
        objective,
        minimized,
        equalities,
        inequalities,
        lower,
        upper,
        start,
        limit: iterations.unwrap_or(MAX_OUTER),
        precision,
    }))
}

/// Each variable's ends: at least 0 when `nonnegative`, each range
/// intersected; `None` for an end that is no real number.
fn bounds(
    ctx: &Context,
    xs: &[Expr],
    ranges: Vec<(&Expr, &Expr, &Expr)>,
    nonnegative: bool,
) -> Option<(Ends, Ends)> {
    let n = xs.len();
    let zero = End {
        exact: Expr::from(0),
        value: 0.0,
    };
    let mut lower: Vec<Option<End>> = vec![nonnegative.then(|| zero.clone()); n];
    let mut upper: Vec<Option<End>> = vec![None; n];
    for (x, a, b) in ranges {
        let j = xs.iter().position(|v| v == x).expect("collected");
        let (Ok(low), Ok(high)) = (end(ctx, a), end(ctx, b)) else {
            return None;
        };
        if let Some(low) = low
            && lower[j].as_ref().is_none_or(|l| l.value < low.value)
        {
            lower[j] = Some(low);
        }
        if let Some(high) = high
            && upper[j].as_ref().is_none_or(|h| h.value > high.value)
        {
            upper[j] = Some(high);
        }
    }
    Some((lower, upper))
}

/// `nlp_initialpoint`'s values, a region's (all ranges) midpoints; `None`
/// for anything else, giac's several points too.
fn initial_point(
    ev: &mut Eval<'_>,
    xs: &[Expr],
    initial: &[Expr],
) -> Result<Option<Vec<Option<Expr>>>, EvalError> {
    let regions = initial
        .iter()
        .filter(|p| p.args().get(1).and_then(head) == Some(Symbol::RANGE))
        .count();
    if regions != 0 && regions != initial.len() {
        return Ok(None);
    }
    let mut start = vec![None; xs.len()];
    for p in initial {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, v]) = (head(p), p.args()) else {
            return Ok(None);
        };
        let Some(j) = xs.iter().position(|y| y == x) else {
            return Ok(None);
        };
        start[j] = Some(match (head(v), v.args()) {
            (Some(Symbol::RANGE), [a, b]) => ev.eval(&((a.clone() + b.clone()) / Expr::from(2)))?,
            (Some(Symbol::LIST), _) => return Ok(None),
            _ => v.clone(),
        });
    }
    Ok(Some(start))
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
    /// Stopped by `nlp_iterationlimit` before it converged.
    Unconverged,
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

/// The augmented Lagrangian's outer loop from `x0`, at most `limit`
/// iterations (and [`MAX_OUTER`]).
fn lagrangian(
    c: &Compiled,
    x0: Vec<f64>,
    lo: &[f64],
    hi: &[f64],
    limit: usize,
) -> Result<Search, Failure> {
    let mut m = Multipliers {
        lambda: vec![0.0; c.equalities.len()],
        mu: vec![0.0; c.inequalities.len()],
        rho: 10.0,
    };
    let mut x = x0;
    let mut previous = f64::INFINITY;
    let mut tol = 1e-6;
    let mut converged = false;
    for _ in 0..limit.min(MAX_OUTER) {
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
            converged = true;
            break;
        }
        if violation > previous / 4.0 && m.rho < MAX_PENALTY {
            m.rho *= 10.0;
        }
        previous = violation;
        tol = (tol / 10.0).max(1e-11);
    }
    if !converged && limit < MAX_OUTER {
        return Err(Failure::Unconverged);
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

/// The search's point certified as a Karush–Kuhn–Tucker point at `bits`
/// bits: the answer, or `None`.
fn certify(
    ev: &mut Eval<'_>,
    p: &Program,
    c: &Compiled,
    s: &Search,
    lo: &[f64],
    hi: &[f64],
    bits: u64,
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
        match certified_root_at(ev, equations, &unknowns, guess, bits)? {
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
