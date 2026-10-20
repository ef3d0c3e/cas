//! [`Trace::verify`]: each step claims `before = after` under the context's
//! assumptions, checked at sample points (CAS-DOC-4 §9).
//!
//! A point gives each free symbol of a step an exact value: an integer
//! where the facts say integer, a rational that is not an integer in the
//! Real domain or where the facts say real, a Gaussian rational otherwise.
//! It must satisfy each assumption that mentions its symbols, evaluated at
//! it; a point where one is not `True` is drawn again, and the symbols such
//! an assumption ties to the step's get values too. Both sides are evaluated
//! at the point. Exact values (rationals, Gaussian rationals, `True` and
//! `False`, lists of them) are compared; other values agree when they are
//! equal and tell nothing otherwise, nor do infinities and failures. A step
//! whose sides differ at a point is wrong; one where no point told is
//! unchecked. So rational expressions are checked now, and those whose
//! values need balls (`Sin[1/2]`) when Balls lands.

use foldhash::fast::FixedState;
use hashbrown::HashSet;

use super::{Step, StepKind, Trace};
use crate::assume::Domain;
use crate::context::Context;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::{Complex, Rational};

/// The points a step must agree at.
const POINTS: usize = 4;

/// The points drawn for a step at most, counting those rejected.
const DRAWS: usize = 32;

/// What [`Trace::verify`] found.
#[derive(Clone, Debug, PartialEq)]
pub struct Verification<'t> {
    /// The first step, in the order taken, whose sides differ at a point;
    /// checking stops there.
    pub wrong: Option<Counterexample<'t>>,
    /// The steps no point told about: their sides are not exact values at
    /// the points drawn, or no point drawn satisfies the assumptions.
    pub unchecked: Vec<&'t Step>,
    /// How many steps agreed at every point that told.
    pub passed: usize,
}

/// A step whose sides differ at a point.
#[derive(Clone, Debug, PartialEq)]
pub struct Counterexample<'t> {
    /// The step.
    pub step: &'t Step,
    /// The value of each free symbol.
    pub point: Vec<(Symbol, Expr)>,
    /// `before` at the point, evaluated.
    pub before: Expr,
    /// `after` at the point, evaluated.
    pub after: Expr,
}

impl Trace {
    /// Checks each step at exact sample points that satisfy the
    /// assumptions of `ctx`, in the order the steps were taken, up to the
    /// first wrong one: "the result is wrong" becomes "this rule, at this
    /// position, is wrong" (CAS-DOC-4 §9). The points are drawn the same
    /// way each time.
    ///
    /// ```
    /// use cas_core::{Context, EvalOptions, Expr, Symbol, TraceLevel};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let mut ctx = Context::new();
    /// // A deliberately wrong rule: `x^2` is `2 x`.
    /// ctx.define_up_delayed(Symbol::new("x"), &read("Power[x, 2]"), &read("Times[2, x]"))
    ///     .unwrap();
    /// let options = EvalOptions { trace: TraceLevel::Rules };
    /// let (result, report) = ctx.eval_with(&read("Plus[Power[x, 2], 1]"), &options);
    /// assert_eq!(result.unwrap(), read("Plus[1, Times[2, x]]"));
    /// let wrong = report.trace.verify(&ctx).wrong.unwrap();
    /// assert_eq!(wrong.step.before, read("Power[x, 2]"));
    /// assert_eq!(wrong.step.position, [1]);
    /// ```
    #[must_use]
    pub fn verify(&self, ctx: &Context) -> Verification<'_> {
        let mut verification = Verification {
            wrong: None,
            unchecked: Vec::new(),
            passed: 0,
        };
        for step in &self.steps {
            if step.before == step.after || matches!(step.kind, StepKind::Attempts(_)) {
                continue;
            }
            match check(ctx, step) {
                Outcome::Agrees => verification.passed += 1,
                Outcome::Unchecked => verification.unchecked.push(step),
                Outcome::Differs(point, before, after) => {
                    verification.wrong = Some(Counterexample {
                        step,
                        point,
                        before,
                        after,
                    });
                    break;
                }
            }
        }
        verification
    }
}

enum Outcome {
    Agrees,
    Unchecked,
    Differs(Vec<(Symbol, Expr)>, Expr, Expr),
}

/// What a value at a point tells.
enum Value {
    /// A number, a truth value, or a list of them.
    Exact(Expr),
    Other(Expr),
    /// An infinity, `Indeterminate`, or a failure.
    Nothing,
}

/// What each kind of free symbol is drawn from.
#[derive(Clone, Copy)]
enum Range {
    Integers,
    Reals,
    Complexes,
}

fn check(ctx: &Context, step: &Step) -> Outcome {
    let Some(symbols) = free_symbols(ctx, &[&step.before, &step.after]) else {
        return Outcome::Unchecked;
    };
    let (symbols, assumptions) = with_assumptions(ctx, symbols);
    let ranges: Vec<Range> = symbols.iter().map(|&s| range(ctx, s)).collect();
    let mut rng = SplitMix(0x9e37_79b9_7f4a_7c15 ^ step.id as u64);
    let mut told = 0;
    let draws = if symbols.is_empty() { 1 } else { DRAWS };
    for _ in 0..draws {
        let point: Vec<(Symbol, Expr)> = symbols
            .iter()
            .zip(&ranges)
            .map(|(&s, &r)| (s, rng.value(r)))
            .collect();
        let rules: Vec<(Expr, Expr)> = point
            .iter()
            .map(|(s, v)| (Expr::from(*s), v.clone()))
            .collect();
        let holds = assumptions.iter().all(|a| {
            ctx.eval(&a.substitute(&rules))
                .is_ok_and(|t| t.as_symbol() == Some(Symbol::TRUE))
        });
        if !holds {
            continue;
        }
        let at = |side: &Expr| value(ctx.eval(&side.substitute(&rules)));
        match (at(&step.before), at(&step.after)) {
            (Value::Exact(b), Value::Exact(a)) if a != b => {
                return Outcome::Differs(point, b, a);
            }
            (Value::Exact(b) | Value::Other(b), Value::Exact(a) | Value::Other(a)) if a == b => {
                told += 1;
            }
            _ => {}
        }
        if told == POINTS || symbols.is_empty() {
            break;
        }
    }
    if told > 0 {
        Outcome::Agrees
    } else {
        Outcome::Unchecked
    }
}

/// The symbols of `sides` that stand for values: not builtins, without a
/// value in `ctx`, and never a head. `None` when a side binds or matches
/// (`Function`, a pattern), where values cannot be put in.
fn free_symbols(ctx: &Context, sides: &[&Expr]) -> Option<Vec<Symbol>> {
    let mut heads: HashSet<Symbol, FixedState> = HashSet::default();
    let mut symbols: Vec<Symbol> = Vec::new();
    for side in sides {
        if side.has_pattern() {
            return None;
        }
        for part in side.distinct() {
            match part.view() {
                View::Apply { head, .. } => {
                    if let Some(h) = head.as_symbol() {
                        if h == Symbol::FUNCTION {
                            return None;
                        }
                        heads.insert(h);
                    }
                }
                View::Symbol(s)
                    if !s.is_builtin() && ctx.value(s).is_none() && !symbols.contains(&s) =>
                {
                    symbols.push(s);
                }
                _ => {}
            }
        }
    }
    symbols.retain(|s| !heads.contains(s));
    symbols.sort_by_cached_key(|&s| Expr::from(s));
    Some(symbols)
}

/// `symbols` and the assumptions of `ctx` that mention them, with the free
/// symbols those mention, until none is added.
fn with_assumptions(ctx: &Context, mut symbols: Vec<Symbol>) -> (Vec<Symbol>, Vec<Expr>) {
    let mut relevant = vec![false; ctx.assumptions().len()];
    loop {
        let mut added = false;
        for (a, taken) in ctx.assumptions().iter().zip(&mut relevant) {
            if *taken
                || !symbols
                    .iter()
                    .any(|&s| a.may_contain(s) && a.contains(&s.into()))
            {
                continue;
            }
            *taken = true;
            for s in free_symbols(ctx, &[a]).unwrap_or_default() {
                if !symbols.contains(&s) {
                    symbols.push(s);
                    added = true;
                }
            }
        }
        if !added {
            break;
        }
    }
    let assumptions = ctx
        .assumptions()
        .iter()
        .zip(&relevant)
        .filter(|(_, taken)| **taken)
        .map(|(a, _)| a.clone())
        .collect();
    (symbols, assumptions)
}

/// What the facts and the domain say `s` ranges over.
fn range(ctx: &Context, s: Symbol) -> Range {
    let element = |set: Symbol| Expr::apply(&Symbol::ELEMENT.into(), &[s.into(), set.into()]);
    if ctx.decide(&element(Symbol::INTEGERS)) == crate::Truth::True {
        Range::Integers
    } else if ctx.domain() == Domain::Real
        || ctx.decide(&element(Symbol::REALS)) == crate::Truth::True
    {
        Range::Reals
    } else {
        Range::Complexes
    }
}

fn value(result: Result<Expr, crate::EvalError>) -> Value {
    let Ok(v) = result else {
        return Value::Nothing;
    };
    if [
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
        Symbol::INFINITY,
        Symbol::INDETERMINATE,
    ]
    .into_iter()
    .any(|s| v.may_contain(s) && v.contains(&s.into()))
    {
        return Value::Nothing;
    }
    if is_exact(&v) {
        Value::Exact(v)
    } else {
        Value::Other(v)
    }
}

/// Whether `v` is a number, `True`, `False`, or a list of such.
fn is_exact(v: &Expr) -> bool {
    match v.view() {
        View::Symbol(s) => s == Symbol::TRUE || s == Symbol::FALSE,
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::LIST) => {
            args.iter().all(is_exact)
        }
        _ => v.is_number(),
    }
}

/// `SplitMix64`: the same points each time, for a step.
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// An integer in `[lo, hi]`.
    fn between(&mut self, lo: i64, hi: i64) -> i64 {
        let width = hi.abs_diff(lo) + 1;
        lo + i64::try_from(self.next() % width).expect("a small width")
    }

    /// A rational `p/q` with `|p| <= 40` and `2 <= q <= 12` that is not an
    /// integer: integers, among them 0 and ±1, where identities that do not
    /// hold elsewhere might, come only where the facts say integer.
    fn rational(&mut self) -> Rational {
        loop {
            let (p, q) = (self.between(-40, 40), self.between(2, 12));
            if p % q != 0 {
                return Rational::new(p.into(), q.into());
            }
        }
    }

    fn value(&mut self, range: Range) -> Expr {
        match range {
            Range::Integers => loop {
                let k = self.between(-30, 30);
                if k.abs() > 1 {
                    return Expr::from(k);
                }
            },
            Range::Reals => Expr::from(self.rational()),
            Range::Complexes => Expr::from(Complex::new(self.rational(), self.rational())),
        }
    }
}
