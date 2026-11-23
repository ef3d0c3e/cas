//! `Replace`, `ReplaceAll`, `ReplaceRepeated` and `MatchQ` (CAS-DOC-4 §5.3).
//!
//! Rules are `Rule[lhs, rhs]` (the right-hand side evaluated with the rule),
//! `RuleDelayed[lhs, rhs]` (held), or a list of them, compiled on each call
//! with the attributes of the evaluation's context.
//! A malformed rule or pattern leaves the call unevaluated. What a rule
//! rewrites is evaluated after, by the loop.
//!
//! `AlgebraicSubstitute[l = r, e]`, Xcas's `algsubs` (CAS-566), giac's
//! `_algsubs` (solve.cc): `l = r` (`Equation` or `Equal`) used as an
//! algebraic relation, not a pattern. With `n/d` the `Together` form of `e`
//! and `ln/ld` that of `l`, `Eliminate[{ln − ld r, n − d T}, ids]` for a
//! fresh `T`, `ids` the symbols of the equation in order of first
//! appearance (left side, then right, depth first), dropping the last id
//! until a relation is left; `Solve` of it for `T`, its first root under
//! `Together`: `algsubs(x^2=u,1+x^2+x^4)` is `u^2+u+1`. No relation for any
//! prefix, or no root: `e`. A `T` left in the root, an `Eliminate` or
//! `Solve` that stays or fails: the call stays.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::pattern::{self, Pattern};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `ReplaceAll[e, rules]`: every part of `e` that a rule rewrites, from the
/// top; a rewritten part is not looked into again.
pub(super) fn replace_all(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [subject, rules] = e.args() else {
        return Ok(None);
    };
    let Some(rules) = pattern::rules(rules, ev.context()) else {
        return Ok(None);
    };
    pattern::replace_all(ev, &rules, subject).map(Some)
}

/// `Replace[e, rules]`: `e` itself, by the first rule that rewrites it.
pub(super) fn replace(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [subject, rules] = e.args() else {
        return Ok(None);
    };
    let Some(rules) = pattern::rules(rules, ev.context()) else {
        return Ok(None);
    };
    Ok(Some(
        pattern::rewrite(ev, &rules, subject)?.unwrap_or_else(|| subject.clone()),
    ))
}

/// `ReplaceRepeated[e, rules]`: `ReplaceAll`, evaluated, again until nothing
/// changes; each round counts against `Limits::iterations`.
pub(super) fn replace_repeated(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [subject, rules] = e.args() else {
        return Ok(None);
    };
    let Some(rules) = pattern::rules(rules, ev.context()) else {
        return Ok(None);
    };
    let mut current = subject.clone();
    for _ in 0..ev.context().limits().iterations {
        let replaced = pattern::replace_all(ev, &rules, &current)?;
        let next = ev.eval(&replaced)?;
        if next == current {
            return Ok(Some(next));
        }
        current = next;
    }
    Err(EvalError::IterationLimit)
}

/// `MatchQ[e, pattern]`.
pub(super) fn match_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [subject, p] = e.args() else {
        return Ok(None);
    };
    let Ok(p) = Pattern::with_context(p, ev.context()) else {
        return Ok(None);
    };
    let matched = p.matches(ev, subject)?;
    Ok(Some(Expr::from(if matched {
        Symbol::TRUE
    } else {
        Symbol::FALSE
    })))
}

/// The symbols of `e` that are no numeric constant, heads aside, in order
/// of first appearance, depth first.
fn symbols(e: &Expr, out: &mut Vec<Expr>) {
    if e.as_symbol().is_some() && !e.is_numeric() {
        if !out.contains(e) {
            out.push(e.clone());
        }
    } else {
        for a in e.args() {
            symbols(a, out);
        }
    }
}

/// `Numerator` and `Denominator` of `Together[e]`; `None` when either stays.
fn fraction(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<[Expr; 2]>, EvalError> {
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [e.clone()]))?;
    let n = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let d = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    let stays = |p: &Expr, s| p.head().and_then(Expr::as_symbol) == Some(s);
    if stays(&n, Symbol::NUMERATOR) || stays(&d, Symbol::DENOMINATOR) {
        return Ok(None);
    }
    Ok(Some([n, d]))
}

/// The elements of a list, `None` for anything else.
fn elements(e: &Expr) -> Option<&[Expr]> {
    (e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)).then(|| e.args())
}

/// `AlgebraicSubstitute[l = r, e]`.
pub(super) fn algebraic_substitute(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [eq, subject] = e.args() else {
        return Ok(None);
    };
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) =
        (eq.head().and_then(Expr::as_symbol), eq.args())
    else {
        return Ok(None);
    };
    // Elimination and solving are heavy and can fail on what the fuzzer
    // builds: an error other than an interrupt leaves the call.
    match substitute(ev, eq, l, r, subject) {
        Err(EvalError::Interrupted) => Err(EvalError::Interrupted),
        Err(_) => Ok(None),
        ok => ok,
    }
}

/// giac's loop: eliminate, dropping the last id until a relation is left,
/// then solve it for the fresh `T`.
fn substitute(ev: &mut Eval<'_>, eq: &Expr, lhs: &Expr, rhs: &Expr, subject: &Expr) -> Native {
    let t = Expr::from(Symbol::new("AlgebraicSubstitute$T"));
    let Some([num, den]) = fraction(ev, subject)? else {
        return Ok(None);
    };
    let Some([ln, ld]) = fraction(ev, lhs)? else {
        return Ok(None);
    };
    let ee = num - den * t.clone();
    let eq1 = ln - ld * rhs.clone();
    let mut ids = Vec::new();
    symbols(eq, &mut ids);
    let mut relation = None;
    while !ids.is_empty() {
        let polys = Expr::call(Symbol::LIST, [eq1.clone(), ee.clone()]);
        let vars = Expr::call(Symbol::LIST, ids.clone());
        let sol = ev.eval(&Expr::call(Symbol::ELIMINATE, [polys, vars]))?;
        let Some(found) = elements(&sol) else {
            return Ok(None);
        };
        if !found.is_empty() {
            relation = Some(sol);
            break;
        }
        ids.pop();
    }
    let Some(relation) = relation else {
        return Ok(Some(subject.clone()));
    };
    let roots = ev.eval(&Expr::call(Symbol::SOLVE, [relation, t.clone()]))?;
    let Some(roots) = elements(&roots) else {
        return Ok(None);
    };
    let Some(root) = roots.first() else {
        return Ok(Some(subject.clone()));
    };
    if root.contains(&t) {
        return Ok(None);
    }
    ev.eval(&Expr::call(Symbol::TOGETHER, [root.clone()]))
        .map(Some)
}
