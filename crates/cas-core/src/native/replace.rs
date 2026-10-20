//! `Replace`, `ReplaceAll`, `ReplaceRepeated` and `MatchQ` (CAS-DOC-4 §5.3).
//!
//! Rules are `Rule[lhs, rhs]` (the right-hand side evaluated with the rule),
//! `RuleDelayed[lhs, rhs]` (held), or a list of them, compiled on each call
//! with the attributes of the evaluation's context.
//! A malformed rule or pattern leaves the call unevaluated. What a rule
//! rewrites is evaluated after, by the loop.

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
