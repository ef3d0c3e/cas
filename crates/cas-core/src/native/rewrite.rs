//! Rewriting forms (CAS-DOC-36): Xcas's commands that turn exponentials,
//! logarithms, powers, trigonometric, hyperbolic and inverse trigonometric
//! functions into one another.
//!
//! Each builtin takes one expression and rewrites every part of it, lists
//! and equations among them, in one pass from the leaves up: a part a rule
//! rewrote is not visited again, and what it does not match stays (§1).
//! The answer is evaluated, as every native's is.
//!
//! - [`exponential`]: `HyperbolicToExp`, `ExpToPower`, `PowerToExp`,
//!   `ExpExpand`, `LogExpand` and `LogCombine` (§2).

mod exponential;
#[cfg(test)]
mod testing;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// A rule: the rewrite of one part, whose own parts it has rewritten
/// already, or `None` to keep it.
type Rule = fn(&mut Eval<'_>, &Expr) -> Result<Option<Expr>, EvalError>;

/// `e = f[x]` with `rule` applied to each part of `x` from the leaves up,
/// the step noted as `name`.
fn rewrite(ev: &mut Eval<'_>, e: &Expr, name: &'static str, rule: Rule) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let mut failure = None;
    let result = x.map_bottom_up(|part| {
        if failure.is_some() {
            return part.clone();
        }
        match rule(ev, part) {
            Ok(Some(new)) => new,
            Ok(None) => part.clone(),
            Err(error) => {
                failure = Some(error);
                part.clone()
            }
        }
    });
    if let Some(error) = failure {
        return Err(error);
    }
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// The head of `e` if it is a symbol.
fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// The terms of a sum, or `u` alone.
fn terms(u: &Expr) -> &[Expr] {
    match head(u) {
        Some(Symbol::PLUS) => u.args(),
        _ => std::slice::from_ref(u),
    }
}

/// The factors of a product, or `u` alone.
fn factors(u: &Expr) -> &[Expr] {
    match head(u) {
        Some(Symbol::TIMES) => u.args(),
        _ => std::slice::from_ref(u),
    }
}

/// `HyperbolicToExp[x]`, Xcas's `hyp2exp`.
pub(super) fn hyperbolic_to_exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "HyperbolicToExp", exponential::hyperbolic)
}

/// `ExpToPower[x]`, Xcas's `exp2pow`.
pub(super) fn exp_to_power(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ExpToPower", exponential::exp_to_power)
}

/// `PowerToExp[x]`, Xcas's `pow2exp`.
pub(super) fn power_to_exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "PowerToExp", exponential::power_to_exp)
}

/// `ExpExpand[x]`, Xcas's `expexpand`.
pub(super) fn exp_expand(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ExpExpand", exponential::exp_expand)
}

/// `LogExpand[x]`, Xcas's `lnexpand`.
pub(super) fn log_expand(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "LogExpand", exponential::log_expand)
}

/// `LogCombine[x]`, Xcas's `lncollect`.
pub(super) fn log_combine(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "LogCombine", exponential::log_combine)
}
