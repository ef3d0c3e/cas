//! Xcas's session commands (CAS-619, CAS-620), as natives.
//!
//! - **The history**: `Ans[n]` (`ans`, giac prog.cc:8611) answers the
//!   answer of an earlier statement of the session, `Quest[n]` (`quest`,
//!   `entry`, prog.cc:8643) its input, as a string where giac keeps the
//!   expression. `n ≥ 0` is the `(n+1)`th entry the history keeps, `n < 0`
//!   the `(−n)`th from the end, `−1` by default. An empty history answers
//!   `undef`, as giac does; an index out of range, not an integer, or more
//!   than one argument leaves the call with a warning. The session records
//!   each statement ([`Context::record`](crate::Context::record)); reading
//!   the history makes the evaluation volatile ([`Eval::volatile`]), so
//!   that nothing it gave is stamped or memoised.
//! - `Version[]` (`version`, usual.cc:8153): this CAS's name and version.
//! - `NoDisplay[e]` (`nodisp`, prog.cc:10341): `"Done"`, as giac's.
//! - `Restart[]` (`restart`, maple.cc:492), `RemoveAllVariables[]`
//!   (`rm_all_vars`, usual.cc:251), `RemoveAToZ[]` (`rm_a_z`, usual.cc:231),
//!   `Variables[]` (`VARS`, rpn.cc:754), `Time[e]` (`time`, maple.cc:545)
//!   and `Monotonic[]` (`monotonic`, prog.cc:11072) change or read the
//!   session, which runs them as statements; nested in an expression, they
//!   stay with a warning, where giac answers.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// An entry of the history, by the index of `e`, through `pick`.
fn entry(ev: &mut Eval<'_>, e: &Expr, pick: impl FnOnce(&(String, Expr)) -> Expr) -> Native {
    let n = match e.args() {
        [] => -1,
        [n] => match n.as_i64() {
            Some(n) => n,
            None => return refuse(ev, e, "the index of the history is an integer"),
        },
        _ => return refuse(ev, e, "the history takes one index at most"),
    };
    ev.volatile();
    let found = {
        let history = ev.context().history();
        if history.is_empty() {
            return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
        }
        history.get(n).map(pick)
    };
    match found {
        Some(answer) => Ok(Some(answer)),
        None => refuse(ev, e, "the history has no entry at this index"),
    }
}

/// `Ans[]`, `Ans[n]`.
pub(super) fn ans(ev: &mut Eval<'_>, e: &Expr) -> Native {
    entry(ev, e, |(_, answer)| answer.clone())
}

/// `Quest[]`, `Quest[n]`.
pub(super) fn quest(ev: &mut Eval<'_>, e: &Expr) -> Native {
    entry(ev, e, |(input, _)| Expr::string(input))
}

/// `Version[]`.
pub(super) fn version(_: &mut Eval<'_>, e: &Expr) -> Native {
    if !e.args().is_empty() {
        return Ok(None);
    }
    Ok(Some(Expr::string(concat!(
        "cas ",
        env!("CARGO_PKG_VERSION")
    ))))
}

/// `NoDisplay[e]`.
pub(super) fn no_display(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [_] => Ok(Some(Expr::string("Done"))),
        _ => Ok(None),
    }
}

/// A command the session runs as a statement, nested: with no arguments,
/// it stays with a warning; otherwise it stays.
fn statement_only(ev: &mut Eval<'_>, e: &Expr, arity: usize, why: &'static str) -> Native {
    if e.args().len() != arity {
        return Ok(None);
    }
    refuse(ev, e, why)
}

/// `Restart[]`, nested.
pub(super) fn restart(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 0, "restart runs as a statement of the session")
}

/// `RemoveAllVariables[]`, nested.
pub(super) fn remove_all_variables(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 0, "rm_all_vars runs as a statement of the session")
}

/// `RemoveAToZ[]`, nested.
pub(super) fn remove_a_to_z(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 0, "rm_a_z runs as a statement of the session")
}

/// `Variables[]`, nested.
pub(super) fn variables(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 0, "VARS runs as a statement of the session")
}

/// `Time[e]`, nested.
pub(super) fn time(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 1, "time runs as a statement of the session")
}

/// `Monotonic[]`, nested.
pub(super) fn monotonic(ev: &mut Eval<'_>, e: &Expr) -> Native {
    statement_only(ev, e, 0, "monotonic runs as a statement of the session")
}
