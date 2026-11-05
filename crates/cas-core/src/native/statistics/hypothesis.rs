//! Hypothesis tests (CAS-DOC-52): Xcas's `normalt`, `studentt`,
//! `chisquaret`, `kolmogorovt` and `wilcoxont`, with the laws and
//! statistics they take.
//!
//! - The alternative is one of the strings `"!="` (the default), `"<"` and
//!   `">"`, where giac takes its quoted operators; the level `α`, 0.05 by
//!   default, a real number in `(0, 1)`.
//! - A test answers giac's number: 1 when the null hypothesis is not
//!   rejected at `α`, 0 when it is, by the sign of `p − α`, decided exactly
//!   or by balls; undecided, the call stays with a warning. The statistic
//!   and the p-value go into the evaluation step's note, where giac prints
//!   its report.
//! - Location tests on a mean ([`location`], §2, CAS-503).
//! - Pearson's χ² ([`chisquare`], §3, CAS-504).
//! - Kolmogorov's distribution and the Kolmogorov–Smirnov test
//!   ([`kolmogorov`], §4, CAS-505, CAS-506).

pub(in crate::native) mod chisquare;
pub(in crate::native) mod kolmogorov;
pub(in crate::native) mod location;

use std::cmp::Ordering;
use std::fmt::Write as _;

use cas_num::Rational;

use super::refuse;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::rounding::real_sign;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The alternative hypothesis.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Alternative {
    /// `"!="`: two-sided.
    Unequal,
    /// `"<"`.
    Less,
    /// `">"`.
    Greater,
}

impl Alternative {
    /// The alternative written as the string `x`; `None` for another.
    fn of(x: &Expr) -> Option<Alternative> {
        match x.view() {
            View::String("!=" | "<>") => Some(Alternative::Unequal),
            View::String("<") => Some(Alternative::Less),
            View::String(">") => Some(Alternative::Greater),
            _ => None,
        }
    }

    fn text(self) -> &'static str {
        match self {
            Alternative::Unequal => "!=",
            Alternative::Less => "<",
            Alternative::Greater => ">",
        }
    }
}

/// The trailing `[alt], [α]` of a test; `None` after a warning.
pub(super) fn options(
    ev: &mut Eval<'_>,
    e: &Expr,
    rest: &[Expr],
) -> Result<Option<(Alternative, Expr)>, EvalError> {
    let (alternative, level) = match rest {
        [] => (None, None),
        [a] if is_string(a) => (Some(a), None),
        [a, alpha] => (Some(a), Some(alpha)),
        _ => {
            refuse(
                ev,
                e,
                "expected an alternative \"!=\", \"<\" or \">\" and a level",
            )?;
            return Ok(None);
        }
    };
    let alternative = match alternative.map(Alternative::of) {
        None => Alternative::Unequal,
        Some(Some(a)) => a,
        Some(None) => {
            refuse(ev, e, "the alternative is one of \"!=\", \"<\" and \">\"")?;
            return Ok(None);
        }
    };
    let alpha = match level {
        None => Expr::from(Rational::new(1.into(), 20.into())),
        Some(alpha) => {
            let inside = positive(ev, alpha)
                && real_sign(ev.context(), &ev.eval(&(Expr::from(1) - alpha.clone()))?)
                    == Some(Ordering::Greater);
            if !inside {
                refuse(ev, e, "the level is a number between 0 and 1")?;
                return Ok(None);
            }
            alpha.clone()
        }
    };
    Ok(Some((alternative, alpha)))
}

/// Whether a ball occurs in `x`.
pub(super) fn has_ball(x: &Expr) -> bool {
    matches!(x.view(), View::Ball(_) | View::ComplexBall(_)) || x.args().iter().any(has_ball)
}

/// `x` as one ball at the working precision when `inexact`, as `N` makes
/// it: an expression of balls left by exact steps.
pub(super) fn settled(ev: &mut Eval<'_>, x: Expr, inexact: bool) -> Result<Expr, EvalError> {
    if !inexact {
        return Ok(x);
    }
    let bits = ev.context().working_precision();
    crate::native::approximation(ev.context(), &x, bits)
}

/// Whether `x` is a string: an alternative, or a wrong one.
pub(super) fn is_string(x: &Expr) -> bool {
    matches!(x.view(), View::String(_))
}

/// Whether `x` is a real number, shown positive.
pub(super) fn positive(ev: &Eval<'_>, x: &Expr) -> bool {
    real_sign(ev.context(), x) == Some(Ordering::Greater)
}

/// Whether `x` is a real number, its sign decided.
pub(super) fn real(ev: &Eval<'_>, x: &Expr) -> bool {
    real_sign(ev.context(), x).is_some()
}

/// The p-value of a statistic `s` with the cumulative function `cdf`, for
/// a law symmetric about 0: `F(s)`, `1 − F(s)`, or `2 F(−|s|)`.
pub(super) fn symmetric_p(alternative: Alternative, s: &Expr, cdf: &dyn Fn(Expr) -> Expr) -> Expr {
    match alternative {
        Alternative::Less => cdf(s.clone()),
        Alternative::Greater => Expr::from(1) - cdf(s.clone()),
        Alternative::Unequal => {
            Expr::from(2) * cdf(Expr::from(-1) * Expr::call(Symbol::ABS, [s.clone()]))
        }
    }
}

/// The answer of a test of p-value `p` at level `alpha`: 0 when `p < α`,
/// 1 otherwise; the step notes `what`, the statistics `shown` and `p`.
pub(super) fn decide(
    ev: &mut Eval<'_>,
    e: &Expr,
    what: &str,
    shown: &[(&str, Expr)],
    p: &Expr,
    alpha: &Expr,
) -> Native {
    let p = ev.eval(p)?;
    let gap = ev.eval(&(p.clone() - alpha.clone()))?;
    let answer = match real_sign(ev.context(), &gap) {
        Some(Ordering::Less) => 0,
        Some(_) => 1,
        None => {
            return refuse(ev, e, "the p-value is not told apart from the level");
        }
    };
    let answer = Expr::from(answer);
    if ev.tracing() {
        let mut note = format!("{what}:");
        for (name, value) in shown.iter().chain(std::iter::once(&("p", p))) {
            let value = ev.eval(&Expr::call(Symbol::N, [value.clone(), Expr::from(10)]))?;
            let _ = write!(note, " {name} = {value},");
        }
        let _ = write!(note, " α = {alpha}");
        ev.step(note, e, &answer);
    }
    Ok(Some(answer))
}

/// The alternative's text, for notes.
pub(super) fn written(alternative: Alternative) -> &'static str {
    alternative.text()
}

#[cfg(test)]
mod tests;
