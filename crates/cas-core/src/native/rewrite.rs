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
//! - [`substitution`]: `HalfTan`, `HalfTanHyperbolicToExp`, `TanToSinCos`,
//!   `SinToCosTan`, `CosToSinTan`, `TanToSinCos2` and `TanToCosSin2` (§3).
//! - [`pythagorean`]: `TrigToCos`, `TrigToSin` and `TrigToTan` (§4), which
//!   take passes of their own, with evaluation between.
//! - [`collect`]: `TrigCollect` (§5), after `TrigReduce`.
//! - [`inverse`]: `ArcCosToArcSin`, `ArcSinToArcCos`, `ArcCosToArcTan`,
//!   `ArcSinToArcTan`, `ArcTanToArcCos`, `ArcTanToArcSin` and
//!   `ArcTrigToLog` (§6).
//! - [`complex`]: `ComplexExpand`, `MultiplyConjugate` and
//!   `MultiplyComplexConjugate` (§7), on the whole expression.

mod collect;
mod complex;
mod exponential;
mod inverse;
mod pythagorean;
mod substitution;
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

/// `x` with `form` applied to it, or to each entry of a list or relation,
/// for the forms that take one expression whole.
fn entrywise(
    ev: &mut Eval<'_>,
    x: &Expr,
    form: &mut dyn FnMut(&mut Eval<'_>, &Expr) -> Result<Expr, EvalError>,
) -> Result<Expr, EvalError> {
    if let (
        Some(
            Symbol::LIST
            | Symbol::EQUAL
            | Symbol::UNEQUAL
            | Symbol::LESS
            | Symbol::LESS_EQUAL
            | Symbol::GREATER
            | Symbol::GREATER_EQUAL,
        ),
        Some(h),
    ) = (head(x), x.head())
    {
        let entries = x
            .args()
            .iter()
            .map(|a| entrywise(ev, a, form))
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(Expr::apply(h, &entries));
    }
    form(ev, x)
}

/// `e = f[x]` with `form` applied entrywise, the step noted as `name`.
fn whole(
    ev: &mut Eval<'_>,
    e: &Expr,
    name: &'static str,
    form: &mut dyn FnMut(&mut Eval<'_>, &Expr) -> Result<Expr, EvalError>,
) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = entrywise(ev, x, form)?;
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

/// `HalfTan[x]`, Xcas's `halftan`.
pub(super) fn half_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "HalfTan", substitution::half_tan)
}

/// `HalfTanHyperbolicToExp[x]`, Xcas's `halftan_hyp2exp`.
pub(super) fn half_tan_hyperbolic_to_exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(
        ev,
        e,
        "HalfTanHyperbolicToExp",
        substitution::half_tan_hyperbolic,
    )
}

/// `TanToSinCos[x]`, Xcas's `tan2sincos`.
pub(super) fn tan_to_sin_cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "TanToSinCos", substitution::tan_to_sin_cos)
}

/// `SinToCosTan[x]`, Xcas's `sin2costan`.
pub(super) fn sin_to_cos_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "SinToCosTan", substitution::sin_to_cos_tan)
}

/// `CosToSinTan[x]`, Xcas's `cos2sintan`.
pub(super) fn cos_to_sin_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "CosToSinTan", substitution::cos_to_sin_tan)
}

/// `TanToSinCos2[x]`, Xcas's `tan2sincos2`.
pub(super) fn tan_to_sin_cos_2(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "TanToSinCos2", substitution::tan_to_sin_cos_2)
}

/// `TanToCosSin2[x]`, Xcas's `tan2cossin2`.
pub(super) fn tan_to_cos_sin_2(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "TanToCosSin2", substitution::tan_to_cos_sin_2)
}

/// `e = f[x]` in the Pythagorean form favouring `favoured`, the step noted
/// as `name`.
fn pythagorean_form(
    ev: &mut Eval<'_>,
    e: &Expr,
    name: &'static str,
    favoured: pythagorean::Favoured,
) -> Native {
    whole(ev, e, name, &mut |ev, x| pythagorean::form(ev, x, favoured))
}

/// `TrigToCos[x]`, Xcas's `trigcos`.
pub(super) fn trig_to_cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    pythagorean_form(ev, e, "TrigToCos", pythagorean::Favoured::Cos)
}

/// `TrigToSin[x]`, Xcas's `trigsin`.
pub(super) fn trig_to_sin(ev: &mut Eval<'_>, e: &Expr) -> Native {
    pythagorean_form(ev, e, "TrigToSin", pythagorean::Favoured::Sin)
}

/// `TrigToTan[x]`, Xcas's `trigtan`.
pub(super) fn trig_to_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    pythagorean_form(ev, e, "TrigToTan", pythagorean::Favoured::Tan)
}

/// `TrigCollect[x]`, Xcas's `tcollect`.
pub(super) fn trig_collect(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = collect::collect(ev, x)?;
    ev.step("TrigCollect", e, &result);
    Ok(Some(result))
}

/// `ArcCosToArcSin[x]`, Xcas's `acos2asin`.
pub(super) fn arc_cos_to_arc_sin(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcCosToArcSin", inverse::arc_cos_to_arc_sin)
}

/// `ArcSinToArcCos[x]`, Xcas's `asin2acos`.
pub(super) fn arc_sin_to_arc_cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcSinToArcCos", inverse::arc_sin_to_arc_cos)
}

/// `ArcCosToArcTan[x]`, Xcas's `acos2atan`.
pub(super) fn arc_cos_to_arc_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcCosToArcTan", inverse::arc_cos_to_arc_tan)
}

/// `ArcSinToArcTan[x]`, Xcas's `asin2atan`.
pub(super) fn arc_sin_to_arc_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcSinToArcTan", inverse::arc_sin_to_arc_tan)
}

/// `ArcTanToArcCos[x]`, Xcas's `atan2acos`.
pub(super) fn arc_tan_to_arc_cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcTanToArcCos", inverse::arc_tan_to_arc_cos)
}

/// `ArcTanToArcSin[x]`, Xcas's `atan2asin`.
pub(super) fn arc_tan_to_arc_sin(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcTanToArcSin", inverse::arc_tan_to_arc_sin)
}

/// `ArcTrigToLog[x]`, Xcas's `atrig2ln`.
pub(super) fn arc_trig_to_log(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rewrite(ev, e, "ArcTrigToLog", inverse::arc_trig_to_log)
}

/// `ComplexExpand[x]`, Xcas's `evalc`.
pub(super) fn complex_expand(ev: &mut Eval<'_>, e: &Expr) -> Native {
    whole(ev, e, "ComplexExpand", &mut complex::complex_expand)
}

/// `MultiplyConjugate[x]`, Xcas's `mult_conjugate`.
pub(super) fn multiply_conjugate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    whole(ev, e, "MultiplyConjugate", &mut |_, x| {
        Ok(complex::multiply_conjugate(x))
    })
}

/// `MultiplyComplexConjugate[x]`, Xcas's `mult_c_conjugate`.
pub(super) fn multiply_complex_conjugate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    whole(
        ev,
        e,
        "MultiplyComplexConjugate",
        &mut complex::multiply_complex_conjugate,
    )
}
