//! giac's expression, calculus, integer and box leftovers of the gap ledger
//! (CAS-640 to CAS-643), under giac's names. Each call is rewritten into
//! builtins for the loop to evaluate; where giac errs, the call stays with
//! a warning:
//!
//! - `NthRoot[a, b]` (`root`, CAS-640): `b^(1/a)`; `Indeterminate` for
//!   `a = 0`, giac's `undef`;
//! - `ArcTan2[y, x]` (`atan2`, CAS-640): the angle of the point `(x, y)`
//!   in `(−π, π]`. With the signs of `x` and `y` known: `0` at the origin
//!   and on the positive axis, `π` on the negative one, `±π/2` on the
//!   other axis, `atan(y/x) + (1 − sign x) sign y π/2` off the axes.
//!   Otherwise that expression, giac's symbolic form, or `sign y π/2` for
//!   `x = 0`;
//! - `Has[e, v]` (`has`, CAS-640): 1 when the user symbol `v` occurs in
//!   `e`, as a head too, else 0. giac looks for names only, so a builtin
//!   symbol or anything else gives 0;
//! - `FloorMod[a, b]` (`fmod`, CAS-640): `a − b ⌊a/b⌋`; `Indeterminate`
//!   for `b = 0`;
//! - `FractionalIntegerParts[a]` (`modf`, CAS-640): of a real number, the
//!   sequence of `a − trunc a` and `trunc a`, as decimals as giac's;
//! - `MultinomialProbability[n, [p…], [k…]]` (`multinomial`, CAS-642):
//!   `n!/∏ kᵢ! ∏ pᵢ^kᵢ` for nonnegative integers `kᵢ` summing to `n`;
//! - `ModularOrder[a % n]` (`order`, CAS-642): the order of `a` in
//!   `(ℤ/nℤ)*`, `MultiplicativeOrder[a, n]`; for `a` not invertible the
//!   call stays;
//! - `SignedFactorList[n]` (`maple_ifactors`, CAS-642): Maple's
//!   `[sign n, [[p, e], …]]`, `[0, []]` for 0;
//! - `AverageRateOfChange[f, x, h]` (`avgRC`, CAS-641): `(f(x + h) −
//!   f(x))/h`, `h` the decimal 0.001 by default as giac's; `Indeterminate`
//!   for `h = 0`;
//! - `BoxConstraints[[x…], [[a, b]…]]` (`box_constraints`, CAS-643): the
//!   sequence `x = a..b, …`.
//!
//! giac's `addtable` (`AddTable`, CAS-641) is a statement the session runs
//! (`Context::add_transform`); `additionally` (CAS-640) is `assume`, which
//! adds to the assumptions already.

use std::cmp::Ordering;

use cas_num::Integer;

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The call stays, with giac's error as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// Whether `x` is 0, exact or a decimal.
fn is_zero(ev: &mut Eval<'_>, x: &Expr) -> Result<bool, EvalError> {
    let equal = ev.eval(&Expr::call(Symbol::EQUAL, [x.clone(), Expr::from(0)]))?;
    Ok(equal.as_symbol() == Some(Symbol::TRUE))
}

/// Whether `x` is a list.
fn is_list(x: &Expr) -> bool {
    x.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// `Indeterminate`, giac's `undef`.
fn undefined() -> Native {
    Ok(Some(Expr::from(Symbol::INDETERMINATE)))
}

/// `NthRoot[a, b]`, giac's `root`.
pub(super) fn nth_root(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if is_zero(ev, a)? {
        return undefined();
    }
    Ok(Some(b.clone().pow(a.clone().pow(-1))))
}

/// `ArcTan2[y, x]`, giac's `atan2`.
pub(super) fn arc_tan2(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [y, x] = e.args() else {
        return Ok(None);
    };
    let half_pi = || Expr::from(Symbol::PI) / Expr::from(2);
    let ordinal = |o: Ordering| match o {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    };
    let sy = real_sign(ev.context(), y).map(ordinal);
    let sx = real_sign(ev.context(), x).map(ordinal);
    let result = match (sy, sx) {
        (Some(0), Some(0 | 1)) => Expr::from(0),
        (Some(0), Some(_)) => Expr::from(Symbol::PI),
        (Some(sy), Some(0)) => Expr::from(sy) * half_pi(),
        (Some(sy), Some(sx)) => {
            let atan = Expr::call(Symbol::ARC_TAN, [y.clone() / x.clone()]);
            atan + Expr::from((1 - sx) * sy) * half_pi()
        }
        (_, Some(0)) => Expr::call(Symbol::SIGN, [y.clone()]) * half_pi(),
        _ => {
            let atan = Expr::call(Symbol::ARC_TAN, [y.clone() / x.clone()]);
            let sign = |v: &Expr| Expr::call(Symbol::SIGN, [v.clone()]);
            atan + (Expr::from(1) - sign(x)) * sign(y) * half_pi()
        }
    };
    Ok(Some(result))
}

/// `Has[e, v]`, giac's `has`.
pub(super) fn has(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, v] = e.args() else {
        return Ok(None);
    };
    let name = v.as_symbol().is_some_and(|s| !s.is_builtin());
    Ok(Some(Expr::from(i64::from(name && x.contains(v)))))
}

/// `FloorMod[a, b]`, giac's `fmod`.
pub(super) fn floor_mod(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if is_zero(ev, b)? {
        return undefined();
    }
    let floor = Expr::call(Symbol::FLOOR, [a.clone() / b.clone()]);
    Ok(Some(a.clone() - b.clone() * floor))
}

/// `FractionalIntegerParts[a]`, giac's `modf`.
pub(super) fn fractional_integer_parts(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    if a.to_rational().is_none() && !matches!(a.view(), View::Ball(_)) {
        return refuse(ev, e, "a real number is needed");
    }
    let whole = Expr::call(Symbol::INTEGER_PART, [a.clone()]);
    let decimal = |x: Expr| Expr::call(Symbol::N, [x]);
    Ok(Some(Expr::call(
        Symbol::SEQUENCE,
        [decimal(a.clone() - whole.clone()), decimal(whole)],
    )))
}

/// `MultinomialProbability[n, p, k]`, giac's `multinomial`.
pub(super) fn multinomial_probability(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, p, k] = e.args() else {
        return Ok(None);
    };
    if !is_list(p) || !is_list(k) {
        return Ok(None);
    }
    if p.args().len() != k.args().len() {
        return refuse(ev, e, "as many probabilities as counts are needed");
    }
    let Some(counts) = k
        .args()
        .iter()
        .map(Expr::to_integer)
        .collect::<Option<Vec<Integer>>>()
        .filter(|ks| ks.iter().all(|k| !k.is_negative()))
    else {
        return refuse(ev, e, "nonnegative integer counts are needed");
    };
    let total = counts.iter().fold(Integer::ZERO, |sum, k| &sum + k);
    if !is_zero(ev, &(n.clone() - Expr::from(total)))? {
        return refuse(ev, e, "the counts must sum to n");
    }
    let factorial = |x: &Expr| Expr::call(Symbol::FACTORIAL, [x.clone()]);
    let mut result = factorial(n);
    for (p, k) in p.args().iter().zip(k.args()) {
        result = result * p.clone().pow(k.clone()) / factorial(k);
    }
    Ok(Some(result))
}

/// `ModularOrder[a % n]`, giac's `order`.
pub(super) fn modular_order(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [modular] = e.args() else {
        return Ok(None);
    };
    if modular.head().and_then(Expr::as_symbol) != Some(Symbol::MODULAR) {
        return Ok(None);
    }
    let [a, n] = modular.args() else {
        return Ok(None);
    };
    let (Some(residue), Some(modulus)) = (a.to_integer(), n.to_integer()) else {
        return Ok(None);
    };
    if residue.gcd(&modulus) != Integer::ONE {
        return refuse(ev, e, "not invertible");
    }
    Ok(Some(Expr::call(
        Symbol::MULTIPLICATIVE_ORDER,
        [a.clone(), n.clone()],
    )))
}

/// `SignedFactorList[n]`, giac's `maple_ifactors`.
pub(super) fn signed_factor_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let Some(n) = x.to_integer() else {
        return refuse(ev, e, "an integer is needed");
    };
    let list = |items: Vec<Expr>| Expr::call(Symbol::LIST, items);
    let sign = Expr::from(n.sign() as i64);
    let magnitude = Integer::from(n.unsigned_abs());
    if n.is_zero() || magnitude == Integer::ONE {
        return Ok(Some(list(vec![sign, list(vec![])])));
    }
    let flat = ev.eval(&Expr::call(
        Symbol::INTEGER_FACTOR_LIST,
        [Expr::from(magnitude)],
    ))?;
    if !is_list(&flat) {
        return Ok(None);
    }
    let pairs = flat.args().chunks(2).map(|pe| list(pe.to_vec())).collect();
    Ok(Some(list(vec![sign, list(pairs)])))
}

/// `AverageRateOfChange[f, x]` and `[f, x, h]`, giac's `avgRC`.
pub(super) fn average_rate_of_change(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, x, h) = match e.args() {
        [f, x] => (
            f,
            x,
            ev.eval(&Expr::call(Symbol::N, [Expr::from(1) / Expr::from(1000)]))?,
        ),
        [f, x, h] => (f, x, h.clone()),
        _ => return Ok(None),
    };
    if x.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    if is_zero(ev, &h)? {
        return undefined();
    }
    let moved = f.substitute(&[(x.clone(), x.clone() + h.clone())]);
    Ok(Some((moved - f.clone()) / h))
}

/// `BoxConstraints[[x…], [[a, b]…]]`, giac's `box_constraints`.
pub(super) fn box_constraints(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [vars, bounds] = e.args() else {
        return Ok(None);
    };
    if !is_list(vars) || !is_list(bounds) {
        return Ok(None);
    }
    let pair = |b: &Expr| is_list(b) && b.args().len() == 2;
    if vars.args().len() != bounds.args().len() || !bounds.args().iter().all(pair) {
        return refuse(ev, e, "invalid list of bounds");
    }
    let equations = vars.args().iter().zip(bounds.args()).map(|(v, b)| {
        let range = Expr::call(Symbol::RANGE, b.args().iter().cloned());
        Expr::call(Symbol::EQUATION, [v.clone(), range])
    });
    Ok(Some(Expr::call(Symbol::SEQUENCE, equations)))
}

#[cfg(test)]
mod tests;
