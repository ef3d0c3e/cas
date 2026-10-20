//! Integer division and parity (CAS-122): Xcas's `irem`, `iquo`, `even` and
//! `odd`, which procedures need (CAS-DOC-7 §8).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `Remainder[a, b]` and `Quotient[a, b]` of integers: Euclidean division,
/// `a = q b + r` with `0 <= r < |b|`, as giac for a positive `b`. A
/// remainder by 0 is `Indeterminate` and a quotient `ComplexInfinity`, with
/// a warning; anything else stays.
fn divide(ev: &mut Eval<'_>, e: &Expr, quotient: bool) -> Native {
    let [dividend, divisor] = e.args() else {
        return Ok(None);
    };
    let (Some(dividend), Some(divisor)) = (dividend.to_integer(), divisor.to_integer()) else {
        return Ok(None);
    };
    if divisor.is_zero() {
        return Ok(Some(if quotient {
            ev.warn(Warning::Infinite(e.clone()));
            Symbol::COMPLEX_INFINITY.into()
        } else {
            ev.warn(Warning::Indeterminate(e.clone()));
            Symbol::INDETERMINATE.into()
        }));
    }
    let (whole, rest) = dividend.div_rem_euclid(&divisor);
    Ok(Some(Expr::from(if quotient { whole } else { rest })))
}

pub(super) fn remainder(ev: &mut Eval<'_>, e: &Expr) -> Native {
    divide(ev, e, false)
}

pub(super) fn quotient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    divide(ev, e, true)
}

/// `EvenQ[n]` and `OddQ[n]`: for a number, whether it is an even or an odd
/// integer; other expressions stay, so a parity that is not known is not
/// taken as false.
fn parity(e: &Expr, odd: bool) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    if !n.is_number() {
        return Ok(None);
    }
    let holds = n
        .to_integer()
        .is_some_and(|n| n.unsigned_abs().bit(0) == odd);
    Ok(Some(Expr::from(if holds {
        Symbol::TRUE
    } else {
        Symbol::FALSE
    })))
}

pub(super) fn even_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    parity(e, false)
}

pub(super) fn odd_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    parity(e, true)
}
