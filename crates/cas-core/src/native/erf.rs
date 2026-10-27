//! `Erf` and `Erfc` (CAS-DOC-15 §8, CAS-247), Xcas's `erf` and `erfc`: the
//! error function and its complement `1 − Erf[x]`, which giac writes as
//! that difference and this CAS keeps, so that a ball's `Erfc` keeps its
//! relative accuracy however small.
//!
//! - `Erf[0] = 0` and `Erfc[0] = 1`; `Erf[±Infinity] = ±1`, `Erfc[Infinity]
//!   = 0` and `Erfc[-Infinity] = 2`.
//! - A real ball by `Ball::erf` and `Ball::erfc`.
//! - The sign of a negative number or coefficient comes out: `Erf[−x] =
//!   −Erf[x]` and `Erfc[−x] = 2 − Erfc[x]`.
//! - Anything else stays.

use super::number::{Infinity, Num, Precision};
use super::special::negated;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `Erf[x]` and `Erfc[x]`.
pub(super) fn erf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(f), [x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let complement = f == Symbol::ERFC;
    if *x == Expr::from(0) {
        return Ok(Some(Expr::from(i64::from(complement))));
    }
    if let Some(Infinity::Directed(d)) = Infinity::of(x) {
        let sign = if d.is_one() {
            1
        } else if d.negated().is_one() {
            -1
        } else {
            return Ok(None);
        };
        return Ok(Some(Expr::from(if complement { 1 - sign } else { sign })));
    }
    if let Some(ref number @ Num::Ball(ref ball)) = Num::of(x) {
        let bits = number.precision(Precision::of(ev.context()));
        let value = if complement {
            ball.erfc(bits + 16)
        } else {
            ball.erf(bits + 16)
        };
        return Ok(Some(
            Num::Ball(Box::new(value.round_to_accuracy(bits))).to_expr(),
        ));
    }
    let Some(y) = negated(x) else {
        return Ok(None);
    };
    let value = Expr::call(f, [y]);
    Ok(Some(if complement {
        Expr::from(2) - value
    } else {
        -value
    }))
}
