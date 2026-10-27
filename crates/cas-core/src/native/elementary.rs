//! The elementary functions of balls (CAS-DOC-5 §9.2, CAS-84): `Exp`,
//! `Log`, `Sin`, `Cos`, `Tan`, `ArcSin`, `ArcCos` and `ArcTan` of a ball,
//! real or complex, are the balls of their values, as arithmetic makes balls
//! of exact numbers met with balls (§9.1). The real function where it has a
//! value, else the principal complex one, as `Power` does; where neither has
//! one, `Indeterminate` with a warning, and `-Infinity` for the logarithm of
//! `0.0`. Exact arguments have their special values first (CAS-200,
//! `special`), `Sin[Pi/6]` is `1/2`; others stay as they are.

use cas_num::ball::{Ball, ComplexBall};

use super::number::{Infinity, Num, Precision};
use super::special;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `f[x]` at a special value of an exact `x` ([`special::value`]), of an
/// inverse function's value or the other way ([`special::composition`]), or for a
/// ball `x`: `real` of a real ball, `complex` of a complex one or where the
/// real function has no value.
fn elementary(
    ev: &mut Eval<'_>,
    e: &Expr,
    real: fn(&Ball, u64) -> Ball,
    complex: fn(&ComplexBall, u64) -> ComplexBall,
) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if let Some(f) = e.head().and_then(Expr::as_symbol)
        && let Some(v) = special::value(f, x).or_else(|| special::composition(ev.context(), f, x))
    {
        return Ok(Some(v));
    }
    let Some(n) = Num::of(x).filter(Num::is_ball) else {
        return Ok(None);
    };
    let working = Precision::of(ev.context());
    let value = n
        .ball_function(working, real, complex)
        .expect("a ball")
        .to_expr();
    if value.as_symbol() == Some(Symbol::INDETERMINATE) {
        ev.warn(Warning::Indeterminate(e.clone()));
    }
    Ok(Some(value))
}

pub(super) fn exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::exp, ComplexBall::exp)
}

/// `Log`, and `-Infinity` for `Log[0.0]`, as for an exact 0 the limit is.
pub(super) fn log(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if let [x] = e.args() {
        if let Some(n @ Num::Ball(_)) = Num::of(x) {
            if n.is_zero() {
                ev.warn(Warning::Infinite(e.clone()));
                return Ok(Some(Infinity::Directed(Num::one().negated()).to_expr()));
            }
        }
    }
    elementary(ev, e, Ball::log, ComplexBall::log)
}

pub(super) fn sin(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::sin, ComplexBall::sin)
}

pub(super) fn cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::cos, ComplexBall::cos)
}

pub(super) fn tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::tan, ComplexBall::tan)
}

pub(super) fn arc_sin(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::asin, ComplexBall::asin)
}

pub(super) fn arc_cos(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::acos, ComplexBall::acos)
}

pub(super) fn arc_tan(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::atan, ComplexBall::atan)
}
