//! The elementary functions of balls (CAS-DOC-5 §9.2, CAS-84): `Exp`,
//! `Log`, `Sin`, `Cos`, `Tan`, `ArcSin`, `ArcCos` and `ArcTan` of a ball,
//! and the hyperbolic `Sinh` … `ArcTanh` (CAS-DOC-15 §2, CAS-238),
//! real or complex, are the balls of their values, as arithmetic makes balls
//! of exact numbers met with balls (§9.1). The real function where it has a
//! value, else the principal complex one, as `Power` does; where neither has
//! one, `Indeterminate` with a warning, and `-Infinity` for the logarithm of
//! `0.0`. Exact arguments have their special values first (CAS-200,
//! `special`), `Sin[Pi/6]` is `1/2`; others stay as they are.
//!
//! Xcas's logarithms to a base and angle conversions, rewritten exactly,
//! listable:
//!
//! - `LogB[x, b]` (`logb`, CAS-565), giac's: `ln(x)/ln(b)`, and for an
//!   integer `x ≠ 0` and an integer `b > 1` the `n` factors `b` of `x`
//!   taken out first, `ln(x/b^n)/ln(b) + n`, so `logb(24, 2)` is
//!   `ln(3)/ln(2) + 3`; `LogB[x]` is `Log[x]`. `Log2[x]` (`log2`) and
//!   `Log10[x]` (`log10`) are `LogB[x, 2]` and `LogB[x, 10]`; `Exp10[x]`
//!   (`alog10`) is `10^x`.
//! - `Radians[x]` (`radians`, CAS-564) is `π x/180` and `Degrees[x]`
//!   (`degrees`) `180 x/π`, exact where giac's are decimals.

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Integer, Rational};

use super::number::{Infinity, Num, Precision};
use super::special;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `f[x]` at a special value of an exact `x` ([`special::value`]), of an
/// inverse function's value or the other way ([`special::composition`]), at
/// integer multiples of `π` ([`special::periodic`]), or for a ball `x`:
/// `real` of a real ball, `complex` of a complex one or where the real
/// function has no value.
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
        && let Some(v) = special::value(f, x)
            .or_else(|| special::composition(ev.context(), f, x))
            .or_else(|| special::periodic(ev.context(), f, x))
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
    // `Exp[k W(z)] = (z/W(z))^k` (CAS-DOC-31 §3.3).
    if let [x] = e.args()
        && let Some(value) = super::lambert::exponential(x)
    {
        return Ok(Some(value));
    }
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

pub(super) fn sinh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::sinh, ComplexBall::sinh)
}

pub(super) fn cosh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::cosh, ComplexBall::cosh)
}

pub(super) fn tanh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::tanh, ComplexBall::tanh)
}

pub(super) fn arc_sinh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::asinh, ComplexBall::asinh)
}

pub(super) fn arc_cosh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::acosh, ComplexBall::acosh)
}

pub(super) fn arc_tanh(ev: &mut Eval<'_>, e: &Expr) -> Native {
    elementary(ev, e, Ball::atanh, ComplexBall::atanh)
}

/// The reciprocal trigonometric functions, rewritten as giac does
/// (CAS-DOC-15 §3, CAS-239): `Sec[x] = 1/Cos[x]`, `Csc[x] = 1/Sin[x]`,
/// `Cot[x] = Cos[x]/Sin[x]`, `ArcSec[x] = ArcCos[1/x]`, `ArcCsc[x] =
/// ArcSin[1/x]` and `ArcCot[x] = ArcTan[1/x]`, `ArcCot[0] = Pi/2`.
pub(super) fn reciprocal(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (Some(f), [x]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return Ok(None);
    };
    let call = |g: Symbol, u: Expr| Expr::call(g, [u]);
    let inverse = || x.pow(-1);
    Ok(Some(match f {
        Symbol::SEC => call(Symbol::COS, x.clone()).pow(-1),
        Symbol::CSC => call(Symbol::SIN, x.clone()).pow(-1),
        Symbol::COT => call(Symbol::COS, x.clone()) * call(Symbol::SIN, x.clone()).pow(-1),
        Symbol::ARC_SEC => call(Symbol::ARC_COS, inverse()),
        Symbol::ARC_CSC => call(Symbol::ARC_SIN, inverse()),
        Symbol::ARC_COT if *x == Expr::from(0) => {
            Expr::from(Rational::new(Integer::ONE, Integer::from(2))) * Expr::from(Symbol::PI)
        }
        Symbol::ARC_COT => call(Symbol::ARC_TAN, inverse()),
        _ => return Ok(None),
    }))
}

/// `LogB[x, b]`, `LogB[x]`, `Log2[x]` and `Log10[x]`.
pub(super) fn log_b(_: &mut Eval<'_>, e: &Expr) -> Native {
    let base = match e.head().and_then(Expr::as_symbol) {
        Some(Symbol::LOG2) => Some(Expr::from(2)),
        Some(Symbol::LOG10) => Some(Expr::from(10)),
        _ => None,
    };
    let (x, b) = match (base, e.args()) {
        (Some(b), [x]) => (x.clone(), b),
        (None, [x]) => return Ok(Some(Expr::call(Symbol::LOG, [x.clone()]))),
        (None, [x, b]) => (x.clone(), b.clone()),
        _ => return Ok(None),
    };
    let (mut x, mut n) = (x, 0_u64);
    if let (Some(mut m), Some(base)) = (x.to_integer(), b.to_integer())
        && base > Integer::ONE
        && !m.is_zero()
    {
        while m.is_divisible_by(&base) {
            m = m.div_exact(&base);
            n += 1;
        }
        x = Expr::from(m);
    }
    let log = |u: Expr| Expr::call(Symbol::LOG, [u]);
    Ok(Some(log(x) / log(b) + Expr::from(Integer::from(n))))
}

/// `Exp10[x]`: `10^x`.
pub(super) fn exp10(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(Expr::from(10).pow(x.clone())))
}

/// `Radians[x]` and `Degrees[x]`.
pub(super) fn angle(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let ratio = Expr::from(Rational::new(Integer::ONE, Integer::from(180))) * Symbol::PI;
    Ok(Some(match e.head().and_then(Expr::as_symbol) {
        Some(Symbol::RADIANS) => x * ratio,
        Some(Symbol::DEGREES) => x / ratio,
        _ => return Ok(None),
    }))
}
