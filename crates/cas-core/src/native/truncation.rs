//! `IntegerPart` (`iPart`), `FractionalPart` (`fPart`, `frac`) and
//! `Truncate` (`trunc`) (CAS-622), Xcas's (giac ti89.cc `_iPart` and
//! `fPart`, maple.cc `_trunc`), on [`rounding`](super::rounding)'s `decide`:
//! a rational exactly, a ball by its ends, a numeric expression by its
//! approximation until the ends agree.
//!
//! - Toward zero: `iPart(-7/2)` is −3, `trunc(-4.7)` −4, and `fPart(x)` is
//!   `x − iPart(x)`, so `fPart(-1/2)` is −1/2.
//! - Exact on exact numbers, where giac approximates: `iPart(7/2)` is 3
//!   (giac 3.0), `fPart(sqrt(2))` `sqrt(2) − 1`, `trunc(sqrt(2), 3)` 707/500
//!   (giac 1.414); a decimal gives a decimal for `iPart`, `fPart` and
//!   `trunc(x, n)`, and `trunc(x)` an integer, as giac.
//! - `trunc(x, n)`, `n` an integer: `x` to `b = n` decimals for `n ≥ 0`;
//!   for `n < 0` giac's `b = −1 − n − floor(log10 |x|)`, `−n` significant
//!   digits (`trunc(12345.678, -2)` is 12000.0); toward zero where giac
//!   floors (`trunc(-sqrt(2), 3)`: giac −1.415, ours −1.414); a `10^b`
//!   past `Limits::number_bits` stays.
//! - A non-real number part by part; `iPart(x, y)` and `fPart(x, y)` the
//!   sequence of each, as giac; lists thread (`Listable`); what does not
//!   decide stays, free symbols too (giac's `fPart(x)` is
//!   `x-sign(x)*floor(x*sign(x))`, its `trunc(x)` `floor(x)`).

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::rounding::decide;
use crate::build;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// Which part.
#[derive(Clone, Copy)]
enum Part {
    Integer,
    Fractional,
    Truncate,
    /// `trunc(x, n)`.
    Digits(i64),
}

/// `10^b`.
fn power_of_ten(b: i64) -> Rational {
    let p = Rational::from(Integer::from(10_u64).pow(b.unsigned_abs()));
    if b < 0 { p.recip() } else { p }
}

/// `floor(log10 q)` of a positive `q`: an estimate from the bit lengths,
/// within 2 of it, then corrected.
fn floor_log10(q: &Rational) -> Option<i64> {
    let bits = |n: &Integer| i64::try_from(n.bit_length()).ok();
    let mut e = (bits(q.numer())? - bits(q.denom())?) * 30103 / 100_000;
    while power_of_ten(e) > *q {
        e -= 1;
    }
    while power_of_ten(e + 1) <= *q {
        e += 1;
    }
    Some(e)
}

/// `floor(log10 |x|)`, where the ends of `x`'s value agree on it.
fn magnitude(ctx: &Context, x: &Expr) -> Option<i64> {
    decide(ctx, x, |lo, hi| {
        let same_sign =
            (lo.is_positive() && hi.is_positive()) || (lo.is_negative() && hi.is_negative());
        if !same_sign {
            return None;
        }
        let (a, b) = (floor_log10(&lo.abs())?, floor_log10(&hi.abs())?);
        (a == b).then_some(a)
    })
}

/// `x s` truncated toward zero, over `s`, where the ends of `x`'s value
/// agree on it.
fn truncated(ctx: &Context, x: &Expr, s: &Rational) -> Option<Rational> {
    let t = decide(ctx, x, |lo, hi| {
        let (a, b) = ((lo * s).trunc(), (hi * s).trunc());
        (a == b).then_some(a)
    })?;
    Some(&Rational::from(t) / s)
}

/// `part` of a real `x`; `None` where its value does not decide.
fn real(ev: &mut Eval<'_>, part: Part, x: &Expr) -> Native {
    let ctx = ev.context();
    let decimal = matches!(x.view(), View::Ball(_));
    let ball = |q: &Rational| Expr::from(Ball::from_rational(q, ctx.working_precision()));
    let b = match part {
        Part::Digits(n) if n >= 0 => n,
        Part::Digits(n) => {
            let Some(b) = magnitude(ctx, x).and_then(|m| (-1_i64).checked_sub(n)?.checked_sub(m))
            else {
                return Ok(None);
            };
            b
        }
        _ => 0,
    };
    // `10^b` has about `3.33 |b|` bits.
    if b.unsigned_abs().saturating_mul(4) > ctx.limits().number_bits {
        return Ok(None);
    }
    let Some(q) = truncated(ctx, x, &power_of_ten(b)) else {
        return Ok(None);
    };
    Ok(Some(match part {
        Part::Integer | Part::Digits(_) if decimal => ball(&q),
        Part::Fractional => ev.eval(&(x.clone() - Expr::from(q)))?,
        _ => Expr::from(q),
    }))
}

/// `part` of `x`, a non-real number part by part.
fn part_of(ev: &mut Eval<'_>, part: Part, x: &Expr) -> Native {
    if let Some(r) = real(ev, part, x)? {
        return Ok(Some(r));
    }
    if !x.is_numeric() {
        return Ok(None);
    }
    let im = ev.eval(&Expr::call(Symbol::IM, [x.clone()]))?;
    if im.to_rational().is_some_and(|q| q.is_zero()) {
        return Ok(None);
    }
    let re = ev.eval(&Expr::call(Symbol::RE, [x.clone()]))?;
    let (Some(a), Some(b)) = (real(ev, part, &re)?, real(ev, part, &im)?) else {
        return Ok(None);
    };
    Ok(Some(ev.eval(&(a + build::i() * b))?))
}

/// `IntegerPart[x]` and `FractionalPart[x]`; of several arguments the
/// sequence of each.
pub(super) fn integer_part(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(f) = e.head().and_then(Expr::as_symbol) else {
        return Ok(None);
    };
    let part = if f == Symbol::INTEGER_PART {
        Part::Integer
    } else {
        Part::Fractional
    };
    match e.args() {
        [] => Ok(None),
        [x] => part_of(ev, part, x),
        args => Ok(Some(Expr::call(
            Symbol::SEQUENCE,
            args.iter().map(|a| Expr::call(f, [a.clone()])),
        ))),
    }
}

/// `Truncate[x]` and `Truncate[x, n]`.
pub(super) fn truncate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [x] => part_of(ev, Part::Truncate, x),
        [x, n] => {
            let View::Int(n) = n.view() else {
                return Ok(None);
            };
            part_of(ev, Part::Digits(n), x)
        }
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests;
