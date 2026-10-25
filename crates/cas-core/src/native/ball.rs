//! Balls in the language (CAS-DOC-5 §6.4, §6.5): `Ball[c, r]`, its
//! `ball(c, r)`, the ball of midpoint `c` and radius `r`; `Decimal[n, e]`,
//! a decimal literal; and `Exact[e]`, its `exact`, back to fractions.

use std::cmp::Ordering;

use cas_num::ball::{Ball, Float, Mag, Round};
use cas_num::{Complex, Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::print::ball_decimal;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// An upper bound of `r >= 0` as a radius; `None` for a negative or
/// symbolic `r`.
fn radius(r: &Expr) -> Option<Mag> {
    if let View::Ball(b) = r.view() {
        return (b.sign() != Some(Ordering::Less)).then(|| Mag::from_float(b.mid()).add(b.rad()));
    }
    let q = r.to_rational()?;
    if q.is_negative() {
        return None;
    }
    let (up, _) = Float::from_rational(&q, 30, Round::Up);
    Some(Mag::from_float(&up))
}

/// `Ball[c, r]`: an exact real `c` rounded to the working precision, its
/// radius grown by `r`, which may be exact or a ball (its upper bound); a
/// ball `c` with `r` more radius. It threads over a list `c`, with a
/// scalar `r` for every element or a list `r` of the same length paired with
/// it, other shapes being [`EvalError::Shape`] (CAS-DOC-5 §6.5). Symbols, a
/// negative `r` and a complex `c` leave it unevaluated.
pub(super) fn ball(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [c, r] = e.args() else {
        return Ok(None);
    };
    if is_list(c) {
        let items: Vec<Expr> = if is_list(r) {
            if r.args().len() != c.args().len() {
                return Err(EvalError::Shape);
            }
            c.args()
                .iter()
                .zip(r.args())
                .map(|(c, r)| Expr::apply(&Symbol::BALL.into(), &[c.clone(), r.clone()]))
                .collect()
        } else {
            c.args()
                .iter()
                .map(|c| Expr::apply(&Symbol::BALL.into(), &[c.clone(), r.clone()]))
                .collect()
        };
        return Ok(Some(Expr::apply(&Symbol::LIST.into(), &items)));
    }
    if is_list(r) {
        return Err(EvalError::Shape);
    }
    let Some(r) = radius(r) else {
        return Ok(None);
    };
    let center = match c.view() {
        View::Ball(b) => b.clone(),
        _ => match c.to_rational() {
            Some(q) => Ball::from_rational(&q, ev.context().working_precision()),
            None => return Ok(None),
        },
    };
    Ok(Some(Expr::from(center.add_error(&r))))
}

/// `Decimal[n, e]`, the literal `n · 10^e` with the digits of `n` as
/// written: its ball at the working precision, or at more bits for a
/// literal with more digits than that precision reads back
/// ([`Ball::from_literal`], CAS-148). Exact when binary holds the value,
/// as for `0.5` and `2.0`.
pub(super) fn decimal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, exponent] = e.args() else {
        return Ok(None);
    };
    let (Some(n), Some(exponent)) = (
        n.to_integer(),
        exponent.to_integer().and_then(|x| x.to_i64()),
    ) else {
        return Ok(None);
    };
    let working = ev.context().working_precision();
    Ok(Some(Expr::from(Ball::from_literal(&n, exponent, working))))
}

/// How many digits fewer than its precision shows a fraction must match for
/// [`exact`]: 12 at 53 bits, giac's tolerance of `10^-12`.
const SLACK_DIGITS: u64 = 3;

/// `Exact[e]`: `e` with each ball the simplest fraction that matches the
/// digits it shows (D4, CAS-148, ours). Exact numbers and other parts stay.
///
/// - An exact ball is its midpoint: `exact(0.5)` is `1/2`.
/// - Otherwise the ball shows the decimal `s`, of `d` significant digits, at
///   the working precision (D2). Take the simplest fraction, of least
///   denominator, that rounds to `s` at `max(d, D - 3)` digits, `D` being
///   the digits the precision carries (15 at 53 bits). It is the result
///   when it is far from a coincidence: a fraction of denominator `q` falls
///   in an interval of half-width `w` by chance with a probability of about
///   `q² w`, so it must have `q² w <= 10^-3`. Otherwise `s` is. (No other
///   fraction of `s`'s denominator fits an interval that narrow, so the
///   result is `s` or simpler.)
/// - So `exact(0.1)` is `1/10` and `exact(0.333333333333)` is `1/3`, while
///   `exact(0.3333)` stays `3333/10000`, `exact(0.7)` stays `7/10`, and the
///   15 digits of `3.14159265358979` stay as they are: no fraction near them
///   is that much simpler.
/// - A ball that shows no digit is the simplest fraction inside it, 0 when it
///   holds 0.
/// - A complex ball, part by part.
///
/// A fraction above `Limits::number_bits` bits leaves the ball as it is,
/// and so does a ball that shows no digit, does not hold 0, and whose
/// midpoint or radius is above it.
pub(super) fn exact(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    let prec = ev.context().working_precision();
    let max_bits = ev.context().limits().number_bits;
    Ok(Some(arg.replace(|part| match part.view() {
        View::Ball(b) => exact_ball(b, prec, max_bits).map(Expr::from),
        View::ComplexBall(z) => {
            let (re, im) = (
                exact_ball(z.re(), prec, max_bits)?,
                exact_ball(z.im(), prec, max_bits)?,
            );
            Some(Expr::from(Complex::new(re, im)))
        }
        _ => None,
    })))
}

/// The fraction [`exact`] makes of `b`; `None` when it would take more than
/// `max_bits` bits.
///
/// The midpoint and the radius become fractions only within `max_bits`
/// bits, which their exponents tell first: a ball near `2^(-10^9)`, as
/// `(4/9)^Ball[9×10^8, 1]` is, would otherwise make fractions of 10^9 bits
/// to refuse them (CAS-170).
fn exact_ball(b: &Ball, prec: u64, max_bits: u64) -> Option<Rational> {
    let fits = |q: &Rational| q.numer().bit_length().max(q.denom().bit_length()) <= max_bits;
    if b.is_exact() {
        return b.mid().to_rational_within(max_bits);
    }
    let Some((negative, digits, exponent)) = shown_decimal(&ball_decimal(b, prec)) else {
        let rad = b.rad().to_float()?;
        if b.sign().is_none() {
            // The ball holds 0, the simplest fraction of all.
            return Some(Rational::from(0));
        }
        let mid = b.mid().to_rational_within(max_bits)?;
        let rad = rad.to_rational_within(max_bits)?;
        return Rational::simplest_between(&(&mid - &rad), &(&mid + &rad)).filter(fits);
    };
    let accuracy = u64::try_from(b.rel_accuracy_bits()).unwrap_or(0);
    let carried = prec.max(accuracy) * 30_103 / 100_000;
    let shown = u64::try_from(digits.to_string().len()).ok()?;
    let matched = shown.max(carried.saturating_sub(SLACK_DIGITS)).max(1);
    // The decimal exponent of the first digit, which also bounds the size
    // of the result.
    let e10 = i64::try_from(shown).ok()? - 1 + exponent;
    if e10
        .unsigned_abs()
        .max(exponent.unsigned_abs())
        .saturating_mul(3322)
        / 1000
        > max_bits
    {
        return None;
    }
    let power = |k: i64| {
        let p = Rational::from(Integer::from(10u64).pow(k.unsigned_abs()));
        if k >= 0 { p } else { p.recip() }
    };
    let value = &Rational::from(if negative { -digits } else { digits }) * &power(exponent);
    // Half a unit in the `matched`-th digit.
    let half = &power(e10 - i64::try_from(matched).ok()? + 1) / &Rational::from(2);
    let simplest = Rational::simplest_between(&(&value - &half), &(&value + &half))?;
    let q = Rational::from(simplest.denom().clone());
    let no_chance = &(&q * &q) * &half <= Rational::new(Integer::from(1), Integer::from(1000));
    let chosen = if no_chance { simplest } else { value };
    Some(chosen).filter(fits)
}

/// The decimal `ball_decimal` writes, `-1.25e-5` as `(true, 125, -7)`: its
/// sign, its digits without leading zeros, and the power of 10 they are
/// scaled by. `None` for a ball that shows no digit (`[± 1.2e-16]`).
pub(super) fn shown_decimal(s: &str) -> Option<(bool, Integer, i64)> {
    let (mantissa, exponent) = match s.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i64>().ok()?),
        None => (s, 0),
    };
    let (negative, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => (true, m),
        None => (false, mantissa),
    };
    let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if int.is_empty() || !int.bytes().chain(frac.bytes()).all(|c| c.is_ascii_digit()) {
        return None;
    }
    let all = format!("{int}{frac}");
    let digits = all.trim_start_matches('0');
    let digits = Integer::from_str_radix(if digits.is_empty() { "0" } else { digits }, 10).ok()?;
    let shift = i64::try_from(frac.len()).ok()?;
    Some((negative, digits, exponent - shift))
}

#[cfg(test)]
mod tests;
