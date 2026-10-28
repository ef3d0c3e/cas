//! `Factorial`, the language's `n!` (CAS-DOC-7 §9).

use cas_num::{Integer, Natural, Rational};

use super::number::{Infinity, Num, Precision};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// An upper bound on the bits of `n!`: `n` factors of at most the bits of
/// `n` each. Within a sixth of the truth from `n = 10^5` on.
fn factorial_bits(n: u64) -> u64 {
    n.saturating_mul(u64::from(u64::BITS - n.leading_zeros()))
}

/// `Factorial[n]`, `Listable`, exactly, as in Xcas:
///
/// - an integer from 0: `n!` ([`Natural::factorial`]);
/// - a half-integer: `Γ(n + 1)`, a rational times `√π`, such as
///   `(1/2)! = √π/2` and `(-1/2)! = √π`;
/// - a negative integer: a pole of `Γ`, `ComplexInfinity`, with a warning,
///   as `1/0`;
/// - a result estimated above `Limits::number_bits` bits stays
///   unevaluated, as a large power of numbers does.
///
/// A real ball's factorial is `Γ(x + 1)` (CAS-DOC-15 §6, CAS-243), so that
/// `N` of any rational's has a value. Anything else stays.
pub(super) fn factorial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if matches!(Num::of(arg), Some(Num::Ball(_))) {
        let shifted = ev.eval(&(arg.clone() + Expr::from(1)))?;
        return Ok(ball_gamma(ev, e, &shifted));
    }
    let Some(value) = arg.to_rational() else {
        return Ok(None);
    };
    let max_bits = ev.context().limits().number_bits;
    if value.is_integer() {
        if value.is_negative() {
            ev.warn(Warning::Infinite(e.clone()));
            return Ok(Some(Symbol::COMPLEX_INFINITY.into()));
        }
        let Ok(n) = u64::try_from(value.numer()) else {
            return Ok(None);
        };
        if factorial_bits(n) > max_bits {
            return Ok(None);
        }
        return Ok(Some(Expr::from(Integer::from(Natural::factorial(n)))));
    }
    if value.denom().to_i64() != Some(2) {
        return Ok(None);
    }
    // `value + 1 = k + 1/2`, and `Γ(k + 1/2) = c √π` with
    // `c = (2k)! / (4^k k!)` for `k >= 0`, `(-4)^j j! / (2j)!` for `k = -j`.
    let Some(numer) = value.numer().to_i64() else {
        return Ok(None);
    };
    // `(numer + 1) / 2`, for an odd `numer`.
    let k = numer.div_euclid(2) + 1;
    let j = k.unsigned_abs();
    if factorial_bits(2 * j) > max_bits {
        return Ok(None);
    }
    let ratio = Rational::new(
        Integer::from(Natural::factorial(2 * j)),
        Integer::from(Natural::factorial(j) << (2 * j)),
    );
    let coefficient = if k >= 0 {
        ratio
    } else if j % 2 == 0 {
        ratio.recip()
    } else {
        -ratio.recip()
    };
    let half = Expr::from(Rational::new(1.into(), 2.into()));
    let sqrt_pi = Expr::apply(&Symbol::POWER.into(), &[Symbol::PI.into(), half]);
    Ok(Some(Expr::apply(
        &Symbol::TIMES.into(),
        &[Expr::from(coefficient), sqrt_pi],
    )))
}

/// The most factors `Binomial[n, k]` multiplies out for a symbolic `n`.
const MAX_BINOMIAL_FACTORS: u64 = 1000;

/// `Binomial[n, k]`, Xcas's `binomial` and `comb` (CAS-DOC-13 §1,
/// CAS-225): for an integer `k`, `n (n − 1) ⋯ (n − k + 1) / k!`, 0 for a
/// negative `k`:
///
/// - a rational `n`, exactly: `binomial(5,2)` is `10`, and the product
///   vanishes for an integer `0 ≤ n < k`; an estimate past
///   `Limits::number_bits` stays;
/// - another `n`, the product built, up to 1000 factors, as giac's
///   `binomial(n,2)` is `n!/(2*(n-2)!)`, here `n (n − 1)/2`.
///
/// Anything else stays.
pub(super) fn binomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, k] = e.args() else {
        // Xcas's `binomial(n, p, k)`, the binomial distribution's mass.
        return super::statistics::distributions::discrete::binomial_mass(ev, e);
    };
    let Some(k) = k.as_i64() else {
        return Ok(None);
    };
    if k < 0 {
        return Ok(Some(Expr::from(0)));
    }
    let k = k.unsigned_abs();
    if let Some(n) = n.to_rational() {
        // About k times the bits of n, and those of k! over it.
        let bits = n.numer().bit_length().max(n.denom().bit_length()).max(1);
        if k.saturating_mul(bits) > ev.context().limits().number_bits {
            return Ok(None);
        }
        let mut value = Rational::ONE;
        for i in 0..k {
            let factor =
                &(&n - &Rational::from(Integer::from(i))) / &Rational::from(Integer::from(i + 1));
            value = &value * &factor;
            if value.is_zero() {
                break;
            }
        }
        return Ok(Some(Expr::from(value)));
    }
    if k > MAX_BINOMIAL_FACTORS {
        return Ok(None);
    }
    let factors: Vec<Expr> = (0..k)
        .map(|i| n.clone() - Expr::from(i64::try_from(i).expect("a small count")))
        .collect();
    let denominator = Expr::from(Integer::from(Natural::factorial(k)));
    Ok(Some(
        Expr::apply(&Symbol::TIMES.into(), &factors) / denominator,
    ))
}

/// `Γ` of a real ball at its precision, `Indeterminate` with a warning
/// where the ball meets a pole; `None` for another argument.
fn ball_gamma(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Option<Expr> {
    let number = Num::of(x)?;
    let Num::Ball(ball) = &number else {
        return None;
    };
    let bits = number.precision(Precision::of(ev.context()));
    let value = ball.gamma(bits + 16).round_to_accuracy(bits);
    if !value.is_finite() {
        ev.warn(Warning::Indeterminate(e.clone()));
        return Some(Symbol::INDETERMINATE.into());
    }
    Some(Num::Ball(Box::new(value)).to_expr())
}

/// `Gamma[z]`, Xcas's `Gamma` (CAS-DOC-15 §6, CAS-243): `(z − 1)!` at the
/// positive integers and the half-integers (`Gamma[1/2] = √π`), through
/// `Factorial`; `ComplexInfinity` with a warning at the poles, the integers
/// `≤ 0`; a real ball's by `Ball::gamma`. Anything else stays.
pub(super) fn gamma(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [z] = e.args() else {
        return Ok(None);
    };
    if let Some(q) = z.to_rational() {
        if q.is_integer() && !q.is_positive() {
            ev.warn(Warning::Infinite(e.clone()));
            return Ok(Some(Symbol::COMPLEX_INFINITY.into()));
        }
        if q.is_integer() || q.denom().to_i64() == Some(2) {
            let below = &q - &Rational::ONE;
            return Ok(Some(Expr::call(Symbol::FACTORIAL, [Expr::from(below)])));
        }
        return Ok(None);
    }
    Ok(ball_gamma(ev, e, z))
}

/// `Beta[a, b]`, Xcas's `Beta` (CAS-DOC-15 §6, CAS-243): `Gamma[a]
/// Gamma[b] / Gamma[a + b]` where each has a value, finite; otherwise it
/// stays.
pub(super) fn beta(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let mut values = Vec::with_capacity(3);
    for x in [a.clone(), b.clone(), a.clone() + b.clone()] {
        let g = ev.eval(&Expr::call(Symbol::GAMMA, [x]))?;
        let stays = g.contains(&Expr::from(Symbol::GAMMA))
            || [
                Symbol::COMPLEX_INFINITY,
                Symbol::INDETERMINATE,
                Symbol::INFINITY,
            ]
            .into_iter()
            .any(|s| g.contains(&Expr::from(s)))
            || Infinity::of(&g).is_some();
        if stays {
            return Ok(None);
        }
        values.push(g);
    }
    let [ga, gb, gab]: [Expr; 3] = values.try_into().expect("three values");
    Ok(Some(ga * gb / gab))
}
