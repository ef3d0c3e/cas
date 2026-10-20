//! `Factorial`, the language's `n!` (CAS-DOC-7 §9).

use cas_num::{Integer, Natural, Rational};

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
/// Anything else stays, waiting for the Gamma function.
pub(super) fn factorial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
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
