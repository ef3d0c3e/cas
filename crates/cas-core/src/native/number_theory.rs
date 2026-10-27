//! Number theory (CAS-DOC-17): Xcas's integer functions over cas-num's
//! kernels.
//!
//! - `PrimeQ[n]` (`isprime`, CAS-267): `True` for a prime, proven below 2^64
//!   and a Baillie–PSW probable prime above, `False` for every other integer
//!   and for a number that is not one; a symbol stays.
//! - `NextPrime[n]` (`nextprime`): the least prime above `n`, 2 below 2.
//!   `PreviousPrime[n]` (`prevprime`): the greatest below `n`; `n ≤ 2` stays,
//!   with a warning (giac's `prevprime(2)` is an error).
//! - Factorisation (CAS-268), by trial division, perfect powers and
//!   Pollard's rho, whose effort shrinks as `(128/bits)²` above 128 bits so
//!   that a call stays within about a second; a composite it leaves
//!   unsplit leaves the call, with a warning. `IntegerFactor[n]`
//!   (`ifactor`) is the product that does not multiply out, `Hold[2^3 3^2
//!   5]` for 360 (giac's `2^3*3^2*5` stays, a product of integers evaluates
//!   here), the sign in front; a prime, 0 and ±1 are themselves.
//!   `IntegerFactorList[n]` (`ifactors`): giac's flat list `[p₁, e₁, …]`,
//!   `[-1, 1, …]` first for a negative `n`. `Divisors[n]` (`idivis`): the
//!   positive divisors of `|n|`, ascending, at most 2^20 of them.
//!   `EulerPhi[n]` (`euler`): `n Π (1 − 1/p)`, `φ(0) = 0`; negatives stay.

use cas_num::{Factorization, Integer, Natural};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The one argument of `e`, an integer.
fn integer_argument(e: &Expr) -> Option<Integer> {
    match e.args() {
        [n] => n.to_integer(),
        _ => None,
    }
}

/// `n`'s magnitude.
fn magnitude(n: &Integer) -> Natural {
    n.clone().into_parts().1
}

/// The call stays, with giac's error as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// `PrimeQ[n]`.
pub(super) fn prime_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.to_integer() else {
        return Ok(n.is_number().then(|| Symbol::FALSE.into()));
    };
    let prime = !n.is_negative() && magnitude(&n).primality().is_probable_prime();
    Ok(Some(
        if prime { Symbol::TRUE } else { Symbol::FALSE }.into(),
    ))
}

/// `NextPrime[n]`.
pub(super) fn next_prime(_: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    let p = if n.is_negative() {
        Natural::from(2u64)
    } else {
        magnitude(&n).next_prime()
    };
    Ok(Some(Expr::from(Integer::from(p))))
}

/// `PreviousPrime[n]`.
pub(super) fn previous_prime(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    let found = if n.is_negative() {
        None
    } else {
        magnitude(&n).previous_prime()
    };
    match found {
        Some(p) => Ok(Some(Expr::from(Integer::from(p)))),
        None => refuse(ev, e, "no prime below"),
    }
}

/// Pollard rho's iterations per composite up to 128 bits, cas-num's
/// default, enough for factors up to about 2^40.
const RHO_ITERATIONS: u64 = 1 << 22;

/// The most divisors `Divisors` lists.
const MAX_DIVISORS: u64 = 1 << 20;

/// The factorisation of `|n|`, `n ≠ 0`, with an effort that shrinks with
/// its size; `None` when a composite is left unsplit.
fn factorization(n: &Integer) -> Option<Factorization> {
    let m = magnitude(n);
    let bits = u128::from(m.bit_length().max(1));
    let effort = if bits <= 128 {
        RHO_ITERATIONS
    } else {
        let scaled = u128::from(RHO_ITERATIONS) * 128 * 128 / (bits * bits);
        u64::try_from(scaled).unwrap_or(RHO_ITERATIONS).max(1 << 12)
    };
    let f = m.factor_with_effort(effort);
    f.is_complete().then_some(f)
}

/// `IntegerFactor[n]`.
pub(super) fn integer_factor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    if n.is_zero() || magnitude(&n).is_one() {
        return Ok(Some(Expr::from(n)));
    }
    let Some(f) = factorization(&n) else {
        return refuse(ev, e, "a composite factor is not split");
    };
    if let [only] = f.factors()
        && only.exponent() == 1
    {
        return Ok(Some(Expr::from(n)));
    }
    let mut factors = Vec::with_capacity(f.factors().len() + 1);
    if n.is_negative() {
        factors.push(Expr::from(-1));
    }
    for x in f.factors() {
        let p = Expr::from(Integer::from(x.base().clone()));
        factors.push(if x.exponent() == 1 {
            p
        } else {
            Expr::call(Symbol::POWER, [p, Expr::from(Integer::from(x.exponent()))])
        });
    }
    let product = match factors.len() {
        1 => factors.pop().expect("one"),
        _ => Expr::call(Symbol::TIMES, factors),
    };
    Ok(Some(Expr::call(Symbol::HOLD, [product])))
}

/// `IntegerFactorList[n]`.
pub(super) fn integer_factor_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    if n.is_zero() {
        return refuse(ev, e, "factorisation of zero");
    }
    let Some(f) = factorization(&n) else {
        return refuse(ev, e, "a composite factor is not split");
    };
    let mut items = Vec::with_capacity(2 * f.factors().len() + 2);
    if n.is_negative() {
        items.extend([Expr::from(-1), Expr::from(1)]);
    }
    for x in f.factors() {
        items.push(Expr::from(Integer::from(x.base().clone())));
        items.push(Expr::from(Integer::from(x.exponent())));
    }
    Ok(Some(Expr::call(Symbol::LIST, items)))
}

/// `Divisors[n]`.
pub(super) fn divisors(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    if n.is_zero() {
        return refuse(ev, e, "every integer divides 0");
    }
    let Some(f) = factorization(&n) else {
        return refuse(ev, e, "a composite factor is not split");
    };
    let count = f
        .factors()
        .iter()
        .try_fold(1u64, |c, x| c.checked_mul(x.exponent().saturating_add(1)));
    if count.is_none_or(|c| c > MAX_DIVISORS) {
        return refuse(ev, e, "too many divisors to list");
    }
    let mut all = vec![Natural::one()];
    for x in f.factors() {
        let mut next =
            Vec::with_capacity(all.len() * usize::try_from(x.exponent() + 1).unwrap_or(1));
        for d in &all {
            let mut power = d.clone();
            next.push(power.clone());
            for _ in 0..x.exponent() {
                power = &power * x.base();
                next.push(power.clone());
            }
        }
        all = next;
    }
    all.sort();
    Ok(Some(Expr::call(
        Symbol::LIST,
        all.into_iter().map(|d| Expr::from(Integer::from(d))),
    )))
}

/// `EulerPhi[n]`.
pub(super) fn euler_phi(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(n) = integer_argument(e) else {
        return Ok(None);
    };
    if n.is_negative() {
        return refuse(ev, e, "Euler's totient of a negative integer");
    }
    if n.is_zero() || n == Integer::ONE {
        return Ok(Some(Expr::from(n)));
    }
    let Some(f) = factorization(&n) else {
        return refuse(ev, e, "a composite factor is not split");
    };
    let phi = f.factors().iter().fold(Natural::one(), |acc, x| {
        acc * &x.base().pow(x.exponent() - 1) * &(x.base() - 1u64)
    });
    Ok(Some(Expr::from(Integer::from(phi))))
}

#[cfg(test)]
mod tests;
