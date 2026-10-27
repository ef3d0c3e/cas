//! Number theory (CAS-DOC-17): Xcas's integer functions over cas-num's
//! kernels.
//!
//! - `PrimeQ[n]` (`isprime`, CAS-267): `True` for a prime, proven below 2^64
//!   and a Baillie–PSW probable prime above, `False` for every other integer
//!   and for a number that is not one; a symbol stays.
//! - `NextPrime[n]` (`nextprime`): the least prime above `n`, 2 below 2.
//!   `PreviousPrime[n]` (`prevprime`): the greatest below `n`; `n ≤ 2` stays,
//!   with a warning (giac's `prevprime(2)` is an error).

use cas_num::{Integer, Natural};

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

#[cfg(test)]
mod tests;
