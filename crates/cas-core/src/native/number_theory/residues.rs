//! Division, pseudoprimes and cyclotomic polynomials (CAS-DOC-39 §2,
//! CAS-422).
//!
//! - `QuotientRemainder[a, b]` (`iquorem`): `{q, r}` with `a = b q + r`, `0
//!   ≤ r < |b|`, `Quotient` and `Remainder` together: `iquorem(17,-5)` is
//!   `[-3,2]`. A divisor 0 stays, with a warning.
//! - `SymmetricMod[a, n]` (`smod`): `a mod n` in `(−n/2, n/2]`, `n` by its
//!   size: `smod(8,5)` is −2.
//! - `PseudoprimeQ[n]` (`is_pseudoprime`), `|n|`: 2 for a prime proven by
//!   cas-num (below 2^64), 1 for a probable prime that passed Baillie and
//!   PSW's test, 0 otherwise. giac answers 1 for `2^61-1`, which cas-num
//!   proves prime.
//! - `Cyclotomic[n]` (`cyclotomic`), `1 ≤ n ≤` [`MAX_CYCLOTOMIC`]: `Φ_n =
//!   Π_(d | n) (x^d − 1)^μ(n/d)`, the products and exact quotients by `x^d
//!   − 1` in ℤ[x], each in one pass over the coefficients; the coefficients
//!   from the highest, giac's `[1,0,-1,0,1]` for `n = 12`. `Cyclotomic[n,
//!   x]` is `Φ_n(x)`.
//! - `MultiplicativeOrder[a, n]` (`znorder`, CAS-423): the least `k ≥ 1`
//!   with `a^k ≡ 1 (mod n)`, a divisor of Carmichael's `λ(n)`: each prime
//!   `q` of `λ(n)` divided out while the power stays 1. 0 when `gcd(a, n) >
//!   1`, as giac's; `n ≥ 1`.
//! - `PrimitiveRoot[n]` (`znprimroot`): for `n` 2, 4, `p^k` or `2 p^k`, `p`
//!   an odd prime, the least `g` coprime to `n` with `g^(φ(n)/q) ≢ 1` for
//!   each prime `q` of `φ(n)`; 0 for 1, as giac's; `Indeterminate` for
//!   another `n`, giac's `undef`. giac answers `g % n`, a residue.
//!
//! A number cas-num does not factor fully leaves the call.

use cas_num::{Integer, Natural};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The largest index of a cyclotomic polynomial.
const MAX_CYCLOTOMIC: i64 = 1 << 14;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `QuotientRemainder[a, b]`.
pub(in crate::native) fn quotient_remainder(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let (Some(a), Some(b)) = (a.to_integer(), b.to_integer()) else {
        return Ok(None);
    };
    if b.is_zero() {
        ev.warn(Warning::Native {
            message: "division by 0".into(),
            expr: e.clone(),
        });
        return Ok(None);
    }
    let (q, r) = a.div_rem_euclid(&b);
    Ok(Some(list(vec![Expr::from(q), Expr::from(r)])))
}

/// `SymmetricMod[a, n]`.
pub(in crate::native) fn symmetric_mod(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, n] = e.args() else {
        return Ok(None);
    };
    let (Some(a), Some(n)) = (a.to_integer(), n.to_integer()) else {
        return Ok(None);
    };
    if n.is_zero() {
        return Ok(None);
    }
    let n = Integer::from(n.unsigned_abs());
    let (_, r) = a.div_rem_euclid(&n);
    // r in [0, n): past n/2 it is r − n.
    let twice = &r + &r;
    Ok(Some(Expr::from(if twice > n { &r - &n } else { r })))
}

/// `PseudoprimeQ[n]`.
pub(in crate::native) fn pseudoprime_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.to_integer() else {
        return Ok(None);
    };
    let level = match n.unsigned_abs().primality() {
        cas_num::Primality::Prime => 2,
        cas_num::Primality::ProbablePrime => 1,
        cas_num::Primality::NotPrime => 0,
    };
    Ok(Some(Expr::from(level)))
}

/// Möbius's `μ(m)` from `m`'s prime factors.
fn mobius(m: u64) -> i64 {
    let mut m = m;
    let mut sign = 1;
    let mut p = 2;
    while p * p <= m {
        if m % p == 0 {
            m /= p;
            if m % p == 0 {
                return 0;
            }
            sign = -sign;
        }
        p += 1;
    }
    if m > 1 { -sign } else { sign }
}

/// `p (x^d − 1)`, coefficients from the constant.
fn times(p: &[Integer], d: usize) -> Vec<Integer> {
    let mut q = vec![Integer::ZERO; p.len() + d];
    for (k, c) in p.iter().enumerate() {
        q[k + d] = &q[k + d] + c;
        q[k] = &q[k] - c;
    }
    q
}

/// `p/(x^d − 1)`, exact: `p[k] = q[k − d] − q[k]`.
fn divided(p: &[Integer], d: usize) -> Vec<Integer> {
    let n = p.len() - d;
    let mut q: Vec<Integer> = Vec::with_capacity(n);
    for k in 0..n {
        let before = if k >= d {
            q[k - d].clone()
        } else {
            Integer::ZERO
        };
        q.push(&before - &p[k]);
    }
    q
}

/// `Φ_n`'s coefficients, from the constant.
pub(super) fn cyclotomic_coefficients(n: u64) -> Vec<Integer> {
    let divisors: Vec<u64> = (1..=n).filter(|d| n % d == 0).collect();
    let mut p = vec![Integer::ONE];
    for &d in &divisors {
        if mobius(n / d) == 1 {
            p = times(&p, usize::try_from(d).expect("small"));
        }
    }
    for &d in &divisors {
        if mobius(n / d) == -1 {
            p = divided(&p, usize::try_from(d).expect("small"));
        }
    }
    p
}

/// `Cyclotomic[n]` and `Cyclotomic[n, x]`.
pub(in crate::native) fn cyclotomic(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (n, x) = match e.args() {
        [n] => (n, None),
        [n, x] => (n, Some(x)),
        _ => return Ok(None),
    };
    let Some(n) = n.as_i64().filter(|n| (1..=MAX_CYCLOTOMIC).contains(n)) else {
        return Ok(None);
    };
    let p = cyclotomic_coefficients(n.unsigned_abs());
    Ok(Some(match x {
        None => list(p.into_iter().rev().map(Expr::from).collect()),
        Some(x) => Expr::call(
            Symbol::PLUS,
            p.into_iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(|(k, c)| Expr::from(c) * x.pow(i64::try_from(k).expect("small"))),
        ),
    }))
}

/// `n`'s factorization, primes and exponents; `None` when not complete.
fn factorization(n: &Natural) -> Option<Vec<(Integer, u64)>> {
    let f = n.factor();
    f.is_complete().then(|| {
        f.factors()
            .iter()
            .map(|p| (Integer::from(p.base().clone()), p.exponent()))
            .collect()
    })
}

/// Carmichael's `λ(n)` from `n`'s factorization.
fn carmichael(factors: &[(Integer, u64)]) -> Integer {
    let two = Integer::from(2);
    factors.iter().fold(Integer::ONE, |acc, (p, k)| {
        let k = *k;
        let part = if *p == two {
            match k {
                1 => Integer::ONE,
                2 => two.clone(),
                _ => two.pow(k - 2),
            }
        } else {
            &(p - &Integer::ONE) * &p.pow(k - 1)
        };
        acc.lcm(&part)
    })
}

/// `a^e mod n`, `n ≥ 1`, `e ≥ 0`.
fn power(a: &Integer, e: &Integer, n: &Integer) -> Integer {
    a.pow_mod(e, n).expect("a nonnegative exponent")
}

/// `MultiplicativeOrder[a, n]`.
pub(in crate::native) fn multiplicative_order(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, n] = e.args() else {
        return Ok(None);
    };
    let (Some(a), Some(n)) = (a.to_integer(), n.to_integer()) else {
        return Ok(None);
    };
    if n < Integer::ONE {
        return Ok(None);
    }
    if a.gcd(&n) != Integer::ONE {
        return Ok(Some(Expr::from(0)));
    }
    let Some(factors) = factorization(&n.unsigned_abs()) else {
        return Ok(None);
    };
    let lambda = carmichael(&factors);
    let Some(primes) = factorization(&lambda.unsigned_abs()) else {
        return Ok(None);
    };
    let mut order = lambda;
    for (q, _) in primes {
        loop {
            let (part, rest) = order.div_rem_euclid(&q);
            if !rest.is_zero() || power(&a, &part, &n) != Integer::ONE % &n {
                break;
            }
            order = part;
        }
    }
    Ok(Some(Expr::from(order)))
}

/// `PrimitiveRoot[n]`.
pub(in crate::native) fn primitive_root(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.to_integer() else {
        return Ok(None);
    };
    if n < Integer::ONE {
        return Ok(None);
    }
    let small = n.to_i64();
    match small {
        Some(1) => return Ok(Some(Expr::from(0))),
        Some(2) => return Ok(Some(Expr::from(1))),
        Some(4) => return Ok(Some(Expr::from(3))),
        _ => {}
    }
    let Some(factors) = factorization(&n.unsigned_abs()) else {
        return Ok(None);
    };
    // p^k or 2 p^k, p odd.
    let two = Integer::from(2);
    let odd: Vec<&(Integer, u64)> = factors.iter().filter(|(p, _)| *p != two).collect();
    let twos = factors
        .iter()
        .find(|(p, _)| *p == two)
        .map_or(0, |(_, k)| *k);
    let [(p, k)] = odd.as_slice() else {
        return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
    };
    if twos > 1 {
        return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
    }
    let phi = &(p - &Integer::ONE) * &p.pow(*k - 1);
    let Some(primes) = factorization(&phi.unsigned_abs()) else {
        return Ok(None);
    };
    let one = Integer::ONE;
    let mut g = two;
    while g < n {
        let primitive = g.gcd(&n) == one
            && primes
                .iter()
                .all(|(q, _)| power(&g, &phi.div_exact(q), &n) != one);
        if primitive {
            return Ok(Some(Expr::from(g)));
        }
        g = &g + &one;
    }
    Ok(None)
}

#[cfg(test)]
mod tests;
