//! Primality proofs: the `n − 1` test of Brillhart, Lehmer and Selfridge
//! \[BLS75\] over a partial factorisation of `n − 1`.
#![allow(
    clippy::many_single_char_names,
    reason = "the textbook names (n, m, f, r, q, s, a, y, g)"
)]

use super::Primality;
use crate::{Integer, Natural};

#[cfg(test)]
mod tests;

/// The most bases tried for one prime `q` of `n − 1`: for a prime `n`, a
/// base fails with probability `1/q ≤ 1/2`.
const MAX_BASES: u64 = 1000;

/// The outcome of [`Natural::prove_prime`].
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Proof {
    /// Proven prime.
    Proven,
    /// Composite, or below 2.
    Composite,
    /// A probable prime the limits did not prove: `n − 1` was not factored
    /// far enough.
    Unknown,
}

impl Natural {
    /// Whether `self` is prime, with a proof.
    ///
    /// Below 2^64, and for a composite, [`Natural::primality`] decides. A
    /// Baillie–PSW probable prime `n` above is proven by the `n − 1` test
    /// \[BLS75\]: `n − 1 = F R` with `F` the product of the full powers of
    /// the primes `q` of `n − 1` that are found (trial division, then
    /// `rho_iterations` of Pollard's rho and `curves` curves of the
    /// elliptic curve method at stage-1 bound `b1`, as
    /// [`Natural::factor_with_ecm`]) and proven (by this test, with the
    /// same limits). Each `q` needs a base `a` with `a^(n−1) ≡ 1` and
    /// `gcd(a^((n−1)/q) − 1, n) = 1` (Pocklington); then every prime factor
    /// of `n` is `1 mod F`. `n` is prime when `F² > n`, or by \[BLS75\]
    /// Theorem 5 when `R = 2 F s + r`, `1 ≤ r < 2F`,
    /// `(F + 1)(2F² + (r − 1) F + 1) > n` and `s = 0` or `r² − 8s` is not
    /// a square. A base that shows `n` composite gives
    /// [`Proof::Composite`]; an `F` too small, or no base among the first
    /// thousand, [`Proof::Unknown`].
    #[must_use]
    pub fn prove_prime(&self, rho_iterations: u64, b1: u64, curves: u32) -> Proof {
        match self.primality() {
            Primality::NotPrime => Proof::Composite,
            Primality::Prime => Proof::Proven,
            Primality::ProbablePrime => n_minus_one(self, rho_iterations, b1, curves),
        }
    }
}

/// The `n − 1` test of [`Natural::prove_prime`] for an odd `n > 2^64`.
fn n_minus_one(n: &Natural, rho_iterations: u64, b1: u64, curves: u32) -> Proof {
    let m = n - 1;
    let mut f = Natural::one();
    let mut r = m.clone();
    for factor in m.factor_with_ecm(rho_iterations, b1, curves).factors() {
        let q = factor.base();
        let proven = match factor.primality() {
            Primality::Prime => true,
            Primality::ProbablePrime => q.prove_prime(rho_iterations, b1, curves) == Proof::Proven,
            Primality::NotPrime => false,
        };
        if !proven {
            continue;
        }
        // The full power of q: an unsplit composite may hold more of it.
        loop {
            let (quotient, rest) = r.div_rem(q);
            if !rest.is_zero() {
                break;
            }
            r = quotient;
            f = &f * q;
        }
        match pocklington(n, &m, q) {
            Proof::Proven => {}
            other => return other,
        }
    }
    if &f * &f > *n {
        return Proof::Proven;
    }
    // F is even (2 divides n − 1 and is proven) and prime to R, so r is odd.
    let two_f = &f + &f;
    let (s, r) = r.div_rem(&two_f);
    let bound = (&f + 1) * &(&(&(&two_f * &f) + &(&(&r - 1) * &f)) + 1);
    if bound <= *n {
        return Proof::Unknown;
    }
    if s.is_zero() {
        return Proof::Proven;
    }
    let d = Integer::from(&r * &r) - Integer::from(&s * 8);
    if d.is_negative() || d.isqrt().pow(2) != d {
        Proof::Proven
    } else {
        // [BLS75] makes n composite here; we only claim what we checked.
        Proof::Unknown
    }
}

/// Pocklington's condition for the prime `q` of `m = n − 1`: the least
/// base `a ≥ 2` with `a^m ≡ 1` and `gcd(a^(m/q) − 1, n) = 1` gives
/// [`Proof::Proven`]; a base that fails Fermat's test or meets a proper
/// factor gives [`Proof::Composite`]; none among [`MAX_BASES`],
/// [`Proof::Unknown`].
fn pocklington(n: &Natural, m: &Natural, q: &Natural) -> Proof {
    let cofactor = m.div_exact(q);
    for a in 2..2 + MAX_BASES {
        let a = Natural::from(a);
        if !a.pow_mod(m, n).is_one() {
            return Proof::Composite;
        }
        let y = a.pow_mod(&cofactor, n);
        if y.is_zero() {
            return Proof::Composite;
        }
        let g = (y - 1).gcd(n);
        if g.is_one() {
            return Proof::Proven;
        }
        if g != *n {
            return Proof::Composite;
        }
    }
    Proof::Unknown
}
