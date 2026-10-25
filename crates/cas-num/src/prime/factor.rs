//! Integer factorisation: trial division, perfect powers and Pollard's rho.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, c, x, y, q, g, k, r)"
)]

use super::Primality;
use super::small::{SIEVE_LIMIT, groups, primes};
use crate::kernel::gcd_1;
use crate::modular::{ModRing, ModWord, Residue};
use crate::{Limb, Natural};

#[cfg(test)]
mod tests;

/// The rho iterations [`Natural::factor`] spends on each composite before
/// giving up on it; enough for factors up to about 2^40.
pub(crate) const DEFAULT_RHO_ITERATIONS: u64 = 1 << 22;

/// Products of this many rho differences share one GCD (Brent).
const RHO_BATCH: u64 = 128;

/// A factorisation `n = ∏ base_i^exponent_i`, sorted by increasing base.
///
/// Every base is labelled: [`Primality::Prime`] (proven),
/// [`Primality::ProbablePrime`] (passed Baillie–PSW), or
/// [`Primality::NotPrime`], a composite the search gave up on. Bases are
/// distinct, but an unsplit composite may share prime factors with another.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Factorization {
    factors: Vec<Factor>,
}

/// A base of a [`Factorization`] with its exponent and label.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Factor {
    base: Natural,
    exponent: u64,
    primality: Primality,
}

impl Factor {
    #[must_use]
    pub fn base(&self) -> &Natural {
        &self.base
    }

    #[must_use]
    pub fn exponent(&self) -> u64 {
        self.exponent
    }

    #[must_use]
    pub fn primality(&self) -> Primality {
        self.primality
    }
}

impl Factorization {
    /// The factors, by increasing base; empty for 1.
    #[must_use]
    pub fn factors(&self) -> &[Factor] {
        &self.factors
    }

    /// Whether every base is prime or a probable prime.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.factors.iter().all(|f| f.primality.is_probable_prime())
    }

    /// Whether every base is proven prime.
    #[must_use]
    pub fn is_proven(&self) -> bool {
        self.factors.iter().all(|f| f.primality == Primality::Prime)
    }

    /// The product of the factors.
    #[must_use]
    pub fn value(&self) -> Natural {
        self.factors
            .iter()
            .fold(Natural::one(), |p, f| p * &f.base.pow(f.exponent))
    }
}

impl Natural {
    /// The factorisation of `self` with the default effort: factors up to
    /// about 2^40 are found reliably, larger ones by chance.
    ///
    /// # Panics
    ///
    /// If `self` is zero.
    #[must_use]
    pub fn factor(&self) -> Factorization {
        self.factor_with_effort(DEFAULT_RHO_ITERATIONS)
    }

    /// The factorisation of `self`, spending about `rho_iterations` steps
    /// of Pollard's rho (at most twice that) on each composite.
    ///
    /// Primes below 2^16 are removed by trial division; a perfect power is
    /// factored through its base; the remaining composites are split by
    /// Pollard's rho with Brent's cycle detection, which finds a prime
    /// factor `p` in about `√p` steps. Composites that resist the given
    /// effort are kept, labelled [`Primality::NotPrime`].
    ///
    /// # Panics
    ///
    /// If `self` is zero.
    #[must_use]
    pub fn factor_with_effort(&self, rho_iterations: u64) -> Factorization {
        assert!(!self.is_zero(), "factorisation of zero");
        let mut found = Vec::new();
        let rest = trial_division(self, &mut found);
        if !rest.is_one() {
            split(rest, 1, rho_iterations, &mut found);
        }
        found.sort_by(|a, b| a.base.cmp(&b.base));
        let mut factors: Vec<Factor> = Vec::with_capacity(found.len());
        for f in found {
            match factors.last_mut() {
                Some(last) if last.base == f.base => last.exponent += f.exponent,
                _ => factors.push(f),
            }
        }
        Factorization { factors }
    }
}

/// Removes the primes below 2^16 from `n` into `found`; returns the rest,
/// which is 1, a prime, or has only prime factors above 2^16.
fn trial_division(n: &Natural, found: &mut Vec<Factor>) -> Natural {
    let primes = primes();
    let mut n = n.clone();
    for (product, range) in groups() {
        let r = n.rem_limb(*product);
        for &p in &primes[range.clone()] {
            let p = Limb::from(p);
            if r % p != 0 {
                continue;
            }
            let mut exponent = 0;
            while n.is_divisible_by_limb(p) {
                n = n.div_exact_limb(p);
                exponent += 1;
            }
            found.push(Factor {
                base: Natural::from(p),
                exponent,
                primality: Primality::Prime,
            });
        }
        // Once n < p^2 for the next prime p, it is 1 or prime.
        let next = primes
            .get(range.end)
            .map_or(Limb::from(SIEVE_LIMIT), |&p| Limb::from(p));
        if n < Natural::from(next) * next {
            break;
        }
    }
    let bound = Natural::from(SIEVE_LIMIT).square();
    if !n.is_one() && n < bound {
        found.push(Factor {
            base: n,
            exponent: 1,
            primality: Primality::Prime,
        });
        return Natural::one();
    }
    n
}

/// Factors `n > 1` without prime factors below 2^16, each found base
/// raised to `multiplicity`.
fn split(n: Natural, multiplicity: u64, effort: u64, found: &mut Vec<Factor>) {
    let primality = n.primality();
    if primality.is_probable_prime() {
        found.push(Factor {
            base: n,
            exponent: multiplicity,
            primality,
        });
        return;
    }
    if let Some((base, k)) = n.perfect_power() {
        split(base, multiplicity * k, effort, found);
        return;
    }
    match rho(&n, effort) {
        Some(d) => {
            let e = n.div_exact(&d);
            split(d, multiplicity, effort, found);
            split(e, multiplicity, effort, found);
        }
        None => found.push(Factor {
            base: n,
            exponent: multiplicity,
            primality: Primality::NotPrime,
        }),
    }
}

/// The arithmetic Pollard's rho needs modulo the composite `n`.
trait RhoRing {
    type Elem: Clone + PartialEq;

    fn elem(&self, x: Limb) -> Self::Elem;
    /// `x^2 + c`.
    fn step(&self, x: &Self::Elem, c: &Self::Elem) -> Self::Elem;
    /// `q (x - y)`.
    fn accumulate(&self, q: &Self::Elem, x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    /// `gcd(q, n)`.
    fn gcd(&self, q: &Self::Elem) -> Natural;
}

impl RhoRing for ModWord {
    type Elem = Limb;

    fn elem(&self, x: Limb) -> Limb {
        self.reduce(x)
    }

    fn step(&self, x: &Limb, c: &Limb) -> Limb {
        self.add(self.mul(*x, *x), *c)
    }

    fn accumulate(&self, q: &Limb, x: &Limb, y: &Limb) -> Limb {
        self.mul(*q, self.sub(*x, *y))
    }

    fn gcd(&self, q: &Limb) -> Natural {
        Natural::from(gcd_1(*q, self.modulus()))
    }
}

impl RhoRing for ModRing {
    type Elem = Residue;

    fn elem(&self, x: Limb) -> Residue {
        ModRing::elem(self, &Natural::from(x))
    }

    fn step(&self, x: &Residue, c: &Residue) -> Residue {
        self.add(&self.square(x), c)
    }

    fn accumulate(&self, q: &Residue, x: &Residue, y: &Residue) -> Residue {
        self.mul(q, &self.sub(x, y))
    }

    fn gcd(&self, q: &Residue) -> Natural {
        // n is odd, so the Montgomery factor R does not change the GCD.
        self.to_natural(q).gcd(self.modulus())
    }
}

/// A non-trivial factor of the odd composite `n` (not a perfect power), or
/// `None` after `effort` iterations.
fn rho(n: &Natural, effort: u64) -> Option<Natural> {
    match *n.limbs() {
        [w] => rho_in(&ModWord::new(w), n, effort),
        _ => rho_in(&ModRing::new(n.clone()), n, effort),
    }
}

/// Brent's variant of Pollard's rho (Brent 1980) on `x -> x^2 + c`,
/// retrying with the next `c` when a cycle closes without a factor.
fn rho_in<R: RhoRing>(ring: &R, n: &Natural, effort: u64) -> Option<Natural> {
    let mut spent = 0u64;
    for c in 1.. {
        let c = ring.elem(c);
        let mut y = ring.elem(2);
        let mut q = ring.elem(1);
        let mut g = Natural::one();
        let (mut x, mut ys) = (y.clone(), y.clone());
        let mut r = 1u64;
        while g.is_one() {
            x = y.clone();
            for _ in 0..r {
                y = ring.step(&y, &c);
            }
            let mut k = 0;
            while k < r && g.is_one() {
                ys = y.clone();
                let batch = RHO_BATCH.min(r - k);
                for _ in 0..batch {
                    y = ring.step(&y, &c);
                    q = ring.accumulate(&q, &x, &y);
                }
                g = ring.gcd(&q);
                k += batch;
            }
            spent += 2 * r;
            if spent > effort {
                return None;
            }
            r *= 2;
        }
        if g == *n {
            // The batch overshot: replay it one step at a time.
            #[cfg(test)]
            hits::BACKTRACK.with(|c| c.set(c.get() + 1));
            loop {
                ys = ring.step(&ys, &c);
                g = ring.gcd(&ring.accumulate(&ring.elem(1), &x, &ys));
                if !g.is_one() {
                    break;
                }
            }
        }
        if g != *n {
            return Some(g);
        }
    }
    unreachable!("the loop over c only ends by returning")
}

/// Counts how often a batch overshoots and is replayed step by step, so
/// tests can check that their inputs reach the replay.
///
/// Last in the file: clippy's `items_after_test_module` takes an inline
/// `#[cfg(test)]` module for the tests (CAS-174).
#[cfg(test)]
pub(crate) mod hits {
    use std::cell::Cell;

    thread_local! {
        pub(crate) static BACKTRACK: Cell<u64> = const { Cell::new(0) };
    }
}
