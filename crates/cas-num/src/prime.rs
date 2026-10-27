//! Primality testing and integer factorisation.
//!
//! - [`Natural::primality`]: trial division, then Miller–Rabin with a set of
//!   bases that is deterministic below 2^64, and the Baillie–PSW test above.
//!   Results are labelled: [`Primality::Prime`] is proven,
//!   [`Primality::ProbablePrime`] passed BPSW, for which no counterexample is
//!   known.
//! - [`Natural::is_strong_probable_prime`], [`Natural::next_prime`],
//!   [`Integer::jacobi`](crate::Integer::jacobi).
//! - [`Natural::factor`]: trial division, perfect powers and Pollard's rho
//!   (Brent's variant), with an effort bound; factors it cannot split stay
//!   labelled [`Primality::NotPrime`].
//! - [`prime_pi`] and [`nth_prime`]: a segmented sieve of Eratosthenes, to
//!   [`PRIME_PI_LIMIT`] and [`NTH_PRIME_LIMIT`].

mod factor;
mod jacobi;
mod probable;
mod sieve;
mod small;

pub use factor::{Factor, Factorization};
pub(crate) use probable::is_prime_u64;
pub use sieve::{NTH_PRIME_LIMIT, PRIME_PI_LIMIT, nth_prime, prime_pi};

#[cfg(doc)]
use crate::Natural;

/// The outcome of a primality test.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Primality {
    /// Not prime: composite, or below 2.
    NotPrime,
    /// Passed the Baillie–PSW test. No composite is known to pass it, but
    /// none is proven not to exist above 2^64.
    ProbablePrime,
    /// Proven prime.
    Prime,
}

impl Primality {
    /// Whether the value is prime or a probable prime.
    #[must_use]
    pub fn is_probable_prime(self) -> bool {
        self != Primality::NotPrime
    }
}
