//! The primes below 2^16, and groups of them multiplied into one limb so
//! that trial division costs one multi-precision remainder per group.

use std::ops::Range;
use std::sync::OnceLock;

/// The table holds the primes below this.
pub(crate) const SIEVE_LIMIT: u32 = 1 << 16;

/// [`Natural::primality`](crate::Natural::primality) trial-divides by the
/// primes below this before the probabilistic tests.
pub(crate) const TRIAL_LIMIT: u32 = 1 << 12;

/// The primes below [`SIEVE_LIMIT`], in increasing order (sieve of
/// Eratosthenes).
pub(crate) fn primes() -> &'static [u32] {
    static PRIMES: OnceLock<Vec<u32>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let n = SIEVE_LIMIT as usize;
        let mut composite = vec![false; n];
        let mut primes = Vec::new();
        for i in 2..n {
            if composite[i] {
                continue;
            }
            primes.push(u32::try_from(i).expect("below 2^16"));
            for j in (i * i..n).step_by(i) {
                composite[j] = true;
            }
        }
        primes
    })
}

/// Consecutive runs of [`primes`] with their product, each product below
/// 2^64, as `(product, index range)`.
pub(crate) fn groups() -> &'static [(u64, Range<usize>)] {
    static GROUPS: OnceLock<Vec<(u64, Range<usize>)>> = OnceLock::new();
    GROUPS.get_or_init(|| {
        let primes = primes();
        let mut groups = Vec::new();
        let mut start = 0;
        let mut product = 1u64;
        for (i, &p) in primes.iter().enumerate() {
            if let Some(q) = product.checked_mul(u64::from(p)) {
                product = q;
            } else {
                groups.push((product, start..i));
                start = i;
                product = u64::from(p);
            }
        }
        groups.push((product, start..primes.len()));
        groups
    })
}

/// Whether `n < SIEVE_LIMIT` is prime.
pub(crate) fn is_small_prime(n: u32) -> bool {
    debug_assert!(n < SIEVE_LIMIT);
    primes().binary_search(&n).is_ok()
}
