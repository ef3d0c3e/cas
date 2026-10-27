//! The sieves against known values and the primality test.

use super::{NTH_PRIME_LIMIT, PRIME_PI_LIMIT, nth_prime, prime_pi};
use crate::Natural;

#[test]
fn known_counts() {
    for (x, pi) in [
        (0, 0),
        (1, 0),
        (2, 1),
        (3, 2),
        (4, 2),
        (10, 4),
        (100, 25),
        (1000, 168),
        (10_000, 1229),
        (100_000, 9592),
        (1_000_000, 78_498),
    ] {
        assert_eq!(prime_pi(x), Some(pi), "π({x})");
    }
    assert_eq!(prime_pi(PRIME_PI_LIMIT + 1), None);
}

#[test]
fn known_primes() {
    for (n, p) in [
        (1, 2),
        (2, 3),
        (3, 5),
        (5, 11),
        (6, 13),
        (10, 29),
        (100, 541),
        (1000, 7919),
        (10_000, 104_729),
        (100_000, 1_299_709),
    ] {
        assert_eq!(nth_prime(n), Some(p), "p_{n}");
    }
    assert_eq!(nth_prime(0), None);
    assert_eq!(nth_prime(NTH_PRIME_LIMIT + 1), None);
}

/// Each count against the primality test, across segment boundaries.
#[test]
fn against_primality() {
    let mut count = 0;
    for x in 0..=600_000u64 {
        if Natural::from(x).primality().is_probable_prime() {
            count += 1;
        }
        if x % 9973 == 0 || (524_280..=524_300).contains(&x) {
            assert_eq!(prime_pi(x), Some(count), "π({x})");
        }
    }
}
