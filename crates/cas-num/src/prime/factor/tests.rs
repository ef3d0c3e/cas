//! Factorisation of products of known primes (found with GMP), compared
//! exactly, including the labels.

use std::collections::BTreeMap;

use proptest::prelude::*;
use rug::Integer;
use rug::ops::Pow;

use super::super::small::{groups, primes};
use super::{hits, rho};
use crate::testkit::{Mag, ToRug, strategy};
use crate::{Natural, Primality};

/// A prime of the given size class, from GMP's `next_prime`.
fn prime(bits: std::ops::Range<u32>) -> BoxedStrategy<Integer> {
    (bits, any::<u64>(), strategy::natural_up_to(3))
        .prop_map(|(b, x, big)| {
            let start = if b <= 64 {
                Integer::from(x) >> (64 - b)
            } else {
                big.to_rug() | (Integer::from(1) << (b - 1))
            };
            start.next_prime()
        })
        .boxed()
}

/// Up to six prime powers: tiny primes (trial division), primes between 2^16
/// and 2^34 (rho), and at most one prime of up to 3 limbs (BPSW).
fn product() -> BoxedStrategy<Vec<(Integer, u32)>> {
    let small = (prime(2..17), 1..6u32);
    let medium = (prime(17..35), 1..3u32);
    let factors = prop::collection::vec(prop_oneof![3 => small, 2 => medium], 0..5);
    let large = prop::option::weighted(0.4, (prime(65..190), 1..3u32));
    (factors, large)
        .prop_map(|(mut f, large)| {
            f.extend(large);
            f
        })
        .boxed()
}

fn expected(f: &[(Integer, u32)]) -> Vec<(Integer, u64, Primality)> {
    let mut merged: BTreeMap<Integer, u64> = BTreeMap::new();
    for (p, e) in f {
        *merged.entry(p.clone()).or_default() += u64::from(*e);
    }
    merged
        .into_iter()
        .map(|(p, e)| {
            let label = if p.significant_bits() <= 64 {
                Primality::Prime
            } else {
                Primality::ProbablePrime
            };
            (p, e, label)
        })
        .collect()
}

fn ours(n: &Natural, effort: Option<u64>) -> Vec<(Integer, u64, Primality)> {
    let f = match effort {
        None => n.factor(),
        Some(e) => n.factor_with_effort(e),
    };
    assert_eq!(f.value(), *n, "product of the factors");
    f.factors()
        .iter()
        .map(|f| (f.base().to_rug(), f.exponent(), f.primality()))
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn factors_products_of_known_primes(f in product()) {
        let n = f.iter().fold(Integer::from(1), |n, (p, e)| n * p.clone().pow(*e));
        let n = Mag::from_rug(&n).to_natural();
        prop_assert_eq!(ours(&n, None), expected(&f));
    }

    /// Any value: the factors multiply back, increase, and carry the right
    /// labels.
    #[test]
    fn factors_any_value(n in strategy::natural_up_to(2).prop_filter("non-zero", |n| !n.is_zero())) {
        let n = n.to_natural();
        let f = n.factor_with_effort(1 << 14);
        prop_assert_eq!(f.value(), n);
        for w in f.factors().windows(2) {
            prop_assert!(w[0].base() < w[1].base());
        }
        for x in f.factors() {
            let gmp = x.base().to_rug().is_probably_prime(40) != rug::integer::IsPrime::No;
            prop_assert_eq!(x.primality().is_probable_prime(), gmp);
            prop_assert!(x.exponent() >= 1);
        }
    }
}

#[test]
fn edge_cases() {
    let one = Natural::one().factor();
    assert!(one.factors().is_empty() && one.is_proven() && one.value().is_one());
    let f = Natural::from(2u64).pow(100).factor();
    assert_eq!(ours(&Natural::from(2u64).pow(100), None).len(), 1);
    assert_eq!(f.factors()[0].exponent(), 100);
    // 2^64 - 59 is prime; 2^64 + 1 = 274177 · 67280421310721.
    assert_eq!(
        ours(&Natural::from(u64::MAX - 58), None),
        vec![(Integer::from(u64::MAX - 58), 1, Primality::Prime)]
    );
    assert_eq!(
        ours(&(Natural::power_of_two(64) + 1u64), None),
        vec![
            (Integer::from(274_177), 1, Primality::Prime),
            (Integer::from(67_280_421_310_721u64), 1, Primality::Prime),
        ]
    );
    // A perfect power of a probable prime above 2^64.
    let p = Natural::power_of_two(127) - 1u64;
    assert_eq!(
        ours(&p.pow(3), None),
        vec![(p.to_rug(), 3, Primality::ProbablePrime)]
    );
    // F_6 = 2^64 + 1 squared times F_5's factors.
    let f5 = Natural::power_of_two(32) + 1u64;
    let got = ours(&(&f5 * &(Natural::power_of_two(64) + 1u64)), None);
    assert_eq!(got.len(), 4);
    assert!(Natural::from(641u64).factor().is_proven());
}

#[test]
fn hard_composites_stay_unsplit() {
    // Two primes near 2^50: far beyond 1000 rho steps.
    let p = Integer::from(1u64 << 50).next_prime();
    let q = (p.clone() + 1000u32).next_prime();
    let n = Mag::from_rug(&(p * q)).to_natural();
    let f = n.factor_with_effort(1000);
    assert!(!f.is_complete());
    assert_eq!(
        ours(&n, Some(1000)),
        vec![(n.to_rug(), 1, Primality::NotPrime)]
    );
    // With an unsplit square of it, the exponent is kept.
    let f = n.square().factor_with_effort(1000);
    assert_eq!(f.factors().len(), 1);
    assert_eq!(f.factors()[0].exponent(), 2);
}

#[test]
#[should_panic(expected = "factorisation of zero")]
fn zero_panics() {
    let _ = Natural::ZERO.factor();
}

/// Squares and products of the primes at each end of every trial-division
/// group: the early exit after a group must not call such a composite
/// prime.
#[test]
fn group_boundaries() {
    let primes = primes();
    for (_, range) in groups() {
        for i in [range.start, range.end - 1] {
            let p = Integer::from(primes[i]);
            let q = p.clone().next_prime();
            let label = |x: &Integer| (x.clone(), 1, Primality::Prime);
            let square = Mag::from_rug(&p.clone().square()).to_natural();
            assert_eq!(
                ours(&square, None),
                vec![(p.clone(), 2, Primality::Prime)],
                "{p}^2"
            );
            let product = Mag::from_rug(&(p.clone() * &q)).to_natural();
            assert_eq!(
                ours(&product, None),
                vec![label(&p), label(&q)],
                "{p} * {q}"
            );
        }
    }
}

/// Small semiprimes make a batch of rho differences overshoot to a multiple
/// of n; the replay must still find a proper factor. Both the word and the
/// multi-limb ring are used.
#[test]
fn rho_replays_overshooting_batches() {
    hits::BACKTRACK.with(|c| c.set(0));
    let mut p = Integer::from(1 << 16);
    for i in 0..300 {
        p = p.next_prime();
        let q = if i % 2 == 0 {
            Integer::from(&p + 2000).next_prime()
        } else {
            (Integer::from(1) << 70u32) + Integer::from(&p * 7919)
        }
        .next_prime();
        let n = Mag::from_rug(&(p.clone() * &q)).to_natural();
        let d = rho(&n, 1 << 24).expect("small factor").to_rug();
        assert!(d == p || d == q, "{d} is not a prime factor of {p} * {q}");
    }
    assert!(
        hits::BACKTRACK.with(std::cell::Cell::get) >= 10,
        "replay not reached"
    );
}

/// The groups cover every prime of the table once, in order, and the
/// largest primes below 2^16 are removed by trial division.
#[test]
fn groups_cover_the_table() {
    let primes = primes();
    let mut next = 0;
    for (product, range) in groups() {
        assert_eq!(range.start, next);
        let p = primes[range.clone()]
            .iter()
            .fold(1u64, |p, &q| p.checked_mul(u64::from(q)).unwrap());
        assert_eq!(*product, p);
        next = range.end;
    }
    assert_eq!(next, primes.len());
    assert_eq!(primes.last(), Some(&65521));
    for p in [65521u32, 65519, 65497] {
        let p = Integer::from(p);
        let n = Mag::from_rug(&p.clone().square()).to_natural();
        assert_eq!(ours(&n, None), vec![(p, 2, Primality::Prime)]);
    }
}

/// The effort bound holds: a factor near 2^30 needs about 2^15 rho steps.
#[test]
fn effort_is_bounded() {
    let p = Integer::from(1u64 << 30).next_prime();
    let q = Integer::from(1u64 << 33).next_prime();
    let n = Mag::from_rug(&(p.clone() * &q)).to_natural();
    let needed = (0..=24)
        .map(|k| 1u64 << k)
        .find(|&e| rho(&n, e).is_some())
        .expect("found within 2^24 steps");
    assert!(needed >= 1 << 12, "found after {needed} steps");
}
