//! `ModWord` checked against GMP.

use proptest::prelude::*;
use rug::Integer;

use super::ModWord;
use crate::testkit::strategy;
use crate::testkit::{ToRug, differential};

/// Largest prime below 2^64.
const LARGEST_PRIME: u64 = u64::MAX - 58;

/// A word modulus: 1, 2, small, prime near 2^64, at the normalisation edge
/// (top bit set or not), a power of two, or any.
fn modulus() -> BoxedStrategy<u64> {
    prop_oneof![
        1 => prop::sample::select(vec![
            1,
            2,
            3,
            LARGEST_PRIME,
            (1 << 61) - 1,
            u64::MAX,
            1 << 63,
            (1 << 63) - 1,
            (1 << 63) + 1,
        ]),
        1 => (0..64u32).prop_map(|k| 1u64 << k),
        2 => 1..1000u64,
        4 => strategy::top_limb(),
    ]
    .boxed()
}

/// A modulus and two reduced elements, weighted towards 0, 1 and `n - 1`.
fn operands() -> BoxedStrategy<(u64, u64, u64)> {
    modulus()
        .prop_flat_map(|n| (Just(n), element(n), element(n)))
        .boxed()
}

fn element(n: u64) -> BoxedStrategy<u64> {
    prop_oneof![
        1 => Just(0),
        1 => Just(1 % n),
        1 => Just(n - 1),
        4 => strategy::limb().prop_map(move |x| x % n),
    ]
    .boxed()
}

fn rug_mod(x: Integer, n: u64) -> Integer {
    x.modulo(&Integer::from(n))
}

differential! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    fn add((n, a, b) in operands()) {
        ours: ModWord::new(n).add(a, b),
        gmp: rug_mod(Integer::from(a) + b, n),
    }

    fn sub((n, a, b) in operands()) {
        ours: ModWord::new(n).sub(a, b),
        gmp: rug_mod(Integer::from(a) - b, n),
    }

    fn neg((n, a, _b) in operands()) {
        ours: ModWord::new(n).neg(a),
        gmp: rug_mod(-Integer::from(a), n),
    }

    fn mul((n, a, b) in operands()) {
        ours: ModWord::new(n).mul(a, b),
        gmp: rug_mod(Integer::from(a) * b, n),
    }

    fn reduce(n in modulus(), a in strategy::limb()) {
        ours: ModWord::new(n).reduce(a),
        gmp: Integer::from(a % n),
    }

    fn reduce_u128(n in modulus(), hi in strategy::limb(), lo in strategy::limb()) {
        ours: ModWord::new(n).reduce_u128((u128::from(hi) << 64) | u128::from(lo)),
        gmp: rug_mod((Integer::from(hi) << 64) + lo, n),
    }

    fn reduce_natural(n in modulus(), a in strategy::natural()) {
        ours: ModWord::new(n).reduce_natural(&a.to_natural()),
        gmp: rug_mod(a.to_rug(), n),
    }

    fn reduce_integer(n in modulus(), a in strategy::integer()) {
        ours: ModWord::new(n).reduce_integer(&a.to_integer()),
        gmp: rug_mod(a.to_rug(), n),
    }

    fn pow((n, a, _b) in operands(), e in strategy::limb()) {
        ours: ModWord::new(n).pow(a, e),
        gmp: Integer::from(a).pow_mod(&Integer::from(e), &Integer::from(n)).unwrap(),
    }

    fn pow_natural((n, a, _b) in operands(), e in strategy::natural_up_to(3)) {
        ours: ModWord::new(n).pow_natural(a, &e.to_natural()),
        gmp: Integer::from(a).pow_mod(&e.to_rug(), &Integer::from(n)).unwrap(),
    }

    /// Exponents long enough for every window width (8 from 4609 bits).
    fn pow_natural_long((n, a, _b) in operands(), e in strategy::natural_up_to(130)) {
        ours: ModWord::new(n).pow_natural(a, &e.to_natural()),
        gmp: Integer::from(a).pow_mod(&e.to_rug(), &Integer::from(n)).unwrap(),
    }

    fn inv((n, a, _b) in operands()) {
        ours: ModWord::new(n).inv(a),
        gmp: Integer::from(a).invert(&Integer::from(n)).ok(),
    }
}

#[test]
fn edge_cases() {
    let one = ModWord::new(1);
    assert_eq!(one.reduce(u64::MAX), 0);
    assert_eq!(one.pow(0, 0), 0);
    assert_eq!(one.inv(0), Some(0));
    let w = ModWord::new(LARGEST_PRIME);
    assert_eq!(w.pow(0, 0), 1);
    assert_eq!(w.pow(5, LARGEST_PRIME - 1), 1);
    assert_eq!(w.inv(0), None);
    assert_eq!(w.modulus(), LARGEST_PRIME);
    assert_eq!(ModWord::new(6).inv(4), None);
}

#[test]
#[should_panic(expected = "zero modulus")]
fn zero_modulus_panics() {
    let _ = ModWord::new(0);
}
