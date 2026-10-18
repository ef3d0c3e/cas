//! `Natural` powers and modular powers checked against GMP.

use proptest::prelude::*;
use rug::Integer;
use rug::ops::Pow;

use crate::Natural;
use crate::testkit::strategy;
use crate::testkit::{Mag, ToRug, differential};

/// A base and exponent with a result of at most a few thousand limbs: large
/// exponents only on small bases.
fn base_exp() -> BoxedStrategy<(Mag, u32)> {
    prop_oneof![
        3 => (strategy::natural(), 0..=12u32),
        2 => (strategy::natural_up_to(2), 0..=300u32),
        1 => (prop::sample::select(vec![0u64, 1, 2, 3, 4, 6, 1 << 40]), 0..=5000u32)
            .prop_map(|(b, e)| (Mag::from_u64(b), e)),
    ]
    .boxed()
}

/// A modulus: one, a power of two, a single limb, an odd one of three to
/// eight limbs (stack arithmetic), or any.
fn modulus() -> BoxedStrategy<Mag> {
    prop_oneof![
        1 => Just(Mag::from_u64(1)),
        1 => (0..300u32).prop_map(|k| Mag::from_rug(&(Integer::from(1) << k))),
        2 => strategy::natural_up_to(1),
        2 => strategy::natural_with_len(3..=8usize).prop_map(|m| Mag::from_rug(&(m.to_rug() | 1))),
        4 => strategy::natural_up_to(12),
    ]
    .prop_filter("non-zero", |m| !m.is_zero())
    .boxed()
}

differential! {
    fn pow((a, e) in base_exp()) {
        ours: a.to_natural().pow(u64::from(e)),
        gmp: a.to_rug().pow(e),
    }

    fn pow_mod(a in strategy::natural_up_to(16), e in strategy::natural_up_to(4), m in modulus()) {
        ours: a.to_natural().pow_mod(&e.to_natural(), &m.to_natural()),
        gmp: a.to_rug().pow_mod(&e.to_rug(), &m.to_rug()).unwrap(),
    }
}

#[test]
fn pow_edge_cases() {
    assert_eq!(Natural::ZERO.pow(0), Natural::one());
    assert_eq!(Natural::ZERO.pow(1), Natural::ZERO);
    assert_eq!(
        Natural::from(2u64).pow(100_000),
        Natural::power_of_two(100_000)
    );
    assert_eq!(Natural::one().pow(u64::MAX), Natural::one());
    let m = Natural::from(7u64);
    assert_eq!(
        Natural::from(3u64).pow_mod(&Natural::ZERO, &m),
        Natural::one()
    );
    assert_eq!(Natural::ZERO.pow_mod(&Natural::ZERO, &m), Natural::one());
    assert_eq!(
        Natural::from(3u64).pow_mod(&Natural::ZERO, &Natural::one()),
        Natural::ZERO
    );
}

#[test]
#[should_panic(expected = "power too large")]
fn huge_power_panics() {
    let _ = Natural::from(3u64).pow(u64::MAX);
}

#[test]
#[should_panic(expected = "power too large")]
fn huge_power_of_two_panics() {
    let _ = Natural::from(4u64).pow(u64::MAX);
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn pow_mod_zero_panics() {
    let _ = Natural::one().pow_mod(&Natural::one(), &Natural::ZERO);
}
