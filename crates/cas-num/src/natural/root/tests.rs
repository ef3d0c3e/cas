//! `Natural` roots and perfect powers checked against GMP, with the root
//! bounds checked directly.

use proptest::prelude::*;
use rug::Integer;
use rug::ops::Pow;

use super::is_prime_u64;
use crate::Natural;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

/// `s^n + d` for `d` in `-1..=1`, or one of the extremes `s^n + 2s` and
/// `(s + 1)^n` of a square root's remainder: the values where a root or its
/// correction step is off by one.
fn near_power(max_limbs: usize, n: u32) -> BoxedStrategy<Mag> {
    (strategy::natural_up_to(max_limbs), 0..5u8)
        .prop_map(move |(s, shape)| {
            let s = s.to_rug();
            let p = Integer::from((&s).pow(n));
            let x = match shape {
                0 => p - 1u32,
                1 => p,
                2 => p + 1u32,
                3 => p + Integer::from(&s * 2u32),
                _ => (s + 1u32).pow(n),
            };
            Mag::from_rug(&x.max(Integer::new()))
        })
        .boxed()
}

/// `s² + 64 j P`, with `P` the product of the residue primes: congruent to a
/// square modulo everything `is_perfect_square` filters by, but rarely a
/// square, so only the root itself can reject it.
fn filter_passing() -> BoxedStrategy<Mag> {
    (strategy::natural_up_to(DEFAULT_MAX_LIMBS / 2), 1..=1000u32)
        .prop_map(|(s, j)| {
            let p = Integer::from(super::RESIDUE_PRODUCT) * 64u32 * j;
            Mag::from_rug(&(s.to_rug().square() + p))
        })
        .boxed()
}

fn square_input() -> BoxedStrategy<Mag> {
    prop_oneof![
        3 => strategy::natural(),
        3 => near_power(DEFAULT_MAX_LIMBS / 2, 2),
        1 => filter_passing(),
    ]
    .boxed()
}

/// A root degree: small, around the 64-bit word, or large.
fn degree() -> BoxedStrategy<u32> {
    prop_oneof![
        6 => 1..=12u32,
        2 => 13..=200u32,
        1 => prop::sample::select(vec![1000, 4095, 65_536, u32::MAX]),
    ]
    .boxed()
}

/// A degree with a value near its powers; bases shrink as the degree grows.
fn root_input() -> BoxedStrategy<(Mag, u32)> {
    let near = (2..=40u32).prop_flat_map(|n| {
        let limbs = (60 / n as usize).max(1);
        (near_power(limbs, n), Just(n))
    });
    prop_oneof![
        2 => (strategy::natural(), degree()),
        2 => near,
    ]
    .boxed()
}

/// `a^k` with a small base and exponent, near misses, and powers with a
/// factor of two split off, which exercises the trailing-zero rule.
fn power_input() -> BoxedStrategy<Mag> {
    let power =
        (strategy::natural_up_to(4), 2..=24u32, -1i32..=1, 0..3u32).prop_map(|(a, k, d, twos)| {
            let a = a.to_rug() << twos;
            let x = a.pow(k) + d;
            Mag::from_rug(&x.max(Integer::new()))
        });
    let twos = (0..2000u32, -1i32..=1)
        .prop_map(|(t, d)| Mag::from_rug(&((Integer::from(1) << t) + d).max(Integer::new())));
    let mixed = (2..1000u32, 1..=6u32, 0..=6u32)
        .prop_map(|(m, k, t)| Mag::from_rug(&(Integer::from(m).pow(k) << (k * t))));
    prop_oneof![
        4 => power,
        2 => strategy::natural(),
        1 => twos,
        1 => mixed,
    ]
    .boxed()
}

differential! {
    fn sqrt_rem(a in square_input()) {
        ours: {
            let x = a.to_natural();
            let (s, r) = x.sqrt_rem();
            prop_assert!(s.square() <= x && x < (&s + 1).square(), "root bounds");
            prop_assert_eq!(&x.isqrt(), &s);
            (s, r)
        },
        gmp: a.to_rug().sqrt_rem(Integer::new()),
    }

    fn root_rem((a, n) in root_input()) {
        ours: {
            let x = a.to_natural();
            let (s, r) = x.root_rem(n);
            prop_assert_eq!(&x.nth_root(n), &s);
            (s, r)
        },
        gmp: a.to_rug().root_rem(Integer::new(), n),
    }

    fn is_perfect_square(a in square_input()) {
        ours: a.to_natural().is_perfect_square(),
        gmp: a.to_rug().is_perfect_square(),
    }

    fn is_perfect_power(a in power_input()) {
        ours: a.to_natural().is_perfect_power(),
        gmp: a.to_rug().is_perfect_power(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// The root bounds on large values (several recursion levels).
    #[test]
    fn sqrt_rem_large(a in strategy::natural_up_to(400)) {
        let x = a.to_natural();
        let (s, r) = x.sqrt_rem();
        prop_assert!(s.square() + &r == x && r <= &s + &s, "root bounds");
    }

    #[test]
    fn root_rem_large(a in strategy::natural_up_to(400), n in 3..=9u32) {
        let x = a.to_natural();
        let (s, r) = x.root_rem(n);
        prop_assert!(r < (&s + 1).pow(u64::from(n)) - s.pow(u64::from(n)), "root bounds");
    }
}

proptest! {
    /// The largest exponent: `a^k` is the value, `a` is not a perfect
    /// power, and there is one exactly when GMP finds a perfect power.
    #[test]
    fn perfect_power(a in power_input()) {
        let x = a.to_natural();
        let found = x.perfect_power();
        let gmp = a.to_rug() > 1 && a.to_rug().is_perfect_power();
        prop_assert_eq!(found.is_some(), gmp);
        if let Some((base, k)) = found {
            prop_assert!(k >= 2 && base.pow(k) == x, "not a power");
            prop_assert!(!base.is_perfect_power(), "exponent not largest");
        }
    }
}

#[test]
fn perfect_power_examples() {
    let pp = |x: u64| {
        Natural::from(x)
            .perfect_power()
            .map(|(a, k)| (u64::try_from(a).unwrap(), k))
    };
    assert_eq!(pp(0), None);
    assert_eq!(pp(1), None);
    assert_eq!(pp(2), None);
    assert_eq!(pp(4), Some((2, 2)));
    assert_eq!(pp(64), Some((2, 6)));
    assert_eq!(pp(72), None);
    assert_eq!(pp(729), Some((3, 6)));
    assert_eq!(pp(1 << 63), Some((2, 63)));
    assert_eq!(pp(3u64.pow(40)), Some((3, 40)));
    assert_eq!(pp(1_000_000), Some((10, 6)));
    assert!(Natural::ZERO.is_perfect_power() && Natural::one().is_perfect_power());
    let big = Natural::from(10u64).pow(1000);
    assert_eq!(big.perfect_power(), Some((Natural::from(10u64), 1000)));
    let prime_power = Natural::from(1_000_003u64).pow(97);
    assert_eq!(
        prime_power.perfect_power(),
        Some((Natural::from(1_000_003u64), 97))
    );
    assert_eq!((prime_power + 1).perfect_power(), None);
}

#[test]
fn root_edge_cases() {
    let two_64 = Natural::power_of_two(64);
    assert_eq!(Natural::ZERO.nth_root(7), Natural::ZERO);
    assert_eq!(two_64.nth_root(64), Natural::from(2u64));
    assert_eq!((&two_64 - 1).nth_root(64), Natural::one());
    assert_eq!(two_64.nth_root(65), Natural::one());
    assert_eq!(two_64.nth_root(u32::MAX), Natural::one());
    assert_eq!(Natural::from(u64::MAX).nth_root(u32::MAX), Natural::one());
    assert_eq!(
        Natural::from(u64::MAX).isqrt(),
        Natural::from(u64::from(u32::MAX))
    );
    let (s, r) = Natural::power_of_two(1001).sqrt_rem();
    assert_eq!(s.square() + r, Natural::power_of_two(1001));
}

#[test]
#[should_panic(expected = "zeroth root")]
fn zeroth_root_panics() {
    let _ = Natural::one().nth_root(0);
}

#[test]
fn prime_test_matches_gmp() {
    let check = |n: u64| {
        let gmp = Integer::from(n).is_probably_prime(30) != rug::integer::IsPrime::No;
        assert_eq!(is_prime_u64(n), gmp, "{n}");
    };
    (0..20_000).for_each(check);
    (u64::MAX - 20_000..=u64::MAX).for_each(check);
    // Strong pseudoprimes to several small bases.
    [3_215_031_751, 3_825_123_056_546_413_051]
        .into_iter()
        .for_each(check);
}
