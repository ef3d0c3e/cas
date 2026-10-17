//! `Integer` powers, roots and perfect powers checked against GMP.

use proptest::prelude::*;
use rug::Integer as Gmp;
use rug::ops::Pow;

use super::Integer;
use crate::testkit::strategy;
use crate::testkit::{Signed, ToRug, differential};

/// Values near the `i64` bounds, any `i64`, larger values, and small powers
/// with either sign.
fn value() -> BoxedStrategy<Signed> {
    let power = (-40i64..=40, 2..=12u32, -1i64..=1)
        .prop_map(|(a, k, d)| Signed::from_rug(&(Gmp::from(a).pow(k) + d)));
    prop_oneof![
        3 => strategy::integer(),
        2 => power,
    ]
    .boxed()
}

fn exponent() -> BoxedStrategy<u32> {
    prop_oneof![3 => 0..=8u32, 1 => 9..=130u32].boxed()
}

/// A degree; even degrees are kept for non-negative values in the tests.
fn degree() -> BoxedStrategy<u32> {
    prop_oneof![4 => 1..=9u32, 1 => 10..=130u32].boxed()
}

differential! {
    fn pow(a in strategy::integer(), e in exponent()) {
        ours: { let x = a.to_integer(); if x.bit_length() * u64::from(e) > 100_000 { x } else { x.pow(u64::from(e)) } },
        gmp: { let x = a.to_rug(); if u64::from(x.significant_bits()) * u64::from(e) > 100_000 { x } else { x.pow(e) } },
    }

    fn pow_mod(a in strategy::integer(), e in strategy::signed_up_to(3), m in strategy::integer()) {
        ours: a.to_integer().pow_mod(&e.to_integer(), &m.to_integer()),
        gmp: {
            let m = m.to_rug();
            if m == 0 { None } else { a.to_rug().pow_mod(&e.to_rug(), &m).ok() }
        },
    }

    fn sqrt_rem(a in value()) {
        ours: {
            let x = a.to_integer().abs();
            prop_assert_eq!(&x.isqrt(), &x.sqrt_rem().0);
            x.sqrt_rem()
        },
        gmp: a.to_rug().abs().sqrt_rem(Gmp::new()),
    }

    fn root_rem(a in value(), n in degree()) {
        ours: {
            let x = a.to_integer();
            let x = if n % 2 == 0 { x.abs() } else { x };
            prop_assert_eq!(&x.nth_root(n), &x.root_rem(n).0);
            x.root_rem(n)
        },
        gmp: {
            let x = a.to_rug();
            let x = if n % 2 == 0 { x.abs() } else { x };
            x.root_rem(Gmp::new(), n)
        },
    }

    fn is_perfect_square(a in value()) {
        ours: a.to_integer().is_perfect_square(),
        gmp: a.to_rug().is_perfect_square(),
    }

    fn is_perfect_power(a in value()) {
        ours: a.to_integer().is_perfect_power(),
        gmp: a.to_rug().is_perfect_power(),
    }
}

#[test]
fn edge_cases() {
    let i = Integer::from;
    assert_eq!(i(-2).pow(63), i(i64::MIN));
    assert_eq!(i(2).pow(63).to_string(), "9223372036854775808");
    assert_eq!(i(-1).pow(u64::MAX), i(-1));
    assert_eq!(i(-27).nth_root(3), i(-3));
    assert_eq!(i(-26).root_rem(3), (i(-2), i(-18)));
    assert!(i(-1).is_perfect_power() && i(-8).is_perfect_power() && i(0).is_perfect_power());
    assert!(!i(-4).is_perfect_power() && !i(-16).is_perfect_power());
    assert!(i(-64).is_perfect_power() && !i(-4).is_perfect_square());
    assert_eq!(i(3).pow_mod(&i(-1), &i(7)), Some(i(5)));
    assert_eq!(i(2).pow_mod(&i(-1), &i(4)), None);
    assert_eq!(i(2).pow_mod(&i(1), &i(0)), None);
    assert_eq!(i(-2).pow_mod(&i(3), &i(-5)), Some(i(2)));
}

#[test]
#[should_panic(expected = "argument of integer square root cannot be negative")]
fn negative_sqrt_panics() {
    let _ = Integer::from(-1).isqrt();
}

#[test]
#[should_panic(expected = "even root of a negative integer")]
fn negative_even_root_panics() {
    let _ = Integer::from(-1).nth_root(4);
}
