//! Two-limb Montgomery powers checked against GMP.

use proptest::prelude::*;
use rug::Integer;

use super::{ModDouble, mul_wide, pow_mod_even};
use crate::testkit::{ToRug, strategy};

/// Odd moduli in `(2^64, 2^128)`: random, just above `2^64`, just below
/// `2^128` (where the REDC sum overflows), and `2^127 ± 1`.
fn modulus() -> BoxedStrategy<u128> {
    prop_oneof![
        any::<u128>().prop_map(|n| (n | 1) | 1 << 64),
        (0..1000u128).prop_map(|k| (1 << 64) + 2 * k + 1),
        (0..1000u128).prop_map(|k| u128::MAX - 2 * k),
        Just((1 << 127) + 1),
        Just((1 << 127) - 1),
    ]
    .boxed()
}

proptest! {
    #[test]
    fn mul_wide_matches(a in any::<u128>(), b in any::<u128>()) {
        let (lo, hi) = mul_wide(a, b);
        let full = (Integer::from(hi) << 128) + lo;
        prop_assert_eq!(full, Integer::from(a) * b);
    }

    #[test]
    fn pow_matches(
        n in modulus(),
        a in any::<u128>(),
        e in prop_oneof![Just(crate::testkit::Mag::default()), strategy::natural_up_to(130)],
    ) {
        let a = a % n;
        let ours = ModDouble::new(n).pow(a, &e.to_natural());
        let gmp = Integer::from(a).pow_mod(&e.to_rug(), &Integer::from(n)).unwrap();
        prop_assert_eq!(Integer::from(ours), gmp);
    }
}

/// Bases at the edges: 0, 1, `n - 1` (which squares to 1), with exponents
/// 0, 1, 2 and a long one.
#[test]
fn pow_edges() {
    for n in [(1u128 << 64) + 1, u128::MAX, (1 << 127) - 1] {
        let d = ModDouble::new(n);
        for a in [0, 1, 2, n - 1, n - 2] {
            for e in [0u64, 1, 2, 3, u64::MAX] {
                let ours = d.pow(a, &crate::Natural::from(e));
                let gmp = Integer::from(a)
                    .pow_mod(&Integer::from(e), &Integer::from(n))
                    .unwrap();
                assert_eq!(Integer::from(ours), gmp, "n = {n}, a = {a}, e = {e}");
            }
        }
    }
}

proptest! {
    /// The Montgomery-form element operations, mapped back by REDC.
    #[test]
    fn element_ops(n in modulus(), x in any::<u128>(), y in any::<u128>()) {
        let d = ModDouble::new(n);
        let (x, y) = (x % n, y % n);
        let (mx, my) = (d.to_mont(x), d.to_mont(y));
        let back = |m: u128| Integer::from(d.mul(m, 1));
        let gn = Integer::from(n);
        let g = |v: Integer| v.modulo(&gn);
        prop_assert_eq!(back(mx), Integer::from(x));
        prop_assert_eq!(back(d.add(mx, my)), g(Integer::from(x) + y));
        prop_assert_eq!(back(d.sub(mx, my)), g(Integer::from(x) - y));
        prop_assert_eq!(back(d.neg(mx)), g(-Integer::from(x)));
        prop_assert_eq!(back(d.mul(mx, my)), g(Integer::from(x) * y));
        prop_assert_eq!(back(d.one()), Integer::from(1));
    }
}

/// Even moduli in `(2^64, 2^128)`: odd parts of one and two limbs, 1 (a
/// power of two), and the largest powers of two dividing them.
fn even_modulus() -> BoxedStrategy<u128> {
    prop_oneof![
        (any::<u128>(), 1..27u32).prop_map(|(n, s)| ((n | 1) << s) | 1 << 100),
        (1..64u32, any::<u64>()).prop_map(|(s, o)| u128::from(o | 1) << (64 - s + 1).min(64)),
        (65..128u32).prop_map(|s| 1u128 << s),
        Just(u128::MAX - 1),
    ]
    .prop_filter("even, above 2^64", |&m| m >> 64 != 0 && m & 1 == 0)
    .boxed()
}

proptest! {
    #[test]
    fn pow_mod_even_matches(
        m in even_modulus(),
        a in any::<u128>(),
        e in prop_oneof![Just(crate::testkit::Mag::default()), strategy::natural_up_to(4)],
    ) {
        let a = a % m;
        let ours = pow_mod_even(a, &e.to_natural(), m);
        let gmp = Integer::from(a).pow_mod(&e.to_rug(), &Integer::from(m)).unwrap();
        prop_assert_eq!(Integer::from(ours), gmp);
    }
}
