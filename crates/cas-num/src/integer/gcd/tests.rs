#![allow(
    clippy::many_single_char_names,
    reason = "a, b, g, s, t mirror the notation"
)]
//! `Integer` GCD, extended GCD, LCM and modular inverse checked against GMP.

use proptest::prelude::*;
use rug::Integer as Gmp;

use super::Integer;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Signed, ToRug, differential};

/// Pairs of signed values: GCD-shaped magnitudes with random signs, values at
/// the promotion edges of the inline representation, and small multiples of a
/// common value, which hit every exceptional cofactor case (`|a| = |b|`,
/// `|a| = 2g`, `|b| = 2g`, zero).
fn pair() -> BoxedStrategy<(Signed, Signed)> {
    let signs = || (any::<bool>(), any::<bool>());
    let small_multiples = (strategy::natural(), -6i64..=6, -6i64..=6).prop_map(|(g, i, j)| {
        let g = g.to_rug();
        let f = |k: i64| Signed::from_rug(&(Gmp::from(&g * k)));
        (f(i), f(j))
    });
    prop_oneof![
        4 => (strategy::gcd_pair(DEFAULT_MAX_LIMBS), signs()).prop_map(|((a, b), (sa, sb))| {
            (Signed::new(sa, a), Signed::new(sb, b))
        }),
        2 => (strategy::integer(), strategy::integer()),
        2 => small_multiples,
    ]
    .boxed()
}

fn ours(a: &Signed, b: &Signed) -> (Integer, Integer) {
    (a.to_integer(), b.to_integer())
}

differential! {
    fn gcd((a, b) in pair()) {
        ours: {
            let (x, y) = ours(&a, &b);
            let g = x.gcd(&y);
            prop_assert_eq!(&g, &y.gcd(&x));
            g
        },
        gmp: a.to_rug().gcd(&b.to_rug()),
    }

    fn lcm((a, b) in pair()) {
        ours: { let (x, y) = ours(&a, &b); x.lcm(&y) },
        gmp: a.to_rug().lcm(&b.to_rug()),
    }

    /// The exact cofactors GMP returns, and Bézout's identity checked
    /// directly.
    fn gcd_ext((a, b) in pair()) {
        ours: {
            let (x, y) = ours(&a, &b);
            let (g, s, t) = x.gcd_ext(&y);
            prop_assert_eq!(&(&s * &x + &t * &y), &g, "Bézout identity fails");
            (g, s, t)
        },
        gmp: a.to_rug().extended_gcd(b.to_rug(), Gmp::new()),
    }

    fn mod_inverse((a, m) in pair()) {
        ours: {
            let (x, y) = ours(&a, &m);
            let inv = x.mod_inverse(&y);
            if let Some(inv) = &inv {
                let one = (&(inv * &x) - &Integer::ONE).is_divisible_by(&y);
                prop_assert!(one, "not an inverse");
                prop_assert!(!inv.is_negative() && inv < &y.abs());
            }
            inv
        },
        gmp: {
            let m = m.to_rug();
            if m == 0 { None } else { a.to_rug().invert(&m).ok() }
        },
    }
}

#[test]
fn exceptional_cofactors() {
    let ext = |a: i64, b: i64| {
        let (g, s, t) = Integer::from(a).gcd_ext(&Integer::from(b));
        (
            g.to_i64().unwrap(),
            s.to_i64().unwrap(),
            t.to_i64().unwrap(),
        )
    };
    assert_eq!(ext(0, 0), (0, 0, 0));
    assert_eq!(ext(-5, 0), (5, -1, 0));
    assert_eq!(ext(0, -5), (5, 0, -1));
    assert_eq!(ext(7, -7), (7, 0, -1));
    assert_eq!(ext(21, 7), (7, 0, 1));
    assert_eq!(ext(-3, 6), (3, -1, 0));
    assert_eq!(ext(6, 4), (2, 1, -1));
    assert_eq!(ext(240, 46), (2, -9, 47));
}

#[test]
fn inverse_edge_cases() {
    let inv = |a: i64, m: i64| {
        Integer::from(a)
            .mod_inverse(&Integer::from(m))
            .map(|x| x.to_i64().unwrap())
    };
    assert_eq!(inv(3, 0), None);
    assert_eq!(inv(3, 1), Some(0));
    assert_eq!(inv(0, -1), Some(0));
    assert_eq!(inv(0, 7), None);
    assert_eq!(inv(3, 7), Some(5));
    assert_eq!(inv(-3, -7), Some(2));
    assert_eq!(inv(4, 6), None);
}

#[test]
fn gcd_of_i64_min_is_boxed() {
    let min = Integer::from(i64::MIN);
    let g = min.gcd(&Integer::ZERO);
    assert_eq!(g.to_rug(), Gmp::from(1u64 << 63));
    assert_eq!(min.gcd(&min).to_rug(), Gmp::from(1u64 << 63));
    let (g, s, t) = min.gcd_ext(&min);
    assert_eq!(
        (g.to_rug(), s, t),
        (Gmp::from(1u64 << 63), Integer::ZERO, Integer::from(-1))
    );
}
