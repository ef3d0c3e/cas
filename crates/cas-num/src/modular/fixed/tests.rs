//! Fixed-size Montgomery arithmetic checked against GMP, for every length
//! from 2 to `FIXED_MAX` limbs.

use proptest::prelude::*;
use rug::Integer;

use super::{FIXED_MAX, Fixed, ModFixed};
use crate::Natural;
use crate::testkit::{Mag, ToRug, strategy};

/// An odd modulus of exactly `k` limbs: random, `B^k - 1` (where the
/// products' carries are largest), `B^(k-1) + 1`, or `2^(64k - 1) ± 1`.
fn modulus() -> BoxedStrategy<Integer> {
    (2..=FIXED_MAX)
        .prop_flat_map(|k| {
            let bits = u32::try_from(64 * k).unwrap();
            let random = strategy::natural_with_len(Just(k))
                .prop_map(|m| m.to_rug() | 1)
                .boxed();
            prop_oneof![
                4 => random,
                1 => Just((Integer::from(1) << bits) - 1),
                1 => Just((Integer::from(1) << (bits - 64)) + 1),
                1 => Just((Integer::from(1) << (bits - 1)) + 1),
                1 => Just((Integer::from(1) << bits) - 3),
            ]
        })
        .boxed()
}

/// A modulus and two operands below it: 0, 1, `n - 1`, `n - 2` or random.
fn operands() -> BoxedStrategy<(Integer, Integer, Integer)> {
    modulus()
        .prop_flat_map(|n| {
            let near = |n: &Integer| {
                let random = strategy::natural_up_to(FIXED_MAX)
                    .prop_map({
                        let n = n.clone();
                        move |a| a.to_rug() % &n
                    })
                    .boxed();
                prop_oneof![
                    1 => Just(Integer::new()),
                    1 => Just(Integer::from(1)),
                    1 => Just(n.clone() - 1),
                    1 => Just(n.clone() - 2),
                    4 => random,
                ]
            };
            (Just(n.clone()), near(&n), near(&n))
        })
        .boxed()
}

fn nat(a: &Integer) -> Natural {
    Mag::from_rug(a).to_natural()
}

/// Every operation of `ModFixed::<K>` on `a` and `b` modulo `n` (of `K`
/// limbs), against GMP.
fn check<const K: usize>(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    e: &Natural,
) -> Result<(), TestCaseError> {
    let ring = ModFixed::<K>::new(&nat(n));
    let (x, y) = (ring.elem(&nat(a)), ring.elem(&nat(b)));
    let reduce = |v: Integer| {
        let r = v % n;
        if r < 0 { r + n } else { r }
    };
    let gmp_pow = a.clone().pow_mod(&e.to_rug(), n).unwrap();
    for (op, z, want) in [
        ("elem", x, a.clone()),
        ("one", ring.one(), reduce(Integer::from(1))),
        ("add", ring.add(x, y), reduce(Integer::from(a + b))),
        ("sub", ring.sub(x, y), reduce(Integer::from(a - b))),
        ("neg", ring.neg(x), reduce(Integer::from(-a))),
        ("mul", ring.mul(x, y), reduce(Integer::from(a * b))),
        ("square", ring.square(x), reduce(Integer::from(a * a))),
        ("pow", ring.pow(x, e), gmp_pow),
    ] {
        // The representation is canonical (below n), not only congruent.
        prop_assert!(limbs(z) < *n, "{} not reduced", op);
        prop_assert_eq!(ring.value(z).to_rug(), want, "{}", op);
    }
    Ok(())
}

/// The limbs of `z`, as a number.
fn limbs<const K: usize>(z: Fixed<K>) -> Integer {
    Integer::from_digits(&z.0, rug::integer::Order::Lsf)
}

/// `a mod n` through `ModFixed::<K>`'s elements, for any `a`: the value,
/// and the representation, which must be that of the reduced value.
fn reduced<const K: usize>(n: &Integer, a: &Natural) -> Result<Natural, TestCaseError> {
    let ring = ModFixed::<K>::new(&nat(n));
    let x = ring.elem(a);
    let value = ring.value(x);
    prop_assert_eq!(x, ring.elem(&value));
    prop_assert!(limbs(x) < *n);
    Ok(value)
}

fn dispatch(n: &Integer, a: &Integer, b: &Integer, e: &Natural) -> Result<(), TestCaseError> {
    match n.significant_digits::<u64>() {
        2 => check::<2>(n, a, b, e),
        3 => check::<3>(n, a, b, e),
        4 => check::<4>(n, a, b, e),
        5 => check::<5>(n, a, b, e),
        6 => check::<6>(n, a, b, e),
        7 => check::<7>(n, a, b, e),
        8 => check::<8>(n, a, b, e),
        k => unreachable!("modulus of {k} limbs"),
    }
}

proptest! {
    #[test]
    fn ops_match((n, a, b) in operands(), e in strategy::natural_up_to(3)) {
        dispatch(&n, &a, &b, &e.to_natural())?;
    }

    #[test]
    fn elem_reduces((n, a) in modulus().prop_flat_map(|n| {
        let near = prop_oneof![
            Just(n.clone()),
            Just(Integer::from(&n + 1)),
            Just(Integer::from(&n * 2u32)),
            Just(Integer::from(&n << 64)),
            Just((Integer::from(1) << (64 * n.significant_digits::<u64>())) - 1),
            strategy::natural_up_to(2 * FIXED_MAX).prop_map(|a| a.to_rug()),
        ];
        (Just(n), near)
    })) {
        let a_nat = nat(&a);
        let got = match n.significant_digits::<u64>() {
            2 => reduced::<2>(&n, &a_nat)?,
            3 => reduced::<3>(&n, &a_nat)?,
            4 => reduced::<4>(&n, &a_nat)?,
            5 => reduced::<5>(&n, &a_nat)?,
            6 => reduced::<6>(&n, &a_nat)?,
            7 => reduced::<7>(&n, &a_nat)?,
            8 => reduced::<8>(&n, &a_nat)?,
            k => unreachable!("modulus of {k} limbs"),
        };
        prop_assert_eq!(got.to_rug(), a % &n);
    }
}

/// Products of the largest residue, `n - 1`, modulo the largest modulus of
/// each length, `B^k - 1`, and powers to `n - 1` and `n - 2`.
#[test]
fn extremes() {
    for k in 2..=FIXED_MAX {
        let n = (Integer::from(1) << u32::try_from(64 * k).unwrap()) - 1;
        let top = Integer::from(&n - 1);
        for e in [
            Integer::from(&n - 1),
            Integer::from(&n - 2),
            Integer::from(1) << 200,
        ] {
            dispatch(&n, &top, &top, &nat(&e)).unwrap();
        }
    }
}
