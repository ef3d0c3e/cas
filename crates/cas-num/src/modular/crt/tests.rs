//! Chinese remaindering checked against GMP: a solution is unique in
//! `0..lcm`, so the congruences, the range and the modulus pin it down.

use proptest::prelude::*;
use rug::Integer;

use super::{CrtBasis, crt};
use crate::Natural;
use crate::testkit::{Mag, ToRug, strategy};

/// Non-zero moduli, often sharing a large factor.
fn moduli() -> BoxedStrategy<(Mag, Mag)> {
    strategy::gcd_pair(6)
        .prop_filter("non-zero", |(m, n)| !m.is_zero() && !n.is_zero())
        .boxed()
}

/// Word moduli made pairwise coprime by dropping any that share a factor
/// with an earlier one; up to `max` of them.
fn coprime_moduli(max: usize) -> BoxedStrategy<Vec<u64>> {
    prop::collection::vec(
        prop_oneof![1 => 1..100u64, 3 => strategy::top_limb()],
        0..=max,
    )
    .prop_map(|ms| {
        let mut kept: Vec<u64> = Vec::new();
        for m in ms {
            if kept
                .iter()
                .all(|&k| Integer::from(k).gcd(&Integer::from(m)) == 1)
            {
                kept.push(m);
            }
        }
        kept
    })
    .boxed()
}

fn product(ms: &[u64]) -> Integer {
    ms.iter().fold(Integer::from(1), |p, &m| p * m)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn crt_pair((m, n) in moduli(), a in strategy::natural_up_to(7), b in strategy::natural_up_to(7), same in prop::bool::weighted(0.4)) {
        // Equal inputs are always compatible, whatever the moduli share.
        let b = if same { a.clone() } else { b };
        let (mr, nr, ar, br) = (m.to_rug(), n.to_rug(), a.to_rug(), b.to_rug());
        let g = mr.clone().gcd(&nr);
        let compatible = (ar.clone() - &br).is_divisible(&g);
        match crt(&a.to_natural(), &m.to_natural(), &b.to_natural(), &n.to_natural()) {
            None => prop_assert!(!compatible, "compatible congruences rejected"),
            Some((x, l)) => {
                prop_assert!(compatible, "incompatible congruences solved");
                let (x, l) = (x.to_rug(), l.to_rug());
                prop_assert_eq!(&l, &mr.clone().lcm(&nr));
                prop_assert!(x < l, "solution not reduced");
                prop_assert_eq!(x.clone().modulo(&mr), ar.modulo(&mr));
                prop_assert_eq!(x.modulo(&nr), br.modulo(&nr));
            }
        }
    }

    #[test]
    fn basis_combine(ms in coprime_moduli(12), seed in prop::collection::vec(strategy::limb(), 12)) {
        let basis = CrtBasis::new(&ms).expect("pairwise coprime");
        let big_m = product(&ms);
        prop_assert_eq!(basis.modulus().to_rug(), big_m.clone());
        prop_assert_eq!(basis.moduli().collect::<Vec<_>>(), ms.clone());
        prop_assert_eq!(basis.len(), ms.len());
        let residues = &seed[..ms.len()];
        let x = basis.combine(residues).to_rug();
        prop_assert!(x < big_m, "not reduced");
        for (&m, &r) in ms.iter().zip(residues) {
            prop_assert_eq!(x.clone() % m, r % m);
        }
        let y = basis.combine_symmetric(residues).to_rug();
        prop_assert!(Integer::from(&y * 2) <= big_m && Integer::from(&y * 2) > -big_m.clone(), "not symmetric");
        prop_assert_eq!(y.modulo(&big_m), x);
    }

    /// A value below the product comes back from its residues.
    #[test]
    fn basis_round_trip(ms in coprime_moduli(40), x in strategy::natural_up_to(40)) {
        let basis = CrtBasis::new(&ms).expect("pairwise coprime");
        let x = x.to_rug().modulo(&product(&ms));
        let residues: Vec<u64> = ms.iter().map(|&m| (x.clone() % m).to_u64().unwrap()).collect();
        prop_assert_eq!(basis.combine(&residues).to_rug(), x);
    }

    /// Moduli sharing a factor are refused.
    #[test]
    fn basis_rejects_common_factor(ms in coprime_moduli(6), m in 2..u64::MAX, k in 2..1000u64, at in any::<prop::sample::Index>()) {
        let mut ms = ms;
        ms.insert(at.index(ms.len() + 1), m);
        let shared = m.checked_mul(k).unwrap_or(m);
        ms.push(shared);
        prop_assert!(CrtBasis::new(&ms).is_none());
    }
}

#[test]
fn edge_cases() {
    let empty = CrtBasis::new(&[]).unwrap();
    assert!(empty.is_empty());
    assert_eq!(empty.modulus(), &Natural::one());
    assert_eq!(empty.combine(&[]), Natural::ZERO);
    assert_eq!(empty.combine_symmetric(&[]), crate::Integer::ZERO);
    let ones = CrtBasis::new(&[1, 1, 7]).unwrap();
    assert_eq!(ones.combine(&[5, 9, 10]), Natural::from(3u64));
    let even = CrtBasis::new(&[4, 3]).unwrap();
    // 6 = M/2 is the top of the symmetric range, 7 wraps to -5.
    assert_eq!(even.combine_symmetric(&[2, 0]), crate::Integer::from(6));
    assert_eq!(even.combine_symmetric(&[3, 1]), crate::Integer::from(-5));
    let (x, l) = crt(
        &Natural::from(3u64),
        &Natural::from(4u64),
        &Natural::from(5u64),
        &Natural::from(6u64),
    )
    .unwrap();
    assert_eq!((x, l), (Natural::from(11u64), Natural::from(12u64)));
    assert!(
        crt(
            &Natural::from(0u64),
            &Natural::from(4u64),
            &Natural::from(1u64),
            &Natural::from(6u64)
        )
        .is_none()
    );
}

#[test]
#[should_panic(expected = "one residue per modulus")]
fn residue_count_mismatch_panics() {
    let _ = CrtBasis::new(&[3, 5]).unwrap().combine(&[1]);
}

#[test]
#[should_panic(expected = "zero modulus")]
fn zero_word_modulus_panics() {
    let _ = CrtBasis::new(&[3, 0]);
}

#[test]
#[should_panic(expected = "zero modulus")]
fn zero_crt_modulus_panics() {
    let _ = crt(
        &Natural::one(),
        &Natural::ZERO,
        &Natural::one(),
        &Natural::one(),
    );
}
