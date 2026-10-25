//! Tests of multifactor Hensel lifting (CAS-181).

use cas_num::algebra::{EuclideanDomain, Integers, Ring, WordPrimeField};
use cas_num::{Integer, Natural};
use proptest::prelude::*;

use crate::cancel::Never;
use crate::dense::{Dense, DenseRing};
use crate::factor::finite::factor_squarefree;
use crate::factor::hensel::lift;

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

/// `p` modulo the prime of `fp`.
fn image(fp: &DenseRing<WordPrimeField>, p: &Dense<Integer>) -> Dense<u64> {
    fp.poly(
        p.coeffs()
            .iter()
            .map(|c| fp.base().elem_integer(c))
            .collect(),
    )
}

/// `p` with its coefficients reduced into `[0, m)`.
fn reduce(p: &Dense<Integer>, m: &Natural) -> Dense<Integer> {
    let m = Integer::from(m);
    zx().poly(p.coeffs().iter().map(|c| c.div_rem_euclid(&m).1).collect())
}

/// The lifts of the factors of `f` to `pᵏ`, checked: monic, coefficients in
/// `[0, pᵏ)`, each congruent to its factor modulo `p`, and `lc(f)` times
/// their product congruent to `f` modulo `pᵏ`.
fn check_lift(
    fp: &DenseRing<WordPrimeField>,
    f: &Dense<Integer>,
    factors: &[Dense<u64>],
    k: u64,
) -> Result<Vec<Dense<Integer>>, TestCaseError> {
    let zx = zx();
    let m = fp.characteristic().pow(k);
    let lifted = lift(f, factors, fp, k, &Never).unwrap();
    prop_assert_eq!(lifted.len(), factors.len());
    for (big, small) in lifted.iter().zip(factors) {
        prop_assert_eq!(big.lead(), Some(&Integer::ONE), "a monic lift");
        prop_assert_eq!(&reduce(big, &m), big, "coefficients below p^k");
        prop_assert_eq!(&image(fp, big), small, "congruent modulo p");
    }
    let product = lifted
        .iter()
        .fold(zx.constant(f.lead().unwrap().clone()), |acc, g| {
            zx.mul(&acc, g)
        });
    prop_assert_eq!(
        reduce(&product, &m),
        reduce(f, &m),
        "the product modulo p^k"
    );
    Ok(lifted)
}

/// `x² + 1` is `(x − 2)(x + 2)` modulo 5, and lifts to `(x − a)(x + a)` for
/// the 5-adic square root `a ≡ 2` of −1.
#[test]
fn square_roots_of_minus_one() {
    let fp = DenseRing::new(WordPrimeField::new(5).unwrap());
    let f = zx().from_i64s(&[1, 0, 1]);
    let factors = [fp.from_i64s(&[-2, 1]), fp.from_i64s(&[2, 1])];
    for k in [1, 2, 3, 7, 20] {
        let m = Integer::from(Natural::from(5u64).pow(k));
        let lifted = lift(&f, &factors, &fp, k, &Never).unwrap();
        let a = (-&lifted[0].coeffs()[0]).div_rem_euclid(&m).1;
        assert_eq!(
            (&(&a * &a) + &Integer::ONE).div_rem_euclid(&m).1,
            Integer::ZERO
        );
        assert_eq!(a.div_rem_euclid(&Integer::from(5)).1, Integer::from(2));
        assert_eq!(lifted[1].coeffs()[0], a);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(150))]

    /// The local factors of products of random polynomials, square-free
    /// modulo a small prime, lift to several exponents; the lifts to `pᵏ`
    /// reduce to those to `pʲ`, `j < k`; and products of groups of local
    /// factors lift as well.
    #[test]
    fn lifts_of_local_factors(
        parts in prop::collection::vec(prop::collection::vec(-20i64..20, 2..5), 1..5),
        lead in prop::sample::select(vec![1i64, -1, 2, 6, 35]),
        p in prop::sample::select(vec![3u64, 11, 101, 65_537]),
        k in 1u64..12,
        groups in prop::collection::vec(any::<bool>(), 12),
    ) {
        let zx = zx();
        let fp = DenseRing::new(WordPrimeField::new(p).unwrap());
        let f = parts.iter().fold(zx.constant(Integer::from(lead)), |acc, c| zx.mul(&acc, &zx.from_i64s(c)));
        let fbar = image(&fp, &f);
        prop_assume!(fbar.degree() == f.degree() && fbar.degree() > Some(0));
        prop_assume!(EuclideanDomain::gcd(&fp, &fbar, &fp.derivative(&fbar)).degree() == Some(0));
        let monic = EuclideanDomain::normal_form(&fp, &fbar).0;
        let factors = factor_squarefree(&fp, &monic, &Never).unwrap();
        let lifted = check_lift(&fp, &f, &factors, k)?;
        for j in 1..k {
            let m = Natural::from(p).pow(j);
            let lower = check_lift(&fp, &f, &factors, j)?;
            for (high, low) in lifted.iter().zip(&lower) {
                prop_assert_eq!(&reduce(high, &m), low, "lifts are unique");
            }
        }
        // Two groups of the local factors, when both are non-empty.
        let (first, second): (Vec<_>, Vec<_>) = factors.iter().zip(groups.iter().cycle()).partition(|(_, g)| **g);
        if !first.is_empty() && !second.is_empty() {
            let product = |v: &[(&Dense<u64>, &bool)]| v.iter().fold(fp.one(), |acc, (g, _)| fp.mul(&acc, g));
            check_lift(&fp, &f, &[product(&first), product(&second)], k)?;
        }
    }
}
