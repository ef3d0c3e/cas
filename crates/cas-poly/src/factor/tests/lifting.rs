//! Tests of multivariate Hensel lifting (CAS-183).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use cas_num::algebra::{Integers, Ring};
use cas_num::{Integer, Natural};
use proptest::prelude::*;

use crate::cancel::Never;
use crate::dense::{Dense, DenseRing};
use crate::factor::lifting::lift_variables;
use crate::sparse::{Sparse, SparseRing};

/// A prime above every coefficient here.
fn prime() -> Natural {
    Natural::power_of_two(80).next_prime()
}

/// The inputs of the lifting for known factors: their images at the point,
/// univariate in `main`, and their leading coefficients in `main`.
fn inputs(
    ring: &SparseRing<Integers>,
    factors: &[Sparse<Integer>],
    main: usize,
    points: &[(usize, Integer)],
) -> (Vec<Dense<Integer>>, Vec<Sparse<Integer>>) {
    let univariate = factors
        .iter()
        .map(|g| {
            let at = points
                .iter()
                .fold(g.clone(), |p, (v, a)| ring.eval_var(&p, *v, a));
            ring.to_dense(&at, main).unwrap()
        })
        .collect();
    let lcs = factors.iter().map(|g| ring.lead_coeff(g, main)).collect();
    (univariate, lcs)
}

/// Whether the images are pairwise coprime and keep their degrees in
/// `main`.
fn good_point(
    ring: &SparseRing<Integers>,
    factors: &[Sparse<Integer>],
    univariate: &[Dense<Integer>],
    main: usize,
) -> bool {
    let zx = DenseRing::new(Integers);
    factors
        .iter()
        .zip(univariate)
        .all(|(g, u)| u.degree().map(|d| d as u64) == ring.degree(g, main))
        && univariate.iter().enumerate().all(|(i, a)| {
            univariate[i + 1..]
                .iter()
                .all(|b| zx.gcd(a, b).degree() == Some(0))
        })
}

/// Lifts the images of known factors and expects them back.
fn lifts_back(
    ring: &SparseRing<Integers>,
    factors: &[Sparse<Integer>],
    main: usize,
    points: &[(usize, Integer)],
) -> Option<Vec<Sparse<Integer>>> {
    let f = factors.iter().fold(ring.one(), |acc, g| ring.mul(&acc, g));
    let (univariate, lcs) = inputs(ring, factors, main, points);
    lift_variables(&f, main, points, &univariate, &lcs, &prime(), &Never).unwrap()
}

#[test]
fn examples() {
    let ring = SparseRing::new(Integers, 3);
    let (x, y, z) = (ring.var(0), ring.var(1), ring.var(2));
    let int = |n: i64| ring.constant(Integer::from(n));
    let point = |v: usize, a: i64| (v, Integer::from(a));
    // (x + y)(x - y) at y = 2.
    let factors = [ring.add(&x, &y), ring.sub(&x, &y)];
    assert_eq!(
        lifts_back(&ring, &factors, 0, &[point(1, 2)]),
        Some(factors.to_vec())
    );
    // (y x + 1)(x + y^2) at y = 2: a leading coefficient y.
    let factors = [
        ring.add(&ring.mul(&y, &x), &int(1)),
        ring.add(&x, &ring.mul(&y, &y)),
    ];
    assert_eq!(
        lifts_back(&ring, &factors, 0, &[point(1, 2)]),
        Some(factors.to_vec())
    );
    // (x + y + z)(x y - z + 1)(x^2 + z y + 3) at y = 1, z = 2, then in the
    // other order, and with z as the main variable.
    let factors = [
        ring.add(&ring.add(&x, &y), &z),
        ring.add(&ring.sub(&ring.mul(&x, &y), &z), &int(1)),
        ring.add(&ring.add(&ring.mul(&x, &x), &ring.mul(&z, &y)), &int(3)),
    ];
    for points in [[point(1, 1), point(2, 2)], [point(2, 2), point(1, 1)]] {
        assert_eq!(
            lifts_back(&ring, &factors, 0, &points),
            Some(factors.to_vec())
        );
    }
    let factors = [
        ring.add(&ring.add(&x, &y), &z),
        ring.add(&ring.sub(&ring.mul(&z, &y), &x), &int(2)),
    ];
    assert_eq!(
        lifts_back(&ring, &factors, 2, &[point(0, 3), point(1, -1)]),
        Some(factors.to_vec())
    );
}

/// Leading coefficients that agree with the images but that no
/// factorisation has: the lifts do not multiply back.
#[test]
fn wrong_leading_coefficients() {
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let one = ring.one();
    // (y x + 1)(x + y^2) at y = 2 is (2x + 1)(x + 4), whose leading
    // coefficients 2 and 1 are not y and 1.
    let factors = [
        ring.add(&ring.mul(&y, &x), &one),
        ring.add(&x, &ring.mul(&y, &y)),
    ];
    let f = ring.mul(&factors[0], &factors[1]);
    let points = [(1, Integer::from(2))];
    let (univariate, _) = inputs(&ring, &factors, 0, &points);
    let lcs = [ring.constant(Integer::from(2)), one];
    let result = lift_variables(&f, 0, &points, &univariate, &lcs, &prime(), &Never).unwrap();
    assert_eq!(result, None);
}

/// Terms as exponents and coefficients.
type Terms = Vec<(Vec<u64>, i64)>;

/// Factors `c xᵏ + …` in `n` variables: `k`, the terms of the leading
/// coefficient in the other variables, and terms of lower degree in `x`.
fn factor(n: usize) -> impl Strategy<Value = (u64, Terms, Terms)> {
    (
        1u64..3,
        prop::collection::vec((prop::collection::vec(0u64..2, n - 1), 1i64..4), 1..3),
        prop::collection::vec((prop::collection::vec(0u64..3, n), -4i64..5), 0..4),
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    /// Products of two or three random factors in two or three variables
    /// lift back from their images at random points where those are
    /// coprime.
    #[test]
    fn random_factors_lift_back(
        (n, factors) in (2usize..=3).prop_flat_map(|n| (Just(n), prop::collection::vec(factor(n), 2..4))),
        values in prop::collection::vec(-3i64..4, 2),
    ) {
        let ring = SparseRing::new(Integers, n);
        let factors: Vec<Sparse<Integer>> = factors
            .iter()
            .map(|(k, lead, lower)| {
                // The leading coefficient in x₀, without x₀, times x₀ᵏ, plus
                // terms of lower degree in x₀.
                let mut terms: Vec<(Vec<u64>, Integer)> = lead
                    .iter()
                    .map(|(e, c)| {
                        let mut exps = vec![*k];
                        exps.extend(e);
                        (exps, Integer::from(*c))
                    })
                    .collect();
                terms.extend(lower.iter().filter(|(e, _)| e[0] < *k).map(|(e, c)| (e.clone(), Integer::from(*c))));
                ring.from_terms(terms)
            })
            .collect();
        let points: Vec<(usize, Integer)> = (1..n).map(|v| (v, Integer::from(values[v - 1]))).collect();
        let (univariate, _) = inputs(&ring, &factors, 0, &points);
        prop_assume!(good_point(&ring, &factors, &univariate, 0));
        prop_assert_eq!(lifts_back(&ring, &factors, 0, &points), Some(factors.clone()));
    }
}
