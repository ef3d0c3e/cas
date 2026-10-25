//! Tests of factoring multivariate polynomials over ℤ and ℚ (CAS-184).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use std::sync::atomic::AtomicBool;

use cas_num::algebra::{Integers, Rationals, Ring};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use crate::cancel::Cancelled;
use crate::dense::DenseRing;
use crate::factor::Factors;
use crate::sparse::{Sparse, SparseRing};

/// `c ∏ fᵢ^eᵢ`.
fn product(ring: &SparseRing<Integers>, f: &Factors<Integer, Sparse<Integer>>) -> Sparse<Integer> {
    f.factors()
        .iter()
        .fold(ring.constant(f.constant().clone()), |acc, (g, e)| {
            ring.mul(&acc, &ring.pow(g, *e))
        })
}

fn total_degree<E>(p: &Sparse<E>) -> u64 {
    (0..p.term_count())
        .map(|k| p.exponents(k).iter().sum())
        .max()
        .unwrap_or(0)
}

/// The order of `factor`: by total degree, then by terms from the leading
/// one, then by exponent.
fn sort(factors: &mut [(Sparse<Integer>, u64)]) {
    let key = |p: &Sparse<Integer>| {
        let terms: Vec<(Vec<u64>, Integer)> = p.terms().map(|(e, c)| (e, c.clone())).collect();
        (total_degree(p), terms)
    };
    factors.sort_by(|(a, e), (b, x)| (key(a), e).cmp(&(key(b), x)));
}

/// The factorisation of `p`, checked. It multiplies back. Its factors have
/// variables, are primitive with positive leading coefficients, distinct
/// and sorted, and each factors as itself.
fn check(
    ring: &SparseRing<Integers>,
    p: &Sparse<Integer>,
) -> Result<Factors<Integer, Sparse<Integer>>, TestCaseError> {
    let f = ring.factor(p);
    prop_assert_eq!(&product(ring, &f), p);
    let mut sorted = f.factors().to_vec();
    sort(&mut sorted);
    prop_assert_eq!(f.factors(), sorted.as_slice(), "sorted factors");
    for (i, (g, e)) in f.factors().iter().enumerate() {
        prop_assert!(*e >= 1);
        prop_assert!(total_degree(g) > 0, "a factor with a variable");
        prop_assert_eq!(ring.content(g), Integer::ONE, "a primitive factor");
        prop_assert!(g.coeffs()[0].is_positive());
        prop_assert!(
            f.factors()[i + 1..].iter().all(|(h, _)| h != g),
            "distinct factors"
        );
        let again = ring.factor(g);
        prop_assert_eq!(again.constant(), &Integer::ONE);
        prop_assert_eq!(
            again.factors(),
            &[(g.clone(), 1)][..],
            "an irreducible factor"
        );
    }
    Ok(f)
}

/// The factorisation of `∏ aᵢ^eᵢ` from those of the `aᵢ`.
fn union(
    ring: &SparseRing<Integers>,
    parts: &[(Sparse<Integer>, u64)],
) -> (Integer, Vec<(Sparse<Integer>, u64)>) {
    let mut constant = Integer::ONE;
    let mut factors: Vec<(Sparse<Integer>, u64)> = Vec::new();
    for (a, e) in parts {
        let f = ring.factor(a);
        constant = &constant * &f.constant().pow(*e);
        for (g, k) in f.factors() {
            match factors.iter_mut().find(|(h, _)| h == g) {
                Some((_, x)) => *x += k * e,
                None => factors.push((g.clone(), k * e)),
            }
        }
    }
    sort(&mut factors);
    (constant, factors)
}

/// Expects the factorisation `constant ∏ fᵢ^eᵢ` of their product.
fn expect(ring: &SparseRing<Integers>, constant: i64, factors: &[(Sparse<Integer>, u64)]) {
    let mut sorted = factors.to_vec();
    sort(&mut sorted);
    let p = sorted
        .iter()
        .fold(ring.constant(Integer::from(constant)), |acc, (g, e)| {
            ring.mul(&acc, &ring.pow(g, *e))
        });
    let f = ring.factor(&p);
    assert_eq!(f.constant(), &Integer::from(constant), "{p:?}");
    assert_eq!(f.factors(), sorted.as_slice(), "{p:?}");
}

/// giac's answers, leading coefficients, contents, powers and signs.
#[test]
fn examples() {
    let ring = SparseRing::new(Integers, 3);
    let (x, y, z) = (ring.var(0), ring.var(1), ring.var(2));
    let int = |n: i64| ring.constant(Integer::from(n));
    let sum = |ps: &[&Sparse<Integer>]| ps.iter().fold(ring.zero(), |acc, p| ring.add(&acc, p));
    let times = |a: &Sparse<Integer>, b: &Sparse<Integer>| ring.mul(a, b);
    let (xx, yy) = (times(&x, &x), times(&y, &y));
    // factor(x^4-y^4) = (x-y)*(x+y)*(x^2+y^2)
    expect(
        &ring,
        1,
        &[
            (ring.sub(&x, &y), 1),
            (ring.add(&x, &y), 1),
            (ring.add(&xx, &yy), 1),
        ],
    );
    // factor(x^2*y+x*y^2) = y*x*(x+y)
    expect(
        &ring,
        1,
        &[(x.clone(), 1), (y.clone(), 1), (ring.add(&x, &y), 1)],
    );
    // x^2 - 2 y^2 and x^2 + y^2 + 1 are irreducible over ℚ.
    expect(&ring, 1, &[(ring.sub(&xx, &times(&int(2), &yy)), 1)]);
    expect(&ring, 1, &[(sum(&[&xx, &yy, &int(1)]), 1)]);
    // Three variables.
    expect(
        &ring,
        1,
        &[
            (sum(&[&x, &y, &z]), 1),
            (ring.add(&ring.sub(&times(&x, &y), &z), &int(1)), 1),
            (sum(&[&xx, &times(&z, &y), &int(3)]), 1),
        ],
    );
    // Leading coefficients y and y^2 in x.
    expect(
        &ring,
        1,
        &[
            (ring.add(&times(&y, &x), &int(1)), 1),
            (sum(&[&times(&yy, &x), &y, &int(2)]), 1),
        ],
    );
    // A content y + 1 in x, powers, and a constant with a sign.
    expect(
        &ring,
        -6,
        &[
            (ring.add(&y, &int(1)), 1),
            (ring.add(&x, &y), 2),
            (ring.sub(&x, &times(&int(2), &y)), 3),
            (sum(&[&xx, &y, &int(1)]), 1),
        ],
    );
    // Constants.
    let f = ring.factor(&int(-12));
    assert_eq!(f.constant(), &Integer::from(-12));
    assert!(f.factors().is_empty());
    let f = ring.factor(&ring.zero());
    assert_eq!(f.constant(), &Integer::ZERO);
}

#[test]
fn cancellation() {
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let p = ring.sub(&ring.pow(&x, 4), &ring.pow(&y, 4));
    let stop = AtomicBool::new(true);
    assert_eq!(ring.factor_with(&p, &stop), Err(Cancelled));
}

/// Terms of `n` variables, small exponents and coefficients.
fn terms(n: usize, len: usize) -> impl Strategy<Value = Vec<(Vec<u64>, i64)>> {
    prop::collection::vec((prop::collection::vec(0u64..3, n), -5i64..6), 1..=len)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Products of powers of random polynomials in one to three variables:
    /// the checks, the union of the parts' factorisations, in one variable
    /// the dense factorisation, and over ℚ the monic associates.
    #[test]
    fn random_products(
        (n, parts) in (1usize..=3).prop_flat_map(|n| (Just(n), prop::collection::vec((terms(n, 4), 1u64..3), 1..4))),
    ) {
        let ring = SparseRing::new(Integers, n);
        let parts: Vec<(Sparse<Integer>, u64)> = parts
            .iter()
            .map(|(t, e)| (ring.from_terms(t.iter().map(|(x, c)| (x.clone(), Integer::from(*c)))), *e))
            .collect();
        prop_assume!(parts.iter().all(|(a, _)| !a.is_zero()));
        let p = parts.iter().fold(ring.one(), |acc, (a, e)| ring.mul(&acc, &ring.pow(a, *e)));
        let f = check(&ring, &p)?;
        let (constant, factors) = union(&ring, &parts);
        prop_assert_eq!(f.constant(), &constant);
        prop_assert_eq!(f.factors(), factors.as_slice());
        if n == 1 {
            let dense = DenseRing::new(Integers).factor(&ring.to_dense(&p, 0).unwrap());
            prop_assert_eq!(dense.constant(), f.constant());
            prop_assert_eq!(dense.factors().len(), f.factors().len());
            for (g, e) in dense.factors() {
                prop_assert!(f.factors().contains(&(ring.from_dense(g, 0), *e)));
            }
        }
        let qs = SparseRing::new(Rationals, n);
        let pq = ring.map_into(&p, &qs, |c| Rational::new(c.clone(), Integer::from(4)));
        let fq = qs.factor(&pq);
        prop_assert_eq!(fq.constant(), &pq.coeffs()[0]);
        let monic: Vec<(Sparse<Rational>, u64)> = f
            .factors()
            .iter()
            .map(|(g, e)| {
                let lead = g.coeffs()[0].clone();
                (ring.map_into(g, &qs, |c| Rational::new(c.clone(), lead.clone())), *e)
            })
            .collect();
        prop_assert_eq!(fq.factors(), monic.as_slice());
    }
}
