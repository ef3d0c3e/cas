//! Tests of factoring univariate polynomials over ℤ and ℚ (CAS-182).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use std::sync::atomic::AtomicBool;

use cas_num::algebra::{Integers, Rationals, Ring};
use cas_num::{Integer, Natural, Rational};
use proptest::prelude::*;

use crate::cancel::Cancelled;
use crate::dense::{Dense, DenseRing};
use crate::factor::Factors;

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

/// `c ∏ fᵢ^eᵢ`.
fn product(ring: &DenseRing<Integers>, f: &Factors<Integer, Dense<Integer>>) -> Dense<Integer> {
    f.factors()
        .iter()
        .fold(ring.constant(f.constant().clone()), |acc, (g, e)| {
            ring.mul(&acc, &ring.pow(g, *e))
        })
}

/// The order of `factor`: by degree, then by coefficients from the
/// constant term up, then by exponent.
fn sort(factors: &mut [(Dense<Integer>, u64)]) {
    factors.sort_by(|(a, e), (b, x)| (a.degree(), a.coeffs(), e).cmp(&(b.degree(), b.coeffs(), x)));
}

/// The factorisation of `p`, checked. It multiplies back. Its factors have
/// positive degrees, are primitive with positive leading coefficients, are
/// distinct and sorted, and each factors as itself.
fn check(p: &Dense<Integer>) -> Result<Factors<Integer, Dense<Integer>>, TestCaseError> {
    let ring = zx();
    let f = ring.factor(p);
    prop_assert_eq!(&product(&ring, &f), p);
    for pair in f.factors().windows(2) {
        let ((a, _), (b, _)) = (&pair[0], &pair[1]);
        prop_assert!(
            (a.degree(), a.coeffs()) < (b.degree(), b.coeffs()),
            "sorted, distinct factors"
        );
    }
    for (g, e) in f.factors() {
        prop_assert!(*e >= 1);
        prop_assert!(g.degree() > Some(0), "a factor of positive degree");
        prop_assert_eq!(ring.content(g), Integer::ONE, "a primitive factor");
        prop_assert!(g.lead().unwrap().is_positive());
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

/// The factorisation of `∏ aᵢ^eᵢ` from those of the `aᵢ`: the product of
/// their constants, and their factors with the exponents summed.
fn union(parts: &[(Dense<Integer>, u64)]) -> (Integer, Vec<(Dense<Integer>, u64)>) {
    let ring = zx();
    let mut constant = Integer::ONE;
    let mut factors: Vec<(Dense<Integer>, u64)> = Vec::new();
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

/// The cyclotomic polynomial `Φₙ`: `xⁿ − 1` over the `Φ_d` of the proper
/// divisors `d` of `n`.
fn cyclotomic(n: usize) -> Dense<Integer> {
    let ring = zx();
    let mut p = ring.sub(&ring.monomial(Integer::ONE, n), &ring.one());
    for d in (1..n).filter(|d| n % d == 0) {
        p = ring.div_exact(&p, &cyclotomic(d)).unwrap();
    }
    p
}

/// A case: coefficients from the constant term up, the constant, and the
/// factors with their exponents.
type Case = (Vec<i64>, i64, Vec<(Vec<i64>, u64)>);

/// giac's answers, and polynomials that split modulo every prime.
#[test]
fn examples() {
    let ring = zx();
    let cases: Vec<Case> = vec![
        // factor(x^6-1) = (x-1)*(x+1)*(x^2-x+1)*(x^2+x+1)
        (
            vec![-1, 0, 0, 0, 0, 0, 1],
            1,
            vec![
                (vec![-1, 1], 1),
                (vec![1, 1], 1),
                (vec![1, -1, 1], 1),
                (vec![1, 1, 1], 1),
            ],
        ),
        // factor(2*x^2-2) = 2*(x-1)*(x+1)
        (vec![-2, 0, 2], 2, vec![(vec![-1, 1], 1), (vec![1, 1], 1)]),
        // factor(-x^2+1) = -(x-1)*(x+1)
        (vec![1, 0, -1], -1, vec![(vec![-1, 1], 1), (vec![1, 1], 1)]),
        // factor(6*x+12) = 6*(x+2)
        (vec![12, 6], 6, vec![(vec![2, 1], 1)]),
        // x^4 - y^4 at y = 1: (x-1)*(x+1)*(x^2+1)
        (
            vec![-1, 0, 0, 0, 1],
            1,
            vec![(vec![-1, 1], 1), (vec![1, 1], 1), (vec![1, 0, 1], 1)],
        ),
        // x^4 + 4 = (x^2-2x+2)(x^2+2x+2), Sophie Germain's identity.
        (
            vec![4, 0, 0, 0, 1],
            1,
            vec![(vec![2, -2, 1], 1), (vec![2, 2, 1], 1)],
        ),
        // x^4 + x^2 + 1 = (x^2-x+1)(x^2+x+1)
        (
            vec![1, 0, 1, 0, 1],
            1,
            vec![(vec![1, -1, 1], 1), (vec![1, 1, 1], 1)],
        ),
        // x^3 - x^2 = x^2 (x - 1)
        (
            vec![0, 0, -1, 1],
            1,
            vec![(vec![-1, 1], 1), (vec![0, 1], 2)],
        ),
        // (x^2 + 1)^3 (x - 1)^2
        (
            vec![1, -2, 4, -6, 6, -6, 4, -2, 1],
            1,
            vec![(vec![-1, 1], 2), (vec![1, 0, 1], 3)],
        ),
        // Irreducible: x^2 - 2, 4x^2 - 8x + 1, x^4 + 1, and the
        // Swinnerton-Dyer polynomials of 2 and 3, and of 2, 3 and 5, which
        // split into factors of degree at most 2 modulo every prime.
        (vec![-2, 0, 1], 1, vec![(vec![-2, 0, 1], 1)]),
        (vec![1, -8, 4], 1, vec![(vec![1, -8, 4], 1)]),
        (vec![1, 0, 0, 0, 1], 1, vec![(vec![1, 0, 0, 0, 1], 1)]),
        (vec![1, 0, -10, 0, 1], 1, vec![(vec![1, 0, -10, 0, 1], 1)]),
        (
            vec![576, 0, -960, 0, 352, 0, -40, 0, 1],
            1,
            vec![(vec![576, 0, -960, 0, 352, 0, -40, 0, 1], 1)],
        ),
        // Constants.
        (vec![-5], -5, vec![]),
        (vec![], 0, vec![]),
    ];
    for (input, constant, factors) in cases {
        let f = ring.factor(&ring.from_i64s(&input));
        assert_eq!(f.constant(), &Integer::from(constant), "{input:?}");
        let expected: Vec<(Dense<Integer>, u64)> = factors
            .iter()
            .map(|(g, e)| (ring.from_i64s(g), *e))
            .collect();
        assert_eq!(f.factors(), expected.as_slice(), "{input:?}");
    }
}

/// The Swinnerton-Dyer polynomial of 2, 3, 5 and 7, of degree 16, is
/// irreducible, and splits into 8 factors modulo every prime.
#[test]
fn swinnerton_dyer() {
    let ring = zx();
    let s = ring.from_i64s(&[
        46225, 0, -5_596_840, 0, 13_950_764, 0, -7_453_176, 0, 1_513_334, 0, -141_912, 0, 6476, 0,
        -136, 0, 1,
    ]);
    let f = ring.factor(&s);
    assert_eq!(f.factors(), &[(s.clone(), 1)]);
    // Times x^2 - 2, whose roots are among its local factors' too.
    let q = ring.from_i64s(&[-2, 0, 1]);
    let f = ring.factor(&ring.mul(&s, &q));
    assert_eq!(f.factors(), &[(q, 1), (s, 1)]);
}

/// `xⁿ − 1` is the product of the `Φ_d` for the divisors `d` of `n`.
#[test]
fn cyclotomic_polynomials() {
    let ring = zx();
    for n in 1..=30 {
        let p = ring.sub(&ring.monomial(Integer::ONE, n), &ring.one());
        let mut expected: Vec<(Dense<Integer>, u64)> = (1..=n)
            .filter(|d| n % d == 0)
            .map(|d| (cyclotomic(d), 1))
            .collect();
        sort(&mut expected);
        let f = ring.factor(&p);
        assert_eq!(f.constant(), &Integer::ONE);
        assert_eq!(f.factors(), expected.as_slice(), "x^{n} - 1");
    }
}

/// Factors with coefficients of hundreds of bits, and many linear factors.
#[test]
fn large_coefficients_and_many_factors() {
    let ring = zx();
    let big = Integer::from(Natural::power_of_two(100));
    let a = ring.poly(vec![-big, Integer::ONE]);
    let b = ring.poly(vec![Integer::ONE, Integer::from(3).pow(50)]);
    let c = ring.from_i64s(&[7, 0, 0, 1]);
    let f = ring.factor(&ring.mul(&ring.mul(&a, &b), &c));
    assert_eq!(f.factors(), &[(a, 1), (b, 1), (c, 1)]);
    let linear: Vec<Dense<Integer>> = (1..=12).map(|i| ring.from_i64s(&[-i, 1])).collect();
    let p = linear.iter().fold(ring.one(), |acc, l| ring.mul(&acc, l));
    let mut expected: Vec<(Dense<Integer>, u64)> = linear.into_iter().map(|l| (l, 1)).collect();
    sort(&mut expected);
    assert_eq!(ring.factor(&p).factors(), expected.as_slice());
}

#[test]
fn cancellation() {
    let ring = zx();
    let p = ring.from_i64s(&[576, 0, -960, 0, 352, 0, -40, 0, 1]);
    let stop = AtomicBool::new(true);
    assert_eq!(ring.factor_with(&p, &stop), Err(Cancelled));
}

/// Polynomials irreducible by Eisenstein's criterion for a prime `q`: `q`
/// divides every coefficient but the leading one, and `q²` does not divide
/// the constant term. Their primitive parts are.
fn eisenstein() -> impl Strategy<Value = Vec<i64>> {
    (
        prop::sample::select(vec![2i64, 3, 5, 7]),
        prop::collection::vec(-4i64..4, 0..5),
        1i64..4,
        1i64..6,
    )
        .prop_map(|(q, middle, constant, lead)| {
            let constant = q * (q * constant + 1);
            let lead = q * lead + 1;
            let mut c = vec![constant];
            c.extend(middle.iter().map(|m| q * m));
            c.push(lead);
            c
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(150))]

    /// Products of powers of random polynomials: the checks, and the union
    /// of the parts' factorisations; over ℚ, the monic associates with the
    /// leading coefficient.
    #[test]
    fn random_products(parts in prop::collection::vec((prop::collection::vec(-30i64..30, 2..6), 1u64..3), 1..4), d in 1i64..5) {
        let ring = zx();
        let parts: Vec<(Dense<Integer>, u64)> = parts.iter().map(|(c, e)| (ring.from_i64s(c), *e)).collect();
        prop_assume!(parts.iter().all(|(a, _)| !a.is_zero()));
        let p = parts.iter().fold(ring.one(), |acc, (a, e)| ring.mul(&acc, &ring.pow(a, *e)));
        let f = check(&p)?;
        let (constant, factors) = union(&parts);
        prop_assert_eq!(f.constant(), &constant);
        prop_assert_eq!(f.factors(), factors.as_slice());
        // Over ℚ: p / d.
        let qx = DenseRing::new(Rationals);
        let pq = qx.poly(p.coeffs().iter().map(|c| Rational::new(c.clone(), Integer::from(d))).collect());
        let fq = qx.factor(&pq);
        prop_assert_eq!(fq.constant(), pq.lead().unwrap());
        let monic: Vec<(Dense<Rational>, u64)> = f.factors().iter().map(|(g, e)| {
            let lead = g.lead().unwrap();
            (qx.poly(g.coeffs().iter().map(|c| Rational::new(c.clone(), lead.clone())).collect()), *e)
        }).collect();
        prop_assert_eq!(fq.factors(), monic.as_slice());
    }

    /// Products of distinct Eisenstein polynomials factor into exactly
    /// them.
    #[test]
    fn eisenstein_products(parts in prop::collection::vec(eisenstein(), 1..4), k in 1i64..9) {
        let ring = zx();
        let mut expected: Vec<(Dense<Integer>, u64)> = Vec::new();
        for c in &parts {
            let g = ring.primitive_part(&ring.from_i64s(c));
            if !expected.iter().any(|(h, _)| *h == g) {
                expected.push((g, 1));
            }
        }
        let p = expected.iter().fold(ring.constant(Integer::from(k)), |acc, (g, _)| ring.mul(&acc, g));
        sort(&mut expected);
        let f = check(&p)?;
        prop_assert_eq!(f.constant(), &Integer::from(k));
        prop_assert_eq!(f.factors(), expected.as_slice());
    }
}
