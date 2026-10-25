//! Tests of the integration of rational functions (CAS-196): examples,
//! Sturm's count of real roots, and properties: the derivative of the
//! result and the part left add up to the integrand, and denominators that
//! split over ℚ leave nothing.

use cas_num::algebra::{Rationals, Ring};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use super::{RationalIntegral, real_roots};
use crate::dense::{Dense, DenseRing};

fn qx() -> DenseRing<Rationals> {
    DenseRing::new(Rationals)
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn poly(coeffs: &[i64]) -> Dense<Rational> {
    qx().poly(coeffs.iter().map(|&c| q(c, 1)).collect())
}

/// `n / d` as a fraction plus `m / e`.
fn add(
    ring: &DenseRing<Rationals>,
    (n, d): (Dense<Rational>, Dense<Rational>),
    (m, e): (Dense<Rational>, Dense<Rational>),
) -> (Dense<Rational>, Dense<Rational>) {
    (
        ring.add(&ring.mul(&n, &e), &ring.mul(&m, &d)),
        ring.mul(&d, &e),
    )
}

/// The derivative of the result, plus the part left: a fraction.
fn derivative(
    ring: &DenseRing<Rationals>,
    r: &RationalIntegral,
) -> (Dense<Rational>, Dense<Rational>) {
    let mut sum = (ring.derivative(&r.polynomial), ring.one());
    // (n / d)′ = (n′ d − n d′) / d².
    let (n, d) = &r.rational;
    let numerator = ring.sub(
        &ring.mul(&ring.derivative(n), d),
        &ring.mul(n, &ring.derivative(d)),
    );
    sum = add(ring, sum, (numerator, ring.mul(d, d)));
    for log in &r.logarithms {
        let numerator = ring.scale(&ring.derivative(&log.argument), &log.coefficient);
        sum = add(ring, sum, (numerator, log.argument.clone()));
    }
    if let Some((p, e)) = &r.rest {
        sum = add(ring, sum, (p.clone(), e.clone()));
    }
    sum
}

/// The derivative of `∫ a / d` and the part left add up to `a / d`.
fn check(a: &Dense<Rational>, d: &Dense<Rational>) -> Result<RationalIntegral, TestCaseError> {
    let ring = qx();
    let r = ring.integrate_rational(a, d);
    let (n, e) = derivative(&ring, &r);
    prop_assert_eq!(
        ring.mul(&n, d),
        ring.mul(a, &e),
        "{:?} / {:?}: {:?}",
        a,
        d,
        r
    );
    // The rational part in lowest terms over a monic denominator; the
    // arguments monic and square-free; the part left proper.
    let (gn, gd) = &r.rational;
    prop_assert!(ring.is_one(&ring.gcd(gn, gd)) || gn.is_zero());
    prop_assert_eq!(gd.lead(), Some(&Rational::ONE));
    for log in &r.logarithms {
        prop_assert_eq!(log.argument.lead(), Some(&Rational::ONE));
        let s = &log.argument;
        prop_assert!(ring.is_one(&ring.gcd(s, &ring.derivative(s))));
        prop_assert_eq!(log.real_root, real_roots(&ring, s) > 0);
    }
    if let Some((p, e)) = &r.rest {
        prop_assert!(p.degree() < e.degree());
    }
    Ok(r)
}

#[test]
fn examples() {
    let ring = qx();
    // ∫ 1/x² = −1/x.
    let r = check(&poly(&[1]), &poly(&[0, 0, 1])).unwrap();
    assert_eq!(r.rational, (poly(&[-1]), poly(&[0, 1])));
    assert!(r.logarithms.is_empty() && r.rest.is_none());
    // ∫ 1/(x² (x + 1)) = −1/x − ln|x| + ln|x + 1|.
    let r = check(&poly(&[1]), &poly(&[0, 0, 1, 1])).unwrap();
    assert_eq!(r.rational, (poly(&[-1]), poly(&[0, 1])));
    let logs: Vec<(Rational, Dense<Rational>, bool)> = r
        .logarithms
        .iter()
        .map(|l| (l.coefficient.clone(), l.argument.clone(), l.real_root))
        .collect();
    assert!(logs.contains(&(q(-1, 1), poly(&[0, 1]), true)));
    assert!(logs.contains(&(q(1, 1), poly(&[1, 1]), true)));
    // ∫ 1/(x² − 1) = ln|x − 1|/2 − ln|x + 1|/2.
    let r = check(&poly(&[1]), &poly(&[-1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 2);
    assert!(r.rest.is_none());
    // ∫ 2x/(x² + 1) = ln(x² + 1): no real root, no absolute value.
    let r = check(&poly(&[0, 2]), &poly(&[1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 1);
    assert_eq!(r.logarithms[0].argument, poly(&[1, 0, 1]));
    assert!(!r.logarithms[0].real_root);
    // ∫ 1/(x² + 1): residues ±i/2, left for CAS-197.
    let r = check(&poly(&[1]), &poly(&[1, 0, 1])).unwrap();
    assert!(r.logarithms.is_empty());
    assert_eq!(r.rest, Some((poly(&[1]), poly(&[1, 0, 1]))));
    // ∫ (x³ + 1)/(x − 1) = x³/3 + x²/2 + x + 2 ln|x − 1|.
    let r = check(&poly(&[1, 0, 0, 1]), &poly(&[-1, 1])).unwrap();
    assert_eq!(
        r.polynomial,
        ring.poly(vec![q(0, 1), q(1, 1), q(1, 2), q(1, 3)])
    );
    // ∫ 1/(x (x² + 1)): ln|x| and −ln(x² + 1)/2 together, residues 1 and
    // −1/2 twice.
    let r = check(&poly(&[1]), &poly(&[0, 1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 2);
    assert!(r.rest.is_none());
    // ∫ 1/((x − 1)(x² + 1)): ln|x − 1|/2 and the rest with x² + 1.
    let r = check(&poly(&[1]), &poly(&[-1, 1, -1, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 2);
    let r = check(&poly(&[1]), &poly(&[-1, 2, -1, 0, 1])).unwrap();
    assert!(r.rest.is_some());
    // ∫ 1/(x² + 1)² = x/(2(x² + 1)) + ∫ 1/(2(x² + 1)).
    let r = check(&poly(&[1]), &poly(&[1, 0, 2, 0, 1])).unwrap();
    assert_eq!(
        r.rational,
        (ring.poly(vec![q(0, 1), q(1, 2)]), poly(&[1, 0, 1]))
    );
    assert_eq!(r.rest, Some((ring.poly(vec![q(1, 2)]), poly(&[1, 0, 1]))));
}

#[test]
fn sturm_counts_real_roots() {
    let ring = qx();
    for (p, n) in [
        (poly(&[1, 0, 1]), 0),
        (poly(&[-1, 0, 1]), 2),
        (poly(&[1, 0, -1, 1]), 1),
        (poly(&[0, -1, 0, 1]), 3),
        (poly(&[1, 0, 0, 0, 1]), 0),
        (poly(&[-2, 0, 0, 0, 1]), 2),
        (poly(&[5]), 0),
        (poly(&[3, 1]), 1),
    ] {
        assert_eq!(real_roots(&ring, &p), n, "{p:?}");
    }
}

/// A product of `(x − r)^e` and of `x² + c` for `c > 0`, times a constant.
fn denominator() -> impl Strategy<Value = (Vec<(i64, u32)>, Vec<i64>, i64)> {
    (
        prop::collection::vec((-6i64..7, 1u32..4), 0..4),
        prop::collection::vec(1i64..6, 0..2),
        prop_oneof![Just(1i64), Just(-2), Just(3)],
    )
}

fn build(roots: &[(i64, u32)], quadratics: &[i64], k: i64) -> Dense<Rational> {
    let ring = qx();
    let mut d = poly(&[k]);
    for &(r, e) in roots {
        for _ in 0..e {
            d = ring.mul(&d, &poly(&[-r, 1]));
        }
    }
    for &c in quadratics {
        d = ring.mul(&d, &poly(&[c, 0, 1]));
    }
    d
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    /// Any numerator over products of linear factors with multiplicities
    /// and positive quadratics: the derivative checks, and without the
    /// quadratics nothing is left.
    #[test]
    fn derivatives_add_up(
        a in prop::collection::vec(-9i64..10, 0..8),
        (roots, quadratics, k) in denominator(),
    ) {
        let d = build(&roots, &quadratics, k);
        let r = check(&poly(&a), &d)?;
        if quadratics.is_empty() {
            prop_assert!(r.rest.is_none(), "{:?}", r.rest);
        }
    }
}
