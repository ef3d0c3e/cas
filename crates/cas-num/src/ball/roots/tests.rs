//! Tests of [`polynomial_roots`]: polynomials with known roots, each ball
//! holding one root and `p` vanishing on it, the real ones found real, the
//! others apart from the real line; and the polynomials it refuses.

use proptest::prelude::*;

use super::{polynomial_roots, value};
use crate::ball::ComplexBall;
use crate::{Complex, Integer, Rational};

fn c(n: i64) -> Complex {
    Complex::from(n)
}

/// The product of `p` and `q`, coefficients from the lowest.
fn times(p: &[Complex], q: &[Complex]) -> Vec<Complex> {
    let mut r = vec![c(0); p.len() + q.len() - 1];
    for (i, a) in p.iter().enumerate() {
        for (j, b) in q.iter().enumerate() {
            r[i + j] = &r[i + j] + &(a * b);
        }
    }
    r
}

/// `∏ (x − r)`.
fn with_roots(roots: &[i64]) -> Vec<Complex> {
    roots
        .iter()
        .fold(vec![c(1)], |p, &r| times(&p, &[c(-r), c(1)]))
}

/// Whether `p` holds 0 on the ball `z`, at 256 bits.
fn vanishes(coeffs: &[Complex], z: &ComplexBall) -> bool {
    let p: Vec<ComplexBall> = coeffs
        .iter()
        .map(|c| ComplexBall::from_complex(c, 256))
        .collect();
    value(&p, z, 256).may_be_zero()
}

/// The roots at 53 bits, at most 4096, and that `p` vanishes on each.
fn roots(coeffs: &[Complex]) -> Vec<ComplexBall> {
    let roots = polynomial_roots(coeffs, 53, 4096).expect("roots");
    assert_eq!(roots.len(), coeffs.len() - 1);
    for z in &roots {
        assert!(vanishes(coeffs, z), "{z:?}");
    }
    roots
}

fn real_count(roots: &[ComplexBall]) -> usize {
    roots.iter().filter(|z| z.is_real()).count()
}

#[test]
fn quadratics() {
    // x² − 2: ±√2, real, each squared holding 2.
    let found = roots(&[c(-2), c(0), c(1)]);
    assert_eq!(real_count(&found), 2);
    for z in &found {
        assert!(z.mul(z, 256).contains_complex(&c(2)));
    }
    // x² + 1: ±i.
    let found = roots(&[c(1), c(0), c(1)]);
    assert_eq!(real_count(&found), 0);
    assert!(found.iter().any(|z| z.contains_complex(&Complex::I)));
    assert!(found.iter().any(|z| z.contains_complex(&-Complex::I)));
}

#[test]
fn integer_roots() {
    let found = roots(&with_roots(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]));
    assert_eq!(real_count(&found), 12);
    for r in 1..=12 {
        let r = c(r);
        assert_eq!(found.iter().filter(|z| z.contains_complex(&r)).count(), 1);
    }
}

#[test]
fn a_quintic() {
    // x⁵ − x − 1: one real root, 1.1673039782614187, and two pairs.
    let found = roots(&[c(-1), c(-1), c(0), c(0), c(0), c(1)]);
    assert_eq!(real_count(&found), 1);
    let real = found.iter().find(|z| z.is_real()).expect("a real root");
    assert!((real.re().mid().to_f64() - 1.167_303_978_261_418_7).abs() < 1e-12);
    assert!(real.re().rel_accuracy_bits() >= 50);
}

#[test]
fn quartics() {
    // x⁴ − 10x² + 1: ±√2 ± √3, real; x⁴ + 1: none real.
    assert_eq!(real_count(&roots(&[c(1), c(0), c(-10), c(0), c(1)])), 4);
    assert_eq!(real_count(&roots(&[c(1), c(0), c(0), c(0), c(1)])), 0);
}

#[test]
fn complex_coefficients() {
    // x² − i: ±(1 + i)/√2, each squared holding i; nothing is known real.
    let found = roots(&[-Complex::I, c(0), c(1)]);
    assert_eq!(real_count(&found), 0);
    for z in &found {
        assert!(z.mul(z, 256).contains_complex(&Complex::I));
    }
}

#[test]
fn a_root_at_zero() {
    // x³ + x: 0 exactly, ±i.
    let found = roots(&[c(0), c(1), c(0), c(1)]);
    assert!(found.contains(&ComplexBall::ZERO));
    assert_eq!(real_count(&found), 1);
}

#[test]
fn far_and_near_roots() {
    // (x − 10⁶)(x − 1/1000) times 1000, and x³ − 10⁻⁹.
    let found = roots(&times(&[c(-1_000_000), c(1)], &[c(-1), c(1000)]));
    assert_eq!(real_count(&found), 2);
    let small = Complex::from(Rational::new(
        Integer::from(-1),
        Integer::from(1_000_000_000),
    ));
    let found = roots(&[small, c(0), c(0), c(1)]);
    assert_eq!(real_count(&found), 1);
}

#[test]
fn what_is_refused() {
    assert!(polynomial_roots(&[], 53, 1024).is_none());
    assert!(polynomial_roots(&[c(3)], 53, 1024).is_none());
    assert!(polynomial_roots(&[c(1), c(1), c(0)], 53, 1024).is_none());
    // (x − 1)²: a double root is never apart.
    assert!(polynomial_roots(&[c(1), c(-2), c(1)], 53, 256).is_none());
    // x²: a double root at 0.
    assert!(polynomial_roots(&[c(0), c(0), c(1)], 53, 256).is_none());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Distinct integer roots and quadratics without real roots: each
    /// integer in one real root's ball, the quadratics' roots not real.
    #[test]
    fn known_roots(
        integers in proptest::collection::btree_set(-30i64..30, 1..6),
        quadratics in proptest::collection::vec((-6i64..6, 1i64..20), 0..3),
    ) {
        let integers: Vec<i64> = integers.into_iter().collect();
        let mut p = with_roots(&integers);
        // x² + b x + c with b² < 4c, distinct.
        let mut seen = Vec::new();
        for (b, c0) in quadratics {
            if b * b < 4 * c0 && !seen.contains(&(b, c0)) {
                seen.push((b, c0));
                p = times(&p, &[c(c0), c(b), c(1)]);
            }
        }
        let found = roots(&p);
        prop_assert_eq!(real_count(&found), integers.len());
        for r in &integers {
            let r = c(*r);
            prop_assert_eq!(
                found.iter().filter(|z| z.is_real() && z.contains_complex(&r)).count(),
                1
            );
        }
    }
}
