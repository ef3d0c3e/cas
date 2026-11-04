//! Tests of the certified eigenpairs (CAS-498): the balls hold the exact
//! eigenpair, are narrow, and a double eigenvalue fails.

use cas_num::ball::{Ball, ComplexBall, Mag};
use cas_num::{Complex, Integer, Rational};

use super::{attempt, eigenpair};

const PREC: u64 = 128;

fn matrix(rows: &[&[i64]]) -> Vec<Vec<ComplexBall>> {
    rows.iter()
        .map(|r| {
            r.iter()
                .map(|&x| ComplexBall::from_real(Ball::from_integer(&Integer::from(x))))
                .collect()
        })
        .collect()
}

/// `p/q` as a real ball.
fn real(p: i64, q: i64) -> ComplexBall {
    ComplexBall::from_real(Ball::from_rational(
        &Rational::new(p.into(), q.into()),
        PREC,
    ))
}

fn narrow(z: &ComplexBall) -> bool {
    let limit = Mag::power_of_two(&Integer::from(-90));
    *z.re().rad() < limit && *z.im().rad() < limit
}

#[test]
fn symmetric() {
    // [[2, 1], [1, 3]]: λ = (5 ± √5)/2, x = (1, λ − 2) scaled.
    let a = matrix(&[&[2, 1], &[1, 3]]);
    let (lambda, x) = attempt(&a, &real(3618, 1000), PREC).unwrap_or_else(|e| panic!("{e}"));
    assert!(lambda.is_real() && narrow(&lambda), "{lambda:?}");
    // λ² − 5λ + 5 holds 0.
    let p = lambda
        .mul(&lambda, PREC)
        .sub(&lambda.mul(&real(5, 1), PREC), PREC)
        .add(&real(5, 1), PREC);
    assert!(p.may_be_zero(), "{p:?}");
    // The largest component is 1: x = ((λ − 2)⁻¹... ) with x₂ = 1 here.
    assert!(
        x[1].contains_complex(&Complex::from(Rational::ONE)),
        "{x:?}"
    );
    assert!(x.iter().all(narrow), "{x:?}");
    // A x − λ x holds 0.
    for i in 0..2 {
        let mut r = a[i][0]
            .mul(&x[0], PREC)
            .add(&a[i][1].mul(&x[1], PREC), PREC);
        r = r.sub(&lambda.mul(&x[i], PREC), PREC);
        assert!(r.may_be_zero(), "{r:?}");
    }
}

#[test]
fn rotation() {
    // [[0, −1], [1, 0]]: λ = i, x = (i, 1) up to scale.
    let a = matrix(&[&[0, -1], &[1, 0]]);
    let guess = ComplexBall::new(
        Ball::from_rational(&Rational::new(1.into(), 100.into()), PREC),
        Ball::from_rational(&Rational::new(99.into(), 100.into()), PREC),
    );
    let (lambda, x) = attempt(&a, &guess, PREC).unwrap_or_else(|e| panic!("{e}"));
    let i = Complex::new(Rational::ZERO, Rational::ONE);
    assert!(lambda.contains_complex(&i) && narrow(&lambda), "{lambda:?}");
    let ratio = x[0].div(&x[1], PREC);
    assert!(ratio.contains_complex(&i), "{x:?}");
}

#[test]
fn double_eigenvalue_fails() {
    let a = matrix(&[&[1, 0], &[0, 1]]);
    assert!(eigenpair(&a, &real(1, 1), PREC).is_none());
    let jordan = matrix(&[&[2, 1], &[0, 2]]);
    assert!(eigenpair(&jordan, &real(2, 1), PREC).is_none());
}
