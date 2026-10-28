//! Tests of Taylor's method (CAS-392): the coefficients of `y′ = y` and of
//! `y′ = t`, and the orders.

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::super::rhs::Node;
use super::{coefficients, order};
use crate::{Context, Expr};

/// `y′ = y` from `y(0) = 1`: `1/k!`; `y′ = t` from `y(1/2) = 0`: `1/2 + s`
/// gives `0, 1/2, 1/2`.
#[test]
fn coefficient_recurrence() {
    let ctx = Context::new();
    let variables = [Expr::symbol("t"), Expr::symbol("y")];
    let f = Node::compile(&ctx, &Expr::symbol("y"), &variables, 128).unwrap();
    let one = Ball::from_integer(&Integer::ONE);
    let jets = coefficients(&[f], &Ball::ZERO, std::slice::from_ref(&one), 7, 128).unwrap();
    let mut factorial = one.clone();
    for (k, c) in jets[0].coeffs().iter().enumerate() {
        if k > 0 {
            factorial = factorial.mul(&Ball::from_integer(&Integer::from(k as u64)), 128);
        }
        assert!(
            c.overlaps(&one.div(&factorial, 128)),
            "coefficient {k}: {c:?}"
        );
    }
    let f = Node::compile(&ctx, &Expr::symbol("t"), &variables, 128).unwrap();
    let half = Ball::from_rational(&Rational::new(Integer::ONE, Integer::from(2)), 128);
    let jets = coefficients(&[f], &half, &[Ball::ZERO], 4, 128).unwrap();
    let expected = [Ball::ZERO, half.clone(), half, Ball::ZERO];
    for (c, e) in jets[0].coeffs().iter().zip(&expected) {
        assert!(c.overlaps(e), "{c:?} against {e:?}");
    }
}

/// The order grows with the precision, within `[8, 60]`.
#[test]
fn orders() {
    assert_eq!(order(8), 8);
    assert_eq!(order(53), 17);
    assert_eq!(order(1000), 60);
}
