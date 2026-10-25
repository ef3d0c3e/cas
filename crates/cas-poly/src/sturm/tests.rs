//! Tests of Sturm's count of real roots (CAS-197, CAS-198): examples on the
//! line and between bounds, roots at the bounds, repeated roots; and a
//! proptest against products of linear factors with known roots.

use cas_num::algebra::{Rationals, Ring};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use crate::dense::{Dense, DenseRing};

fn qx() -> DenseRing<Rationals> {
    DenseRing::new(Rationals)
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn poly(coeffs: &[i64]) -> Dense<Rational> {
    qx().from_i64s(coeffs)
}

#[test]
fn on_the_line() {
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
        // (x − 1)² (x + 1): distinct roots.
        (poly(&[1, -1, -1, 1]), 2),
    ] {
        assert_eq!(ring.count_real_roots(&p, None, None), n, "{p:?}");
    }
}

#[test]
fn between_bounds() {
    let ring = qx();
    let two = poly(&[-2, 0, 1]);
    let count = |p: &Dense<Rational>, lo: Option<Rational>, hi: Option<Rational>| {
        ring.count_real_roots(p, lo.as_ref(), hi.as_ref())
    };
    // ±√2.
    assert_eq!(count(&two, Some(q(0, 1)), Some(q(1, 1))), 0);
    assert_eq!(count(&two, Some(q(1, 1)), Some(q(2, 1))), 1);
    assert_eq!(count(&two, Some(q(-2, 1)), Some(q(2, 1))), 2);
    assert_eq!(count(&two, Some(q(141, 100)), Some(q(142, 100))), 1);
    assert_eq!(count(&two, Some(q(142, 100)), None), 0);
    assert_eq!(count(&two, None, Some(q(0, 1))), 1);
    assert_eq!(count(&two, Some(q(2, 1)), Some(q(1, 1))), 0);
    // Roots at the bounds count: x − 1 and x (x − 1).
    let one = poly(&[-1, 1]);
    assert_eq!(count(&one, Some(q(1, 1)), Some(q(2, 1))), 1);
    assert_eq!(count(&one, Some(q(0, 1)), Some(q(1, 1))), 1);
    assert_eq!(count(&one, Some(q(1, 1)), Some(q(1, 1))), 1);
    assert_eq!(count(&one, None, Some(q(1, 1))), 1);
    assert_eq!(count(&one, Some(q(1, 1)), None), 1);
    assert_eq!(count(&one, Some(q(3, 2)), Some(q(5, 2))), 0);
    let both = poly(&[0, -1, 1]);
    assert_eq!(count(&both, Some(q(0, 1)), Some(q(1, 1))), 2);
    assert_eq!(count(&both, Some(q(1, 2)), Some(q(1, 1))), 1);
    // (x − 1)³ (x + 2)²: one root in [0, 2], two on the line.
    let mut repeated = poly(&[1]);
    for factor in [
        poly(&[-1, 1]),
        poly(&[-1, 1]),
        poly(&[-1, 1]),
        poly(&[2, 1]),
        poly(&[2, 1]),
    ] {
        repeated = ring.mul(&repeated, &factor);
    }
    assert_eq!(count(&repeated, Some(q(0, 1)), Some(q(2, 1))), 1);
    assert_eq!(count(&repeated, None, None), 2);
    // A constant.
    assert_eq!(count(&poly(&[7]), None, None), 0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// For a product of `(x − rᵢ)^eᵢ`, times `x² + 1` or not, the count in
    /// `[lo, hi]` is the number of distinct `rᵢ` there.
    #[test]
    fn counts_known_roots(
        roots in prop::collection::vec(((-12i64..13, 1i64..5), 1u32..4), 0..6),
        positive in any::<bool>(),
        lo in prop::option::of((-15i64..16, 1i64..5)),
        hi in prop::option::of((-15i64..16, 1i64..5)),
    ) {
        let ring = qx();
        let mut p = poly(&[3]);
        let mut distinct: Vec<Rational> = Vec::new();
        for &((n, d), e) in &roots {
            let r = q(n, d);
            for _ in 0..e {
                p = ring.mul(&p, &ring.poly(vec![-&r, Rational::ONE]));
            }
            if !distinct.contains(&r) {
                distinct.push(r);
            }
        }
        if positive {
            p = ring.mul(&p, &poly(&[1, 0, 1]));
        }
        let (lo, hi) = (lo.map(|(n, d)| q(n, d)), hi.map(|(n, d)| q(n, d)));
        let expected = if lo.as_ref().zip(hi.as_ref()).is_some_and(|(lo, hi)| lo > hi) {
            0
        } else {
            distinct
                .iter()
                .filter(|r| lo.as_ref().is_none_or(|lo| *r >= lo) && hi.as_ref().is_none_or(|hi| *r <= hi))
                .count()
        };
        prop_assert_eq!(ring.count_real_roots(&p, lo.as_ref(), hi.as_ref()), expected, "{:?}", p);
    }
}
