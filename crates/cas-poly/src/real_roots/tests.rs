//! Tests of the real root isolation (CAS-538): giac's LMQ bounds and VAS
//! intervals on its own examples, and that each interval VAS gives holds
//! exactly one root, by Sturm's count.

use cas_num::algebra::{Rationals, Ring};
use cas_num::{Integer, Rational};

use super::{Isolated, LowerBound, poslbd_lmq, posub_lmq, vas, vas_positive};
use crate::dense::{Dense, DenseRing};

fn qx() -> DenseRing<Rationals> {
    DenseRing::new(Rationals)
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// A polynomial from its coefficients, the constant first.
fn poly(coeffs: &[i64]) -> Dense<Rational> {
    qx().from_i64s(coeffs)
}

fn root(n: i64, d: i64) -> Isolated {
    Isolated::Root(q(n, d))
}

fn between(a: (i64, i64), b: (i64, i64)) -> Isolated {
    Isolated::Interval(q(a.0, a.1), q(b.0, b.1))
}

/// giac's `posubLMQ` answers.
#[test]
fn upper_bounds() {
    for (p, bound) in [
        (poly(&[7, -7, 0, 1]), 4),
        (poly(&[-7, 7, 0, -1]), 4),
        (poly(&[-2, 0, 1]), 3),
        (poly(&[1, 0, 1]), 0),
        (poly(&[1, 0, -7, 7]), 3),
        (poly(&[-1, 0, 2]), 2),
        (poly(&[-4, 0, 1]), 3),
        (poly(&[-3, 0, 1]), 3),
        (poly(&[-5, 0, 1]), 4),
        (poly(&[-2, 1]), 5),
        (poly(&[-5, 2]), 6),
        (poly(&[-8, 0, 0, 1]), 3),
        (poly(&[-100, 0, 1]), 15),
        (poly(&[1, 0, -10, 0, 1]), 5),
        (poly(&[2, -3, 1]), 7),
        (poly(&[-1, 1, -3, 0, 0, 1]), 2),
        (poly(&[-7, 1, 0, -2, 3, 1]), 2),
        (poly(&[0, -1, 1]), 3),
        (poly(&[0, 0, 1]), 0),
        (poly(&[5]), 0),
        (qx().poly(vec![q(-3, 1), q(1, 2)]), 13),
    ] {
        assert_eq!(posub_lmq(&p), Integer::from(bound), "{p:?}");
    }
}

/// giac's `poslbdLMQ` answers.
#[test]
fn lower_bounds() {
    for (p, bound) in [
        (poly(&[7, -7, 0, 1]), LowerBound::Power(-1)),
        (poly(&[-2, 0, 1]), LowerBound::Power(0)),
        (poly(&[1, 1, 0, 1]), LowerBound::Infinity),
        (poly(&[1, 0, -10, 0, 1]), LowerBound::Power(-3)),
        (poly(&[-2, 1]), LowerBound::Power(0)),
        (poly(&[-3, 1]), LowerBound::Power(0)),
        (poly(&[-100, 1]), LowerBound::Power(5)),
        (poly(&[2, -3, 1]), LowerBound::Power(-2)),
        (poly(&[-1, 1, -3, 0, 0, 1]), LowerBound::Power(-1)),
        (poly(&[-7, 1, 0, -2, 3, 1]), LowerBound::Power(-1)),
        (poly(&[5, -1]), LowerBound::Power(1)),
        (poly(&[0, 0, 1]), LowerBound::Infinity),
        (poly(&[0, -1, 1]), LowerBound::Zero),
        (poly(&[5]), LowerBound::Zero),
        (
            qx().poly(vec![q(1, 5), q(-7, 1), q(1, 3)]),
            LowerBound::Power(-7),
        ),
    ] {
        assert_eq!(poslbd_lmq(&p), bound, "{p:?}");
    }
    assert_eq!(LowerBound::Power(-1).value(), Some(q(1, 2)));
    assert_eq!(LowerBound::Infinity.value(), None);
}

/// giac's `VAS` and `VAS_positive` answers.
#[test]
fn isolation() {
    let cubic = poly(&[7, -7, 0, 1]);
    assert_eq!(
        vas(&cubic),
        [
            between((-4, 1), (0, 1)),
            between((1, 1), (3, 2)),
            between((3, 2), (2, 1)),
        ]
    );
    assert_eq!(
        vas(&poly(&[-7, 7, 0, -1])),
        [
            between((-4, 1), (0, 1)),
            between((1, 1), (3, 2)),
            between((3, 2), (2, 1)),
        ]
    );
    assert_eq!(
        vas_positive(&cubic),
        [between((1, 1), (3, 2)), between((3, 2), (2, 1))]
    );
    for (p, expected) in [
        (poly(&[2, -3, 1]), vec![root(1, 1), between((1, 1), (4, 1))]),
        (
            poly(&[0, 2, -3, 1]),
            vec![root(1, 1), between((1, 1), (3, 1))],
        ),
        (poly(&[0, -1, 1]), vec![between((0, 1), (3, 1))]),
        (poly(&[1, 0, 1]), vec![]),
    ] {
        assert_eq!(vas_positive(&p), expected, "{p:?}");
    }
    for p in [poly(&[1, 0, 1]), poly(&[5]), poly(&[0])] {
        assert!(vas(&p).is_empty(), "{p:?}");
    }
}

/// giac's `VAS` answers on more polynomials: exact roots met by the
/// search, the root 0, even polynomials mirrored, large roots scaled.
#[test]
fn isolation_examples() {
    for (p, expected) in [
        (
            poly(&[-2, 0, 0, 0, 1]),
            vec![between((-2, 1), (0, 1)), between((0, 1), (2, 1))],
        ),
        (poly(&[0, 0, 1]), vec![root(0, 1)]),
        (poly(&[0, 0, 0, 0, 1]), vec![root(0, 1)]),
        // (x − 1)²(x − 2).
        (
            poly(&[-2, 5, -4, 1]),
            vec![root(1, 1), between((1, 1), (4, 1))],
        ),
        (poly(&[-1, 1, -3, 0, 0, 1]), vec![between((1, 1), (2, 1))]),
        (
            poly(&[1, 0, -10, 0, 1]),
            vec![
                between((-6, 1), (-1, 1)),
                between((-1, 1), (0, 1)),
                between((0, 1), (1, 1)),
                between((1, 1), (6, 1)),
            ],
        ),
        (
            qx().poly(vec![q(-1, 1), q(0, 1), q(1, 2)]),
            vec![between((-3, 1), (0, 1)), between((0, 1), (3, 1))],
        ),
        (poly(&[-100, 1]), vec![between((0, 1), (204, 1))]),
        (
            poly(&[20000, -300, 1]),
            vec![between((96, 1), (128, 1)), between((128, 1), (224, 1))],
        ),
        (poly(&[-1, 3]), vec![between((0, 1), (1, 1))]),
        (
            poly(&[5, -16, 3]),
            vec![between((0, 1), (1, 1)), between((1, 1), (8, 1))],
        ),
        (
            poly(&[0, -1, 0, 1]),
            vec![
                between((-2, 1), (0, 1)),
                root(0, 1),
                between((0, 1), (2, 1)),
            ],
        ),
        // (x² − 2)(x − 1)².
        (
            poly(&[-2, 4, -1, -2, 1]),
            vec![
                between((-3, 1), (0, 1)),
                root(1, 1),
                between((1, 1), (3, 1)),
            ],
        ),
        (
            poly(&[-1, 0, 100]),
            vec![between((-1, 1), (0, 1)), between((0, 1), (1, 1))],
        ),
        (
            poly(&[0, -1, 0, 0, 0, 1]),
            vec![
                between((-2, 1), (0, 1)),
                root(0, 1),
                between((0, 1), (2, 1)),
            ],
        ),
        (
            poly(&[4, 0, -5, 0, 1]),
            vec![
                between((-3, 1), (-1, 1)),
                root(-1, 1),
                root(1, 1),
                between((1, 1), (3, 1)),
            ],
        ),
        (
            poly(&[0, 2, -3, 1]),
            vec![root(0, 1), root(1, 1), between((1, 1), (4, 1))],
        ),
        (
            poly(&[-36, 0, 49, 0, -14, 0, 1]),
            vec![
                between((-4, 1), (-2, 1)),
                root(-2, 1),
                root(-1, 1),
                root(1, 1),
                root(2, 1),
                between((2, 1), (4, 1)),
            ],
        ),
        (
            poly(&[6, -11, -3, 2]),
            vec![
                between((-4, 1), (0, 1)),
                between((0, 1), (1, 1)),
                between((1, 1), (5, 1)),
            ],
        ),
    ] {
        assert_eq!(vas(&p), expected, "{p:?}");
    }
}

/// Each interval of `vas(p)` holds exactly one root in its interior, each
/// exact entry is a root, and the entries count the distinct real roots.
#[test]
fn intervals_hold_one_root() {
    let ring = qx();
    let wilkinson = (1..=8).fold(poly(&[1]), |p, k| ring.mul(&p, &poly(&[-k, 1])));
    let mixed = ring.mul(
        &ring.mul(&poly(&[-2, 0, 1]), &poly(&[-1, 3])),
        &poly(&[1, 1, 1]),
    );
    for p in [
        poly(&[7, -7, 0, 1]),
        poly(&[1, 0, -10, 0, 1]),
        poly(&[-2, 5, -4, 1]),
        poly(&[-1, 1, -3, 0, 0, 1]),
        poly(&[20000, -300, 1]),
        poly(&[-36, 0, 49, 0, -14, 0, 1]),
        // Chebyshev's T₈: eight roots in (−1, 1).
        poly(&[1, 0, -32, 0, 160, 0, -256, 0, 128]),
        poly(&[-2, 0, 0, 0, 0, 0, 0, 1]),
        // Mignotte's x¹⁰ − 2(10x − 1)²: two roots close to 1/10.
        poly(&[-2, 40, -200, 0, 0, 0, 0, 0, 0, 0, 1]),
        wilkinson,
        mixed,
    ] {
        let found = vas(&p);
        assert_eq!(found.len(), ring.count_real_roots(&p, None, None), "{p:?}");
        for entry in &found {
            match entry {
                Isolated::Root(r) => assert!(ring.eval(&p, r).is_zero(), "{p:?}: {r:?}"),
                Isolated::Interval(lo, hi) => {
                    let at_ends = [lo, hi]
                        .into_iter()
                        .filter(|x| ring.eval(&p, x).is_zero())
                        .count();
                    let inside = ring.count_real_roots(&p, Some(lo), Some(hi)) - at_ends;
                    assert!(lo < hi && inside == 1, "{p:?}: [{lo:?}, {hi:?}]");
                }
            }
        }
    }
}
