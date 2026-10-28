//! Tests of the Gauss–Legendre rules: the small rules exactly, the sums
//! of the weights, exactness for odd and even powers up to `x^(2n−1)`,
//! order and symmetry, and the cache.

use super::{gauss_legendre, legendre};
use crate::ball::Ball;
use crate::{Integer, Rational};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// Whether `b` contains `x`, or the ball of `x` overlaps it.
fn holds(b: &Ball, x: &Ball) -> bool {
    b.overlaps(x)
}

#[test]
fn small_rules() {
    let prec = 128;
    let one = gauss_legendre(1, prec);
    assert_eq!(one.len(), 1);
    assert!(one[0].0.contains_rational(&q(0, 1)));
    assert!(one[0].1.contains_rational(&q(2, 1)));
    // ±1/√3, weights 1.
    let two = gauss_legendre(2, prec);
    let third = Ball::from_rational(&q(1, 3), 4 * prec).sqrt(4 * prec);
    assert!(holds(&two[0].0, &third) && holds(&two[1].0, &third.neg()));
    assert!(two.iter().all(|(_, w)| w.contains_rational(&q(1, 1))));
    // 0, ±√(3/5); 8/9 and 5/9.
    let three = gauss_legendre(3, prec);
    let root = Ball::from_rational(&q(3, 5), 4 * prec).sqrt(4 * prec);
    assert!(holds(&three[0].0, &root) && holds(&three[2].0, &root.neg()));
    assert!(three[1].0.contains_rational(&q(0, 1)));
    assert!(three[1].1.contains_rational(&q(8, 9)));
    assert!(three[0].1.contains_rational(&q(5, 9)) && three[2].1.contains_rational(&q(5, 9)));
    // Tight: radii near 2^-prec.
    for (x, w) in three.iter() {
        assert!(x.rel_accuracy_bits() >= i64::try_from(prec).unwrap() - 8 || x.is_exact());
        assert!(w.rel_accuracy_bits() >= i64::try_from(prec).unwrap() - 8);
    }
}

/// `Σ w_k x_k^j = ∫_{−1}^{1} x^j = 2/(j + 1)` for even `j`, 0 for odd `j`,
/// for `j ≤ 2n − 1`; the nodes descend in `(−1, 1)`, apart, symmetric.
#[test]
fn exactness_and_order() {
    for (n, prec) in [
        (4, 64),
        (7, 100),
        (16, 200),
        (33, 128),
        (64, 256),
        (100, 160),
    ] {
        let rule = gauss_legendre(n, prec);
        assert_eq!(rule.len(), n);
        let wp = prec + 64;
        for j in 0..2 * n {
            let mut sum = Ball::ZERO;
            for (x, w) in rule.iter() {
                let mut power = Ball::from_integer(&Integer::from(1));
                for _ in 0..j {
                    power = power.mul(x, wp);
                }
                sum = sum.add(&w.mul(&power, wp), wp);
            }
            let exact = if j % 2 == 1 {
                q(0, 1)
            } else {
                q(2, i64::try_from(j).unwrap() + 1)
            };
            assert!(sum.contains_rational(&exact), "n = {n}, x^{j}: {sum}");
        }
        for pair in rule.windows(2) {
            let (a, b) = (&pair[0].0, &pair[1].0);
            assert_eq!(
                b.sub(a, wp).sign(),
                Some(std::cmp::Ordering::Less),
                "n = {n}"
            );
        }
        for (k, (x, w)) in rule.iter().enumerate() {
            let (mirror, v) = &rule[n - 1 - k];
            assert!(x.overlaps(&mirror.neg()) && w.overlaps(v), "n = {n}");
            let (p, _) = legendre(n, x, wp);
            assert!(p.contains_rational(&q(0, 1)), "n = {n}: P_n(x_{k}) = {p}");
        }
    }
}

#[test]
fn the_cache() {
    let a = gauss_legendre(12, 96);
    let b = gauss_legendre(12, 96);
    assert!(std::sync::Arc::ptr_eq(&a, &b));
    super::clear_rules();
    let c = gauss_legendre(12, 96);
    assert!(!std::sync::Arc::ptr_eq(&a, &c));
    assert_eq!(a.len(), c.len());
}
