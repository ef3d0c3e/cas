//! Rational reconstruction: exhaustively against a brute-force search for
//! small moduli, and on known fractions (built with GMP) for large ones.

use proptest::prelude::*;
use rug::Integer;

use super::{rational_reconstruction, rational_reconstruction_with_bounds};
use crate::testkit::{Mag, ToRug, strategy};
use crate::{Natural, Rational};

/// Every `n/d` with `n ≡ a d (mod m)`, `|n| <= nb`, `0 < d <= db` and
/// `gcd(d, m) = 1`, by trying each `d`.
fn brute_force(a: u64, m: u64, nb: u64, db: u64) -> Vec<Rational> {
    let mut found = Vec::new();
    for d in 1..=db {
        if Integer::from(d).gcd(&Integer::from(m)) != 1 {
            continue;
        }
        let n = u64::try_from(u128::from(a) * u128::from(d) % u128::from(m)).unwrap();
        let n = if n <= nb {
            i128::from(n)
        } else if m - n <= nb {
            -i128::from(m - n)
        } else {
            continue;
        };
        found.push(Rational::new(n.into(), d.into()));
    }
    found
}

/// The search finds at most one value (uniqueness when `2 nb db < m`), and
/// reconstruction agrees with it.
fn check(a: u64, m: u64, nb: u64, db: u64, ours: Option<&Rational>) {
    let found = brute_force(a, m, nb, db);
    assert!(
        found.windows(2).all(|w| w[0] == w[1]),
        "several fractions for a = {a}, m = {m}"
    );
    assert_eq!(ours, found.first(), "a = {a}, m = {m}, bounds {nb}, {db}");
}

#[test]
fn small_moduli_exhaustive() {
    for m in 1..=300u64 {
        let nm = Natural::from(m);
        let bound = (m - 1) / 2;
        let bound = (0..=bound).take_while(|b| 2 * b * b < m).last().unwrap();
        for a in 0..m {
            check(
                a,
                m,
                bound,
                bound,
                rational_reconstruction(&Natural::from(a), &nm).as_ref(),
            );
        }
    }
}

#[test]
fn small_moduli_with_bounds() {
    for m in 1..=120u64 {
        let nm = Natural::from(m);
        for nb in 0..m {
            for db in 0..m {
                if 2 * nb * db >= m {
                    break;
                }
                for a in 0..m {
                    let ours = rational_reconstruction_with_bounds(
                        &Natural::from(a),
                        &nm,
                        &Natural::from(nb),
                        &Natural::from(db),
                    );
                    check(a, m, nb, db, ours.as_ref());
                }
            }
        }
    }
}

/// `n d^(-1) mod m` with `m` above `2 max(|n|, d)^2` (or `2 |n| d` for the
/// bounded form), coprime to `d`.
fn known_fraction(n: &Integer, d: &Integer, symmetric: bool, pad: &Integer) -> (Integer, Integer) {
    let big = Integer::from(n.abs_ref()).max(d.clone());
    let floor = if symmetric {
        Integer::from(big.square_ref()) * 2
    } else {
        Integer::from(n.abs_ref()) * d * 2
    };
    let mut m: Integer = floor + pad + 1;
    while m.clone().gcd(d) != 1 {
        m += 1;
    }
    let a = (n.clone() * d.clone().invert(&m).unwrap()).modulo(&m);
    (a, m)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    #[test]
    fn recovers_known_fraction(
        n in strategy::signed_up_to(12),
        d in strategy::natural_up_to(12),
        pad in strategy::natural_up_to(4),
    ) {
        let (n, d) = (n.to_rug(), d.to_rug().max(Integer::from(1)));
        let expected = rug::Rational::from((n.clone(), d.clone()));
        let pad = pad.to_rug();

        let to_nat = |x: &Integer| Mag::from_rug(x).to_natural();

        let (a, m) = known_fraction(&n, &d, true, &pad);
        let ours = rational_reconstruction(&to_nat(&a), &to_nat(&m));
        prop_assert_eq!(ours.map(|q| q.to_rug()), Some(expected.clone()));

        let (a, m) = known_fraction(&n, &d, false, &pad);
        let ours = rational_reconstruction_with_bounds(
            &to_nat(&a),
            &to_nat(&m),
            &to_nat(&Integer::from(n.abs_ref())),
            &to_nat(&d),
        );
        prop_assert_eq!(ours.map(|q| q.to_rug()), Some(expected));
    }
}

#[test]
#[should_panic(expected = "zero modulus")]
fn zero_modulus_panics() {
    let _ = rational_reconstruction(&Natural::one(), &Natural::ZERO);
}

#[test]
#[should_panic(expected = "bounds too large")]
fn large_bounds_panic() {
    let five = Natural::from(5u64);
    let _ =
        rational_reconstruction_with_bounds(&Natural::one(), &Natural::from(50u64), &five, &five);
}
