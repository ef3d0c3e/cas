//! The elliptic curve method on products of known primes (found with GMP).

use proptest::prelude::*;
use rug::Integer;

use crate::Natural;
use crate::testkit::{Mag, ToRug};

fn natural(x: &Integer) -> Natural {
    Mag::from_rug(x).to_natural()
}

/// `ecm`'s answer, checked to be a proper divisor.
fn proper(n: &Natural, b1: u64, curves: u32, seed: u64) -> Option<Natural> {
    let d = n.ecm(b1, curves, seed)?;
    let (n, d) = (n.to_rug(), d.to_rug());
    assert!(d > 1 && d < n && n.is_divisible(&d), "{d} of {n}");
    Some(natural(&d))
}

#[test]
fn nothing_to_split() {
    for x in [0u64, 1, 2, 3, 5, 17, 65_537, (1 << 61) - 1] {
        assert_eq!(Natural::from(x).ecm(2000, 5, 0), None, "{x}");
    }
    let p = natural(&(Integer::from(1) << 200u32).next_prime());
    assert_eq!(p.ecm(2000, 5, 0), None);
    assert_eq!(
        Natural::from(360u64).ecm(2000, 5, 0),
        Some(Natural::from(2u64))
    );
    assert_eq!(
        Natural::from(4u64).ecm(2000, 5, 0),
        Some(Natural::from(2u64))
    );
}

/// 2^128 + 1 = 59649589127497217 · 5704689200685129054721, the 17-digit
/// factor some 2^28 rho steps away.
#[test]
fn fermat_seven() {
    let n = natural(&((Integer::from(1) << 128u32) + 1));
    let d = proper(&n, 11_000, 50, 0).expect("a factor");
    assert_eq!(d.to_rug(), 59_649_589_127_497_217u64);
    let f = n.factor_with_ecm(1000, 11_000, 50);
    assert!(f.is_complete());
    assert_eq!(f.factors().len(), 2);
    assert_eq!(f.value(), n);
    // Without curves, the same rho effort leaves it whole.
    assert!(!n.factor_with_effort(1000).is_complete());
}

/// Factors of 16 and 18 digits: B1 = 2000 needs stage 2 to find the
/// smaller within 40 curves, and the factorisation completes.
#[test]
fn two_stages() {
    let p = Integer::from(10u64.pow(15) + 37);
    let q = Integer::from(10u64.pow(17) + 3);
    let n = natural(&(p.clone() * &q));
    assert_eq!(proper(&n, 2000, 40, 0).map(|d| d.to_rug()), Some(p.clone()));
    let r = Integer::from(10u128.pow(20) + 39);
    let m = natural(&(p * q * r));
    let f = m.factor_with_ecm(1 << 12, 11_000, 100);
    assert!(f.is_complete());
    assert_eq!(f.factors().len(), 3);
    assert_eq!(f.value(), m);
}

/// One-limb moduli, where both factors are near 2^31.
#[test]
fn one_limb() {
    let p = Integer::from(1u64 << 31).next_prime();
    let q = Integer::from(3u64 << 30).next_prime();
    let n = natural(&(p * q));
    assert!(proper(&n, 500, 50, 0).is_some());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Whatever a curve finds is a proper divisor, for semiprimes and for
    /// prime powers times a prime.
    #[test]
    fn divisors(a in 20..44u32, b in 20..64u32, x in any::<u64>(), y in any::<u64>(),
                square in any::<bool>(), seed in 0..1000u64) {
        let p = (Integer::from(x >> (64 - a)) | (Integer::from(1) << (a - 1))).next_prime();
        let q = (Integer::from(y >> (64 - b)) | (Integer::from(1) << (b - 1))).next_prime();
        let p = if square { p.clone() * &p } else { p };
        let n = natural(&(p * q));
        let _ = proper(&n, 300, 8, seed);
    }
}
