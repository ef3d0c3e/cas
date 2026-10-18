//! Primality tests checked against GMP, against naive definitions, and on
//! known pseudoprimes.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, d, p, q, u, v, k, s)"
)]

use proptest::prelude::*;
use rug::Integer;
use rug::integer::IsPrime;

use super::strong_probable_prime as spp;
use super::{TestRing, U64_BASES, is_prime_u64, strong_lucas};
use crate::modular::{ModDouble, ModRing};
use crate::testkit::{Mag, ToRug, strategy};
use crate::{Natural, Primality};

fn gmp_is_prime(n: &Integer) -> bool {
    n.is_probably_prime(40) != IsPrime::No
}

/// The label we must give: proven below 2^64, probable above.
fn expected(n: &Integer) -> Primality {
    match (gmp_is_prime(n), n.significant_bits() <= 64) {
        (false, _) => Primality::NotPrime,
        (true, true) => Primality::Prime,
        (true, false) => Primality::ProbablePrime,
    }
}

/// Strong probable-prime test written directly with GMP.
fn gmp_strong_probable_prime(n: &Integer, base: &Integer) -> bool {
    if *n < 2 || n.is_even() {
        return *n == 2;
    }
    let n1 = Integer::from(n - 1);
    let s = n1.find_one(0).unwrap();
    if base.is_divisible(n) {
        return true;
    }
    let mut y = base.clone().pow_mod(&(n1.clone() >> s), n).unwrap();
    if y == 1 || y == n1 {
        return true;
    }
    for _ in 1..s {
        y = y.square().modulo(n);
        if y == n1 {
            return true;
        }
    }
    false
}

/// Primes, composites with large factors and plain random values.
fn prime_input() -> BoxedStrategy<Mag> {
    let next = |m: Mag| Mag::from_rug(&m.to_rug().next_prime());
    prop_oneof![
        2 => strategy::natural_up_to(6),
        3 => strategy::natural_up_to(6).prop_map(next),
        1 => strategy::natural_up_to(1).prop_map(next),
        2 => (strategy::natural_up_to(3), strategy::natural_up_to(3))
            .prop_map(move |(a, b)| Mag::from_rug(&(a.to_rug().next_prime() * b.to_rug().next_prime()))),
        1 => strategy::natural_up_to(3).prop_map(|a| Mag::from_rug(&a.to_rug().next_prime().square())),
    ]
    .boxed()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn primality(n in prime_input()) {
        prop_assert_eq!(n.to_natural().primality(), expected(&n.to_rug()));
    }

    #[test]
    fn strong_probable_prime(n in strategy::natural_up_to(4), base in strategy::natural_up_to(5)) {
        prop_assert_eq!(
            n.to_natural().is_strong_probable_prime(&base.to_natural()),
            gmp_strong_probable_prime(&n.to_rug(), &base.to_rug())
        );
    }

    #[test]
    fn next_prime(n in strategy::natural_up_to(3)) {
        prop_assert_eq!(n.to_natural().next_prime().to_rug(), n.to_rug().next_prime());
    }
}

#[test]
fn small_values_exhaustive() {
    for n in 0..200_000u64 {
        assert_eq!(
            Natural::from(n).primality(),
            expected(&Integer::from(n)),
            "{n}"
        );
    }
}

#[test]
fn word_boundaries() {
    let ranges = [
        (1 << 16) - 100..(1 << 16) + 100,
        (1 << 32) - 200..(1 << 32) + 200,
        u64::MAX - 20_000..u64::MAX,
    ];
    for n in ranges.into_iter().flatten().chain([u64::MAX]) {
        assert_eq!(is_prime_u64(n), gmp_is_prime(&Integer::from(n)), "{n}");
    }
    // Just above 2^64: the BPSW path.
    for k in 0..400u64 {
        let n = Natural::power_of_two(64) + k;
        assert_eq!(n.primality(), expected(&n.to_rug()), "2^64 + {k}");
    }
}

#[test]
fn next_prime_edges() {
    assert_eq!(Natural::ZERO.next_prime(), Natural::from(2u64));
    assert_eq!(Natural::from(2u64).next_prime(), Natural::from(3u64));
    // The last prime below 2^64 is 2^64 - 59; the first above is 2^64 + 13.
    let last = Natural::from(u64::MAX - 58);
    assert_eq!((&last - 1u64).next_prime(), last);
    assert_eq!(last.next_prime(), Natural::power_of_two(64) + 13u64);
    assert_eq!(
        Natural::from(u64::MAX).next_prime(),
        Natural::power_of_two(64) + 13u64
    );
}

/// `(n, bases)`: strong pseudoprimes to every listed base.
const STRONG_PSEUDOPRIMES: &[(&str, &[u64])] = &[
    ("2047", &[2]),
    ("3277", &[2]),
    ("4033", &[2]),
    ("4681", &[2]),
    ("8321", &[2]),
    ("1373653", &[2, 3]),
    ("25326001", &[2, 3, 5]),
    ("3215031751", &[2, 3, 5, 7]),
    ("2152302898747", &[2, 3, 5, 7, 11]),
    ("3474749660383", &[2, 3, 5, 7, 11, 13]),
    ("341550071728321", &[2, 3, 5, 7, 11, 13, 17]),
    ("3825123056546413051", &[2, 3, 5, 7, 11, 13, 17, 19, 23]),
    (
        "318665857834031151167461",
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37],
    ),
    (
        "3317044064679887385961981",
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41],
    ),
];

#[test]
fn strong_pseudoprimes_are_caught() {
    for &(n, bases) in STRONG_PSEUDOPRIMES {
        let n: Natural = n.parse().unwrap();
        for &b in bases {
            assert!(
                n.is_strong_probable_prime(&Natural::from(b)),
                "{n} base {b}"
            );
        }
        assert!(!gmp_is_prime(&n.to_rug()), "{n} is composite");
        assert_eq!(n.primality(), Primality::NotPrime, "{n}");
    }
}

/// Carmichael numbers fool the Fermat test to every coprime base.
#[test]
fn carmichael_numbers_are_caught() {
    for n in [
        561u64,
        1105,
        1729,
        2465,
        2821,
        6601,
        8911,
        41041,
        825_265,
        321_197_185,
        5_394_826_801,
        232_250_619_601,
        9_746_347_772_161,
    ] {
        assert_eq!(Natural::from(n).primality(), Primality::NotPrime, "{n}");
    }
}

/// `U_k` and `V_k` by the recurrences `W_(k+1) = P W_k - Q W_(k-1)`.
fn naive_strong_lucas(n: i128) -> bool {
    let jacobi = |d: i128| Integer::from(d).jacobi(&Integer::from(n));
    let mut d = 5i128;
    loop {
        match jacobi(d) {
            -1 => break,
            0 => return false,
            _ => d = if d > 0 { -d - 2 } else { -d + 2 },
        }
    }
    let (p, q) = (1i128, (1 - d) / 4);
    let n1 = n + 1;
    let s = n1.trailing_zeros();
    let k = n1 >> s;
    let (mut u, mut u_prev) = (1i128, 0i128);
    let (mut v, mut v_prev) = (p, 2i128);
    let mut us = vec![0, 1];
    let mut vs = vec![2, p];
    for _ in 2..=n1 {
        (u, u_prev) = ((p * u - q * u_prev).rem_euclid(n), u);
        (v, v_prev) = ((p * v - q * v_prev).rem_euclid(n), v);
        us.push(u);
        vs.push(v);
    }
    let idx = |j: i128| usize::try_from(j).unwrap();
    us[idx(k)] == 0 || (0..s).any(|r| vs[idx(k << r)] == 0)
}

#[test]
fn strong_lucas_matches_definition() {
    for n in (1001..6001u64).step_by(2) {
        let nat = Natural::from(n);
        if nat.is_perfect_square() {
            continue;
        }
        let ours = strong_lucas(&ModRing::new(nat), &Natural::from(n));
        assert_eq!(ours, naive_strong_lucas(i128::from(n)), "{n}");
    }
}

/// Strong Lucas pseudoprimes with Selfridge's parameters (OEIS A217255),
/// which pass the Lucas half of BPSW but are composite.
#[test]
fn strong_lucas_pseudoprimes() {
    for n in [
        5459u64, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519,
    ] {
        let nat = Natural::from(n);
        assert!(strong_lucas(&ModRing::new(nat.clone()), &nat), "{n}");
        assert!(!gmp_is_prime(&Integer::from(n)), "{n} is composite");
        assert_eq!(nat.primality(), Primality::NotPrime);
    }
}

#[test]
fn known_primes() {
    for p in [61u64, 89, 107, 127, 521, 607, 1279] {
        let m = Natural::power_of_two(p) - 1u64;
        let label = if p <= 64 {
            Primality::Prime
        } else {
            Primality::ProbablePrime
        };
        assert_eq!(m.primality(), label, "2^{p} - 1");
    }
    for p in [67u64, 101, 257] {
        let m = Natural::power_of_two(p) - 1u64;
        assert_eq!(m.primality(), Primality::NotPrime, "2^{p} - 1");
    }
}

#[test]
fn u64_bases_are_sinclairs() {
    assert_eq!(
        U64_BASES,
        [2, 325, 9375, 28178, 450_775, 9_780_504, 1_795_265_022]
    );
}

/// A prime dividing one of the bases sees that base as 0 and must skip it
/// (for example 299210837, a factor of 1795265022).
#[test]
fn primes_dividing_a_base() {
    let mut checked = 0;
    for b in U64_BASES {
        let mut m = b;
        let mut p = 2;
        while m > 1 {
            if m % p == 0 {
                assert!(is_prime_u64(p), "{p} divides base {b}");
                while m % p == 0 {
                    m /= p;
                }
                checked += u32::from(p >= 1 << 16);
            }
            p += 1;
            if p * p > m && m > 1 {
                assert!(is_prime_u64(m), "{m} divides base {b}");
                checked += u32::from(m >= 1 << 16);
                break;
            }
        }
    }
    assert!(checked >= 1, "no base has a prime factor past the table");
}

/// Odd two-limb values: primes, products of two primes near `2^64`, near
/// squares and plain random ones.
fn two_limb_odd() -> BoxedStrategy<Natural> {
    let two_limbs = (any::<u64>(), 1..u64::MAX)
        .prop_map(|(lo, hi)| Integer::from(lo) | (Integer::from(hi) << 64));
    prop_oneof![
        two_limbs.clone(),
        two_limbs.clone().prop_map(|n: Integer| n.next_prime()),
        (any::<u32>(), any::<u32>()).prop_map(|(a, b)| {
            let p = Integer::from((1u64 << 32) + u64::from(a)).next_prime();
            let q = Integer::from((1u64 << 40) + u64::from(b)).next_prime();
            p * q
        }),
        two_limbs.prop_map(|n: Integer| {
            let r = n.sqrt();
            r.clone() * (r + 2)
        }),
    ]
    .prop_map(|n| Mag::from_rug(&(n | 1)).to_natural())
    .prop_filter("two limbs", |n| n.limb_len() == 2)
    .boxed()
}

proptest! {
    /// The `u128` ring gives the same Miller–Rabin and Lucas answers as
    /// `ModRing`, base by base.
    #[test]
    fn two_limb_ring_matches(n in two_limb_odd(), base in 2..1000u64) {
        let [lo, hi] = *n.limbs() else { unreachable!() };
        let double = ModDouble::new(u128::from(lo) | u128::from(hi) << 64);
        let ring = ModRing::new(n.clone());
        let b = i64::try_from(base).unwrap();
        prop_assert_eq!(
            spp(&double, &n, &double.small(b)),
            spp(&ring, &n, &ring.small(b))
        );
        if !n.is_perfect_square() {
            prop_assert_eq!(strong_lucas(&double, &n), strong_lucas(&ring, &n));
        }
        prop_assert_eq!(n.primality(), expected(&n.to_rug()));
    }
}
