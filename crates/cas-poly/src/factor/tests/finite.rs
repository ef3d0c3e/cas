//! Tests of factoring over word-size prime fields (CAS-180).

#![allow(clippy::many_single_char_names, reason = "the names of the algorithms")]

use std::sync::LazyLock;
use std::sync::atomic::AtomicBool;

use cas_num::algebra::{EuclideanDomain, Ring, WordPrimeField};
use proptest::prelude::*;

use crate::cancel::{Cancelled, Never};
use crate::dense::{Dense, DenseRing};
use crate::factor::finite::{Random, distinct_degree, equal_degree, factor_squarefree};

fn fp(p: u64) -> DenseRing<WordPrimeField> {
    DenseRing::new(WordPrimeField::new(p).unwrap())
}

/// The monic polynomials of degree `d` over 𝔽ₚ, all `p^d` of them.
fn monics(p: u64, d: usize) -> Vec<Vec<u64>> {
    let count = p.pow(u32::try_from(d).unwrap());
    (0..count)
        .map(|mut k| {
            let mut c: Vec<u64> = (0..d)
                .map(|_| {
                    let digit = k % p;
                    k /= p;
                    digit
                })
                .collect();
            c.push(1);
            c
        })
        .collect()
}

/// The monic irreducible polynomials of degrees 1 to 4 over 𝔽ₚ, by trial
/// division by every monic polynomial of at most half their degree.
fn irreducibles(p: u64) -> Vec<Dense<u64>> {
    let ring = fp(p);
    let mut out = Vec::new();
    for d in 1..=4 {
        let divisors: Vec<Dense<u64>> = (1..=d / 2)
            .flat_map(|e| monics(p, e))
            .map(|c| ring.poly(c))
            .collect();
        for c in monics(p, d) {
            let f = ring.poly(c);
            if divisors.iter().all(|g| !ring.div_rem(&f, g).1.is_zero()) {
                out.push(f);
            }
        }
    }
    out
}

static IRREDUCIBLES: LazyLock<Vec<(u64, Vec<Dense<u64>>)>> = LazyLock::new(|| {
    [3, 5, 7]
        .into_iter()
        .map(|p| (p, irreducibles(p)))
        .collect()
});

/// The order of [`factor_squarefree`]'s results: by degree, then by
/// coefficients from the constant term up.
fn sorted(mut v: Vec<Dense<u64>>) -> Vec<Dense<u64>> {
    v.sort_by(|a, b| {
        a.degree()
            .cmp(&b.degree())
            .then_with(|| a.coeffs().cmp(b.coeffs()))
    });
    v
}

fn product(ring: &DenseRing<WordPrimeField>, fs: &[Dense<u64>]) -> Dense<u64> {
    fs.iter().fold(ring.one(), |acc, f| ring.mul(&acc, f))
}

/// `x^(p^k) mod f` by `k` powers by `p`, each by square-and-multiply on the
/// bits of `p`.
fn frobenius(ring: &DenseRing<WordPrimeField>, f: &Dense<u64>, k: usize) -> Dense<u64> {
    let p = u64::try_from(&ring.characteristic()).unwrap();
    let mut h = ring.div_rem(&ring.x(), f).1;
    for _ in 0..k {
        let mut r = ring.one();
        for i in (0..u64::BITS - p.leading_zeros()).rev() {
            r = ring.div_rem(&ring.mul(&r, &r), f).1;
            if (p >> i) & 1 == 1 {
                r = ring.div_rem(&ring.mul(&r, &h), f).1;
            }
        }
        h = r;
    }
    h
}

/// Rabin's test: monic `f` of degree `n` is irreducible over 𝔽ₚ exactly
/// when `x^(p^n) ≡ x (mod f)` and `gcd(x^(p^(n/r)) − x, f) = 1` for each
/// prime `r` dividing `n`.
fn rabin(ring: &DenseRing<WordPrimeField>, f: &Dense<u64>) -> bool {
    let n = f.degree().unwrap();
    let x = ring.div_rem(&ring.x(), f).1;
    if frobenius(ring, f, n) != x {
        return false;
    }
    (2..=n)
        .filter(|&r| n % r == 0 && (2..r).all(|s| r % s != 0))
        .all(|r| {
            let h = ring.sub(&frobenius(ring, f, n / r), &x);
            EuclideanDomain::gcd(ring, &h, f).degree() == Some(0)
        })
}

#[test]
fn splitmix() {
    let mut random = Random::new(0);
    assert_eq!(random.next(), 0xE220_A839_7B1D_CDAF);
    assert_eq!(random.next(), 0x6E78_9E6A_A1B9_65F4);
    assert_eq!(random.next(), 0x06C4_5D18_8009_454F);
}

/// The counts of Gauss's formula, `(1/d) Σ_(e | d) μ(e) p^(d/e)`, check the
/// brute-force reference.
#[test]
fn counts_of_irreducibles() {
    let expected = [
        (3, [3, 3, 8, 18]),
        (5, [5, 10, 40, 150]),
        (7, [7, 21, 112, 588]),
    ];
    for ((p, list), (q, counts)) in IRREDUCIBLES.iter().zip(expected) {
        assert_eq!(*p, q);
        for (d, count) in (1..=4).zip(counts) {
            assert_eq!(
                list.iter().filter(|f| f.degree() == Some(d)).count(),
                count,
                "degree {d} over F_{p}"
            );
        }
    }
}

#[test]
fn small_cases() {
    let ring = fp(5);
    assert!(
        factor_squarefree(&ring, &ring.one(), &Never)
            .unwrap()
            .is_empty()
    );
    assert!(
        distinct_degree(&ring, &ring.one(), &Never)
            .unwrap()
            .is_empty()
    );
    let x = ring.x();
    assert_eq!(
        factor_squarefree(&ring, &x, &Never).unwrap(),
        vec![x.clone()]
    );
    // x^p - x is the product of the x - a.
    for p in [3, 5, 7, 11] {
        let ring = fp(p);
        let f = ring.sub(&ring.pow(&ring.x(), p), &ring.x());
        let expected: Vec<Dense<u64>> = (0..p).map(|a| ring.poly(vec![a, 1])).collect();
        assert_eq!(factor_squarefree(&ring, &f, &Never).unwrap(), expected);
    }
    // x^9 - x over F_3: the three linear and three quadratic irreducibles.
    let ring = fp(3);
    let f = ring.sub(&ring.pow(&ring.x(), 9), &ring.x());
    let ddf = distinct_degree(&ring, &f, &Never).unwrap();
    let linear = ring.sub(&ring.pow(&ring.x(), 3), &ring.x());
    let quadratic = ring.div_rem(&f, &linear).0;
    assert_eq!(ddf, vec![(linear, 1), (quadratic, 2)]);
    let expected: Vec<Dense<u64>> = IRREDUCIBLES[0]
        .1
        .iter()
        .filter(|g| g.degree() <= Some(2))
        .cloned()
        .collect();
    assert_eq!(
        factor_squarefree(&ring, &f, &Never).unwrap(),
        sorted(expected)
    );
}

#[test]
fn cancellation() {
    let ring = fp(7);
    let f = ring.sub(&ring.pow(&ring.x(), 7), &ring.x());
    let stop = AtomicBool::new(true);
    assert_eq!(factor_squarefree(&ring, &f, &stop), Err(Cancelled));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Products of distinct irreducibles of degrees 1 to 4 over 𝔽₃, 𝔽₅ and
    /// 𝔽₇ split back into them, with their products by degree from the
    /// distinct-degree step.
    #[test]
    fn products_of_irreducibles_split_back((field, picks) in (0usize..3).prop_flat_map(|i| (Just(i), prop::sample::subsequence((0..IRREDUCIBLES[i].1.len()).collect::<Vec<_>>(), 1..=6)))) {
        let (p, list) = &IRREDUCIBLES[field];
        let ring = fp(*p);
        let factors: Vec<Dense<u64>> = picks.iter().map(|&k| list[k].clone()).collect();
        let f = product(&ring, &factors);
        prop_assert_eq!(factor_squarefree(&ring, &f, &Never).unwrap(), sorted(factors.clone()));
        let mut degrees: Vec<usize> = factors.iter().map(|g| g.degree().unwrap()).collect();
        degrees.sort_unstable();
        degrees.dedup();
        let by_degree: Vec<(Dense<u64>, usize)> = degrees
            .iter()
            .map(|&d| {
                let of_d: Vec<Dense<u64>> = factors.iter().filter(|g| g.degree() == Some(d)).cloned().collect();
                (product(&ring, &of_d), d)
            })
            .collect();
        prop_assert_eq!(distinct_degree(&ring, &f, &Never).unwrap(), by_degree.clone());
        let mut random = Random::new(1);
        for (g, d) in by_degree {
            let split = equal_degree(&ring, &g, d, &mut random, &Never).unwrap();
            let expected: Vec<Dense<u64>> = factors.iter().filter(|h| h.degree() == Some(d)).cloned().collect();
            prop_assert_eq!(sorted(split), sorted(expected));
        }
    }

    /// Random monic square-free polynomials over large word primes: the
    /// factors multiply back, are monic, and pass Rabin's test.
    #[test]
    fn large_primes(p in prop::sample::select(vec![65_537u64, (1 << 61) - 1, (1 << 62) - 57]), seed in any::<u64>(), n in 1usize..9) {
        let ring = fp(p);
        let mut random = Random::new(seed);
        let mut coeffs: Vec<u64> = (0..n).map(|_| random.next() % p).collect();
        coeffs.push(1);
        let f = ring.poly(coeffs);
        prop_assume!(EuclideanDomain::gcd(&ring, &f, &ring.derivative(&f)).degree() == Some(0));
        let factors = factor_squarefree(&ring, &f, &Never).unwrap();
        prop_assert_eq!(product(&ring, &factors), f);
        prop_assert_eq!(&sorted(factors.clone()), &factors);
        for g in &factors {
            prop_assert_eq!(g.lead(), Some(&1));
            prop_assert!(rabin(&ring, g), "an irreducible factor");
        }
    }
}
