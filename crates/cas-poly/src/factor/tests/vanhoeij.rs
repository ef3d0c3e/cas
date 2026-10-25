//! Tests of van Hoeij's recombination (CAS-188), through `factor`, which
//! takes it for more than `THRESHOLD` local factors.

use cas_num::Integer;
use cas_num::algebra::{Integers, Ring};
use proptest::prelude::*;

use crate::dense::{Dense, DenseRing};

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

/// The order of `factor`: by degree, then by coefficients from the
/// constant term up.
fn sorted(mut v: Vec<(Dense<Integer>, u64)>) -> Vec<(Dense<Integer>, u64)> {
    v.sort_by(|(a, e), (b, x)| (a.degree(), a.coeffs(), e).cmp(&(b.degree(), b.coeffs(), x)));
    v
}

/// The cyclotomic polynomial `Φₙ`, by division.
fn cyclotomic(n: usize) -> Dense<Integer> {
    let ring = zx();
    let mut p = ring.sub(&ring.monomial(Integer::ONE, n), &ring.one());
    for d in (1..n).filter(|d| n % d == 0) {
        p = ring.div_exact(&p, &cyclotomic(d)).unwrap();
    }
    p
}

/// The Swinnerton-Dyer polynomial of √2, √3, √5, √7 and √11, of degree 32:
/// irreducible, with at least 16 local factors modulo every prime.
const SWINNERTON_DYER: [i64; 33] = [
    2_000_989_041_197_056,
    0,
    -44_660_812_492_570_624,
    0,
    183_876_928_237_731_840,
    0,
    -255_690_851_718_529_024,
    0,
    172_580_952_324_702_208,
    0,
    -65_892_492_886_671_360,
    0,
    15_459_151_516_270_592,
    0,
    -2_349_014_746_136_576,
    0,
    239_210_760_462_336,
    0,
    -16_665_641_517_056,
    0,
    801_918_722_048,
    0,
    -26_625_650_688,
    0,
    602_397_952,
    0,
    -9_028_096,
    0,
    84864,
    0,
    -448,
    0,
    1,
];

#[test]
fn swinnerton_dyer() {
    let ring = zx();
    let s = ring.from_i64s(&SWINNERTON_DYER);
    assert_eq!(ring.factor(&s).factors(), &[(s.clone(), 1)]);
    // With x^2 - 3, whose roots are among the conjugates' halves.
    let q = ring.from_i64s(&[-3, 0, 1]);
    let f = ring.factor(&ring.mul(&s, &q));
    assert_eq!(f.factors(), &[(q, 1), (s, 1)]);
}

/// `xⁿ − 1` for `n` with many divisors, whose cyclotomic factors split into
/// many local factors each.
#[test]
fn cyclotomic_polynomials() {
    let ring = zx();
    for n in [36, 48, 60, 72, 84, 90, 105, 120] {
        let p = ring.sub(&ring.monomial(Integer::ONE, n), &ring.one());
        let expected = sorted(
            (1..=n)
                .filter(|d| n % d == 0)
                .map(|d| (cyclotomic(d), 1))
                .collect(),
        );
        let f = ring.factor(&p);
        assert_eq!(f.constant(), &Integer::ONE);
        assert_eq!(f.factors(), expected.as_slice(), "x^{n} - 1");
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// Products of 9 to 16 distinct linear factors, with powers and a
    /// constant: more local factors than the threshold, each a true factor.
    #[test]
    fn many_linear_factors(roots in prop::collection::btree_map(-40i64..40, 1u64..3, 9..17), k in 1i64..20) {
        let ring = zx();
        let expected: Vec<(Dense<Integer>, u64)> = roots
            .iter()
            .map(|(&a, &e)| (ring.from_i64s(&[-a, 1]), e))
            .collect();
        let p = expected
            .iter()
            .fold(ring.constant(Integer::from(k)), |acc, (g, e)| ring.mul(&acc, &ring.pow(g, *e)));
        let f = ring.factor(&p);
        prop_assert_eq!(f.constant(), &Integer::from(k));
        let expected = sorted(expected);
        prop_assert_eq!(f.factors(), expected.as_slice());
    }

    /// Products of a few irreducible quadratics `x² + b x + c`, each with
    /// as many local factors as the prime splits it into, and linear ones.
    #[test]
    fn quadratics_and_lines(
        quadratics in prop::collection::btree_set((-9i64..10, 1i64..30), 3..7),
        lines in prop::collection::btree_set(-9i64..10, 3..6),
    ) {
        let ring = zx();
        let mut expected: Vec<(Dense<Integer>, u64)> = Vec::new();
        for &(b, c) in &quadratics {
            // Irreducible: a negative discriminant.
            if b * b - 4 * c < 0 {
                expected.push((ring.from_i64s(&[c, b, 1]), 1));
            }
        }
        expected.extend(lines.iter().map(|&a| (ring.from_i64s(&[a, 1]), 1)));
        let p = expected.iter().fold(ring.one(), |acc, (g, _)| ring.mul(&acc, g));
        let f = ring.factor(&p);
        let expected = sorted(expected);
        prop_assert_eq!(f.factors(), expected.as_slice());
    }
}
