//! Factoring polynomials over ℤ: univariate (CAS-182) and in two or three
//! variables (CAS-184). Products of powers of random polynomials factor
//! into distinct, primitive, normalised factors that multiply back. The
//! factorisation of the product is the union of the parts', and each factor
//! factors as itself.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Integer;
use cas_num::algebra::{Integers, Ring};
use cas_num_fuzz::integer;
use cas_poly::{Dense, DenseRing, Factors, Sparse, SparseRing};
use libfuzzer_sys::fuzz_target;

/// Univariate inputs past this total degree are skipped. Zassenhaus's
/// recombination tries subsets of the local factors, exponentially many in
/// the worst case, until van Hoeij's method (CAS-188).
const MAX_DEGREE: usize = 20;

/// The same for the total degree of multivariate inputs, whose lifting
/// costs grow with each variable.
const MAX_TOTAL_DEGREE: u64 = 12;

/// A coefficient: mostly small, sometimes of several limbs.
fn coefficient(u: &mut Unstructured<'_>) -> arbitrary::Result<Integer> {
    Ok(if u.ratio(3u8, 4)? {
        Integer::from(u.int_in_range(-40i64..=40)?)
    } else {
        integer(u, 3)?.0
    })
}

/// The order of the dense `factor`: by degree, then by coefficients from
/// the constant term up, then by exponent.
fn sort_dense(factors: &mut [(Dense<Integer>, u64)]) {
    factors.sort_by(|(a, e), (b, x)| (a.degree(), a.coeffs(), e).cmp(&(b.degree(), b.coeffs(), x)));
}

fn total_degree<E>(p: &Sparse<E>) -> u64 {
    (0..p.term_count()).map(|k| p.exponents(k).iter().sum()).max().unwrap_or(0)
}

/// The order of the sparse `factor`: by total degree, then by terms from the
/// leading one, then by exponent.
fn sort_sparse(factors: &mut [(Sparse<Integer>, u64)]) {
    let key = |p: &Sparse<Integer>| {
        let terms: Vec<(Vec<u64>, Integer)> = p.terms().map(|(e, c)| (e, c.clone())).collect();
        (total_degree(p), terms)
    };
    factors.sort_by(|(a, e), (b, x)| (key(a), e).cmp(&(key(b), x)));
}

/// The factorisation of `∏ aᵢ^eᵢ` from those of the `aᵢ`.
fn union<P: Clone + PartialEq>(
    parts: &[(P, u64)],
    factor: impl Fn(&P) -> Factors<Integer, P>,
) -> (Integer, Vec<(P, u64)>) {
    let mut constant = Integer::ONE;
    let mut union: Vec<(P, u64)> = Vec::new();
    for (a, e) in parts {
        let fa = factor(a);
        constant = &constant * &fa.constant().pow(*e);
        for (g, k) in fa.factors() {
            match union.iter_mut().find(|(h, _)| h == g) {
                Some((_, x)) => *x += k * e,
                None => union.push((g.clone(), k * e)),
            }
        }
    }
    (constant, union)
}

fn univariate(u: &mut Unstructured<'_>) -> arbitrary::Result<()> {
    let zx = DenseRing::new(Integers);
    let count = u.int_in_range(1u8..=4)?;
    let mut parts: Vec<(Dense<Integer>, u64)> = Vec::new();
    let mut degree = 0;
    for _ in 0..count {
        let (d, e) = (u.int_in_range(0usize..=8)?, u.int_in_range(1u64..=3)?);
        let mut coeffs = Vec::with_capacity(d + 1);
        for _ in 0..=d {
            coeffs.push(coefficient(u)?);
        }
        let a = zx.poly(coeffs);
        if a.is_zero() {
            return Ok(());
        }
        degree += a.degree().unwrap_or(0) * e as usize;
        parts.push((a, e));
    }
    if degree > MAX_DEGREE {
        return Ok(());
    }
    let p = parts.iter().fold(zx.one(), |acc, (a, e)| zx.mul(&acc, &zx.pow(a, *e)));
    let f = zx.factor(&p);
    let back = f
        .factors()
        .iter()
        .fold(zx.constant(f.constant().clone()), |acc, (g, e)| zx.mul(&acc, &zx.pow(g, *e)));
    assert_eq!(back, p, "the factors multiply back");
    for pair in f.factors().windows(2) {
        let ((a, _), (b, _)) = (&pair[0], &pair[1]);
        assert!((a.degree(), a.coeffs()) < (b.degree(), b.coeffs()), "sorted, distinct factors");
    }
    for (g, e) in f.factors() {
        assert!(*e >= 1 && g.degree() > Some(0), "a factor of positive degree");
        assert_eq!(zx.content(g), Integer::ONE, "a primitive factor");
        assert!(g.lead().is_some_and(Integer::is_positive), "a positive leading coefficient");
        assert_eq!(zx.factor(g).factors(), &[(g.clone(), 1)][..], "a factor factors as itself");
    }
    let (constant, mut factors) = union(&parts, |a| zx.factor(a));
    sort_dense(&mut factors);
    assert_eq!(f.constant(), &constant, "the product of the parts' constants");
    assert_eq!(f.factors(), factors.as_slice(), "the union of the parts' factors");
    Ok(())
}

fn multivariate(u: &mut Unstructured<'_>, n: usize) -> arbitrary::Result<()> {
    let ring = SparseRing::new(Integers, n);
    let count = u.int_in_range(1u8..=3)?;
    let mut parts: Vec<(Sparse<Integer>, u64)> = Vec::new();
    let mut degree = 0;
    for _ in 0..count {
        let (terms, e) = (u.int_in_range(1usize..=4)?, u.int_in_range(1u64..=2)?);
        let mut t = Vec::with_capacity(terms);
        for _ in 0..terms {
            let mut exps = Vec::with_capacity(n);
            for _ in 0..n {
                exps.push(u.int_in_range(0u64..=3)?);
            }
            t.push((exps, coefficient(u)?));
        }
        let a = ring.from_terms(t);
        if a.is_zero() {
            return Ok(());
        }
        degree += total_degree(&a) * e;
        parts.push((a, e));
    }
    if degree > MAX_TOTAL_DEGREE {
        return Ok(());
    }
    let p = parts.iter().fold(ring.one(), |acc, (a, e)| ring.mul(&acc, &ring.pow(a, *e)));
    let f = ring.factor(&p);
    let back = f
        .factors()
        .iter()
        .fold(ring.constant(f.constant().clone()), |acc, (g, e)| ring.mul(&acc, &ring.pow(g, *e)));
    assert_eq!(back, p, "the factors multiply back");
    let mut sorted = f.factors().to_vec();
    sort_sparse(&mut sorted);
    assert_eq!(f.factors(), sorted.as_slice(), "sorted factors");
    for (i, (g, e)) in f.factors().iter().enumerate() {
        assert!(*e >= 1 && total_degree(g) > 0, "a factor with a variable");
        assert_eq!(ring.content(g), Integer::ONE, "a primitive factor");
        assert!(g.coeffs()[0].is_positive(), "a positive leading coefficient");
        assert!(f.factors()[i + 1..].iter().all(|(h, _)| h != g), "distinct factors");
        assert_eq!(ring.factor(g).factors(), &[(g.clone(), 1)][..], "a factor factors as itself");
    }
    let (constant, mut factors) = union(&parts, |a| ring.factor(a));
    sort_sparse(&mut factors);
    assert_eq!(f.constant(), &constant, "the product of the parts' constants");
    assert_eq!(f.factors(), factors.as_slice(), "the union of the parts' factors");
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(n) = u.int_in_range(1usize..=3) else {
        return;
    };
    let _ = if n == 1 { univariate(&mut u) } else { multivariate(&mut u, n) };
});
