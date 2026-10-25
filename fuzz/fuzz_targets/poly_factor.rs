//! Factoring univariate polynomials over ℤ (CAS-182). Products of powers
//! of random polynomials factor into distinct, primitive, normalised factors
//! that multiply back. The factorisation of the product is the union of the
//! parts', and each factor factors as itself.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Integer;
use cas_num::algebra::{Integers, Ring};
use cas_num_fuzz::integer;
use cas_poly::{Dense, DenseRing};
use libfuzzer_sys::fuzz_target;

/// Inputs past this total degree are skipped. Zassenhaus's recombination
/// tries subsets of the local factors, exponentially many in the worst case,
/// until van Hoeij's method (CAS-188).
const MAX_DEGREE: usize = 20;

/// A coefficient: mostly small, sometimes of several limbs.
fn coefficient(u: &mut Unstructured<'_>) -> arbitrary::Result<Integer> {
    Ok(if u.ratio(3u8, 4)? {
        Integer::from(u.int_in_range(-40i64..=40)?)
    } else {
        integer(u, 3)?.0
    })
}

/// The order of `factor`: by degree, then by coefficients from the constant
/// term up, then by exponent.
fn sort(factors: &mut [(Dense<Integer>, u64)]) {
    factors.sort_by(|(a, e), (b, x)| (a.degree(), a.coeffs(), e).cmp(&(b.degree(), b.coeffs(), x)));
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let zx = DenseRing::new(Integers);
    let Ok(count) = u.int_in_range(1u8..=4) else {
        return;
    };
    let mut parts: Vec<(Dense<Integer>, u64)> = Vec::new();
    let mut degree = 0;
    for _ in 0..count {
        let (Ok(d), Ok(e)) = (u.int_in_range(0usize..=8), u.int_in_range(1u64..=3)) else {
            return;
        };
        let mut coeffs = Vec::with_capacity(d + 1);
        for _ in 0..=d {
            let Ok(c) = coefficient(&mut u) else {
                return;
            };
            coeffs.push(c);
        }
        let a = zx.poly(coeffs);
        if a.is_zero() {
            return;
        }
        degree += a.degree().unwrap_or(0) * e as usize;
        parts.push((a, e));
    }
    if degree > MAX_DEGREE {
        return;
    }
    let p = parts
        .iter()
        .fold(zx.one(), |acc, (a, e)| zx.mul(&acc, &zx.pow(a, *e)));
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
    // The union of the parts' factorisations.
    let mut constant = Integer::ONE;
    let mut union: Vec<(Dense<Integer>, u64)> = Vec::new();
    for (a, e) in &parts {
        let fa = zx.factor(a);
        constant = &constant * &fa.constant().pow(*e);
        for (g, k) in fa.factors() {
            match union.iter_mut().find(|(h, _)| h == g) {
                Some((_, x)) => *x += k * e,
                None => union.push((g.clone(), k * e)),
            }
        }
    }
    sort(&mut union);
    assert_eq!(f.constant(), &constant, "the product of the parts' constants");
    assert_eq!(f.factors(), union.as_slice(), "the union of the parts' factors");
});
