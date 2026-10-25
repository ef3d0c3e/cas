//! Tests of factoring over ℚ(i) (CAS-186).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use std::cmp::Ordering;

use cas_num::algebra::{EuclideanDomain, Field, GaussianIntegers, GaussianRationals, Ring};
use cas_num::{Complex, Rational};
use proptest::prelude::*;

use crate::factor::Factors;
use crate::sparse::{Sparse, SparseRing};

fn gauss(re: i64, im: i64) -> Complex {
    Complex::new(Rational::from(re), Rational::from(im))
}

/// The polynomial of terms `(exponents, re, im)`.
fn poly(ring: &SparseRing<GaussianIntegers>, terms: &[(&[u64], i64, i64)]) -> Sparse<Complex> {
    ring.from_terms(terms.iter().map(|&(e, re, im)| (e.to_vec(), gauss(re, im))))
}

/// `c ∏ fᵢ^eᵢ`.
fn product(
    ring: &SparseRing<GaussianIntegers>,
    f: &Factors<Complex, Sparse<Complex>>,
) -> Sparse<Complex> {
    f.factors()
        .iter()
        .fold(ring.constant(f.constant().clone()), |acc, (g, e)| {
            ring.mul(&acc, &ring.pow(g, *e))
        })
}

fn total_degree<E>(p: &Sparse<E>) -> u64 {
    (0..p.term_count())
        .map(|k| p.exponents(k).iter().sum())
        .max()
        .unwrap_or(0)
}

/// The order of `factor`: by total degree, then by terms from the leading
/// one, each by exponents, then real and imaginary parts, then by exponent.
fn sort(factors: &mut [(Sparse<Complex>, u64)]) {
    let compare = |a: &Sparse<Complex>, b: &Sparse<Complex>| -> Ordering {
        total_degree(a).cmp(&total_degree(b)).then_with(|| {
            a.terms()
                .zip(b.terms())
                .map(|((e, c), (f, d))| {
                    e.cmp(&f)
                        .then_with(|| c.re().cmp(d.re()))
                        .then_with(|| c.im().cmp(d.im()))
                })
                .find(|o| o.is_ne())
                .unwrap_or_else(|| a.term_count().cmp(&b.term_count()))
        })
    };
    factors.sort_by(|(a, e), (b, x)| compare(a, b).then(e.cmp(x)));
}

/// The factorisation of `p`, checked: it multiplies back, and its factors
/// have variables, are primitive over ℤ\[i\] with their leading coefficients
/// in the first quadrant, are distinct and sorted, and factor as
/// themselves.
fn check(
    ring: &SparseRing<GaussianIntegers>,
    p: &Sparse<Complex>,
) -> Result<Factors<Complex, Sparse<Complex>>, TestCaseError> {
    let f = ring.factor(p);
    prop_assert_eq!(&product(ring, &f), p);
    let mut sorted = f.factors().to_vec();
    sort(&mut sorted);
    prop_assert_eq!(f.factors(), sorted.as_slice(), "sorted factors");
    let base = ring.base();
    for (i, (g, e)) in f.factors().iter().enumerate() {
        prop_assert!(*e >= 1);
        prop_assert!(total_degree(g) > 0, "a factor with a variable");
        prop_assert!(base.is_one(&ring.content(g)), "a primitive factor");
        prop_assert!(
            base.is_one(&base.normal_form(&g.coeffs()[0]).1),
            "a normalised factor"
        );
        prop_assert!(
            f.factors()[i + 1..].iter().all(|(h, _)| h != g),
            "distinct factors"
        );
        let again = ring.factor(g);
        prop_assert!(base.is_one(again.constant()));
        prop_assert_eq!(
            again.factors(),
            &[(g.clone(), 1)][..],
            "an irreducible factor"
        );
    }
    Ok(f)
}

/// Expects the factorisation of the product.
fn expect(
    ring: &SparseRing<GaussianIntegers>,
    constant: &Complex,
    factors: &[(Sparse<Complex>, u64)],
) {
    let mut sorted = factors.to_vec();
    sort(&mut sorted);
    let p = sorted
        .iter()
        .fold(ring.constant(constant.clone()), |acc, (g, e)| {
            ring.mul(&acc, &ring.pow(g, *e))
        });
    let f = ring.factor(&p);
    assert_eq!(f.constant(), constant, "{p:?}");
    assert_eq!(f.factors(), sorted.as_slice(), "{p:?}");
}

/// giac's answers, and polynomials over ℤ that split over ℚ(i).
#[test]
fn examples() {
    let ring = SparseRing::new(GaussianIntegers, 2);
    let x = |c: (i64, i64)| poly(&ring, &[(&[1, 0], 1, 0), (&[0, 0], c.0, c.1)]);
    // factor(x^2+(1+2*i)*x+2*i) = (x+1)*(x+2*i)
    expect(&ring, &gauss(1, 0), &[(x((1, 0)), 1), (x((0, 2)), 1)]);
    // factor(i*x^2+i) = i*(x+i)*(x-i)
    expect(&ring, &gauss(0, 1), &[(x((0, 1)), 1), (x((0, -1)), 1)]);
    // factor((2+i)*x+(1-2*i)) = (2+i)*(x-i)
    expect(&ring, &gauss(2, 1), &[(x((0, -1)), 1)]);
    // factor(x^2+i*x+2) = (x+2*i)*(x-i), factor(2*i*x+2) = 2*i*(x-i)
    expect(&ring, &gauss(1, 0), &[(x((0, 2)), 1), (x((0, -1)), 1)]);
    expect(&ring, &gauss(0, 2), &[(x((0, -1)), 1)]);
    // factor(x^2+2*x+1+i*(x+1)) = (x+1)*(x+1+i), and powers:
    // factor((x+i)^2*(x-i)) = (x+i)^2*(x-i).
    expect(&ring, &gauss(1, 0), &[(x((1, 0)), 1), (x((1, 1)), 1)]);
    expect(&ring, &gauss(1, 0), &[(x((0, 1)), 2), (x((0, -1)), 1)]);
    // x^4 + 4 i stays, as in giac; x^2 + 2 too.
    let quartic = poly(&ring, &[(&[4, 0], 1, 0), (&[0, 0], 0, 4)]);
    expect(&ring, &gauss(1, 0), &[(quartic, 1)]);
    let quadratic = poly(&ring, &[(&[2, 0], 1, 0), (&[0, 0], 2, 0)]);
    expect(&ring, &gauss(1, 0), &[(quadratic, 1)]);
    // Over ℤ, x^4 + 1 = (x^2 - i)(x^2 + i) and x^2 + y^2 = (x - i y)(x + i y).
    let x2 = |im: i64| poly(&ring, &[(&[2, 0], 1, 0), (&[0, 0], 0, im)]);
    expect(&ring, &gauss(1, 0), &[(x2(-1), 1), (x2(1), 1)]);
    let linear = |im: i64| poly(&ring, &[(&[1, 0], 1, 0), (&[0, 1], 0, im)]);
    expect(&ring, &gauss(1, 0), &[(linear(-1), 1), (linear(1), 1)]);
    // A content y + i in x, and the norm's shift past a real factor.
    let y_plus_i = poly(&ring, &[(&[0, 1], 1, 0), (&[0, 0], 0, 1)]);
    let x_plus_y = poly(&ring, &[(&[1, 0], 1, 0), (&[0, 1], 1, 0)]);
    let x2_plus_3 = poly(&ring, &[(&[2, 0], 1, 0), (&[0, 0], 3, 0)]);
    expect(
        &ring,
        &gauss(3, 0),
        &[(y_plus_i, 1), (x_plus_y, 2), (x2_plus_3, 1)],
    );
    // Constants.
    let f = ring.factor(&ring.constant(gauss(0, 5)));
    assert_eq!(f.constant(), &gauss(0, 5));
    assert!(f.factors().is_empty());
    // Over ℚ(i): monic factors, the leading coefficient as the constant.
    let fields = SparseRing::new(GaussianRationals, 2);
    let half = GaussianRationals.inv(&gauss(2, 0)).unwrap();
    let p = fields.scale(&poly(&ring, &[(&[2, 0], 0, 1), (&[0, 0], 0, 1)]), &half);
    let f = fields.factor(&p);
    assert_eq!(
        f.constant(),
        &Complex::new(Rational::ZERO, Rational::new(1.into(), 2.into()))
    );
    assert_eq!(f.factors(), &[(x((0, -1)), 1), (x((0, 1)), 1)]);
}

/// Terms of `n` variables with small Gaussian coefficients.
fn terms(n: usize, len: usize) -> impl Strategy<Value = Vec<(Vec<u64>, i64, i64)>> {
    prop::collection::vec(
        (prop::collection::vec(0u64..3, n), -3i64..4, -3i64..4),
        1..=len,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Products of powers of random Gaussian polynomials in one or two
    /// variables: the checks, and the union of the parts' factorisations.
    #[test]
    fn random_products(
        (n, parts) in (1usize..=2).prop_flat_map(|n| (Just(n), prop::collection::vec((terms(n, 3), 1u64..3), 1..4))),
    ) {
        let ring = SparseRing::new(GaussianIntegers, n);
        let parts: Vec<(Sparse<Complex>, u64)> = parts
            .iter()
            .map(|(t, e)| (ring.from_terms(t.iter().map(|(x, re, im)| (x.clone(), gauss(*re, *im)))), *e))
            .collect();
        prop_assume!(parts.iter().all(|(a, _)| !a.is_zero()));
        let p = parts.iter().fold(ring.one(), |acc, (a, e)| ring.mul(&acc, &ring.pow(a, *e)));
        let f = check(&ring, &p)?;
        let mut constant = Complex::from(1);
        let mut union: Vec<(Sparse<Complex>, u64)> = Vec::new();
        for (a, e) in &parts {
            let fa = ring.factor(a);
            constant = &constant * &GaussianIntegers.pow(fa.constant(), *e);
            for (g, k) in fa.factors() {
                match union.iter_mut().find(|(h, _)| h == g) {
                    Some((_, x)) => *x += k * e,
                    None => union.push((g.clone(), k * e)),
                }
            }
        }
        sort(&mut union);
        prop_assert_eq!(f.constant(), &constant);
        prop_assert_eq!(f.factors(), union.as_slice());
    }
}
