//! Tests of square-free decomposition (CAS-179).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use cas_num::algebra::{Field, GaussianIntegers, GaussianRationals, Integers, Rationals, Ring};
use cas_num::{Complex, Integer, Rational};
use proptest::prelude::*;

use crate::dense::{Dense, DenseRing};
use crate::factor::{CharacteristicZero, Factors};
use crate::sparse::{Sparse, SparseRing};

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

fn small_coeffs(len: usize) -> impl Strategy<Value = Vec<Integer>> {
    prop::collection::vec((-9i64..9).prop_map(Integer::from), 1..=len)
}

/// `c ∏ fᵢ^eᵢ`.
fn dense_product<R: CharacteristicZero>(
    ring: &DenseRing<R>,
    s: &Factors<R::Elem, Dense<R::Elem>>,
) -> Dense<R::Elem> {
    s.factors()
        .iter()
        .fold(ring.constant(s.constant().clone()), |acc, (f, e)| {
            ring.mul(&acc, &ring.pow(f, *e))
        })
}

fn sparse_product<R: CharacteristicZero>(
    ring: &SparseRing<R>,
    s: &Factors<R::Elem, Sparse<R::Elem>>,
) -> Sparse<R::Elem> {
    s.factors()
        .iter()
        .fold(ring.constant(s.constant().clone()), |acc, (f, e)| {
            ring.mul(&acc, &ring.pow(f, *e))
        })
}

/// The decomposition of `p`, checked: it multiplies back, its exponents
/// increase, and its factors have positive degrees, are primitive and
/// normalised, square-free, and pairwise coprime.
fn check_dense<R: CharacteristicZero>(
    ring: &DenseRing<R>,
    p: &Dense<R::Elem>,
) -> Result<Factors<R::Elem, Dense<R::Elem>>, TestCaseError> {
    let s = ring.squarefree(p);
    prop_assert_eq!(&dense_product(ring, &s), p);
    let base = ring.base();
    for pair in s.factors().windows(2) {
        prop_assert!(pair[0].1 < pair[1].1, "increasing exponents");
    }
    for (i, (f, e)) in s.factors().iter().enumerate() {
        prop_assert!(*e >= 1);
        prop_assert!(f.degree() > Some(0), "a factor of positive degree");
        prop_assert!(base.is_one(&ring.content(f)), "a primitive factor");
        prop_assert!(
            base.is_one(&base.normal_form(f.lead().unwrap()).1),
            "a normalised factor"
        );
        prop_assert_eq!(
            ring.gcd(f, &ring.derivative(f)).degree(),
            Some(0),
            "a square-free factor"
        );
        for (g, _) in &s.factors()[i + 1..] {
            prop_assert_eq!(ring.gcd(f, g).degree(), Some(0), "coprime factors");
        }
    }
    Ok(s)
}

/// Whether `p` has no variable.
fn is_constant<E>(p: &Sparse<E>) -> bool {
    p.degrees().iter().all(|&d| d == 0)
}

/// As [`check_dense`], for sparse polynomials, whose exponents can repeat:
/// square-free means coprime to the derivative in each variable it has, up
/// to factors without it.
fn check_sparse<R: CharacteristicZero>(
    ring: &SparseRing<R>,
    p: &Sparse<R::Elem>,
) -> Result<Factors<R::Elem, Sparse<R::Elem>>, TestCaseError> {
    let s = ring.squarefree(p);
    prop_assert_eq!(&sparse_product(ring, &s), p);
    let base = ring.base();
    for pair in s.factors().windows(2) {
        prop_assert!(pair[0].1 <= pair[1].1, "sorted exponents");
    }
    for (i, (f, e)) in s.factors().iter().enumerate() {
        prop_assert!(*e >= 1);
        prop_assert!(!is_constant(f), "a factor with a variable");
        prop_assert!(base.is_one(&ring.content(f)), "a primitive factor");
        prop_assert!(
            base.is_one(&base.normal_form(&f.coeffs()[0]).1),
            "a normalised factor"
        );
        for (var, &d) in f.degrees().iter().enumerate() {
            if d > 0 {
                let g = ring.gcd(f, &ring.derivative(f, var));
                prop_assert_eq!(ring.degree(&g, var), Some(0), "a square-free factor");
            }
        }
        for (g, _) in &s.factors()[i + 1..] {
            prop_assert!(is_constant(&ring.gcd(f, g)), "coprime factors");
        }
    }
    Ok(s)
}

/// Musser's square-free decomposition, a reference independent of Yun's.
/// From `a = gcd(q, q')` and `b = q / a`, each round's `c = gcd(a, b)`
/// leaves the part `b / c` and moves on to `a / c` and `c`.
fn musser(ring: &DenseRing<Integers>, p: &Dense<Integer>) -> Vec<(Dense<Integer>, u64)> {
    let q = ring.primitive_part(p);
    let mut out = Vec::new();
    if q.degree().unwrap_or(0) == 0 {
        return out;
    }
    let mut a = ring.gcd(&q, &ring.derivative(&q));
    let mut b = ring.div_exact(&q, &a).unwrap();
    let mut k = 1;
    while b.degree() > Some(0) {
        let c = ring.gcd(&a, &b);
        let f = ring.div_exact(&b, &c).unwrap();
        if f.degree() > Some(0) {
            out.push((ring.primitive_part(&f), k));
        }
        a = ring.div_exact(&a, &c).unwrap();
        b = c;
        k += 1;
    }
    out
}

/// The monic associate over ℚ of an integer polynomial.
fn monic(p: &Dense<Integer>) -> Dense<Rational> {
    let qx = DenseRing::new(Rationals);
    let lead = Rational::from(p.lead().unwrap().clone());
    qx.poly(
        p.coeffs()
            .iter()
            .map(|c| Rational::from(c.clone()) / &lead)
            .collect(),
    )
}

/// Parts grouped by exponent: the product of each group, in increasing
/// order of exponents.
fn grouped<P: Ring>(ring: &P, parts: &[(P::Elem, u64)]) -> Vec<(P::Elem, u64)> {
    let mut exponents: Vec<u64> = parts.iter().map(|&(_, e)| e).collect();
    exponents.sort_unstable();
    exponents.dedup();
    exponents
        .into_iter()
        .map(|e| {
            let product = parts
                .iter()
                .filter(|&&(_, x)| x == e)
                .fold(ring.one(), |acc, (f, _)| ring.mul(&acc, f));
            (product, e)
        })
        .collect()
}

/// `k ∏ fᵢ^eᵢ`.
fn power_product<P: Ring>(ring: &P, k: P::Elem, parts: &[(P::Elem, u64)]) -> P::Elem {
    parts
        .iter()
        .fold(k, |acc, (f, e)| ring.mul(&acc, &ring.pow(f, *e)))
}

/// Distinct roots `a / d`, `d > 0` and coprime to `a`, each with an
/// exponent.
fn roots() -> impl Strategy<Value = Vec<((i64, i64), u64)>> {
    prop::collection::vec(((-12i64..12, 1i64..4), 1u64..5), 1..6).prop_map(|v| {
        let mut out: Vec<((i64, i64), u64)> = Vec::new();
        for ((a, d), e) in v {
            let g = Integer::from(a).gcd(&Integer::from(d)).to_i64().unwrap();
            let root = (a / g, d / g);
            if !out.iter().any(|(r, _)| *r == root) {
                out.push((root, e));
            }
        }
        out
    })
}

#[test]
fn zero_and_constants() {
    let ring = zx();
    let s = ring.squarefree(&ring.zero());
    assert_eq!(s.constant(), &Integer::ZERO);
    assert!(s.factors().is_empty());
    let s = ring.squarefree(&ring.from_i64s(&[-7]));
    assert_eq!(s.constant(), &Integer::from(-7));
    assert!(s.factors().is_empty());
    let xy = SparseRing::new(Integers, 2);
    let s = xy.squarefree(&xy.constant(Integer::from(6)));
    assert_eq!(s.constant(), &Integer::from(6));
    assert!(s.factors().is_empty());
}

/// A polynomial, its constant, and its parts with their exponents.
type Case = (Vec<i64>, i64, Vec<(Vec<i64>, u64)>);

/// giac's answers, and signs and contents.
#[test]
fn dense_examples() {
    let ring = zx();
    let cases: Vec<Case> = vec![
        // sqrfree(12*x^4-24*x^2+12) = 12*(x^2-1)^2
        (vec![12, 0, -24, 0, 12], 12, vec![(vec![-1, 0, 1], 2)]),
        // sqrfree(x^3-x^2-x+1) = (x+1)*(x-1)^2
        (
            vec![1, -1, -1, 1],
            1,
            vec![(vec![1, 1], 1), (vec![-1, 1], 2)],
        ),
        // -(x - 1)^2
        (vec![-1, 2, -1], -1, vec![(vec![-1, 1], 2)]),
        // -6 x^3 (2x + 1): the content and the sign.
        (
            vec![0, 0, 0, -6, -12],
            -6,
            vec![(vec![1, 2], 1), (vec![0, 1], 3)],
        ),
        // x^2 + 1 is square-free.
        (vec![1, 0, 1], 1, vec![(vec![1, 0, 1], 1)]),
        // (x + 1)^10
        (
            vec![1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1],
            1,
            vec![(vec![1, 1], 10)],
        ),
    ];
    for (input, constant, factors) in cases {
        let s = ring.squarefree(&ring.from_i64s(&input));
        assert_eq!(s.constant(), &Integer::from(constant), "{input:?}");
        let expected: Vec<(Dense<Integer>, u64)> = factors
            .iter()
            .map(|(f, e)| (ring.from_i64s(f), *e))
            .collect();
        assert_eq!(s.factors(), expected.as_slice(), "{input:?}");
    }
}

#[test]
fn sparse_examples() {
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let int = |n: i64| ring.constant(Integer::from(n));
    // (x - y)^3 (x + y) y^2: the content y^2 in x, decomposed in y.
    let (xmy, xpy) = (ring.sub(&x, &y), ring.add(&x, &y));
    let p = ring.mul(&ring.mul(&ring.pow(&xmy, 3), &xpy), &ring.pow(&y, 2));
    let s = ring.squarefree(&p);
    assert_eq!(s.constant(), &Integer::ONE);
    assert_eq!(
        s.factors(),
        &[(xpy.clone(), 1), (y.clone(), 2), (xmy.clone(), 3)]
    );
    // 4 x^2 y - 4 y: the content y and x^2 - 1 share the exponent 1, and
    // stay apart.
    let xx1 = ring.sub(&ring.mul(&x, &x), &int(1));
    let p = ring.mul(&int(4), &ring.mul(&y, &xx1));
    let s = ring.squarefree(&p);
    assert_eq!(s.constant(), &Integer::from(4));
    assert_eq!(s.factors(), &[(y.clone(), 1), (xx1, 1)]);
    // x^2 y^2 and (y + 1)^2 x^2: the content's part first, as giac's
    // sqrfree(x^2*y^2+2*x*y^2+y^2) is y^2*(x+1)^2.
    let p = ring.mul(&ring.mul(&x, &x), &ring.mul(&y, &y));
    let s = ring.squarefree(&p);
    assert_eq!(s.factors(), &[(y.clone(), 2), (x.clone(), 2)]);
    let p = ring.mul(&ring.pow(&ring.add(&y, &int(1)), 2), &ring.mul(&x, &x));
    let s = ring.squarefree(&p);
    assert_eq!(s.factors(), &[(ring.add(&y, &int(1)), 2), (x.clone(), 2)]);
    // -(y + 1)^2 (x - 1): the sign of the lexicographic lead.
    let p = ring.neg(&ring.mul(&ring.pow(&ring.add(&y, &int(1)), 2), &ring.sub(&x, &int(1))));
    let s = ring.squarefree(&p);
    assert_eq!(s.constant(), &Integer::from(-1));
    assert_eq!(
        s.factors(),
        &[(ring.sub(&x, &int(1)), 1), (ring.add(&y, &int(1)), 2)]
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// `k f g² h³` over ℤ: the checks, Musser's decomposition, and over ℚ
    /// the monic associates with the leading coefficient.
    #[test]
    fn dense_decompositions(f in small_coeffs(4), g in small_coeffs(4), h in small_coeffs(3), k in -20i64..20, d in 1i64..6) {
        let ring = zx();
        let (f, g, h) = (ring.poly(f), ring.poly(g), ring.poly(h));
        let p = ring.scale(
            &ring.mul(&f, &ring.mul(&ring.pow(&g, 2), &ring.pow(&h, 3))),
            &Integer::from(k),
        );
        let s = check_dense(&ring, &p)?;
        let reference = musser(&ring, &p);
        prop_assert_eq!(s.factors(), reference.as_slice());
        // Over ℚ: p / d.
        let qx = DenseRing::new(Rationals);
        let pq = qx.poly(p.coeffs().iter().map(|c| Rational::new(c.clone(), Integer::from(d))).collect());
        let sq = check_dense(&qx, &pq)?;
        let expected: Vec<(Dense<Rational>, u64)> = s.factors().iter().map(|(f, e)| (monic(f), *e)).collect();
        prop_assert_eq!(sq.factors(), expected.as_slice());
        if let Some(lead) = pq.lead() {
            prop_assert_eq!(sq.constant(), lead);
        }
    }

    /// Products of powers of distinct linear factors and of `x² + b`, which
    /// are irreducible and pairwise coprime, decompose into exactly those
    /// parts.
    #[test]
    fn dense_decompositions_of_known_parts(roots in roots(), quadratics in prop::collection::btree_map(1i64..6, 1u64..4, 0..3), k in 1i64..30) {
        let ring = zx();
        let mut parts: Vec<(Dense<Integer>, u64)> = roots.iter().map(|&((a, d), e)| (ring.from_i64s(&[-a, d]), e)).collect();
        parts.extend(quadratics.iter().map(|(&b, &e)| (ring.from_i64s(&[b, 0, 1]), e)));
        let p = power_product(&ring, ring.constant(Integer::from(k)), &parts);
        let s = check_dense(&ring, &p)?;
        prop_assert_eq!(s.constant(), &Integer::from(k));
        let expected = grouped(&ring, &parts);
        prop_assert_eq!(s.factors(), expected.as_slice());
    }

    /// Over ℤ[i], linear factors with distinct Gaussian roots; over ℚ(i)
    /// the same over 2, monic.
    #[test]
    fn gaussian_decompositions(roots in prop::collection::btree_map((-4i64..4, -4i64..4), 1u64..4, 1..5), k in (-3i64..3, -3i64..3), f in prop::collection::vec((-3i64..3, -3i64..3), 1..4)) {
        prop_assume!(k != (0, 0));
        let gauss = |(a, b): (i64, i64)| Complex::new(Rational::from(a), Rational::from(b));
        let ring = DenseRing::new(GaussianIntegers);
        let parts: Vec<(Dense<Complex>, u64)> = roots.iter().map(|(&r, &e)| (ring.poly(vec![-gauss(r), Complex::ONE]), e)).collect();
        let p = power_product(&ring, ring.constant(gauss(k)), &parts);
        let s = check_dense(&ring, &p)?;
        prop_assert_eq!(s.constant(), &gauss(k));
        let expected = grouped(&ring, &parts);
        prop_assert_eq!(s.factors(), expected.as_slice());
        // Times the square of an arbitrary polynomial.
        let f = ring.poly(f.into_iter().map(gauss).collect());
        check_dense(&ring, &ring.mul(&p, &ring.square(&f)))?;
        // Over ℚ(i).
        let field = DenseRing::new(GaussianRationals);
        let half = GaussianRationals.inv(&Complex::new(Rational::from(2), Rational::ZERO)).unwrap();
        let pq = field.scale(&field.poly(p.coeffs().to_vec()), &half);
        let sq = check_dense(&field, &pq)?;
        let expected = grouped(&field, &parts);
        prop_assert_eq!(sq.factors(), expected.as_slice());
    }
}

/// Terms of `n` variables, small exponents and coefficients.
fn terms(n: usize, len: usize) -> impl Strategy<Value = Vec<(Vec<u64>, i64)>> {
    prop::collection::vec((prop::collection::vec(0u64..3, n), -5i64..5), 1..=len)
}

fn build(ring: &SparseRing<Integers>, t: &[(Vec<u64>, i64)]) -> Sparse<Integer> {
    ring.from_terms(t.iter().map(|(e, c)| (e.clone(), Integer::from(*c))))
}

/// Linear forms `a₀ + a₁ x₀ + … + aₙ xₙ₋₁` with at least one variable.
fn forms(n: usize) -> impl Strategy<Value = Vec<(Vec<i64>, u64)>> {
    prop::collection::vec((prop::collection::vec(-3i64..4, n + 1), 1u64..4), 1..5)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    /// `k f g² h³` in one to three variables over ℤ: the checks, the dense
    /// decomposition in one variable, and over ℚ the monic associates.
    #[test]
    fn sparse_decompositions((n, f, g, h) in (1usize..=3).prop_flat_map(|n| (Just(n), terms(n, 3), terms(n, 3), terms(n, 2))), k in -6i64..6) {
        let ring = SparseRing::new(Integers, n);
        let (f, g, h) = (build(&ring, &f), build(&ring, &g), build(&ring, &h));
        let p = ring.scale(
            &ring.mul(&f, &ring.mul(&ring.pow(&g, 2), &ring.pow(&h, 3))),
            &Integer::from(k),
        );
        let s = check_sparse(&ring, &p)?;
        if n == 1 {
            let zx = zx();
            let dense = zx.squarefree(&ring.to_dense(&p, 0).unwrap());
            prop_assert_eq!(s.constant(), dense.constant());
            let expected: Vec<(Sparse<Integer>, u64)> = dense.factors().iter().map(|(f, e)| (ring.from_dense(f, 0), *e)).collect();
            prop_assert_eq!(s.factors(), expected.as_slice());
        }
        let qs = SparseRing::new(Rationals, n);
        let pq = ring.map_into(&p, &qs, |c| Rational::new(c.clone(), Integer::from(3)));
        let sq = check_sparse(&qs, &pq)?;
        prop_assert_eq!(sq.factors().len(), s.factors().len());
        for ((a, e), (b, x)) in s.factors().iter().zip(sq.factors()) {
            prop_assert_eq!(e, x);
            let lead = Rational::from(a.coeffs()[0].clone());
            prop_assert_eq!(&ring.map_into(a, &qs, |c| Rational::from(c.clone()) / &lead), b);
        }
    }

    /// Products of powers of distinct linear forms in one to four
    /// variables, irreducible and pairwise coprime, decompose into parts
    /// that, multiplied together by exponent, are those of the forms: the
    /// unique decomposition with one part an exponent. Forms without the
    /// first variable come from contents.
    #[test]
    fn sparse_decompositions_of_known_parts((n, forms) in (1usize..=4).prop_flat_map(|n| (Just(n), forms(n))), k in 1i64..10) {
        let ring = SparseRing::new(Integers, n);
        let mut parts: Vec<(Sparse<Integer>, u64)> = Vec::new();
        for (a, e) in forms {
            let form = ring.from_terms(a.iter().enumerate().map(|(i, &c)| {
                let mut x = vec![0; n];
                if i > 0 {
                    x[i - 1] = 1;
                }
                (x, Integer::from(c))
            }));
            if is_constant(&form) {
                continue;
            }
            let form = ring.primitive_part(&form);
            if !parts.iter().any(|(f, _)| *f == form) {
                parts.push((form, e));
            }
        }
        prop_assume!(!parts.is_empty());
        let p = power_product(&ring, ring.constant(Integer::from(k)), &parts);
        let s = check_sparse(&ring, &p)?;
        prop_assert_eq!(s.constant(), &Integer::from(k));
        prop_assert_eq!(grouped(&ring, s.factors()), grouped(&ring, &parts));
    }
}
