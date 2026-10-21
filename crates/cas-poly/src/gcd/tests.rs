//! Tests of gcds, resultants and discriminants.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use cas_num::algebra::{
    EuclideanDomain, Field, GaussianIntegers, GaussianRationals, Integers, Rationals, Ring,
    WordPrimeField,
};
use cas_num::{Complex, Integer, Natural, Rational};
use proptest::prelude::*;

use super::modular::{self, Primes};
use super::{GcdDomain, GcdMethod, heuristic, resultant, subres};
use std::sync::atomic::AtomicBool;

use crate::cancel::{Cancelled, Never};
use crate::dense::{Dense, DenseRing};
use crate::shared::RationalPolys;
use crate::sparse::{Sparse, SparseRing};

/// The determinant of the Sylvester matrix of `a` and `b` over ℚ, by
/// Gaussian elimination: the definition of the resultant.
fn sylvester(a: &[Integer], b: &[Integer]) -> Rational {
    let q =
        |p: &[Integer]| -> Vec<Rational> { p.iter().map(|c| Rational::from(c.clone())).collect() };
    sylvester_q(&q(a), &q(b))
}

fn sylvester_q(a: &[Rational], b: &[Rational]) -> Rational {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let size = m + n;
    let mut rows: Vec<Vec<Rational>> = Vec::with_capacity(size);
    for (count, p) in [(n, a), (m, b)] {
        for shift in 0..count {
            let mut row = vec![Rational::ZERO; size];
            for (k, c) in p.iter().rev().enumerate() {
                row[shift + k] = c.clone();
            }
            rows.push(row);
        }
    }
    let q = Rationals;
    let mut det = Rational::ONE;
    for col in 0..size {
        let Some(pivot) = (col..size).find(|&r| !rows[r][col].is_zero()) else {
            return Rational::ZERO;
        };
        if pivot != col {
            rows.swap(pivot, col);
            det = -det;
        }
        let p = rows[col][col].clone();
        det = q.mul(&det, &p);
        let pivot_row = rows[col].clone();
        for row in &mut rows[col + 1..] {
            let f = q.div(&row[col], &p).unwrap();
            for (x, y) in row.iter_mut().zip(&pivot_row).skip(col) {
                *x = q.sub(x, &q.mul(&f, y));
            }
        }
    }
    det
}

fn small_coeffs(len: usize) -> impl Strategy<Value = Vec<Integer>> {
    prop::collection::vec((-9i64..9).prop_map(Integer::from), 1..=len)
}

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

/// Terms as exponent vectors and small integer coefficients.
type TermList = Vec<(Vec<u64>, i64)>;

/// Terms of `n` variables, small exponents and coefficients.
fn terms(n: usize, len: usize) -> impl Strategy<Value = TermList> {
    prop::collection::vec((prop::collection::vec(0u64..4, n), -9i64..9), 1..=len)
}

fn build<R: Ring>(ring: &SparseRing<R>, t: &[(Vec<u64>, i64)]) -> Sparse<R::Elem> {
    ring.from_terms(
        t.iter()
            .map(|(e, c)| (e.clone(), ring.base().elem_integer(&Integer::from(*c)))),
    )
}

/// Three polynomials of 1 to 3 variables: a common factor and two
/// cofactors.
fn triples() -> impl Strategy<Value = (usize, [TermList; 3])> {
    (1usize..=3).prop_flat_map(|n| (Just(n), [terms(n, 4), terms(n, 5), terms(n, 5)]))
}

/// Whether `a` and `b` are coprime, without a gcd: no common factor with
/// the first variable they have, by a resultant, and coprime contents in
/// it, recursively.
fn coprime<R: GcdDomain>(ring: &SparseRing<R>, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> bool {
    let (da, db) = (a.degrees(), b.degrees());
    let Some(var) = (0..ring.nvars()).find(|&v| da[v] > 0 || db[v] > 0) else {
        // Constants: their gcd is a unit.
        let base = ring.base();
        let value =
            |p: &Sparse<R::Elem>| p.coeffs().first().cloned().unwrap_or_else(|| base.zero());
        return base.is_one(&base.gcd(&value(a), &value(b)));
    };
    if da[var] > 0 && db[var] > 0 && ring.resultant(a, b, var).is_zero() {
        return false;
    }
    coprime(ring, &ring.content_in(a, var), &ring.content_in(b, var))
}

/// The gcd of `f g` and `f h` by `method`, if it applies: it and its
/// cofactors multiply back, the cofactors are coprime, `f` divides it, and
/// it is normalised.
fn check_gcd<R: GcdDomain>(
    ring: &SparseRing<R>,
    f: &Sparse<R::Elem>,
    g: &Sparse<R::Elem>,
    h: &Sparse<R::Elem>,
    method: GcdMethod,
) -> Result<Option<Sparse<R::Elem>>, TestCaseError> {
    let (a, b) = (ring.mul(f, g), ring.mul(f, h));
    let Some((d, x, y)) = ring.gcd_by(&a, &b, method, &Never).unwrap() else {
        return Ok(None);
    };
    prop_assert_eq!(&ring.mul(&d, &x), &a);
    prop_assert_eq!(&ring.mul(&d, &y), &b);
    if !f.is_zero() {
        prop_assert!(
            ring.div_exact(&d, f).is_some(),
            "the common factor divides the gcd"
        );
    }
    if !d.is_zero() {
        let base = ring.base();
        let lead = &d.coeffs()[0];
        prop_assert!(base.is_one(&base.normal_form(lead).1), "a normalised gcd");
        prop_assert!(coprime(ring, &x, &y), "coprime cofactors");
    }
    Ok(Some(d))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// The modular resultant over ℤ, subresultants, and Euclid's over 𝔽ₚ
    /// and ℚ, against the Sylvester determinant.
    #[test]
    fn dense_resultants_are_sylvester_determinants(a in small_coeffs(7), b in small_coeffs(7), d in 1i64..6) {
        let ring = zx();
        let (pa, pb) = (ring.poly(a), ring.poly(b));
        prop_assume!(!pa.is_zero() && !pb.is_zero());
        let expected = sylvester(pa.coeffs(), pb.coeffs());
        let res = ring.resultant(&pa, &pb);
        prop_assert_eq!(Rational::from(res.clone()), expected.clone());
        prop_assert_eq!(subres::resultant(&ring, &pa, &pb, &Never).unwrap(), res.clone());
        for p in [2u64, 7, 65_537, (1 << 61) - 1] {
            let f = WordPrimeField::new(p).unwrap();
            let fp = DenseRing::new(f.clone());
            let image = |q: &Dense<Integer>| fp.poly(q.coeffs().iter().map(|c| f.elem_integer(c)).collect());
            let (ia, ib) = (image(&pa), image(&pb));
            if ia.degree() == pa.degree() && ib.degree() == pb.degree() {
                prop_assert_eq!(resultant::field(&fp, &ia, &ib, &Never).unwrap(), f.elem_integer(&res));
            }
            if !ia.is_zero() && !ib.is_zero() {
                prop_assert_eq!(fp.resultant(&ia, &ib), subres::resultant(&fp, &ia, &ib, &Never).unwrap());
            }
        }
        // Over ℚ: b over d.
        let qx = DenseRing::new(Rationals);
        let qa: Vec<Rational> = pa.coeffs().iter().map(|c| Rational::from(c.clone())).collect();
        let qb: Vec<Rational> = pb.coeffs().iter().map(|c| Rational::new(c.clone(), Integer::from(d))).collect();
        prop_assert_eq!(qx.resultant(&qx.poly(qa.clone()), &qx.poly(qb.clone())), sylvester_q(&qa, &qb));
    }

    /// Hadamard's bound with coefficients of several limbs.
    #[test]
    fn modular_resultants_of_large_coefficients(
        a in prop::collection::vec((any::<bool>(), prop::collection::vec(any::<u64>(), 1..4)), 1..9),
        b in prop::collection::vec((any::<bool>(), prop::collection::vec(any::<u64>(), 1..4)), 1..9),
    ) {
        let ring = zx();
        let big = |v: &[(bool, Vec<u64>)]| ring.poly(v.iter().map(|(s, l)| Integer::from_parts(*s, Natural::from_vec(l.clone()))).collect());
        let (pa, pb) = (big(&a), big(&b));
        prop_assume!(!pa.is_zero() && !pb.is_zero());
        prop_assert_eq!(ring.resultant(&pa, &pb), subres::resultant(&ring, &pa, &pb, &Never).unwrap());
    }

    /// The resultant in a variable is 0 exactly when the gcd has that
    /// variable.
    #[test]
    fn resultants_vanish_with_common_factors((n, [f, g, h]) in triples(), var in 0usize..3) {
        let ring = SparseRing::new(Integers, n);
        let var = var % n;
        let (f, g, h) = (build(&ring, &f), build(&ring, &g), build(&ring, &h));
        // Often a common factor, sometimes none.
        let f = if ring.degree(&f, var).unwrap_or(0) % 2 == 0 { f } else { ring.one() };
        let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
        prop_assume!(!a.is_zero() && !b.is_zero());
        let common = ring.degree(&ring.gcd(&a, &b), var).unwrap_or(0) > 0;
        prop_assert_eq!(ring.resultant(&a, &b, var).is_zero(), common);
    }

    #[test]
    fn dense_gcds_over_the_integers(f in small_coeffs(4), g in small_coeffs(5), h in small_coeffs(5)) {
        let ring = zx();
        let (f, g, h) = (ring.poly(f), ring.poly(g), ring.poly(h));
        prop_assume!(!f.is_zero() && !g.is_zero() && !h.is_zero());
        let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
        let (d, x, y) = ring.gcd_cofactors(&a, &b);
        prop_assert_eq!(ring.mul(&d, &x), a.clone());
        prop_assert_eq!(ring.mul(&d, &y), b.clone());
        prop_assert!(d.lead().unwrap().is_positive());
        prop_assert!(ring.div_exact(&d, &f).is_some(), "the common factor divides the gcd");
        // Against Euclid over ℚ: the monic associate, and the contents.
        let qx = DenseRing::new(Rationals);
        let to_q = |p: &Dense<Integer>| qx.poly(p.coeffs().iter().map(|c| Rational::from(c.clone())).collect());
        let monic = EuclideanDomain::gcd(&qx, &to_q(&a), &to_q(&b));
        let lead = Rational::from(d.lead().unwrap().clone());
        prop_assert_eq!(qx.scale(&to_q(&d), &Rationals.inv(&lead).unwrap()), monic);
        prop_assert_eq!(ring.content(&d), Integers.gcd(&ring.content(&a), &ring.content(&b)));
    }

    #[test]
    fn sparse_gcds_by_subresultants((n, [f, g, h]) in triples()) {
        let ring = SparseRing::new(Integers, n);
        let (f, g, h) = (build(&ring, &f), build(&ring, &g), build(&ring, &h));
        let d = check_gcd(&ring, &f, &g, &h, GcdMethod::Subresultant)?.expect("always applies");
        // One variable: the dense gcd.
        if n == 1 {
            let zx = zx();
            let dense = |p: &Sparse<Integer>| ring.to_dense(p, 0).unwrap();
            let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
            prop_assert_eq!(dense(&d), zx.gcd_by(&dense(&a), &dense(&b), GcdMethod::Subresultant, &Never).unwrap().unwrap().0);
        }
        prop_assert_eq!(check_gcd(&ring, &f, &g, &h, GcdMethod::Modular)?, Some(d.clone()));
        prop_assert_eq!(ring.gcd(&ring.mul(&f, &g), &ring.mul(&f, &h)), d);
        for p in [2u64, 3, 65_537, (1 << 61) - 1] {
            let fp = SparseRing::new(WordPrimeField::new(p).unwrap(), n);
            let image = |q: &Sparse<Integer>| ring.map_into(q, &fp, |c| fp.base().elem_integer(c));
            let (f, g, h) = (image(&f), image(&g), image(&h));
            let reference = check_gcd(&fp, &f, &g, &h, GcdMethod::Subresultant)?;
            // Brown's evaluation can run out of points in a small field.
            if let Some(x) = check_gcd(&fp, &f, &g, &h, GcdMethod::Modular)? {
                prop_assert_eq!(Some(x), reference.clone());
            }
            prop_assert_eq!(Some(fp.gcd(&fp.mul(&f, &g), &fp.mul(&f, &h))), reference);
        }
    }

    /// The heuristic, where it succeeds, against subresultants, over ℤ,
    /// ℚ, ℤ[i] and ℚ(i).
    #[test]
    fn heuristic_gcds((n, [f, g, h]) in triples(), im in prop::collection::vec(-3i64..3, 5), d in 1i64..5) {
        let ring = SparseRing::new(Integers, n);
        let (f, g, h) = (build(&ring, &f), build(&ring, &g), build(&ring, &h));
        let reference = check_gcd(&ring, &f, &g, &h, GcdMethod::Subresultant)?;
        if let Some(d) = check_gcd(&ring, &f, &g, &h, GcdMethod::Heuristic)? {
            prop_assert_eq!(Some(d), reference.clone());
        }
        prop_assert_eq!(Some(check_gcd(&ring, &f, &g, &h, GcdMethod::Subresultant)?.unwrap()), reference);
        // ℤ[i]: the same polynomials, each coefficient times 1 + k i.
        let zi = SparseRing::new(GaussianIntegers, n);
        let twist = |p: &Sparse<Integer>, k: i64| ring.map_into(p, &zi, |c| Complex::new(Rational::from(c.clone()), Rational::from(c * &Integer::from(k))));
        let (fi, gi, hi) = (twist(&f, im[0]), twist(&g, im[1]), twist(&h, im[2]));
        let reference = check_gcd(&zi, &fi, &gi, &hi, GcdMethod::Subresultant)?;
        if let Some(x) = check_gcd(&zi, &fi, &gi, &hi, GcdMethod::Heuristic)? {
            prop_assert_eq!(Some(x), reference);
        }
        // ℚ and ℚ(i): over a denominator.
        let q = SparseRing::new(Rationals, n);
        let over = |p: &Sparse<Integer>| ring.map_into(p, &q, |c| Rational::new(c.clone(), Integer::from(d)));
        let reference = check_gcd(&q, &over(&f), &over(&g), &over(&h), GcdMethod::Subresultant)?;
        prop_assert_eq!(check_gcd(&q, &over(&f), &over(&g), &over(&h), GcdMethod::Heuristic)?, reference);
        let qi = SparseRing::new(GaussianRationals, n);
        let over_i = |p: &Sparse<Complex>| zi.map_into(p, &qi, |c| Complex::new(c.re() / &Rational::from(Integer::from(d)), c.im().clone()));
        check_gcd(&qi, &over_i(&fi), &over_i(&gi), &over_i(&hi), GcdMethod::Subresultant)?;
    }

    #[test]
    fn dense_heuristic_gcds(f in small_coeffs(5), g in small_coeffs(6), h in small_coeffs(6), k in 1i64..1000) {
        let ring = zx();
        let (f, g, h) = (ring.poly(f), ring.poly(g), ring.poly(h));
        prop_assume!(!f.is_zero() && !g.is_zero() && !h.is_zero());
        // Larger coefficients: f times k.
        let f = ring.scale(&f, &Integer::from(k));
        let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
        let reference = ring.gcd_by(&a, &b, GcdMethod::Subresultant, &Never).unwrap().unwrap();
        if let Some(parts) = ring.gcd_by(&a, &b, GcdMethod::Heuristic, &Never).unwrap() {
            prop_assert_eq!(&parts, &reference);
        }
        prop_assert_eq!(ring.gcd_by(&a, &b, GcdMethod::Modular, &Never).unwrap(), Some(reference.clone()));
        prop_assert_eq!(ring.gcd_cofactors(&a, &b), reference);
    }

    /// Dense gcds over ℤ[i], ℚ, ℚ(i), 𝔽ₚ and ℚ[x] with one denominator,
    /// against subresultants and Euclid.
    #[test]
    fn dense_gcds_over_every_domain(f in small_coeffs(4), g in small_coeffs(5), h in small_coeffs(5), k in -3i64..3, d in 1i64..5) {
        let z = zx();
        let (f, g, h) = (z.poly(f), z.poly(g), z.poly(h));
        prop_assume!(!f.is_zero() && !g.is_zero() && !h.is_zero());
        let (a, b) = (z.mul(&f, &g), z.mul(&f, &h));
        // ℤ[i]: each coefficient c becomes c (1 + k i).
        let zi = DenseRing::new(GaussianIntegers);
        let twist = |p: &Dense<Integer>| zi.poly(p.coeffs().iter().map(|c| Complex::new(Rational::from(c.clone()), Rational::from(c * &Integer::from(k)))).collect());
        let (ai, bi) = (twist(&a), twist(&b));
        let reference = zi.gcd_by(&ai, &bi, GcdMethod::Subresultant, &Never).unwrap().unwrap();
        prop_assert_eq!(zi.gcd_cofactors(&ai, &bi), reference.clone());
        prop_assert_eq!(zi.mul(&reference.0, &reference.1), ai.clone());
        prop_assert!(zi.div_exact(&reference.0, &twist(&f)).is_some());
        // ℚ, against Euclid; and ℚ[x] over one denominator.
        let qx = DenseRing::new(Rationals);
        let to_q = |p: &Dense<Integer>| qx.poly(p.coeffs().iter().map(|c| Rational::new(c.clone(), Integer::from(d))).collect());
        let (aq, bq) = (to_q(&a), to_q(&b));
        let euclid = EuclideanDomain::gcd(&qx, &aq, &bq);
        let (gq, xq, yq) = qx.gcd_cofactors(&aq, &bq);
        prop_assert_eq!(&gq, &euclid);
        prop_assert_eq!(qx.mul(&gq, &xq), aq.clone());
        prop_assert_eq!(qx.mul(&gq, &yq), bq.clone());
        let shared = RationalPolys::new();
        let (sa, sb) = (shared.from_fractions(aq.coeffs()), shared.from_fractions(bq.coeffs()));
        prop_assert_eq!(shared.to_fractions(&EuclideanDomain::gcd(&shared, &sa, &sb)), euclid.coeffs().to_vec());
        prop_assert_eq!(shared.resultant(&sa, &sb), qx.resultant(&aq, &bq));
        prop_assert_eq!(shared.discriminant(&sa), qx.discriminant(&aq));
        // ℚ(i).
        let qi = DenseRing::new(GaussianRationals);
        let over = |p: &Dense<Complex>| qi.poly(p.coeffs().iter().map(|c| Complex::new(c.re() / &Rational::from(Integer::from(d)), c.im().clone())).collect());
        let (aqi, bqi) = (over(&ai), over(&bi));
        prop_assert_eq!(qi.gcd(&aqi, &bqi), EuclideanDomain::gcd(&qi, &aqi, &bqi));
        // 𝔽ₚ: Euclid against subresultants.
        for p in [2u64, 5, 65_537] {
            let fp = DenseRing::new(WordPrimeField::new(p).unwrap());
            let image = |q: &Dense<Integer>| fp.poly(q.coeffs().iter().map(|c| fp.base().elem_integer(c)).collect());
            let (ap, bp) = (image(&a), image(&b));
            if ap.is_zero() || bp.is_zero() {
                continue;
            }
            prop_assert_eq!(fp.gcd_cofactors(&ap, &bp), fp.gcd_by(&ap, &bp, GcdMethod::Subresultant, &Never).unwrap().unwrap());
        }
    }

    /// Coefficients of several limbs, so that Brown's algorithm combines
    /// several primes, in one and in two variables.
    #[test]
    fn gcds_of_large_coefficients(
        f in prop::collection::vec((any::<bool>(), prop::collection::vec(any::<u64>(), 1..4)), 1..5),
        g in small_coeffs(5),
        h in small_coeffs(5),
    ) {
        let z = zx();
        let f = z.poly(f.iter().map(|(s, l)| Integer::from_parts(*s, Natural::from_vec(l.clone()))).collect());
        let (g, h) = (z.poly(g), z.poly(h));
        prop_assume!(!f.is_zero() && !g.is_zero() && !h.is_zero());
        let (a, b) = (z.mul(&f, &g), z.mul(&f, &h));
        let reference = z.gcd_by(&a, &b, GcdMethod::Subresultant, &Never).unwrap().unwrap();
        prop_assert_eq!(z.gcd_by(&a, &b, GcdMethod::Modular, &Never).unwrap(), Some(reference.clone()));
        if let Some(parts) = z.gcd_by(&a, &b, GcdMethod::Heuristic, &Never).unwrap() {
            prop_assert_eq!(&parts, &reference);
        }
        // Two variables: x for the dense variable, times (x + y) in f.
        let ring = SparseRing::new(Integers, 2);
        let lift = |p: &Dense<Integer>| ring.from_dense(p, 0);
        let common = ring.mul(&lift(&f), &ring.add(&ring.var(0), &ring.var(1)));
        let (gs, hs) = (ring.add(&lift(&g), &ring.var(1)), lift(&h));
        let d = check_gcd(&ring, &common, &gs, &hs, GcdMethod::Subresultant)?.unwrap();
        prop_assert_eq!(check_gcd(&ring, &common, &gs, &hs, GcdMethod::Modular)?, Some(d.clone()));
        if let Some(x) = check_gcd(&ring, &common, &gs, &hs, GcdMethod::Heuristic)? {
            prop_assert_eq!(x, d);
        }
    }

    /// Four variables, for the recursion of every method.
    #[test]
    fn gcds_in_four_variables(f in terms(4, 3), g in terms(4, 4), h in terms(4, 4)) {
        let ring = SparseRing::new(Integers, 4);
        let (f, g, h) = (build(&ring, &f), build(&ring, &g), build(&ring, &h));
        let d = check_gcd(&ring, &f, &g, &h, GcdMethod::Subresultant)?.unwrap();
        prop_assert_eq!(check_gcd(&ring, &f, &g, &h, GcdMethod::Modular)?, Some(d.clone()));
        if let Some(x) = check_gcd(&ring, &f, &g, &h, GcdMethod::Heuristic)? {
            prop_assert_eq!(x, d.clone());
        }
        let fp = SparseRing::new(WordPrimeField::new((1 << 61) - 1).unwrap(), 4);
        let image = |q: &Sparse<Integer>| ring.map_into(q, &fp, |c| fp.base().elem_integer(c));
        let reference = check_gcd(&fp, &image(&f), &image(&g), &image(&h), GcdMethod::Subresultant)?;
        prop_assert_eq!(check_gcd(&fp, &image(&f), &image(&g), &image(&h), GcdMethod::Modular)?, reference);
    }

    /// The resultant in the first variable against the dense resultant at
    /// points of the others, where neither leading coefficient vanishes.
    #[test]
    fn sparse_resultants_at_points(
        (n, [a, b]) in (2usize..=3).prop_flat_map(|n| (Just(n), [terms(n, 5), terms(n, 5)])),
        point in prop::collection::vec(-4i64..4, 3),
    ) {
        let ring = SparseRing::new(Integers, n);
        let (a, b) = (build(&ring, &a), build(&ring, &b));
        let res = ring.resultant(&a, &b, 0);
        prop_assert_eq!(ring.degree(&res, 0).unwrap_or(0), 0);
        let values: Vec<Integer> = point[..n].iter().map(|&v| Integer::from(v)).collect();
        let at = |p: &Sparse<Integer>| (1..n).fold(p.clone(), |q, v| ring.eval_var(&q, v, &values[v]));
        let (pa, pb) = (at(&a), at(&b));
        if ring.degree(&pa, 0) != ring.degree(&a, 0) || ring.degree(&pb, 0) != ring.degree(&b, 0) {
            // A leading coefficient vanishes at the point.
            return Ok(());
        }
        let zx = zx();
        let dense = |p: &Sparse<Integer>| ring.to_dense(p, 0).unwrap();
        let expected = zx.resultant(&dense(&pa), &dense(&pb));
        prop_assert_eq!(ring.eval(&res, &values), expected);
    }
}

/// `b² - 4ac` and `-4p³ - 27q²`, and the edge cases of degree 0 and 1.
#[test]
fn discriminants() {
    let ring = SparseRing::new(Integers, 4);
    let var = |v: usize| ring.var(v);
    let int = |n: i64| ring.constant(Integer::from(n));
    let (x, a, b, c) = (var(0), var(1), var(2), var(3));
    let quadratic = ring.add(
        &ring.add(&ring.mul(&a, &ring.mul(&x, &x)), &ring.mul(&b, &x)),
        &c,
    );
    let expected = ring.sub(&ring.mul(&b, &b), &ring.mul(&int(4), &ring.mul(&a, &c)));
    assert_eq!(ring.discriminant(&quadratic, 0), Some(expected));
    // x³ + p x + q, with p and q the variables 1 and 2.
    let (p, q) = (var(1), var(2));
    let cubic = ring.add(&ring.add(&ring.pow(&x, 3), &ring.mul(&p, &x)), &q);
    let expected = ring.sub(
        &ring.mul(&int(-4), &ring.pow(&p, 3)),
        &ring.mul(&int(27), &ring.pow(&q, 2)),
    );
    assert_eq!(ring.discriminant(&cubic, 0), Some(expected));
    // Degree 1 gives 1; constants and zero have none.
    assert_eq!(
        ring.discriminant(&ring.add(&ring.mul(&int(3), &x), &a), 0),
        Some(int(1))
    );
    assert_eq!(ring.discriminant(&a, 0), None);
    assert_eq!(ring.discriminant(&ring.zero(), 0), None);
    // Dense: x⁴ + 1 has 256, 2x² + 3x + 1 has 1.
    let zx = zx();
    assert_eq!(
        zx.discriminant(&zx.from_i64s(&[1, 0, 0, 0, 1])),
        Some(Integer::from(256))
    );
    assert_eq!(
        zx.discriminant(&zx.from_i64s(&[1, 3, 2])),
        Some(Integer::ONE)
    );
    assert_eq!(zx.discriminant(&zx.from_i64s(&[5])), None);
}

/// giac's resultants (CAS-137), where giac follows the definition.
#[test]
fn resultants_of_examples() {
    let zx = zx();
    let r = |a: &[i64], b: &[i64]| zx.resultant(&zx.from_i64s(a), &zx.from_i64s(b));
    // resultant(x^2-2, x^3-3, x) = 1; (x^3+x, x^2+3) = 12 both ways;
    // (2x^2+3, 4x^3-1) = 440.
    assert_eq!(r(&[-2, 0, 1], &[-3, 0, 0, 1]), Integer::ONE);
    assert_eq!(r(&[0, 1, 0, 1], &[3, 0, 1]), Integer::from(12));
    assert_eq!(r(&[3, 0, 1], &[0, 1, 0, 1]), Integer::from(12));
    assert_eq!(r(&[3, 0, 2], &[-1, 0, 0, 4]), Integer::from(440));
    // A constant to the other's degree, where giac gives 1; zero gives 0.
    assert_eq!(r(&[2], &[1, 0, 1]), Integer::from(4));
    assert_eq!(r(&[1, 0, 1], &[2]), Integer::from(4));
    assert_eq!(r(&[2], &[3]), Integer::ONE);
    assert_eq!(r(&[], &[1, 1]), Integer::ZERO);
    // resultant(x y - 1, x^2 + y^2 - 4, x) = y^4 - 4 y^2 + 1.
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let one = ring.one();
    let a = ring.sub(&ring.mul(&x, &y), &one);
    let b = ring.sub(
        &ring.add(&ring.mul(&x, &x), &ring.mul(&y, &y)),
        &ring.constant(Integer::from(4)),
    );
    let expected = ring.add(
        &ring.sub(
            &ring.pow(&y, 4),
            &ring.mul(&ring.constant(Integer::from(4)), &ring.pow(&y, 2)),
        ),
        &one,
    );
    assert_eq!(ring.resultant(&a, &b, 0), expected);
}

/// `x + 2` and `x - 2` at the first point, 6: their values 8 and 4 have the
/// gcd 4, whose digits make `x - 2`, which does not divide `x + 2`. The
/// next point, 16, gives the gcd 1.
#[test]
fn a_heuristic_point_that_fails() {
    let ring = zx();
    let (a, b) = (ring.from_i64s(&[2, 1]), ring.from_i64s(&[-2, 1]));
    let xi = heuristic::first_point::<Integers>(&Natural::from(2u64), &Natural::from(2u64));
    assert_eq!(xi, Integer::from(6));
    let gamma = Integers.gcd(&ring.eval(&a, &xi), &ring.eval(&b, &xi));
    let candidate = heuristic::digits_dense(&ring, &gamma, &xi);
    assert_eq!(candidate, b);
    assert_eq!(ring.div_exact(&a, &candidate), None);
    assert_eq!(heuristic::next_point(&xi), Integer::from(16));
    let (g, x, y) = ring
        .gcd_by(&a, &b, GcdMethod::Heuristic, &Never)
        .unwrap()
        .unwrap();
    assert_eq!((g, x, y), (ring.one(), a, b));
}

/// Past its size limit the heuristic gives up, and `gcd` goes on to
/// Brown's algorithm.
#[test]
fn the_heuristic_gives_up_on_large_values() {
    let ring = zx();
    // Degree 1200 with coefficients of 240 bits: ξ has some 242 bits, and
    // its 1202nd power far more than the limit.
    const { assert!(242 * 1202 > heuristic::MAX_BITS) };
    let big = Integer::from(Natural::power_of_two(239));
    let coeffs = |shift: i64| -> Vec<Integer> {
        (0..=1200i64)
            .map(|k| &big + &Integer::from(k * shift))
            .collect()
    };
    let f = ring.poly(coeffs(3));
    let (a, b) = (
        ring.mul(&f, &ring.from_i64s(&[1, 1])),
        ring.mul(&f, &ring.from_i64s(&[-1, 1])),
    );
    assert!(
        ring.gcd_by(&a, &b, GcdMethod::Heuristic, &Never)
            .unwrap()
            .is_none()
    );
    let (g, x, y) = ring.gcd_cofactors(&a, &b);
    assert!(g == f, "the gcd is f");
    assert_eq!((x, y), (ring.from_i64s(&[1, 1]), ring.from_i64s(&[-1, 1])));
}

/// The first word prime `p₀` is unlucky for `(x + p₀)(x - 1)` and `x (x -
/// 1)`, whose images are equal: the second prime shows it and the
/// combination starts over. With `p₁` instead, the second prime is the
/// unlucky one, and is dropped.
#[test]
fn unlucky_primes() {
    let ring = zx();
    let primes: Vec<u64> = Primes::new().take(2).map(|(p, _)| p).collect();
    let b = ring.from_i64s(&[0, -1, 1]);
    for &p in &primes {
        let shifted = ring.poly(vec![Integer::from(p), Integer::ONE]);
        let a = ring.mul(&shifted, &ring.from_i64s(&[-1, 1]));
        let (g, x, y) = ring
            .gcd_by(&a, &b, GcdMethod::Modular, &Never)
            .unwrap()
            .unwrap();
        assert_eq!(g, ring.from_i64s(&[-1, 1]), "unlucky {p}");
        assert_eq!((x, y), (shifted, ring.x()));
        // Several variables: (x + p y)(x + y) and x (x + y).
        let sparse = SparseRing::new(Integers, 2);
        let (x, y) = (sparse.var(0), sparse.var(1));
        let common = sparse.add(&x, &y);
        let a = sparse.mul(
            &sparse.add(&x, &sparse.scale(&y, &Integer::from(p))),
            &common,
        );
        let b = sparse.mul(&x, &common);
        let (g, _, _) = sparse
            .gcd_by(&a, &b, GcdMethod::Modular, &Never)
            .unwrap()
            .unwrap();
        assert_eq!(g, common, "unlucky {p}, two variables");
    }
}

/// Bad points of evaluation, from the points 1, 2, …: `x + y - 1` and `x`
/// are coprime, but at `y = 1` both are `x`. `(x + y - 2)(x + y)` and
/// `x (x + y)` have the gcd `x + y`, but at `y = 2` both are `x (x + 2)`;
/// the point is dropped and `y = 1, 3` interpolate `x + y`.
#[test]
fn bad_points() {
    let f = WordPrimeField::new(65_537).unwrap();
    let ring = SparseRing::new(f, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let int = |n: i64| ring.constant(ring.base().elem_integer(&Integer::from(n)));
    let a = ring.sub(&ring.add(&x, &y), &int(1));
    let g = modular::pgcd_from(&ring, &a, &x, 1, &Never).unwrap();
    assert_eq!(g, Some(ring.one()));
    let common = ring.add(&x, &y);
    let a = ring.mul(&ring.sub(&common, &int(2)), &common);
    let b = ring.mul(&x, &common);
    assert_eq!(
        modular::pgcd_from(&ring, &a, &b, 1, &Never).unwrap(),
        Some(common)
    );
}

/// Over 𝔽₂ the leading coefficients `y (y + 1)` vanish at both points: the
/// evaluation runs out of points and `gcd` takes subresultants.
#[test]
fn a_field_too_small_for_its_points() {
    let ring = SparseRing::new(WordPrimeField::new(2).unwrap(), 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let lead = ring.mul(&y, &ring.add(&y, &ring.one()));
    let common = ring.add(&ring.mul(&lead, &x), &ring.one());
    let a = ring.mul(&common, &ring.add(&x, &y));
    let b = ring.mul(&common, &ring.add(&x, &ring.one()));
    assert_eq!(
        ring.gcd_by(&a, &b, GcdMethod::Modular, &Never).unwrap(),
        None
    );
    assert_eq!(ring.gcd(&a, &b), common);
}

/// giac's answers where they follow from the normal forms (CAS-DOC-8 §15,
/// detail 5), and zeros.
#[test]
fn normal_forms() {
    let z = zx();
    let p = |c: &[i64]| z.from_i64s(c);
    // gcd(6x^2 - 6, 4x - 4) = 2x - 2, with the signs flipped too.
    assert_eq!(z.gcd(&p(&[-6, 0, 6]), &p(&[-4, 4])), p(&[-2, 2]));
    assert_eq!(z.gcd(&p(&[6, 0, -6]), &p(&[4, -4])), p(&[-2, 2]));
    // lcm(6x^2 - 6, 4x - 4) = 12x^2 - 12; lcm(x^2 - 1, x^2 + 2x + 1).
    assert_eq!(z.lcm(&p(&[-6, 0, 6]), &p(&[-4, 4])), p(&[-12, 0, 12]));
    assert_eq!(z.lcm(&p(&[-1, 0, 1]), &p(&[1, 2, 1])), p(&[-1, -1, 1, 1]));
    // gcd(4, 6x) = 2; zeros.
    assert_eq!(z.gcd(&p(&[4]), &p(&[0, 6])), p(&[2]));
    assert_eq!(z.gcd(&p(&[]), &p(&[2, -2])), p(&[-2, 2]));
    assert_eq!(
        z.gcd_cofactors(&p(&[]), &p(&[2, -2])),
        (p(&[-2, 2]), p(&[]), p(&[-1]))
    );
    assert_eq!(z.gcd_cofactors(&p(&[]), &p(&[])), (p(&[]), p(&[]), p(&[])));
    assert_eq!(z.lcm(&p(&[1, -1]), &p(&[])), p(&[]));
    // Over ℚ the gcd is monic: gcd(x^2 - 1/100, x - 1/10) = x - 1/10, and
    // giac's 10x - 1 is the gcd of the numerators.
    let q = RationalPolys::new();
    let r = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
    let a = q.from_fractions(&[r(-1, 100), r(0, 1), r(1, 1)]);
    let b = q.from_fractions(&[r(-1, 10), r(1, 1)]);
    assert_eq!(
        q.to_fractions(&EuclideanDomain::gcd(&q, &a, &b)),
        vec![r(-1, 10), r(1, 1)]
    );
    assert_eq!(z.gcd(a.numerator(), b.numerator()), p(&[-1, 10]));
    // Several variables: the sign of the leading term, x first.
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let a = ring.sub(&ring.mul(&y, &y), &ring.mul(&x, &x));
    let b = ring.sub(&y, &x);
    assert_eq!(ring.gcd(&a, &b), ring.sub(&x, &y));
    assert_eq!(
        ring.content(&ring.scale(&a, &Integer::from(-6))),
        Integer::from(6)
    );
    assert_eq!(
        ring.primitive_part(&ring.scale(&a, &Integer::from(-6))),
        ring.neg(&a)
    );
    // In x: x y + y = y (x + 1).
    let xy = ring.add(&ring.mul(&x, &y), &y);
    assert_eq!(ring.content_in(&xy, 0), y);
    assert_eq!(ring.primitive_part_in(&xy, 0), ring.add(&x, &ring.one()));
    assert_eq!(
        ring.content_in(&ring.neg(&xy), 1),
        ring.add(&x, &ring.one())
    );
}

/// A raised flag stops every method.
#[test]
fn gcds_stop_when_cancelled() {
    let stop = AtomicBool::new(true);
    let z = zx();
    let (a, b) = (z.from_i64s(&[-1, 0, 1]), z.from_i64s(&[1, 2, 1]));
    for method in [
        GcdMethod::Heuristic,
        GcdMethod::Modular,
        GcdMethod::Subresultant,
    ] {
        assert_eq!(
            z.gcd_by(&a, &b, method, &stop),
            Err(Cancelled),
            "{method:?}"
        );
    }
    assert_eq!(z.gcd_with(&a, &b, &stop), Err(Cancelled));
    assert_eq!(z.resultant_with(&a, &b, &stop), Err(Cancelled));
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    let a = ring.sub(&ring.mul(&x, &x), &ring.mul(&y, &y));
    let b = ring.sub(&x, &y);
    for method in [
        GcdMethod::Heuristic,
        GcdMethod::Modular,
        GcdMethod::Subresultant,
    ] {
        assert_eq!(
            ring.gcd_by(&a, &b, method, &stop),
            Err(Cancelled),
            "{method:?}"
        );
    }
    assert_eq!(ring.resultant_with(&a, &b, 0, &stop), Err(Cancelled));
}

/// Remainder sequences with gaps: degrees that drop by more than one,
/// against the Sylvester determinant and Euclid over ℚ.
#[test]
fn abnormal_remainder_sequences() {
    let z = zx();
    let qx = DenseRing::new(Rationals);
    let to_q = |p: &Dense<Integer>| {
        qx.poly(
            p.coeffs()
                .iter()
                .map(|c| Rational::from(c.clone()))
                .collect(),
        )
    };
    let cases: [(&[i64], &[i64]); 5] = [
        // x^8 + x^6 - 3x^4 - 3x^3 + 8x^2 + 2x - 5 and 3x^6 + 5x^4 - 4x^2 - 9x
        // + 21: Knuth's example (TAOCP 2, §4.6.1), whose remainders skip
        // degrees.
        (&[-5, 2, 8, -3, -3, 0, 1, 0, 1], &[21, -9, -4, 0, 5, 0, 3]),
        (&[1, 0, 0, 0, 0, 1], &[1, 0, 1]),
        (&[1, 0, 0, 1, 0, 0, 1], &[0, 1, 0, 1]),
        (&[2, 0, 0, 0, 3, 0, 0, 1], &[1, 0, 0, 0, 1]),
        (&[-1, 0, 0, 0, 0, 0, 0, 0, 1], &[1, 0, 0, 0, 0, 0, 1]),
    ];
    // Its subresultants end in 260708, the resultant up to sign.
    let knuth = z.resultant(&z.from_i64s(cases[0].0), &z.from_i64s(cases[0].1));
    assert_eq!(knuth.abs(), Integer::from(260_708));
    for (a, b) in cases {
        let (pa, pb) = (z.from_i64s(a), z.from_i64s(b));
        let expected = sylvester(pa.coeffs(), pb.coeffs());
        assert_eq!(
            Rational::from(subres::resultant(&z, &pa, &pb, &Never).unwrap()),
            expected,
            "{a:?}, {b:?}"
        );
        assert_eq!(Rational::from(z.resultant(&pa, &pb)), expected);
        let (g, _, _) = z
            .gcd_by(&pa, &pb, GcdMethod::Subresultant, &Never)
            .unwrap()
            .unwrap();
        let monic = EuclideanDomain::gcd(&qx, &to_q(&pa), &to_q(&pb));
        assert_eq!(g.degree(), monic.degree());
        assert!(z.div_exact(&pa, &g).is_some() && z.div_exact(&pb, &g).is_some());
    }
}

/// Pseudo-random words, the same in every run (`SplitMix64`).
struct Words(u64);

impl Words {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// An integer of `bits` bits, of either sign.
    fn integer(&mut self, bits: u64) -> Integer {
        let limbs: Vec<u64> = (0..bits.div_ceil(64)).map(|_| self.next()).collect();
        let n = Natural::from_vec(limbs) >> (bits.div_ceil(64) * 64 - bits);
        Integer::from_parts(self.next() & 1 == 1, n)
    }
}

/// A workload for instruction counts (CAS-137), as `products_workload`:
/// `CAS_POLY_GCD` is `dense:METHOD:DEGREE:BITS`, the gcd of `f g` and `f h`
/// for pseudo-random `f`, `g`, `h` of that degree and coefficient size, or
/// `sparse:METHOD:N`, the same for `f = (1 + x + y + z)^N + x`, `g` and `h`
/// likewise with `y` and `z`. `METHOD` is `heuristic`, `modular`,
/// `subresultant`, `auto`, or `inputs` to build the inputs alone.
#[test]
#[ignore = "a workload for instruction counts"]
fn gcd_workload() {
    let spec = std::env::var("CAS_POLY_GCD").unwrap_or_else(|_| "dense:auto:100:64".into());
    let parts: Vec<&str> = spec.split(':').collect();
    let method = match parts[1] {
        "heuristic" => Some(GcdMethod::Heuristic),
        "modular" => Some(GcdMethod::Modular),
        "subresultant" => Some(GcdMethod::Subresultant),
        "auto" | "inputs" => None,
        other => panic!("unknown method {other}"),
    };
    let inputs_only = parts[1] == "inputs";
    let mut words = Words(137);
    match parts[0] {
        "dense" => {
            let degree: usize = parts[2].parse().unwrap();
            let bits: u64 = parts[3].parse().unwrap();
            let z = zx();
            let mut random = || z.poly((0..=degree).map(|_| words.integer(bits)).collect());
            let (f, g, h) = (random(), random(), random());
            let (a, b) = (z.mul(&f, &g), z.mul(&f, &h));
            if inputs_only {
                return;
            }
            let (d, _, _) = match method {
                Some(m) => z
                    .gcd_by(&a, &b, m, &Never)
                    .unwrap()
                    .expect("the method applies"),
                None => z.gcd_cofactors(&a, &b),
            };
            assert!(z.div_exact(&d, &f).is_some());
        }
        "sparse" => {
            let n: u64 = parts[2].parse().unwrap();
            let ring = SparseRing::new(Integers, 3);
            let (x, y, w) = (ring.var(0), ring.var(1), ring.var(2));
            let base = ring.pow(&ring.add(&ring.add(&ring.add(&ring.one(), &x), &y), &w), n);
            let (f, g, h) = (
                ring.add(&base, &x),
                ring.add(&base, &y),
                ring.add(&base, &w),
            );
            let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
            if inputs_only {
                return;
            }
            let (d, _, _) = match method {
                Some(m) => ring
                    .gcd_by(&a, &b, m, &Never)
                    .unwrap()
                    .expect("the method applies"),
                None => ring.gcd_cofactors(&a, &b),
            };
            assert_eq!(d, f);
        }
        other => panic!("unknown workload {other}"),
    }
}

/// The first word prime divides the leading coefficient of `p₀ x² + x +
/// 1`: its image drops a degree, and so a resultant modulo `p₀` that is not
/// the image of `res(p₀ x² + x + 1, x - 2) = 4 p₀ + 3`. The modular
/// resultant skips that prime.
#[test]
fn resultants_skip_primes_of_the_leading_coefficients() {
    let z = zx();
    let p0 = Primes::new().next().unwrap().0;
    let a = z.poly(vec![Integer::ONE, Integer::ONE, Integer::from(p0)]);
    let b = z.from_i64s(&[-2, 1]);
    let expected = &(&Integer::from(p0) * &Integer::from(4)) + &Integer::from(3);
    assert_eq!(z.resultant(&a, &b), expected);
    assert_eq!(z.resultant(&b, &a), expected);
    assert_eq!(resultant::integers(&a, &b, &Never).unwrap(), expected);
}
