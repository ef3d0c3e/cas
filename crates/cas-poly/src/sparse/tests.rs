//! Tests of sparse polynomials against a naive reference, a map from
//! exponent vectors to coefficients, over ℤ, ℚ and 𝔽ₚ; the product's
//! kernels and accumulators; the packed fields at their boundaries; the
//! benchmark products of CAS-DOC-8 §13.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use std::collections::BTreeMap;
use std::sync::atomic::AtomicBool;

use cas_num::algebra::{EuclideanDomain, Integers, Rationals, Ring, WordPrimeField};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use super::heap::Overflow;
use super::layout::Layout;
use super::product::{Accumulator, Wide, heap_terms, mul_by, through_arrays};
use super::{Sparse, SparseRing};
use crate::cancel::{Cancelled, Never};
use crate::dense::DenseRing;

/// A polynomial as a map from exponent vectors to non-zero coefficients.
type Naive<E> = BTreeMap<Vec<u64>, E>;

fn naive<R: Ring>(p: &Sparse<R::Elem>) -> Naive<R::Elem> {
    p.terms().map(|(e, c)| (e, c.clone())).collect()
}

fn add_term<R: Ring>(r: &R, p: &mut Naive<R::Elem>, e: Vec<u64>, c: R::Elem) {
    let sum = match p.remove(&e) {
        Some(d) => r.add(&d, &c),
        None => c,
    };
    if !r.is_zero(&sum) {
        p.insert(e, sum);
    }
}

fn naive_add<R: Ring>(r: &R, a: &Naive<R::Elem>, b: &Naive<R::Elem>) -> Naive<R::Elem> {
    let mut s = a.clone();
    for (e, c) in b {
        add_term(r, &mut s, e.clone(), c.clone());
    }
    s
}

fn naive_mul<R: Ring>(r: &R, a: &Naive<R::Elem>, b: &Naive<R::Elem>) -> Naive<R::Elem> {
    let mut p = Naive::new();
    for (e, c) in a {
        for (f, d) in b {
            let m: Vec<u64> = e.iter().zip(f).map(|(x, y)| x + y).collect();
            add_term(r, &mut p, m, r.mul(c, d));
        }
    }
    p
}

fn naive_pow<R: Ring>(r: &R, a: &Naive<R::Elem>, n: u64, nvars: usize) -> Naive<R::Elem> {
    let mut p = Naive::new();
    add_term(r, &mut p, vec![0; nvars], r.one());
    for _ in 0..n {
        p = naive_mul(r, &p, a);
    }
    p
}

/// Division by the leading term in lexicographic order, the largest term
/// first, as the heap does it.
fn naive_div<R: EuclideanDomain>(
    r: &R,
    a: &Naive<R::Elem>,
    b: &Naive<R::Elem>,
) -> (Naive<R::Elem>, Naive<R::Elem>) {
    let (lm, lc) = b.last_key_value().expect("a divisor");
    let mut p = a.clone();
    let (mut q, mut rem) = (Naive::new(), Naive::new());
    while let Some((m, c)) = p.pop_last() {
        let divides = m.iter().zip(lm).all(|(x, y)| x >= y);
        let (t, rest) = if divides {
            r.div_rem(&c, lc)
        } else {
            (r.zero(), c.clone())
        };
        if divides && r.is_zero(&rest) {
            let tm: Vec<u64> = m.iter().zip(lm).map(|(x, y)| x - y).collect();
            for (f, d) in b {
                if f == lm {
                    continue;
                }
                let pm: Vec<u64> = tm.iter().zip(f).map(|(x, y)| x + y).collect();
                add_term(r, &mut p, pm, r.neg(&r.mul(&t, d)));
            }
            add_term(r, &mut q, tm, t);
        } else {
            add_term(r, &mut rem, m, c);
        }
    }
    (q, rem)
}

fn naive_eval<R: Ring>(r: &R, a: &Naive<R::Elem>, point: &[R::Elem]) -> R::Elem {
    a.iter().fold(r.zero(), |s, (e, c)| {
        let t = e
            .iter()
            .zip(point)
            .fold(c.clone(), |t, (&k, x)| r.mul(&t, &r.pow(x, k)));
        r.add(&s, &t)
    })
}

/// The operations whose cost does not grow with the exponents, against the
/// naive reference: sums, products, scaling, division, the derivative,
/// degrees, coefficients, the recursive view, the leading term.
fn check_arithmetic<R: EuclideanDomain + Clone>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
) -> Result<(), TestCaseError> {
    let r = ring.base();
    let n = ring.nvars();
    let (na, nb) = (naive::<R>(a), naive::<R>(b));
    let from = |p: &Naive<R::Elem>| ring.from_terms(p.iter().map(|(e, c)| (e.clone(), c.clone())));
    // The construction round trip and canonical packing.
    prop_assert_eq!(&from(&na), a);
    prop_assert_eq!(ring.add(a, b), from(&naive_add(r, &na, &nb)));
    prop_assert_eq!(ring.add(b, a), from(&naive_add(r, &na, &nb)));
    let neg_b: Naive<R::Elem> = nb.iter().map(|(e, c)| (e.clone(), r.neg(c))).collect();
    prop_assert_eq!(ring.sub(a, b), from(&naive_add(r, &na, &neg_b)));
    prop_assert_eq!(ring.mul(a, b), from(&naive_mul(r, &na, &nb)));
    prop_assert_eq!(ring.mul(b, a), from(&naive_mul(r, &na, &nb)));
    let k = r.elem_integer(&Integer::from(-3));
    let scaled: Naive<R::Elem> = na
        .iter()
        .map(|(e, x)| (e.clone(), r.mul(&k, x)))
        .filter(|(_, x)| !r.is_zero(x))
        .collect();
    prop_assert_eq!(ring.scale(a, &k), from(&scaled));
    // Division, with remainder and exact.
    if !b.is_zero() {
        let (q, rem) = ring.div_rem(a, b);
        let (nq, nrem) = naive_div(r, &na, &nb);
        prop_assert_eq!(&q, &from(&nq));
        prop_assert_eq!(&rem, &from(&nrem));
        prop_assert_eq!(ring.add(&ring.mul(&q, b), &rem), a.clone());
        let ab = ring.mul(a, b);
        prop_assert_eq!(ring.div_exact(&ab, b), Some(a.clone()));
        prop_assert_eq!(ring.div_exact(a, b), nrem.is_empty().then(|| from(&nq)));
    }
    for var in 0..n {
        // The derivative.
        let mut derivative = Naive::new();
        for (e, coeff) in &na {
            if e[var] > 0 {
                let mut f = e.clone();
                f[var] -= 1;
                add_term(
                    r,
                    &mut derivative,
                    f,
                    r.mul(&r.elem_integer(&Integer::from(e[var])), coeff),
                );
            }
        }
        prop_assert_eq!(ring.derivative(a, var), from(&derivative));
        // Degrees, coefficients, the recursive view.
        let degree = na.keys().map(|e| e[var]).max();
        prop_assert_eq!(ring.degree(a, var), degree);
        let parts = ring.to_recursive(a, var);
        prop_assert!(parts.windows(2).all(|w| w[0].0 > w[1].0));
        for (d, coeff) in &parts {
            prop_assert_eq!(ring.degree(coeff, var).unwrap_or(0), 0);
            prop_assert_eq!(coeff, &ring.coeff(a, var, *d));
        }
        prop_assert_eq!(&ring.from_recursive(var, &parts), a);
        if let Some(d) = degree {
            prop_assert_eq!(ring.lead_coeff(a, var), ring.coeff(a, var, d));
        }
    }
    // The leading term is the largest.
    prop_assert_eq!(
        ring.lead_term(a).map(|(e, c)| (e, c.clone())),
        na.last_key_value().map(|(e, c)| (e.clone(), c.clone()))
    );
    Ok(())
}

/// Every operation against the naive reference: [`check_arithmetic`], and
/// powers, evaluation at a point, partial evaluation and substitution, on
/// small exponents.
fn check_ops<R: EuclideanDomain + Clone>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    c: &Sparse<R::Elem>,
    point: &[R::Elem],
) -> Result<(), TestCaseError> {
    check_arithmetic(ring, a, b)?;
    let r = ring.base();
    let n = ring.nvars();
    let (na, nc) = (naive::<R>(a), naive::<R>(c));
    let from = |p: &Naive<R::Elem>| ring.from_terms(p.iter().map(|(e, c)| (e.clone(), c.clone())));
    prop_assert_eq!(ring.pow(c, 3), from(&naive_pow(r, &nc, 3, n)));
    prop_assert_eq!(ring.eval(a, point), naive_eval(r, &na, point));
    for (var, x) in point.iter().enumerate() {
        let mut partial = Naive::new();
        for (e, coeff) in &na {
            let mut f = e.clone();
            f[var] = 0;
            add_term(r, &mut partial, f, r.mul(coeff, &r.pow(x, e[var])));
        }
        prop_assert_eq!(ring.eval_var(a, var, x), from(&partial));
        // x_var -> c.
        let substituted = na.iter().fold(Naive::new(), |s, (e, coeff)| {
            let mut f = e.clone();
            f[var] = 0;
            let mut m = Naive::new();
            m.insert(f, coeff.clone());
            naive_add(r, &s, &naive_mul(r, &m, &naive_pow(r, &nc, e[var], n)))
        });
        prop_assert_eq!(ring.substitute(a, var, c), from(&substituted));
    }
    Ok(())
}

/// Terms as exponent vectors and small integer coefficients.
type TermList = Vec<(Vec<u64>, i64)>;

/// A ring's number of variables, three polynomials' terms, and a point.
type Case = (usize, [TermList; 3], Vec<i64>);

/// Terms of `n` variables, small exponents and coefficients, often
/// cancelling.
fn terms(n: usize) -> impl Strategy<Value = TermList> {
    prop::collection::vec((prop::collection::vec(0u64..5, n), -9i64..9), 0..8)
}

/// Terms with exponents around the field boundaries: 127, the largest of 8
/// bits, and 255, of 9; so operands come in different layouts, and their
/// products cross into wider ones.
fn wide_terms(n: usize) -> impl Strategy<Value = TermList> {
    let exponent = prop_oneof![4 => 0u64..5, 2 => 120u64..140, 1 => 250u64..260];
    prop::collection::vec((prop::collection::vec(exponent, n), -9i64..9), 0..6)
}

fn cases() -> impl Strategy<Value = Case> {
    (1usize..=4).prop_flat_map(|n| {
        (
            Just(n),
            [terms(n), terms(n), terms(n)],
            prop::collection::vec(-5i64..5, n),
        )
    })
}

fn build<R: Ring>(ring: &SparseRing<R>, t: &[(Vec<u64>, i64)]) -> Sparse<R::Elem> {
    ring.from_terms(
        t.iter()
            .map(|(e, c)| (e.clone(), ring.base().elem_integer(&Integer::from(*c)))),
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(160))]

    #[test]
    fn integers((n, [a, b, c], x) in cases()) {
        let ring = SparseRing::new(Integers, n);
        let point: Vec<Integer> = x.iter().map(|&v| Integer::from(v)).collect();
        check_ops(&ring, &build(&ring, &a), &build(&ring, &b), &build(&ring, &c), &point)?;
    }

    #[test]
    fn rationals((n, [a, b, c], x) in cases(), d in 1i64..7) {
        let ring = SparseRing::new(Rationals, n);
        let over = |t: &[(Vec<u64>, i64)]| ring.from_terms(t.iter().map(|(e, c)| (e.clone(), Rational::new(Integer::from(*c), Integer::from(d)))));
        let point: Vec<Rational> = x.iter().map(|&v| Rational::new(Integer::from(v), Integer::from(2))).collect();
        check_ops(&ring, &over(&a), &over(&b), &over(&c), &point)?;
    }

    /// Mixed layouts: operands of 8, 9 and 10 bits, products wider.
    #[test]
    fn exponents_across_fields((n, a, b) in (1usize..=3).prop_flat_map(|n| (Just(n), wide_terms(n), wide_terms(n)))) {
        let ring = SparseRing::new(Integers, n);
        check_arithmetic(&ring, &build(&ring, &a), &build(&ring, &b))?;
        let f = SparseRing::new(WordPrimeField::new(7).unwrap(), n);
        check_arithmetic(&f, &build(&f, &a), &build(&f, &b))?;
    }

    #[test]
    fn prime_field((n, [a, b, c], x) in cases(), p in prop::sample::select(vec![2u64, 3, 7, 65_537, (1 << 61) - 1])) {
        let ring = SparseRing::new(WordPrimeField::new(p).unwrap(), n);
        let point: Vec<u64> = x.iter().map(|&v| ring.base().elem_integer(&Integer::from(v))).collect();
        check_ops(&ring, &build(&ring, &a), &build(&ring, &b), &build(&ring, &c), &point)?;
    }
}

/// Exponents at the edges of the fields: 8-bit fields hold 127, so `x^127 ·
/// x` overflows them, which the guards report; the product repacks into 9
/// bits. Wider fields and words likewise, up to 2^62.
#[test]
fn packed_fields_at_their_boundaries() {
    let ring = SparseRing::new(Integers, 3);
    let mono = |e: &[u64]| ring.monomial(Integer::ONE, e);
    // The guards report an overflow in a layout too narrow.
    let layout = Layout::new(3, 8);
    let (a, b) = (
        mono(&[127, 0, 1]).repacked(layout),
        mono(&[1, 0, 0]).repacked(layout),
    );
    assert_eq!(heap_terms(&ring, &a, &b, layout), Err(Overflow));
    let a = mono(&[126, 0, 1]).repacked(layout);
    assert_eq!(heap_terms(&ring, &a, &b, layout), Ok(mono(&[127, 0, 1])));
    for var in 0..3 {
        let mut e = [0; 3];
        e[var] = 127;
        let mut one = [0; 3];
        one[var] = 1;
        let (a, b) = (mono(&e).repacked(layout), mono(&one).repacked(layout));
        assert_eq!(heap_terms(&ring, &a, &b, layout), Err(Overflow), "x{var}");
    }
    // The public product sizes its fields first.
    for (top, bits) in [
        (127u64, 8),
        (255, 9),
        (1 << 15, 17),
        ((1 << 31) - 1, 33),
        (1 << 61, 64),
    ] {
        let p = ring.add(&mono(&[top, 0, 1]), &mono(&[0, 1, 0]));
        let product = ring.mul(&p, &mono(&[top, 1, 0]));
        let expected = ring.from_terms([
            (vec![2 * top, 1, 1], Integer::ONE),
            (vec![top, 2, 0], Integer::ONE),
        ]);
        assert_eq!(product, expected, "{top}");
        assert_eq!(
            product.layout.bits(),
            bits.max(Layout::holding(3, 2 * top).bits()),
            "{top}"
        );
        // Repacked back down when the large exponents cancel.
        let back = ring.sub(&product, &mono(&[2 * top, 1, 1]));
        assert_eq!(back, mono(&[top, 2, 0]));
        assert_eq!(back.layout, Layout::holding(3, top));
    }
    // Squaring across a field: (x^100 + y)^2 has x^200.
    let p = ring.add(&mono(&[100, 0, 0]), &mono(&[0, 1, 0]));
    let square = ring.pow(&p, 2);
    assert_eq!(ring.degree(&square, 0), Some(200));
    assert_eq!(square.term_count(), 3);
    // Division keeps the dividend's fields.
    let (q, r) = ring.div_rem(&square, &p);
    assert_eq!((q, r), (p.clone(), ring.zero()));
    // Many variables over several words.
    let ring = SparseRing::new(Integers, 17);
    let mut e = vec![0; 17];
    e[16] = 300;
    e[0] = 1;
    let x = ring.monomial(Integer::from(2), &e);
    let y = ring.var(16);
    let xy = ring.mul(&x, &y);
    assert_eq!(xy.exponents(0)[16], 301);
    assert_eq!(ring.div_exact(&xy, &y), Some(x));
}

#[test]
fn examples() {
    let ring = SparseRing::new(Integers, 2);
    let (x, y) = (ring.var(0), ring.var(1));
    // (x + y)(x - y) = x^2 - y^2, in lexicographic order.
    let p = ring.mul(&ring.add(&x, &y), &ring.sub(&x, &y));
    let terms: Vec<(Vec<u64>, Integer)> = p.terms().map(|(e, c)| (e, c.clone())).collect();
    assert_eq!(
        terms,
        [(vec![2, 0], Integer::ONE), (vec![0, 2], Integer::from(-1))]
    );
    // x^2 - y^2 = (x + y)(x - y), and 2x is not divisible by 3x over ℤ.
    assert_eq!(
        ring.div_exact(&p, &ring.add(&x, &y)),
        Some(ring.sub(&x, &y))
    );
    let three_x = ring.scale(&x, &Integer::from(3));
    assert_eq!(
        ring.div_exact(&ring.scale(&x, &Integer::from(2)), &three_x),
        None
    );
    // 2x divided by 3x over ℤ: the term goes to the remainder.
    let (q, r) = ring.div_rem(&ring.scale(&x, &Integer::from(2)), &three_x);
    assert_eq!((q, r), (ring.zero(), ring.scale(&x, &Integer::from(2))));
    // Dense round trip in y.
    let dense = DenseRing::new(Integers);
    let q = ring.add(&ring.pow(&y, 3), &ring.scale(&y, &Integer::from(-2)));
    let d = ring.to_dense(&q, 1).unwrap();
    assert_eq!(d, dense.from_i64s(&[0, -2, 0, 1]));
    assert_eq!(ring.from_dense(&d, 1), q);
    assert_eq!(ring.to_dense(&p, 1), None);
    // x := y + 1 in x^2 - y^2 gives 2y + 1.
    let s = ring.substitute(&p, 0, &ring.add(&y, &ring.one()));
    assert_eq!(s, ring.add(&ring.scale(&y, &Integer::from(2)), &ring.one()));
    // Lexicographic division can raise later exponents past every input's:
    // x^10 / (x - y^200) leaves y^2000, and x / (x + y^200) leaves -y^200,
    // beyond the 8-bit fields of x and x^10.
    let ten = ring.pow(&x, 10);
    let b = ring.sub(&x, &ring.pow(&y, 200));
    let (q, r) = ring.div_rem(&ten, &b);
    assert_eq!(r, ring.pow(&y, 2000));
    assert_eq!(ring.add(&ring.mul(&q, &b), &r), ten);
    let (q, r) = ring.div_rem(&x, &ring.add(&x, &ring.pow(&y, 200)));
    assert_eq!((q, r), (ring.one(), ring.neg(&ring.pow(&y, 200))));
    // An exact division stops at the degrees of the quotient (CAS-137):
    // x - y^200 has more y than x^10; x^2 y + 1 over x - y gives x y first,
    // past the y-degree 0 of a quotient; x^3 y + y^4 = y (x + y)(x^2 - x y
    // + y^2) divides.
    assert_eq!(ring.div_exact(&ten, &b), None);
    let xy = ring.mul(&x, &y);
    let x_minus_y = ring.sub(&x, &y);
    assert_eq!(
        ring.div_exact(&ring.add(&ring.mul(&x, &xy), &ring.one()), &x_minus_y),
        None
    );
    let a = ring.add(&ring.mul(&ring.pow(&x, 3), &y), &ring.pow(&y, 4));
    assert_eq!(
        ring.div_exact(&a, &ring.add(&x, &y))
            .map(|q| ring.mul(&q, &ring.add(&x, &y))),
        Some(a)
    );
    // Constants and zero.
    let none = SparseRing::new(Integers, 0);
    assert_eq!(
        none.mul(
            &none.constant(Integer::from(3)),
            &none.constant(Integer::from(4))
        ),
        none.constant(Integer::from(12))
    );
    assert!(ring.zero().is_zero());
}

/// Every kernel and accumulator that can take `a b` gives the naive
/// product, as does the product's own choice.
fn check_kernels<R: Ring>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
) -> Result<(), TestCaseError> {
    let expected = ring.from_terms(naive_mul(ring.base(), &naive::<R>(a), &naive::<R>(b)));
    prop_assert_eq!(&ring.mul(a, b), &expected);
    for arrays in [false, true] {
        for accumulator in [Accumulator::Words, Accumulator::Mixed, Accumulator::General] {
            if let Some(p) = mul_by(ring, a, b, arrays, accumulator, &Never) {
                prop_assert_eq!(
                    &p.unwrap(),
                    &expected,
                    "arrays {}, {:?}",
                    arrays,
                    accumulator
                );
            }
        }
    }
    Ok(())
}

/// Coefficients for the accumulators: small ones, the edges of `i64`, and
/// some past them, so that products and sums cross `2^126`, `2^127` and
/// `2^128`, and some coefficients have no small preimage.
fn coefficient() -> impl Strategy<Value = Integer> {
    prop_oneof![
        6 => (-9i64..9).prop_map(Integer::from),
        2 => prop::sample::select(vec![i64::MIN, i64::MIN + 1, i64::MAX - 1, i64::MAX])
            .prop_map(Integer::from),
        1 => any::<i64>().prop_map(Integer::from),
        1 => (any::<i64>(), any::<u64>())
            .prop_map(|(h, l)| Integer::from(i128::from(h) << 64 | i128::from(l))),
    ]
}

/// Terms of `n` variables for the kernels: exponents small enough for the
/// arrays to take every product, coefficients from [`coefficient`].
fn kernel_terms(n: usize) -> impl Strategy<Value = Vec<(Vec<u64>, Integer)>> {
    prop::collection::vec((prop::collection::vec(0u64..4, n), coefficient()), 0..14)
}

type KernelCase = (usize, Vec<(Vec<u64>, Integer)>, Vec<(Vec<u64>, Integer)>);

fn kernel_cases() -> impl Strategy<Value = KernelCase> {
    (1usize..=4).prop_flat_map(|n| (Just(n), kernel_terms(n), kernel_terms(n)))
}

/// Products of two `i64`s at their edges, or anywhere.
fn edge() -> impl Strategy<Value = i64> {
    prop_oneof![
        3 => prop::sample::select(vec![i64::MIN, i64::MIN + 1, -1, 0, 1, i64::MAX]),
        1 => any::<i64>(),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(160))]

    /// Over ℤ, coefficients in words and past them.
    #[test]
    fn kernels_over_integers((n, a, b) in kernel_cases()) {
        let ring = SparseRing::new(Integers, n);
        check_kernels(&ring, &ring.from_terms(a), &ring.from_terms(b))?;
    }

    /// Over ℚ, integers and fractions together; over word prime fields,
    /// residues below `2^63`, which have small preimages, and above.
    #[test]
    fn kernels_over_rationals_and_fields(
        (n, a, b) in kernel_cases(),
        d in 1i64..4,
        p in prop::sample::select(vec![7u64, (1 << 61) - 1, u64::MAX - 58]),
    ) {
        let q = SparseRing::new(Rationals, n);
        let over = |t: &[(Vec<u64>, Integer)]| {
            q.from_terms(t.iter().enumerate().map(|(k, (e, c))| {
                let denominator = Integer::from(if k % 2 == 0 { d } else { 1 });
                (e.clone(), Rational::new(c.clone(), denominator))
            }))
        };
        check_kernels(&q, &over(&a), &over(&b))?;
        let f = SparseRing::new(WordPrimeField::new(p).unwrap(), n);
        let into = |t: &[(Vec<u64>, Integer)]| {
            f.from_terms(t.iter().map(|(e, c)| (e.clone(), f.base().elem_integer(c))))
        };
        check_kernels(&f, &into(&a), &into(&b))?;
    }

    /// Sums of products in three words, against the integers.
    #[test]
    fn wide_sums(products in prop::collection::vec((edge(), edge()), 0..48)) {
        let mut wide = Wide::default();
        let mut sum = Integer::ZERO;
        for (x, y) in products {
            wide.add_product(x, y);
            sum = &sum + &(&Integer::from(x) * &Integer::from(y));
            prop_assert_eq!(wide.to_integer(), sum.clone(), "after {} × {}", x, y);
        }
    }
}

/// Sums in three words past `2^128` and down below `−2^128`: `(−2^63)²`
/// sixty-four times, then `−2^63 (2^63 − 1)` two hundred times.
#[test]
fn wide_sums_past_two_words() {
    let mut wide = Wide::default();
    let mut sum = Integer::ZERO;
    let up = &Integer::from(i64::MIN) * &Integer::from(i64::MIN);
    for _ in 0..64 {
        wide.add_product(i64::MIN, i64::MIN);
        sum = &sum + &up;
        assert_eq!(wide.to_integer(), sum);
    }
    assert_eq!(sum, Integer::from(2).pow(132));
    let down = &Integer::from(i64::MIN) * &Integer::from(i64::MAX);
    for _ in 0..200 {
        wide.add_product(i64::MIN, i64::MAX);
        sum = &sum + &down;
        assert_eq!(wide.to_integer(), sum);
    }
    assert!(sum < -Integer::from(2).pow(133));
}

/// `n!`.
fn factorial(n: u64) -> Integer {
    (1..=n).fold(Integer::ONE, |f, k| &f * &Integer::from(k))
}

/// `(1 + x + y + z)^7` squared by every kernel and accumulator, against the
/// multinomial coefficients of the 14th power; the heap's products chain
/// the most in such a dense product.
#[test]
fn dense_products_by_every_kernel() {
    let ring = SparseRing::new(Integers, 3);
    let base = (0..3).fold(ring.one(), |s, v| ring.add(&s, &ring.var(v)));
    let f = ring.pow(&base, 7);
    let mut terms = Vec::new();
    for x in 0..=14u64 {
        for y in 0..=14 - x {
            for z in 0..=14 - x - y {
                let parts = factorial(x) * factorial(y) * factorial(z) * factorial(14 - x - y - z);
                terms.push((vec![x, y, z], factorial(14).div_exact(&parts)));
            }
        }
    }
    let expected = ring.from_terms(terms);
    assert_eq!(ring.mul(&f, &f), expected);
    for arrays in [false, true] {
        for accumulator in [Accumulator::Words, Accumulator::Mixed, Accumulator::General] {
            let p = mul_by(&ring, &f, &f, arrays, accumulator, &Never);
            assert_eq!(p, Some(Ok(expected.clone())), "{arrays} {accumulator:?}");
        }
    }
}

/// Products dense in their degree box go through the arrays, sparse ones
/// through the heap; a chunk past its limit cannot take a product.
#[test]
fn the_kernel_suits_the_product() {
    // Fateman's product at the sixth power: 13³ sums a chunk, 13 chunks,
    // for 210 × 210 products of terms.
    let ring = SparseRing::new(Integers, 4);
    let base = (0..4).fold(ring.one(), |s, v| ring.add(&s, &ring.var(v)));
    let f = ring.pow(&base, 6);
    assert!(through_arrays(&f, &ring.add(&f, &ring.one())));
    // Monagan and Pearce's at the fourth power: 25 chunks of 17³ 25 sums
    // for 126 × 126 products.
    let ring = SparseRing::new(Integers, 5);
    let m = |c: i64, v: usize, e: u64| {
        let mut x = vec![0; 5];
        x[v] = e;
        ring.monomial(Integer::from(c), &x)
    };
    let sum = |parts: &[Sparse<Integer>]| parts.iter().fold(ring.one(), |s, p| ring.add(&s, p));
    let f = ring.pow(
        &sum(&[m(1, 0, 1), m(1, 1, 1), m(2, 2, 2), m(3, 3, 3), m(5, 4, 5)]),
        4,
    );
    let g = ring.pow(
        &sum(&[m(1, 4, 1), m(1, 3, 1), m(2, 2, 2), m(3, 1, 3), m(5, 0, 5)]),
        4,
    );
    assert!(!through_arrays(&f, &g));
    // Dense in one variable, and with a gap.
    let ring = SparseRing::new(Integers, 1);
    let dense = ring.from_terms((0..20).map(|k| (vec![k], Integer::from(k + 1))));
    assert!(through_arrays(&dense, &dense));
    let gap = ring.add(&ring.one(), &ring.monomial(Integer::ONE, &[1000]));
    assert!(!through_arrays(&gap, &gap));
    // A chunk of 300001 sums is past the limit.
    let ring = SparseRing::new(Integers, 2);
    let (a, b) = (
        ring.monomial(Integer::ONE, &[0, 200_000]),
        ring.monomial(Integer::ONE, &[0, 100_000]),
    );
    assert_eq!(
        mul_by(&ring, &a, &b, true, Accumulator::General, &Never),
        None
    );
    assert_eq!(ring.mul(&a, &b), ring.monomial(Integer::ONE, &[0, 300_000]));
}

/// Both kernels stop when cancelled, whatever they accumulate in.
#[test]
fn cancelled_products() {
    let ring = SparseRing::new(Integers, 3);
    let base = (0..3).fold(ring.one(), |s, v| ring.add(&s, &ring.var(v)));
    // 455 terms: 207025 products of terms, 2925 monomials.
    let f = ring.pow(&base, 12);
    let stop = AtomicBool::new(true);
    for arrays in [false, true] {
        for accumulator in [Accumulator::Words, Accumulator::Mixed, Accumulator::General] {
            let p = mul_by(&ring, &f, &f, arrays, accumulator, &stop);
            assert_eq!(p, Some(Err(Cancelled)), "{arrays} {accumulator:?}");
        }
    }
    assert_eq!(ring.mul_with(&f, &f, &stop), Err(Cancelled));
}

/// A workload for instruction counts (CAS-136), as `products_workload`:
/// `CAS_POLY_SPARSE` is `fateman` or `pearce`, the products of CAS-DOC-8
/// §13 over ℤ, or `fateman-inputs` and `pearce-inputs`, their factors
/// alone, to subtract.
#[test]
#[ignore = "a workload for instruction counts"]
fn sparse_workload() {
    let Ok(which) = std::env::var("CAS_POLY_SPARSE") else {
        return;
    };
    let (name, inputs) = match which.strip_suffix("-inputs") {
        Some(name) => (name, true),
        None => (which.as_str(), false),
    };
    let (n, f, g) = match name {
        // f = (1 + x + y + z + t)^20, f (f + 1).
        "fateman" => {
            let ring = SparseRing::new(Integers, 4);
            let base = (0..4).fold(ring.one(), |s, v| ring.add(&s, &ring.var(v)));
            let f = ring.pow(&base, 20);
            let g = ring.add(&f, &ring.one());
            (4, f, g)
        }
        // (1 + x + y + 2z^2 + 3t^3 + 5u^5)^12 (1 + u + t + 2z^2 + 3y^3 + 5x^5)^12.
        "pearce" => {
            let ring = SparseRing::new(Integers, 5);
            let m = |c: i64, v: usize, e: u64| {
                let mut x = vec![0; 5];
                x[v] = e;
                ring.monomial(Integer::from(c), &x)
            };
            let sum =
                |parts: &[Sparse<Integer>]| parts.iter().fold(ring.one(), |s, p| ring.add(&s, p));
            let f = ring.pow(
                &sum(&[m(1, 0, 1), m(1, 1, 1), m(2, 2, 2), m(3, 3, 3), m(5, 4, 5)]),
                12,
            );
            let g = ring.pow(
                &sum(&[m(1, 4, 1), m(1, 3, 1), m(2, 2, 2), m(3, 1, 3), m(5, 0, 5)]),
                12,
            );
            (5, f, g)
        }
        _ => panic!("CAS_POLY_SPARSE=fateman or pearce"),
    };
    if inputs {
        return;
    }
    let ring = SparseRing::new(Integers, n);
    let p = ring.mul(&f, &g);
    eprintln!(
        "{which}: {} x {} terms -> {}",
        f.term_count(),
        g.term_count(),
        p.term_count()
    );
}
