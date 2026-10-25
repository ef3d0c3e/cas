//! Sparse products (CAS-136, CAS-155) against the naive product, over ℤ, ℚ
//! and word prime fields, in one to four variables: factors dense in a
//! small degree box, which the arrays take, or with some exponents spread
//! far, which the heap takes; coefficients small, at the edges of `i64`, or
//! of up to three limbs, which the three-word, mixed and ring accumulators
//! take. The product is the naive one either way round, and divides back
//! exactly, through the division's heap.
#![no_main]

use std::collections::BTreeMap;

use arbitrary::Unstructured;
use cas_num::algebra::{EuclideanDomain, Integers, Rationals, Ring, WordPrimeField};
use cas_num::{Integer, Rational};
use cas_num_fuzz::integer;
use cas_poly::{Sparse, SparseRing};
use libfuzzer_sys::fuzz_target;

/// At most this many terms a factor, for the naive product.
const MAX_TERMS: usize = 64;

/// A coefficient: small, at the edges of `i64`, any `i64`, or of up to
/// three limbs.
fn coefficient(u: &mut Unstructured<'_>) -> arbitrary::Result<Integer> {
    Ok(match u.int_in_range(0u8..=7)? {
        0..=3 => Integer::from(u.int_in_range(-20i64..=20)?),
        4 => Integer::from(*u.choose(&[i64::MIN, i64::MIN + 1, i64::MAX - 1, i64::MAX])?),
        5 => Integer::from(u.arbitrary::<i64>()?),
        _ => integer(u, 3)?.0,
    })
}

type Terms = Vec<(Vec<u64>, Integer)>;

/// Terms in `n` variables, exponents below `width`, now and then one up to
/// `2^20`.
fn terms(u: &mut Unstructured<'_>, n: usize, width: u64) -> arbitrary::Result<Terms> {
    let count = u.int_in_range(0..=MAX_TERMS)?;
    let mut terms = Vec::with_capacity(count);
    for _ in 0..count {
        let mut e = Vec::with_capacity(n);
        for _ in 0..n {
            e.push(if u.ratio(1u8, 32)? {
                u.int_in_range(0..=1 << 20)?
            } else {
                u.int_in_range(0..=width - 1)?
            });
        }
        terms.push((e, coefficient(u)?));
    }
    Ok(terms)
}

/// `a b` is the naive product, as is `b a`, and `a b / b = a`.
fn check<R: EuclideanDomain>(ring: &SparseRing<R>, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) {
    let r = ring.base();
    let mut naive: BTreeMap<Vec<u64>, R::Elem> = BTreeMap::new();
    for (e, x) in a.terms() {
        for (f, y) in b.terms() {
            let m: Vec<u64> = e.iter().zip(&f).map(|(p, q)| p + q).collect();
            let c = r.mul(x, y);
            let sum = match naive.remove(&m) {
                Some(d) => r.add(&d, &c),
                None => c,
            };
            if !r.is_zero(&sum) {
                naive.insert(m, sum);
            }
        }
    }
    let expected = ring.from_terms(naive);
    let p = ring.mul(a, b);
    assert_eq!(p, expected, "a b");
    assert_eq!(ring.mul(b, a), expected, "b a");
    if !b.is_zero() {
        assert_eq!(ring.div_exact(&p, b).as_ref(), Some(a), "a b / b");
    }
}

fn run(u: &mut Unstructured<'_>) -> arbitrary::Result<()> {
    let n = u.int_in_range(1usize..=4)?;
    let width = u.int_in_range(1u64..=6)?;
    let (ta, tb) = (terms(u, n, width)?, terms(u, n, width)?);
    match u.int_in_range(0u8..=2)? {
        0 => {
            let ring = SparseRing::new(Integers, n);
            check(&ring, &ring.from_terms(ta), &ring.from_terms(tb));
        }
        1 => {
            // Every other term over a denominator: integers and fractions.
            let d = Integer::from(u.int_in_range(1i64..=5)?);
            let ring = SparseRing::new(Rationals, n);
            let over = |t: Terms| {
                ring.from_terms(t.into_iter().enumerate().map(|(k, (e, c))| {
                    let denominator = if k % 2 == 0 { d.clone() } else { Integer::ONE };
                    (e, Rational::new(c, denominator))
                }))
            };
            check(&ring, &over(ta), &over(tb));
        }
        _ => {
            // Residues below 2^63, which have small preimages, and above.
            let p = *u.choose(&[2u64, 7, 65_537, (1 << 61) - 1, u64::MAX - 58])?;
            let ring = SparseRing::new(WordPrimeField::new(p).expect("a prime"), n);
            let into = |t: Terms| {
                ring.from_terms(
                    t.into_iter()
                        .map(|(e, c)| (e, ring.base().elem_integer(&c))),
                )
            };
            check(&ring, &into(ta), &into(tb));
        }
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
