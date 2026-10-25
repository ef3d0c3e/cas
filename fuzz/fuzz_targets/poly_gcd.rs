//! Polynomial gcds (CAS-137, CAS-152) of the products `f g` and `f h` of
//! random sparse polynomials in one to three variables, over ℤ and over the
//! word prime `2^61 − 1`: every method that applies (the heuristic, Brown's,
//! Zippel's) gives the gcd and cofactors of subresultants, the reference,
//! and so does the default; the gcd times each cofactor is the product.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Integer;
use cas_num::algebra::{Integers, Ring, WordPrimeField};
use cas_num_fuzz::integer;
use cas_poly::{GcdDomain, GcdMethod, Never, Sparse, SparseRing};
use libfuzzer_sys::fuzz_target;

/// At most this many terms a polynomial, each exponent at most 3: the
/// subresultants stay quick.
const MAX_TERMS: usize = 8;

/// A coefficient: mostly small, sometimes of two limbs.
fn coefficient(u: &mut Unstructured<'_>) -> arbitrary::Result<Integer> {
    Ok(if u.ratio(7u8, 8)? {
        Integer::from(u.int_in_range(-20i64..=20)?)
    } else {
        integer(u, 2)?.0
    })
}

fn poly(u: &mut Unstructured<'_>, ring: &SparseRing<Integers>) -> arbitrary::Result<Sparse<Integer>> {
    let n = ring.nvars();
    let count = u.int_in_range(1..=MAX_TERMS)?;
    let mut terms = Vec::with_capacity(count);
    for _ in 0..count {
        let e: Vec<u64> = (0..n)
            .map(|_| u.int_in_range(0..=3))
            .collect::<arbitrary::Result<_>>()?;
        terms.push((e, coefficient(u)?));
    }
    Ok(ring.from_terms(terms))
}

/// Every method that applies, and the default, against subresultants.
fn check<R: GcdDomain>(ring: &SparseRing<R>, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) {
    let reference = ring
        .gcd_by(a, b, GcdMethod::Subresultant, &Never)
        .unwrap()
        .expect("subresultants always apply");
    let (g, x, y) = &reference;
    assert_eq!(&ring.mul(g, x), a, "g (a / g)");
    assert_eq!(&ring.mul(g, y), b, "g (b / g)");
    for method in [GcdMethod::Heuristic, GcdMethod::Modular, GcdMethod::Zippel] {
        if let Some(parts) = ring.gcd_by(a, b, method, &Never).unwrap() {
            assert_eq!(parts, reference, "{method:?}");
        }
    }
    assert_eq!(ring.gcd_cofactors(a, b), reference, "the default");
}

fn run(u: &mut Unstructured<'_>) -> arbitrary::Result<()> {
    let n = u.int_in_range(1usize..=3)?;
    let ring = SparseRing::new(Integers, n);
    let (f, g, h) = (poly(u, &ring)?, poly(u, &ring)?, poly(u, &ring)?);
    let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
    if a.is_zero() || b.is_zero() {
        return Ok(());
    }
    check(&ring, &a, &b);
    let field = WordPrimeField::new((1 << 61) - 1).expect("a prime");
    let fp = SparseRing::new(field, n);
    let image = |p: &Sparse<Integer>| {
        fp.from_terms(p.terms().map(|(e, c)| (e, fp.base().elem_integer(c))))
    };
    let (ap, bp) = (image(&a), image(&b));
    if !ap.is_zero() && !bp.is_zero() {
        check(&fp, &ap, &bp);
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
