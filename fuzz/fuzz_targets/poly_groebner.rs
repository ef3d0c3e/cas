//! Gröbner bases (CAS-287) of small random ideals in one to three
//! variables, under each monomial order: the reduced basis has the defining
//! properties (every generator reduces to 0 by it, and so does every
//! S-polynomial of two of its elements; each element is monic, and no term
//! of one is divisible by another's leading monomial), and the basis of the
//! basis is itself, the reduced basis being unique. Limits on the work and
//! on the coefficients bound each input (CAS-291).
#![no_main]

use arbitrary::Unstructured;
use cas_num::{Integer, Rational};
use cas_poly::{
    Distributed, GroebnerLimits, MonomialOrder, Never, groebner_basis, groebner_reduce,
    s_polynomial,
};
use libfuzzer_sys::fuzz_target;

/// At most this many generators, of this many terms, each exponent at most
/// 3: the bases stay small.
const MAX_GENERATORS: usize = 4;
const MAX_TERMS: usize = 4;

/// What a computation may do: a step budget let one input run for minutes,
/// a division step merging thousands of terms with growing coefficients.
const LIMITS: GroebnerLimits = GroebnerLimits {
    work: 1 << 20,
    coefficient_bits: 1 << 12,
};

fn poly(
    u: &mut Unstructured<'_>,
    n: usize,
    order: MonomialOrder,
) -> arbitrary::Result<Distributed> {
    let count = u.int_in_range(1..=MAX_TERMS)?;
    let mut terms = Vec::with_capacity(count);
    for _ in 0..count {
        let e: Vec<u64> = (0..n)
            .map(|_| u.int_in_range(0..=3))
            .collect::<arbitrary::Result<_>>()?;
        let num = u.int_in_range(-9i64..=9)?;
        let den = u.int_in_range(1i64..=4)?;
        terms.push((e, Rational::new(Integer::from(num), Integer::from(den))));
    }
    Ok(Distributed::new(terms, order))
}

fn run(u: &mut Unstructured<'_>) -> arbitrary::Result<()> {
    let n = u.int_in_range(1usize..=3)?;
    let order = *u.choose(&[
        MonomialOrder::Lex,
        MonomialOrder::DegLex,
        MonomialOrder::DegRevLex,
    ])?;
    let count = u.int_in_range(1..=MAX_GENERATORS)?;
    let generators = (0..count)
        .map(|_| poly(u, n, order))
        .collect::<arbitrary::Result<Vec<_>>>()?;
    // Past the limits or the bounds: nothing to check.
    let Ok(basis) = groebner_basis(&generators, order, LIMITS, &Never) else {
        return Ok(());
    };
    // The divisions share the work of one computation.
    let checks = generators.len() + basis.len() * basis.len() / 2;
    let share = GroebnerLimits {
        work: LIMITS.work / u64::try_from(checks.max(1)).unwrap_or(u64::MAX),
        ..LIMITS
    };
    let reduces_to_zero =
        |p: &Distributed| groebner_reduce(p, &basis, order, share, &Never).map(|r| r.is_zero());
    for g in &generators {
        match reduces_to_zero(g) {
            Ok(zero) => assert!(zero, "{g:?} by {basis:?}"),
            Err(_) => return Ok(()),
        }
    }
    for (i, a) in basis.iter().enumerate() {
        assert!(a.terms()[0].1.is_one(), "{a:?} monic in {basis:?}");
        for (j, b) in basis.iter().enumerate() {
            if i < j {
                let s = s_polynomial(a, b, order);
                match reduces_to_zero(&s) {
                    Ok(zero) => assert!(zero, "S({a:?}, {b:?}) in {basis:?}"),
                    Err(_) => return Ok(()),
                }
            }
            if i != j {
                let lead = b.leading_monomial().expect("non-zero");
                let divisible = |m: &[u64]| m.iter().zip(lead).all(|(x, y)| y <= x);
                assert!(
                    a.terms().iter().all(|(m, _)| !divisible(m)),
                    "{a:?} reduced by {b:?}"
                );
            }
        }
    }
    if let Ok(again) = groebner_basis(&basis, order, LIMITS, &Never) {
        assert_eq!(again, basis, "the basis of {basis:?}");
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
