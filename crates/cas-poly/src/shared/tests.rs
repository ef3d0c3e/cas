//! Tests of the shared-denominator polynomials over ℚ and ℚ(i), against
//! coefficient-wise polynomials over [`Rationals`] and
//! [`GaussianRationals`].

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use cas_num::algebra::{EuclideanDomain, GaussianRationals, Rationals, Ring};
use cas_num::{Complex, Integer, Natural, Rational};
use proptest::prelude::*;

use super::{GaussianRationalPolys, Numerators, RationalPolys, Shared, SharedRing};
use crate::dense::{Coefficients, Dense, DenseRing};

/// A rational of small or large parts.
fn rational() -> impl Strategy<Value = Rational> {
    let integer = prop_oneof![
        (-30i64..30).prop_map(Integer::from),
        (any::<bool>(), prop::collection::vec(any::<u64>(), 1..=2))
            .prop_map(|(negative, l)| Integer::from_parts(negative, Natural::from_vec(l))),
    ];
    (integer, 1u64..1000).prop_map(|(n, d)| Rational::new(n, Integer::from(d)))
}

fn rationals(len: usize) -> impl Strategy<Value = Vec<Rational>> {
    prop::collection::vec(rational(), 0..=len)
}

fn complexes(len: usize) -> impl Strategy<Value = Vec<Complex>> {
    prop::collection::vec(
        (rational(), rational()).prop_map(|(re, im)| Complex::new(re, im)),
        0..=len,
    )
}

/// Reduced: a positive denominator coprime to the integer parts of the
/// numerator, and 0 over 1.
fn reduced<N: Numerators>(p: &Shared<N::Elem>) -> Result<(), TestCaseError> {
    prop_assert!(!p.denominator().is_zero());
    if p.is_zero() {
        prop_assert!(p.denominator().is_one());
    } else {
        let g = p
            .numerator()
            .coeffs()
            .iter()
            .fold(p.denominator().clone(), |g, c| N::gcd_parts(&g, c));
        prop_assert!(g.is_one());
    }
    Ok(())
}

/// Every operation of `SharedRing<N>` against coefficient-wise arithmetic
/// in `reference`, from the same coefficients.
fn agrees<N, F>(
    shared: &SharedRing<N>,
    reference: &DenseRing<F>,
    a: &[F::Elem],
    b: &[F::Elem],
    x: &F::Elem,
) -> Result<(), TestCaseError>
where
    N: Numerators<Fractions = F>,
    F: Coefficients + cas_num::algebra::Field + Default,
{
    let (sa, sb) = (shared.from_fractions(a), shared.from_fractions(b));
    let (ra, rb) = (reference.poly(a.to_vec()), reference.poly(b.to_vec()));
    let same = |s: &Shared<N::Elem>, r: &Dense<F::Elem>| -> Result<(), TestCaseError> {
        reduced::<N>(s)?;
        prop_assert_eq!(&shared.to_fractions(s), r.coeffs());
        Ok(())
    };
    same(&sa, &ra)?;
    same(&shared.add(&sa, &sb), &reference.add(&ra, &rb))?;
    same(&shared.sub(&sa, &sb), &reference.sub(&ra, &rb))?;
    same(&shared.neg(&sa), &reference.neg(&ra))?;
    same(&shared.mul(&sa, &sb), &reference.mul(&ra, &rb))?;
    same(&shared.square(&sa), &reference.square(&ra))?;
    same(&shared.derivative(&sa), &reference.derivative(&ra))?;
    same(&shared.scale(&sa, x), &reference.scale(&ra, x))?;
    prop_assert_eq!(shared.eval(&sa, x), reference.eval(&ra, x));
    if !rb.is_zero() {
        let (q, r) = shared.div_rem(&sa, &sb);
        let (rq, rr) = reference.div_rem(&ra, &rb);
        same(&q, &rq)?;
        same(&r, &rr)?;
        let (monic, unit) = shared.normal_form(&sb);
        let (rmonic, runit) = reference.normal_form(&rb);
        same(&monic, &rmonic)?;
        same(&unit, &runit)?;
        // (a b) / b = a.
        let (q, r) = shared.div_rem(&shared.mul(&sa, &sb), &sb);
        prop_assert_eq!(q, sa.clone());
        prop_assert!(r.is_zero());
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn rational_polynomials(a in rationals(14), b in rationals(9), x in rational()) {
        agrees(&RationalPolys::new(), &DenseRing::new(Rationals), &a, &b, &x)?;
    }

    #[test]
    fn gaussian_rational_polynomials(a in complexes(10), b in complexes(6), x in (rational(), rational())) {
        let x = Complex::new(x.0, x.1);
        agrees(&GaussianRationalPolys::new(), &DenseRing::new(GaussianRationals), &a, &b, &x)?;
    }

    /// The ring laws, and the gcd of the Euclidean algorithm, on ℚ[x].
    #[test]
    fn rational_ring_laws(a in rationals(8), b in rationals(8), c in rationals(4)) {
        let ring = RationalPolys::new();
        let (a, b, c) = (ring.from_fractions(&a), ring.from_fractions(&b), ring.from_fractions(&c));
        prop_assert_eq!(ring.mul(&a, &ring.add(&b, &c)), ring.add(&ring.mul(&a, &b), &ring.mul(&a, &c)));
        prop_assert_eq!(ring.mul(&ring.mul(&a, &b), &c), ring.mul(&a, &ring.mul(&b, &c)));
        prop_assert_eq!(ring.add(&a, &ring.zero()), a.clone());
        prop_assert_eq!(ring.mul(&a, &ring.one()), a.clone());
        if !c.is_zero() {
            let g = ring.gcd(&ring.mul(&a, &c), &ring.mul(&b, &c));
            prop_assert!(ring.div_rem(&g, &c).1.is_zero() || g.is_zero());
        }
    }
}

#[test]
fn examples() {
    let ring = RationalPolys::new();
    let half = Rational::new(Integer::ONE, Integer::from(2));
    // x/2 + 1/3 is (3x + 2) / 6.
    let p = ring.from_fractions(&[Rational::new(Integer::ONE, Integer::from(3)), half.clone()]);
    assert_eq!(
        p.numerator().coeffs(),
        &[Integer::from(2), Integer::from(3)]
    );
    assert_eq!(p.denominator(), &Natural::from(6u64));
    // (x/2 + 1/3) 6 = 3x + 2, over 1.
    let six = ring.scale(&p, &Rational::from(6));
    assert_eq!(six, ring.from_i64s(&[2, 3]));
    assert!(six.denominator().is_one());
    // x^2 - 1/4 = (x - 1/2)(x + 1/2).
    let q = ring.from_fractions(&[
        Rational::new(Integer::from(-1), Integer::from(4)),
        Rational::from(0),
        Rational::from(1),
    ]);
    let (quot, rem) = ring.div_rem(
        &q,
        &ring.from_fractions(&[-half.clone(), Rational::from(1)]),
    );
    assert_eq!(ring.to_fractions(&quot), [half, Rational::from(1)]);
    assert!(rem.is_zero());
    // Zero is 0 / 1.
    assert_eq!(ring.sub(&p, &p), ring.zero());
    assert!(ring.zero().denominator().is_one());
    // Over ℚ(i): (x + i)^2 = x^2 + 2i x - 1.
    let ring = GaussianRationalPolys::new();
    let i = Complex::new(Rational::from(0), Rational::from(1));
    let p = ring.from_fractions(&[
        i.clone(),
        Complex::new(Rational::from(1), Rational::from(0)),
    ]);
    let two_i = Complex::new(Rational::from(0), Rational::from(2));
    assert_eq!(
        ring.to_fractions(&ring.square(&p)),
        [
            Complex::new(Rational::from(-1), Rational::from(0)),
            two_i,
            Complex::new(Rational::from(1), Rational::from(0))
        ]
    );
    // Monic over ℚ(i): (2i x + 1) / (2i) = x - i/2.
    let p = ring.from_fractions(&[
        Complex::new(Rational::from(1), Rational::from(0)),
        Complex::new(Rational::from(0), Rational::from(2)),
    ]);
    let (monic, _) = ring.normal_form(&p);
    assert_eq!(
        ring.to_fractions(&monic),
        [
            Complex::new(
                Rational::from(0),
                Rational::new(Integer::from(-1), Integer::from(2))
            ),
            Complex::new(Rational::from(1), Rational::from(0))
        ]
    );
}
