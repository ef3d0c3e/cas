//! Tests of dense polynomials: ring laws, division, and every product
//! algorithm against schoolbook, over ℤ, ℚ, ℤ\[i\], ℚ(i), 𝔽ₚ and ℤ/nℤ.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities"
)]

use std::cell::Cell;
use std::sync::atomic::AtomicBool;

use cas_num::algebra::{
    BigPrimeField, EuclideanDomain, GaussianIntegers, GaussianRationals, Integers, Rationals, Ring,
    WordPrimeField,
};
use cas_num::ball::RealBalls;
use cas_num::modular::ModWord;
use cas_num::{Complex, Integer, Natural, Rational};
use proptest::prelude::*;

use super::{Coefficients, Dense, DenseRing, kronecker, mul};
use crate::cancel::{Cancelled, Never};

/// An integer of up to `limbs` limbs, of either sign, often small.
fn integer(limbs: usize) -> impl Strategy<Value = Integer> {
    prop_oneof![
        (-20i64..20).prop_map(Integer::from),
        (
            any::<bool>(),
            prop::collection::vec(any::<u64>(), 1..=limbs)
        )
            .prop_map(|(negative, l)| Integer::from_parts(negative, Natural::from_vec(l))),
    ]
}

/// Integer coefficients, `len` of them at most, of up to `limbs` limbs.
fn coeffs(len: usize, limbs: usize) -> impl Strategy<Value = Vec<Integer>> {
    prop::collection::vec(integer(limbs), 0..=len)
}

/// Primes that fit a word: small ones, where products often vanish, and
/// the largest.
fn word_prime() -> impl Strategy<Value = u64> {
    prop::sample::select(vec![2, 3, 7, 101, 65_537, (1 << 61) - 1, u64::MAX - 58])
}

fn zx() -> DenseRing<Integers> {
    DenseRing::new(Integers)
}

/// The image of integer coefficients in `R`.
fn image<R: Coefficients>(ring: &DenseRing<R>, cs: &[Integer]) -> Dense<R::Elem> {
    ring.poly(cs.iter().map(|c| ring.base().elem_integer(c)).collect())
}

/// Gaussian integers from pairs of integer coefficients.
fn gaussian(re: &[Integer], im: &[Integer]) -> Vec<Complex> {
    let n = re.len().max(im.len());
    let at = |v: &[Integer], i: usize| Rational::from(v.get(i).cloned().unwrap_or(Integer::ZERO));
    (0..n).map(|i| Complex::new(at(re, i), at(im, i))).collect()
}

/// The laws of a commutative ring with identity on `x`, `y`, `z`.
fn ring_laws<R: Ring>(r: &R, x: &R::Elem, y: &R::Elem, z: &R::Elem) -> Result<(), TestCaseError> {
    prop_assert_eq!(r.add(x, y), r.add(y, x));
    prop_assert_eq!(r.mul(x, y), r.mul(y, x));
    prop_assert_eq!(r.add(&r.add(x, y), z), r.add(x, &r.add(y, z)));
    prop_assert_eq!(r.mul(&r.mul(x, y), z), r.mul(x, &r.mul(y, z)));
    prop_assert_eq!(r.mul(x, &r.add(y, z)), r.add(&r.mul(x, y), &r.mul(x, z)));
    prop_assert_eq!(r.add(x, &r.zero()), x.clone());
    prop_assert_eq!(r.mul(x, &r.one()), x.clone());
    prop_assert!(r.is_zero(&r.add(x, &r.neg(x))));
    prop_assert_eq!(r.sub(x, y), r.add(x, &r.neg(y)));
    prop_assert_eq!(r.square(x), r.mul(x, x));
    prop_assert_eq!(r.pow(x, 3), r.mul(x, &r.mul(x, x)));
    Ok(())
}

/// Division over a field: `a = q b + r` with `deg r < deg b`, and `(a b) / b
/// = a` exactly; pseudo-division and exact division agree; the content of a
/// non-zero polynomial is 1 and its primitive part monic.
fn division_laws<F: Coefficients + cas_num::algebra::Field>(
    ring: &DenseRing<F>,
    a: &Dense<F::Elem>,
    b: &Dense<F::Elem>,
) -> Result<(), TestCaseError> {
    if b.is_zero() {
        return Ok(());
    }
    let f = ring.base();
    let (q, r) = ring.pseudo_div_rem(a, b);
    let e = (a.coeffs().len() + 1).saturating_sub(b.coeffs().len());
    let power = f.pow(b.lead().unwrap(), u64::try_from(e).unwrap());
    prop_assert_eq!(ring.scale(a, &power), ring.add(&ring.mul(&q, b), &r));
    prop_assert_eq!(ring.div_exact(&ring.mul(a, b), b), Some(a.clone()));
    if let Some(q) = ring.div_exact(a, b) {
        prop_assert_eq!(ring.mul(&q, b), a.clone());
    }
    prop_assert!(f.is_one(&ring.content(b)));
    prop_assert_eq!(ring.primitive_part(b), ring.normal_form(b).0);
    let (q, r) = ring.div_rem(a, b);
    prop_assert_eq!(ring.add(&ring.mul(&q, b), &r), a.clone());
    prop_assert!(r.degree() < b.degree());
    let (q, r) = ring.div_rem(&ring.mul(a, b), b);
    prop_assert_eq!(q, a.clone());
    prop_assert!(r.is_zero());
    let (monic, unit) = ring.normal_form(b);
    prop_assert_eq!(ring.mul(b, &unit), monic.clone());
    prop_assert!(ring.base().is_one(monic.lead().unwrap()));
    Ok(())
}

/// Horner's rule and the derivative agree with products: `(a b)(x) =
/// a(x) b(x)` and `(a b)' = a' b + a b'`.
fn eval_laws<R: Coefficients>(
    ring: &DenseRing<R>,
    a: &Dense<R::Elem>,
    b: &Dense<R::Elem>,
    x: &R::Elem,
) -> Result<(), TestCaseError> {
    let r = ring.base();
    let ab = ring.mul(a, b);
    prop_assert_eq!(ring.eval(&ab, x), r.mul(&ring.eval(a, x), &ring.eval(b, x)));
    prop_assert_eq!(
        ring.derivative(&ab),
        ring.add(
            &ring.mul(&ring.derivative(a), b),
            &ring.mul(a, &ring.derivative(b))
        )
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn integers_form_a_ring(x in coeffs(20, 2), y in coeffs(20, 2), z in coeffs(12, 2)) {
        let ring = zx();
        let (x, y, z) = (ring.poly(x), ring.poly(y), ring.poly(z));
        ring_laws(&ring, &x, &y, &z)?;
        eval_laws(&ring, &x, &y, &Integer::from(-3))?;
    }

    /// Kronecker substitution against schoolbook, with negative coefficients
    /// of up to 256 bits, across the threshold.
    #[test]
    fn kronecker_is_schoolbook(a in coeffs(40, 4), b in coeffs(40, 4)) {
        prop_assume!(!a.is_empty() && !b.is_empty());
        let expected = mul::schoolbook(&Integers, &a, &b, &Never).unwrap();
        prop_assert_eq!(kronecker::mul(&a, &b), expected);
        prop_assert_eq!(kronecker::square(&a), mul::schoolbook(&Integers, &a, &a, &Never).unwrap());
    }

    /// Karatsuba against schoolbook, over ℚ and a big prime field.
    #[test]
    fn karatsuba_is_schoolbook(a in coeffs(70, 1), b in coeffs(70, 1)) {
        prop_assume!(!a.is_empty() && !b.is_empty());
        let q = |v: &[Integer]| -> Vec<Rational> {
            v.iter().map(|c| Rational::new(c.clone(), Integer::from(7))).collect()
        };
        let (qa, qb) = (q(&a), q(&b));
        prop_assert_eq!(
            mul::karatsuba(&Rationals, &qa, &qb, &Never).unwrap(),
            mul::schoolbook(&Rationals, &qa, &qb, &Never).unwrap()
        );
        // The Mersenne prime 2^127 - 1, of two limbs.
        let f = BigPrimeField::new(Natural::power_of_two(127) - 1).unwrap();
        let (fa, fb): (Vec<_>, Vec<_>) = (
            a.iter().map(|c| f.elem_integer(c)).collect(),
            b.iter().map(|c| f.elem_integer(c)).collect(),
        );
        prop_assert_eq!(
            mul::karatsuba(&f, &fa, &fb, &Never).unwrap(),
            mul::schoolbook(&f, &fa, &fb, &Never).unwrap()
        );
    }

    /// Pseudo-division and exact division over ℤ.
    #[test]
    fn integer_division(a in coeffs(16, 2), b in coeffs(8, 2)) {
        let ring = zx();
        let (a, b) = (ring.poly(a), ring.poly(b));
        prop_assume!(!b.is_zero());
        let (q, r) = ring.pseudo_div_rem(&a, &b);
        let e = (a.coeffs().len() + 1).saturating_sub(b.coeffs().len());
        let power = b.lead().unwrap().pow(u64::try_from(e).unwrap());
        prop_assert_eq!(
            ring.scale(&a, &power),
            ring.add(&ring.mul(&q, &b), &r)
        );
        prop_assert!(r.degree() < b.degree());
        prop_assert_eq!(ring.div_exact(&ring.mul(&a, &b), &b), Some(a.clone()));
        // Not a multiple: a b + 1 over b of degree 1 or more.
        if b.degree() >= Some(1) {
            let off = ring.add(&ring.mul(&a, &b), &ring.one());
            prop_assert_eq!(ring.div_exact(&off, &b), None);
        }
        // Whatever it says divides does.
        if let Some(q) = ring.div_exact(&a, &b) {
            prop_assert_eq!(ring.mul(&q, &b), a.clone());
        }
    }

    /// The content and primitive part over ℤ: `p = ± content · pp`, with a
    /// primitive `pp` of positive leading coefficient.
    #[test]
    fn content_and_primitive_part(p in coeffs(12, 2), k in integer(1)) {
        let ring = zx();
        let p = ring.scale(&ring.poly(p), &k);
        let (c, pp) = (ring.content(&p), ring.primitive_part(&p));
        if p.is_zero() {
            prop_assert!(c.is_zero() && pp.is_zero());
            return Ok(());
        }
        prop_assert!(c.is_positive());
        prop_assert!(pp.lead().unwrap().is_positive());
        prop_assert!(ring.content(&pp).is_one());
        let back = ring.scale(&pp, &c);
        prop_assert!(back == p || ring.neg(&back) == p);
        prop_assert!(Integer::from(&k.unsigned_abs()).is_zero() || c.is_divisible_by(&k.abs()));
    }

    /// 𝔽ₚ for word primes: the ring laws, Kronecker products against
    /// schoolbook, division and evaluation.
    #[test]
    fn word_prime_fields(p in word_prime(), a in coeffs(80, 2), b in coeffs(80, 2), c in coeffs(10, 1)) {
        let ring = DenseRing::new(WordPrimeField::new(p).unwrap());
        let (a, b, c) = (image(&ring, &a), image(&ring, &b), image(&ring, &c));
        ring_laws(&ring, &a, &b, &c)?;
        division_laws(&ring, &a, &b)?;
        division_laws(&ring, &b, &c)?;
        eval_laws(&ring, &a, &b, &ring.base().elem_integer(&Integer::from(5)))?;
        if !a.is_zero() && !b.is_zero() {
            let f = ring.base();
            prop_assert_eq!(
                f.mul_coeffs(a.coeffs(), b.coeffs(), &Never).unwrap(),
                mul::schoolbook(f, a.coeffs(), b.coeffs(), &Never).unwrap()
            );
        }
        // The Euclidean gcd: a common factor stays.
        if !c.is_zero() {
            let g = ring.gcd(&ring.mul(&a, &c), &ring.mul(&b, &c));
            if !g.is_zero() {
                prop_assert!(ring.div_rem(&g, &ring.normal_form(&c).0).1.is_zero());
            }
        }
    }

    /// ℤ/nℤ for a word modulus, zero divisors included: Kronecker products
    /// against schoolbook.
    #[test]
    fn word_moduli(n in 2u64..1000, a in coeffs(60, 1), b in coeffs(60, 1)) {
        let m = ModWord::new(n);
        let (a, b): (Vec<u64>, Vec<u64>) = (
            a.iter().map(|c| m.elem_integer(c)).collect(),
            b.iter().map(|c| m.elem_integer(c)).collect(),
        );
        prop_assume!(!a.is_empty() && !b.is_empty());
        prop_assert_eq!(
            m.mul_coeffs(&a, &b, &Never).unwrap(),
            mul::schoolbook(&m, &a, &b, &Never).unwrap()
        );
    }

    /// ℚ: the ring laws and division.
    #[test]
    fn rationals(a in coeffs(12, 1), b in coeffs(8, 1), c in coeffs(6, 1), d in 1i64..50) {
        let ring = DenseRing::new(Rationals);
        let q = |v: &[Integer]| ring.poly(v.iter().map(|c| Rational::new(c.clone(), Integer::from(d))).collect());
        let (a, b, c) = (q(&a), q(&b), q(&c));
        ring_laws(&ring, &a, &b, &c)?;
        division_laws(&ring, &a, &b)?;
        eval_laws(&ring, &a, &b, &Rational::new(Integer::from(2), Integer::from(3)))?;
    }

    /// ℤ\[i\]: the ring laws, and Kronecker products by parts against
    /// schoolbook, past the threshold.
    #[test]
    fn gaussian_integers(re in coeffs(30, 3), im in coeffs(30, 3), re2 in coeffs(30, 3), im2 in coeffs(30, 3)) {
        let ring = DenseRing::new(GaussianIntegers);
        let (a, b) = (gaussian(&re, &im), gaussian(&re2, &im2));
        if !a.is_empty() && !b.is_empty() {
            prop_assert_eq!(
                GaussianIntegers.mul_coeffs(&a, &b, &Never).unwrap(),
                mul::schoolbook(&GaussianIntegers, &a, &b, &Never).unwrap()
            );
        }
        let (a, b) = (ring.poly(a), ring.poly(b));
        ring_laws(&ring, &a, &b, &ring.from_i64s(&[1, -2]))?;
        let i = Complex::new(Rational::from(0), Rational::from(1));
        eval_laws(&ring, &a, &b, &i)?;
        if !b.is_zero() {
            prop_assert_eq!(ring.div_exact(&ring.mul(&a, &b), &b), Some(a.clone()));
            let (q, r) = ring.pseudo_div_rem(&a, &b);
            let e = (a.coeffs().len() + 1).saturating_sub(b.coeffs().len());
            let power = GaussianIntegers.pow(b.lead().unwrap(), u64::try_from(e).unwrap());
            prop_assert_eq!(ring.scale(&a, &power), ring.add(&ring.mul(&q, &b), &r));
        }
    }

    /// ℚ(i): the ring laws and division.
    #[test]
    fn gaussian_rationals(re in coeffs(10, 1), im in coeffs(10, 1), re2 in coeffs(6, 1), im2 in coeffs(6, 1)) {
        let ring = DenseRing::new(GaussianRationals);
        let half = |v: Vec<Complex>| -> Vec<Complex> {
            let two = Rational::from(2);
            v.into_iter().map(|c| Complex::new(c.re() / &two, c.im().clone())).collect()
        };
        let (a, b) = (ring.poly(half(gaussian(&re, &im))), ring.poly(gaussian(&re2, &im2)));
        ring_laws(&ring, &a, &b, &ring.from_i64s(&[3, 0, 1]))?;
        division_laws(&ring, &a, &b)?;
        eval_laws(&ring, &a, &b, &Complex::new(Rational::from(1), Rational::new(Integer::from(1), Integer::from(2))))?;
    }
}

/// Every coefficient of the product at its largest: factors of the most
/// magnitude their bits allow, all of one sign, as many as fill the bits
/// counted for the sum. A digit a bit short wraps.
#[test]
fn kronecker_at_the_bounds() {
    for bits in [1, 2, 63, 64, 65, 200] {
        let top = Integer::from(Natural::power_of_two(bits)) - Integer::ONE;
        for terms in [1, 2, 3, 4, 7, 8, 15, 16, 31, 33] {
            for sign in [1i64, -1] {
                let a = vec![&top * &Integer::from(sign); terms];
                let b = vec![top.clone(); terms];
                let expected = mul::schoolbook(&Integers, &a, &b, &Never).unwrap();
                assert_eq!(
                    kronecker::mul(&a, &b),
                    expected,
                    "{bits} bits, {terms} terms"
                );
                assert_eq!(
                    kronecker::square(&a),
                    mul::schoolbook(&Integers, &a, &a, &Never).unwrap()
                );
            }
        }
    }
    // The largest word prime, every residue p - 1.
    let p = u64::MAX - 58;
    let m = ModWord::new(p);
    for terms in [33, 64, 100, 257] {
        let a = vec![p - 1; terms];
        assert_eq!(
            kronecker::mul_mod(&a, &a, &m),
            mul::schoolbook(&m, &a, &a, &Never).unwrap(),
            "{terms} terms"
        );
    }
}

#[test]
fn exact_division_checks_every_coefficient() {
    let ring = zx();
    // 3x / 2x: the leading coefficient does not divide.
    assert_eq!(
        ring.div_exact(&ring.from_i64s(&[0, 3]), &ring.from_i64s(&[0, 2])),
        None
    );
    assert_eq!(
        ring.div_exact(&ring.from_i64s(&[1, 3, 2]), &ring.from_i64s(&[1, 2])),
        Some(ring.from_i64s(&[1, 1]))
    );
    assert_eq!(
        ring.div_exact(&ring.from_i64s(&[1, 3]), &ring.from_i64s(&[1, 2, 1])),
        None
    );
    assert_eq!(
        ring.div_exact(&ring.zero(), &ring.from_i64s(&[5])),
        Some(ring.zero())
    );
}

#[test]
fn small_examples() {
    let ring = zx();
    // (x + 1)^5.
    let p = ring.pow(&ring.from_i64s(&[1, 1]), 5);
    assert_eq!(p, ring.from_i64s(&[1, 5, 10, 10, 5, 1]));
    assert_eq!(ring.derivative(&p), ring.from_i64s(&[5, 20, 30, 20, 5]));
    assert_eq!(ring.eval(&p, &Integer::from(1)), Integer::from(32));
    assert_eq!(p.degree(), Some(5));
    assert_eq!(ring.zero().degree(), None);
    assert_eq!(
        ring.monomial(Integer::from(3), 2),
        ring.from_i64s(&[0, 0, 3])
    );
    assert_eq!(ring.shift(&ring.x(), 2), ring.from_i64s(&[0, 0, 0, 1]));
    // Trailing zeros go.
    assert_eq!(ring.from_i64s(&[1, 0, 0]), ring.one());
    // Content and primitive part: -6x^2 - 4 = -2 (3x^2 + 2).
    let p = ring.from_i64s(&[-4, 0, -6]);
    assert_eq!(ring.content(&p), Integer::from(2));
    assert_eq!(ring.primitive_part(&p), ring.from_i64s(&[2, 0, 3]));
    // Pseudo-division: 2^2 (x^2 + 1) = (2x - 1)(2x + 1) + 5.
    let (q, r) = ring.pseudo_div_rem(&ring.from_i64s(&[1, 0, 1]), &ring.from_i64s(&[-1, 2]));
    assert_eq!((q, r), (ring.from_i64s(&[1, 2]), ring.from_i64s(&[5])));
    // Over 𝔽_7, x^2 - 1 = (x + 1)(x - 1).
    let f7 = DenseRing::new(WordPrimeField::new(7).unwrap());
    let (q, r) = f7.div_rem(&f7.from_i64s(&[-1, 0, 1]), &f7.from_i64s(&[1, 1]));
    assert_eq!((q, r), (f7.from_i64s(&[-1, 1]), f7.zero()));
    // Over ℤ/4ℤ, (2x + 1)(2x + 1) = 4x^2 + 4x + 1 = 1: the leading product
    // vanishes.
    let z4 = DenseRing::new(ModWord::new(4));
    let p = z4.from_i64s(&[1, 2]);
    assert_eq!(z4.mul(&p, &p), z4.one());
}

/// A ring that counts its products, to bound an algorithm's work.
struct Counting<R> {
    ring: R,
    products: Cell<u64>,
}

impl<R: Ring> Ring for Counting<R> {
    type Elem = R::Elem;

    fn zero(&self) -> R::Elem {
        self.ring.zero()
    }

    fn one(&self) -> R::Elem {
        self.ring.one()
    }

    fn add(&self, a: &R::Elem, b: &R::Elem) -> R::Elem {
        self.ring.add(a, b)
    }

    fn sub(&self, a: &R::Elem, b: &R::Elem) -> R::Elem {
        self.ring.sub(a, b)
    }

    fn neg(&self, a: &R::Elem) -> R::Elem {
        self.ring.neg(a)
    }

    fn mul(&self, a: &R::Elem, b: &R::Elem) -> R::Elem {
        self.products.set(self.products.get() + 1);
        self.ring.mul(a, b)
    }

    fn elem_integer(&self, n: &Integer) -> R::Elem {
        self.ring.elem_integer(n)
    }

    fn characteristic(&self) -> Natural {
        self.ring.characteristic()
    }
}

impl<R: Ring> Coefficients for Counting<R> {}

/// Pseudo-division by a divisor of degree `n` takes `O(s n)` products for
/// `s = deg a - n + 1` rounds, and none by a leading coefficient of 1,
/// where Algorithm R multiplied the whole remainder and quotient by it at
/// each round: `s^2`, 4 million, for `x^2000 + …` by `x + 1` (CAS-172). The
/// results are Algorithm R's.
#[test]
fn pseudo_division_scales_lazily() {
    let counting = DenseRing::new(Counting {
        ring: Integers,
        products: Cell::new(0),
    });
    let zx = zx();
    let mut a = vec![0; 2001];
    (a[0], a[7], a[1999], a[2000]) = (-1, 3, 5, 1);
    for (b, monic) in [
        (&[1, 1][..], true),
        (&[1, 0, 1][..], true),
        (&[1, -1][..], false),
        (&[3, 1, 2][..], false),
        (&[-4, 0, 0, 0, 0, 0, 0, 0, 0, 7][..], false),
        (&[2][..], false),
    ] {
        let (a, b) = (counting.from_i64s(&a), counting.from_i64s(b));
        counting.base().products.set(0);
        let (q, r) = counting.pseudo_div_rem(&a, &b);
        let products = counting.base().products.get();
        let n = b.degree().unwrap() as u64;
        let s = 2000 - n + 1;
        let bound = if monic { s * n } else { s * (2 * n + 3) };
        assert!(products <= bound, "{products} products by {b:?}");
        let power = b.lead().unwrap().pow(s);
        assert_eq!(
            zx.scale(&a, &power),
            zx.add(&zx.mul(&q, &b), &r),
            "by {b:?}"
        );
        assert!(r.degree() < b.degree());
    }
}

#[test]
fn polynomials_over_polynomials() {
    // ℤ[x][y]: (y + x)(y - x) = y^2 - x^2.
    let zx = zx();
    let zxy = DenseRing::new(zx.clone());
    let x = zx.x();
    let a = zxy.poly(vec![x.clone(), zx.one()]);
    let b = zxy.poly(vec![zx.neg(&x), zx.one()]);
    let p = zxy.mul(&a, &b);
    assert_eq!(
        p,
        zxy.poly(vec![zx.neg(&zx.mul(&x, &x)), zx.zero(), zx.one()])
    );
}

#[test]
fn balls_contain_the_product() {
    // Exact integer balls multiply exactly.
    let ring = DenseRing::new(RealBalls::new(64));
    let a = ring.from_i64s(&[1, 2, 3]);
    let b = ring.from_i64s(&[-1, 5]);
    assert_eq!(ring.mul(&a, &b), ring.from_i64s(&[-1, 3, 7, 15]));
}

#[test]
fn products_stop_when_cancelled() {
    let ring = zx();
    let a = ring.poly((0..200).map(Integer::from).collect());
    let stop = AtomicBool::new(true);
    assert_eq!(ring.mul_with(&a, &a, &stop), Err(Cancelled));
    let f = DenseRing::new(Rationals);
    let q = f.poly((0..100).map(Rational::from).collect());
    assert_eq!(f.mul_with(&q, &q, &stop), Err(Cancelled));
    assert_eq!(f.div_rem_with(&q, &f.x(), &stop), Err(Cancelled));
    let go = AtomicBool::new(false);
    assert_eq!(ring.mul_with(&a, &a, &go), Ok(ring.mul(&a, &a)));
}

/// A workload for instruction counts (CAS-135), run by `perf stat` on the
/// test binary with `--ignored --exact`: `CAS_POLY_WORK` is
/// `algorithm:length:bits:repeats`, the algorithm one of `schoolbook`,
/// `karatsuba`, `kronecker`, `word-schoolbook`, `word-kronecker` or `mul`
/// (the product `DenseRing` picks), the factors of `length` coefficients
/// of `bits` bits, of either sign over ℤ, and modulo a prime of about
/// `bits` bits for the word algorithms. Without it, nothing: the cost of
/// the harness alone.
#[test]
#[ignore = "a workload for instruction counts"]
fn products_workload() {
    let Ok(spec) = std::env::var("CAS_POLY_WORK") else {
        return;
    };
    let parts: Vec<&str> = spec.split(':').collect();
    let [algorithm, length, bits, repeats] = parts[..] else {
        panic!("CAS_POLY_WORK=algorithm:length:bits:repeats");
    };
    let (length, bits, repeats): (usize, u64, usize) = (
        length.parse().unwrap(),
        bits.parse().unwrap(),
        repeats.parse().unwrap(),
    );
    // xorshift64*: the same factors in every run.
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state.wrapping_mul(0x2545_f491_4f6c_dd1d)
    };
    let mut factor = || -> Vec<Integer> {
        (0..length)
            .map(|_| {
                let limbs: Vec<u64> = (0..bits.div_ceil(64)).map(|_| next()).collect();
                let magnitude = Natural::from_vec(limbs) >> (bits.div_ceil(64) * 64 - bits);
                Integer::from_parts(next() & 1 == 1, magnitude)
            })
            .collect()
    };
    let (a, b) = (factor(), factor());
    // The word algorithms work modulo a prime of about `bits` bits.
    let p = match bits {
        ..=17 => 65_537,
        18..=31 => (1 << 31) - 1,
        32..=61 => (1 << 61) - 1,
        _ => u64::MAX - 58,
    };
    let m = ModWord::new(p);
    let (wa, wb): (Vec<u64>, Vec<u64>) = (
        a.iter().map(|c| m.elem_integer(c)).collect(),
        b.iter().map(|c| m.elem_integer(c)).collect(),
    );
    let ring = zx();
    let (pa, pb) = (ring.poly(a.clone()), ring.poly(b.clone()));
    for _ in 0..repeats {
        let n = match algorithm {
            "schoolbook" => mul::schoolbook(&Integers, &a, &b, &Never).unwrap().len(),
            "karatsuba" => mul::karatsuba(&Integers, &a, &b, &Never).unwrap().len(),
            "kronecker" => kronecker::mul(&a, &b).len(),
            "word-schoolbook" => mul::schoolbook(&m, &wa, &wb, &Never).unwrap().len(),
            "word-kronecker" => kronecker::mul_mod(&wa, &wb, &m).len(),
            "mul" => ring.mul(&pa, &pb).coeffs().len(),
            _ => panic!("unknown algorithm {algorithm}"),
        };
        std::hint::black_box(n);
    }
}
