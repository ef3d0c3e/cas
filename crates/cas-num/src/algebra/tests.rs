//! The ring laws for every instance, and the generic Euclidean algorithm on
//! ℤ (against GMP), on prime fields, and on a test-only polynomial ring
//! 𝔽ₚ[x] / ℚ[x] built on the traits.
#![allow(
    clippy::many_single_char_names,
    reason = "the laws use the textbook names (r, a, b, c, g, s, t, q)"
)]

use proptest::prelude::*;

use super::{
    BigPrimeField, EuclideanDomain, Field, GaussianIntegers, GaussianRationals, Integers,
    Rationals, Ring, WordPrimeField,
};
use crate::modular::{ModRing, ModWord};
use crate::testkit::{ToRug, strategy};
use crate::{Complex, Integer, Natural, Rational};

/// Fails unless `x == y`, printing both.
fn same<T: PartialEq + std::fmt::Debug>(x: &T, y: &T, law: &str) -> Result<(), TestCaseError> {
    prop_assert_eq!(x, y, "{}", law);
    Ok(())
}

/// The commutative ring axioms on `a, b, c`, plus the provided methods and
/// the canonical map from ℤ on `m, n`.
fn ring_laws<R: Ring>(
    r: &R,
    a: &R::Elem,
    b: &R::Elem,
    c: &R::Elem,
    m: &Integer,
    n: &Integer,
) -> Result<(), TestCaseError> {
    same(
        &r.add(&r.add(a, b), c),
        &r.add(a, &r.add(b, c)),
        "+ associative",
    )?;
    same(&r.add(a, b), &r.add(b, a), "+ commutative")?;
    same(&r.add(a, &r.zero()), a, "0 neutral")?;
    same(&r.add(a, &r.neg(a)), &r.zero(), "negation")?;
    same(&r.sub(a, b), &r.add(a, &r.neg(b)), "subtraction")?;
    same(
        &r.mul(&r.mul(a, b), c),
        &r.mul(a, &r.mul(b, c)),
        "* associative",
    )?;
    same(&r.mul(a, b), &r.mul(b, a), "* commutative")?;
    same(&r.mul(a, &r.one()), a, "1 neutral")?;
    same(
        &r.mul(a, &r.add(b, c)),
        &r.add(&r.mul(a, b), &r.mul(a, c)),
        "distributive",
    )?;
    same(&r.square(a), &r.mul(a, a), "square")?;
    prop_assert_eq!(r.is_zero(a), *a == r.zero());
    prop_assert_eq!(r.is_one(a), *a == r.one());
    let cube = r.mul(&r.square(a), a);
    same(&r.pow(a, 3), &cube, "pow 3")?;
    same(&r.pow(a, 0), &r.one(), "pow 0")?;
    same(&r.pow(a, 13), &r.mul(&r.pow(a, 6), &r.pow(a, 7)), "pow 13")?;
    same(
        &r.elem_integer(&(m + n)),
        &r.add(&r.elem_integer(m), &r.elem_integer(n)),
        "hom +",
    )?;
    same(
        &r.elem_integer(&(m * n)),
        &r.mul(&r.elem_integer(m), &r.elem_integer(n)),
        "hom *",
    )?;
    same(&r.elem_integer(&Integer::ONE), &r.one(), "hom 1")?;
    let ch = Integer::from(r.characteristic());
    same(&r.elem_integer(&ch), &r.zero(), "characteristic")?;
    for x in [a, b, c] {
        if let Some(k) = r.small_preimage(x) {
            same(&r.elem_integer(&Integer::from(k)), x, "small preimage")?;
        }
    }
    Ok(())
}

/// The field laws: `a / b * b = a` and inverses.
fn field_laws<F: Field>(f: &F, a: &F::Elem, b: &F::Elem) -> Result<(), TestCaseError> {
    match f.inv(b) {
        None => prop_assert!(f.is_zero(b), "non-zero element without inverse"),
        Some(i) => {
            same(&f.mul(b, &i), &f.one(), "inverse")?;
            same(&f.mul(&f.div(a, b).unwrap(), b), a, "division")?;
            let (q, r) = f.div_rem(a, b);
            prop_assert!(f.is_zero(&r));
            same(&f.mul(&q, b), a, "div_rem")?;
            same(&f.normal_form(b).0, &f.one(), "normal form")?;
        }
    }
    prop_assert!(f.div(a, &f.zero()).is_none());
    same(
        &f.normal_form(&f.zero()),
        &(f.zero(), f.one()),
        "normal form of 0",
    )?;
    Ok(())
}

/// `gcd` and `gcd_ext` on `a, b`: Bézout, the GCD divides both, it is in
/// normal form, and `gcd` agrees with `gcd_ext`.
fn gcd_laws<R: EuclideanDomain>(r: &R, a: &R::Elem, b: &R::Elem) -> Result<(), TestCaseError> {
    let (g, s, t) = r.gcd_ext(a, b);
    same(&r.add(&r.mul(&s, a), &r.mul(&t, b)), &g, "Bézout")?;
    same(&r.gcd(a, b), &g, "gcd = gcd_ext")?;
    same(&r.normal_form(&g).0, &g, "normal form")?;
    if r.is_zero(&g) {
        prop_assert!(r.is_zero(a) && r.is_zero(b));
    } else {
        prop_assert!(r.is_zero(&r.div_rem(a, &g).1), "gcd divides a");
        prop_assert!(r.is_zero(&r.div_rem(b, &g).1), "gcd divides b");
    }
    Ok(())
}

/// ℤ with the trait's own Euclidean algorithm (no Lehmer override).
struct PlainIntegers;

impl Ring for PlainIntegers {
    type Elem = Integer;
    fn zero(&self) -> Integer {
        Integer::ZERO
    }
    fn one(&self) -> Integer {
        Integer::ONE
    }
    fn add(&self, a: &Integer, b: &Integer) -> Integer {
        a + b
    }
    fn sub(&self, a: &Integer, b: &Integer) -> Integer {
        a - b
    }
    fn neg(&self, a: &Integer) -> Integer {
        -a
    }
    fn mul(&self, a: &Integer, b: &Integer) -> Integer {
        a * b
    }
    fn elem_integer(&self, n: &Integer) -> Integer {
        n.clone()
    }
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl EuclideanDomain for PlainIntegers {
    fn div_rem(&self, a: &Integer, b: &Integer) -> (Integer, Integer) {
        Integers.div_rem(a, b)
    }
    fn normal_form(&self, a: &Integer) -> (Integer, Integer) {
        Integers.normal_form(a)
    }
}

/// Dense polynomials over a field, lowest coefficient first, without
/// trailing zeros: enough of 𝔽[x] to show that polynomial rings fit the
/// traits.
struct Poly<F>(F);

impl<F: Field> Poly<F> {
    fn trim(&self, mut v: Vec<F::Elem>) -> Vec<F::Elem> {
        while v.last().is_some_and(|c| self.0.is_zero(c)) {
            v.pop();
        }
        v
    }

    fn poly(&self, cs: &[i64]) -> Vec<F::Elem> {
        self.trim(
            cs.iter()
                .map(|&c| self.0.elem_integer(&Integer::from(c)))
                .collect(),
        )
    }
}

impl<F: Field> Ring for Poly<F> {
    type Elem = Vec<F::Elem>;

    fn zero(&self) -> Self::Elem {
        Vec::new()
    }

    fn one(&self) -> Self::Elem {
        self.trim(vec![self.0.one()])
    }

    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        let n = a.len().max(b.len());
        let z = self.0.zero();
        self.trim(
            (0..n)
                .map(|i| self.0.add(a.get(i).unwrap_or(&z), b.get(i).unwrap_or(&z)))
                .collect(),
        )
    }

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.add(a, &self.neg(b))
    }

    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        a.iter().map(|c| self.0.neg(c)).collect()
    }

    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut r = vec![self.0.zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                r[i + j] = self.0.add(&r[i + j], &self.0.mul(x, y));
            }
        }
        self.trim(r)
    }

    fn elem_integer(&self, n: &Integer) -> Self::Elem {
        self.trim(vec![self.0.elem_integer(n)])
    }

    fn characteristic(&self) -> Natural {
        self.0.characteristic()
    }
}

impl<F: Field> EuclideanDomain for Poly<F> {
    /// Long division by the leading coefficient's inverse.
    fn div_rem(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem) {
        let lead = self
            .0
            .inv(b.last().expect("attempt to divide by zero"))
            .unwrap();
        let mut r = a.clone();
        let mut q = vec![self.0.zero(); a.len().saturating_sub(b.len() - 1)];
        while r.len() >= b.len() {
            let shift = r.len() - b.len();
            let c = self.0.mul(r.last().unwrap(), &lead);
            for (i, y) in b.iter().enumerate() {
                r[shift + i] = self.0.sub(&r[shift + i], &self.0.mul(&c, y));
            }
            q[shift] = c;
            r = self.trim(r);
        }
        (self.trim(q), r)
    }

    /// Monic.
    fn normal_form(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        match a.last() {
            None => (Vec::new(), self.one()),
            Some(lead) => {
                let u = self.0.inv(lead).unwrap();
                (a.iter().map(|c| self.0.mul(c, &u)).collect(), vec![u])
            }
        }
    }
}

const WORD_PRIMES: [u64; 6] = [2, 3, 7, 65521, (1 << 61) - 1, u64::MAX - 58];

fn word_prime() -> impl Strategy<Value = u64> {
    prop::sample::select(&WORD_PRIMES[..])
}

fn big_prime() -> impl Strategy<Value = Natural> {
    prop::sample::select(vec![89u64, 127, 521]).prop_map(|k| Natural::power_of_two(k) - 1u64)
}

fn rational() -> impl Strategy<Value = Rational> {
    strategy::rational().prop_map(|(n, d)| Rational::new(n.to_integer(), d.to_integer()))
}

/// Complex rationals (or Gaussian integers) with parts of at most two limbs,
/// so that `pow(13)` in the ring laws stays small.
fn small_complex(integral: bool) -> impl Strategy<Value = Complex> {
    let part = move || {
        (strategy::signed_up_to(2), strategy::signed_up_to(1)).prop_filter_map(
            "zero denominator",
            move |(n, d)| {
                let d = if integral {
                    Integer::ONE
                } else {
                    d.to_integer()
                };
                Rational::checked_new(n.to_integer(), d)
            },
        )
    };
    (part(), part()).prop_map(|(re, im)| Complex::new(re, im))
}

fn small_poly() -> impl Strategy<Value = Vec<i64>> {
    prop::collection::vec(-50i64..50, 0..6)
}

proptest! {
    #[test]
    fn integers(a in strategy::integer(), b in strategy::integer(), c in strategy::integer()) {
        let (a, b, c) = (a.to_integer(), b.to_integer(), c.to_integer());
        ring_laws(&Integers, &a, &b, &c, &a, &b)?;
        gcd_laws(&Integers, &a, &b)?;
        gcd_laws(&PlainIntegers, &a, &b)?;
        if !b.is_zero() {
            let (q, r) = Integers.div_rem(&a, &b);
            prop_assert!(!r.is_negative() && r < b.abs() && q * &b + &r == a, "Euclidean division");
        }
        prop_assert_eq!(Integers.normal_form(&Integer::ZERO), (Integer::ZERO, Integer::ONE));
        // The provided Euclidean algorithm agrees with GMP.
        prop_assert_eq!(PlainIntegers.gcd(&a, &b).to_rug(), a.to_rug().gcd(&b.to_rug()));
    }

    #[test]
    fn gaussian_rationals(a in small_complex(false), b in small_complex(false), c in small_complex(false), m in strategy::integer(), n in strategy::integer()) {
        let (m, n) = (m.to_integer(), n.to_integer());
        ring_laws(&GaussianRationals, &a, &b, &c, &m, &n)?;
        field_laws(&GaussianRationals, &a, &b)?;
        gcd_laws(&GaussianRationals, &a, &b)?;
    }

    #[test]
    fn gaussian_integers(a in small_complex(true), b in small_complex(true), c in small_complex(true), m in strategy::integer(), n in strategy::integer()) {
        let (m, n) = (m.to_integer(), n.to_integer());
        ring_laws(&GaussianIntegers, &a, &b, &c, &m, &n)?;
        gcd_laws(&GaussianIntegers, &a, &b)?;
    }

    #[test]
    fn rationals(a in rational(), b in rational(), c in rational(), m in strategy::integer(), n in strategy::integer()) {
        let (m, n) = (m.to_integer(), n.to_integer());
        ring_laws(&Rationals, &a, &b, &c, &m, &n)?;
        field_laws(&Rationals, &a, &b)?;
        gcd_laws(&Rationals, &a, &b)?;
    }

    #[test]
    fn mod_word(n in strategy::top_limb(), a in strategy::integer(), b in strategy::integer(), c in strategy::integer()) {
        let r = ModWord::new(n);
        let (a, b, c) = (a.to_integer(), b.to_integer(), c.to_integer());
        let (x, y, z) = (r.elem_integer(&a), r.elem_integer(&b), r.elem_integer(&c));
        ring_laws(&r, &x, &y, &z, &a, &b)?;
    }

    #[test]
    fn mod_ring(n in strategy::natural_up_to(4).prop_filter("non-zero", |n| !n.is_zero()), a in strategy::integer(), b in strategy::integer(), c in strategy::integer()) {
        let r = ModRing::new(n.to_natural());
        let (a, b, c) = (a.to_integer(), b.to_integer(), c.to_integer());
        let (x, y, z) = (r.elem_integer(&a), r.elem_integer(&b), r.elem_integer(&c));
        ring_laws(&r, &x, &y, &z, &a, &b)?;
    }

    #[test]
    fn word_prime_field(p in word_prime(), a in strategy::integer(), b in strategy::integer(), c in strategy::integer()) {
        let f = WordPrimeField::new(p).unwrap();
        let (a, b, c) = (a.to_integer(), b.to_integer(), c.to_integer());
        let (x, y, z) = (f.elem_integer(&a), f.elem_integer(&b), f.elem_integer(&c));
        ring_laws(&f, &x, &y, &z, &a, &b)?;
        field_laws(&f, &x, &y)?;
        gcd_laws(&f, &x, &y)?;
        prop_assert_eq!(f.ring().modulus(), p);
    }

    #[test]
    fn big_prime_field(p in big_prime(), a in strategy::integer(), b in strategy::integer(), c in strategy::integer()) {
        let f = BigPrimeField::new(p.clone()).unwrap();
        let (a, b, c) = (a.to_integer(), b.to_integer(), c.to_integer());
        let (x, y, z) = (f.elem_integer(&a), f.elem_integer(&b), f.elem_integer(&c));
        ring_laws(&f, &x, &y, &z, &a, &b)?;
        field_laws(&f, &x, &y)?;
        gcd_laws(&f, &x, &y)?;
        prop_assert_eq!(f.characteristic(), p);
    }

    /// 𝔽ₚ[x]: the GCD of `g u` and `g v` is a multiple of `g`, and Bézout
    /// holds. Small primes make common factors beyond `g` frequent.
    #[test]
    fn polynomials_over_prime_field(p in word_prime(), g in small_poly(), u in small_poly(), v in small_poly(), c in small_poly()) {
        let ring = Poly(WordPrimeField::new(p).unwrap());
        let (g, u, v, c) = (ring.poly(&g), ring.poly(&u), ring.poly(&v), ring.poly(&c));
        let (a, b) = (ring.mul(&g, &u), ring.mul(&g, &v));
        ring_laws(&ring, &a, &b, &c, &Integer::from(3), &Integer::from(-5))?;
        gcd_laws(&ring, &a, &b)?;
        let d = ring.gcd(&a, &b);
        if !ring.is_zero(&d) {
            prop_assert!(ring.is_zero(&ring.div_rem(&d, &ring.normal_form(&g).0).1) || ring.is_zero(&g));
        }
    }

    /// ℚ[x] over the same inputs.
    #[test]
    fn polynomials_over_rationals(g in small_poly(), u in small_poly(), v in small_poly()) {
        let ring = Poly(Rationals);
        let (g, u, v) = (ring.poly(&g), ring.poly(&u), ring.poly(&v));
        let (a, b) = (ring.mul(&g, &u), ring.mul(&g, &v));
        gcd_laws(&ring, &a, &b)?;
        let d = ring.gcd(&a, &b);
        if !ring.is_zero(&g) {
            prop_assert!(ring.is_zero(&ring.div_rem(&d, &ring.normal_form(&g).0).1));
        }
    }
}

#[test]
fn known_polynomial_gcd() {
    // Over 𝔽_7: x^2 - 1 = (x - 1)(x + 1) and x^2 + 2x + 1 = (x + 1)^2.
    let ring = Poly(WordPrimeField::new(7).unwrap());
    let a = ring.poly(&[-1, 0, 1]);
    let b = ring.poly(&[1, 2, 1]);
    assert_eq!(ring.gcd(&a, &b), ring.poly(&[1, 1]));
    // Over 𝔽_2, x^2 + 1 = (x + 1)^2 is a multiple of x + 1.
    let ring = Poly(WordPrimeField::new(2).unwrap());
    let a = ring.poly(&[1, 0, 1]);
    assert_eq!(ring.gcd(&a, &ring.poly(&[1, 1])), ring.poly(&[1, 1]));
}

#[test]
fn prime_field_construction() {
    assert!(WordPrimeField::new(1).is_none());
    assert!(WordPrimeField::new(91).is_none());
    assert!(WordPrimeField::new(u64::MAX).is_none());
    assert!(BigPrimeField::new(Natural::power_of_two(128) + 1u64).is_none());
    assert!(BigPrimeField::new(Natural::power_of_two(127) - 1u64).is_some());
    // A small prime also works over the multi-precision ring.
    let f = BigPrimeField::new(Natural::from(5u64)).unwrap();
    let three = f.elem_integer(&Integer::from(3));
    assert_eq!(
        f.ring().to_natural(&f.inv(&three).unwrap()),
        Natural::from(2u64)
    );
}

/// Which elements have a small preimage: integers in `i64`, integral
/// rationals likewise, word residues below `2^63`; none over the rings
/// without one at hand.
#[test]
fn small_preimages() {
    let big = Integer::from(i64::MAX) + Integer::ONE;
    assert_eq!(Integers.small_preimage(&Integer::from(-5)), Some(-5));
    assert_eq!(
        Integers.small_preimage(&Integer::from(i64::MIN)),
        Some(i64::MIN)
    );
    assert_eq!(Integers.small_preimage(&big), None);
    assert_eq!(Rationals.small_preimage(&Rational::from(7)), Some(7));
    let half = Rational::new(Integer::ONE, Integer::from(2));
    assert_eq!(Rationals.small_preimage(&half), None);
    assert_eq!(Rationals.small_preimage(&Rational::from(big.clone())), None);
    let r = ModWord::new(u64::MAX);
    assert_eq!(r.small_preimage(&r.elem_integer(&Integer::from(-1))), None);
    assert_eq!(
        r.small_preimage(&r.elem_integer(&Integer::from(9))),
        Some(9)
    );
    let f = WordPrimeField::new(u64::MAX - 58).unwrap();
    assert_eq!(f.small_preimage(&f.elem_integer(&Integer::from(-1))), None);
    assert_eq!(
        f.small_preimage(&f.elem_integer(&Integer::from(i64::MAX))),
        Some(i64::MAX)
    );
    let g = BigPrimeField::new(Natural::from(5u64)).unwrap();
    assert_eq!(g.small_preimage(&g.elem_integer(&Integer::from(3))), None);
    let i = Complex::new(Rational::ZERO, Rational::ONE);
    assert_eq!(GaussianIntegers.small_preimage(&i), None);
}

#[test]
fn characteristics() {
    assert_eq!(Integers.characteristic(), Natural::ZERO);
    assert_eq!(Rationals.characteristic(), Natural::ZERO);
    assert_eq!(ModWord::new(12).characteristic(), Natural::from(12u64));
    assert_eq!(Ring::one(&ModWord::new(1)), 0);
}
