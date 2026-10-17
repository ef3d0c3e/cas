//! Algebraic structures as ring objects.
//!
//! A ring is a value: the traits are implemented by a type describing the
//! ring (`R: Ring`), its elements are `R::Elem`, and the operations are
//! methods on `&R`. Rings fixed at run time, such as ℤ/nℤ for a modulus read
//! from input or a polynomial ring over it, are then ordinary values, and
//! generic algorithms take the ring as an argument:
//!
//! ```
//! use cas_num::algebra::{EuclideanDomain, Integers, WordPrimeField};
//! use cas_num::Integer;
//!
//! fn gcd_of_all<R: EuclideanDomain>(ring: &R, xs: &[R::Elem]) -> R::Elem {
//!     xs.iter().fold(ring.zero(), |g, x| ring.gcd(&g, x))
//! }
//!
//! let zz = Integers;
//! let xs = [Integer::from(12), Integer::from(-18), Integer::from(30)];
//! assert_eq!(gcd_of_all(&zz, &xs), Integer::from(6));
//!
//! let f7 = WordPrimeField::new(7).unwrap();
//! assert_eq!(gcd_of_all(&f7, &[3, 5]), 1);
//! ```
//!
//! Instances: [`Integers`] (ℤ), [`Rationals`] (ℚ), [`ModWord`] and
//! [`ModRing`] (ℤ/nℤ, rings only), and [`PrimeField`] over either (𝔽ₚ).

use crate::modular::{ModRing, ModWord, Residue};
use crate::{Integer, Natural, Rational};

#[cfg(test)]
mod tests;

/// A commutative ring with identity.
pub trait Ring {
    /// The elements; equality is equality in the ring.
    type Elem: Clone + PartialEq + std::fmt::Debug;

    fn zero(&self) -> Self::Elem;

    fn one(&self) -> Self::Elem;

    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    fn neg(&self, a: &Self::Elem) -> Self::Elem;

    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    /// The image of `n` under the canonical map ℤ → R.
    fn elem_integer(&self, n: &Integer) -> Self::Elem;

    /// The characteristic: 0 for ℤ and ℚ, `n` for ℤ/nℤ.
    fn characteristic(&self) -> Natural;

    fn is_zero(&self, a: &Self::Elem) -> bool {
        *a == self.zero()
    }

    fn is_one(&self, a: &Self::Elem) -> bool {
        *a == self.one()
    }

    fn square(&self, a: &Self::Elem) -> Self::Elem {
        self.mul(a, a)
    }

    /// `a^e`, with `a^0 = 1`, by binary exponentiation.
    fn pow(&self, a: &Self::Elem, e: u64) -> Self::Elem {
        let mut r = self.one();
        for i in (0..u64::BITS - e.leading_zeros()).rev() {
            r = self.square(&r);
            if (e >> i) & 1 == 1 {
                r = self.mul(&r, a);
            }
        }
        r
    }
}

/// An integral domain with division with remainder.
///
/// The provided `gcd` and `gcd_ext` run the Euclidean algorithm on
/// [`EuclideanDomain::div_rem`] and normalise the result with
/// [`EuclideanDomain::normal_form`]; instances with a faster algorithm
/// override them.
pub trait EuclideanDomain: Ring {
    /// `(q, r)` with `a = q b + r` and `r` smaller than `b` in the domain's
    /// Euclidean function (so the Euclidean algorithm terminates).
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    fn div_rem(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem);

    /// `(n, u)` with `n = a u` for a unit `u`, where `n` is the canonical
    /// representative of the associates of `a`: `|a|` in ℤ, 1 for a non-zero
    /// element of a field, a monic polynomial. Zero gives `(0, 1)`.
    fn normal_form(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem);

    /// The greatest common divisor, in normal form; `gcd(0, 0) = 0`.
    fn gcd(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        let (mut a, mut b) = (a.clone(), b.clone());
        while !self.is_zero(&b) {
            let r = self.div_rem(&a, &b).1;
            (a, b) = (b, r);
        }
        self.normal_form(&a).0
    }

    /// `(g, s, t)` with `g = gcd(a, b) = s a + t b`, `g` in normal form.
    fn gcd_ext(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem, Self::Elem) {
        // Invariants: r0 = s0 a + t0 b and r1 = s1 a + t1 b.
        let (mut r0, mut r1) = (a.clone(), b.clone());
        let (mut s0, mut s1) = (self.one(), self.zero());
        let (mut t0, mut t1) = (self.zero(), self.one());
        while !self.is_zero(&r1) {
            let (q, r) = self.div_rem(&r0, &r1);
            let s = self.sub(&s0, &self.mul(&q, &s1));
            let t = self.sub(&t0, &self.mul(&q, &t1));
            (r0, r1) = (r1, r);
            (s0, s1) = (s1, s);
            (t0, t1) = (t1, t);
        }
        let (g, u) = self.normal_form(&r0);
        (g, self.mul(&s0, &u), self.mul(&t0, &u))
    }
}

/// A field.
pub trait Field: EuclideanDomain {
    /// `1 / a`, or `None` for zero.
    fn inv(&self, a: &Self::Elem) -> Option<Self::Elem>;

    /// `a / b`, or `None` if `b` is zero.
    fn div(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.mul(a, &self.inv(b)?))
    }
}

/// The ring ℤ, with [`Integer`] elements.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct Integers;

/// The field ℚ, with [`Rational`] elements.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct Rationals;

impl Ring for Integers {
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

    fn is_zero(&self, a: &Integer) -> bool {
        a.is_zero()
    }

    fn pow(&self, a: &Integer, e: u64) -> Integer {
        a.pow(e)
    }
}

impl EuclideanDomain for Integers {
    /// The Euclidean division: `0 <= r < |b|`.
    fn div_rem(&self, a: &Integer, b: &Integer) -> (Integer, Integer) {
        a.div_rem_euclid(b)
    }

    fn normal_form(&self, a: &Integer) -> (Integer, Integer) {
        if a.is_negative() {
            (-a, Integer::from(-1))
        } else {
            (a.clone(), Integer::ONE)
        }
    }

    /// Lehmer's algorithm ([`Integer::gcd`]).
    fn gcd(&self, a: &Integer, b: &Integer) -> Integer {
        a.gcd(b)
    }

    /// [`Integer::gcd_ext`], whose cofactors are GMP's.
    fn gcd_ext(&self, a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
        a.gcd_ext(b)
    }
}

impl Ring for Rationals {
    type Elem = Rational;

    fn zero(&self) -> Rational {
        Rational::ZERO
    }

    fn one(&self) -> Rational {
        Rational::ONE
    }

    fn add(&self, a: &Rational, b: &Rational) -> Rational {
        a + b
    }

    fn sub(&self, a: &Rational, b: &Rational) -> Rational {
        a - b
    }

    fn neg(&self, a: &Rational) -> Rational {
        -a
    }

    fn mul(&self, a: &Rational, b: &Rational) -> Rational {
        a * b
    }

    fn elem_integer(&self, n: &Integer) -> Rational {
        Rational::from(n)
    }

    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }

    fn is_zero(&self, a: &Rational) -> bool {
        a.is_zero()
    }
}

impl Field for Rationals {
    fn inv(&self, a: &Rational) -> Option<Rational> {
        a.checked_recip()
    }
}

/// Euclidean structure of a field: exact division, zero remainder, every
/// non-zero element associate to 1.
macro_rules! field_euclidean {
    ($t:ty) => {
        impl EuclideanDomain for $t {
            fn div_rem(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem) {
                let q = self.div(a, b).expect("attempt to divide by zero");
                (q, self.zero())
            }

            fn normal_form(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
                match self.inv(a) {
                    Some(u) => (self.one(), u),
                    None => (self.zero(), self.one()),
                }
            }
        }
    };
}

field_euclidean!(Rationals);

impl Ring for ModWord {
    type Elem = u64;

    fn zero(&self) -> u64 {
        0
    }

    fn one(&self) -> u64 {
        self.reduce(1)
    }

    fn add(&self, a: &u64, b: &u64) -> u64 {
        ModWord::add(self, *a, *b)
    }

    fn sub(&self, a: &u64, b: &u64) -> u64 {
        ModWord::sub(self, *a, *b)
    }

    fn neg(&self, a: &u64) -> u64 {
        ModWord::neg(self, *a)
    }

    fn mul(&self, a: &u64, b: &u64) -> u64 {
        ModWord::mul(self, *a, *b)
    }

    fn elem_integer(&self, n: &Integer) -> u64 {
        self.reduce_integer(n)
    }

    fn characteristic(&self) -> Natural {
        Natural::from(self.modulus())
    }

    fn pow(&self, a: &u64, e: u64) -> u64 {
        ModWord::pow(self, *a, e)
    }
}

impl Ring for ModRing {
    type Elem = Residue;

    fn zero(&self) -> Residue {
        ModRing::zero(self)
    }

    fn one(&self) -> Residue {
        ModRing::one(self)
    }

    fn add(&self, a: &Residue, b: &Residue) -> Residue {
        ModRing::add(self, a, b)
    }

    fn sub(&self, a: &Residue, b: &Residue) -> Residue {
        ModRing::sub(self, a, b)
    }

    fn neg(&self, a: &Residue) -> Residue {
        ModRing::neg(self, a)
    }

    fn mul(&self, a: &Residue, b: &Residue) -> Residue {
        ModRing::mul(self, a, b)
    }

    fn elem_integer(&self, n: &Integer) -> Residue {
        ModRing::elem_integer(self, n)
    }

    fn characteristic(&self) -> Natural {
        self.modulus().clone()
    }

    fn is_zero(&self, a: &Residue) -> bool {
        ModRing::is_zero(self, a)
    }

    fn square(&self, a: &Residue) -> Residue {
        ModRing::square(self, a)
    }

    fn pow(&self, a: &Residue, e: u64) -> Residue {
        ModRing::pow(self, a, &Natural::from(e))
    }
}

/// The prime field 𝔽ₚ = ℤ/pℤ over a [`ModWord`] or a [`ModRing`], with
/// the same elements. Construction checks that `p` is prime.
#[derive(Clone, Debug)]
pub struct PrimeField<M> {
    ring: M,
}

/// 𝔽ₚ for a prime that fits a limb, with `u64` elements.
pub type WordPrimeField = PrimeField<ModWord>;

/// 𝔽ₚ for a prime of any size, with [`Residue`] elements.
pub type BigPrimeField = PrimeField<ModRing>;

impl<M> PrimeField<M> {
    /// The underlying ring ℤ/pℤ, for conversions.
    pub fn ring(&self) -> &M {
        &self.ring
    }
}

impl PrimeField<ModWord> {
    /// 𝔽ₚ, or `None` unless `p` is prime (proven, as every word is).
    #[must_use]
    pub fn new(p: u64) -> Option<Self> {
        crate::prime::is_prime_u64(p).then(|| PrimeField {
            ring: ModWord::new(p),
        })
    }
}

impl PrimeField<ModRing> {
    /// 𝔽ₚ, or `None` if `p` is not prime. Above 2^64 a probable prime (in
    /// the sense of [`Natural::primality`]) is accepted.
    #[must_use]
    pub fn new(p: Natural) -> Option<Self> {
        p.primality().is_probable_prime().then(|| PrimeField {
            ring: ModRing::new(p),
        })
    }
}

impl<M: Ring> Ring for PrimeField<M> {
    type Elem = M::Elem;

    fn zero(&self) -> M::Elem {
        self.ring.zero()
    }

    fn one(&self) -> M::Elem {
        self.ring.one()
    }

    fn add(&self, a: &M::Elem, b: &M::Elem) -> M::Elem {
        self.ring.add(a, b)
    }

    fn sub(&self, a: &M::Elem, b: &M::Elem) -> M::Elem {
        self.ring.sub(a, b)
    }

    fn neg(&self, a: &M::Elem) -> M::Elem {
        self.ring.neg(a)
    }

    fn mul(&self, a: &M::Elem, b: &M::Elem) -> M::Elem {
        self.ring.mul(a, b)
    }

    fn elem_integer(&self, n: &Integer) -> M::Elem {
        self.ring.elem_integer(n)
    }

    fn characteristic(&self) -> Natural {
        self.ring.characteristic()
    }

    fn is_zero(&self, a: &M::Elem) -> bool {
        self.ring.is_zero(a)
    }

    fn square(&self, a: &M::Elem) -> M::Elem {
        self.ring.square(a)
    }

    fn pow(&self, a: &M::Elem, e: u64) -> M::Elem {
        self.ring.pow(a, e)
    }
}

impl Field for PrimeField<ModWord> {
    fn inv(&self, a: &u64) -> Option<u64> {
        // p is prime, so only 0 has no inverse.
        self.ring.inv(*a)
    }
}

impl Field for PrimeField<ModRing> {
    fn inv(&self, a: &Residue) -> Option<Residue> {
        self.ring.inv(a)
    }
}

field_euclidean!(PrimeField<ModWord>);
field_euclidean!(PrimeField<ModRing>);
