//! ℚ\[x\] and ℚ(i)\[x\] as integer, or Gaussian integer, polynomials over
//! one shared natural denominator (CAS-DOC-8 §11.2), as FLINT's
//! `fmpq_poly`: sums, products and divisions run on the numerators,
//! products by Kronecker substitution, and one gcd with the denominator
//! reduces the result.

use cas_num::algebra::{
    EuclideanDomain, Field, GaussianIntegers, GaussianRationals, Integers, Rationals, Ring,
};
use cas_num::{Complex, Integer, Natural, Rational};

use crate::cancel::{Cancel, Cancelled, Never};
use crate::dense::{Coefficients, Dense, DenseRing};
use crate::gcd::GcdDomain;

#[cfg(test)]
mod tests;

/// The numerator rings of shared-denominator polynomials, ℤ and ℤ\[i\]:
/// what reducing against a natural denominator takes, and the field of
/// fractions the polynomials are over.
pub trait Numerators: GcdDomain + Default {
    /// ℚ for ℤ, ℚ(i) for ℤ\[i\].
    type Fractions: Field + Default;

    /// `gcd(g, the integer parts of c)`.
    fn gcd_parts(g: &Natural, c: &Self::Elem) -> Natural;

    /// `c / n`, for an `n` that divides the integer parts of `c`.
    fn div_natural(c: &Self::Elem, n: &Natural) -> Self::Elem;

    /// `c · n`.
    fn mul_natural(c: &Self::Elem, n: &Natural) -> Self::Elem;

    /// `(u, n)` with `d u = n` a positive natural, for a non-zero `d`:
    /// `(sign d, |d|)` over ℤ, `(conj d, norm d)` over ℤ\[i\].
    fn to_natural(d: &Self::Elem) -> (Self::Elem, Natural);

    /// `c / den` in the field of fractions.
    fn fraction(c: &Self::Elem, den: &Natural) -> <Self::Fractions as Ring>::Elem;

    /// `(c, den)` with `x = c / den` and `den` positive.
    fn split(x: &<Self::Fractions as Ring>::Elem) -> (Self::Elem, Natural);
}

impl Numerators for Integers {
    type Fractions = Rationals;

    fn gcd_parts(g: &Natural, c: &Integer) -> Natural {
        g.gcd(&c.unsigned_abs())
    }

    fn div_natural(c: &Integer, n: &Natural) -> Integer {
        c.div_exact(&Integer::from(n))
    }

    fn mul_natural(c: &Integer, n: &Natural) -> Integer {
        c * &Integer::from(n)
    }

    fn to_natural(d: &Integer) -> (Integer, Natural) {
        let sign = if d.is_negative() {
            Integer::from(-1)
        } else {
            Integer::ONE
        };
        (sign, d.unsigned_abs())
    }

    fn fraction(c: &Integer, den: &Natural) -> Rational {
        Rational::new(c.clone(), Integer::from(den))
    }

    fn split(x: &Rational) -> (Integer, Natural) {
        (x.numer().clone(), x.denom().unsigned_abs())
    }
}

impl Numerators for GaussianIntegers {
    type Fractions = GaussianRationals;

    fn gcd_parts(g: &Natural, c: &Complex) -> Natural {
        g.gcd(&c.re().numer().unsigned_abs())
            .gcd(&c.im().numer().unsigned_abs())
    }

    fn div_natural(c: &Complex, n: &Natural) -> Complex {
        let n = Integer::from(n);
        Complex::new(
            Rational::from(c.re().numer().div_exact(&n)),
            Rational::from(c.im().numer().div_exact(&n)),
        )
    }

    fn mul_natural(c: &Complex, n: &Natural) -> Complex {
        let n = Rational::from(Integer::from(n));
        Complex::new(c.re() * &n, c.im() * &n)
    }

    fn to_natural(d: &Complex) -> (Complex, Natural) {
        (d.conj(), d.norm().numer().unsigned_abs())
    }

    fn fraction(c: &Complex, den: &Natural) -> Complex {
        let den = Rational::from(Integer::from(den));
        Complex::new(c.re() / &den, c.im() / &den)
    }

    fn split(x: &Complex) -> (Complex, Natural) {
        let den = x
            .re()
            .denom()
            .unsigned_abs()
            .lcm(&x.im().denom().unsigned_abs());
        let scale = Rational::from(Integer::from(&den));
        (Complex::new(x.re() * &scale, x.im() * &scale), den)
    }
}

/// A polynomial over the field of fractions of `E`'s ring: `num / den`, an
/// integer (or Gaussian integer) polynomial over a positive natural,
/// reduced: `den` and the integer parts of `num`'s coefficients are coprime,
/// and zero is `0 / 1`. Its operations are those of a [`SharedRing`].
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Shared<E> {
    num: Dense<E>,
    den: Natural,
}

/// A polynomial over ℚ.
pub type RationalPoly = Shared<Integer>;

/// A polynomial over ℚ(i).
pub type GaussianRationalPoly = Shared<Complex>;

impl<E> Shared<E> {
    /// The numerator, an integer (or Gaussian integer) polynomial.
    #[must_use]
    pub fn numerator(&self) -> &Dense<E> {
        &self.num
    }

    /// The positive denominator.
    #[must_use]
    pub fn denominator(&self) -> &Natural {
        &self.den
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    /// The degree; `None` for the zero polynomial.
    #[must_use]
    pub fn degree(&self) -> Option<usize> {
        self.num.degree()
    }
}

/// The ring of [`Shared`] polynomials over the field of fractions of `N`:
/// ℚ\[x\] for ℤ, ℚ(i)\[x\] for ℤ\[i\]. A [`Ring`] and a
/// [`EuclideanDomain`], with monic polynomials as the normal forms.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct SharedRing<N> {
    num: DenseRing<N>,
}

/// ℚ\[x\].
pub type RationalPolys = SharedRing<Integers>;

/// ℚ(i)\[x\].
pub type GaussianRationalPolys = SharedRing<GaussianIntegers>;

/// The elements of the field of fractions of `N`.
pub(crate) type Fraction<N> = <<N as Numerators>::Fractions as Ring>::Elem;

/// A quotient and a remainder.
type SharedQuoRem<E> = (Shared<E>, Shared<E>);

impl<N: Numerators> SharedRing<N> {
    #[must_use]
    pub fn new() -> Self {
        SharedRing {
            num: DenseRing::new(N::default()),
        }
    }

    /// The ring of the numerators.
    #[must_use]
    pub fn numerators(&self) -> &DenseRing<N> {
        &self.num
    }

    /// `num / den`, reduced.
    ///
    /// # Panics
    ///
    /// If `den` is zero.
    #[must_use]
    pub fn from_parts(&self, num: Dense<N::Elem>, den: Natural) -> Shared<N::Elem> {
        assert!(!den.is_zero(), "a zero denominator");
        if num.is_zero() {
            return self.zero();
        }
        let mut g = den.clone();
        for c in num.coeffs() {
            g = N::gcd_parts(&g, c);
            if g.is_one() {
                return Shared { num, den };
            }
        }
        let num = self
            .num
            .poly(num.coeffs().iter().map(|c| N::div_natural(c, &g)).collect());
        Shared {
            num,
            den: den.div_exact(&g),
        }
    }

    /// The polynomial of `coeffs` in the field of fractions, lowest degree
    /// first, over the least common denominator.
    #[must_use]
    pub fn from_fractions(&self, coeffs: &[Fraction<N>]) -> Shared<N::Elem> {
        let parts: Vec<(N::Elem, Natural)> = coeffs.iter().map(N::split).collect();
        let den = parts.iter().fold(Natural::one(), |l, (_, d)| l.lcm(d));
        let num = parts
            .iter()
            .map(|(c, d)| N::mul_natural(c, &den.div_exact(d)))
            .collect();
        self.from_parts(self.num.poly(num), den)
    }

    /// The polynomial of integer coefficients, lowest degree first.
    #[must_use]
    pub fn from_i64s(&self, coeffs: &[i64]) -> Shared<N::Elem> {
        self.from_parts(self.num.from_i64s(coeffs), Natural::one())
    }

    /// The coefficients in the field of fractions, lowest degree first.
    #[must_use]
    pub fn to_fractions(&self, p: &Shared<N::Elem>) -> Vec<Fraction<N>> {
        p.num
            .coeffs()
            .iter()
            .map(|c| N::fraction(c, &p.den))
            .collect()
    }

    /// `c p`.
    #[must_use]
    pub fn scale(&self, p: &Shared<N::Elem>, c: &Fraction<N>) -> Shared<N::Elem> {
        let (cn, cd) = N::split(c);
        self.from_parts(self.num.scale(&p.num, &cn), &p.den * &cd)
    }

    /// `p(x)`, by Horner's rule in the field of fractions.
    #[must_use]
    pub fn eval(&self, p: &Shared<N::Elem>, x: &Fraction<N>) -> Fraction<N> {
        let f = N::Fractions::default();
        self.to_fractions(p)
            .iter()
            .rev()
            .fold(f.zero(), |acc, c| f.add(&f.mul(&acc, x), c))
    }

    /// The derivative `p'`.
    #[must_use]
    pub fn derivative(&self, p: &Shared<N::Elem>) -> Shared<N::Elem> {
        self.from_parts(self.num.derivative(&p.num), p.den.clone())
    }

    /// `p · n` for a natural `n`, coefficient by coefficient.
    fn mul_natural(&self, p: &Dense<N::Elem>, n: &Natural) -> Dense<N::Elem> {
        if n.is_one() {
            return p.clone();
        }
        self.num
            .poly(p.coeffs().iter().map(|c| N::mul_natural(c, n)).collect())
    }

    /// `a b`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn mul_with(
        &self,
        a: &Shared<N::Elem>,
        b: &Shared<N::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Shared<N::Elem>, Cancelled> {
        let num = self.num.mul_with(&a.num, &b.num, cancel)?;
        Ok(self.from_parts(num, &a.den * &b.den))
    }

    /// `(q, r)` with `a = q b + r` and `deg r < deg b`, through the
    /// pseudo-division of the numerators: with `d` the leading coefficient
    /// of `b`'s numerator and `e = deg a - deg b + 1`, `d^e num(a) = Q
    /// num(b) + R`, so `q = Q den(b) / (d^e den(a))` and `r = R / (d^e
    /// den(a))`, `1 / d^e` made a natural's inverse by [`Numerators::to_natural`].
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[allow(
        clippy::many_single_char_names,
        reason = "the names of a = q b + r, d^e for the leading coefficient's power"
    )]
    pub fn div_rem_with(
        &self,
        a: &Shared<N::Elem>,
        b: &Shared<N::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<SharedQuoRem<N::Elem>, Cancelled> {
        let d = b.num.lead().expect("attempt to divide by zero");
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let (q, r) = self.num.pseudo_div_rem(&a.num, &b.num);
        let e = (a.num.coeffs().len() + 1).saturating_sub(b.num.coeffs().len());
        let (u, n) = N::to_natural(&self.num.base().pow(d, u64::try_from(e).expect("fits")));
        let den = &n * &a.den;
        let q = self.mul_natural(&self.num.scale(&q, &u), &b.den);
        let r = self.num.scale(&r, &u);
        Ok((self.from_parts(q, den.clone()), self.from_parts(r, den)))
    }
}

impl<N: Numerators> Ring for SharedRing<N> {
    type Elem = Shared<N::Elem>;

    fn zero(&self) -> Self::Elem {
        Shared {
            num: Dense::default(),
            den: Natural::one(),
        }
    }

    fn one(&self) -> Self::Elem {
        self.from_i64s(&[1])
    }

    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        if a.den == b.den {
            return self.from_parts(self.num.add(&a.num, &b.num), a.den.clone());
        }
        let g = a.den.gcd(&b.den);
        let (ca, cb) = (b.den.div_exact(&g), a.den.div_exact(&g));
        let num = self.num.add(
            &self.mul_natural(&a.num, &ca),
            &self.mul_natural(&b.num, &cb),
        );
        self.from_parts(num, &a.den * &ca)
    }

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.add(a, &self.neg(b))
    }

    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        Shared {
            num: self.num.neg(&a.num),
            den: a.den.clone(),
        }
    }

    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        match self.mul_with(a, b, &Never) {
            Ok(p) => p,
            Err(Cancelled) => unreachable!("never cancelled"),
        }
    }

    fn elem_integer(&self, n: &Integer) -> Self::Elem {
        self.from_parts(self.num.elem_integer(n), Natural::one())
    }

    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }

    fn is_zero(&self, a: &Self::Elem) -> bool {
        a.is_zero()
    }
}

impl<N: Numerators> EuclideanDomain for SharedRing<N> {
    fn div_rem(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem) {
        match self.div_rem_with(a, b, &Never) {
            Ok(qr) => qr,
            Err(Cancelled) => unreachable!("never cancelled"),
        }
    }

    /// The gcd of the numerators, [`DenseRing::gcd`], made monic.
    fn gcd(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        let g = self.num.gcd(&a.num, &b.num);
        self.normal_form(&self.from_parts(g, Natural::one())).0
    }

    /// The monic associate `a / lc(a)`, with `1 / lc(a)` as the unit.
    fn normal_form(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        let Some(lead) = a.num.lead() else {
            return (self.zero(), self.one());
        };
        let (u, n) = N::to_natural(lead);
        let monic = self.from_parts(self.num.scale(&a.num, &u), n.clone());
        let unit = self.from_parts(self.num.constant(N::mul_natural(&u, &a.den)), n);
        (monic, unit)
    }
}

/// Polynomials over ℚ\[x\] or ℚ(i)\[x\] multiply by Karatsuba's method.
impl<N: Numerators> Coefficients for SharedRing<N> {}
