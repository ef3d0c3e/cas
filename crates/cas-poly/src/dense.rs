//! Dense univariate polynomials (CAS-DOC-8 §11.2): coefficients lowest
//! degree first, without trailing zeros. `doc/algorithms/dense.md`
//! describes the algorithms.

mod div;
mod kronecker;
mod mul;

#[cfg(test)]
mod tests;

pub use mul::Coefficients;

use cas_num::algebra::{EuclideanDomain, Field, Ring};
use cas_num::{Integer, Natural};

use crate::cancel::{Cancel, Cancelled, Never, uncancelled};

/// A polynomial `c₀ + c₁x + … + cₙxⁿ` with `cₙ ≠ 0`: its coefficients
/// lowest degree first, without trailing zeros, so that the zero polynomial
/// has none and equal polynomials are equal values. The operations are
/// those of a [`DenseRing`] over the ring of the coefficients.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Dense<E> {
    coeffs: Vec<E>,
}

impl<E> Default for Dense<E> {
    fn default() -> Self {
        Dense { coeffs: Vec::new() }
    }
}

impl<E> Dense<E> {
    /// The coefficients, lowest degree first.
    #[must_use]
    pub fn coeffs(&self) -> &[E] {
        &self.coeffs
    }

    #[must_use]
    pub fn into_coeffs(self) -> Vec<E> {
        self.coeffs
    }

    /// Whether this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The degree; `None` for the zero polynomial.
    #[must_use]
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// The leading coefficient; `None` for the zero polynomial.
    #[must_use]
    pub fn lead(&self) -> Option<&E> {
        self.coeffs.last()
    }

    /// The coefficient of `x^k`; `None` past the degree, where it is 0.
    #[must_use]
    pub fn coeff(&self, k: usize) -> Option<&E> {
        self.coeffs.get(k)
    }
}

/// The ring `R[x]` of dense polynomials over the ring object `R`, itself a
/// ring object: a [`Ring`] whose elements are [`Dense`] polynomials, and a
/// [`EuclideanDomain`] when `R` is a field, with monic polynomials as the
/// normal forms.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct DenseRing<R> {
    base: R,
}

impl<R: Coefficients> DenseRing<R> {
    #[must_use]
    pub fn new(base: R) -> Self {
        DenseRing { base }
    }

    /// The ring of the coefficients.
    #[must_use]
    pub fn base(&self) -> &R {
        &self.base
    }

    /// The polynomial of `coeffs`, lowest degree first: trailing zeros are
    /// dropped.
    #[must_use]
    pub fn poly(&self, mut coeffs: Vec<R::Elem>) -> Dense<R::Elem> {
        while coeffs.last().is_some_and(|c| self.base.is_zero(c)) {
            coeffs.pop();
        }
        Dense { coeffs }
    }

    /// The polynomial of integer coefficients, lowest degree first, through
    /// the canonical map ℤ → R.
    #[must_use]
    pub fn from_i64s(&self, coeffs: &[i64]) -> Dense<R::Elem> {
        self.poly(
            coeffs
                .iter()
                .map(|&c| self.base.elem_integer(&Integer::from(c)))
                .collect(),
        )
    }

    /// The constant `c`.
    #[must_use]
    pub fn constant(&self, c: R::Elem) -> Dense<R::Elem> {
        self.poly(vec![c])
    }

    /// `c x^k`.
    #[must_use]
    pub fn monomial(&self, c: R::Elem, k: usize) -> Dense<R::Elem> {
        if self.base.is_zero(&c) {
            return Dense::default();
        }
        let mut coeffs = vec![self.base.zero(); k];
        coeffs.push(c);
        Dense { coeffs }
    }

    /// The variable `x`.
    #[must_use]
    pub fn x(&self) -> Dense<R::Elem> {
        self.monomial(self.base.one(), 1)
    }

    /// `c p`.
    #[must_use]
    pub fn scale(&self, p: &Dense<R::Elem>, c: &R::Elem) -> Dense<R::Elem> {
        self.poly(p.coeffs.iter().map(|a| self.base.mul(c, a)).collect())
    }

    /// `p x^k`.
    #[must_use]
    pub fn shift(&self, p: &Dense<R::Elem>, k: usize) -> Dense<R::Elem> {
        if p.is_zero() {
            return Dense::default();
        }
        let mut coeffs = vec![self.base.zero(); k];
        coeffs.extend(p.coeffs.iter().cloned());
        Dense { coeffs }
    }

    /// `p(x)`, by Horner's rule.
    #[must_use]
    pub fn eval(&self, p: &Dense<R::Elem>, x: &R::Elem) -> R::Elem {
        p.coeffs.iter().rev().fold(self.base.zero(), |acc, c| {
            self.base.add(&self.base.mul(&acc, x), c)
        })
    }

    /// The derivative `p'`.
    #[must_use]
    pub fn derivative(&self, p: &Dense<R::Elem>) -> Dense<R::Elem> {
        self.poly(
            p.coeffs
                .iter()
                .enumerate()
                .skip(1)
                .map(|(k, c)| self.base.mul(&self.base.elem_integer(&Integer::from(k)), c))
                .collect(),
        )
    }

    /// `a b`, polling `cancel` in the long loops.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn mul_with(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Dense<R::Elem>, Cancelled> {
        if a.is_zero() || b.is_zero() {
            return Ok(Dense::default());
        }
        Ok(self.poly(self.base.mul_coeffs(&a.coeffs, &b.coeffs, cancel)?))
    }

    /// `a²`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn square_with(
        &self,
        a: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Dense<R::Elem>, Cancelled> {
        if a.is_zero() {
            return Ok(Dense::default());
        }
        Ok(self.poly(self.base.square_coeffs(&a.coeffs, cancel)?))
    }

    /// `a^n`, with `a^0 = 1`, by binary powering, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn pow_with(
        &self,
        a: &Dense<R::Elem>,
        n: u64,
        cancel: &dyn Cancel,
    ) -> Result<Dense<R::Elem>, Cancelled> {
        let mut r = self.one();
        for i in (0..u64::BITS - n.leading_zeros()).rev() {
            r = self.square_with(&r, cancel)?;
            if (n >> i) & 1 == 1 {
                r = self.mul_with(&r, a, cancel)?;
            }
        }
        Ok(r)
    }

    /// The coefficient-wise combination of `a` and `b` by `f`, a missing
    /// coefficient being 0.
    fn zip(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        f: impl Fn(&R::Elem, &R::Elem) -> R::Elem,
    ) -> Dense<R::Elem> {
        let zero = self.base.zero();
        let n = a.coeffs.len().max(b.coeffs.len());
        self.poly(
            (0..n)
                .map(|i| {
                    f(
                        a.coeffs.get(i).unwrap_or(&zero),
                        b.coeffs.get(i).unwrap_or(&zero),
                    )
                })
                .collect(),
        )
    }
}

impl<R: Coefficients> Ring for DenseRing<R> {
    type Elem = Dense<R::Elem>;

    fn zero(&self) -> Self::Elem {
        Dense::default()
    }

    fn one(&self) -> Self::Elem {
        self.constant(self.base.one())
    }

    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.zip(a, b, |x, y| self.base.add(x, y))
    }

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.zip(a, b, |x, y| self.base.sub(x, y))
    }

    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.poly(a.coeffs.iter().map(|c| self.base.neg(c)).collect())
    }

    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        uncancelled(self.mul_with(a, b, &Never))
    }

    fn elem_integer(&self, n: &Integer) -> Self::Elem {
        self.constant(self.base.elem_integer(n))
    }

    fn characteristic(&self) -> Natural {
        self.base.characteristic()
    }

    fn is_zero(&self, a: &Self::Elem) -> bool {
        a.is_zero()
    }

    fn square(&self, a: &Self::Elem) -> Self::Elem {
        uncancelled(self.square_with(a, &Never))
    }

    fn pow(&self, a: &Self::Elem, e: u64) -> Self::Elem {
        uncancelled(self.pow_with(a, e, &Never))
    }
}

/// Polynomials over polynomials multiply by Karatsuba's method.
impl<R: Coefficients> Coefficients for DenseRing<R> {}

impl<F: Coefficients + Field> EuclideanDomain for DenseRing<F> {
    /// Long division: `a = q b + r` with `deg r < deg b`.
    fn div_rem(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem) {
        uncancelled(self.div_rem_with(a, b, &Never))
    }

    /// The monic associate, with the inverse of the leading coefficient as
    /// the unit.
    fn normal_form(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        let Some(lead) = a.lead() else {
            return (self.zero(), self.one());
        };
        let inv = self.base.inv(lead).expect("a non-zero leading coefficient");
        (self.scale(a, &inv), self.constant(inv))
    }
}

impl<R: Coefficients + EuclideanDomain> DenseRing<R> {
    /// The content: the gcd of the coefficients, in the normal form of `R`
    /// (non-negative over ℤ); 0 for the zero polynomial.
    #[must_use]
    pub fn content(&self, p: &Dense<R::Elem>) -> R::Elem {
        let mut g = self.base.zero();
        for c in &p.coeffs {
            g = self.base.gcd(&g, c);
            if self.base.is_one(&g) {
                break;
            }
        }
        g
    }

    /// The primitive part: `p` divided by its content, then by the unit
    /// that puts its leading coefficient in normal form, so that over ℤ it
    /// has a positive leading coefficient. `p` is a unit times
    /// `content(p) · primitive_part(p)`. Zero for zero.
    #[must_use]
    pub fn primitive_part(&self, p: &Dense<R::Elem>) -> Dense<R::Elem> {
        let Some(lead) = p.lead() else {
            return Dense::default();
        };
        let g = self.content(p);
        let unit = self.base.normal_form(lead).1;
        self.poly(
            p.coeffs
                .iter()
                .map(|c| {
                    let (q, r) = self.base.div_rem(c, &g);
                    debug_assert!(self.base.is_zero(&r), "the content divides");
                    self.base.mul(&q, &unit)
                })
                .collect(),
        )
    }
}
