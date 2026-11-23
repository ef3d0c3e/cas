//! Greatest common divisors, resultants and discriminants (CAS-DOC-8
//! §11.3), over the coefficient domains that have them here: the
//! [`GcdDomain`]s ℤ, ℚ, ℤ\[i\], ℚ(i) and the prime fields.
//! `doc/algorithms/gcd.md` describes the algorithms.
//!
//! - A gcd is normalised: its leading coefficient, the lexicographic one
//!   for sparse polynomials, is in the normal form of the coefficients, and
//!   over ℤ and ℤ\[i\] the gcd of the contents is a factor of it. Over ℤ it
//!   has a positive leading coefficient (CAS-DOC-8 §15, detail 5); over a
//!   field it is monic.
//! - `gcd` chooses the algorithm; `gcd_by` runs one of the [`GcdMethod`]s,
//!   for tests and benchmarks.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the identities: a = g x, b = g y"
)]

mod ezgcd;
mod gaussian;
mod heuristic;
mod modular;
mod resultant;
mod subres;
mod vandermonde;
mod zippel;

#[cfg(test)]
mod tests;

use cas_num::algebra::{
    EuclideanDomain, Field, GaussianIntegers, GaussianRationals, Integers, PrimeField, Rationals,
    Ring,
};
use cas_num::modular::{ModRing, ModWord};
use cas_num::{Complex, Integer, Natural};

use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::dense::{Coefficients, Dense, DenseRing};
use crate::shared::{Fraction, Numerators, Shared, SharedRing};
use crate::sparse::{Sparse, SparseRing};

/// A gcd and its cofactors: `(g, a / g, b / g)`.
pub type Cofactors<P> = (P, P, P);

/// `(g, u, v, c)` of an extended gcd, `u a + v b = c g`
/// ([`SparseRing::extended_gcd`]).
pub type Bezout<P> = (P, P, P, P);

/// `(u, v, d)` of `u a + v b = d c` ([`SparseRing::abcuv`]).
pub type Abcuv<P> = (P, P, P);

/// The gcd algorithms. `gcd` chooses among them; `gcd_by` runs one.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum GcdMethod {
    /// Evaluation at a large integer, the gcd of the values, and the
    /// polynomial read back from its digits, checked by division \[CGG89,
    /// Par05\]. Over ℤ and ℤ\[i\], and ℚ and ℚ(i) through their numerators.
    /// It can give up.
    Heuristic,
    /// Brown's algorithm \[Bro71\]: over ℤ, images modulo word primes,
    /// combined by Chinese remaindering; over ℤ\[i\], sparse, the same
    /// modulo word primes `p ≡ 1 (mod 4)` under both maps `i ↦ ±√−1`
    /// (CAS-494); over a prime field, images at points of the last
    /// variable, combined by interpolation. Over ℚ and ℚ(i) through the
    /// numerators.
    Modular,
    /// Subresultant remainder sequences \[BT71, Bro78, Coh93\], in the first
    /// variable over the gcds of the coefficients in the others: every
    /// domain. The reference.
    Subresultant,
    /// Zippel's sparse interpolation \[Zip79, Zip90\] in de Kleine, Monagan
    /// and Wittkopf's form, LINZIP \[DKMW05\]: the form of the gcd from a
    /// first image, the others by solving for its coefficients. For sparse
    /// polynomials over prime fields of at least `2^31` elements, and over
    /// ℤ, from the form of a first prime's image, and ℚ through the
    /// numerators. It can give up.
    Zippel,
    /// EZ-GCD \[MY73, Wang80, GCL92 §7.6\]: images at integer points of
    /// the variables other than the main one, their univariate gcd, and
    /// multivariate Hensel lifting back to the gcd (CAS-526). Over ℤ, sparse
    /// in several variables, Brown's algorithm in one, and ℚ through the
    /// numerators. It can give up.
    EzGcd,
}

mod sealed {
    pub trait Sealed {}
}

/// The coefficient domains whose polynomials have gcds, resultants and
/// discriminants here: ℤ, ℚ, ℤ\[i\], ℚ(i), and the prime fields of word and
/// of any size. Each chooses its algorithms. The trait is sealed; its
/// methods are the hooks behind [`DenseRing::gcd`], [`SparseRing::gcd`] and
/// their kin, which are the interface.
pub trait GcdDomain: Coefficients + EuclideanDomain + Clone + sealed::Sealed {
    /// The gcd of non-zero `a` and `b` with its cofactors, not normalised,
    /// by the domain's own choice of algorithm.
    #[doc(hidden)]
    fn dense_gcd(
        ring: &DenseRing<Self>,
        a: &Dense<Self::Elem>,
        b: &Dense<Self::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Dense<Self::Elem>>, Cancelled>;

    /// As [`GcdDomain::dense_gcd`], by `method`: `None` when it does not
    /// apply over this domain, or gives up.
    #[doc(hidden)]
    fn dense_gcd_by(
        ring: &DenseRing<Self>,
        a: &Dense<Self::Elem>,
        b: &Dense<Self::Elem>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Dense<Self::Elem>>>, Cancelled>;

    /// As [`GcdDomain::dense_gcd`], for sparse polynomials.
    #[doc(hidden)]
    fn sparse_gcd(
        ring: &SparseRing<Self>,
        a: &Sparse<Self::Elem>,
        b: &Sparse<Self::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Sparse<Self::Elem>>, Cancelled>;

    /// As [`GcdDomain::dense_gcd_by`], for sparse polynomials.
    #[doc(hidden)]
    fn sparse_gcd_by(
        ring: &SparseRing<Self>,
        a: &Sparse<Self::Elem>,
        b: &Sparse<Self::Elem>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Sparse<Self::Elem>>>, Cancelled>;

    /// The resultant of non-zero `a` and `b`.
    #[doc(hidden)]
    fn dense_resultant(
        ring: &DenseRing<Self>,
        a: &Dense<Self::Elem>,
        b: &Dense<Self::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Self::Elem, Cancelled> {
        subres::resultant(ring, a, b, cancel)
    }

    /// The resultant of non-zero `a` and `b` in `x_var`.
    #[doc(hidden)]
    fn sparse_resultant(
        ring: &SparseRing<Self>,
        a: &Sparse<Self::Elem>,
        b: &Sparse<Self::Elem>,
        var: usize,
        cancel: &dyn Cancel,
    ) -> Result<Sparse<Self::Elem>, Cancelled> {
        subres::sparse_resultant(ring, a, b, var, cancel)
    }
}

/// The inverse of a unit of `ring`.
fn unit_inverse<R: EuclideanDomain>(ring: &R, unit: &R::Elem) -> R::Elem {
    ring.div_rem(&ring.one(), unit).0
}

/// `p` times the unit that puts its leading coefficient in normal form.
pub(crate) fn normalized_dense<R: Coefficients + EuclideanDomain>(
    ring: &DenseRing<R>,
    p: &Dense<R::Elem>,
) -> Dense<R::Elem> {
    match p.lead() {
        Some(lead) => {
            let unit = ring.base().normal_form(lead).1;
            if ring.base().is_one(&unit) {
                p.clone()
            } else {
                ring.scale(p, &unit)
            }
        }
        None => p.clone(),
    }
}

/// `p` times the unit that puts its leading coefficient, in lexicographic
/// order, in normal form.
pub(crate) fn normalized_sparse<R: EuclideanDomain>(
    ring: &SparseRing<R>,
    p: &Sparse<R::Elem>,
) -> Sparse<R::Elem> {
    match p.coeffs().first() {
        Some(lead) => {
            let unit = ring.base().normal_form(lead).1;
            if ring.base().is_one(&unit) {
                p.clone()
            } else {
                ring.scale(p, &unit)
            }
        }
        None => p.clone(),
    }
}

/// The gcd scaled into normal form, the cofactors by the inverse unit.
fn normal_cofactors_dense<R: Coefficients + EuclideanDomain>(
    ring: &DenseRing<R>,
    (g, a, b): Cofactors<Dense<R::Elem>>,
) -> Cofactors<Dense<R::Elem>> {
    let Some(lead) = g.lead() else {
        return (g, a, b);
    };
    let base = ring.base();
    let unit = base.normal_form(lead).1;
    if base.is_one(&unit) {
        return (g, a, b);
    }
    let inverse = unit_inverse(base, &unit);
    (
        ring.scale(&g, &unit),
        ring.scale(&a, &inverse),
        ring.scale(&b, &inverse),
    )
}

/// As [`normal_cofactors_dense`], for sparse polynomials.
fn normal_cofactors_sparse<R: EuclideanDomain>(
    ring: &SparseRing<R>,
    (g, a, b): Cofactors<Sparse<R::Elem>>,
) -> Cofactors<Sparse<R::Elem>> {
    let Some(lead) = g.coeffs().first() else {
        return (g, a, b);
    };
    let base = ring.base();
    let unit = base.normal_form(lead).1;
    if base.is_one(&unit) {
        return (g, a, b);
    }
    let inverse = unit_inverse(base, &unit);
    (
        ring.scale(&g, &unit),
        ring.scale(&a, &inverse),
        ring.scale(&b, &inverse),
    )
}

/// Whether `(n (n - 1) / 2)` is odd: the sign of a discriminant.
fn discriminant_negates(n: u64) -> bool {
    (n % 4) >= 2
}

impl<R: GcdDomain> DenseRing<R> {
    /// The greatest common divisor, normalised: its leading coefficient in
    /// the normal form of `R`, so monic over a field and positive over ℤ,
    /// and over ℤ and ℤ\[i\] the gcd of the contents a factor of it.
    /// `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, a: &Dense<R::Elem>, b: &Dense<R::Elem>) -> Dense<R::Elem> {
        uncancelled(self.gcd_with(a, b, &Never))
    }

    /// As [`DenseRing::gcd`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_with(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Dense<R::Elem>, Cancelled> {
        Ok(self.gcd_cofactors_with(a, b, cancel)?.0)
    }

    /// `(g, a / g, b / g)` for `g` the gcd, normalised as by
    /// [`DenseRing::gcd`]; `(0, 0, 0)` for two zeros.
    #[must_use]
    pub fn gcd_cofactors(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
    ) -> Cofactors<Dense<R::Elem>> {
        uncancelled(self.gcd_cofactors_with(a, b, &Never))
    }

    /// As [`DenseRing::gcd_cofactors`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_cofactors_with(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Dense<R::Elem>>, Cancelled> {
        let parts = match self.with_zero(a, b) {
            Some(parts) => parts,
            None => R::dense_gcd(self, a, b, cancel)?,
        };
        Ok(normal_cofactors_dense(self, parts))
    }

    /// The gcd and its cofactors by `method`, normalised as by
    /// [`DenseRing::gcd`], polling `cancel`. `None` when the method does
    /// not apply over `R`, or gives up, as the heuristic can.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_by(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Dense<R::Elem>>>, Cancelled> {
        let parts = match self.with_zero(a, b) {
            Some(parts) => Some(parts),
            None => R::dense_gcd_by(self, a, b, method, cancel)?,
        };
        Ok(parts.map(|parts| normal_cofactors_dense(self, parts)))
    }

    /// The gcd and cofactors when `a` or `b` is zero.
    fn with_zero(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
    ) -> Option<Cofactors<Dense<R::Elem>>> {
        match (a.is_zero(), b.is_zero()) {
            (true, true) => Some((self.zero(), self.zero(), self.zero())),
            (true, false) => Some((b.clone(), self.zero(), self.one())),
            (false, true) => Some((a.clone(), self.one(), self.zero())),
            (false, false) => None,
        }
    }

    /// The least common multiple `a (b / gcd(a, b))`, normalised as the
    /// gcd; 0 when either is 0.
    #[must_use]
    pub fn lcm(&self, a: &Dense<R::Elem>, b: &Dense<R::Elem>) -> Dense<R::Elem> {
        let (_, _, b_cofactor) = self.gcd_cofactors(a, b);
        normalized_dense(self, &self.mul(a, &b_cofactor))
    }

    /// The resultant `lc(a)^(deg b) ∏ b(α)` over the roots `α` of `a`, the
    /// determinant of the Sylvester matrix: 0 exactly when `a` and `b` have
    /// a common factor, `c^(deg b)` for a constant `a = c`, and 0 when
    /// either is 0.
    #[must_use]
    pub fn resultant(&self, a: &Dense<R::Elem>, b: &Dense<R::Elem>) -> R::Elem {
        uncancelled(self.resultant_with(a, b, &Never))
    }

    /// As [`DenseRing::resultant`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn resultant_with(
        &self,
        a: &Dense<R::Elem>,
        b: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<R::Elem, Cancelled> {
        if a.is_zero() || b.is_zero() {
            return Ok(self.base().zero());
        }
        R::dense_resultant(self, a, b, cancel)
    }

    /// The discriminant `(-1)^(n (n - 1) / 2) res(p, p') / lc(p)` of a
    /// polynomial of degree `n ≥ 1`, `lc(p)^(2n - 2) ∏ (αᵢ - αⱼ)²` over the
    /// pairs of its roots: 0 exactly when it has a repeated root, 1 in
    /// degree 1. `None` for a constant or zero.
    #[must_use]
    pub fn discriminant(&self, p: &Dense<R::Elem>) -> Option<R::Elem> {
        uncancelled(self.discriminant_with(p, &Never))
    }

    /// As [`DenseRing::discriminant`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn discriminant_with(
        &self,
        p: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Option<R::Elem>, Cancelled> {
        let (Some(n), Some(lead)) = (p.degree(), p.lead()) else {
            return Ok(None);
        };
        if n == 0 {
            return Ok(None);
        }
        let base = self.base();
        let res = self.resultant_with(p, &self.derivative(p), cancel)?;
        let (d, r) = base.div_rem(&res, lead);
        debug_assert!(base.is_zero(&r), "the leading coefficient divides");
        Ok(Some(if discriminant_negates(n as u64) {
            base.neg(&d)
        } else {
            d
        }))
    }
}

impl<R: GcdDomain> SparseRing<R> {
    /// The greatest common divisor, normalised: its leading coefficient in
    /// lexicographic order, the first variable most significant, in the
    /// normal form of `R`, so monic over a field and positive over ℤ, and
    /// over ℤ and ℤ\[i\] the gcd of the contents a factor of it.
    /// `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> Sparse<R::Elem> {
        uncancelled(self.gcd_with(a, b, &Never))
    }

    /// As [`SparseRing::gcd`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Sparse<R::Elem>, Cancelled> {
        Ok(self.gcd_cofactors_with(a, b, cancel)?.0)
    }

    /// `(g, a / g, b / g)` for `g` the gcd, normalised as by
    /// [`SparseRing::gcd`]; `(0, 0, 0)` for two zeros.
    #[must_use]
    pub fn gcd_cofactors(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
    ) -> Cofactors<Sparse<R::Elem>> {
        uncancelled(self.gcd_cofactors_with(a, b, &Never))
    }

    /// As [`SparseRing::gcd_cofactors`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_cofactors_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Sparse<R::Elem>>, Cancelled> {
        let parts = match self.with_zero(a, b) {
            Some(parts) => parts,
            None => R::sparse_gcd(self, a, b, cancel)?,
        };
        Ok(normal_cofactors_sparse(self, parts))
    }

    /// The gcd and its cofactors by `method`, normalised as by
    /// [`SparseRing::gcd`], polling `cancel`. `None` when the method does
    /// not apply over `R`, or gives up, as the heuristic can.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn gcd_by(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Sparse<R::Elem>>>, Cancelled> {
        let parts = match self.with_zero(a, b) {
            Some(parts) => Some(parts),
            None => R::sparse_gcd_by(self, a, b, method, cancel)?,
        };
        Ok(parts.map(|parts| normal_cofactors_sparse(self, parts)))
    }

    /// The gcd and cofactors when `a` or `b` is zero.
    fn with_zero(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
    ) -> Option<Cofactors<Sparse<R::Elem>>> {
        match (a.is_zero(), b.is_zero()) {
            (true, true) => Some((self.zero(), self.zero(), self.zero())),
            (true, false) => Some((b.clone(), self.zero(), self.one())),
            (false, true) => Some((a.clone(), self.one(), self.zero())),
            (false, false) => None,
        }
    }

    /// The least common multiple `a (b / gcd(a, b))`, normalised as the
    /// gcd; 0 when either is 0.
    #[must_use]
    pub fn lcm(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> Sparse<R::Elem> {
        let (_, _, b_cofactor) = self.gcd_cofactors(a, b);
        normalized_sparse(self, &self.mul(a, &b_cofactor))
    }

    /// The content: the gcd of the coefficients, in the normal form of
    /// `R` (non-negative over ℤ, 1 over a field); 0 for zero.
    #[must_use]
    pub fn content(&self, p: &Sparse<R::Elem>) -> R::Elem {
        let base = self.base();
        let mut g = base.zero();
        for c in p.coeffs() {
            g = base.gcd(&g, c);
            if base.is_one(&g) {
                break;
            }
        }
        g
    }

    /// The primitive part: `p` over its content, times the unit that puts
    /// its leading coefficient in normal form, so that `p` is a unit times
    /// `content(p) · primitive_part(p)`. Zero for zero.
    #[must_use]
    pub fn primitive_part(&self, p: &Sparse<R::Elem>) -> Sparse<R::Elem> {
        let c = self.content(p);
        if p.is_zero() || self.base().is_one(&c) {
            return normalized_sparse(self, p);
        }
        let base = self.base();
        normalized_sparse(self, &self.map_into(p, self, |x| base.div_rem(x, &c).0))
    }

    /// The content in `x_var`: the gcd of the coefficients of `p` as a
    /// polynomial in `x_var`, a polynomial without it, normalised as the
    /// gcd; 0 for zero.
    #[must_use]
    pub fn content_in(&self, p: &Sparse<R::Elem>, var: usize) -> Sparse<R::Elem> {
        let mut g = self.zero();
        for (_, c) in self.to_recursive(p, var) {
            g = self.gcd(&g, &c);
            if self.is_one(&g) {
                break;
            }
        }
        g
    }

    /// The primitive part in `x_var`: `p` over its content in `x_var`,
    /// normalised as the gcd. Zero for zero.
    #[must_use]
    pub fn primitive_part_in(&self, p: &Sparse<R::Elem>, var: usize) -> Sparse<R::Elem> {
        if p.is_zero() {
            return self.zero();
        }
        let c = self.content_in(p, var);
        match self.div_exact(p, &c) {
            Some(q) => normalized_sparse(self, &q),
            None => unreachable!("the content divides"),
        }
    }

    /// The resultant in `x_var`, a polynomial in the other variables: the
    /// determinant of the Sylvester matrix of `a` and `b` as polynomials in
    /// `x_var`, 0 exactly when they have a common factor that has `x_var`.
    /// 0 when either is 0, `c^(deg b)` for an `a = c` without `x_var`.
    ///
    /// # Panics
    ///
    /// If `var` is not a variable of the ring.
    #[must_use]
    pub fn resultant(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        var: usize,
    ) -> Sparse<R::Elem> {
        uncancelled(self.resultant_with(a, b, var, &Never))
    }

    /// As [`SparseRing::resultant`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `var` is not a variable of the ring.
    pub fn resultant_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        var: usize,
        cancel: &dyn Cancel,
    ) -> Result<Sparse<R::Elem>, Cancelled> {
        assert!(var < self.nvars(), "variable {var} of {}", self.nvars());
        if a.is_zero() || b.is_zero() {
            return Ok(self.zero());
        }
        R::sparse_resultant(self, a, b, var, cancel)
    }

    /// The discriminant in `x_var`, `(-1)^(n (n - 1) / 2) res(p, p') /
    /// lc(p)` for `p` of degree `n ≥ 1` in `x_var`, the derivative and the
    /// leading coefficient in `x_var`: a polynomial in the other variables,
    /// `b² - 4ac` for `a x² + b x + c`. `None` when `p` is zero or has no
    /// `x_var`.
    ///
    /// # Panics
    ///
    /// If `var` is not a variable of the ring.
    #[must_use]
    pub fn discriminant(&self, p: &Sparse<R::Elem>, var: usize) -> Option<Sparse<R::Elem>> {
        uncancelled(self.discriminant_with(p, var, &Never))
    }

    /// As [`SparseRing::discriminant`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `var` is not a variable of the ring.
    pub fn discriminant_with(
        &self,
        p: &Sparse<R::Elem>,
        var: usize,
        cancel: &dyn Cancel,
    ) -> Result<Option<Sparse<R::Elem>>, Cancelled> {
        assert!(var < self.nvars(), "variable {var} of {}", self.nvars());
        let n = match self.degree(p, var) {
            None | Some(0) => return Ok(None),
            Some(n) => n,
        };
        let res = self.resultant_with(p, &self.derivative(p, var), var, cancel)?;
        let d = self
            .div_exact_with(&res, &self.lead_coeff(p, var), cancel)?
            .expect("the leading coefficient divides");
        Ok(Some(if discriminant_negates(n) {
            self.neg(&d)
        } else {
            d
        }))
    }

    /// The extended gcd in `x_var`: `(g, u, v, c)` with `u a + v b = c g`,
    /// `g` the gcd of [`SparseRing::gcd`] and `c` a non-zero polynomial
    /// without `x_var`, giac's `Tegcdpsr` and `Tegcdtoabcuv` (poly.h): the
    /// extended Euclidean algorithm on pseudo-remainders
    /// ([`SparseRing::pseudo_div_rem`]), each remainder and its two
    /// cofactors divided by the gcd of their contents in `x_var`. `(0, 0,
    /// 0, 1)` for two zeros.
    ///
    /// # Panics
    ///
    /// Never: the gcd divides the last remainder.
    #[must_use]
    pub fn extended_gcd(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        var: usize,
    ) -> Bezout<Sparse<R::Elem>> {
        if a.is_zero() && b.is_zero() {
            return (self.zero(), self.zero(), self.zero(), self.one());
        }
        let (mut r0, mut s0, mut t0) = (a.clone(), self.one(), self.zero());
        let (mut r1, mut s1, mut t1) = (b.clone(), self.zero(), self.one());
        while !r1.is_zero() {
            // `l^δ r0 = q r1 + r`, so `r = (l^δ s0 − q s1) a + (l^δ t0 − q t1) b`.
            let delta = pseudo_exponent(self.degree(&r0, var), self.degree(&r1, var));
            let f = self.pow(&self.lead_coeff(&r1, var), delta);
            let (q, r) = self.pseudo_div_rem(&r0, &r1, var);
            let s = self.sub(&self.mul(&f, &s0), &self.mul(&q, &s1));
            let t = self.sub(&self.mul(&f, &t0), &self.mul(&q, &t1));
            let [r, s, t] = self.without_common_content([r, s, t], var);
            r0 = std::mem::replace(&mut r1, r);
            s0 = std::mem::replace(&mut s1, s);
            t0 = std::mem::replace(&mut t1, t);
        }
        // `s0 a + t0 b = r0 = h g/k` for `g = k g₁`, `g₁` primitive in
        // `x_var` and `h`, `k` without it.
        let g = self.gcd(a, b);
        let g1 = self.primitive_part_in(&g, var);
        let h = self
            .div_exact(&r0, &g1)
            .expect("the gcd divides the last remainder");
        let k = self.div_exact(&g, &g1).expect("a primitive part divides");
        let [u, v, c] = self.without_common_content([self.mul(&s0, &k), self.mul(&t0, &k), h], var);
        (g, u, v, c)
    }

    /// `(u, v, d)` with `u a + v b = d c`, `deg u < deg b` in `x_var` when
    /// `b` has `x_var`, and `d` a non-zero polynomial without `x_var`:
    /// giac's `Tabcuv` (poly.h), from [`SparseRing::extended_gcd`]. `None`
    /// when the gcd of `a` and `b` does not divide `c` as polynomials in
    /// `x_var`.
    #[must_use]
    pub fn abcuv(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        c: &Sparse<R::Elem>,
        var: usize,
    ) -> Option<Abcuv<Sparse<R::Elem>>> {
        if c.is_zero() {
            return Some((self.zero(), self.zero(), self.one()));
        }
        let (g, u0, v0, c0) = self.extended_gcd(a, b, var);
        if g.is_zero() {
            return None;
        }
        // `l^δ c = w g`, so `u0 w a + v0 w b = c0 l^δ c`.
        let delta = pseudo_exponent(self.degree(c, var), self.degree(&g, var));
        let (w, rem) = self.pseudo_div_rem(c, &g, var);
        if !rem.is_zero() {
            return None;
        }
        let mut d = self.mul(&c0, &self.pow(&self.lead_coeff(&g, var), delta));
        let (mut u, mut v) = (self.mul(&u0, &w), self.mul(&v0, &w));
        if self.degree(b, var).is_some_and(|db| db > 0) {
            // `l^δ u = q b + u′`: `u′ a + (l^δ v + q a) b = l^δ d c`.
            let delta = pseudo_exponent(self.degree(&u, var), self.degree(b, var));
            let f = self.pow(&self.lead_coeff(b, var), delta);
            let (q, r) = self.pseudo_div_rem(&u, b, var);
            u = r;
            v = self.add(&self.mul(&f, &v), &self.mul(&q, a));
            d = self.mul(&f, &d);
        }
        let [u, v, d] = self.without_common_content([u, v, d], var);
        Some((u, v, d))
    }

    /// The partial fractions of `num / den` in `x_var` over a square-free
    /// decomposition `den = κ f₁^m₁ ⋯ fₖ^mₖ` of the denominator, `κ`
    /// without `x_var` and the `fᵢ` pairwise coprime, giac's `Tpartfrac`
    /// (poly.h): `num / den = q / d + Σ nᵢ / (dᵢ fᵢ^mᵢ)` with
    /// `deg nᵢ < deg fᵢ^mᵢ` in `x_var` and `d`, `dᵢ` without it. The
    /// quotient comes from [`SparseRing::pseudo_div_rem`], each piece from
    /// [`SparseRing::abcuv`] against the product of the factors after it.
    /// `None` when the factors do not divide `den` with a cofactor `κ`
    /// without `x_var`, or two of them are not coprime.
    ///
    /// # Panics
    ///
    /// If `den` is zero.
    #[must_use]
    pub fn partial_fractions(
        &self,
        num: &Sparse<R::Elem>,
        den: &Sparse<R::Elem>,
        factors: &[(Sparse<R::Elem>, u64)],
        var: usize,
    ) -> Option<PartialFractions<R::Elem>> {
        let powers: Vec<Sparse<R::Elem>> = factors.iter().map(|(f, m)| self.pow(f, *m)).collect();
        let product = powers.iter().fold(self.one(), |p, f| self.mul(&p, f));
        let kappa = self.div_exact(den, &product)?;
        if self.degree(&kappa, var).is_some_and(|d| d > 0) {
            return None;
        }
        let delta = pseudo_exponent(self.degree(num, var), self.degree(den, var));
        let quotient_den = self.pow(&self.lead_coeff(den, var), delta);
        let (quotient, mut rest) = self.pseudo_div_rem(num, den, var);
        // `rest / (c · Πⱼ≥ᵢ fⱼ^mⱼ)` remains to split.
        let mut c = self.mul(&quotient_den, &kappa);
        let mut after = product;
        let mut pieces = Vec::with_capacity(factors.len());
        for (k, p) in powers.iter().enumerate() {
            if k + 1 == powers.len() {
                pieces.push(self.reduced_piece(rest, c, var));
                break;
            }
            after = self.div_exact(&after, p).expect("a factor of the product");
            if self
                .degree(&self.gcd(&after, p), var)
                .is_some_and(|d| d > 0)
            {
                return None;
            }
            let (u, v, d) = self.abcuv(&after, p, &rest, var)?;
            // `d rest = u Q + v P`: `rest / (c P Q) = u / (c d P) + v / (c d Q)`.
            c = self.mul(&c, &d);
            pieces.push(self.reduced_piece(u, c.clone(), var));
            rest = v;
        }
        Some(PartialFractions {
            quotient,
            quotient_den,
            pieces,
        })
    }

    /// `n / d` for `d` without `x_var`, both divided by their common
    /// content.
    fn reduced_piece(
        &self,
        n: Sparse<R::Elem>,
        d: Sparse<R::Elem>,
        var: usize,
    ) -> (Sparse<R::Elem>, Sparse<R::Elem>) {
        let [n, d] = self.without_common_content([n, d], var);
        (n, d)
    }

    /// The polynomials `ps` divided by the gcd of their contents in
    /// `x_var`; unchanged when all are zero.
    fn without_common_content<const N: usize>(
        &self,
        ps: [Sparse<R::Elem>; N],
        var: usize,
    ) -> [Sparse<R::Elem>; N] {
        let mut g = self.zero();
        for p in &ps {
            g = self.gcd(&g, &self.content_in(p, var));
            if self.is_one(&g) {
                return ps;
            }
        }
        if g.is_zero() {
            return ps;
        }
        ps.map(|p| self.div_exact(&p, &g).expect("the content divides"))
    }
}

/// `max(m − n + 1, 0)` for the degrees `m` and `n` of a pseudo-division, 0
/// for a zero dividend.
fn pseudo_exponent(m: Option<u64>, n: Option<u64>) -> u64 {
    match (m, n) {
        (Some(m), Some(n)) if m >= n => m - n + 1,
        _ => 0,
    }
}

/// The partial fractions of [`SparseRing::partial_fractions`]: `num / den =
/// quotient / quotient_den + Σ n / d (f^m)` over its `pieces` `(n, d)`, one a
/// factor `f^m` in their order.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct PartialFractions<E> {
    /// The polynomial part's numerator.
    pub quotient: Sparse<E>,
    /// Its denominator, without the variable.
    pub quotient_den: Sparse<E>,
    /// `(nᵢ, dᵢ)`: the numerator and the denominator without the variable
    /// of each factor's piece.
    pub pieces: Vec<(Sparse<E>, Sparse<E>)>,
}

impl<N: Numerators> SharedRing<N> {
    /// The resultant over the field of fractions, from that of the
    /// numerators: `res(n / d, n' / d') = res(n, n') / (d^(deg n') d'^(deg
    /// n))`. 0 when either is 0.
    #[must_use]
    pub fn resultant(
        &self,
        a: &Shared<N::Elem>,
        b: &Shared<N::Elem>,
    ) -> <N::Fractions as Ring>::Elem {
        let (Some(m), Some(n)) = (a.degree(), b.degree()) else {
            return N::fraction(&self.numerators().base().zero(), &Natural::one());
        };
        let res = self.numerators().resultant(a.numerator(), b.numerator());
        let den = &a.denominator().pow(n as u64) * &b.denominator().pow(m as u64);
        N::fraction(&res, &den)
    }

    /// The discriminant, from that of the numerator: `disc(n / d) = disc(n)
    /// / d^(2 deg n - 2)`. `None` for a constant or zero.
    #[must_use]
    pub fn discriminant(&self, p: &Shared<N::Elem>) -> Option<<N::Fractions as Ring>::Elem> {
        let d = self.numerators().discriminant(p.numerator())?;
        let n = p.degree()? as u64;
        Some(N::fraction(&d, &p.denominator().pow(2 * n - 2)))
    }
}

/// The gcd over a field by Euclid's algorithm, monic, with its cofactors.
fn euclid<F: Field + Coefficients>(
    ring: &DenseRing<F>,
    a: &Dense<F::Elem>,
    b: &Dense<F::Elem>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Dense<F::Elem>>, Cancelled> {
    let (mut x, mut y) = (a.clone(), b.clone());
    while !y.is_zero() {
        let r = ring.div_rem_with(&x, &y, cancel)?.1;
        (x, y) = (y, r);
    }
    let g = normalized_dense(ring, &x);
    let a_cofactor = ring.div_rem_with(a, &g, cancel)?.0;
    let b_cofactor = ring.div_rem_with(b, &g, cancel)?.0;
    Ok((g, a_cofactor, b_cofactor))
}

/// `(g, a / g, b / g)` for a gcd `g` of `a` and `b`.
fn cofactors_of<R: GcdDomain>(
    ring: &SparseRing<R>,
    g: Sparse<R::Elem>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Sparse<R::Elem>>, Cancelled> {
    let quotient = |p: &Sparse<R::Elem>| match ring.div_exact_with(p, &g, cancel)? {
        Some(q) => Ok(q),
        None => unreachable!("a gcd divides"),
    };
    let (a_cofactor, b_cofactor) = (quotient(a)?, quotient(b)?);
    Ok((g, a_cofactor, b_cofactor))
}

/// Whether `a` and `b` have two variables or more between them.
fn several_variables<E>(a: &Sparse<E>, b: &Sparse<E>) -> bool {
    let (da, db) = (a.degrees(), b.degrees());
    da.iter()
        .zip(&db)
        .filter(|(x, y)| **x > 0 || **y > 0)
        .count()
        > 1
}

/// Points of Zippel's algorithm are drawn from fields of at least this many
/// elements: fewer make its random points collide.
const ZIPPEL_FIELD: u64 = 1 << 31;

/// The gcd over a prime field by Zippel's algorithm, with its cofactors,
/// which test it: `None` when the field is small, the algorithm gives up, or
/// its gcd does not divide.
fn zippel_cofactors<F: GcdDomain + Field>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Sparse<F::Elem>>>, Cancelled> {
    if u64::try_from(&ring.base().characteristic()).is_ok_and(|p| p < ZIPPEL_FIELD) {
        return Ok(None);
    }
    let Some(g) = zippel::field_gcd(ring, a, b, cancel)? else {
        return Ok(None);
    };
    let Some(a_cofactor) = ring.div_exact_with(a, &g, cancel)? else {
        return Ok(None);
    };
    let Some(b_cofactor) = ring.div_exact_with(b, &g, cancel)? else {
        return Ok(None);
    };
    Ok(Some((g, a_cofactor, b_cofactor)))
}

macro_rules! sealed {
    ($($t:ty),*) => {$(
        impl sealed::Sealed for $t {}
    )*};
}

sealed!(
    Integers,
    Rationals,
    GaussianIntegers,
    GaussianRationals,
    PrimeField<ModWord>,
    PrimeField<ModRing>
);

/// ℤ: the heuristic, which gives up on large problems, then Brown's
/// algorithm, or in several variables Zippel's first; resultants modulo
/// primes in one variable.
impl GcdDomain for Integers {
    fn dense_gcd(
        ring: &DenseRing<Integers>,
        a: &Dense<Integer>,
        b: &Dense<Integer>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Dense<Integer>>, Cancelled> {
        match heuristic::dense(ring, a, b, cancel)? {
            Some(parts) => Ok(parts),
            None => modular::dense(ring, a, b, cancel),
        }
    }

    fn dense_gcd_by(
        ring: &DenseRing<Integers>,
        a: &Dense<Integer>,
        b: &Dense<Integer>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Dense<Integer>>>, Cancelled> {
        match method {
            GcdMethod::Heuristic => heuristic::dense(ring, a, b, cancel),
            GcdMethod::Modular | GcdMethod::EzGcd => Ok(Some(modular::dense(ring, a, b, cancel)?)),
            GcdMethod::Subresultant => Ok(Some(subres::dense_gcd(ring, a, b, cancel)?)),
            GcdMethod::Zippel => Ok(None),
        }
    }

    fn sparse_gcd(
        ring: &SparseRing<Integers>,
        a: &Sparse<Integer>,
        b: &Sparse<Integer>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Sparse<Integer>>, Cancelled> {
        if let Some(parts) = heuristic::sparse(ring, a, b, cancel)? {
            return Ok(parts);
        }
        if several_variables(a, b)
            && let Some(parts) = zippel::integers(ring, a, b, cancel)?
        {
            return Ok(parts);
        }
        modular::sparse(ring, a, b, cancel)
    }

    fn sparse_gcd_by(
        ring: &SparseRing<Integers>,
        a: &Sparse<Integer>,
        b: &Sparse<Integer>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Sparse<Integer>>>, Cancelled> {
        match method {
            GcdMethod::Heuristic => heuristic::sparse(ring, a, b, cancel),
            GcdMethod::EzGcd if several_variables(a, b) => ezgcd::integers(ring, a, b, cancel),
            // EZ-GCD in one variable: Brown's algorithm.
            GcdMethod::Modular | GcdMethod::EzGcd => Ok(Some(modular::sparse(ring, a, b, cancel)?)),
            GcdMethod::Subresultant => Ok(Some(subres::sparse_gcd(ring, a, b, cancel)?)),
            GcdMethod::Zippel => zippel::integers(ring, a, b, cancel),
        }
    }

    fn dense_resultant(
        _: &DenseRing<Integers>,
        a: &Dense<Integer>,
        b: &Dense<Integer>,
        cancel: &dyn Cancel,
    ) -> Result<Integer, Cancelled> {
        resultant::integers(a, b, cancel)
    }

    fn sparse_resultant(
        ring: &SparseRing<Integers>,
        a: &Sparse<Integer>,
        b: &Sparse<Integer>,
        var: usize,
        cancel: &dyn Cancel,
    ) -> Result<Sparse<Integer>, Cancelled> {
        // In one variable, the modular resultant.
        if let (Some(da), Some(db)) = (ring.to_dense(a, var), ring.to_dense(b, var)) {
            return Ok(ring.constant(resultant::integers(&da, &db, cancel)?));
        }
        subres::sparse_resultant(ring, a, b, var, cancel)
    }
}

/// ℤ\[i\]: the heuristic, then subresultants.
impl GcdDomain for GaussianIntegers {
    fn dense_gcd(
        ring: &DenseRing<GaussianIntegers>,
        a: &Dense<Complex>,
        b: &Dense<Complex>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Dense<Complex>>, Cancelled> {
        match heuristic::dense(ring, a, b, cancel)? {
            Some(parts) => Ok(parts),
            None => subres::dense_gcd(ring, a, b, cancel),
        }
    }

    fn dense_gcd_by(
        ring: &DenseRing<GaussianIntegers>,
        a: &Dense<Complex>,
        b: &Dense<Complex>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Dense<Complex>>>, Cancelled> {
        match method {
            GcdMethod::Heuristic => heuristic::dense(ring, a, b, cancel),
            GcdMethod::Modular | GcdMethod::Zippel | GcdMethod::EzGcd => Ok(None),
            GcdMethod::Subresultant => Ok(Some(subres::dense_gcd(ring, a, b, cancel)?)),
        }
    }

    fn sparse_gcd(
        ring: &SparseRing<GaussianIntegers>,
        a: &Sparse<Complex>,
        b: &Sparse<Complex>,
        cancel: &dyn Cancel,
    ) -> Result<Cofactors<Sparse<Complex>>, Cancelled> {
        // Brown's algorithm modulo Gaussian primes where the heuristic gives
        // up (CAS-494): subresultants took minutes in two variables of
        // degree 150.
        match heuristic::sparse(ring, a, b, cancel)? {
            Some(parts) => Ok(parts),
            None => gaussian::sparse(ring, a, b, cancel),
        }
    }

    fn sparse_gcd_by(
        ring: &SparseRing<GaussianIntegers>,
        a: &Sparse<Complex>,
        b: &Sparse<Complex>,
        method: GcdMethod,
        cancel: &dyn Cancel,
    ) -> Result<Option<Cofactors<Sparse<Complex>>>, Cancelled> {
        match method {
            GcdMethod::Heuristic => heuristic::sparse(ring, a, b, cancel),
            GcdMethod::Modular => Ok(Some(gaussian::sparse(ring, a, b, cancel)?)),
            GcdMethod::Zippel | GcdMethod::EzGcd => Ok(None),
            GcdMethod::Subresultant => Ok(Some(subres::sparse_gcd(ring, a, b, cancel)?)),
        }
    }
}

/// Coefficients over one natural denominator: `(n, d)` with `p = n / d`.
pub(crate) fn split<N: Numerators>(coeffs: &[Fraction<N>]) -> (Vec<N::Elem>, Natural) {
    let parts: Vec<(N::Elem, Natural)> = coeffs.iter().map(N::split).collect();
    let den = parts.iter().fold(Natural::one(), |l, (_, d)| l.lcm(d));
    let nums = parts
        .into_iter()
        .map(|(c, d)| N::mul_natural(&c, &den.div_exact(&d)))
        .collect();
    (nums, den)
}

/// `p = n / d` over ℚ or ℚ(i), `n` over ℤ or ℤ\[i\].
pub(crate) fn split_sparse<N: Numerators>(
    ring: &SparseRing<N::Fractions>,
    numerators: &SparseRing<N>,
    p: &Sparse<Fraction<N>>,
) -> (Sparse<N::Elem>, Natural)
where
    N::Fractions: Clone,
{
    let den = p
        .coeffs()
        .iter()
        .fold(Natural::one(), |l, c| l.lcm(&N::split(c).1));
    let num = ring.map_into(p, numerators, |c| {
        let (n, d) = N::split(c);
        N::mul_natural(&n, &den.div_exact(&d))
    });
    (num, den)
}

/// The numerators of `a` and `b` over ℤ or ℤ\[i\], and their denominators.
struct Split<P> {
    a: P,
    b: P,
    da: Natural,
    db: Natural,
}

fn split_dense<N: Numerators>(
    numerators: &DenseRing<N>,
    a: &Dense<Fraction<N>>,
    b: &Dense<Fraction<N>>,
) -> Split<Dense<N::Elem>> {
    let ((na, da), (nb, db)) = (split::<N>(a.coeffs()), split::<N>(b.coeffs()));
    Split {
        a: numerators.poly(na),
        b: numerators.poly(nb),
        da,
        db,
    }
}

/// Over the fractions: `a = n / d` and `n = g q` give the gcd `g` and the
/// cofactor `q / d`.
fn lift_dense<N: Numerators>(
    ring: &DenseRing<N::Fractions>,
    split: &Split<Dense<N::Elem>>,
    (g, qa, qb): Cofactors<Dense<N::Elem>>,
) -> Cofactors<Dense<Fraction<N>>>
where
    N::Fractions: Coefficients,
{
    let lift = |p: &Dense<N::Elem>, den: &Natural| {
        ring.poly(p.coeffs().iter().map(|c| N::fraction(c, den)).collect())
    };
    (
        lift(&g, &Natural::one()),
        lift(&qa, &split.da),
        lift(&qb, &split.db),
    )
}

fn split_sparse_pair<N: Numerators>(
    ring: &SparseRing<N::Fractions>,
    numerators: &SparseRing<N>,
    a: &Sparse<Fraction<N>>,
    b: &Sparse<Fraction<N>>,
) -> Split<Sparse<N::Elem>>
where
    N::Fractions: Clone,
{
    let ((na, da), (nb, db)) = (
        split_sparse(ring, numerators, a),
        split_sparse(ring, numerators, b),
    );
    Split {
        a: na,
        b: nb,
        da,
        db,
    }
}

/// As [`lift_dense`], for sparse polynomials.
fn lift_sparse<N: Numerators>(
    ring: &SparseRing<N::Fractions>,
    numerators: &SparseRing<N>,
    split: &Split<Sparse<N::Elem>>,
    (g, qa, qb): Cofactors<Sparse<N::Elem>>,
) -> Cofactors<Sparse<Fraction<N>>>
where
    N::Fractions: Clone,
{
    let lift =
        |p: &Sparse<N::Elem>, den: &Natural| numerators.map_into(p, ring, |c| N::fraction(c, den));
    (
        lift(&g, &Natural::one()),
        lift(&qa, &split.da),
        lift(&qb, &split.db),
    )
}

/// ℚ and ℚ(i): the gcd of the numerators over ℤ and ℤ\[i\], by the same
/// method.
macro_rules! over_numerators {
    ($t:ty, $n:ty) => {
        impl GcdDomain for $t {
            fn dense_gcd(
                ring: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Cofactors<Dense<Self::Elem>>, Cancelled> {
                let numerators = DenseRing::new(<$n>::default());
                let split = split_dense(&numerators, a, b);
                let parts = <$n>::dense_gcd(&numerators, &split.a, &split.b, cancel)?;
                Ok(lift_dense::<$n>(ring, &split, parts))
            }

            fn dense_gcd_by(
                ring: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                method: GcdMethod,
                cancel: &dyn Cancel,
            ) -> Result<Option<Cofactors<Dense<Self::Elem>>>, Cancelled> {
                let numerators = DenseRing::new(<$n>::default());
                let split = split_dense(&numerators, a, b);
                let parts = <$n>::dense_gcd_by(&numerators, &split.a, &split.b, method, cancel)?;
                Ok(parts.map(|parts| lift_dense::<$n>(ring, &split, parts)))
            }

            fn sparse_gcd(
                ring: &SparseRing<Self>,
                a: &Sparse<Self::Elem>,
                b: &Sparse<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Cofactors<Sparse<Self::Elem>>, Cancelled> {
                let numerators = SparseRing::new(<$n>::default(), ring.nvars());
                let split = split_sparse_pair(ring, &numerators, a, b);
                let parts = <$n>::sparse_gcd(&numerators, &split.a, &split.b, cancel)?;
                Ok(lift_sparse::<$n>(ring, &numerators, &split, parts))
            }

            fn sparse_gcd_by(
                ring: &SparseRing<Self>,
                a: &Sparse<Self::Elem>,
                b: &Sparse<Self::Elem>,
                method: GcdMethod,
                cancel: &dyn Cancel,
            ) -> Result<Option<Cofactors<Sparse<Self::Elem>>>, Cancelled> {
                let numerators = SparseRing::new(<$n>::default(), ring.nvars());
                let split = split_sparse_pair(ring, &numerators, a, b);
                let parts = <$n>::sparse_gcd_by(&numerators, &split.a, &split.b, method, cancel)?;
                Ok(parts.map(|parts| lift_sparse::<$n>(ring, &numerators, &split, parts)))
            }

            /// `res(n / d, n' / d') = res(n, n') / (d^(deg n') d'^(deg n))`.
            fn dense_resultant(
                _: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Self::Elem, Cancelled> {
                let numerators = DenseRing::new(<$n>::default());
                let split = split_dense(&numerators, a, b);
                let res = <$n>::dense_resultant(&numerators, &split.a, &split.b, cancel)?;
                let degree = |p: &Dense<_>| p.degree().unwrap_or(0) as u64;
                let den = &split.da.pow(degree(&split.b)) * &split.db.pow(degree(&split.a));
                Ok(<$n>::fraction(&res, &den))
            }

            fn sparse_resultant(
                ring: &SparseRing<Self>,
                a: &Sparse<Self::Elem>,
                b: &Sparse<Self::Elem>,
                var: usize,
                cancel: &dyn Cancel,
            ) -> Result<Sparse<Self::Elem>, Cancelled> {
                let numerators = SparseRing::new(<$n>::default(), ring.nvars());
                let split = split_sparse_pair(ring, &numerators, a, b);
                let res = <$n>::sparse_resultant(&numerators, &split.a, &split.b, var, cancel)?;
                let degree = |p| numerators.degree(p, var).unwrap_or(0);
                let den = &split.da.pow(degree(&split.b)) * &split.db.pow(degree(&split.a));
                Ok(numerators.map_into(&res, ring, |c| <$n>::fraction(c, &den)))
            }
        }
    };
}

over_numerators!(Rationals, Integers);
over_numerators!(GaussianRationals, GaussianIntegers);

/// Prime fields: Euclid's algorithm in one variable, Brown's evaluation in
/// several, and subresultants where the field is too small for it.
macro_rules! prime_field {
    ($($t:ty),*) => {$(
        impl GcdDomain for $t {
            fn dense_gcd(
                ring: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Cofactors<Dense<Self::Elem>>, Cancelled> {
                euclid(ring, a, b, cancel)
            }

            fn dense_gcd_by(
                ring: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                method: GcdMethod,
                cancel: &dyn Cancel,
            ) -> Result<Option<Cofactors<Dense<Self::Elem>>>, Cancelled> {
                match method {
                    GcdMethod::Subresultant => Ok(Some(subres::dense_gcd(ring, a, b, cancel)?)),
                    GcdMethod::Heuristic
                    | GcdMethod::Modular
                    | GcdMethod::Zippel
                    | GcdMethod::EzGcd => Ok(None),
                }
            }

            fn sparse_gcd(
                ring: &SparseRing<Self>,
                a: &Sparse<Self::Elem>,
                b: &Sparse<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Cofactors<Sparse<Self::Elem>>, Cancelled> {
                match modular::pgcd(ring, a, b, cancel)? {
                    Some(g) => cofactors_of(ring, g, a, b, cancel),
                    None => subres::sparse_gcd(ring, a, b, cancel),
                }
            }

            fn sparse_gcd_by(
                ring: &SparseRing<Self>,
                a: &Sparse<Self::Elem>,
                b: &Sparse<Self::Elem>,
                method: GcdMethod,
                cancel: &dyn Cancel,
            ) -> Result<Option<Cofactors<Sparse<Self::Elem>>>, Cancelled> {
                match method {
                    GcdMethod::Modular => match modular::pgcd(ring, a, b, cancel)? {
                        Some(g) => Ok(Some(cofactors_of(ring, g, a, b, cancel)?)),
                        None => Ok(None),
                    },
                    GcdMethod::Subresultant => Ok(Some(subres::sparse_gcd(ring, a, b, cancel)?)),
                    GcdMethod::Zippel => zippel_cofactors(ring, a, b, cancel),
                    GcdMethod::Heuristic | GcdMethod::EzGcd => Ok(None),
                }
            }

            fn dense_resultant(
                ring: &DenseRing<Self>,
                a: &Dense<Self::Elem>,
                b: &Dense<Self::Elem>,
                cancel: &dyn Cancel,
            ) -> Result<Self::Elem, Cancelled> {
                resultant::field(ring, a, b, cancel)
            }
        }
    )*};
}

prime_field!(PrimeField<ModWord>, PrimeField<ModRing>);
