//! Factoring (CAS-DOC-8 §4, CAS-142, planned in CAS-PLAN-1): the
//! square-free decomposition, the factors of square-free polynomials over
//! word-size prime fields, and their Hensel lifts to prime powers.
//! `doc/algorithms/factor.md` describes the algorithms.

mod finite;
mod hensel;
mod squarefree;

#[cfg(test)]
mod tests;

use cas_num::algebra::{GaussianIntegers, GaussianRationals, Integers, Rationals};

use crate::gcd::GcdDomain;

/// A polynomial as a constant times powers of polynomials, `c f₁^e₁ ⋯
/// fₖ^eₖ`. Each `fᵢ` has a positive degree, is primitive, and has its
/// leading coefficient in the normal form of the coefficients: positive over
/// ℤ, monic over a field, in the first quadrant over ℤ\[i\]. For sparse
/// polynomials that is the lexicographic leading coefficient. The exponents
/// are at least 1.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Factors<C, P> {
    constant: C,
    factors: Vec<(P, u64)>,
}

impl<C, P> Factors<C, P> {
    /// The constant `c`. Over ℤ it is the content with the sign of the
    /// leading coefficient, and over a field the leading coefficient. It is 0
    /// for the zero polynomial.
    #[must_use]
    pub fn constant(&self) -> &C {
        &self.constant
    }

    /// The factors with their exponents.
    #[must_use]
    pub fn factors(&self) -> &[(P, u64)] {
        &self.factors
    }

    #[must_use]
    pub fn into_parts(self) -> (C, Vec<(P, u64)>) {
        (self.constant, self.factors)
    }
}

/// The [`GcdDomain`]s of characteristic 0: ℤ, ℚ, ℤ\[i\] and ℚ(i). Over these
/// domains a factor that a polynomial has more than once divides its
/// derivative, and square-free decomposition relies on that. The prime
/// fields do not implement it.
pub trait CharacteristicZero: GcdDomain {}

impl CharacteristicZero for Integers {}

impl CharacteristicZero for Rationals {}

impl CharacteristicZero for GaussianIntegers {}

impl CharacteristicZero for GaussianRationals {}
