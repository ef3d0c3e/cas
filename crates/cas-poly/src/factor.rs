//! Factoring (CAS-DOC-8 §4, CAS-142, planned in CAS-PLAN-1): the
//! square-free decomposition; univariate factoring over ℤ and ℚ by
//! Zassenhaus's algorithm, from the factors of square-free polynomials over
//! word-size prime fields and their Hensel lifts to prime powers;
//! multivariate factoring over ℤ and ℚ by Wang's algorithm, with its
//! Hensel lifting; factoring over ℚ(i) by Trager's norms; and van Hoeij's
//! recombination by lattice reduction.
//! `doc/algorithms/factor.md` describes the algorithms.

mod finite;
mod hensel;
mod lifting;
mod lll;
mod squarefree;
mod trager;
mod vanhoeij;
mod wang;
mod zassenhaus;

#[cfg(test)]
mod tests;

use std::cell::Cell;
use std::fmt;

use cas_num::algebra::{GaussianIntegers, GaussianRationals, Integers, Rationals};

use crate::cancel::{Cancel, Cancelled};
use crate::gcd::GcdDomain;

/// Why a bounded factorisation stopped.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum FactorError {
    /// The [`Cancel`] said so.
    Cancelled,
    /// The searches would pass their bound: Zassenhaus's subsets of local
    /// factors, van Hoeij's columns, and Wang's points.
    Search,
}

impl fmt::Display for FactorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FactorError::Cancelled => f.write_str("cancelled"),
            FactorError::Search => f.write_str("too long a search"),
        }
    }
}

impl std::error::Error for FactorError {}

impl From<Cancelled> for FactorError {
    fn from(Cancelled: Cancelled) -> FactorError {
        FactorError::Cancelled
    }
}

/// The searches of a factorisation, counted against their bound, and what
/// the long loops poll.
pub(crate) struct Search<'a> {
    cancel: &'a dyn Cancel,
    left: Cell<u64>,
}

impl<'a> Search<'a> {
    pub(crate) fn new(cancel: &'a dyn Cancel, bound: u64) -> Search<'a> {
        Search {
            cancel,
            left: Cell::new(bound),
        }
    }

    /// The [`Cancel`] of the loops.
    pub(crate) fn cancel(&self) -> &'a dyn Cancel {
        self.cancel
    }

    /// One more subset or point, polling the [`Cancel`].
    pub(crate) fn step(&self) -> Result<(), FactorError> {
        if self.cancel.cancelled() {
            return Err(FactorError::Cancelled);
        }
        let left = self.left.get();
        if left == 0 {
            return Err(FactorError::Search);
        }
        self.left.set(left - 1);
        Ok(())
    }
}

/// The value of an unbounded search, which only a cancellation stops.
pub(crate) fn unbounded<T>(result: Result<T, FactorError>) -> Result<T, Cancelled> {
    match result {
        Ok(value) => Ok(value),
        Err(FactorError::Cancelled) => Err(Cancelled),
        Err(FactorError::Search) => unreachable!("an unbounded search"),
    }
}

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
