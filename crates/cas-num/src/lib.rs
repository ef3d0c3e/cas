//! Arbitrary-precision numerics for the CAS engine.
//!
//! Correctness is the primary goal: every operation is differentially tested
//! against GMP (test-only dependency). See the "Numerics design" page in Lific.
//!
//! - [`kernel`]: allocation-free limb kernels on slices, the building blocks of
//!   every algorithm.
//! - [`Natural`]: unsigned arbitrary-precision integers.
//! - [`Integer`]: signed arbitrary-precision integers; values that fit an `i64`
//!   never allocate.
//! - [`Rational`]: rational numbers, always reduced with a positive
//!   denominator.
//! - [`modular`]: residues modulo word-size and multi-precision moduli,
//!   Chinese remaindering and rational reconstruction.

mod error;
mod integer;
pub mod kernel;
pub mod modular;
mod natural;
mod rational;

pub use error::{OutOfRangeError, ParseErrorKind, ParseNumberError};
pub use integer::{Integer, Rounding};
pub use natural::Natural;
pub use rational::Rational;

/// A machine word of a multi-precision number. Numbers are little-endian
/// sequences of limbs.
pub type Limb = u64;

#[cfg(test)]
mod testkit;
