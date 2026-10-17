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

mod error;
mod integer;
pub mod kernel;
mod natural;

pub use error::OutOfRangeError;
pub use integer::{Integer, Rounding};
pub use natural::Natural;

/// A machine word of a multi-precision number. Numbers are little-endian
/// sequences of limbs.
pub type Limb = u64;

#[cfg(test)]
mod testkit;
