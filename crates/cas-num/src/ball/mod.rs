//! Ball arithmetic (CAS-DOC-5): arbitrary-precision binary floats as
//! midpoints, and the types built on them.
//!
//! - [`Float`]: an exact binary float `±0.m × 2^exp` with an `Integer`
//!   exponent, exact operations, and correctly rounded ones in every
//!   [`Round`] mode (§3).

mod float;

pub use float::{Float, Round};
