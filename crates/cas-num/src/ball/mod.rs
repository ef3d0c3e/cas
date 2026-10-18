//! Ball arithmetic (CAS-DOC-5): arbitrary-precision binary floats as
//! midpoints, and the types built on them.
//!
//! - [`Float`]: an exact binary float `±0.m × 2^exp` with an `Integer`
//!   exponent, exact operations, and correctly rounded ones in every
//!   [`Round`] mode (§3).
//! - [`Mag`]: a 30-bit upper (or lower) bound on a magnitude with an
//!   `Integer` exponent, the radius type (§4).
//! - [`Ball`]: midpoint ± radius, with arithmetic whose results contain every
//!   exact result for points of the inputs, three-valued predicates, and the
//!   [`RealBalls`] ring object (§5).

mod float;
mod mag;
mod real;

pub use float::{APPROX_ULPS, Float, Round};
pub use mag::Mag;
pub use real::{Ball, RealBalls};
