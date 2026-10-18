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
//! - Constants: [`Ball::const_pi`], [`Ball::const_log2`] and [`Ball::const_e`]
//!   by binary splitting, cached process-wide (§7.1).

mod bsplit;
mod consts;
mod float;
mod mag;
mod real;

pub use consts::clear_constant_cache;
pub use float::{APPROX_ULPS, Float, Round};
pub use mag::Mag;
pub use real::{Ball, RealBalls};
