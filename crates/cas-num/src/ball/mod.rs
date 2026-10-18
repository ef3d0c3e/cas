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
//! - Elementary functions: [`Ball::exp`], [`Ball::expm1`], [`Ball::log`],
//!   [`Ball::log1p`] (§7.2).
//! - Decimal printing with Arb's digit rule ([`Ball::to_decimal`], `Display`)
//!   and decimal parsing ([`Ball::from_decimal`]); exact hexadecimal floats
//!   for `FullForm` ([`Float::to_hex`], [`Mag::to_hex`]) (§6).
//! - [`MachineBall`]: `f64` midpoint and radius with proven bounds, the fast
//!   path up to 53 bits; `None` sends the caller back to [`Ball`] (§5.4).

mod bsplit;
mod consts;
mod decimal;
mod float;
mod fns;
mod hex;
mod machine;
mod mag;
mod real;

pub use consts::clear_constant_cache;
pub use decimal::DecimalStyle;
pub use float::{APPROX_ULPS, Float, Round};
pub use machine::MachineBall;
pub use mag::Mag;
pub use real::{Ball, RealBalls};
