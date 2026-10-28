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
//! - Constants: [`Ball::const_pi`], [`Ball::const_log2`], [`Ball::const_e`],
//!   [`Ball::const_euler`] and [`Ball::const_catalan`] by binary splitting,
//!   cached process-wide (§7.1).
//! - Elementary functions: [`Ball::exp`], [`Ball::expm1`], [`Ball::log`],
//!   [`Ball::log1p`], [`Ball::sin_cos`], [`Ball::tan`], [`Ball::sin_pi`],
//!   [`Ball::atan`], [`Ball::asin`], [`Ball::acos`], [`Ball::atan2`] (§7.2);
//!   [`Ball::pow_ball`], [`Ball::root`], [`Ball::sinh`] and the other
//!   hyperbolic functions (§7.3).
//! - [`ComplexBall`]: rectangular complex balls with principal branches
//!   (§8).
//! - Decimal printing with Arb's digit rule ([`Ball::to_decimal`], `Display`)
//!   and decimal parsing ([`Ball::from_decimal`]); exact hexadecimal floats
//!   for `FullForm` ([`Float::to_hex`], [`Mag::to_hex`]) (§6).
//! - [`MachineBall`]: `f64` midpoint and radius with proven bounds, the fast
//!   path up to 53 bits; `None` sends the caller back to [`Ball`] (§5.4).
//! - [`polynomial_roots`]: the roots of a square-free polynomial with exact
//!   coefficients, each in a validated ball (CAS-211).
//! - [`gauss_legendre`]: the nodes and weights of Gauss–Legendre rules in
//!   validated balls, cached (CAS-DOC-26 §2).
//!
//! The algorithms, the tricks on top of them and their sources are
//! documented in this crate's `doc/algorithms/`, after GMP's "Algorithms"
//! chapter and MPFR's *Algorithms and Proofs*. Keys in brackets in these doc
//! comments, such as \[BZ10\] or \[Joh17\], refer to its `references.md`.

mod bsplit;
mod complex;
mod consts;
mod decimal;
mod float;
mod fns;
mod hex;
#[cfg(test)]
mod huge_tests;
mod machine;
mod mag;
mod quadrature;
mod real;
mod roots;

pub use complex::ComplexBall;
pub use consts::clear_constant_cache;
pub use decimal::DecimalStyle;
pub use float::{APPROX_ULPS, EXACT_BITS_MAX, Float, Round};
pub use fns::bernoulli;
pub use machine::MachineBall;
pub use mag::Mag;
pub use quadrature::{MAX_CACHED, Rule, gauss_legendre};
pub use real::{Ball, RealBalls};
pub use roots::polynomial_roots;
