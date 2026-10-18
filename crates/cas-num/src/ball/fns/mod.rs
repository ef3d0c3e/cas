//! Elementary functions on balls (CAS-DOC-5 §7.2): the baselines.
//!
//! Each function computes its value at the exact midpoint as a ball, at a
//! few guard bits above the requested precision, then widens it by a
//! derivative bound for the input radius. Only the series sums are proved
//! by hand ([`fixed`]); argument reductions and reconstructions run on ball
//! arithmetic, which tracks their errors itself.
#![allow(
    clippy::many_single_char_names,
    reason = "x, t, m, k, h, n, w, e as in the error analyses of the comments"
)]

mod exp;
mod fixed;
mod log;

#[cfg(test)]
mod tests;

use super::{Ball, Mag};

/// The working precision for a result at `prec` bits: guard bits for the
/// few ball operations after the core evaluation and its final rounding.
fn working(prec: u64) -> u64 {
    prec + 20 + u64::from(u64::BITS - prec.leading_zeros())
}

/// An upper bound of `|x|` over the ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// The exact ball 1.
fn one() -> Ball {
    Ball::from_integer(&crate::Integer::ONE)
}
