//! Elementary functions on balls (CAS-DOC-5 §7.2).
//!
//! Each function computes its value at the exact midpoint as a ball, at a
//! few guard bits above the requested precision, then widens it by a
//! derivative bound for the input radius. Only the series sums are proved
//! by hand ([`fixed`], [`bsplit`](super::bsplit)); argument reductions and
//! reconstructions run on ball arithmetic, which tracks their errors itself.
//!
//! After each function's own reduction (`exp` by `log 2`, `sin`/`cos` by
//! `π/2`, `atan` to `[0, 1]`), the working precision picks the core: lazy
//! tables ([`tables`]) up to [`TABLE_MAX_PREC`], bit-burst ([`bitburst`])
//! from [`BIT_BURST_MIN_PREC`], the baselines (halvings and a series)
//! between, for arguments the tables do not cover, and to build the tables.
//! `log` takes the tables directly; above them it stays "approximate, then
//! certify" on `exp`, so it follows `exp`'s cores.
#![allow(
    clippy::many_single_char_names,
    reason = "x, t, m, k, h, n, w, e as in the error analyses of the comments"
)]

mod atan;
mod bitburst;
mod exp;
mod fixed;
mod hyp;
mod log;
mod pow;
mod tables;
mod trig;

#[cfg(test)]
mod core_tests;
#[cfg(test)]
mod tests;

use super::{Ball, Mag};

pub(super) use tables::clear as clear_tables;

/// The highest working precision served by the tables, and the lowest by
/// bit-burst. Untuned placeholders (CAS-86): Johansson's 4600 bits and the
/// design's 10^4. Measuring them is cas-tune's job (CAS-PLAN-7 step #254,
/// deferred); until then they stay fixed here.
const TABLE_MAX_PREC: u64 = 4600;
const BIT_BURST_MIN_PREC: u64 = 10_000;

/// The working precision for a result at `prec` bits: guard bits for the
/// few ball operations after the core evaluation and its final rounding.
fn working(prec: u64) -> u64 {
    prec + 20 + u64::from(u64::BITS - prec.leading_zeros())
}

/// `x`'s exponent as an `i64`, saturated: an exponent beyond `i64` is far
/// outside every threshold the functions compare it with, on its side.
fn exponent(x: &super::Float) -> i64 {
    let e = x.exponent();
    e.to_i64()
        .unwrap_or(if e.is_negative() { i64::MIN } else { i64::MAX })
}

/// An upper bound of `|x|` over the ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// A lower bound of `n!`.
fn factorial_lower(n: u64) -> Mag {
    (2..=n).fold(Mag::ONE, |acc, k| {
        acc.mul_lower(&Mag::from_integer_lower(&crate::Integer::from(k)))
    })
}

/// The exact ball 1.
fn one() -> Ball {
    Ball::from_integer(&crate::Integer::ONE)
}
