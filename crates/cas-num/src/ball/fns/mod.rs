//! Elementary functions on balls (CAS-DOC-5 §7.2).
//!
//! Each function computes its value at the exact midpoint as a ball, at a few
//! guard bits above the requested precision, then widens it by a derivative
//! bound for the input radius, as Arb does \[Joh17\]. Only the series sums are
//! proved by hand ([`fixed`], [`bsplit`](super::bsplit)); argument reductions
//! and reconstructions run on ball arithmetic, which tracks their errors
//! itself.
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
mod erf;
mod exp;
mod expint;
mod fixed;
mod gamma;
mod hyp;
mod lambert;
mod log;
mod polygamma;
mod pow;
mod tables;
mod trig;
mod zeta;

#[cfg(test)]
mod core_tests;
#[cfg(test)]
mod tests;

use super::{Ball, Mag};

pub use gamma::bernoulli;
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

/// `f` at the working precision of `prec` bits, and for an exact `b` again
/// at more while the relative accuracy falls short of `prec`, up to `4 (wp +
/// bits(b)) + 256` more: next to a zero of a function that is not known
/// exactly, a fixed precision is only absolute (CAS-DOC-31 §2).
fn refined(b: &Ball, prec: u64, f: impl Fn(u64) -> Ball) -> Ball {
    let wp = working(prec);
    let mut value = f(wp);
    if !b.is_exact() {
        return value;
    }
    let want = i64::try_from(prec).unwrap_or(i64::MAX).saturating_add(2);
    let cap = 4 * (wp + b.mid().bit_length()) + 256;
    let mut extra = 0u64;
    while value.is_finite() && value.rel_accuracy_bits() < want && extra < cap {
        let known = u64::try_from(value.rel_accuracy_bits().max(0)).unwrap_or(0);
        let short = u64::try_from(want)
            .unwrap_or(u64::MAX)
            .saturating_sub(known);
        extra = (2 * extra).max(short + 32).min(cap);
        value = f(wp + extra);
    }
    value
}
