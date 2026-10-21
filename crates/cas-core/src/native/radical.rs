//! Rational powers of rationals (CAS-DOC-4 §6.1, O4).
//!
//! The canonical form of `b^r`, for a rational `b` and a rational `r` that
//! is not an integer, is a coefficient times powers of positive integers with
//! exponents in `(0, 1)`, one power per exponent:
//!
//! - the integer part of the exponent comes out: `2^(3/2)` is `2·2^(1/2)`,
//!   and `2^(-1/2)` is `2^(1/2)/2`;
//! - perfect powers come out of each base: `12^(1/2)` is `2·3^(1/2)`,
//!   `8^(2/3)` is `4`, `72^(1/3)` is `2·3^(2/3)`. The bases are factored by
//!   trial division below 2¹⁶ and perfect-power tests only, so a huge
//!   radicand never starts a factorisation;
//! - a denominator goes into the numerator, `(2/3)^(1/2)` is `6^(1/2)/3`,
//!   since `D^(-f)` is `D^(1 - f)/D`;
//! - the square root of a negative number is imaginary, `(-4)^(1/2)` is
//!   `2·i`. In the Real domain an odd root of a negative is the real root
//!   (O7): `(-8)^(1/3)` is `-2` and `(-2)^(2/3)` is `2^(2/3)`. Other roots of
//!   negatives, and odd ones in the Complex domain, stay as they are.

use std::collections::BTreeMap;

use cas_num::{Integer, Natural, Rational};

use super::number::{EXACT_ONLY, Num};
use crate::expr::Expr;
use crate::symbol::Symbol;

/// Bases with more bits are not searched for perfect powers.
const MAX_FACTORED_BITS: u64 = 4096;

/// Exponents with a larger denominator are left alone.
const MAX_ROOT: u64 = 1 << 16;

/// The canonical form of `b^r` for rationals `b != 0` and `r` not an
/// integer, or `None` when it is `b^r` itself or cannot be computed, as when
/// its integer part would have more than `max_bits` bits by estimate;
/// `real_domain` takes odd roots of negatives as real ones.
pub(super) fn power(b: &Rational, r: &Rational, real_domain: bool, max_bits: u64) -> Option<Expr> {
    let (whole, numer, denom) = split_exponent(r)?;
    if denom > MAX_ROOT {
        return None;
    }
    // The sign of the root of `b`, the integer part of the exponent taking
    // its own: `(-1)^(p/2)` is `i^p`, and in the Real domain `(-1)^(p/q)` is
    // `(-1)^p` for an odd `q`.
    let mut coefficient = Num::one();
    if b.is_negative() {
        if denom == 2 {
            let i = Num::of(&crate::build::i()).expect("a number");
            coefficient = i.pow(i64::try_from(numer % 4).expect("small"), EXACT_ONLY)?;
        } else if real_domain && denom % 2 == 1 {
            if numer % 2 == 1 {
                coefficient = Num::Real(Rational::from(-1));
            }
        } else {
            return None;
        }
    }
    let magnitude = b.abs();
    // The integer part of the exponent.
    let whole_power = Num::Real(b.clone());
    if whole_power.power_bits(whole) > max_bits {
        return None;
    }
    coefficient = coefficient.mul(&whole_power.pow(whole, EXACT_ONLY)?, EXACT_ONLY);
    // The numerator to the power `numer/denom`, and the denominator `D` as
    // `D^((denom - numer)/denom) / D`.
    let mut roots = Roots::default();
    let (n, d) = (
        magnitude.numer().unsigned_abs(),
        magnitude.denom().unsigned_abs(),
    );
    let out_n = roots.extract(&n, numer, denom);
    let out_d = roots.extract(&d, denom - numer, denom);
    coefficient = coefficient.mul(
        &Num::Real(Rational::new(
            Integer::from(out_n) * Integer::from(out_d),
            Integer::from(d),
        )),
        EXACT_ONLY,
    );
    let mut factors = Vec::new();
    if !coefficient.is_one() {
        factors.push(coefficient.to_expr());
    }
    for (exponent, base) in roots.0 {
        if base != Natural::from(1u64) {
            factors.push(Expr::apply(
                &Symbol::POWER.into(),
                &[Expr::from(Integer::from(base)), Expr::from(exponent)],
            ));
        }
    }
    Some(match factors.len() {
        0 => Expr::from(1),
        1 => factors.pop().expect("one factor"),
        _ => Expr::apply(&Symbol::TIMES.into(), &factors),
    })
}

/// `r` as `whole + numer/denom` with `0 < numer < denom`, if `r` is not an
/// integer and its parts fit.
fn split_exponent(r: &Rational) -> Option<(i64, u64, u64)> {
    let (q, rem) = r.numer().div_rem_floor(r.denom());
    let whole = q.to_i64()?;
    let numer = u64::try_from(rem.to_i64()?).ok()?;
    let denom = u64::try_from(r.denom().to_i64()?).ok()?;
    (numer > 0).then_some((whole, numer, denom))
}

/// Radicals by exponent: the product of the bases raised to each exponent.
#[derive(Default)]
struct Roots(BTreeMap<Rational, Natural>);

impl Roots {
    /// Takes the perfect powers out of `n^(p/q)`: returns what comes out,
    /// and keeps the radicals.
    fn extract(&mut self, n: &Natural, p: u64, q: u64) -> Natural {
        let mut out = Natural::from(1u64);
        if n.bit_length() > MAX_FACTORED_BITS {
            self.add(Rational::new(Integer::from(p), Integer::from(q)), n);
            return out;
        }
        if *n == Natural::from(1u64) {
            return out;
        }
        for f in n.factor_with_effort(0).factors() {
            let total = f.exponent() * p;
            let (whole, rest) = (total / q, total % q);
            out = &out * &f.base().pow(whole);
            if rest > 0 {
                self.add(
                    Rational::new(Integer::from(rest), Integer::from(q)),
                    f.base(),
                );
            }
        }
        out
    }

    fn add(&mut self, exponent: Rational, base: &Natural) {
        let product = self
            .0
            .entry(exponent)
            .or_insert_with(|| Natural::from(1u64));
        *product = &*product * base;
    }
}
