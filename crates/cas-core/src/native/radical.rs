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
//!   `2·i`. Other roots of negative numbers stay as they are until the number
//!   domain (CAS-58) decides between the real root and the principal one.

use std::collections::BTreeMap;

use cas_num::{Integer, Natural, Rational};

use super::number::Num;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// Bases with more bits are not searched for perfect powers.
const MAX_FACTORED_BITS: u64 = 4096;

/// Exponents with a larger denominator are left alone.
const MAX_ROOT: u64 = 1 << 16;

/// Results with more bits than this, by estimate, are left unevaluated.
pub(super) const MAX_POWER_BITS: u64 = 1 << 24;

/// The canonical form of `b^r` for rationals `b != 0` and `r` not an
/// integer, or `None` when it is `b^r` itself or cannot be computed.
pub(super) fn power(b: &Rational, r: &Rational) -> Option<Expr> {
    let (whole, numer, denom) = split_exponent(r)?;
    if denom > MAX_ROOT {
        return None;
    }
    // The sign: `(-1)^(p/2)` is `i^p`; other roots of negatives wait.
    let mut coefficient = Num::one();
    if b.is_negative() {
        if denom != 2 {
            return None;
        }
        let i = Num::of(&crate::build::i()).expect("a number");
        coefficient = i.pow(i64::try_from(numer % 4).expect("small"))?;
    }
    let magnitude = b.abs();
    // The integer part of the exponent.
    let whole_power = Num::Real(b.clone());
    if whole_power.power_bits(whole) > MAX_POWER_BITS {
        return None;
    }
    coefficient = coefficient.mul(&whole_power.pow(whole)?);
    // The numerator to the power `numer/denom`, and the denominator `D` as
    // `D^((denom - numer)/denom) / D`.
    let mut roots = Roots::default();
    let (n, d) = (
        magnitude.numer().unsigned_abs(),
        magnitude.denom().unsigned_abs(),
    );
    let out_n = roots.extract(&n, numer, denom);
    let out_d = roots.extract(&d, denom - numer, denom);
    coefficient = coefficient.mul(&Num::Real(Rational::new(
        Integer::from(out_n) * Integer::from(out_d),
        Integer::from(d),
    )));
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
