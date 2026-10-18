//! Conversions with rounding: from `Rational`, to and from `f64`.

use std::cmp::Ordering;

use super::{Float, LIMB_BITS, Round};
use crate::kernel;
use crate::{Integer, Limb, Natural, Rational};

impl Float {
    /// `r` correctly rounded to `prec` bits, and `rounded.cmp(r)`.
    ///
    /// The quotient of the numerator, scaled to `prec + 2` bits or more, by
    /// the denominator, then rounded with the remainder as the sticky bit.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn from_rational(r: &Rational, prec: u64, rnd: Round) -> (Float, Ordering) {
        assert!(prec >= 1, "Float::from_rational: precision 0");
        if r.is_zero() {
            return (Float::ZERO, Ordering::Equal);
        }
        let (neg, num) = r.numer().clone().into_parts();
        let den = r.denom().unsigned_abs();
        // q = floor(num 2^s / den) has at least prec + 2 bits.
        let s = i128::from(prec + 2) + i128::from(den.bit_length()) - i128::from(num.bit_length());
        let (q, rem) = if s >= 0 {
            let s = u64::try_from(s).expect("shift fits");
            (&num << s).div_rem(&den)
        } else {
            let s = u64::try_from(-s).expect("shift fits");
            num.div_rem(&(&den << s))
        };
        let (mut m, bits) = aligned(&q);
        // value = (q + rem / den) 2^-s, and q has `bits` bits.
        let exp = Integer::from(i128::from(bits) - s);
        Float::finish_round(neg, &mut m, !rem.is_zero(), exp, prec, rnd)
    }

    /// The nearest `f64` (ties to even), with IEEE overflow to ±∞ and gradual
    /// underflow; zero is `+0.0`.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        if self.is_zero() {
            return 0.0;
        }
        let signed = |x: f64| if self.neg { -x } else { x };
        let Some(e) = self.exp.to_i64() else {
            return signed(if self.exp.is_negative() {
                0.0
            } else {
                f64::INFINITY
            });
        };
        if e > 1024 {
            return signed(f64::INFINITY);
        }
        // Below 2^-1075 the nearest double is 0.
        if e < -1074 {
            return signed(0.0);
        }
        // Normal doubles have 53 bits; below 2^-1022 the grid is 2^-1074.
        let prec = if e >= -1021 { 53 } else { e + 1074 };
        if prec == 0 {
            // |x| in [2^-1075, 2^-1074): exactly half rounds to even, 0.
            let half = self.mant().len() == 1 && self.mant()[0] == 1 << (Limb::BITS - 1);
            return signed(if half { 0.0 } else { f64::from_bits(1) });
        }
        let (r, _) = self.round(prec.unsigned_abs(), Round::Nearest);
        // At most one above `e` (a rounding carry), so in `i64`.
        let e = r.exp.to_i64().unwrap_or(i64::MAX);
        if e > 1024 {
            return signed(f64::INFINITY);
        }
        // r = top / 2^64 × 2^e with at most 53 bits in top.
        let top = r.mant()[r.mant().len() - 1];
        let m = top >> 11;
        debug_assert_eq!(m << 11, top);
        #[allow(clippy::cast_precision_loss, reason = "m has at most 53 bits")]
        signed(ldexp(m as f64, e - 53))
    }

    /// The exact value of a finite `f64`; `None` for infinities and NaN.
    /// Both zeros give zero.
    #[must_use]
    pub fn from_f64(x: f64) -> Option<Float> {
        if !x.is_finite() {
            return None;
        }
        let bits = x.to_bits();
        let neg = bits >> 63 == 1;
        #[allow(clippy::cast_possible_wrap, reason = "11 bits")]
        let biased = ((bits >> 52) & 0x7ff) as i64;
        let frac = bits & ((1 << 52) - 1);
        let (m, e) = if biased == 0 {
            (frac, -1074)
        } else {
            (frac | 1 << 52, biased - 1075)
        };
        Some(Float::from_raw(neg, &[m], Integer::from(e)))
    }
}

/// `q` top-aligned (shifted left until its top bit is set) and its bit length.
fn aligned(q: &Natural) -> (Vec<Limb>, u64) {
    let mut m = q.limbs().to_vec();
    let lz = m[m.len() - 1].leading_zeros();
    let out = kernel::shl_assign(&mut m, lz);
    debug_assert_eq!(out, 0);
    let bits = LIMB_BITS * m.len() as u64 - u64::from(lz);
    (m, bits)
}

/// `x 2^k`, exact when the result is a finite double: the scaling goes in
/// steps of normal powers of two, through values that are all multiples of
/// the result's grid.
fn ldexp(mut x: f64, mut k: i64) -> f64 {
    const STEP: i64 = 1000;
    let pow2 = |k: i64| f64::from_bits(u64::try_from(k + 1023).expect("normal exponent") << 52);
    while k > STEP {
        x *= pow2(STEP);
        k -= STEP;
    }
    while k < -STEP {
        x *= pow2(-STEP);
        k += STEP;
    }
    x * pow2(k)
}
