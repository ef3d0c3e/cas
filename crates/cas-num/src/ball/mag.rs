//! `Mag`: a fixed-size bound on a magnitude, the radius type of balls
//! (CAS-DOC-5 §4, Arb's `mag_t`).
//!
//! A value is `man × 2^(exp - 30)` with a 30-bit mantissa `man` in
//! `[2^29, 2^30)` and an [`Integer`] exponent, or zero, or `+∞`. Every
//! operation returns an **upper bound** of the exact result (the `*_lower`
//! forms a lower bound), off by at most `2^-29` relative for one rounding:
//! two 30-bit mantissas multiply exactly in a `u64`, then one shift rounds
//! the product up. In the `i64` exponent range nothing allocates.
//!
//! Saturation can only make a bound looser, never wrong: an upper bound that
//! cannot be represented is `+∞`, and a lower bound of something unknown is 0.
#![allow(
    clippy::many_single_char_names,
    reason = "m, e, k, t, v, w: mantissas, exponents and the series' textbook names"
)]

use std::cmp::Ordering;
use std::fmt;

use super::Float;
use crate::{Integer, Limb};

#[cfg(test)]
mod tests;

/// Mantissa bits.
const BITS: u32 = 30;
/// The smallest normalised mantissa, `2^29`.
const LOW: u32 = 1 << (BITS - 1);
/// One past the largest, `2^30`.
const HIGH: u32 = 1 << BITS;
/// The mantissa that marks `+∞` (with exponent 0).
const INF_MAN: u32 = u32::MAX;

/// `e + d`, through `i64` when it fits.
#[inline]
fn offset(e: &Integer, d: i64) -> Integer {
    match e.to_i64().and_then(|v| v.checked_add(d)) {
        Some(v) => Integer::from(v),
        None => e + &Integer::from(d),
    }
}

/// `a + b + d` (or `a - b + d` with `minus`), through `i64` when it fits.
#[inline]
fn combine(a: &Integer, b: &Integer, minus: bool, d: i64) -> Integer {
    let small = a.to_i64().zip(b.to_i64()).and_then(|(x, y)| {
        let s = if minus {
            x.checked_sub(y)
        } else {
            x.checked_add(y)
        };
        s?.checked_add(d)
    });
    match small {
        Some(v) => Integer::from(v),
        None if minus => a - b + Integer::from(d),
        None => a + b + Integer::from(d),
    }
}

/// `a - b` as a `u32` shift, if it is one.
#[inline]
fn gap(a: &Integer, b: &Integer) -> Option<u32> {
    let d = match a.to_i64().zip(b.to_i64()) {
        Some((x, y)) => x.checked_sub(y)?,
        None => (a - b).to_i64()?,
    };
    u32::try_from(d).ok()
}

/// An upper (or lower) bound on a non-negative real number.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Mag {
    man: u32,
    exp: Integer,
}

impl Mag {
    /// Zero.
    pub const ZERO: Mag = Mag {
        man: 0,
        exp: Integer::ZERO,
    };

    /// `+∞`.
    pub const INF: Mag = Mag {
        man: INF_MAN,
        exp: Integer::ZERO,
    };

    /// One.
    pub const ONE: Mag = Mag {
        man: LOW,
        exp: Integer::ONE,
    };

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.man == 0
    }

    #[must_use]
    pub fn is_inf(&self) -> bool {
        self.man == INF_MAN
    }

    #[must_use]
    pub fn is_finite(&self) -> bool {
        !self.is_inf()
    }

    /// The exponent `e` with `2^(e - 1) <= self < 2^e` for a finite non-zero
    /// value; 0 for zero and `+∞`.
    #[must_use]
    pub fn exponent(&self) -> &Integer {
        &self.exp
    }

    /// `m × 2^(exp - 30)` for any `m` (`0 < m < 2^63`), normalised and rounded
    /// up (or down) to 30 bits.
    #[allow(
        clippy::needless_pass_by_value,
        reason = "callers build the exponent as a temporary"
    )]
    fn from_u64(m: u64, exp: Integer, up: bool) -> Mag {
        debug_assert!(m > 0 && m >> 63 == 0);
        let bits = u64::BITS - m.leading_zeros();
        if bits <= BITS {
            // Exact: shift left.
            let shift = BITS - bits;
            #[allow(clippy::cast_possible_truncation, reason = "at most 30 bits")]
            let man = (m << shift) as u32;
            return Mag::new(man, offset(&exp, -i64::from(shift)));
        }
        let shift = bits - BITS;
        #[allow(clippy::cast_possible_truncation, reason = "30 bits")]
        let mut man = (m >> shift) as u32;
        let lost = m & ((1 << shift) - 1) != 0;
        let mut adj = i64::from(shift);
        if up && lost {
            man += 1;
            if man == HIGH {
                man = LOW;
                adj += 1;
            }
        }
        Mag::new(man, offset(&exp, adj))
    }

    /// Wraps a normalised mantissa.
    fn new(man: u32, exp: Integer) -> Mag {
        debug_assert!(
            (LOW..HIGH).contains(&man),
            "unnormalised Mag mantissa {man:#x}"
        );
        Mag { man, exp }
    }

    /// `|x|`, rounded up.
    #[must_use]
    pub fn from_float(x: &Float) -> Mag {
        Mag::from_float_dir(x, true)
    }

    /// `|x|`, rounded down.
    #[must_use]
    pub fn from_float_lower(x: &Float) -> Mag {
        Mag::from_float_dir(x, false)
    }

    fn from_float_dir(x: &Float, up: bool) -> Mag {
        let m = x.mant();
        let Some(&top) = m.last() else {
            return Mag::ZERO;
        };
        // The top 30 bits of the top-aligned mantissa; the rest decides
        // whether rounding up adds one.
        let shift = Limb::BITS - BITS;
        #[allow(clippy::cast_possible_truncation, reason = "30 bits")]
        let mut man = (top >> shift) as u32;
        let lost = top & ((1 << shift) - 1) != 0 || m.len() > 1;
        let mut exp = x.exponent().clone();
        if up && lost {
            man += 1;
            if man == HIGH {
                man = LOW;
                exp += &Integer::ONE;
            }
        }
        Mag::new(man, exp)
    }

    /// `|n|`, rounded up.
    #[must_use]
    pub fn from_integer(n: &Integer) -> Mag {
        Mag::from_float(&Float::from(n))
    }

    /// `|n|`, rounded down.
    #[must_use]
    pub fn from_integer_lower(n: &Integer) -> Mag {
        Mag::from_float_lower(&Float::from(n))
    }

    /// `2^k`, exactly.
    #[must_use]
    pub fn power_of_two(k: &Integer) -> Mag {
        Mag::new(LOW, k + &Integer::ONE)
    }

    /// The exact value as a `Float`; `None` for `+∞`.
    #[must_use]
    pub fn to_float(&self) -> Option<Float> {
        if self.is_inf() {
            return None;
        }
        if self.is_zero() {
            return Some(Float::ZERO);
        }
        Some(Float::from_integer_exp(
            &Integer::from(i64::from(self.man)),
            &(&self.exp - &Integer::from(BITS)),
        ))
    }

    /// `self × 2^k`, exactly (`+∞` and zero unchanged).
    #[must_use]
    pub fn mul_2exp(&self, k: &Integer) -> Mag {
        if self.is_zero() || self.is_inf() {
            return self.clone();
        }
        Mag::new(self.man, &self.exp + k)
    }

    /// An upper bound of `self + rhs`.
    #[must_use]
    pub fn add(&self, rhs: &Mag) -> Mag {
        self.add_dir(rhs, true)
    }

    /// A lower bound of `self + rhs`.
    #[must_use]
    pub fn add_lower(&self, rhs: &Mag) -> Mag {
        self.add_dir(rhs, false)
    }

    /// `self + rhs` rounded up or down: the smaller operand shifted down to the
    /// larger's scale, its lost bits rounding the sum (at most one unit of a
    /// 31-bit sum, then the final rounding).
    fn add_dir(&self, rhs: &Mag, up: bool) -> Mag {
        if self.is_inf() || rhs.is_inf() {
            return Mag::INF;
        }
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return rhs.clone();
        }
        let (x, y) = if self.exp >= rhs.exp {
            (self, rhs)
        } else {
            (rhs, self)
        };
        // Both at 2 guard bits: x m 2^2, y shifted by the exponent gap.
        let gap = gap(&x.exp, &y.exp);
        let xm = u64::from(x.man) << 2;
        let ym = u64::from(y.man) << 2;
        let (ys, lost) = match gap {
            Some(d) if d < 64 => {
                let ys = ym >> d;
                (ys, ys << d != ym)
            }
            _ => (0, true),
        };
        let mut s = xm + ys;
        if up && lost {
            s += 1;
        }
        Mag::from_u64(s, offset(&x.exp, -2), up)
    }

    /// An upper bound of `max(self - rhs, 0)`.
    #[must_use]
    pub fn sub(&self, rhs: &Mag) -> Mag {
        self.sub_dir(rhs, true)
    }

    /// A lower bound of `max(self - rhs, 0)` (for the `|b| - s` of a
    /// divisor's radius).
    #[must_use]
    pub fn sub_lower(&self, rhs: &Mag) -> Mag {
        self.sub_dir(rhs, false)
    }

    fn sub_dir(&self, rhs: &Mag, up: bool) -> Mag {
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_inf() {
            // ∞ - finite is ∞; ∞ - ∞ is unknown: ∞ above, 0 below.
            return if rhs.is_inf() && !up {
                Mag::ZERO
            } else {
                Mag::INF
            };
        }
        if rhs.is_inf() || self <= rhs {
            return Mag::ZERO;
        }
        // self > rhs > 0 and self.exp >= rhs.exp.
        let gap = gap(&self.exp, &rhs.exp);
        let xm = u64::from(self.man) << 32;
        let ym = u64::from(rhs.man) << 32;
        let d = match gap {
            Some(d) if d < 63 => {
                let ys = ym >> d;
                let lost = ys << d != ym;
                // y's lost bits make the difference smaller: round the
                // subtracted part up for a lower bound, down for an upper.
                xm - ys - u64::from(lost && !up)
            }
            // y is below every bit of x: x itself (above), or just below it.
            _ => xm - u64::from(!up),
        };
        if d == 0 {
            return Mag::ZERO;
        }
        // d 2^(exp - 62), as man 2^(E - 30) form: E = exp - 32.
        Mag::from_u64(d, offset(&self.exp, -32), up)
    }

    /// The next value up (one unit in the last place more).
    fn next_up(&self) -> Mag {
        if self.is_zero() || self.is_inf() {
            return self.clone();
        }
        let man = self.man + 1;
        if man == HIGH {
            Mag::new(LOW, &self.exp + &Integer::ONE)
        } else {
            Mag::new(man, self.exp.clone())
        }
    }

    /// An upper bound of `self × rhs`.
    #[must_use]
    pub fn mul(&self, rhs: &Mag) -> Mag {
        self.mul_dir(rhs, true)
    }

    /// A lower bound of `self × rhs`.
    #[must_use]
    pub fn mul_lower(&self, rhs: &Mag) -> Mag {
        self.mul_dir(rhs, false)
    }

    fn mul_dir(&self, rhs: &Mag, up: bool) -> Mag {
        if self.is_zero() || rhs.is_zero() {
            // 0 × ∞ only arises from a radius that is both: 0 is right for an
            // exact zero factor.
            return Mag::ZERO;
        }
        if self.is_inf() || rhs.is_inf() {
            return Mag::INF;
        }
        // man_a man_b < 2^60, the value (man_a man_b) 2^(ea + eb - 60).
        let p = u64::from(self.man) * u64::from(rhs.man);
        Mag::from_u64(p, combine(&self.exp, &rhs.exp, false, -i64::from(BITS)), up)
    }

    /// An upper bound of `self / rhs` (`+∞` for a zero divisor).
    #[must_use]
    pub fn div(&self, rhs: &Mag) -> Mag {
        self.div_dir(rhs, true)
    }

    /// A lower bound of `self / rhs` (0 for an infinite divisor).
    #[must_use]
    pub fn div_lower(&self, rhs: &Mag) -> Mag {
        self.div_dir(rhs, false)
    }

    fn div_dir(&self, rhs: &Mag, up: bool) -> Mag {
        if self.is_zero() {
            return Mag::ZERO;
        }
        if rhs.is_inf() {
            // finite / ∞ = 0; ∞ / ∞ is unknown.
            return if self.is_inf() && up {
                Mag::INF
            } else {
                Mag::ZERO
            };
        }
        if rhs.is_zero() || self.is_inf() {
            return Mag::INF;
        }
        // (man_a 2^k) / man_b in [2^29, 2^30), one rounding: k = 29 when
        // man_a >= man_b (ratio in [1, 2)), else 30.
        let k = if self.man >= rhs.man { BITS - 1 } else { BITS };
        let num = u64::from(self.man) << k;
        let den = u64::from(rhs.man);
        let q = num / den;
        let exact = q * den == num;
        // a / b = (q / 2^k) 2^(ea - eb) and man 2^(e - 30) form: e = ea - eb - k + 30.
        let exp = combine(&self.exp, &rhs.exp, true, i64::from(BITS) - i64::from(k));
        let m = Mag::from_u64(q, exp, up);
        if up && !exact { m.next_up() } else { m }
    }

    /// An upper bound of `self^n` (`self^0 = 1`), by binary exponentiation,
    /// within a factor `(1 + 2^-29)^(2n)` of the exact power: a squaring
    /// doubles the relative error it receives, so a rounding made at `a^m`
    /// ends up raised to `n / m`.
    #[must_use]
    pub fn pow(&self, n: u64) -> Mag {
        let mut r = Mag::ONE;
        for i in (0..u64::BITS - n.leading_zeros()).rev() {
            r = r.mul(&r);
            if (n >> i) & 1 == 1 {
                r = r.mul(self);
            }
        }
        r
    }

    /// An upper bound of `sqrt(self)`.
    #[must_use]
    pub fn sqrt(&self) -> Mag {
        self.sqrt_dir(true)
    }

    /// A lower bound of `sqrt(self)`.
    #[must_use]
    pub fn sqrt_lower(&self) -> Mag {
        self.sqrt_dir(false)
    }

    /// With `self = man 2^E` (`E = exp - 30`), make `E` even by doubling `man`
    /// if needed, then `s = isqrt(m 2^32)` (31 or 32 bits), so `sqrt(self) =
    /// sqrt(m 2^32) 2^(E / 2 - 16)`: the integer root rounds once (`2^-31`
    /// relative), the conversion to 30 bits once more.
    fn sqrt_dir(&self, up: bool) -> Mag {
        if self.is_zero() || self.is_inf() {
            return self.clone();
        }
        let e = offset(&self.exp, -i64::from(BITS));
        let odd = match e.to_i64() {
            Some(v) => v & 1 == 1,
            None => e.unsigned_abs().limbs()[0] & 1 == 1,
        };
        let (m, e) = if odd {
            (u64::from(self.man) << 1, offset(&e, -1))
        } else {
            (u64::from(self.man), e)
        };
        let n = m << 32;
        let r = n.isqrt();
        let r = if up && r * r != n { r + 1 } else { r };
        let half = match e.to_i64() {
            Some(v) => Integer::from(v / 2),
            None => e.div_exact(&Integer::from(2)),
        };
        // r 2^(E/2 - 16) as m 2^(X - 30): X = E/2 + 14.
        Mag::from_u64(r, offset(&half, 14), up)
    }

    /// The `i64` exponent of a finite non-zero value, if it fits.
    fn exp_i64(&self) -> Option<i64> {
        self.exp.to_i64()
    }
}

/// Fixed-point helpers for the transcendental bounds: `u128` values with
/// `FRAC` fractional bits, every operation rounding up.
mod fixed {
    /// `(hi, lo)` of the 256-bit `a b`.
    pub(super) fn mul_full(a: u128, b: u128) -> (u128, u128) {
        let mask = u128::from(u64::MAX);
        let (a0, a1) = (a & mask, a >> 64);
        let (b0, b1) = (b & mask, b >> 64);
        let p00 = a0 * b0;
        let p01 = a0 * b1;
        let p10 = a1 * b0;
        let p11 = a1 * b1;
        let mid = (p00 >> 64) + (p01 & mask) + (p10 & mask);
        let lo = (p00 & mask) | (mid << 64);
        let hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
        (hi, lo)
    }

    /// `ceil(a b / 2^frac)` for a product below `2^(128 + frac)`.
    pub(super) fn mul_up(a: u128, b: u128, frac: u32) -> u128 {
        let (hi, lo) = mul_full(a, b);
        debug_assert!(frac > 0 && frac < 128 && hi >> frac == 0);
        let q = (hi << (128 - frac)) | (lo >> frac);
        let rest = lo & ((1u128 << frac) - 1);
        q + u128::from(rest != 0)
    }

    /// `ceil(a / d)`.
    pub(super) fn div_up(a: u128, d: u128) -> u128 {
        a.div_ceil(d)
    }
}

/// `m × 2^(e - 128)` with `m`'s top bit set, rounded up: the intermediate of
/// [`Mag::exp`]'s squarings, so that 70 squarings still leave the result well
/// within `2^-28`.
#[derive(Clone, Copy)]
struct Up128 {
    m: u128,
    e: i64,
}

impl Up128 {
    /// An upper bound of `self^2`; `None` if the exponent leaves `i64`.
    fn square(self) -> Option<Up128> {
        let (mut hi, mut lo) = fixed::mul_full(self.m, self.m);
        let mut e = self.e.checked_mul(2)?;
        if hi >> 127 == 0 {
            hi = (hi << 1) | (lo >> 127);
            lo <<= 1;
            e -= 1;
        }
        if lo != 0 {
            if let Some(h) = hi.checked_add(1) {
                hi = h;
            } else {
                hi = 1 << 127;
                e = e.checked_add(1)?;
            }
        }
        Some(Up128 { m: hi, e })
    }

    /// The value as a `Mag`, rounded up.
    fn to_mag(self) -> Mag {
        let top = self.m >> 98;
        let lost = self.m & ((1 << 98) - 1) != 0;
        let man = u64::try_from(top).expect("30 bits") + u64::from(lost);
        Mag::from_u64(man, Integer::from(self.e), true)
    }
}

/// Fractional bits of the `exp` Taylor sum: `1.0 = 2^126`, `e^y < 2`.
const EXP_FRAC: u32 = 126;
/// Terms of the `exp` series for `y < 2^-8`: the tail after them is below
/// `2 y^15 / 15! < 2^-160`, one unit of `2^-126` covers it.
const EXP_TERMS: u32 = 14;

/// `Σ_{k >= first} y^k / k!` for `y = Y 2^-126 < 2^-8`, rounded up, plus one
/// unit for the tail.
fn exp_series(y: u128, first: u32) -> u128 {
    let one = 1u128 << EXP_FRAC;
    let mut term = one;
    let mut sum = if first == 0 { one } else { 0 };
    for k in 1..=EXP_TERMS {
        term = fixed::div_up(fixed::mul_up(term, y, EXP_FRAC), u128::from(k));
        if k >= first {
            sum += term;
        }
    }
    sum + 1
}

impl Mag {
    /// `y = self 2^-s` in `EXP_FRAC` fixed point with `s = max(0, e + 8)`, for
    /// `2^-31 <= self < 2^62`: then `y < 2^-8`, and exact (its lowest bit is
    /// at `2^-60` or above).
    fn exp_reduced(e: i64, man: u32) -> (u128, u32) {
        let s = u32::try_from((e + 8).max(0)).expect("e < 62");
        // self = man 2^(e - 30); y 2^126 = man 2^(96 + e - s).
        let shift = u32::try_from(96 + e - i64::from(s)).expect("non-negative");
        (u128::from(man) << shift, s)
    }

    /// An upper bound of `e^self`.
    ///
    /// - Below `2^-31`: `e^x <= 1 + 2x <= 1 + 2^-29`, the next `Mag` above 1.
    /// - From `2^62` (and `+∞`): `+∞`.
    /// - Otherwise `y = x 2^-s < 2^-8` exactly, `e^y` by 14 Taylor terms in
    ///   126-bit fixed point, each rounded up, plus one unit for the tail; then
    ///   `s <= 70` squarings at 128 bits, rounded up. Each squaring doubles the
    ///   relative error and adds `2^-127`, so it stays below `2^70 × 2^-120`,
    ///   and the final rounding to 30 bits dominates.
    #[must_use]
    pub fn exp(&self) -> Mag {
        if self.is_zero() {
            return Mag::ONE;
        }
        if self.is_inf() {
            return Mag::INF;
        }
        let Some(e) = self.exp_i64() else {
            return if self.exp.is_negative() {
                Mag::ONE.next_up()
            } else {
                Mag::INF
            };
        };
        if e <= -31 {
            return Mag::ONE.next_up();
        }
        if e > 62 {
            return Mag::INF;
        }
        let (y, s) = Mag::exp_reduced(e, self.man);
        let sum = exp_series(y, 0);
        // sum in [2^126, 2^127): m = sum << 1 at e = 1.
        let mut r = Up128 { m: sum << 1, e: 1 };
        for _ in 0..s {
            match r.square() {
                Some(v) => r = v,
                None => return Mag::INF,
            }
        }
        r.to_mag()
    }

    /// An upper bound of `e^self - 1`, without cancellation near 0.
    ///
    /// - Below `2^-31`: `e^x - 1 <= x + x^2 <= x (1 + 2^-31)`, so the next
    ///   `Mag` above `x`.
    /// - Below `2^-8`: the series from `x^1` on in fixed point (relative
    ///   precision `2^-90` or better, as the sum is at least `2^-31`).
    /// - Otherwise [`Mag::exp`]'s 128-bit value minus 1 exactly: `e^x >= 1 +
    ///   2^-8`, so the subtraction magnifies the relative error 257 times at
    ///   most.
    #[must_use]
    pub fn expm1(&self) -> Mag {
        if self.is_zero() || self.is_inf() {
            return self.clone();
        }
        let Some(e) = self.exp_i64() else {
            return if self.exp.is_negative() {
                self.next_up()
            } else {
                Mag::INF
            };
        };
        if e <= -31 {
            return self.next_up();
        }
        if e > 62 {
            return Mag::INF;
        }
        let (y, s) = Mag::exp_reduced(e, self.man);
        if s == 0 {
            let sum = exp_series(y, 1);
            return fixed_to_mag(sum, EXP_FRAC);
        }
        let sum = exp_series(y, 0);
        let mut r = Up128 { m: sum << 1, e: 1 };
        for _ in 0..s {
            match r.square() {
                Some(v) => r = v,
                None => return Mag::INF,
            }
        }
        // r = m 2^(e - 128) >= 1 + 2^-8; subtract 1 = 2^(128 - e) units.
        if r.e > 128 {
            // 1 is below the last bit: e^x - 1 < e^x.
            return r.to_mag();
        }
        let one = 1u128 << (128 - r.e);
        let m = r.m - one;
        let lz = m.leading_zeros();
        Up128 {
            m: m << lz,
            e: r.e - i64::from(lz),
        }
        .to_mag()
    }

    /// An upper bound of `log(1 + self)`.
    ///
    /// - Below `2^-31`: `log(1 + x) <= x`.
    /// - Otherwise `1 + x = 2^k w` with `w` in `[1, 2]` (rounded up), and
    ///   `log(1 + x) = k log 2 + 2 atanh(t)`, `t = (w - 1) / (w + 1) <= 1/3`, in
    ///   64-bit fixed point: `t` and the series terms rounded up until a term
    ///   is a few units, then the rest bounded by `T_j (9/8) / (2j + 1)`;
    ///   `log 2` rounded up.
    ///   From `x >= 2^63` the 1 is dropped: `log(1 + x) <= log x + 1/x`, one
    ///   more unit.
    /// - Exponents beyond `i64`: `log(1 + x) <= exp(x) log 2 + 1` (in `Mag`
    ///   arithmetic, `exp(x)` here the exponent).
    #[must_use]
    pub fn log1p(&self) -> Mag {
        if self.is_zero() || self.is_inf() {
            return self.clone();
        }
        let Some(e) = self.exp_i64() else {
            if self.exp.is_negative() {
                return self.clone();
            }
            let ln2 = fixed_to_mag(u128::from(LN2_UP_64), 64);
            return Mag::from_integer(&self.exp).mul(&ln2).add(&Mag::ONE);
        };
        if e <= -31 {
            return self.clone();
        }
        // v = (1 + x) 2^64, or x 2^64 / 2^(e - 64) when the 1 is dropped.
        let (v, extra_k, tail) = if e < 64 {
            // -30 <= e < 64, so the shift is in 4..98.
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "4..98"
            )]
            let x = u128::from(self.man) << ((e + 34) as u32);
            ((1u128 << 64) + x, 0i64, 0)
        } else {
            // x / 2^(e - 64) = man 2^34, which is man 2^98 in 64-bit fixed;
            // log(1 + x) <= log x + 1/x, and 1/x < 2^-63 is 2 units.
            (u128::from(self.man) << 98, e - 64, 2)
        };
        // 2^k w with w in [1, 2): k = bits(v) - 65.
        let shift = v.ilog2() - 64;
        let k = i64::from(shift);
        let mut k = k;
        let mut w = (v >> shift) + u128::from(v & ((1u128 << shift) - 1) != 0);
        if w == 1 << 65 {
            // Rounded up to 2: w = 1 one binade higher.
            w = 1 << 64;
            k += 1;
        }
        let t = fixed::div_up((w - (1 << 64)) << 64, w + (1 << 64));
        let t2 = fixed::mul_up(t, t, 64);
        let mut term = t;
        let mut series = t;
        for j in 1u32.. {
            // term >= t^(2j + 1); once it is a few units, this term and all
            // later ones sum to at most T_j / (2j + 1) / (1 - t^2), and
            // 1 / (1 - t^2) <= 9/8 for t <= 1/3. Stopping there instead of
            // after a fixed count keeps the rounded-up tiny terms from adding
            // a unit each. t <= 1/3 makes it stop by j = 21.
            term = fixed::mul_up(term, t2, 64);
            let d = u128::from(2 * j + 1);
            if term <= 64 {
                series += fixed::div_up(term * 9, 8 * d);
                break;
            }
            series += fixed::div_up(term, d);
        }
        #[allow(clippy::cast_sign_loss, reason = "k and extra_k are non-negative")]
        let kk = (k + extra_k) as u128;
        let total = kk * u128::from(LN2_UP_64) + 2 * series + tail;
        fixed_to_mag(total, 64)
    }
}

/// `log 2 × 2^64`, rounded up.
const LN2_UP_64: u64 = 0xB172_17F7_D1CF_79AC;

/// `v 2^-frac` as a `Mag`, rounded up (`v > 0`).
fn fixed_to_mag(v: u128, frac: u32) -> Mag {
    let bits = 128 - v.leading_zeros();
    // from_u64(m, E) is m 2^(E - 30).
    let base = i64::from(BITS) - i64::from(frac);
    if bits <= 62 {
        return Mag::from_u64(
            u64::try_from(v).expect("62 bits"),
            Integer::from(base),
            true,
        );
    }
    let shift = bits - 62;
    let top = u64::try_from(v >> shift).expect("62 bits");
    let lost = v & ((1u128 << shift) - 1) != 0;
    Mag::from_u64(
        top + u64::from(lost),
        Integer::from(base + i64::from(shift)),
        true,
    )
}

impl Ord for Mag {
    fn cmp(&self, other: &Mag) -> Ordering {
        match (
            self.is_zero(),
            other.is_zero(),
            self.is_inf(),
            other.is_inf(),
        ) {
            (true, true, ..) | (_, _, true, true) => Ordering::Equal,
            (true, ..) | (_, _, _, true) => Ordering::Less,
            (_, true, ..) | (_, _, true, _) => Ordering::Greater,
            _ => self.exp.cmp(&other.exp).then(self.man.cmp(&other.man)),
        }
    }
}

impl PartialOrd for Mag {
    fn partial_cmp(&self, other: &Mag) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Default for Mag {
    fn default() -> Self {
        Mag::ZERO
    }
}

/// `0x2abcdef0p-30`-style: the mantissa in hex times a power of two, or `0`,
/// `inf`.
impl fmt::Debug for Mag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            f.write_str("0")
        } else if self.is_inf() {
            f.write_str("inf")
        } else {
            write!(f, "{:#x}p{}", self.man, &self.exp - &Integer::from(BITS))
        }
    }
}
