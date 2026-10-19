//! Addition, subtraction and multiplication, exact and correctly rounded.
//!
//! The general path aligns both operands in a buffer (on the stack up to 16
//! limbs), adds or subtracts, and normalises and rounds in place. One-limb
//! operands with `i64` exponents and precisions up to 64 bits take a `u128`
//! path, rounded by a straight-line copy of the shared routine.
#![allow(
    clippy::many_single_char_names,
    reason = "x, y, e, m, n, p, r: the operands and the textbook names"
)]

use std::cmp::Ordering;

use super::round::round_1;
use super::{EXACT_BITS_MAX, Float, LIMB_BITS, Round, with_limbs};
use crate::kernel;
use crate::{Integer, Limb};

/// The top limb of `0.1 × 2^e`, the mantissa of a power of two.
const HALF: Limb = 1 << (Limb::BITS - 1);

impl Float {
    /// `self + rhs`, exactly.
    ///
    /// # Panics
    ///
    /// If the exact sum would have more than [`EXACT_BITS_MAX`] bits (the
    /// span from the higher operand's top bit to the lower one's lowest);
    /// [`Float::checked_add`] returns `None` instead.
    #[must_use]
    pub fn add(&self, rhs: &Float) -> Float {
        self.checked_add(rhs)
            .expect("Float::add: exact sum beyond EXACT_BITS_MAX bits")
    }

    /// `self - rhs`, exactly.
    ///
    /// # Panics
    ///
    /// As [`Float::add`]; [`Float::checked_sub`] returns `None` instead.
    #[must_use]
    pub fn sub(&self, rhs: &Float) -> Float {
        self.checked_sub(rhs)
            .expect("Float::sub: exact difference beyond EXACT_BITS_MAX bits")
    }

    /// `self + rhs`, exactly, or `None` if that would take more than
    /// [`EXACT_BITS_MAX`] bits.
    #[must_use]
    pub fn checked_add(&self, rhs: &Float) -> Option<Float> {
        exact_sum_fits(self, rhs).then(|| sum(self, rhs, rhs.neg, None).0)
    }

    /// `self - rhs`, exactly, or `None` if that would take more than
    /// [`EXACT_BITS_MAX`] bits.
    #[must_use]
    pub fn checked_sub(&self, rhs: &Float) -> Option<Float> {
        exact_sum_fits(self, rhs).then(|| sum(self, rhs, !rhs.neg && !rhs.is_zero(), None).0)
    }

    /// `self + rhs` correctly rounded to `prec` bits, and `rounded.cmp(exact)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    #[inline]
    pub fn add_round(&self, rhs: &Float, prec: u64, rnd: Round) -> (Float, Ordering) {
        assert!(prec >= 1, "Float::add_round: precision 0");
        if let Some(r) = sum_fast(self, rhs, rhs.neg, prec, rnd) {
            return r;
        }
        sum(self, rhs, rhs.neg, Some((prec, rnd)))
    }

    /// `self - rhs` correctly rounded to `prec` bits, and `rounded.cmp(exact)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    #[inline]
    pub fn sub_round(&self, rhs: &Float, prec: u64, rnd: Round) -> (Float, Ordering) {
        assert!(prec >= 1, "Float::sub_round: precision 0");
        if let Some(r) = sum_fast(self, rhs, !rhs.neg, prec, rnd) {
            return r;
        }
        sum(self, rhs, !rhs.neg && !rhs.is_zero(), Some((prec, rnd)))
    }

    /// `self × rhs`, exactly.
    #[must_use]
    pub fn mul(&self, rhs: &Float) -> Float {
        product(self, rhs, None).0
    }

    /// `self × rhs` correctly rounded to `prec` bits, and `rounded.cmp(exact)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    #[inline]
    pub fn mul_round(&self, rhs: &Float, prec: u64, rnd: Round) -> (Float, Ordering) {
        assert!(prec >= 1, "Float::mul_round: precision 0");
        if let ([x], [y]) = (self.mant(), rhs.mant())
            && prec <= LIMB_BITS
            && let Some(r) = product_1(self.neg != rhs.neg, *x, *y, &self.exp, &rhs.exp, prec, rnd)
        {
            return r;
        }
        product(self, rhs, Some((prec, rnd)))
    }
}

/// `a × b`, exact or rounded.
fn product(a: &Float, b: &Float, rnd: Option<(u64, Round)>) -> (Float, Ordering) {
    if a.is_zero() || b.is_zero() {
        return (Float::ZERO, Ordering::Equal);
    }
    let neg = a.neg != b.neg;
    let (x, y) = (a.mant(), b.mant());
    let n = x.len() + y.len();
    with_limbs(n + kernel::mul_scratch_len(x.len(), y.len()), |buf| {
        let (r, scratch) = buf.split_at_mut(n);
        kernel::mul(r, x, y, scratch);
        // 0.x 0.y = (x y) / B^n.
        let bits = LIMB_BITS * n as u64;
        let e = a
            .exp
            .to_i64()
            .zip(b.exp.to_i64())
            .and_then(|(ea, eb)| ea.checked_add(eb)?.checked_sub(i64::try_from(bits).ok()?))
            .map_or_else(|| &a.exp + &b.exp - Integer::from(bits), Integer::from);
        Float::build(neg, r, &e, false, rnd)
    })
}

/// [`product`] of one-limb mantissas through a `u128`; `None` when an
/// exponent leaves `i64`.
#[inline]
fn product_1(
    neg: bool,
    x: Limb,
    y: Limb,
    ex: &Integer,
    ey: &Integer,
    prec: u64,
    rnd: Round,
) -> Option<(Float, Ordering)> {
    let e = ex.to_i64()?.checked_add(ey.to_i64()?)?;
    let p = u128::from(x) * u128::from(y);
    // Both top bits are set, so the product has 127 or 128 bits.
    let lz = p.leading_zeros();
    let p = p << lz;
    finish_2(neg, p, false, e.checked_sub(i64::from(lz))?, prec, rnd)
}

/// Rounds the top-aligned `±0.p × 2^e` (with set bits below `p` if `sticky`)
/// to `prec <= 64` bits.
#[inline]
fn finish_2(
    neg: bool,
    p: u128,
    sticky: bool,
    e: i64,
    prec: u64,
    rnd: Round,
) -> Option<(Float, Ordering)> {
    #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    let (hi, lo) = ((p >> Limb::BITS) as Limb, p as Limb);
    let (m, carry, ord) = round_1(neg, hi, lo, sticky, prec, rnd);
    let e = if carry { e.checked_add(1)? } else { e };
    Some((Float::from_limb(neg, m, e), ord))
}

/// `a + (±|b|)`, where `b_neg` is the sign given to `b`, exactly or rounded.
///
/// With rounding, an operand far below the other is replaced by a single bit
/// that rounds the same way, so the work is bounded by the precision and the
/// larger operand, never by the exponent gap (`1 + 2^-10^6` at 64 bits). Let
/// `x` be the operand with the larger exponent and `L = exp(x) - max(64 kx,
/// prec) - 2`. If `|y| < 2^L`, the exact `x ± y` and `x ± 2^(L - 1)` agree above
/// `2^L` (for a difference, both are `x - 2^L` plus something in `(0, 2^L)`),
/// and both have set bits below it. The result's exponent is at least
/// `exp(x) - 1`, so its last kept bit is at `2^(L + 1)` or above, and no
/// representable value lies strictly between the two: they round alike, with
/// the same ternary value.
fn sum(a: &Float, b: &Float, b_neg: bool, rnd: Option<(u64, Round)>) -> (Float, Ordering) {
    let finish = |x: Float| match rnd {
        Some((prec, mode)) => x.round(prec, mode),
        None => (x, Ordering::Equal),
    };
    if b.is_zero() {
        return finish(a.clone());
    }
    if a.is_zero() {
        return finish(Float {
            neg: b_neg,
            ..b.clone()
        });
    }
    // x has the larger exponent.
    let (x, x_neg, y, y_neg) = if a.exp >= b.exp {
        (a, a.neg, b, b_neg)
    } else {
        (b, b_neg, a, a.neg)
    };
    let prec = rnd.map(|(prec, _)| prec);
    let Layout {
        y_mant,
        ox,
        oy,
        width,
        low,
    } = Layout::small(x, y, prec).unwrap_or_else(|| Layout::big(x, y, prec));
    let n = usize::try_from(width.div_ceil(LIMB_BITS)).expect("Float: sum too large");
    with_limbs(2 * n, |buf| {
        let (bx, by) = buf.split_at_mut(n);
        place(bx, x.mant(), ox);
        place(by, y_mant, oy);
        let (r, neg): (&mut [Limb], bool) = if x_neg == y_neg {
            let carry = kernel::add_n_assign(bx, by);
            debug_assert_eq!(carry, 0);
            (bx, x_neg)
        } else {
            match kernel::cmp_n(bx, by) {
                Ordering::Equal => return (Float::ZERO, Ordering::Equal),
                Ordering::Greater => {
                    kernel::sub_n_assign(bx, by);
                    (bx, x_neg)
                }
                Ordering::Less => {
                    kernel::sub_n_assign(by, bx);
                    (by, y_neg)
                }
            }
        };
        Float::build(neg, r, &low, false, rnd)
    })
}

/// Where the operands of [`sum`] go in its buffer: `y`'s mantissa (or the
/// single bit standing for it when it is far below), the bit offsets of both
/// mantissas above the lowest bit `low`, and the buffer's width in bits.
struct Layout<'a> {
    y_mant: &'a [Limb],
    ox: u64,
    oy: u64,
    width: u64,
    low: Integer,
}

impl<'a> Layout<'a> {
    /// The layout in `i64` arithmetic; `None` if anything leaves `i64`.
    fn small(x: &Float, y: &'a Float, prec: Option<u64>) -> Option<Layout<'a>> {
        let (xe, mut ye) = (x.exp.to_i64()?, y.exp.to_i64()?);
        let kx = i64::try_from(x.mant_bits_stored()).ok()?;
        let x_low = xe.checked_sub(kx)?;
        let mut y_mant = y.mant();
        if let Some(prec) = prec {
            let l = xe.checked_sub(kx.max(i64::try_from(prec).ok()?).checked_add(2)?)?;
            if ye <= l {
                (y_mant, ye) = (&[HALF], l);
            }
        }
        let y_low = ye.checked_sub(i64::try_from(LIMB_BITS * y_mant.len() as u64).ok()?)?;
        let low = x_low.min(y_low);
        let off = |v: i64| v.checked_sub(low).and_then(|d| u64::try_from(d).ok());
        Some(Layout {
            y_mant,
            ox: off(x_low)?,
            oy: off(y_low)?,
            width: off(xe)?.checked_add(1)?,
            low: Integer::from(low),
        })
    }

    /// The layout in `Integer` arithmetic, for exponents beyond `i64`.
    fn big(x: &Float, y: &'a Float, prec: Option<u64>) -> Layout<'a> {
        let x_low = &x.exp - &Integer::from(x.mant_bits_stored());
        let far_below = prec.and_then(|prec| {
            let l = &x.exp - &Integer::from(x.mant_bits_stored().max(prec) + 2);
            (y.exp <= l).then_some(l)
        });
        let (y_mant, y_exp): (&[Limb], Integer) = match far_below {
            Some(l) => (&[HALF], l),
            None => (y.mant(), y.exp.clone()),
        };
        let y_low = &y_exp - &Integer::from(LIMB_BITS * y_mant.len() as u64);
        let low = if x_low <= y_low {
            x_low.clone()
        } else {
            y_low.clone()
        };
        let off = |v: &Integer| -> u64 {
            (v - &low)
                .to_i64()
                .and_then(|d| u64::try_from(d).ok())
                .expect("Float: exponent gap too large for an exact sum")
        };
        Layout {
            y_mant,
            ox: off(&x_low),
            oy: off(&y_low),
            // One bit above x's top for the carry.
            width: off(&x.exp) + 1,
            low,
        }
    }
}

/// [`sum`] rounded, when both mantissas are one limb (so neither operand is
/// zero) and `prec <= 64`; `None` otherwise or when an exponent leaves `i64`.
#[inline]
fn sum_fast(a: &Float, b: &Float, b_neg: bool, prec: u64, rnd: Round) -> Option<(Float, Ordering)> {
    match (a.mant(), b.mant()) {
        ([x], [y]) if prec <= LIMB_BITS => sum_1(a.neg, *x, &a.exp, b_neg, *y, &b.exp, prec, rnd),
        _ => None,
    }
}

/// [`sum`] of one-limb mantissas through a `u128` window; `None` when an
/// exponent leaves `i64`.
///
/// `x`, the operand with the larger exponent, sits at bits 63..127, with bit
/// 127 free for the carry; `y` is shifted right by the exponent gap `d`.
/// Below a gap of 64 nothing is lost (the window has 63 zero bits under `y`).
/// From 64 on, the bits that fall off only matter through a sticky bit: the
/// true value is the integer part of the window plus a fraction in `(0, 1)`
/// (for a difference, `x - y - 1` plus `1 - frac`). Then `y < 2^63` and the
/// result keeps at least 125 bits, so normalising shifts by at most 2 and the
/// rounding bit stays far above the unknown low bits.
#[allow(clippy::too_many_arguments, reason = "two unpacked operands")]
#[inline]
fn sum_1(
    a_neg: bool,
    a: Limb,
    ea: &Integer,
    b_neg: bool,
    b: Limb,
    eb: &Integer,
    prec: u64,
    rnd: Round,
) -> Option<(Float, Ordering)> {
    let (ea, eb) = (ea.to_i64()?, eb.to_i64()?);
    let ((ex, mx, nx), (ey, my, ny)) = if ea >= eb {
        ((ea, a, a_neg), (eb, b, b_neg))
    } else {
        ((eb, b, b_neg), (ea, a, a_neg))
    };
    let d = ex.abs_diff(ey);
    let x = u128::from(mx) << 63;
    let full = u128::from(my) << 63;
    let (y, lost) = if d < 127 {
        let y = full >> d;
        (y, y << d != full)
    } else {
        (0, true)
    };
    let (s, neg) = if nx == ny {
        (x + y, nx)
    } else {
        match x.cmp(&y) {
            Ordering::Greater => (x - y - u128::from(lost), nx),
            // Only at d = 0, where nothing is lost.
            Ordering::Less => (y - x, ny),
            Ordering::Equal => return Some((Float::ZERO, Ordering::Equal)),
        }
    };
    let lz = s.leading_zeros();
    // The window's bit 127 has weight 2^(ex - 1): x's top bit is at 126.
    let e = ex.checked_add(1)?.checked_sub(i64::from(lz))?;
    finish_2(neg, s << lz, lost, e, prec, rnd)
}

/// Writes `mant << offset` bits into the zeroed `v`.
fn place(v: &mut [Limb], mant: &[Limb], offset: u64) {
    let q = usize::try_from(offset / LIMB_BITS).expect("inside the buffer");
    #[allow(clippy::cast_possible_truncation, reason = "below 64")]
    let r = (offset % LIMB_BITS) as u32;
    let k = mant.len();
    let out = kernel::shl(&mut v[q..q + k], mant, r);
    if q + k < v.len() {
        v[q + k] = out;
    } else {
        debug_assert_eq!(out, 0);
    }
}

/// Whether the exact sum of `a` and `b` fits [`EXACT_BITS_MAX`] bits: from
/// the higher top bit to the lower lowest bit, plus a carry.
fn exact_sum_fits(a: &Float, b: &Float) -> bool {
    if a.is_zero() || b.is_zero() {
        return true;
    }
    let low = |x: &Float| x.exponent() - &Integer::from(x.bit_length());
    let top = std::cmp::max(a.exponent(), b.exponent());
    let bottom = std::cmp::min(low(a), low(b));
    (top - &bottom) < Integer::from(EXACT_BITS_MAX)
}
