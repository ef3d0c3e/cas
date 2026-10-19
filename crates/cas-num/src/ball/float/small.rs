//! Fixed-size paths for rounded sums and products of up to four limbs
//! (CAS-91): operands of at most `N` limbs, a precision of at most `64 N`
//! bits and `i64` exponents, on fixed arrays with `N` a constant, so
//! that the loops unroll and nothing is zeroed or allocated beyond the
//! window; rounding and the result's inline mantissa are fixed-size too. The
//! general path in [`super::arith`] stays the reference; the tests compare
//! both with MPFR.
//!
//! **Sum.** As the one-limb path (`sum_1`), generalised: `x`, the operand
//! with the larger exponent, sits top-aligned in a window of `W = N + 1`
//! limbs, and `y` likewise, shifted right by the exponent gap `d`. `y`'s
//! lowest bit stays inside the window while `d <= 64 (W - ky)`, so bits are
//! lost only from `d > 64` on. Then `y < 2^(64W - 65)` while `x >= 2^(64W -
//! 1)`: the true value is the window's integer value plus a fraction in `(0,
//! 1)` (for a difference, `x - y - 1` plus `1 - frac`), normalising shifts by
//! at most one bit either way, and the rounding bit, at `64 W - 1 - 64 N - 1
//! = 62` or above, lies above the unknown fraction and the bit a left shift
//! brings in: the lost bits only count through the sticky bit. A carry out
//! of the window shifts it right by one bit, whose lowest bit joins the
//! sticky bit for the same reason.
//!
//! **Product.** The schoolbook product of the `N`-limb mantissas (shorter
//! ones padded with low zero limbs, which leave the value alone) is in
//! `[B^W / 4, B^W)`, normalised by at most one bit.
#![allow(
    clippy::many_single_char_names,
    reason = "a, b, x, y, w, e, z, p: the operands, window and textbook names"
)]
#![allow(
    clippy::inline_always,
    reason = "measured: without it the 2-limb shapes keep `finish` out of line, 1.3-1.7x slower"
)]

use std::cmp::Ordering;

use super::{Float, LIMB_BITS, Round};
use crate::natural::limb_vec::{INLINE, LimbVec};
use crate::{Integer, Limb};

/// The top bit of a limb.
const TOP: Limb = 1 << (Limb::BITS - 1);

/// A rounded result: sign, inline mantissa, `i64` exponent and ternary
/// value, built into a new value or written into a destination.
pub(super) struct Parts {
    neg: bool,
    mant: LimbVec,
    exp: i64,
    ord: Ordering,
}

impl Parts {
    const ZERO: Parts = Parts {
        neg: false,
        mant: LimbVec::new(),
        exp: 0,
        ord: Ordering::Equal,
    };

    /// The result as a new value.
    #[inline]
    pub(super) fn build(self) -> (Float, Ordering) {
        (
            Float::from_parts(self.neg, self.mant, Integer::from(self.exp)),
            self.ord,
        )
    }

    /// The result written into `out`, field by field (the old exponent is
    /// only freed if it was boxed); its ternary value.
    #[inline]
    pub(super) fn write(self, out: &mut Float) -> Ordering {
        out.mant = self.mant;
        out.exp = Integer::from(self.exp);
        out.neg = self.neg;
        out.debug_check();
        self.ord
    }
}

/// Limbs of `prec` bits.
fn prec_limbs(prec: u64) -> Option<usize> {
    usize::try_from(prec.div_ceil(LIMB_BITS)).ok()
}

/// `a + (±|b|)` rounded, `b_neg` the sign given to `b`, through a fixed
/// window when both operands and the precision take 2 to 4 limbs; `None`
/// otherwise or when an exponent leaves `i64`.
#[inline]
pub(super) fn sum(a: &Float, b: &Float, b_neg: bool, prec: u64, rnd: Round) -> Option<Parts> {
    let (x, y) = (a.mant(), b.mant());
    if x.is_empty() || y.is_empty() {
        return None;
    }
    let (ea, eb) = (a.exp.to_i64()?, b.exp.to_i64()?);
    let (x, y) = if ea >= eb {
        ((ea, x, a.neg), (eb, y, b_neg))
    } else {
        ((eb, y, b_neg), (ea, x, a.neg))
    };
    match x.1.len().max(y.1.len()).max(prec_limbs(prec)?) {
        2 => sum_n::<2, 3>(x, y, prec, rnd),
        3 => sum_n::<3, 4>(x, y, prec, rnd),
        4 => sum_n::<4, 5>(x, y, prec, rnd),
        _ => None,
    }
}

/// `a × b` rounded through a fixed product when both operands and the
/// precision take 2 to 4 limbs; `None` otherwise or when an exponent leaves
/// `i64`.
#[inline]
pub(super) fn product(a: &Float, b: &Float, prec: u64, rnd: Round) -> Option<Parts> {
    let (x, y) = (a.mant(), b.mant());
    if x.is_empty() || y.is_empty() {
        return None;
    }
    let e = a.exp.to_i64()?.checked_add(b.exp.to_i64()?)?;
    let neg = a.neg != b.neg;
    match x.len().max(y.len()).max(prec_limbs(prec)?) {
        2 => product_n::<2, 4>(x, y, e, neg, prec, rnd),
        3 => product_n::<3, 6>(x, y, e, neg, prec, rnd),
        4 => product_n::<4, 8>(x, y, e, neg, prec, rnd),
        _ => None,
    }
}

/// `m` top-aligned in `W` limbs, `m.len() <= W`.
#[inline(always)]
fn top_aligned<const W: usize>(m: &[Limb]) -> [Limb; W] {
    let mut w = [0; W];
    let off = W - m.len();
    for (i, &x) in m.iter().enumerate() {
        w[off + i] = x;
    }
    w
}

/// `w >> k` bits in place; whether set bits fell off. Plain loops over the
/// fixed array: `fill` and `copy_within` of a variable length become
/// `memset` and `memmove` calls, slower here than the few limbs they move.
#[inline(always)]
fn shr_window<const W: usize>(w: &mut [Limb; W], k: u64) -> bool {
    let q = usize::try_from(k / LIMB_BITS).unwrap_or(usize::MAX);
    if q >= W {
        let mut lost = false;
        for x in w.iter_mut() {
            lost |= *x != 0;
            *x = 0;
        }
        return lost;
    }
    #[allow(clippy::cast_possible_truncation, reason = "below 64")]
    let r = (k % LIMB_BITS) as u32;
    let mut lost = false;
    for &x in &w[..q] {
        lost |= x != 0;
    }
    if r != 0 {
        lost |= w[q] << (Limb::BITS - r) != 0;
    }
    for i in 0..W {
        let lo = if i + q < W { w[i + q] } else { 0 };
        let hi = if i + q + 1 < W { w[i + q + 1] } else { 0 };
        w[i] = if r == 0 {
            lo
        } else {
            (lo >> r) | (hi << (Limb::BITS - r))
        };
    }
    lost
}

/// `w << k` bits in place, for `k < 64 W` and no set bit pushed out.
#[inline(always)]
fn shl_window<const W: usize>(w: &mut [Limb; W], k: u32) {
    if k == 0 {
        return;
    }
    let q = usize::try_from(k / Limb::BITS).expect("small");
    let r = k % Limb::BITS;
    for i in (0..W).rev() {
        let hi = if i >= q { w[i - q] } else { 0 };
        let lo = if i > q { w[i - q - 1] } else { 0 };
        w[i] = if r == 0 {
            hi
        } else {
            (hi << r) | (lo >> (Limb::BITS - r))
        };
    }
}

/// Leading zero bits of the `64 W`-bit `w`.
#[inline(always)]
fn leading_zeros<const W: usize>(w: &[Limb; W]) -> u32 {
    let mut lz = 0;
    for &x in w.iter().rev() {
        lz += x.leading_zeros();
        if x != 0 {
            break;
        }
    }
    lz
}

/// `x += y`, returning the carry.
#[inline(always)]
fn add_assign<const W: usize>(x: &mut [Limb; W], y: &[Limb; W]) -> bool {
    let mut carry = false;
    for i in 0..W {
        let (s, c1) = x[i].overflowing_add(y[i]);
        let (s, c2) = s.overflowing_add(Limb::from(carry));
        x[i] = s;
        carry = c1 || c2;
    }
    carry
}

/// `x -= y + borrow_in`, returning the borrow.
#[inline(always)]
fn sub_assign<const W: usize>(x: &mut [Limb; W], y: &[Limb; W], borrow_in: bool) -> bool {
    let mut borrow = borrow_in;
    for i in 0..W {
        let (s, b1) = x[i].overflowing_sub(y[i]);
        let (s, b2) = s.overflowing_sub(Limb::from(borrow));
        x[i] = s;
        borrow = b1 || b2;
    }
    borrow
}

/// `x.cmp(y)` as `64 W`-bit numbers.
#[inline(always)]
fn cmp<const W: usize>(x: &[Limb; W], y: &[Limb; W]) -> Ordering {
    for i in (0..W).rev() {
        match x[i].cmp(&y[i]) {
            Ordering::Equal => {}
            other => return other,
        }
    }
    Ordering::Equal
}

/// Rounds the top-aligned window `w` (top bit set; set bits below it if
/// `sticky`) of `±0.w × 2^e` to `prec` bits, `prec <= 64 (W - 1)`, and builds
/// the result: `round_mant` on a fixed window, its mantissa written inline.
#[inline(always)]
fn finish<const W: usize>(
    neg: bool,
    w: &mut [Limb; W],
    sticky: bool,
    e: i64,
    prec: u64,
    rnd: Round,
) -> Option<Parts> {
    debug_assert!(w[W - 1] & TOP != 0);
    let kept = prec_limbs(prec)?;
    debug_assert!(kept < W && kept <= INLINE);
    let low = W - kept;
    #[allow(clippy::cast_possible_truncation, reason = "below 64")]
    let z = (LIMB_BITS * kept as u64 - prec) as u32;
    let (round_bit, mut rest) = if z > 0 {
        let below = w[low] & ((1 << z) - 1);
        let half = 1 << (z - 1);
        (below & half != 0, below & (half - 1) != 0)
    } else {
        (w[low - 1] & TOP != 0, w[low - 1] << 1 != 0)
    };
    let below_end = if z > 0 { low } else { low - 1 };
    rest |= sticky || w[..below_end].iter().any(|&x| x != 0);
    let last_bit = (w[low] >> z) & 1 == 1;
    w[low] &= !((1 << z) - 1);
    let mut e = e;
    let ord = if round_bit || rest {
        let up = rnd.increments(neg, round_bit, rest, last_bit);
        if up {
            let mut carry = true;
            let mut add = 1 << z;
            for x in &mut w[low..] {
                let (s, c) = x.overflowing_add(add);
                *x = s;
                carry = c;
                if !c {
                    break;
                }
                add = 1;
            }
            if carry {
                w[W - 1] = TOP;
                e = e.checked_add(1)?;
            }
        }
        let mag = if up {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if neg { mag.reverse() } else { mag }
    } else {
        Ordering::Equal
    };
    // The kept limbs without their low zero limbs, written inline.
    let mut first = low;
    while w[first] == 0 {
        first += 1;
    }
    let mut buf = [0; INLINE];
    let len = W - first;
    // An indexed loop: `copy_from_slice`, or a zip, of a variable length
    // becomes a `memcpy` call, measured 1.6x slower here.
    #[allow(clippy::manual_memcpy, reason = "see above")]
    for i in 0..len {
        buf[i] = w[first + i];
    }
    let mant = LimbVec::Inline {
        len: u8::try_from(len).expect("at most INLINE"),
        buf,
    };
    Some(Parts {
        neg,
        mant,
        exp: e,
        ord,
    })
}

/// `sum` on a window of `W = N + 1` limbs; `x` has the larger exponent.
#[inline(always)]
fn sum_n<const N: usize, const W: usize>(
    (ex, x, x_neg): (i64, &[Limb], bool),
    (ey, y, y_neg): (i64, &[Limb], bool),
    prec: u64,
    rnd: Round,
) -> Option<Parts> {
    debug_assert_eq!(W, N + 1);
    let d = ex.abs_diff(ey);
    let mut wx = top_aligned::<W>(x);
    let mut wy = top_aligned::<W>(y);
    let lost = d > 0 && shr_window(&mut wy, d);
    debug_assert!(!lost || d > LIMB_BITS);
    // The window's top bit has weight 2^(ex - 1).
    let (neg, sticky, e) = if x_neg == y_neg {
        if add_assign(&mut wx, &wy) {
            // A carry: shift it in from the top; the bit out joins the sticky
            // bit, far below the rounding bit.
            let out = wx[0] & 1 != 0;
            shr_window(&mut wx, 1);
            wx[W - 1] |= TOP;
            (x_neg, lost || out, ex.checked_add(1)?)
        } else {
            (x_neg, lost, ex)
        }
    } else {
        let neg = match cmp(&wx, &wy) {
            Ordering::Greater => {
                sub_assign(&mut wx, &wy, lost);
                x_neg
            }
            // Only for small gaps, where nothing is lost.
            Ordering::Less => {
                debug_assert!(!lost);
                sub_assign(&mut wy, &wx, false);
                wx = wy;
                y_neg
            }
            Ordering::Equal => {
                debug_assert!(!lost);
                return Some(Parts::ZERO);
            }
        };
        let lz = leading_zeros(&wx);
        debug_assert!(!lost || lz <= 1);
        shl_window(&mut wx, lz);
        (neg, lost, ex.checked_sub(i64::from(lz))?)
    };
    finish(neg, &mut wx, sticky, e, prec, rnd)
}

/// `product` of `N`-limb mantissas into `W = 2N` limbs, at `2^e`.
#[inline(always)]
fn product_n<const N: usize, const W: usize>(
    x: &[Limb],
    y: &[Limb],
    e: i64,
    neg: bool,
    prec: u64,
    rnd: Round,
) -> Option<Parts> {
    debug_assert_eq!(W, 2 * N);
    let x = top_aligned::<N>(x);
    let y = top_aligned::<N>(y);
    let mut p = [0; W];
    for i in 0..N {
        let mut carry: Limb = 0;
        for j in 0..N {
            let t = u128::from(x[i]) * u128::from(y[j]) + u128::from(p[i + j]) + u128::from(carry);
            #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
            {
                p[i + j] = t as Limb;
                carry = (t >> Limb::BITS) as Limb;
            }
        }
        p[i + N] = carry;
    }
    // 0.x 0.y in [1/4, 1): at most one bit of normalisation.
    let e = if p[W - 1] & TOP == 0 {
        shl_window(&mut p, 1);
        e.checked_sub(1)?
    } else {
        e
    };
    finish(neg, &mut p, false, e, prec, rnd)
}
