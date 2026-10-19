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
    clippy::manual_memcpy,
    reason = "indexed copies: `copy_from_slice` of a variable length becomes a `memcpy` call, measured 1.6x slower here"
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

/// `a / b` rounded through a fixed schoolbook division when both operands
/// and the precision take 2 to 4 limbs (CAS-99); `None` otherwise, for a
/// zero operand, or when an exponent leaves `i64`.
#[inline]
pub(super) fn quotient(a: &Float, b: &Float, prec: u64, rnd: Round) -> Option<Parts> {
    let (x, y) = (a.mant(), b.mant());
    if x.is_empty() || y.is_empty() {
        return None;
    }
    let e = a.exp.to_i64()?.checked_sub(b.exp.to_i64()?)?;
    let neg = a.neg != b.neg;
    match x.len().max(y.len()).max(prec_limbs(prec)?).max(2) {
        2 => quotient_n::<2, 4, 5>(x, y, e, neg, prec, rnd),
        3 => quotient_n::<3, 5, 7>(x, y, e, neg, prec, rnd),
        4 => quotient_n::<4, 6, 9>(x, y, e, neg, prec, rnd),
        _ => None,
    }
}

/// `sqrt(a)` rounded through integer Newton steps on fixed arrays when the
/// operand and the precision take up to 4 limbs (`prec + 1` bits then 2 to
/// 5) (CAS-99);
/// `None` otherwise, for zero or a negative operand, or when the exponent
/// leaves `i64`.
#[inline]
pub(super) fn root(a: &Float, prec: u64, rnd: Round) -> Option<Parts> {
    let x = a.mant();
    // The result's mantissa is written inline: at most 4 limbs.
    if x.is_empty() || a.neg || prec_limbs(prec)? > INLINE {
        return None;
    }
    let e = a.exp.to_i64()?;
    match x.len().max(prec_limbs(prec.checked_add(1)?)?).max(2) {
        2 => root_n::<2, 5, 3, 6, 6>(x, e, prec, rnd),
        3 => root_n::<3, 7, 4, 8, 8>(x, e, prec, rnd),
        4 => root_n::<4, 9, 5, 10, 10>(x, e, prec, rnd),
        5 => root_n::<5, 11, 6, 12, 12>(x, e, prec, rnd),
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
    assert!(
        kept < W && kept <= INLINE,
        "small::finish: precision beyond the window"
    );
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
    // An indexed loop, not `copy_from_slice` (see the module's allows).
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

/// `quotient` of `N`-limb mantissas (shorter ones padded with low zero
/// limbs), at `2^e`: `A B^(N+1) / D` by schoolbook division
/// (`kernel::div::sb_div_qr`, Knuth's Algorithm D \[Knu97\] with the 3-by-2
/// reciprocal of \[MG11\]) on `NN = 2N + 1` limbs, the remainder as the
/// sticky bit. As `A / D` lies in `(1/2, 2)`, the quotient `Q` lies in
/// `(B^(N+1) / 2, 2 B^(N+1))`: at least `64 (N + 1) > prec` bits. Placed in
/// a window of `W = N + 2` limbs (its top bit in limb `N + 1` or at the top
/// of limb `N`) and normalised by `lz` bits,
/// `Q / B^(N+1) = 0.w × 2^(64 - lz)`.
#[inline(always)]
fn quotient_n<const N: usize, const W: usize, const NN: usize>(
    x: &[Limb],
    y: &[Limb],
    e: i64,
    neg: bool,
    prec: u64,
    rnd: Round,
) -> Option<Parts> {
    debug_assert!(W == N + 2 && NN == 2 * N + 1);
    let mut num = top_aligned::<NN>(x);
    let d = top_aligned::<N>(y);
    let v = crate::kernel::div::reciprocal_3by2(d[N - 1], d[N - 2]);
    let mut q = [0; W];
    let qh = crate::kernel::div::sb_div_qr(&mut q[..=N], &mut num, &d, v);
    q[N + 1] = qh;
    let mut sticky = false;
    for &r in &num[..N] {
        sticky |= r != 0;
    }
    let lz = leading_zeros(&q);
    debug_assert!(lz == 63 || lz == 64);
    shl_window(&mut q, lz);
    finish(
        neg,
        &mut q,
        sticky,
        e.checked_add(64 - i64::from(lz))?,
        prec,
        rnd,
    )
}

/// `floor(x / d)` into `q`, for `x` of `NX` limbs and `d` of at least two
/// significant limbs: both shifted so that `d`'s top bit is set
/// (`NX1 = NX + 1` limbs for `x`), then [`sb_div_qr`](crate::kernel::div::sb_div_qr).
#[inline(always)]
fn div_floor<const NX: usize, const NX1: usize>(x: &[Limb; NX], d: &[Limb], q: &mut [Limb; NX]) {
    debug_assert_eq!(NX1, NX + 1);
    let mut dn = d.len();
    while d[dn - 1] == 0 {
        dn -= 1;
    }
    debug_assert!(dn >= 2 && dn <= NX);
    let sh = d[dn - 1].leading_zeros();
    let mut dd = [0; NX];
    let mut nn = [0; NX1];
    for i in 0..dn {
        dd[i] = d[i];
    }
    for i in 0..NX {
        nn[i] = x[i];
    }
    shl_window(&mut dd, sh);
    shl_window(&mut nn, sh);
    let v = crate::kernel::div::reciprocal_3by2(dd[dn - 1], dd[dn - 2]);
    *q = [0; NX];
    let qn = NX1 - dn;
    let qh = crate::kernel::div::sb_div_qr(&mut q[..qn], &mut nn, &dd[..dn], v);
    if qn < NX {
        q[qn] = qh;
    } else {
        debug_assert_eq!(qh, 0);
    }
}

/// `x × x` for `x` of `R` limbs into `SQ = 2R` limbs.
#[inline(always)]
fn square<const R: usize, const SQ: usize>(x: &[Limb; R]) -> [Limb; SQ] {
    debug_assert_eq!(SQ, 2 * R);
    let mut p = [0; SQ];
    for i in 0..R {
        let mut carry: Limb = 0;
        for j in 0..R {
            let t = u128::from(x[i]) * u128::from(x[j]) + u128::from(p[i + j]) + u128::from(carry);
            #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
            {
                p[i + j] = t as Limb;
                carry = (t >> Limb::BITS) as Limb;
            }
        }
        p[i + R] = carry;
    }
    p
}

/// `root` with the mantissa padded to `N` limbs, `prec + 1 <= 64 N`.
///
/// `a = M 2^(e - 64 N)` for the `N`-limb `M`, so `a = X 2^(2 h)` with `X = M
/// B^N` and `2 h = e - 128 N`, or `X = 2 M B^N` and `2 h = e - 128 N - 1` when
/// `e` is odd: `X` has `NX = 2N + 1` limbs and `R = floor(sqrt(X))`, between
/// `B^N / sqrt(2)` and `sqrt(2) B^N`, has at least `64 N > prec` bits and
/// `NR = N + 1` limbs.
///
/// Newton's iteration from above (\[BZ10\], chapter 1): `r0 = (isqrt(T) + 1)
/// 2^(t/2)` for the top 128 bits `T = floor(X / 2^t)` (`t` even) exceeds
/// `sqrt(X)` by a factor below `1 + 2^-62`, and `r <- floor((r + floor(X / r))
/// / 2)` never goes below `R` (it is `floor` of the arithmetic mean of `r` and
/// `X / r`, whose geometric mean is `sqrt(X)`) while it roughly squares the
/// relative error: `k` steps with `62 2^k` above `R`'s bits leave `r` within
/// one or two of `R`. Then `r` steps down while `r^2 > X` (at most four times,
/// else the general path takes over), and `X - R^2 != 0` is the sticky bit.
#[inline(always)]
fn root_n<const N: usize, const NX: usize, const NR: usize, const NX1: usize, const SQ: usize>(
    m: &[Limb],
    e: i64,
    prec: u64,
    rnd: Round,
) -> Option<Parts> {
    debug_assert!(NX == 2 * N + 1 && NR == N + 1 && NX1 == NX + 1 && SQ == 2 * NR);
    let odd = e.rem_euclid(2) == 1;
    let mut x = [0; NX];
    for (i, &v) in top_aligned::<N>(m).iter().enumerate() {
        x[N + i] = v;
    }
    if odd {
        shl_window(&mut x, 1);
    }
    let h = (e - 128 * i64::try_from(N).ok()? - i64::from(odd)) / 2;
    // The overestimate from the top 128 bits.
    let bits = 64 * NX as u64 - u64::from(leading_zeros(&x));
    let t = (bits - 128).next_multiple_of(2);
    let mut top = x;
    shr_window(&mut top, t);
    let top = (u128::from(top[1]) << Limb::BITS) | u128::from(top[0]);
    let start = top.isqrt() + 1;
    let mut r = [0; NR];
    #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    {
        r[0] = start as Limb;
        r[1] = (start >> Limb::BITS) as Limb;
    }
    shl_window(&mut r, u32::try_from(t / 2).ok()?);
    // Newton steps: 62 2^k bits against R's at most 64 N + 1.
    let mut good = 62u64;
    while good < 64 * N as u64 + 4 {
        let mut q = [0; NX];
        div_floor::<NX, NX1>(&x, &r, &mut q);
        // r = (r + q) / 2; q <= X / R < sqrt(2) B^N + 1 fits NR limbs.
        let mut carry = false;
        for i in 0..NR {
            let (s, c1) = r[i].overflowing_add(q[i]);
            let (s, c2) = s.overflowing_add(Limb::from(carry));
            r[i] = s;
            carry = c1 || c2;
        }
        shr_window(&mut r, 1);
        if carry {
            r[NR - 1] |= TOP;
        }
        good *= 2;
    }
    // Step down to R, then the remainder's sign.
    let mut x_wide = [0; SQ];
    for i in 0..NX {
        x_wide[i] = x[i];
    }
    let mut steps = 0;
    let rem = loop {
        let mut s = square::<NR, SQ>(&r);
        if cmp(&s, &x_wide) != Ordering::Greater {
            let mut d = x_wide;
            sub_assign(&mut d, &s, false);
            s = d;
            break s;
        }
        steps += 1;
        if steps > 4 {
            debug_assert!(false, "small::root_n: Newton's steps fell short");
            return None;
        }
        let mut borrow = true;
        for v in &mut r {
            let (s, b) = v.overflowing_sub(Limb::from(borrow));
            *v = s;
            borrow = b;
            if !borrow {
                break;
            }
        }
    };
    let sticky = rem.iter().any(|&v| v != 0);
    let lz = leading_zeros(&r);
    debug_assert!(lz == 63 || lz == 64);
    shl_window(&mut r, lz);
    // R 2^h = 0.r × 2^(64 NR - lz + h).
    let exp = h.checked_add(64 * i64::try_from(NR).ok()? - i64::from(lz))?;
    finish(false, &mut r, sticky, exp, prec, rnd)
}
