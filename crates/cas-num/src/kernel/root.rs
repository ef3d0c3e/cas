//! Square roots with remainder: Zimmermann's Karatsuba square root on limbs,
//! as GMP's `mpn_sqrtrem` (Brent and Zimmermann, *Modern Computer
//! Arithmetic*, Algorithm 1.12).
//!
//! [`dc_sqrtrem`] takes the root of a normalised `2n`-limb value in place:
//! the root of its top half, recursively, then one division by that root
//! gives the low half of the root, and one square its remainder. The
//! divisor is the top half of the root, whose top bit is set, so the
//! division needs no shift. [`sqrtrem`] shifts its operand once, by an
//! even number of bits and a zero limb if needed, and shifts the root and
//! the remainder back.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (s, r, a, q, l, h, c, k)"
)]

use super::div::{
    div_2by1, div_qr_normalized, div_qr_normalized_scratch_len, reciprocal_2by1, reciprocal_3by2,
};
use super::{
    add_1_assign, add_n_assign, addmul_1, normalized_len, shl, shr, shr_assign, sqr,
    sqr_scratch_len, sub_1_assign, sub_n_assign, submul_1,
};
use crate::Limb;

// Counters for the rarely taken branches, so tests can check that their
// inputs reach them.
#[cfg(test)]
thread_local! {
    /// The root's estimate was one too large (the correction step).
    pub(crate) static SQRT_CORRECTION: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    /// The low half's estimate was `B^l`, past its limbs.
    pub(crate) static SQRT_CARRY: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

macro_rules! hit {
    ($counter:ident) => {
        #[cfg(test)]
        $counter.with(|c| c.set(c.get() + 1));
    };
}

/// `s = ⌊√a⌋` for `a` of `2n` limbs with a top limb of at least `B/4`, and
/// `s.len() == n`. Leaves the remainder `a - s²` in `a[..n]` and returns
/// its bit `n` (it is at most `2s`, below `2 B^n`); the limbs of `a` above
/// `n` are left unspecified. `scratch.len() >= dc_sqrtrem_scratch_len(n)`.
///
/// With `l = ⌊n/2⌋` and `h = n - l`, the root `s'` of the top `2h` limbs,
/// with remainder `r'`, is the top half of `s`. Dividing `r' B^l` plus the
/// next `l` limbs by `2s'` gives the low half `q` and a remainder `u`;
/// then `s = s' B^l + q` and `a - s² = u B^l + a₀ - q²`, where `a₀` is the
/// low `l` limbs. That is negative at most once, when `s` is one too large.
pub(crate) fn dc_sqrtrem(s: &mut [Limb], a: &mut [Limb], scratch: &mut [Limb]) -> bool {
    let n = s.len();
    debug_assert!(n >= 1 && a.len() == 2 * n && a[2 * n - 1] >= 1 << (Limb::BITS - 2));
    if n == 1 {
        let x = (u128::from(a[1]) << Limb::BITS) | u128::from(a[0]);
        let root = x.isqrt();
        let r = x - root * root;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "a 64-bit root, and the remainder's low limb"
        )]
        let (root, low) = (root as Limb, r as Limb);
        s[0] = root;
        a[0] = low;
        return r >> Limb::BITS != 0;
    }
    let l = n / 2;
    let h = n - l;
    let (s_lo, s_hi) = s.split_at_mut(l);
    // r' up to 2s' has h limbs and a top bit; dividing by s' instead of 2s',
    // r' - s' is below B^h and adds one to the quotient's limb l.
    let mut top = Limb::from(dc_sqrtrem(s_hi, &mut a[2 * l..], scratch));
    if top != 0 {
        let borrow = sub_n_assign(&mut a[2 * l..n + l], s_hi);
        debug_assert_eq!(borrow, 1);
    }
    // s' has its top bit set, as the root of a normalised value.
    let dividend = &mut a[l..n + l];
    top += if h == 1 {
        let d = s_hi[0];
        let qh = dividend[1] >= d;
        if qh {
            dividend[1] -= d;
        }
        (s_lo[0], dividend[0]) = div_2by1(dividend[1], dividend[0], d, reciprocal_2by1(d));
        Limb::from(qh)
    } else {
        let v = reciprocal_3by2(s_hi[h - 1], s_hi[h - 2]);
        div_qr_normalized(s_lo, dividend, s_hi, v, scratch)
    };
    // Halve the quotient: q = Q / 2 for the divisor 2s', whose remainder is
    // the remainder by s', plus s' if Q is odd.
    let odd = s_lo[0] & 1 != 0;
    shr_assign(s_lo, 1);
    s_lo[l - 1] |= top << (Limb::BITS - 1);
    top >>= 1;
    let carry = if odd {
        add_n_assign(&mut a[l..n], s_hi)
    } else {
        0
    };
    // a - s² = u B^l + a₀ - q², with q = top B^l + s_lo; top = 1 leaves
    // s_lo zero and q² = B^(2l).
    let (rem, free) = a.split_at_mut(n);
    sqr(&mut free[..2 * l], s_lo, scratch);
    let borrow = top + sub_n_assign(&mut rem[..2 * l], &free[..2 * l]);
    let borrow = sub_1_assign(&mut rem[2 * l..], borrow);
    if borrow <= carry {
        debug_assert_eq!(top, 0);
        return carry > borrow;
    }
    // s is one too large (as it always is when top = 1, since s' B^l is
    // the root's top half): a - (s - 1)² = a - s² + 2s - 1, where a - s²
    // is rem - B^n.
    hit!(SQRT_CORRECTION);
    if top != 0 {
        hit!(SQRT_CARRY);
    }
    let top = add_1_assign(s_hi, top);
    let up = addmul_1(rem, s, 2) + 2 * top;
    let down = sub_1_assign(rem, 1);
    let out = sub_1_assign(s, 1);
    debug_assert_eq!(top, out);
    debug_assert!(up - down == 1 || up - down == 2);
    up - down == 2
}

/// Scratch limbs [`dc_sqrtrem`] needs for a root of `n` limbs: the largest
/// division or square of its recursion, which halves `n` at each level.
pub(crate) fn dc_sqrtrem_scratch_len(n: usize) -> usize {
    let mut len = 0;
    let mut n = n;
    while n >= 2 {
        let (l, h) = (n / 2, n - n / 2);
        let div = if h >= 2 {
            div_qr_normalized_scratch_len(n, h)
        } else {
            0
        };
        len = len.max(div).max(sqr_scratch_len(l));
        n = h;
    }
    len
}

/// Scratch limbs [`sqrtrem`] needs for an operand of `an` limbs.
#[must_use]
pub fn sqrtrem_scratch_len(an: usize) -> usize {
    let n = an.div_ceil(2);
    2 * n + dc_sqrtrem_scratch_len(n)
}

/// `s = ⌊√a⌋` and `r = a - s²` for `a` with a non-zero top limb, with
/// `s.len() == a.len().div_ceil(2)`, `r.len() == a.len() / 2 + 1` (the
/// remainder is at most `2s`) and `scratch.len() >= sqrtrem_scratch_len(a.len())`.
/// Returns the remainder's length without its zero high limbs, which are
/// written as zeros; it is 0 when `a` is a square.
///
/// The root's top limb is never zero. `a` is shifted left by `2k` bits so
/// that its length is even and its top limb at least `B/4`; the root `S` of
/// that is `s 2^k + s₀`, and `2^(2k) r = R + 2 s₀ S - s₀²` from its
/// remainder `R`.
pub fn sqrtrem(s: &mut [Limb], r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) -> usize {
    let an = a.len();
    assert!(an >= 1 && a[an - 1] != 0, "sqrtrem: operand not normalised");
    let n = an.div_ceil(2);
    assert!(
        s.len() == n && r.len() == an / 2 + 1,
        "sqrtrem: length mismatch"
    );
    assert!(
        scratch.len() >= sqrtrem_scratch_len(an),
        "sqrtrem: scratch too short"
    );
    if an == 1 {
        s[0] = a[0].isqrt();
        r[0] = a[0] - s[0] * s[0];
        return usize::from(r[0] != 0);
    }
    let odd = an & 1;
    let c = a[an - 1].leading_zeros() / 2;
    let (t, rest) = scratch.split_at_mut(2 * n);
    t[0] = 0;
    let out = shl(&mut t[odd..], a, 2 * c);
    debug_assert_eq!(out, 0);
    let mut rl = Limb::from(dc_sqrtrem(s, t, rest));
    // The root of t has k bits more than the root of a.
    let k = c + if odd == 1 { Limb::BITS / 2 } else { 0 };
    if k == 0 {
        r[..n].copy_from_slice(&t[..n]);
        r[n] = rl;
    } else {
        let s0 = s[0] & ((1 << k) - 1);
        rl += addmul_1(&mut t[..n], s, 2 * s0);
        // Nothing is borrowed from rl: t + rl B^n becomes 2^(2k) r, and
        // B^n is a multiple of 2^(2k) too, so before this subtraction the
        // low n limbs held a multiple of 2^(2k) below B^n plus s₀² < 2^(2k).
        let borrow = submul_1(&mut t[..1], &[s0], s0);
        let out = sub_1_assign(&mut t[1..n], borrow);
        debug_assert_eq!(out, 0);
        shr_assign(s, k);
        t[n] = rl;
        // 2k bits come off the remainder: a zero limb first when a was
        // padded with one.
        let (from, shift) = if 2 * k < Limb::BITS {
            (0, 2 * k)
        } else {
            (1, 2 * k - Limb::BITS)
        };
        let out = shr(r, &t[from..=n], shift);
        debug_assert_eq!(out, 0);
    }
    normalized_len(r)
}
