//! Division kernels: by a single limb with a precomputed reciprocal
//! (Möller–Granlund), schoolbook long division (Knuth's algorithm D with a
//! 3-by-2 quotient estimate), divide and conquer (Burnikel–Ziegler) for large
//! divisors, Newton's reciprocal with Barrett's estimates ([`mu`]) for large
//! divisors and quotients, and exact division by a single limb (Hensel).
//!
//! The multi-limb algorithms work on a normalised divisor (top bit set);
//! [`div_rem`] shifts both operands into caller-provided scratch first.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (q, r, n, d, u, v)"
)]

use super::mul::{mul_add, submul_1};
use super::{
    add_1_assign, add_n_assign, cmp_n, mul, mul_scratch_len, shl, shr, sub_1_assign, sub_n_assign,
};
use crate::Limb;

/// Smallest divisor length at which [`div_rem`] uses divide and conquer
/// instead of schoolbook. Set by the tuning configuration (`build.rs`).
pub(crate) const DC_DIV_THRESHOLD: usize = crate::tuning::DC_DIV_THRESHOLD;

pub(crate) mod bdiv;
pub(crate) mod mu;

use mu::{
    MU_DIV_THRESHOLD, mu_div_qr, mu_div_qr_scratch_len, mu_divappr_q, mu_divappr_q_scratch_len,
};

/// Counters for the rarely taken correction branches, so tests can check that
/// their inputs reach them.
#[cfg(test)]
pub(crate) mod hits {
    use std::cell::Cell;

    thread_local! {
        /// Schoolbook: the 3-by-2 estimate was one too large (add-back step).
        pub(crate) static ADD_BACK: Cell<u64> = const { Cell::new(0) };
        /// Schoolbook: the top limbs equal the divisor's, so the estimate is
        /// capped at `B - 1`.
        pub(crate) static CAPPED: Cell<u64> = const { Cell::new(0) };
        /// Divide and conquer: a block quotient estimate was too large.
        pub(crate) static DC_CORRECTION: Cell<u64> = const { Cell::new(0) };
        /// Newton's reciprocal: the recursive reciprocal was decreased.
        pub(crate) static INV_DECREMENT: Cell<u64> = const { Cell::new(0) };
        /// Barrett division: a block quotient estimate was too small.
        pub(crate) static MU_CORRECTION: Cell<u64> = const { Cell::new(0) };
        /// Barrett division: a block remainder reached `B^m - 1`, so its
        /// residue alone was ambiguous.
        pub(crate) static MU_WRAPPED: Cell<u64> = const { Cell::new(0) };
        /// An approximate quotient's last block reached past its limbs and
        /// carried into the ones above.
        pub(crate) static APPR_CARRY: Cell<u64> = const { Cell::new(0) };
    }
}

macro_rules! hit {
    ($counter:ident) => {
        #[cfg(test)]
        hits::$counter.with(|c| c.set(c.get() + 1));
    };
}

/// Möller–Granlund word-level division primitives, generic over the word
/// width so tests can also run them on 8- and 16-bit words. There random
/// inputs reach every boundary case, including equalities that have
/// probability about 2^-64 with 64-bit limbs.
macro_rules! word_division {
    ($(#[$attr:meta])* $vis:vis mod $name:ident: $w:ty, $dw:ty) => {
        $(#[$attr])*
        $vis mod $name {
            #![allow(clippy::cast_possible_truncation, reason = "splitting double words")]

            type W = $w;
            type D = $dw;

            #[inline]
            pub(crate) fn double(hi: W, lo: W) -> D {
                (D::from(hi) << W::BITS) | D::from(lo)
            }

            #[inline]
            pub(crate) fn halves(x: D) -> (W, W) {
                ((x >> W::BITS) as W, x as W)
            }

            /// `floor((B^2 - 1) / d) - B` for a normalised `d` (top bit set):
            /// the reciprocal used by [`div_2by1`].
            pub(crate) fn reciprocal_2by1(d: W) -> W {
                debug_assert!(d >> (W::BITS - 1) == 1, "divisor not normalised");
                // The quotient lies in [B, 2B); dropping the top bit subtracts B.
                halves(D::MAX / D::from(d)).1
            }

            /// `(floor(u / d), u mod d)` for the double word `u = (u1, u0)` with
            /// `u1 < d`, `d` normalised and `v = reciprocal_2by1(d)`
            /// (Möller–Granlund, algorithm 4).
            #[inline]
            pub(crate) fn div_2by1(u1: W, u0: W, d: W, v: W) -> (W, W) {
                debug_assert!(u1 < d);
                let (mut q1, q0) = halves(D::from(v) * D::from(u1) + double(u1, u0));
                q1 = q1.wrapping_add(1);
                let mut r = u0.wrapping_sub(q1.wrapping_mul(d));
                if r > q0 {
                    q1 = q1.wrapping_sub(1);
                    r = r.wrapping_add(d);
                }
                if r >= d {
                    q1 += 1;
                    r -= d;
                }
                (q1, r)
            }

            /// `floor((B^3 - 1) / (d1 B + d0)) - B` for normalised `d1`: the
            /// reciprocal used by [`div_3by2`] (Möller–Granlund, algorithm 6).
            pub(crate) fn reciprocal_3by2(d1: W, d0: W) -> W {
                let mut v = reciprocal_2by1(d1);
                let mut p = d1.wrapping_mul(v).wrapping_add(d0);
                if p < d0 {
                    v = v.wrapping_sub(1);
                    if p >= d1 {
                        v = v.wrapping_sub(1);
                        p = p.wrapping_sub(d1);
                    }
                    p = p.wrapping_sub(d1);
                }
                let (t1, t0) = halves(D::from(v) * D::from(d0));
                p = p.wrapping_add(t1);
                if p < t1 {
                    v = v.wrapping_sub(1);
                    if double(p, t0) >= double(d1, d0) {
                        v = v.wrapping_sub(1);
                    }
                }
                v
            }

            /// `(floor(u / d), u mod d)` for `u = (u2, u1, u0)` and `d = (d1, d0)`
            /// with `(u2, u1) < (d1, d0)`, `d1` normalised and
            /// `v = reciprocal_3by2(d1, d0)` (Möller–Granlund, algorithm 5). The
            /// remainder is a double word.
            #[inline]
            pub(crate) fn div_3by2(u2: W, u1: W, u0: W, d1: W, d0: W, v: W) -> (W, D) {
                let d = double(d1, d0);
                debug_assert!(double(u2, u1) < d);
                let (mut q1, q0) = halves(D::from(v) * D::from(u2) + double(u2, u1));
                let r1 = u1.wrapping_sub(q1.wrapping_mul(d1));
                let t = D::from(d0) * D::from(q1);
                let mut r = double(r1, u0).wrapping_sub(t).wrapping_sub(d);
                q1 = q1.wrapping_add(1);
                if halves(r).0 >= q0 {
                    q1 = q1.wrapping_sub(1);
                    r = r.wrapping_add(d);
                }
                if r >= d {
                    q1 += 1;
                    r -= d;
                }
                (q1, r)
            }
        }
    };
}

word_division!(pub(crate) mod word: u64, u128);
word_division!(#[cfg(test)] pub(crate) mod word8: u8, u16);
word_division!(#[cfg(test)] pub(crate) mod word16: u16, u32);

pub(crate) use word::{div_2by1, div_3by2, reciprocal_2by1, reciprocal_3by2};

/// `q = a / d` for a single non-zero limb `d`, with `q.len() == a.len()`;
/// returns the remainder.
pub fn div_rem_1(q: &mut [Limb], a: &[Limb], d: Limb) -> Limb {
    assert_eq!(q.len(), a.len(), "div_rem_1: length mismatch");
    assert!(d != 0, "div_rem_1: division by zero");
    let shift = d.leading_zeros();
    let d = d << shift;
    let v = reciprocal_2by1(d);
    let mut r = 0;
    if shift == 0 {
        for (q, &a) in q.iter_mut().zip(a).rev() {
            (*q, r) = div_2by1(r, a, d, v);
        }
        return r;
    }
    // Divide a << shift, shifting on the fly; the top limb's high bits start r.
    let back = Limb::BITS - shift;
    if let Some(&top) = a.last() {
        r = top >> back;
    }
    for i in (0..a.len()).rev() {
        let lower = if i > 0 { a[i - 1] >> back } else { 0 };
        (q[i], r) = div_2by1(r, (a[i] << shift) | lower, d, v);
    }
    r >> shift
}

/// `a mod d` for a single non-zero limb `d`.
#[must_use]
pub fn rem_1(a: &[Limb], d: Limb) -> Limb {
    assert!(d != 0, "rem_1: division by zero");
    let shift = d.leading_zeros();
    let d = d << shift;
    let v = reciprocal_2by1(d);
    let back = Limb::BITS - shift;
    let mut r = match (shift, a.last()) {
        (0, _) | (_, None) => 0,
        (_, Some(&top)) => top >> back,
    };
    for i in (0..a.len()).rev() {
        let lower = if shift > 0 && i > 0 {
            a[i - 1] >> back
        } else {
            0
        };
        r = div_2by1(r, (a[i] << shift) | lower, d, v).1;
    }
    r >> shift
}

/// `q = a / d` for a single non-zero limb `d` that divides `a`, by Hensel
/// (2-adic) division, with `q.len() == a.len()`. Returns whether `d` divides
/// `a`; when it does not, `q` holds an unspecified value.
pub fn divexact_1(q: &mut [Limb], a: &[Limb], d: Limb) -> bool {
    assert_eq!(q.len(), a.len(), "divexact_1: length mismatch");
    assert!(d != 0, "divexact_1: division by zero");
    let shift = d.trailing_zeros();
    if shift > 0 && a.first().is_some_and(|&lo| lo & ((1 << shift) - 1) != 0) {
        return false;
    }
    let odd = d >> shift;
    // Newton iteration doubles the correct low bits: 3 (odd^2 = 1 mod 8), 6, …, 96.
    let mut inv = odd;
    for _ in 0..5 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(odd.wrapping_mul(inv)));
    }
    debug_assert_eq!(odd.wrapping_mul(inv), 1);

    // d q_i = s_i + hi_i B with s_i = a_i - c_i (+ B on a wrap), so
    // d Q = A + c_n B^n: the division is exact iff the final c_n is zero.
    let mut c = 0;
    for i in 0..a.len() {
        let x = if shift == 0 {
            a[i]
        } else {
            let upper = a.get(i + 1).map_or(0, |&u| u << (Limb::BITS - shift));
            (a[i] >> shift) | upper
        };
        let (s, wrapped) = x.overflowing_sub(c);
        q[i] = s.wrapping_mul(inv);
        c = mul_add(q[i], odd, 0).1 + Limb::from(wrapped);
    }
    c == 0
}

/// Schoolbook division (Knuth's algorithm D) of `n` by the normalised `d`
/// (`n.len() >= d.len() >= 2`), with `q.len() == n.len() - d.len()` and
/// `v = reciprocal_3by2` of `d`'s top two limbs.
///
/// Leaves the remainder in `n[..d.len()]` and zeroes the limbs above it;
/// returns the top quotient limb (0 or 1), which does not fit in `q`.
pub(crate) fn sb_div_qr(q: &mut [Limb], n: &mut [Limb], d: &[Limb], v: Limb) -> Limb {
    let m = d.len();
    let qn = n.len() - m;
    debug_assert!(m >= 2 && q.len() == qn);
    let qh = cmp_n(&n[qn..], d).is_ge();
    if qh {
        sub_n_assign(&mut n[qn..], d);
    }
    let (d1, d0) = (d[m - 1], d[m - 2]);
    for i in (0..qn).rev() {
        // The window's top m limbs are below d, so its quotient is one limb.
        let w = &mut n[i..=i + m];
        let (u2, u1, u0) = (w[m], w[m - 1], w[m - 2]);
        let mut qhat = if u2 == d1 && u1 == d0 {
            // The estimate would be B; B - 1 is then exact.
            hit!(CAPPED);
            Limb::MAX
        } else {
            // At most one too large.
            div_3by2(u2, u1, u0, d1, d0, v).0
        };
        let borrow = submul_1(&mut w[..m], d, qhat);
        let (top, negative) = w[m].overflowing_sub(borrow);
        w[m] = top;
        if negative {
            hit!(ADD_BACK);
            qhat -= 1;
            let carry = add_n_assign(&mut w[..m], d);
            w[m] = w[m].wrapping_add(carry);
        }
        debug_assert_eq!(w[m], 0);
        q[i] = qhat;
    }
    Limb::from(qh)
}

/// Computes the `k = q.len()` quotient limbs of the window `w` (`d.len() + k`
/// limbs, `1 <= k <= d.len()`) divided by the normalised `d`: estimates them
/// from the top `2k` limbs and the top `k` limbs of `d`, then subtracts the
/// estimate times the rest of `d` and corrects.
///
/// Leaves the remainder in `w[..d.len()]` (limbs above it are unspecified);
/// returns the top quotient limb (0 or 1).
fn div_block(q: &mut [Limb], w: &mut [Limb], d: &[Limb], v: Limb, scratch: &mut [Limb]) -> Limb {
    let (k, dn) = (q.len(), d.len());
    debug_assert!(1 <= k && k <= dn && w.len() == dn + k);
    if k == 1 {
        return sb_div_qr(q, w, d, v);
    }
    let lo = dn - k;
    let mut qh = if k < DC_DIV_THRESHOLD {
        sb_div_qr(q, &mut w[lo..], &d[lo..], v)
    } else {
        dc_div_qr_n(q, &mut w[lo..], &d[lo..], v, scratch)
    };
    if lo == 0 {
        return qh;
    }
    // w[..dn] -= (qh B^k + q) * d[..lo]
    let d_lo = &d[..lo];
    let (tp, rest) = scratch.split_at_mut(dn);
    mul(tp, q, d_lo, rest);
    let mut borrows = sub_n_assign(&mut w[..dn], tp);
    if qh != 0 {
        borrows += sub_n_assign(&mut w[k..dn], d_lo);
    }
    while borrows != 0 {
        hit!(DC_CORRECTION);
        qh -= sub_1_assign(q, 1);
        borrows -= add_n_assign(&mut w[..dn], d);
    }
    qh
}

/// Divide and conquer division of the `2k`-limb `w` by the normalised `k`-limb
/// `d` (`k >= DC_DIV_THRESHOLD`): two [`div_block`] steps of half size each.
/// Leaves the remainder in `w[..k]`; returns the top quotient limb.
pub(crate) fn dc_div_qr_n(
    q: &mut [Limb],
    w: &mut [Limb],
    d: &[Limb],
    v: Limb,
    scratch: &mut [Limb],
) -> Limb {
    let k = d.len();
    debug_assert!(w.len() == 2 * k && q.len() == k);
    let lo = k / 2;
    let (q_lo, q_hi) = q.split_at_mut(lo);
    let qh = div_block(q_hi, &mut w[lo..], d, v, scratch);
    let ql = div_block(q_lo, &mut w[..k + lo], d, v, scratch);
    debug_assert_eq!(ql, 0);
    qh
}

/// [`div_block`]'s quotient alone: the estimate from the top `2k` limbs of
/// `w` and the top `k` limbs of `d`, without the correction by the rest of
/// `d`. That estimate is never below the block's quotient and at most 2
/// above when exact (as [`div_block`] relies on); for `k` from divide and
/// conquer's threshold it comes from [`dc_divappr_n`], which adds its own
/// excess. Returns the top quotient limb; `w` is left unspecified.
fn div_block_appr(
    q: &mut [Limb],
    w: &mut [Limb],
    d: &[Limb],
    v: Limb,
    scratch: &mut [Limb],
) -> Limb {
    let (k, dn) = (q.len(), d.len());
    debug_assert!(1 <= k && k <= dn && w.len() == dn + k);
    if k == 1 {
        return sb_div_qr(q, w, d, v);
    }
    let lo = dn - k;
    if k < DC_DIV_THRESHOLD {
        sb_div_qr(q, &mut w[lo..], &d[lo..], v)
    } else {
        dc_divappr_n(q, &mut w[lo..], &d[lo..], v, scratch)
    }
}

/// [`dc_div_qr_n`]'s quotient alone: the high half exactly, the low half by
/// [`div_block_appr`]. The quotient is never below the true one, and above
/// it by at most 2 per level of the recursion below this one that reaches
/// divide and conquer's threshold, under 2 log₂ k in all. Returns the top
/// quotient limb; `w` is left unspecified.
fn dc_divappr_n(q: &mut [Limb], w: &mut [Limb], d: &[Limb], v: Limb, scratch: &mut [Limb]) -> Limb {
    let k = d.len();
    debug_assert!(w.len() == 2 * k && q.len() == k);
    let lo = k / 2;
    let (q_lo, q_hi) = q.split_at_mut(lo);
    let qh = div_block(q_hi, &mut w[lo..], d, v, scratch);
    let ql = div_block_appr(q_lo, &mut w[..k + lo], d, v, scratch);
    if ql != 0 {
        hit!(APPR_CARRY);
    }
    qh + add_1_assign(q_hi, ql)
}

fn div_block_scratch_len(k: usize, dn: usize) -> usize {
    if k == 1 {
        return 0;
    }
    let lo = dn - k;
    let inner = if k < DC_DIV_THRESHOLD {
        0
    } else {
        dc_div_qr_n_scratch_len(k)
    };
    let correction = if lo == 0 {
        0
    } else {
        dn + mul_scratch_len(k, lo)
    };
    inner.max(correction)
}

/// Scratch limbs [`dc_div_qr_n`] needs for `k` limbs: the most that one
/// correction product of its recursion takes, `s + mul_scratch_len(j, s -
/// j)` for a block of `j` of an `s`-limb division's quotient limbs (none for
/// `j = 1`). At depth `t` the divisions have `k >> t` limbs or one more, and
/// both halves of one take the same scratch, so one pass over the depths
/// stands for the whole tree, whose size grows like `k`; `div_rem` checks
/// its scratch on every call.
pub(crate) fn dc_div_qr_n_scratch_len(k: usize) -> usize {
    let node = |s: usize| {
        if s < 3 {
            0
        } else {
            s + mul_scratch_len(s - s / 2, s / 2)
        }
    };
    let mut len = node(k);
    let (mut lo, mut hi) = (k / 2, k.div_ceil(2));
    while hi >= DC_DIV_THRESHOLD {
        len = len.max(node(hi));
        if lo != hi && lo >= DC_DIV_THRESHOLD {
            len = len.max(node(lo));
        }
        (lo, hi) = (lo / 2, hi.div_ceil(2));
    }
    len
}

/// Whether [`div_qr_normalized`] divides by Newton's reciprocal.
fn uses_mu(nn: usize, m: usize) -> bool {
    m >= MU_DIV_THRESHOLD && nn - m >= MU_DIV_THRESHOLD
}

/// Division of `n` by the normalised `d` (`d.len() >= 2`) with
/// `q.len() == n.len() - d.len()` and `v = reciprocal_3by2` of `d`'s top two
/// limbs: Barrett's estimates from Newton's reciprocal ([`mu::mu_div_qr`])
/// when both the divisor and the quotient are long, else [`bz_div_qr`].
/// Leaves the remainder in `n[..d.len()]`; returns the top quotient limb.
pub(crate) fn div_qr_normalized(
    q: &mut [Limb],
    n: &mut [Limb],
    d: &[Limb],
    v: Limb,
    scratch: &mut [Limb],
) -> Limb {
    if uses_mu(n.len(), d.len()) {
        mu_div_qr(q, n, d, scratch)
    } else {
        bz_div_qr(q, n, d, v, scratch)
    }
}

pub(crate) fn div_qr_normalized_scratch_len(nn: usize, m: usize) -> usize {
    if uses_mu(nn, m) {
        mu_div_qr_scratch_len(nn, m)
    } else {
        bz_div_qr_scratch_len(nn, m)
    }
}

/// [`div_qr_normalized`] without Newton's reciprocal: schoolbook for short
/// divisors, otherwise quotient blocks of `d.len()` limbs from the top, each
/// by divide and conquer.
pub(crate) fn bz_div_qr(
    q: &mut [Limb],
    n: &mut [Limb],
    d: &[Limb],
    v: Limb,
    scratch: &mut [Limb],
) -> Limb {
    let m = d.len();
    let qn = n.len() - m;
    debug_assert!(m >= 2 && q.len() == qn);
    debug_assert_eq!(v, reciprocal_3by2(d[m - 1], d[m - 2]));
    if m < DC_DIV_THRESHOLD || qn < 2 {
        return sb_div_qr(q, n, d, v);
    }
    let qh = cmp_n(&n[qn..], d).is_ge();
    if qh {
        sub_n_assign(&mut n[qn..], d);
    }
    // The partial block first, then full blocks; each window's top m limbs
    // are the previous remainder, below d, so no block has a top limb.
    let mut top = qn;
    let first = qn % m;
    let mut k = if first == 0 { m } else { first };
    while top > 0 {
        let start = top - k;
        let ql = div_block(&mut q[start..top], &mut n[start..top + m], d, v, scratch);
        debug_assert_eq!(ql, 0);
        top = start;
        k = m;
    }
    Limb::from(qh)
}

/// [`bz_div_qr`]'s quotient alone: the blocks as there, the last (lowest)
/// one by [`div_block_appr`]. The quotient is never below the true one and
/// at most `2 log₂ m` above it. Returns the top quotient limb; `n` is left
/// unspecified.
pub(crate) fn bz_divappr_q(
    q: &mut [Limb],
    n: &mut [Limb],
    d: &[Limb],
    v: Limb,
    scratch: &mut [Limb],
) -> Limb {
    let m = d.len();
    let qn = n.len() - m;
    debug_assert!(m >= 2 && q.len() == qn);
    if m < DC_DIV_THRESHOLD || qn < 2 {
        return sb_div_qr(q, n, d, v);
    }
    let mut qh = Limb::from(cmp_n(&n[qn..], d).is_ge());
    if qh != 0 {
        sub_n_assign(&mut n[qn..], d);
    }
    let mut top = qn;
    let first = qn % m;
    let mut k = if first == 0 { m } else { first };
    while top > 0 {
        let start = top - k;
        let (q_block, w) = (&mut q[start..top], &mut n[start..top + m]);
        if start > 0 {
            let ql = div_block(q_block, w, d, v, scratch);
            debug_assert_eq!(ql, 0);
        } else {
            let ql = div_block_appr(q_block, w, d, v, scratch);
            if ql != 0 {
                hit!(APPR_CARRY);
            }
            qh += add_1_assign(&mut q[top..], ql);
        }
        top = start;
        k = m;
    }
    qh
}

pub(crate) fn bz_div_qr_scratch_len(nn: usize, m: usize) -> usize {
    let qn = nn - m;
    if m < DC_DIV_THRESHOLD || qn < 2 {
        return 0;
    }
    let first = qn % m;
    let partial = if first == 0 {
        0
    } else {
        div_block_scratch_len(first, m)
    };
    let full = if qn >= m {
        div_block_scratch_len(m, m)
    } else {
        0
    };
    partial.max(full)
}

/// The lengths [`divappr_q`] works with for `nn` dividend and `dn` divisor
/// limbs: the divisor's limbs it keeps, `dt`, and how many it drops, `cut`.
fn divappr_lens(nn: usize, dn: usize) -> (usize, usize) {
    let dt = dn.min(nn - dn + 2);
    (dt, dn - dt)
}

/// Whether [`divappr_q`] computes a guard limb and divides approximately,
/// rather than dividing its truncated operands exactly.
fn divappr_guards(dt: usize) -> bool {
    dt >= DC_DIV_THRESHOLD
}

/// Scratch limbs [`divappr_q`] needs for `nn` dividend and `dn` divisor limbs.
pub(crate) fn divappr_q_scratch_len(nn: usize, dn: usize) -> usize {
    if dn <= 1 {
        return 0;
    }
    let (dt, cut) = divappr_lens(nn, dn);
    // The guarded dividend: one limb below the truncated one, which has
    // the shift's extra limb on top.
    let gn = nn - cut + 2;
    let inner = if !divappr_guards(dt) {
        div_qr_normalized_scratch_len(gn - 1, dt)
    } else if uses_mu(gn, dt) {
        mu_divappr_q_scratch_len(gn, dt)
    } else {
        bz_div_qr_scratch_len(gn, dt)
    };
    (dt + 1) + gn + (gn - dt) + inner
}

/// An approximate quotient: `q` with `Q <= q <= Q + 1` for `Q = ⌊n / d⌋`,
/// for `d` with a non-zero top limb, `n.len() >= d.len()`,
/// `q.len() == n.len() - d.len() + 1` and
/// `scratch.len() >= divappr_q_scratch_len(n.len(), d.len())`. As GMP's
/// `mpn_divappr_q`, it costs less than the exact quotient by skipping the
/// last part of the remainder.
///
/// Only the top `qn + 2` limbs of `d` matter (`qn = q.len() - 1`): with `j`
/// limbs dropped from `n` and from `d` (normalised), the quotient changes
/// by less than `B^j / d < 2 B^(-qn - 2)` relative, below `4 / B` since
/// `Q < B^(qn + 1)`. Below divide and conquer's threshold, the truncated
/// operands are divided exactly: their quotient is `Q` or `Q + 1`. From it,
/// the quotient `G` of the truncated `n` times `B` is estimated, one guard
/// limb more, without the last block's remainder ([`mu_divappr_q`], at most
/// 5 below, or [`bz_divappr_q`], at most `2 log₂ dn` above). `G` is at least
/// `Q B - 1` and below `Q B + B + 4`, so adding 6 to the estimate and
/// dropping the guard limb gives `Q` or `Q + 1`.
pub(crate) fn divappr_q(q: &mut [Limb], n: &[Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (nn, dn) = (n.len(), d.len());
    assert!(
        dn >= 1 && nn >= dn && q.len() == nn - dn + 1,
        "divappr_q: length mismatch"
    );
    assert!(d[dn - 1] != 0, "divappr_q: divisor has a zero top limb");
    assert!(
        scratch.len() >= divappr_q_scratch_len(nn, dn),
        "divappr_q: scratch too short"
    );
    if dn == 1 {
        div_rem_1(q, n, d[0]);
        return;
    }
    let shift = d[dn - 1].leading_zeros();
    let (dt, cut) = divappr_lens(nn, dn);
    // The truncated divisor, normalised: its top dt limbs shifted, with the
    // bits from the limb below.
    let (ds, rest) = scratch.split_at_mut(dt + 1);
    let from = cut.saturating_sub(1);
    let out = shl(&mut ds[..dn - from], &d[from..], shift);
    debug_assert_eq!(out, 0);
    let ds = &ds[dn - from - dt..dn - from];
    // The truncated dividend, shifted the same way, above one limb that is
    // its guard limb, or a limb below it when truncating.
    let gn = nn - cut + 2;
    let (g, rest) = rest.split_at_mut(gn);
    let start = usize::from(cut == 0);
    g[gn - 1] = shl(&mut g[start..gn - 1], &n[from..], shift);
    g[0] = 0;
    let v = reciprocal_3by2(ds[dt - 1], ds[dt - 2]);
    let (gq, rest) = rest.split_at_mut(gn - dt);
    let top = if divappr_guards(dt) {
        let gq = &mut gq[..gn - dt];
        let mut top = if uses_mu(gn, dt) {
            mu_divappr_q(gq, g, ds, rest)
        } else {
            bz_divappr_q(gq, g, ds, v, rest)
        };
        top += add_1_assign(gq, 6);
        q.copy_from_slice(&gq[1..]);
        top
    } else {
        div_qr_normalized(q, &mut g[1..], ds, v, rest)
    };
    if top != 0 {
        // At most Q + 1 <= B^(qn + 1), so Q = B^(qn + 1) - 1.
        q.fill(Limb::MAX);
    }
}

/// Scratch limbs [`div_rem`] needs for a dividend of `an` and a divisor of
/// `dn` limbs.
#[must_use]
pub fn div_rem_scratch_len(an: usize, dn: usize) -> usize {
    if dn <= 1 || an < dn {
        return 0;
    }
    // Shifted divisor, shifted dividend with one more limb, then the algorithm
    // on all of it or on its low `an` limbs (see `div_rem`).
    let algorithm = if uses_mu(an, dn) {
        div_qr_normalized_scratch_len(an, dn)
    } else {
        div_qr_normalized_scratch_len(an + 1, dn).max(div_qr_normalized_scratch_len(an, dn))
    };
    dn + an + 1 + algorithm
}

/// `q = a / d` and `r = a mod d` for `a.len() >= d.len() >= 1` with a non-zero
/// top limb of `d`, `q.len() == a.len() - d.len() + 1`, `r.len() == d.len()` and
/// `scratch.len() >= div_rem_scratch_len(a.len(), d.len())`.
pub fn div_rem(q: &mut [Limb], r: &mut [Limb], a: &[Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (an, dn) = (a.len(), d.len());
    assert!(
        dn >= 1 && an >= dn && q.len() == an - dn + 1 && r.len() == dn,
        "div_rem: length mismatch"
    );
    assert!(d[dn - 1] != 0, "div_rem: divisor has a zero top limb");
    assert!(
        scratch.len() >= div_rem_scratch_len(an, dn),
        "div_rem: scratch too short"
    );
    if dn == 1 {
        r[0] = div_rem_1(q, a, d[0]);
        return;
    }
    let shift = d[dn - 1].leading_zeros();
    let (ds, rest) = scratch.split_at_mut(dn);
    let (ns, rest) = rest.split_at_mut(an + 1);
    let out = shl(ds, d, shift);
    debug_assert_eq!(out, 0);
    let above = shl(&mut ns[..an], a, shift);
    let v = reciprocal_3by2(ds[dn - 1], ds[dn - 2]);
    // The quotient has qn + 1 limbs, the top one zero unless the shift
    // carries a limb out of the dividend. Below Barrett's range that limb
    // joins the dividend, whose quotient then fits in q. In Barrett's range,
    // as in GMP's `mpn_tdiv_qr`, the top quotient limb comes from one
    // schoolbook step on the top dn + 1 limbs instead: its quotient fits a
    // limb, since `above` is below the divisor's top limb, and it leaves the
    // top dn limbs below the divisor. The rest is then a division of `an`
    // limbs, whose returned high limb is the top quotient limb when nothing
    // was carried out. Dividing an + 1 limbs there would put the blocks of
    // Barrett's division, and their transforms, just past a power of two for
    // a 2n/n division.
    let qn = an - dn;
    let mut top = 0;
    let nn = if above == 0 {
        an
    } else {
        ns[an] = above;
        if uses_mu(an, dn) {
            let mut qt = [0];
            let qh = sb_div_qr(&mut qt, &mut ns[qn..], ds, v);
            debug_assert_eq!(qh, 0);
            top = qt[0];
            an
        } else {
            an + 1
        }
    };
    let qh = div_qr_normalized(&mut q[..nn - dn], &mut ns[..nn], ds, v, rest);
    if nn == an {
        debug_assert!(top == 0 || qh == 0);
        q[qn] = top + qh;
    } else {
        debug_assert_eq!(qh, 0);
    }
    shr(r, &ns[..dn], shift);
}
