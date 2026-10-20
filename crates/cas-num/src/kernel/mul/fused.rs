//! Fused linear passes for the Toom evaluations and interpolations, in the
//! manner of GMP's `addlsh_n`, `rsh1add_n` and `divexact_by3`: each pass
//! combines an addition or subtraction with a shift, a multiple or an exact
//! division, so a value is read and written once instead of once per step.
//!
//! The passes run in blocks of [`BLOCK`] limbs: the block's carry chain is
//! unrolled so the carry stays in the flags (as in the kernel's `chain`),
//! then the shifts or divisions run over the block's results. A shift of the
//! result needs the limb above, so those passes write each block one limb
//! behind the one they compute.
//!
//! Most passes work on `w`-limb two's complement values modulo `B^w`, as the
//! interpolations do: carries and borrows out of the top limb are dropped.
//! Their right shifts are logical, so they are exact when the value shifted
//! lies in `[0, B^w)`; the callers' values are nonnegative there.

use super::super::{adc, add_1_assign, sbb, sub_1_assign};
use super::{mul_add, split};
use crate::Limb;

/// Limbs per block.
const BLOCK: usize = 8;

/// Limb `i` of `y << k` (`0 < k < 64`), from limbs `i` (`cur`) and `i - 1`
/// (`prev`) of `y`.
#[inline]
const fn lsh_limb(cur: Limb, prev: Limb, k: u32) -> Limb {
    (cur << k) | (prev >> (Limb::BITS - k))
}

/// `d[t], c = op(a[t], b[t], c)` for `t < K`, unrolled.
#[allow(
    clippy::inline_always,
    reason = "the carry stays in the flags only once inlined"
)]
#[inline(always)]
fn chain_k<const K: usize, C: Copy>(
    a: &[Limb; K],
    b: &[Limb; K],
    mut c: C,
    op: impl Fn(Limb, Limb, C) -> (Limb, C),
) -> ([Limb; K], C) {
    let mut d = [0; K];
    for t in 0..K {
        (d[t], c) = op(a[t], b[t], c);
    }
    (d, c)
}

/// Limbs `i..i + K` of `y << k`, `0 < k < 64`.
#[allow(clippy::inline_always, reason = "part of an unrolled block")]
#[inline(always)]
fn lsh_block<const K: usize>(y: &[Limb], i: usize, k: u32) -> [Limb; K] {
    let mut prev = if i == 0 { 0 } else { y[i - 1] };
    let ys: &[Limb; K] = y[i..i + K].try_into().expect("block length");
    let mut out = [0; K];
    for t in 0..K {
        out[t] = lsh_limb(ys[t], prev, k);
        prev = ys[t];
    }
    out
}

/// Limbs `i..i + K` of `x >> j`, `0 < j < 64`, for `x` zero-extended.
#[allow(clippy::inline_always, reason = "part of an unrolled block")]
#[inline(always)]
fn rsh_block<const K: usize>(x: &[Limb], i: usize, j: u32) -> [Limb; K] {
    let next = x.get(i + K).copied().unwrap_or(0);
    let xs: &[Limb; K] = x[i..i + K].try_into().expect("block length");
    let mut out = [0; K];
    for t in 0..K {
        let above = if t + 1 < K { xs[t + 1] } else { next };
        out[t] = (xs[t] >> j) | (above << (Limb::BITS - j));
    }
    out
}

/// Limb `i` of `v << k`, `k < 64`, for `v` zero-extended.
#[inline]
fn lsh_limb_at(v: &[Limb], i: usize, k: u32) -> Limb {
    let cur = v.get(i).copied().unwrap_or(0);
    if k == 0 {
        return cur;
    }
    let prev = if i == 0 {
        0
    } else {
        v.get(i - 1).copied().unwrap_or(0)
    };
    lsh_limb(cur, prev, k)
}

/// `r = (a << j) + (b << k)` for `j, k < 64`, with `a` and `b` zero-extended
/// to `r.len() > max(a.len(), b.len())`; the sum must fit.
pub(super) fn lsh_add(r: &mut [Limb], a: &[Limb], j: u32, b: &[Limb], k: u32) {
    debug_assert!(j < Limb::BITS && k < Limb::BITS);
    assert!(a.len().max(b.len()) < r.len(), "lsh_add: length mismatch");
    let block = |v: &[Limb], i: usize, s: u32| {
        if s == 0 {
            v[i..i + BLOCK].try_into().expect("block length")
        } else {
            lsh_block::<BLOCK>(v, i, s)
        }
    };
    let m = a.len().min(b.len());
    let (mut i, mut c) = (0, false);
    while i + BLOCK <= m {
        let (xs, ys) = (block(a, i, j), block(b, i, k));
        let out: &mut [Limb; BLOCK] = (&mut r[i..i + BLOCK]).try_into().expect("block length");
        (*out, c) = chain_k(&xs, &ys, c, adc);
        i += BLOCK;
    }
    for (i, r) in r.iter_mut().enumerate().skip(i) {
        (*r, c) = adc(lsh_limb_at(a, i, j), lsh_limb_at(b, i, k), c);
    }
    debug_assert!(!c, "lsh_add: the sum overflows");
}

/// One block of [`rsh_pass`]: limbs `i..i + K` of `op(x, y)` from carry `c`,
/// written shifted right by `j` into limbs `i - 1..i + K - 1`, after `prev`,
/// the limb of `op(x, y)` at `i - 1`.
#[allow(clippy::inline_always, reason = "part of an unrolled block")]
#[allow(clippy::too_many_arguments, reason = "the state of one pass")]
#[inline(always)]
fn rsh_pass_block<const K: usize, C: Copy>(
    x: &mut [Limb],
    y: &[Limb],
    i: usize,
    j: u32,
    prev: &mut Limb,
    c: &mut C,
    op: impl Fn(Limb, Limb, C) -> (Limb, C),
) {
    let xs: &[Limb; K] = x[i..i + K].try_into().expect("block length");
    let ys: &[Limb; K] = y[i..i + K].try_into().expect("block length");
    let d;
    (d, *c) = chain_k(xs, ys, *c, op);
    let out: &mut [Limb; K] = (&mut x[i - 1..i - 1 + K]).try_into().expect("block length");
    let mut p = *prev;
    for t in 0..K {
        out[t] = (p >> j) | (d[t] << (Limb::BITS - j));
        p = d[t];
    }
    *prev = p;
}

/// `x = (op(x, y) mod B^w) >> j` for `w = x.len() == y.len()` and
/// `0 < j < 64`, where `op` chains a carry from `c` through the limbs.
#[allow(clippy::inline_always, reason = "one copy per operation")]
#[inline(always)]
fn rsh_pass<C: Copy>(
    x: &mut [Limb],
    y: &[Limb],
    j: u32,
    c: C,
    op: impl Fn(Limb, Limb, C) -> (Limb, C) + Copy,
) {
    debug_assert!(0 < j && j < Limb::BITS);
    assert_eq!(x.len(), y.len(), "rsh_pass: length mismatch");
    let w = x.len();
    if w == 0 {
        return;
    }
    let (mut prev, mut c) = op(x[0], y[0], c);
    let mut i = 1;
    while i + BLOCK <= w {
        rsh_pass_block::<BLOCK, C>(x, y, i, j, &mut prev, &mut c, op);
        i += BLOCK;
    }
    while i < w {
        rsh_pass_block::<1, C>(x, y, i, j, &mut prev, &mut c, op);
        i += 1;
    }
    x[w - 1] = prev >> j;
}

/// `x = ((y + x) mod B^w) >> j` for equal lengths `w` and `0 < j < 64`.
pub(super) fn rsh_add_assign(x: &mut [Limb], y: &[Limb], j: u32) {
    rsh_pass(x, y, j, false, |x, y, c| adc(y, x, c));
}

/// `x = ((y - x) mod B^w) >> j` for equal lengths `w` and `0 < j < 64`.
pub(super) fn rsh_rsub_assign(x: &mut [Limb], y: &[Limb], j: u32) {
    rsh_pass(x, y, j, false, |x, y, c| sbb(y, x, c));
}

/// `x = ((x + y m) mod B^w) >> 1` for equal lengths `w`.
pub(super) fn addmul_1_rsh1(x: &mut [Limb], y: &[Limb], m: Limb) {
    rsh_pass(x, y, 1, 0, |x, y, carry: Limb| {
        // x + y m + carry <= (2^64 - 1) + (2^64 - 1)^2 + (2^64 - 1) < 2^128.
        split(u128::from(y) * u128::from(m) + u128::from(x) + u128::from(carry))
    });
}

/// `x += y << k` (`ADD`) or `x -= y << k` over `y.len() <= x.len()` limbs,
/// `0 < k < 64`; returns the carry or borrow out of the top of `x`, which
/// includes the bits of `y << k` above it.
fn lsh_pass<const ADD: bool>(x: &mut [Limb], y: &[Limb], k: u32) -> Limb {
    debug_assert!(0 < k && k < Limb::BITS);
    assert!(y.len() <= x.len(), "lsh_pass: length mismatch");
    let op = |x, y, c| if ADD { adc(x, y, c) } else { sbb(x, y, c) };
    let block = |x: &mut [Limb], i: usize, c_in: bool| {
        let ys = lsh_block::<BLOCK>(y, i, k);
        let xs: &mut [Limb; BLOCK] = (&mut x[i..i + BLOCK]).try_into().expect("block length");
        let c_out;
        (*xs, c_out) = chain_k(xs, &ys, c_in, op);
        c_out
    };
    let ly = y.len();
    let (mut i, mut c) = (0, false);
    while i + BLOCK <= ly {
        c = block(x, i, c);
        i += BLOCK;
    }
    while i < ly {
        let [ys] = lsh_block::<1>(y, i, k);
        (x[i], c) = op(x[i], ys, c);
        i += 1;
    }
    let top = y.last().map_or(0, |&t| t >> (Limb::BITS - k)) + Limb::from(c);
    if x.len() == ly {
        top
    } else if ADD {
        add_1_assign(&mut x[ly..], top)
    } else {
        sub_1_assign(&mut x[ly..], top)
    }
}

/// `x += y << k` for `y.len() <= x.len()` and `0 < k < 64`; the sum must fit.
pub(super) fn addlsh_assign(x: &mut [Limb], y: &[Limb], k: u32) {
    let out = lsh_pass::<true>(x, y, k);
    debug_assert_eq!(out, 0, "addlsh_assign: the sum overflows");
}

/// `x -= y << k` modulo `B^w`, `w = x.len() >= y.len()`, `0 < k < 64`.
pub(super) fn sublsh_assign(x: &mut [Limb], y: &[Limb], k: u32) {
    lsh_pass::<false>(x, y, k);
}

/// Exact division by a divisor `d` of `B - 1`, limb by limb from the bottom
/// (GMP's `mpn_bdiv_dbm1c`). With `m = (B - 1) / d`, the quotient `q = a m /
/// (B - 1)` satisfies `q B = q + a m`, so each quotient limb is the previous
/// one minus the next limb of `a m`. The products are independent of each
/// other, and the recurrence carries only a subtraction, unlike the
/// multiplication by `d^-1` whose borrow waits on two products per limb.
///
/// Fed every limb of `a`, it yields `a / d mod B^len` for any `a` divisible by
/// `d` modulo `B^len`: two's complement values work too.
#[derive(Clone, Copy)]
struct Dbm1 {
    m: Limb,
    h: Limb,
}

impl Dbm1 {
    const fn new(d: Limb) -> Self {
        debug_assert!(Limb::MAX % d == 0);
        Self {
            m: Limb::MAX / d,
            h: 0,
        }
    }

    #[inline]
    fn step(&mut self, a: Limb) -> Limb {
        let (p0, p1) = mul_add(a, self.m, 0);
        let (q, borrow) = self.h.overflowing_sub(p0);
        self.h = q.wrapping_sub(p1).wrapping_sub(Limb::from(borrow));
        q
    }
}

/// The exact divisions' driver: `x = (lhs - rhs) / d` modulo `B^w` for equal
/// lengths `w`, where `|x, y, i| operands` gives limbs `i..i + K` of `lhs`
/// and `rhs` as blocks (expanded for `K` of [`BLOCK`] and 1, so it may read
/// `x` up to `i + K`, the limb above the block); each limb of the
/// difference goes through `div`.
macro_rules! divexact_pass {
    ($x:expr, $y:expr, |$xv:ident, $yv:ident, $iv:ident| $operands:expr, $div:expr) => {{
        let (x, y): (&mut [Limb], &[Limb]) = ($x, $y);
        assert_eq!(x.len(), y.len(), "divexact pass: length mismatch");
        let w = x.len();
        let mut div = $div;
        let (mut i, mut c) = (0, false);
        while i + BLOCK <= w {
            let (lhs, rhs): ([Limb; BLOCK], [Limb; BLOCK]) = {
                let ($xv, $yv, $iv): (&[Limb], &[Limb], usize) = (&*x, y, i);
                $operands
            };
            let d;
            (d, c) = chain_k(&lhs, &rhs, c, sbb);
            for (x, d) in x[i..i + BLOCK].iter_mut().zip(d) {
                *x = div(d);
            }
            i += BLOCK;
        }
        while i < w {
            let ([lhs], [rhs]): ([Limb; 1], [Limb; 1]) = {
                let ($xv, $yv, $iv): (&[Limb], &[Limb], usize) = (&*x, y, i);
                $operands
            };
            let d;
            (d, c) = sbb(lhs, rhs, c);
            x[i] = div(d);
            i += 1;
        }
    }};
}

/// Limbs `i..i + K` of `x` or `y` as a block.
#[allow(clippy::inline_always, reason = "part of an unrolled block")]
#[inline(always)]
fn block<const K: usize>(v: &[Limb], i: usize) -> [Limb; K] {
    v[i..i + K].try_into().expect("block length")
}

/// `x = ((x >> j) - y) / 3` modulo `B^w` for equal lengths `w` and
/// `0 < j < 64`; the difference must be a multiple of 3 modulo `B^w`.
pub(super) fn rsh_sub_divexact_by3(x: &mut [Limb], y: &[Limb], j: u32) {
    debug_assert!(0 < j && j < Limb::BITS);
    let mut d3 = Dbm1::new(3);
    let div = |d| d3.step(d);
    divexact_pass!(x, y, |x, y, i| (rsh_block(x, i, j), block(y, i)), div);
}

/// `x = (y - x) / 15` modulo `B^w` for equal lengths `w`; the difference must
/// be a multiple of 15 modulo `B^w`.
pub(super) fn rsub_divexact_by15(x: &mut [Limb], y: &[Limb]) {
    let mut d15 = Dbm1::new(15);
    let div = |d| d15.step(d);
    divexact_pass!(x, y, |x, y, i| (block(y, i), block(x, i)), div);
}

/// `x = (x - (y << 3)) / 9` modulo `B^w` for equal lengths `w` (the bits of
/// `y << 3` above `B^w` are dropped); the difference must be a multiple of 9
/// modulo `B^w`. 9 does not divide `B - 1`, so this divides by 3 twice, both
/// chains in the same pass.
pub(super) fn sublsh3_divexact_by9(x: &mut [Limb], y: &[Limb]) {
    let (mut first, mut second) = (Dbm1::new(3), Dbm1::new(3));
    let div = |d| second.step(first.step(d));
    divexact_pass!(x, y, |x, y, i| (block(x, i), lsh_block(y, i, 3)), div);
}

#[cfg(test)]
mod tests;
