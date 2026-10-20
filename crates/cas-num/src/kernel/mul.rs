//! Multiplication kernels: single-limb products, schoolbook, Karatsuba,
//! Toom-3, Toom-4 and number-theoretic transforms ([`super::ntt`]), plus
//! dedicated squaring.
//!
//! [`mul`] and [`sqr`] pick the algorithm by operand size. The recursive
//! algorithms work in a caller-provided scratch slice of at least
//! [`mul_scratch_len`] / [`sqr_scratch_len`] limbs, so a whole product costs at
//! most one allocation, made by the caller. Scratch contents on entry are
//! irrelevant and are clobbered.
//!
//! Very unbalanced operands (the shorter at most half the longer) are cut into
//! chunks of the shorter operand's length, each multiplied with a balanced
//! algorithm.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, b, r, n, s, t)"
)]

mod fused;
pub(crate) mod mulhi;
pub(crate) mod mullo;
pub(crate) mod split;

use super::{
    add, add_1_assign, add_assign, add_n, add_n_assign, cmp_n, ntt, shl_assign, shr_assign,
    sub_assign, sub_n, sub_n_assign,
};
use crate::Limb;

/// Smallest shorter-operand length at which [`mul`] uses Karatsuba instead of
/// schoolbook. These thresholds come from the tuning configuration (`build.rs`).
pub(crate) const MUL_KARATSUBA_THRESHOLD: usize = crate::tuning::MUL_KARATSUBA_THRESHOLD;
/// Smallest shorter-operand length at which [`mul`] uses Toom-3.
pub(crate) const MUL_TOOM3_THRESHOLD: usize = crate::tuning::MUL_TOOM3_THRESHOLD;
/// Smallest shorter-operand length at which [`mul`] uses Toom-4.
pub(crate) const MUL_TOOM4_THRESHOLD: usize = crate::tuning::MUL_TOOM4_THRESHOLD;
/// Shorter-operand length from which [`mul`] uses the NTT for products whose
/// transform is as cheap as the size allows; see [`ntt_pays`].
pub(crate) const MUL_NTT_THRESHOLD: usize = crate::tuning::MUL_NTT_THRESHOLD;
/// Shorter operand's length from which unbalanced products use Toom-32.
pub(crate) const MUL_TOOM32_THRESHOLD: usize = crate::tuning::MUL_TOOM32_THRESHOLD;
/// Shorter operand's length from which unbalanced products use Toom-42, and
/// longer ones pieces of twice that operand.
pub(crate) const MUL_TOOM42_THRESHOLD: usize = crate::tuning::MUL_TOOM42_THRESHOLD;
/// Smallest length at which [`sqr`] uses Karatsuba squaring.
pub(crate) const SQR_KARATSUBA_THRESHOLD: usize = crate::tuning::SQR_KARATSUBA_THRESHOLD;
/// Smallest length at which [`sqr`] uses Toom-3 squaring.
pub(crate) const SQR_TOOM3_THRESHOLD: usize = crate::tuning::SQR_TOOM3_THRESHOLD;
/// Smallest length at which [`sqr`] uses Toom-4 squaring.
pub(crate) const SQR_TOOM4_THRESHOLD: usize = crate::tuning::SQR_TOOM4_THRESHOLD;
/// Length from which [`sqr`] uses the NTT, as [`MUL_NTT_THRESHOLD`].
pub(crate) const SQR_NTT_THRESHOLD: usize = crate::tuning::SQR_NTT_THRESHOLD;
/// Ten times the exponent of [`ntt_pays`]'s efficiency rule.
pub(crate) const NTT_EFFICIENCY_EXPONENT: usize = crate::tuning::NTT_EFFICIENCY_EXPONENT;
/// Ten times the exponent of [`ntt_pays_unbalanced`]'s Toom cost.
pub(crate) const NTT_UNBALANCED_EXPONENT: usize = crate::tuning::NTT_UNBALANCED_EXPONENT;
/// Smallest power-of-two modulus length at which [`mul_wrap`] uses a cyclic
/// transform instead of reducing the full product.
pub(crate) const MUL_WRAP_NTT_THRESHOLD: usize = crate::tuning::MUL_WRAP_NTT_THRESHOLD;
/// Smallest even modulus length at which [`mul_wrap`], below the cyclic
/// transform, splits `B^m - 1` into `B^(m/2) - 1` and `B^(m/2) + 1`
/// ([`split`]) instead of reducing the full product.
pub(crate) const MUL_WRAP_SPLIT_THRESHOLD: usize = crate::tuning::MUL_WRAP_SPLIT_THRESHOLD;

/// `a * b + c` as a double limb, returned as `(low, high)`; cannot overflow.
#[inline]
pub(super) fn mul_add(a: Limb, b: Limb, c: Limb) -> (Limb, Limb) {
    let t = u128::from(a) * u128::from(b) + u128::from(c);
    split(t)
}

#[inline]
#[expect(clippy::cast_possible_truncation, reason = "splitting a double limb")]
fn split(t: u128) -> (Limb, Limb) {
    (t as Limb, (t >> Limb::BITS) as Limb)
}

/// `r = a * b` for a single limb `b`; returns the high limb of the product.
pub fn mul_1(r: &mut [Limb], a: &[Limb], b: Limb) -> Limb {
    assert_eq!(r.len(), a.len(), "mul_1: length mismatch");
    let mut carry = 0;
    for (r, &a) in r.iter_mut().zip(a) {
        (*r, carry) = mul_add(a, b, carry);
    }
    carry
}

/// `a *= b` for a single limb `b`; returns the high limb of the product.
pub fn mul_1_assign(a: &mut [Limb], b: Limb) -> Limb {
    let mut carry = 0;
    for a in a {
        (*a, carry) = mul_add(*a, b, carry);
    }
    carry
}

/// `r += a * b` for a single limb `b` and equal lengths; returns the carry limb.
pub fn addmul_1(r: &mut [Limb], a: &[Limb], b: Limb) -> Limb {
    assert_eq!(r.len(), a.len(), "addmul_1: length mismatch");
    let mut carry = 0;
    for (r, &a) in r.iter_mut().zip(a) {
        // a * b + r + carry <= (2^64 - 1)^2 + 2 (2^64 - 1) = 2^128 - 1.
        let t = u128::from(a) * u128::from(b) + u128::from(*r) + u128::from(carry);
        (*r, carry) = split(t);
    }
    carry
}

/// `r -= a * b` for a single limb `b` and equal lengths; returns the borrow limb
/// (the amount still to subtract above the top of `r`).
pub fn submul_1(r: &mut [Limb], a: &[Limb], b: Limb) -> Limb {
    assert_eq!(r.len(), a.len(), "submul_1: length mismatch");
    let mut carry = 0;
    for (r, &a) in r.iter_mut().zip(a) {
        let (lo, hi) = mul_add(a, b, carry);
        let borrow;
        (*r, borrow) = r.overflowing_sub(lo);
        // hi = 2^64 - 1 forces lo = 0 and no borrow, so this cannot overflow.
        carry = hi + Limb::from(borrow);
    }
    carry
}

/// `r[..n] += a * (b0 + b1 B) + c0` for `n = a.len()` and `r.len() == n + 1`,
/// with `r[n]` overwritten by the limb at that position; returns the limb
/// above.
///
/// One pass for two rows: each limb of `r` is loaded and stored once per two
/// multiplier limbs. Limb `i` takes `r[i] + a[i] b0 + c0` and passes its high
/// half on with `a[i] b1 + c1` to limb `i + 1`; both sums fit a double limb.
fn addmul_2(r: &mut [Limb], a: &[Limb], b0: Limb, b1: Limb, c0: Limb) -> Limb {
    debug_assert_eq!(r.len(), a.len() + 1);
    let (head, top) = r.split_at_mut(a.len());
    let (c0, c1) = addmul_2_carries(head, a, b0, b1, c0);
    top[0] = c0;
    c1
}

/// `r += a * (b0 + b1 B) + c0` over `n = r.len() == a.len()` limbs; returns
/// the two limbs above, `(c0, c1)` (the sum is below `B^(n + 2)`).
pub(crate) fn addmul_2_carries(
    r: &mut [Limb],
    a: &[Limb],
    b0: Limb,
    b1: Limb,
    c0: Limb,
) -> (Limb, Limb) {
    debug_assert_eq!(r.len(), a.len());
    let head = r;
    let (mut c0, mut c1) = (c0, 0);
    let step = |r: &mut Limb, a: Limb, c0: Limb, c1: Limb| {
        let x = u128::from(a) * u128::from(b0) + u128::from(*r) + u128::from(c0);
        let (lo, hi) = split(x);
        *r = lo;
        split(u128::from(a) * u128::from(b1) + u128::from(hi) + u128::from(c1))
    };
    let mut rs = head.chunks_exact_mut(2);
    let mut as_ = a.chunks_exact(2);
    for (r, a) in (&mut rs).zip(&mut as_) {
        (c0, c1) = step(&mut r[0], a[0], c0, c1);
        (c0, c1) = step(&mut r[1], a[1], c0, c1);
    }
    if let ([r], [a]) = (rs.into_remainder(), as_.remainder()) {
        (c0, c1) = step(r, *a, c0, c1);
    }
    (c0, c1)
}

/// `r = a * b` by the schoolbook method, for any lengths with
/// `r.len() == a.len() + b.len()`. Quadratic; prefer [`mul`].
///
/// The first row by [`mul_1`], then two rows at a time by [`addmul_2`]; an
/// even count leaves the last row to [`addmul_1`].
pub fn mul_basecase(r: &mut [Limb], a: &[Limb], b: &[Limb]) {
    assert_eq!(r.len(), a.len() + b.len(), "mul_basecase: length mismatch");
    let Some((&b0, b_rest)) = b.split_first() else {
        r.fill(0);
        return;
    };
    let n = a.len();
    if n == 0 {
        r.fill(0);
        return;
    }
    r[n] = mul_1(&mut r[..n], a, b0);
    let mut pairs = b_rest.chunks_exact(2);
    let mut j = 1;
    for pair in &mut pairs {
        r[n + j + 1] = addmul_2(&mut r[j..=j + n], a, pair[0], pair[1], 0);
        j += 2;
    }
    if let [bj] = pairs.remainder() {
        r[n + j] = addmul_1(&mut r[j..j + n], a, *bj);
    }
}

/// `r = a^2` by the schoolbook method with `r.len() == 2 * a.len()`: the
/// off-diagonal products are computed once and doubled. Prefer [`sqr`].
///
/// Up to 6 limbs the lengths are constants, so LLVM unrolls the rows into
/// straight-line code: the loops' setup costs a third to a half of a
/// square of 3 to 6 limbs. Longer arms grow the code faster than they gain.
pub fn sqr_basecase(r: &mut [Limb], a: &[Limb]) {
    let n = a.len();
    assert_eq!(r.len(), 2 * n, "sqr_basecase: length mismatch");
    match n {
        0 => {}
        1 => (r[0], r[1]) = mul_add(a[0], a[0], 0),
        2 => sqr_rows(&mut r[..4], &a[..2]),
        3 => sqr_rows(&mut r[..6], &a[..3]),
        4 => sqr_rows(&mut r[..8], &a[..4]),
        5 => sqr_rows(&mut r[..10], &a[..5]),
        6 => sqr_rows(&mut r[..12], &a[..6]),
        _ => sqr_rows(r, a),
    }
}

/// [`sqr_basecase`] for `n = a.len() >= 2`: the off-diagonal triangle into
/// `r[1..2n - 1]`, then [`sqr_diag_addlsh1`].
///
/// Row `i` holds `a[i] * a[i+1..]` at position `2i + 1`, and its carry the
/// limb `i + n`, which no earlier row reaches. Row 0 writes its limbs, by
/// [`mul_1`]; later rows add to them, two at a time: `a[i] a[i+1]` at
/// `2i + 1`, then `a[i]` and `a[i+1]` times `a[i+2..]` from `2i + 2` by
/// [`addmul_2`], which takes the high half of the first product as its
/// carry in. An even number of rows after the first leaves the last to
/// [`addmul_1`].
#[allow(
    clippy::inline_always,
    reason = "inlined into each arm of sqr_basecase, so that a constant length unrolls"
)]
#[inline(always)]
fn sqr_rows(r: &mut [Limb], a: &[Limb]) {
    let n = a.len();
    r[n] = mul_1(&mut r[1..n], &a[1..], a[0]);
    let rows = n - 1;
    let mut i = 1;
    while i + 1 < rows {
        let (lo, hi) = mul_add(a[i], a[i + 1], r[2 * i + 1]);
        r[2 * i + 1] = lo;
        r[i + n + 1] = addmul_2(&mut r[2 * i + 2..=i + n], &a[i + 2..], a[i], a[i + 1], hi);
        i += 2;
    }
    if i < rows {
        r[i + n] = addmul_1(&mut r[2 * i + 1..i + n], &a[i + 1..], a[i]);
    }
    sqr_diag_addlsh1(r, a);
}

/// `r = 2 t + d` in one pass, as GMP's `mpn_sqr_diag_addlsh1`, for the
/// triangle `t` in `r[1..2n - 1]` (`r[0]` and `r[2n - 1]` are ignored) and
/// the squares `d = sum a[i]^2 B^(2i)`, with `n = a.len() >= 1`.
///
/// Each limb pair `(2i - 1, 2i)` takes its doubled triangle limbs, with the
/// bit shifted out of the limb below, plus the high half of `a[i - 1]^2`
/// and the low half of `a[i]^2`. The total is `a^2 < B^(2n)`, so the top
/// limb takes the last bit and carry without overflow.
#[allow(
    clippy::inline_always,
    reason = "inlined into each arm of sqr_basecase, so that a constant length unrolls"
)]
#[inline(always)]
fn sqr_diag_addlsh1(r: &mut [Limb], a: &[Limb]) {
    let n = a.len();
    debug_assert_eq!(r.len(), 2 * n);
    let (lo, mut hi) = mul_add(a[0], a[0], 0);
    r[0] = lo;
    let (mid, top) = r[1..].split_at_mut(2 * n - 2);
    let (mut bit, mut carry) = (0, false);
    for (t, &x) in mid.chunks_exact_mut(2).zip(&a[1..]) {
        let (t0, t1) = (t[0], t[1]);
        let (lo, next) = mul_add(x, x, 0);
        let c;
        (t[0], c) = super::adc((t0 << 1) | bit, hi, carry);
        (t[1], carry) = super::adc((t1 << 1) | (t0 >> (Limb::BITS - 1)), lo, c);
        bit = t1 >> (Limb::BITS - 1);
        hi = next;
    }
    top[0] = bit + hi + Limb::from(carry);
}

/// Scratch limbs [`mul`] needs for operands of `an` and `bn` limbs.
#[must_use]
pub fn mul_scratch_len(an: usize, bn: usize) -> usize {
    let (an, bn) = if an >= bn { (an, bn) } else { (bn, an) };
    match algorithm(an, bn) {
        Algorithm::Basecase => 0,
        Algorithm::Chunked => {
            let piece = chunk_len(bn);
            let last = an % piece;
            let last = if last == 0 {
                0
            } else {
                mul_scratch_len(last, bn)
            };
            piece + bn + mul_scratch_len(piece, bn).max(last)
        }
        Algorithm::Karatsuba => karatsuba_scratch_len(an, bn, false),
        Algorithm::Toom3 => toom3_scratch_len(an, bn, false),
        Algorithm::Toom4 => toom4_scratch_len(an, bn, false),
        Algorithm::Toom32 => toom32_scratch_len(an, bn),
        Algorithm::Toom42 => toom42_scratch_len(an, bn),
        Algorithm::Ntt => ntt::mul_scratch_len(an, bn),
    }
}

/// Scratch limbs [`sqr`] needs for an operand of `n` limbs.
#[must_use]
pub fn sqr_scratch_len(n: usize) -> usize {
    match sqr_algorithm(n) {
        Algorithm::Basecase | Algorithm::Chunked | Algorithm::Toom32 | Algorithm::Toom42 => 0,
        Algorithm::Karatsuba => karatsuba_scratch_len(n, n, true),
        Algorithm::Toom3 => toom3_scratch_len(n, n, true),
        Algorithm::Toom4 => toom4_scratch_len(n, n, true),
        Algorithm::Ntt => ntt::sqr_scratch_len(n),
    }
}

/// `r = a * b` for any lengths, with `r.len() == a.len() + b.len()` and
/// `scratch.len() >= mul_scratch_len(a.len(), b.len())`.
pub fn mul(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    assert_eq!(r.len(), a.len() + b.len(), "mul: length mismatch");
    assert!(
        scratch.len() >= mul_scratch_len(a.len(), b.len()),
        "mul: scratch too short"
    );
    if a.len() >= b.len() {
        mul_unchecked(r, a, b, scratch);
    } else {
        mul_unchecked(r, b, a, scratch);
    }
}

/// `r = a^2` with `r.len() == 2 * a.len()` and
/// `scratch.len() >= sqr_scratch_len(a.len())`.
pub fn sqr(r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    assert_eq!(r.len(), 2 * a.len(), "sqr: length mismatch");
    assert!(
        scratch.len() >= sqr_scratch_len(a.len()),
        "sqr: scratch too short"
    );
    sqr_unchecked(r, a, scratch);
}

/// Whether [`mul_wrap`] modulo `B^m - 1` uses the cyclic transform: for a
/// transform length (`2^k` or `3 2^k`) from [`MUL_WRAP_NTT_THRESHOLD`] on.
pub(crate) fn wraps_cyclically(m: usize) -> bool {
    m >= MUL_WRAP_NTT_THRESHOLD && ntt::is_transform_len(m)
}

/// The shortest transform length (a power of two or three times one) at
/// least `x >= 1`: the modulus length `m >= x` for products known modulo
/// `B^m - 1`.
pub(crate) fn wrap_len(x: usize) -> usize {
    let pow2 = x.next_power_of_two();
    let three = 3 * x.div_ceil(3).next_power_of_two();
    pow2.min(three)
}

/// Scratch limbs [`mul_wrap`] needs modulo `B^m - 1` for operands of `an`
/// and `bn` limbs, both at least 1.
pub(crate) fn mul_wrap_scratch_len(m: usize, an: usize, bn: usize) -> usize {
    if an + bn <= m {
        mul_scratch_len(an, bn)
    } else if wraps_cyclically(m) {
        let folded = |len: usize| if len > m { m } else { 0 };
        folded(an) + folded(bn) + ntt::cyclic_scratch_len(m)
    } else if split::splits(m, an, bn) {
        split::mul_split_scratch_len(m, an, bn)
    } else {
        an + bn + mul_scratch_len(an, bn)
    }
}

/// `r = a * b mod (B^m - 1)` for `m = r.len() >= 1` and non-empty operands of
/// any length, as the residue in `[0, B^m - 1)`, with `scratch.len() >=
/// mul_wrap_scratch_len(m, a.len(), b.len())`.
///
/// A transform length `m` (see [`wraps_cyclically`]) takes a cyclic
/// transform of `m` coefficients where the full product would need about
/// `a.len() + b.len()`. Below it, an even `m` with long enough operands is
/// split in two halves ([`split::splits`]); any other `m` reduces the full
/// product.
pub(crate) fn mul_wrap(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let (m, an, bn) = (r.len(), a.len(), b.len());
    assert!(m >= 1 && an >= 1 && bn >= 1, "mul_wrap: empty operand");
    assert!(
        scratch.len() >= mul_wrap_scratch_len(m, an, bn),
        "mul_wrap: scratch too short"
    );
    if an + bn <= m {
        mul(&mut r[..an + bn], a, b, scratch);
        r[an + bn..].fill(0);
    } else if wraps_cyclically(m) {
        let (a, scratch) = folded(a, m, scratch);
        let (b, scratch) = folded(b, m, scratch);
        ntt::mul_cyclic(r, a, b, scratch);
    } else if split::splits(m, an, bn) {
        split::mul_split(r, a, b, scratch);
    } else {
        let (t, scratch) = scratch.split_at_mut(an + bn);
        mul(t, a, b, scratch);
        fold(r, t);
    }
    canonical_wrapped(r);
}

/// `x` itself if it has at most `m` limbs, else `x mod (B^m - 1)` in the head
/// of `scratch`; returns it and the rest of `scratch`.
fn folded<'a>(x: &'a [Limb], m: usize, scratch: &'a mut [Limb]) -> (&'a [Limb], &'a mut [Limb]) {
    if x.len() <= m {
        return (x, scratch);
    }
    let (t, rest) = scratch.split_at_mut(m);
    fold(t, x);
    (t, rest)
}

/// `r = x mod (B^m - 1)` for `m = r.len() >= 1`, as a residue in `[0, B^m)`:
/// the sum of `x`'s `m`-limb pieces, each carry out of the top added back in at
/// the bottom (`B^m = 1`).
pub(crate) fn fold(r: &mut [Limb], x: &[Limb]) {
    let m = r.len();
    let head = x.len().min(m);
    r[..head].copy_from_slice(&x[..head]);
    r[head..].fill(0);
    for piece in x[head..].chunks(m) {
        let mut carry = add_assign(r, piece);
        while carry != 0 {
            carry = add_1_assign(r, carry);
        }
    }
}

/// Maps the residue `B^m - 1` (all ones), which stands for zero, to zero.
fn canonical_wrapped(r: &mut [Limb]) {
    if r.iter().all(|&x| x == Limb::MAX) {
        r.fill(0);
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum Algorithm {
    Basecase,
    Chunked,
    Karatsuba,
    Toom3,
    Toom4,
    Toom32,
    Toom42,
    Ntt,
}

/// Whether the NTT beats Toom for `an x bn` limbs (`an >= bn`) given the
/// configured `threshold`, by size: the NTT's cost is a staircase, so a
/// product whose transform is `1 / e` times the ideal one ([`ntt::efficiency`])
/// needs its length scaled up by `e^-x` before the NTT catches up with Toom,
/// for `x = NTT_EFFICIENCY_EXPONENT / 10`. Toom-4's cost grows like `n^1.4`,
/// which alone would give `x = 1 / (1.4 - 1) = 2.5`; `cas-tune` fits it.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
fn ntt_pays(an: usize, bn: usize, threshold: usize) -> bool {
    let x = NTT_EFFICIENCY_EXPONENT as f64 / 10.0;
    bn >= threshold && bn as f64 * ntt::efficiency(an, bn).powf(x) >= threshold as f64
}

/// The algorithm [`mul`] uses for `an >= bn`.
///
/// Up to `an < 4/3 bn` the balanced algorithms take the product. Longer `a`
/// goes to the NTT if its product (one transform of the whole product, or
/// chunks of `a` by `b`'s prepared transform) costs less than Toom
/// ([`ntt_pays_unbalanced`]); else, from their thresholds, to Toom-32
/// below `an = 7/4 bn`, Toom-42 below `11/4 bn`, and beyond to pieces of
/// `2 bn` limbs, each a Toom-42 product. The boundaries are where the costs
/// of the neighbours cross, counting a product of `n` limbs as `n^1.4`.
/// Below the thresholds, Karatsuba takes `an < 2 bn` and pieces of `bn`
/// limbs the rest.
pub(crate) fn algorithm(an: usize, bn: usize) -> Algorithm {
    debug_assert!(an >= bn);
    if bn < MUL_KARATSUBA_THRESHOLD {
        Algorithm::Basecase
    } else if 3 * an >= 4 * bn {
        if ntt_pays_unbalanced(an, bn) {
            Algorithm::Ntt
        } else if bn >= MUL_TOOM32_THRESHOLD && 4 * an < 7 * bn && toom32_fits(an, bn) {
            Algorithm::Toom32
        } else if bn >= MUL_TOOM42_THRESHOLD && 4 * an < 11 * bn && toom42_fits(an, bn) {
            Algorithm::Toom42
        } else if bn > an.div_ceil(2) {
            Algorithm::Karatsuba
        } else if an > chunk_len(bn) {
            Algorithm::Chunked
        } else {
            // Only tiny shapes, which no Toom fits, get here.
            Algorithm::Basecase
        }
    } else if ntt_pays(an, bn, MUL_NTT_THRESHOLD) {
        Algorithm::Ntt
    } else if bn >= MUL_TOOM4_THRESHOLD && bn > 3 * an.div_ceil(4) {
        Algorithm::Toom4
    } else if bn >= MUL_TOOM3_THRESHOLD && bn > 2 * an.div_ceil(3) {
        Algorithm::Toom3
    } else {
        Algorithm::Karatsuba
    }
}

/// Toom's modelled cost for `an x bn` limbs (`an >= bn`), in the units of the
/// NTT's cost model ([`ntt::cost`]): scaled from the balanced crossover `T =
/// MUL_NTT_THRESHOLD`, where it equals the NTT's, as `(an / bn) (bn /
/// T)^x`, `x = NTT_UNBALANCED_EXPONENT / 10`: pieces of a `bn`-limb
/// product, each growing like Toom.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
fn toom_cost(an: usize, bn: usize) -> f64 {
    let t = MUL_NTT_THRESHOLD;
    let x = NTT_UNBALANCED_EXPONENT as f64 / 10.0;
    let at_t = ntt::cost(t, t) as f64;
    at_t * (an as f64 / bn as f64) * (bn as f64 / t as f64).powf(x)
}

/// Whether the NTT's product of an unbalanced `an x bn` (`an >= 4/3 bn`),
/// in one transform or in chunks ([`ntt::cost`]), costs less than Toom on it
/// ([`toom_cost`]).
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
fn ntt_pays_unbalanced(an: usize, bn: usize) -> bool {
    if an + bn < MUL_NTT_THRESHOLD {
        return false;
    }
    (ntt::cost(an, bn) as f64) < toom_cost(an, bn)
}

/// The modelled cost of [`mul`] for operands of `an, bn >= 1` limbs, in the
/// units of the NTT's cost model: its plan where [`algorithm`] takes the
/// NTT, else [`toom_cost`]. For choices between products of a few hundred
/// limbs and more, where Toom-4 and the NTT take them.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
pub(crate) fn cost(an: usize, bn: usize) -> f64 {
    let (an, bn) = if an >= bn { (an, bn) } else { (bn, an) };
    if algorithm(an, bn) == Algorithm::Ntt {
        ntt::cost(an, bn) as f64
    } else {
        toom_cost(an, bn)
    }
}

/// The algorithm [`sqr`] uses for `n` limbs.
pub(crate) fn sqr_algorithm(n: usize) -> Algorithm {
    if n < SQR_KARATSUBA_THRESHOLD {
        Algorithm::Basecase
    } else if n < SQR_TOOM3_THRESHOLD {
        Algorithm::Karatsuba
    } else if ntt_pays(n, n, SQR_NTT_THRESHOLD) {
        Algorithm::Ntt
    } else if n >= SQR_TOOM4_THRESHOLD && n > 3 * n.div_ceil(4) {
        Algorithm::Toom4
    } else {
        Algorithm::Toom3
    }
}

/// [`mul`] for `a.len() >= b.len()`, with lengths already checked.
fn mul_unchecked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    match algorithm(a.len(), b.len()) {
        Algorithm::Basecase => mul_basecase(r, a, b),
        Algorithm::Chunked => mul_chunked(r, a, b, scratch),
        Algorithm::Karatsuba => karatsuba(r, a, b, scratch, false),
        Algorithm::Toom3 => toom3(r, a, b, scratch, false),
        Algorithm::Toom4 => toom4(r, a, b, scratch, false),
        Algorithm::Toom32 => toom32(r, a, b, scratch),
        Algorithm::Toom42 => toom42(r, a, b, scratch),
        Algorithm::Ntt => ntt::mul(r, a, b, scratch),
    }
}

fn sqr_unchecked(r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    match sqr_algorithm(a.len()) {
        Algorithm::Basecase | Algorithm::Chunked | Algorithm::Toom32 | Algorithm::Toom42 => {
            sqr_basecase(r, a);
        }
        Algorithm::Karatsuba => karatsuba(r, a, a, scratch, true),
        Algorithm::Toom3 => toom3(r, a, a, scratch, true),
        Algorithm::Toom4 => toom4(r, a, a, scratch, true),
        Algorithm::Ntt => ntt::sqr(r, a, scratch),
    }
}

/// A recursive product: `x^2` when squaring (then `x` and `y` are the same
/// value), `x * y` otherwise; `x.len() >= y.len()`.
fn product(r: &mut [Limb], x: &[Limb], y: &[Limb], scratch: &mut [Limb], square: bool) {
    if square {
        sqr_unchecked(r, x, scratch);
    } else {
        mul_unchecked(r, x, y, scratch);
    }
}

fn product_scratch_len(xn: usize, yn: usize, square: bool) -> usize {
    if square {
        sqr_scratch_len(xn)
    } else {
        mul_scratch_len(xn, yn)
    }
}

/// The pieces of the longer operand that [`mul_chunked`] multiplies by a
/// `bn`-limb one: `2 bn` limbs, a Toom-42 shape, from its threshold; `bn`
/// limbs below.
fn chunk_len(bn: usize) -> usize {
    if bn >= MUL_TOOM42_THRESHOLD {
        2 * bn
    } else {
        bn
    }
}

/// `r = x * y` for operands in either order.
fn mul_ordered(r: &mut [Limb], x: &[Limb], y: &[Limb], scratch: &mut [Limb]) {
    if x.len() >= y.len() {
        mul_unchecked(r, x, y, scratch);
    } else {
        mul_unchecked(r, y, x, scratch);
    }
}

/// `r = a * b` for `a.len() > chunk_len(b.len())`: `a` is cut into pieces
/// of [`chunk_len`] limbs (the last one shorter), each multiplied by `b` and
/// added in place.
fn mul_chunked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let bn = b.len();
    let piece = chunk_len(bn);
    let (first, rest) = a.split_at(piece);
    mul_unchecked(&mut r[..piece + bn], first, b, scratch);
    r[piece + bn..].fill(0);
    let (tmp, scratch) = scratch.split_at_mut(piece + bn);
    for (i, chunk) in rest.chunks(piece).enumerate() {
        let tmp = &mut tmp[..chunk.len() + bn];
        mul_ordered(tmp, chunk, b, scratch);
        let carry = add_assign(&mut r[(i + 1) * piece..], tmp);
        debug_assert_eq!(carry, 0);
    }
}

/// `r[offset..] += c`, where `c` may extend past the end of `r` with zero
/// limbs. The caller guarantees the sum fits, so there is no carry out.
fn add_at(r: &mut [Limb], offset: usize, c: &[Limb]) {
    let fit = c.len().min(r.len() - offset);
    debug_assert!(c[fit..].iter().all(|&x| x == 0), "add_at: sum overflows");
    let carry = add_assign(&mut r[offset..], &c[..fit]);
    debug_assert_eq!(carry, 0, "add_at: sum overflows");
}

/// `r = |x - y|` for `x.len() >= y.len()` and `r.len() == x.len()`; returns
/// whether `x < y`.
fn abs_diff(r: &mut [Limb], x: &[Limb], y: &[Limb]) -> bool {
    let (x_lo, x_hi) = x.split_at(y.len());
    if x_hi.iter().any(|&l| l != 0) || cmp_n(x_lo, y).is_ge() {
        let borrow = super::sub(r, x, y);
        debug_assert_eq!(borrow, 0);
        false
    } else {
        // x < y, and x's high limbs are all zero.
        let (r_lo, r_hi) = r.split_at_mut(y.len());
        sub_n(r_lo, y, x_lo);
        r_hi.fill(0);
        true
    }
}

fn karatsuba_split(an: usize, bn: usize) -> (usize, usize, usize) {
    let n = an.div_ceil(2);
    (n, an - n, bn - n)
}

pub(crate) fn karatsuba_scratch_len(an: usize, bn: usize, square: bool) -> usize {
    let (n, s, t) = karatsuba_split(an, bn);
    2 * n
        + product_scratch_len(n, n, square)
            .max(product_scratch_len(s, t, square))
            .max(2 * n + 1)
}

/// Karatsuba: with `X = B^n`, `a = a0 + a1 X` and `b = b0 + b1 X`,
///
/// `a b = a0 b0 + (a0 b0 + a1 b1 - (a0 - a1)(b0 - b1)) X + a1 b1 X^2`,
///
/// three half-size products. Needs `a.len() >= b.len() > ceil(a.len() / 2)`.
pub(crate) fn karatsuba(
    r: &mut [Limb],
    a: &[Limb],
    b: &[Limb],
    scratch: &mut [Limb],
    square: bool,
) {
    let (an, bn) = (a.len(), b.len());
    let (n, s, t) = karatsuba_split(an, bn);
    debug_assert!(s >= t && t >= 1 && r.len() == an + bn);
    debug_assert!(scratch.len() >= karatsuba_scratch_len(an, bn, square));
    let (a0, a1) = a.split_at(n);
    let (b0, b1) = b.split_at(n);

    // vm1 = |a0 - a1| |b0 - b1|; the differences sit in r until v0 overwrites them.
    let (vm1, rest) = scratch.split_at_mut(2 * n);
    let negative = {
        let (da, db) = r[..2 * n].split_at_mut(n);
        let neg_a = abs_diff(da, a0, a1);
        if square {
            product(vm1, da, da, rest, true);
            false
        } else {
            let neg_b = abs_diff(db, b0, b1);
            product(vm1, da, db, rest, false);
            neg_a != neg_b
        }
    };

    let (r_lo, r_hi) = r.split_at_mut(2 * n);
    product(r_lo, a0, b0, rest, square);
    product(r_hi, a1, b1, rest, square);

    // middle = v0 + vinf - (a0 - a1)(b0 - b1) >= 0.
    let mid = &mut rest[..=2 * n];
    mid[..2 * n].copy_from_slice(r_lo);
    mid[2 * n] = 0;
    let carry = add_assign(mid, r_hi);
    debug_assert_eq!(carry, 0);
    let c = if negative {
        add_assign(mid, vm1)
    } else {
        sub_assign(mid, vm1)
    };
    debug_assert_eq!(c, 0);
    add_at(r, n, mid);
}

fn toom3_split(an: usize, bn: usize) -> (usize, usize, usize) {
    let n = an.div_ceil(3);
    (n, an - 2 * n, bn - 2 * n)
}

pub(crate) fn toom3_scratch_len(an: usize, bn: usize, square: bool) -> usize {
    let (n, s, t) = toom3_split(an, bn);
    let e = n + 1;
    product_scratch_len(n, n, square)
        .max(product_scratch_len(s, t, square))
        .max(10 * e + product_scratch_len(e, e, square))
}

/// Two's complement negation in place.
fn neg_assign(a: &mut [Limb]) {
    for x in a.iter_mut() {
        *x = !*x;
    }
    add_1_assign(a, 1);
}

fn is_negative(a: &[Limb]) -> bool {
    a.last().is_some_and(|&top| top >> (Limb::BITS - 1) == 1)
}

/// `r = |a|` for a two's complement `a`; returns whether `a < 0`.
fn abs_twos(r: &mut [Limb], a: &[Limb]) -> bool {
    r.copy_from_slice(a);
    let negative = is_negative(a);
    if negative {
        neg_assign(r);
    }
    negative
}

/// Two's complement halving of an even value.
fn half_assign(a: &mut [Limb]) {
    let top = Limb::from(is_negative(a)) << (Limb::BITS - 1);
    let out = shr_assign(a, 1);
    debug_assert_eq!(out, 0, "half_assign: odd value");
    if let Some(last) = a.last_mut() {
        *last |= top;
    }
}

/// `d^-1 mod 2^64` for an odd `d`, by Newton's iteration: `d d = 1 mod 8`
/// gives 3 correct bits and each step doubles them.
pub(crate) const fn limb_inverse(d: Limb) -> Limb {
    let mut inv = d;
    let mut i = 0;
    while i < 5 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(d.wrapping_mul(inv)));
        i += 1;
    }
    inv
}

/// `a /= d` for a multiple of the small odd `d` (`inv = limb_inverse(d)`),
/// modulo `B^len` (so two's complement values work too): multiplies by the
/// inverse of `d`, propagating a borrow.
fn divexact_odd_assign(a: &mut [Limb], d: Limb, inv: Limb) {
    debug_assert_eq!(d.wrapping_mul(inv), 1);
    let mut borrow = 0;
    for x in a {
        let (s, b) = x.overflowing_sub(borrow);
        let q = s.wrapping_mul(inv);
        *x = q;
        // d q = s + hi B, so the next limb owes hi, plus 1 if s wrapped.
        let (_, hi) = mul_add(q, d, 0);
        borrow = hi + Limb::from(b);
    }
}

fn divexact_by3_assign(a: &mut [Limb]) {
    divexact_odd_assign(a, 3, const { limb_inverse(3) });
}

/// `p = x0 + x2` and `sum = x0 + x1 + x2`, as `n + 1`-limb values; `p` is the
/// starting point for the evaluations at -1 and -2 in [`toom3`].
fn toom3_eval_plus(sum: &mut [Limb], p: &mut [Limb], x0: &[Limb], x1: &[Limb], x2: &[Limb]) {
    let n = x0.len();
    p[..n].copy_from_slice(x0);
    p[n] = 0;
    let c = add_assign(p, x2);
    debug_assert_eq!(c, 0);
    sum.copy_from_slice(p);
    let c = add_assign(sum, x1);
    debug_assert_eq!(c, 0);
}

/// Toom-3 (Toom-Cook with points `0, 1, -1, -2, ∞` and Bodrato's
/// interpolation sequence): with `X = B^n`, `a = a0 + a1 X + a2 X^2` and `b`
/// likewise, five third-size products. Needs `b.len() > 2 ceil(a.len() / 3)`.
///
/// Evaluations and the interpolation run in two's complement on fixed-width
/// buffers, so negative intermediate values need no sign bookkeeping.
pub(crate) fn toom3(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb], square: bool) {
    let (an, bn) = (a.len(), b.len());
    let (n, s, t) = toom3_split(an, bn);
    debug_assert!(s >= t && t >= 1 && r.len() == an + bn);
    debug_assert!(scratch.len() >= toom3_scratch_len(an, bn, square));
    let (a0, a12) = a.split_at(n);
    let (a1, a2) = a12.split_at(n);
    let (b0, b12) = b.split_at(n);
    let (b1, b2) = b12.split_at(n);

    // v0 = a0 b0 and vinf = a2 b2 go straight to their final place.
    product(&mut r[..2 * n], a0, b0, scratch, square);
    product(&mut r[4 * n..], a2, b2, scratch, square);
    r[2 * n..4 * n].fill(0);

    let e = n + 1;
    let w = 2 * e;
    let (v1, rest) = scratch.split_at_mut(w);
    let (vm1, rest) = rest.split_at_mut(w);
    let (vm2, rest) = rest.split_at_mut(w);
    let (pa, rest) = rest.split_at_mut(e);
    let (pb, rest) = rest.split_at_mut(e);
    let (ta, rest) = rest.split_at_mut(e);
    let (tb, rest) = rest.split_at_mut(e);

    // v1 = a(1) b(1); afterwards p = x0 + x2.
    toom3_eval_plus(ta, pa, a0, a1, a2);
    if square {
        product(v1, ta, ta, rest, true);
    } else {
        toom3_eval_plus(tb, pb, b0, b1, b2);
        product(v1, ta, tb, rest, false);
    }

    // vm1 = a(-1) b(-1), with p = x0 - x1 + x2.
    let signed_product = |out: &mut [Limb],
                          pa: &[Limb],
                          pb: &[Limb],
                          ta: &mut [Limb],
                          tb: &mut [Limb],
                          rest: &mut [Limb]| {
        let neg_a = abs_twos(ta, pa);
        if square {
            product(out, ta, ta, rest, true);
        } else {
            let neg_b = abs_twos(tb, pb);
            product(out, ta, tb, rest, false);
            if neg_a != neg_b {
                neg_assign(out);
            }
        }
    };
    sub_assign(pa, a1);
    if !square {
        sub_assign(pb, b1);
    }
    signed_product(vm1, pa, pb, ta, tb, rest);

    // vm2 = a(-2) b(-2), with x(-2) = 2 (x(-1) + x2) - x0.
    let to_minus_two = |p: &mut [Limb], x0: &[Limb], x2: &[Limb]| {
        add_assign(p, x2);
        shl_assign(p, 1);
        sub_assign(p, x0);
    };
    to_minus_two(pa, a0, a2);
    if !square {
        to_minus_two(pb, b0, b2);
    }
    signed_product(vm2, pa, pb, ta, tb, rest);

    toom3_interpolate(r, n, v1, vm1, vm2);
}

/// The interpolation of [`toom3`] and [`toom42`], whose products both have
/// degree 4: from `c0 = v(0)` in `r[..2n]`, `c4 = v(∞)` in `r[4n..]` (and
/// zeros between) and the values `v1`, `vm1`, `vm2` at `1, -1, -2` as
/// two's complement values of `w` limbs (overwritten), adds the middle
/// coefficients `c1`, `c2`, `c3` into `r` at `n`, `2n` and `3n` limbs.
fn toom3_interpolate(
    r: &mut [Limb],
    n: usize,
    v1: &mut [Limb],
    vm1: &mut [Limb],
    vm2: &mut [Limb],
) {
    // All mod B^w: afterwards v1 = c1, vm1 = c2, vm2 = c3.
    let (v0, rest_r) = r.split_at(2 * n);
    let vinf = &rest_r[2 * n..];
    sub_assign(vm2, v1); //            (v(-2) - v(1)) / 3
    divexact_by3_assign(vm2);
    sub_assign(v1, vm1); //            c1' = (v(1) - v(-1)) / 2
    half_assign(v1);
    sub_assign(vm1, v0); //            c2' = v(-1) - v0
    neg_assign(vm2); //                c3 = (c2' - vm2) / 2 + 2 vinf
    add_assign(vm2, vm1);
    half_assign(vm2);
    add_assign(vm2, vinf);
    add_assign(vm2, vinf);
    add_assign(vm1, v1); //            c2 = c2' + c1' - vinf
    sub_assign(vm1, vinf);
    sub_assign(v1, vm2); //            c1 = c1' - c3

    add_at(r, n, v1);
    add_at(r, 2 * n, vm1);
    add_at(r, 3 * n, vm2);
}

/// `x * y` into `out` for two's complement `x` and `y` (`e` limbs each),
/// through their absolute values in `tx` and `ty`.
fn signed_product(
    out: &mut [Limb],
    x: &[Limb],
    y: &[Limb],
    tx: &mut [Limb],
    ty: &mut [Limb],
    scratch: &mut [Limb],
) {
    let neg_x = abs_twos(tx, x);
    let neg_y = abs_twos(ty, y);
    mul_unchecked(out, tx, ty, scratch);
    if neg_x != neg_y {
        neg_assign(out);
    }
}

/// `(n, s, t)` for Toom-32: pieces of `n` limbs, `a`'s top piece `s` and
/// `b`'s `t` limbs.
fn toom32_split(an: usize, bn: usize) -> (usize, usize, usize) {
    let n = an.div_ceil(3).max(bn.div_ceil(2));
    (n, an.saturating_sub(2 * n), bn.saturating_sub(n))
}

/// Whether [`toom32`] takes `an x bn`: both top pieces non-empty.
pub(crate) fn toom32_fits(an: usize, bn: usize) -> bool {
    let (_, s, t) = toom32_split(an, bn);
    s >= 1 && t >= 1
}

pub(crate) fn toom32_scratch_len(an: usize, bn: usize) -> usize {
    let (n, s, t) = toom32_split(an, bn);
    let e = n + 1;
    mul_scratch_len(n, n)
        .max(mul_scratch_len(s, t))
        .max(4 * e + 4 * e + mul_scratch_len(e, e))
}

/// Toom-32 (points `0, 1, -1, ∞`): with `X = B^n`, `a = a0 + a1 X + a2 X^2`
/// and `b = b0 + b1 X`, four products of about `n` limbs for `an` about
/// `3/2 bn`. The product `c0 + c1 X + c2 X^2 + c3 X^3` has `c0 = v(0)`,
/// `c3 = v(∞)`, `c1 = (v(1) - v(-1)) / 2 - c3` and
/// `c2 = (v(1) + v(-1)) / 2 - c0`. Needs [`toom32_fits`].
pub(crate) fn toom32(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let (an, bn) = (a.len(), b.len());
    let (n, s, t) = toom32_split(an, bn);
    debug_assert!(s >= 1 && t >= 1 && r.len() == an + bn);
    debug_assert!(scratch.len() >= toom32_scratch_len(an, bn));
    let (a0, a12) = a.split_at(n);
    let (a1, a2) = a12.split_at(n);
    let (b0, b1) = b.split_at(n);

    // v0 = c0 and vinf = c3 go straight to their final place.
    mul_unchecked(&mut r[..2 * n], a0, b0, scratch);
    mul_ordered(&mut r[3 * n..], a2, b1, scratch);
    r[2 * n..3 * n].fill(0);

    let e = n + 1;
    let w = 2 * e;
    let (v1, rest) = scratch.split_at_mut(w);
    let (vm1, rest) = rest.split_at_mut(w);
    let (pa, rest) = rest.split_at_mut(e);
    let (pb, rest) = rest.split_at_mut(e);
    let (ta, rest) = rest.split_at_mut(e);
    let (tb, rest) = rest.split_at_mut(e);

    // a(1) = (a0 + a2) + a1, a(-1) = (a0 + a2) - a1; b(+-1) = b0 +- b1.
    load(pa, a0);
    add_assign(pa, a2);
    ta.copy_from_slice(pa);
    add_assign(ta, a1);
    sub_assign(pa, a1);
    load(pb, b0);
    tb.copy_from_slice(pb);
    add_assign(tb, b1);
    sub_assign(pb, b1);
    mul_unchecked(v1, ta, tb, rest);
    signed_product(vm1, pa, pb, ta, tb, rest);

    // vm1 = (v(1) - v(-1)) / 2 - c3, v1 = (v(1) + v(-1)) / 2 - c0.
    let (v0, rest_r) = r.split_at(2 * n);
    let vinf = &rest_r[n..];
    neg_assign(vm1);
    add_assign(vm1, v1);
    shl_assign(v1, 1);
    sub_assign(v1, vm1);
    half_assign(vm1);
    half_assign(v1);
    sub_assign(vm1, vinf);
    sub_assign(v1, v0);

    add_at(r, n, vm1);
    add_at(r, 2 * n, v1);
}

/// `(n, s, t)` for Toom-42: pieces of `n` limbs, `a`'s top piece `s` and
/// `b`'s `t` limbs.
fn toom42_split(an: usize, bn: usize) -> (usize, usize, usize) {
    let n = an.div_ceil(4).max(bn.div_ceil(2));
    (n, an.saturating_sub(3 * n), bn.saturating_sub(n))
}

/// Whether [`toom42`] takes `an x bn`: both top pieces non-empty.
pub(crate) fn toom42_fits(an: usize, bn: usize) -> bool {
    let (_, s, t) = toom42_split(an, bn);
    s >= 1 && t >= 1
}

pub(crate) fn toom42_scratch_len(an: usize, bn: usize) -> usize {
    let (n, s, t) = toom42_split(an, bn);
    let e = n + 1;
    mul_scratch_len(n, n)
        .max(mul_scratch_len(s, t))
        .max(3 * 2 * e + 9 * e + mul_scratch_len(e, e))
}

/// Toom-42 (points `0, 1, -1, -2, ∞`): with `X = B^n`,
/// `a = a0 + a1 X + a2 X^2 + a3 X^3` and `b = b0 + b1 X`, five products of
/// about `n` limbs for `an` about `2 bn`. The product has degree 4, so
/// [`toom3`]'s interpolation recovers it. Needs [`toom42_fits`].
pub(crate) fn toom42(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let (an, bn) = (a.len(), b.len());
    let (n, s, t) = toom42_split(an, bn);
    debug_assert!(s >= 1 && t >= 1 && r.len() == an + bn);
    debug_assert!(scratch.len() >= toom42_scratch_len(an, bn));
    let pa = quarters(a, n);
    let (b0, b1) = b.split_at(n);

    // v0 = c0 and vinf = c4 go straight to their final place.
    mul_unchecked(&mut r[..2 * n], pa[0], b0, scratch);
    mul_ordered(&mut r[4 * n..], pa[3], b1, scratch);
    r[2 * n..4 * n].fill(0);

    let e = n + 1;
    let w = 2 * e;
    let (v, rest) = scratch.split_at_mut(3 * w);
    let (ea, rest) = rest.split_at_mut(3 * e);
    let (eb, rest) = rest.split_at_mut(3 * e);
    let (ta, rest) = rest.split_at_mut(e);
    let (tb, rest) = rest.split_at_mut(e);
    let (tmp, rest) = rest.split_at_mut(e);

    // a(+-1) = (a0 + a2) +- (a1 + a3), a(-2) = (a0 + 4 a2) - 2 (a1 + 4 a3).
    let [a0, a1, a2, a3] = pa;
    let (p1, m) = ea.split_at_mut(e);
    let (m1, m2) = m.split_at_mut(e);
    load(p1, a0);
    add_assign(p1, a2);
    load(tmp, a1);
    add_assign(tmp, a3);
    m1.copy_from_slice(p1);
    add_assign(p1, tmp);
    sub_assign(m1, tmp);
    load_shifted(m2, a2, 2);
    add_assign(m2, a0);
    load_shifted(tmp, a3, 2);
    add_assign(tmp, a1);
    shl_assign(tmp, 1);
    sub_assign(m2, tmp);
    // b(+-1) = b0 +- b1, b(-2) = b0 - 2 b1.
    let (q1, m) = eb.split_at_mut(e);
    let (n1, n2) = m.split_at_mut(e);
    load(q1, b0);
    n1.copy_from_slice(q1);
    n2.copy_from_slice(q1);
    add_assign(q1, b1);
    sub_assign(n1, b1);
    load_shifted(tmp, b1, 1);
    sub_assign(n2, tmp);

    for (k, out) in v.chunks_exact_mut(w).enumerate() {
        let (xa, xb) = (&ea[k * e..(k + 1) * e], &eb[k * e..(k + 1) * e]);
        signed_product(out, xa, xb, ta, tb, rest);
    }
    let (v1, v) = v.split_at_mut(w);
    let (vm1, vm2) = v.split_at_mut(w);
    toom3_interpolate(r, n, v1, vm1, vm2);
}

fn toom4_split(an: usize, bn: usize) -> (usize, usize, usize) {
    let n = an.div_ceil(4);
    (n, an - 3 * n, bn - 3 * n)
}

pub(crate) fn toom4_scratch_len(an: usize, bn: usize, square: bool) -> usize {
    let (n, s, t) = toom4_split(an, bn);
    let e = n + 1;
    product_scratch_len(n, n, square)
        .max(product_scratch_len(s, t, square))
        .max(4 * 2 * e + 11 * e + product_scratch_len(e, e, square))
}

/// `buf = x` zero-extended (`x.len() <= buf.len()`).
fn load(buf: &mut [Limb], x: &[Limb]) {
    buf[..x.len()].copy_from_slice(x);
    buf[x.len()..].fill(0);
}

/// `buf = x << k`, zero-extended; the value fits.
fn load_shifted(buf: &mut [Limb], x: &[Limb], k: u32) {
    load(buf, x);
    let out = shl_assign(buf, k);
    debug_assert_eq!(out, 0);
}

/// The values of `x0 + x1 X + x2 X^2 + x3 X^3` at `1, -1, 2, -2` and `8 x(1/2)`
/// (`8 x0 + 4 x1 + 2 x2 + x3`) into the five `n + 1`-limb chunks of `v`, the
/// ones at `-1` and `-2` as magnitudes; returns whether those two are
/// negative. `t` is one chunk of work.
fn toom4_eval(v: &mut [Limb], t: &mut [Limb], x: [&[Limb]; 4]) -> (bool, bool) {
    let e = t.len();
    let n = e - 1;
    let [x0, x1, x2, x3] = x;
    let (p1, v) = v.split_at_mut(e);
    let (m1, v) = v.split_at_mut(e);
    let (p2, v) = v.split_at_mut(e);
    let (m2, h) = v.split_at_mut(e);
    // x(+-1) = (x0 + x2) +- (x1 + x3).
    p1[n] = add_n(&mut p1[..n], x0, x2);
    t[n] = add(&mut t[..n], x1, x3);
    let neg1 = abs_diff(m1, p1, t);
    add_n_assign(p1, t);
    // x(+-2) = (x0 + 4 x2) +- (2 x1 + 8 x3).
    fused::lsh_add(p2, x0, 0, x2, 2);
    fused::lsh_add(t, x1, 1, x3, 3);
    let neg2 = abs_diff(m2, p2, t);
    add_n_assign(p2, t);
    // 8 x(1/2) = 2 (4 x0 + x2) + (4 x1 + x3).
    fused::lsh_add(t, x2, 0, x0, 2);
    fused::lsh_add(h, x3, 0, x1, 2);
    fused::addlsh_assign(h, t, 1);
    (neg1, neg2)
}

/// `x`'s pieces of `n`, `n`, `n` and the rest limbs.
fn quarters(x: &[Limb], n: usize) -> [&[Limb]; 4] {
    [&x[..n], &x[n..2 * n], &x[2 * n..3 * n], &x[3 * n..]]
}

/// Toom-4 (Toom-Cook with points `0, 1, -1, 2, -2, 1/2, ∞`): with `X = B^n`,
/// `a = a0 + a1 X + a2 X^2 + a3 X^3` and `b` likewise, seven quarter-size
/// products. Needs `b.len() > 3 ceil(a.len() / 4)`.
///
/// The evaluations at `-1` and `-2` are magnitudes with a sign, so the
/// products are of nonnegative values and the interpolation adds or
/// subtracts them by sign ([`toom4_interpolate`]). `v(0)` and `v(∞)` land in
/// place, and `v(1)` between them, where `c2` ends up.
pub(crate) fn toom4(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb], square: bool) {
    let (an, bn) = (a.len(), b.len());
    let (n, s, t) = toom4_split(an, bn);
    debug_assert!(s >= t && t >= 1 && r.len() == an + bn);
    debug_assert!(scratch.len() >= toom4_scratch_len(an, bn, square));
    let (pa, pb) = (quarters(a, n), quarters(b, n));

    product(&mut r[..2 * n], pa[0], pb[0], scratch, square);
    product(&mut r[6 * n..], pa[3], pb[3], scratch, square);

    let e = n + 1;
    let w = 2 * e;
    let (vm1, rest) = scratch.split_at_mut(w);
    let (v2, rest) = rest.split_at_mut(w);
    let (vm2, rest) = rest.split_at_mut(w);
    let (vh, rest) = rest.split_at_mut(w);
    let (ea, rest) = rest.split_at_mut(5 * e);
    let (eb, rest) = rest.split_at_mut(5 * e);
    let (tmp, rest) = rest.split_at_mut(e);
    let (a1, a2) = toom4_eval(ea, tmp, pa);
    let (b1, b2) = if square {
        (a1, a2)
    } else {
        toom4_eval(eb, tmp, pb)
    };
    let eb: &[Limb] = if square { ea } else { eb };
    let xa = |k: usize| &ea[k * e..(k + 1) * e];
    let xb = |k: usize| &eb[k * e..(k + 1) * e];
    product(&mut r[2 * n..2 * n + w], xa(0), xb(0), rest, square);
    product(vm1, xa(1), xb(1), rest, square);
    product(v2, xa(2), xb(2), rest, square);
    product(vm2, xa(3), xb(3), rest, square);
    product(vh, xa(4), xb(4), rest, square);
    toom4_interpolate(r, n, [vm1, v2, vm2, vh], (a1 != b1, a2 != b2));
}

/// The interpolation of [`toom4`]: from `c0 = v(0)` in `r[..2n]`, `c6 =
/// v(∞)` in `r[6n..]`, `v(1)` in `r[2n..4n + 2]` and, in `v`, `|v(-1)|`,
/// `v(2)`, `|v(-2)|` and `64 v(1/2)` (all `w = 2n + 2` limbs, overwritten),
/// with `neg` telling whether `v(-1)` and `v(-2)` are negative, writes the
/// product `c0 + c1 X + ... + c6 X^6` into `r`.
///
/// Bodrato's sequence, as in GMP's `mpn_toom_interpolate_7pts`: the odd and
/// even parts at `±1` and `±2` by halving sums or differences, `c2 + 4 c4`
/// and `c2 + c4` from those, and the odd coefficients from `64 v(1/2) +
/// v(2) - 65 (c0 + c2 + c4 + c6) + 45 (c2 + c4) = 34 c1 + 16 c3 + 34 c5`.
/// Each step is one pass modulo `B^w`, where negative intermediate values
/// wrap, with the shifts, multiples and exact divisions fused into the
/// additions and subtractions next to them ([`fused`]).
fn toom4_interpolate(r: &mut [Limb], n: usize, v: [&mut [Limb]; 4], neg: (bool, bool)) {
    let [vm1, v2, vm2, vh] = v;
    let w = vm1.len();
    let rsh1 = |x: &mut [Limb], y: &[Limb], negative: bool| {
        // x = (y - v(-k)) / 2 for x = |v(-k)|.
        if negative {
            fused::rsh_add_assign(x, y, 1);
        } else {
            fused::rsh_rsub_assign(x, y, 1);
        }
    };
    {
        let (c0, rest) = r.split_at_mut(2 * n);
        let (w2, rest) = rest.split_at_mut(w);
        let (c0, c6) = (&*c0, &rest[2 * n - 2..]);
        add_n_assign(vh, v2); //                    vh = 64 v(1/2) + v(2)
        rsh1(vm2, v2, neg.1); //                    vm2 = 2 c1 + 8 c3 + 32 c5
        sub_n_assign(v2, vm2);
        sub_assign(v2, c0);
        fused::sublsh_assign(v2, c6, 6); //         v2 = 4 c2 + 16 c4
        rsh1(vm1, w2, neg.0); //                    vm1 = c1 + c3 + c5
        sub_n_assign(w2, vm1); //                   w2 = c0 + c2 + c4 + c6
        submul_1(vh, w2, 65);
        sub_assign(w2, c6);
        sub_assign(w2, c0); //                      w2 = c2 + c4
        fused::addmul_1_rsh1(vh, w2, 45); //        vh = 17 c1 + 8 c3 + 17 c5
        fused::rsh_sub_divexact_by3(v2, w2, 2); //  v2 = c4
        sub_n_assign(w2, v2); //                    w2 = c2
        fused::rsub_divexact_by15(vm2, vh); //      vm2 = c1 - c5
        fused::sublsh3_divexact_by9(vh, vm1); //    vh = c1 + c5
        sub_n_assign(vm1, vh); //                   vm1 = c3
        fused::rsh_add_assign(vm2, vh, 1); //       vm2 = c1
        sub_n_assign(vh, vm2); //                   vh = c5
    }

    // c2 sits in place, its top limb at 4n; c4 goes over the rest.
    debug_assert!(r[4 * n + 1] == 0 && v2[2 * n + 1] == 0);
    let top = r[4 * n];
    r[4 * n..6 * n].copy_from_slice(&v2[..2 * n]);
    let c = add_1_assign(&mut r[4 * n..], top);
    debug_assert_eq!(c, 0);
    let c = add_1_assign(&mut r[6 * n..], v2[2 * n]);
    debug_assert_eq!(c, 0);
    add_at(r, n, vm2);
    add_at(r, 3 * n, vm1);
    add_at(r, 5 * n, vh);
}
