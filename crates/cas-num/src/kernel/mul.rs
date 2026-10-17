//! Multiplication kernels: single-limb products, schoolbook, Karatsuba,
//! Toom-3 and number-theoretic transforms ([`super::ntt`]), plus dedicated
//! squaring.
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

use super::{add_1_assign, add_assign, cmp_n, ntt, shl_assign, shr_assign, sub_assign, sub_n};
use crate::Limb;

/// Smallest shorter-operand length at which [`mul`] uses Karatsuba instead of
/// schoolbook. These thresholds come from the tuning configuration (`build.rs`).
pub(crate) const MUL_KARATSUBA_THRESHOLD: usize = crate::tuning::MUL_KARATSUBA_THRESHOLD;
/// Smallest shorter-operand length at which [`mul`] uses Toom-3.
pub(crate) const MUL_TOOM3_THRESHOLD: usize = crate::tuning::MUL_TOOM3_THRESHOLD;
/// Smallest shorter-operand length at which [`mul`] uses the NTT. The NTT's cost is a staircase (see [`super::ntt`]), so above
/// this it is ahead on average rather than at every length.
pub(crate) const MUL_NTT_THRESHOLD: usize = crate::tuning::MUL_NTT_THRESHOLD;
/// Smallest length at which [`sqr`] uses Karatsuba squaring.
pub(crate) const SQR_KARATSUBA_THRESHOLD: usize = crate::tuning::SQR_KARATSUBA_THRESHOLD;
/// Smallest length at which [`sqr`] uses Toom-3 squaring.
pub(crate) const SQR_TOOM3_THRESHOLD: usize = crate::tuning::SQR_TOOM3_THRESHOLD;
/// Smallest length at which [`sqr`] uses the NTT.
pub(crate) const SQR_NTT_THRESHOLD: usize = crate::tuning::SQR_NTT_THRESHOLD;

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

/// `r = a * b` by the schoolbook method, for any lengths with
/// `r.len() == a.len() + b.len()`. Quadratic; prefer [`mul`].
pub fn mul_basecase(r: &mut [Limb], a: &[Limb], b: &[Limb]) {
    assert_eq!(r.len(), a.len() + b.len(), "mul_basecase: length mismatch");
    let Some((&b0, b_rest)) = b.split_first() else {
        r.fill(0);
        return;
    };
    let n = a.len();
    r[n] = mul_1(&mut r[..n], a, b0);
    for (j, &bj) in b_rest.iter().enumerate() {
        let j = j + 1;
        r[n + j] = addmul_1(&mut r[j..j + n], a, bj);
    }
}

/// `r = a^2` by the schoolbook method with `r.len() == 2 * a.len()`: the
/// off-diagonal products are computed once and doubled. Prefer [`sqr`].
pub fn sqr_basecase(r: &mut [Limb], a: &[Limb]) {
    let n = a.len();
    assert_eq!(r.len(), 2 * n, "sqr_basecase: length mismatch");
    r.fill(0);
    // Row i adds a[i] * a[i+1..] at position 2i + 1; its carry lands on the
    // still untouched limb i + n.
    for i in 0..n.saturating_sub(1) {
        r[i + n] = addmul_1(&mut r[2 * i + 1..i + n], &a[i + 1..], a[i]);
    }
    let out = shl_assign(r, 1);
    debug_assert_eq!(out, 0);
    let mut carry = false;
    for (i, &x) in a.iter().enumerate() {
        let (lo, hi) = mul_add(x, x, 0);
        let c;
        (r[2 * i], c) = super::adc(r[2 * i], lo, carry);
        (r[2 * i + 1], carry) = super::adc(r[2 * i + 1], hi, c);
    }
    debug_assert!(!carry);
}

/// Scratch limbs [`mul`] needs for operands of `an` and `bn` limbs.
#[must_use]
pub fn mul_scratch_len(an: usize, bn: usize) -> usize {
    let (an, bn) = if an >= bn { (an, bn) } else { (bn, an) };
    match algorithm(an, bn) {
        Algorithm::Basecase => 0,
        Algorithm::Chunked => {
            let last = an % bn;
            let last = if last == 0 {
                0
            } else {
                mul_scratch_len(bn, last)
            };
            2 * bn + mul_scratch_len(bn, bn).max(last)
        }
        Algorithm::Karatsuba => karatsuba_scratch_len(an, bn, false),
        Algorithm::Toom3 => toom3_scratch_len(an, bn, false),
        Algorithm::Ntt => ntt::mul_scratch_len(an, bn),
    }
}

/// Scratch limbs [`sqr`] needs for an operand of `n` limbs.
#[must_use]
pub fn sqr_scratch_len(n: usize) -> usize {
    match sqr_algorithm(n) {
        Algorithm::Basecase | Algorithm::Chunked => 0,
        Algorithm::Karatsuba => karatsuba_scratch_len(n, n, true),
        Algorithm::Toom3 => toom3_scratch_len(n, n, true),
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

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum Algorithm {
    Basecase,
    Chunked,
    Karatsuba,
    Toom3,
    Ntt,
}

/// The algorithm [`mul`] uses for `an >= bn`.
pub(crate) fn algorithm(an: usize, bn: usize) -> Algorithm {
    debug_assert!(an >= bn);
    if bn < MUL_KARATSUBA_THRESHOLD {
        Algorithm::Basecase
    } else if bn <= an.div_ceil(2) {
        // Karatsuba needs b's high half to be non-empty.
        Algorithm::Chunked
    } else if bn >= MUL_NTT_THRESHOLD {
        Algorithm::Ntt
    } else if bn >= MUL_TOOM3_THRESHOLD && bn > 2 * an.div_ceil(3) {
        Algorithm::Toom3
    } else {
        Algorithm::Karatsuba
    }
}

/// The algorithm [`sqr`] uses for `n` limbs.
pub(crate) fn sqr_algorithm(n: usize) -> Algorithm {
    if n < SQR_KARATSUBA_THRESHOLD {
        Algorithm::Basecase
    } else if n < SQR_TOOM3_THRESHOLD {
        Algorithm::Karatsuba
    } else if n < SQR_NTT_THRESHOLD {
        Algorithm::Toom3
    } else {
        Algorithm::Ntt
    }
}

/// [`mul`] for `a.len() >= b.len()`, with lengths already checked.
fn mul_unchecked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    match algorithm(a.len(), b.len()) {
        Algorithm::Basecase => mul_basecase(r, a, b),
        Algorithm::Chunked => mul_chunked(r, a, b, scratch),
        Algorithm::Karatsuba => karatsuba(r, a, b, scratch, false),
        Algorithm::Toom3 => toom3(r, a, b, scratch, false),
        Algorithm::Ntt => ntt::mul(r, a, b, scratch),
    }
}

fn sqr_unchecked(r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    match sqr_algorithm(a.len()) {
        Algorithm::Basecase | Algorithm::Chunked => sqr_basecase(r, a),
        Algorithm::Karatsuba => karatsuba(r, a, a, scratch, true),
        Algorithm::Toom3 => toom3(r, a, a, scratch, true),
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

/// `r = a * b` for `b.len() <= a.len() / 2` (rounded up): `a` is cut into
/// chunks of `b.len()` limbs, each multiplied by `b` and added in place.
fn mul_chunked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let bn = b.len();
    let (first, rest) = a.split_at(bn);
    mul_unchecked(&mut r[..2 * bn], first, b, scratch);
    r[2 * bn..].fill(0);
    let (tmp, scratch) = scratch.split_at_mut(2 * bn);
    for (i, chunk) in rest.chunks(bn).enumerate() {
        let tmp = &mut tmp[..chunk.len() + bn];
        mul_unchecked(tmp, b, chunk, scratch);
        let carry = add_assign(&mut r[(i + 1) * bn..], tmp);
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

/// `a /= 3` for a multiple of 3, modulo `B^len` (so two's complement values
/// work too): multiplies by the inverse of 3, propagating a borrow.
fn divexact_by3_assign(a: &mut [Limb]) {
    const INV3: Limb = 0xAAAA_AAAA_AAAA_AAAB; // 3 * INV3 = 1 mod 2^64
    let mut borrow = 0;
    for x in a {
        let (s, b) = x.overflowing_sub(borrow);
        let q = s.wrapping_mul(INV3);
        *x = q;
        // 3 q = s + hi B, so the next limb owes hi, plus 1 if s wrapped.
        let (_, hi) = mul_add(q, 3, 0);
        borrow = hi + Limb::from(b);
    }
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

    // Interpolation (all mod B^w): afterwards v1 = c1, vm1 = c2, vm2 = c3.
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
