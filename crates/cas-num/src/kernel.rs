//! Low-level limb kernels on slices, in the style of GMP's `mpn` layer.
//!
//! Numbers are little-endian slices of [`Limb`]s. Kernels never allocate and
//! never normalise: lengths are exactly what the caller passes, and carries or
//! borrows out of the top limb are returned. Every algorithm in the crate is
//! built from these, so they are tested once and reused.
//!
//! Naming: `op_n` takes operands of equal length; `op_1` takes a single-limb
//! second operand; plain `op` takes a longer first operand; the `_assign`
//! variants write into their first operand.
//!
//! Length preconditions are checked with `assert!`: a mismatch is a caller bug
//! that would otherwise yield a silently wrong result. These assertions are the
//! only panics, and are documented here once rather than on every function.
#![allow(clippy::missing_panics_doc)]

use std::cmp::Ordering;

use crate::Limb;

pub(crate) mod mul;

pub(crate) mod div;

pub use div::{div_rem, div_rem_1, div_rem_scratch_len, divexact_1, rem_1};
pub use mul::{
    addmul_1, mul, mul_1, mul_1_assign, mul_basecase, mul_scratch_len, sqr, sqr_basecase,
    sqr_scratch_len, submul_1,
};

/// `a + b + carry`, returning `(sum, carry_out)`; `carry` must be 0 or 1.
#[inline]
const fn adc(a: Limb, b: Limb, carry: Limb) -> (Limb, Limb) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, (c1 | c2) as Limb)
}

/// `a - b - borrow`, returning `(difference, borrow_out)`; `borrow` must be 0 or 1.
#[inline]
const fn sbb(a: Limb, b: Limb, borrow: Limb) -> (Limb, Limb) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, (b1 | b2) as Limb)
}

/// `r = a + b` for equal lengths; returns the carry (0 or 1).
pub fn add_n(r: &mut [Limb], a: &[Limb], b: &[Limb]) -> Limb {
    assert!(
        r.len() == a.len() && a.len() == b.len(),
        "add_n: length mismatch"
    );
    let mut carry = 0;
    for ((r, &a), &b) in r.iter_mut().zip(a).zip(b) {
        (*r, carry) = adc(a, b, carry);
    }
    carry
}

/// `a += b` for equal lengths; returns the carry (0 or 1).
pub fn add_n_assign(a: &mut [Limb], b: &[Limb]) -> Limb {
    assert_eq!(a.len(), b.len(), "add_n_assign: length mismatch");
    let mut carry = 0;
    for (a, &b) in a.iter_mut().zip(b) {
        (*a, carry) = adc(*a, b, carry);
    }
    carry
}

/// `r = a + b` for a single limb `b`; returns the carry (0 or 1).
///
/// With `a` empty, returns `b` (the whole sum is carried out).
pub fn add_1(r: &mut [Limb], a: &[Limb], b: Limb) -> Limb {
    assert_eq!(r.len(), a.len(), "add_1: length mismatch");
    let mut carry = b;
    for (r, &a) in r.iter_mut().zip(a) {
        let c;
        (*r, c) = a.overflowing_add(carry);
        carry = Limb::from(c);
    }
    carry
}

/// `a += b` for a single limb `b`; returns the carry (0 or 1).
///
/// Stops as soon as the carry dies out, so this is O(1) on average.
pub fn add_1_assign(a: &mut [Limb], b: Limb) -> Limb {
    let mut carry = b;
    for a in a {
        if carry == 0 {
            break;
        }
        let c;
        (*a, c) = a.overflowing_add(carry);
        carry = Limb::from(c);
    }
    carry
}

/// `r = a + b` with `a.len() >= b.len()` and `r.len() == a.len()`; returns the carry.
pub fn add(r: &mut [Limb], a: &[Limb], b: &[Limb]) -> Limb {
    assert!(
        a.len() >= b.len() && r.len() == a.len(),
        "add: length mismatch"
    );
    let (r_lo, r_hi) = r.split_at_mut(b.len());
    let (a_lo, a_hi) = a.split_at(b.len());
    let carry = add_n(r_lo, a_lo, b);
    add_1(r_hi, a_hi, carry)
}

/// `a += b` with `a.len() >= b.len()`; returns the carry.
pub fn add_assign(a: &mut [Limb], b: &[Limb]) -> Limb {
    assert!(a.len() >= b.len(), "add_assign: length mismatch");
    let (lo, hi) = a.split_at_mut(b.len());
    let carry = add_n_assign(lo, b);
    add_1_assign(hi, carry)
}

/// `r = a - b` for equal lengths; returns the borrow (0 or 1).
///
/// A borrow of 1 means `a < b` and `r` holds `a - b + 2^(64 n)`.
pub fn sub_n(r: &mut [Limb], a: &[Limb], b: &[Limb]) -> Limb {
    assert!(
        r.len() == a.len() && a.len() == b.len(),
        "sub_n: length mismatch"
    );
    let mut borrow = 0;
    for ((r, &a), &b) in r.iter_mut().zip(a).zip(b) {
        (*r, borrow) = sbb(a, b, borrow);
    }
    borrow
}

/// `a -= b` for equal lengths; returns the borrow (0 or 1).
pub fn sub_n_assign(a: &mut [Limb], b: &[Limb]) -> Limb {
    assert_eq!(a.len(), b.len(), "sub_n_assign: length mismatch");
    let mut borrow = 0;
    for (a, &b) in a.iter_mut().zip(b) {
        (*a, borrow) = sbb(*a, b, borrow);
    }
    borrow
}

/// `r = a - b` for a single limb `b`; returns the borrow (0 or 1).
///
/// With `a` empty, returns 1 if `b != 0`.
pub fn sub_1(r: &mut [Limb], a: &[Limb], b: Limb) -> Limb {
    assert_eq!(r.len(), a.len(), "sub_1: length mismatch");
    let mut borrow = b;
    for (r, &a) in r.iter_mut().zip(a) {
        let c;
        (*r, c) = a.overflowing_sub(borrow);
        borrow = Limb::from(c);
    }
    Limb::from(borrow != 0)
}

/// `a -= b` for a single limb `b`; returns the borrow (0 or 1).
///
/// Stops as soon as the borrow dies out.
pub fn sub_1_assign(a: &mut [Limb], b: Limb) -> Limb {
    let mut borrow = b;
    for a in a {
        if borrow == 0 {
            break;
        }
        let c;
        (*a, c) = a.overflowing_sub(borrow);
        borrow = Limb::from(c);
    }
    Limb::from(borrow != 0)
}

/// `r = a - b` with `a.len() >= b.len()` and `r.len() == a.len()`; returns the borrow.
pub fn sub(r: &mut [Limb], a: &[Limb], b: &[Limb]) -> Limb {
    assert!(
        a.len() >= b.len() && r.len() == a.len(),
        "sub: length mismatch"
    );
    let (r_lo, r_hi) = r.split_at_mut(b.len());
    let (a_lo, a_hi) = a.split_at(b.len());
    let borrow = sub_n(r_lo, a_lo, b);
    sub_1(r_hi, a_hi, borrow)
}

/// `a -= b` with `a.len() >= b.len()`; returns the borrow.
pub fn sub_assign(a: &mut [Limb], b: &[Limb]) -> Limb {
    assert!(a.len() >= b.len(), "sub_assign: length mismatch");
    let (lo, hi) = a.split_at_mut(b.len());
    let borrow = sub_n_assign(lo, b);
    sub_1_assign(hi, borrow)
}

/// Compares equal-length slices as numbers.
#[must_use]
pub fn cmp_n(a: &[Limb], b: &[Limb]) -> Ordering {
    assert_eq!(a.len(), b.len(), "cmp_n: length mismatch");
    a.iter().rev().cmp(b.iter().rev())
}

/// Compares normalised slices (no zero top limb) of any lengths as numbers.
#[must_use]
pub fn cmp(a: &[Limb], b: &[Limb]) -> Ordering {
    debug_assert!(is_normalized(a) && is_normalized(b));
    a.len().cmp(&b.len()).then_with(|| cmp_n(a, b))
}

/// Length of `a` without its zero high limbs.
#[must_use]
pub fn normalized_len(a: &[Limb]) -> usize {
    a.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1)
}

/// Whether `a` has no zero top limb (the empty slice is normalised zero).
#[must_use]
pub fn is_normalized(a: &[Limb]) -> bool {
    a.last() != Some(&0)
}

/// `r = a << shift` for `shift < 64`; returns the bits shifted out of the top
/// limb, in the low `shift` bits of the result.
pub fn shl(r: &mut [Limb], a: &[Limb], shift: u32) -> Limb {
    assert_eq!(r.len(), a.len(), "shl: length mismatch");
    assert!(shift < Limb::BITS, "shl: shift out of range");
    if shift == 0 {
        r.copy_from_slice(a);
        return 0;
    }
    let back = Limb::BITS - shift;
    let mut prev = 0;
    for (r, &a) in r.iter_mut().zip(a) {
        *r = (a << shift) | prev;
        prev = a >> back;
    }
    prev
}

/// `a <<= shift` for `shift < 64`; returns the bits shifted out (see [`shl`]).
pub fn shl_assign(a: &mut [Limb], shift: u32) -> Limb {
    assert!(shift < Limb::BITS, "shl_assign: shift out of range");
    if shift == 0 {
        return 0;
    }
    let back = Limb::BITS - shift;
    let mut prev = 0;
    for a in a {
        let x = *a;
        *a = (x << shift) | prev;
        prev = x >> back;
    }
    prev
}

/// `r = a >> shift` for `shift < 64`; returns the bits shifted out of the bottom
/// limb, in the high `shift` bits of the result.
pub fn shr(r: &mut [Limb], a: &[Limb], shift: u32) -> Limb {
    assert_eq!(r.len(), a.len(), "shr: length mismatch");
    assert!(shift < Limb::BITS, "shr: shift out of range");
    if shift == 0 {
        r.copy_from_slice(a);
        return 0;
    }
    let back = Limb::BITS - shift;
    let mut prev = 0;
    for (r, &a) in r.iter_mut().zip(a).rev() {
        *r = (a >> shift) | prev;
        prev = a << back;
    }
    prev
}

/// `a >>= shift` for `shift < 64`; returns the bits shifted out (see [`shr`]).
pub fn shr_assign(a: &mut [Limb], shift: u32) -> Limb {
    assert!(shift < Limb::BITS, "shr_assign: shift out of range");
    if shift == 0 {
        return 0;
    }
    let back = Limb::BITS - shift;
    let mut prev = 0;
    for a in a.iter_mut().rev() {
        let x = *a;
        *a = (x >> shift) | prev;
        prev = x << back;
    }
    prev
}

macro_rules! bitwise_n {
    ($($name:ident, $assign:ident, $op:tt, $op_assign:tt, $doc:literal;)*) => {$(
        #[doc = concat!("`r = a ", $doc, " b` limb-wise, for equal lengths.")]
        pub fn $name(r: &mut [Limb], a: &[Limb], b: &[Limb]) {
            assert!(r.len() == a.len() && a.len() == b.len(), concat!(stringify!($name), ": length mismatch"));
            for ((r, &a), &b) in r.iter_mut().zip(a).zip(b) {
                *r = a $op b;
            }
        }

        #[doc = concat!("`a ", $doc, "= b` limb-wise, for equal lengths.")]
        pub fn $assign(a: &mut [Limb], b: &[Limb]) {
            assert_eq!(a.len(), b.len(), concat!(stringify!($assign), ": length mismatch"));
            for (a, &b) in a.iter_mut().zip(b) {
                *a $op_assign b;
            }
        }
    )*};
}

bitwise_n! {
    and_n, and_n_assign, &, &=, "&";
    or_n, or_n_assign, |, |=, "|";
    xor_n, xor_n_assign, ^, ^=, "^";
}

/// Number of one bits.
#[must_use]
pub fn count_ones(a: &[Limb]) -> u64 {
    a.iter().map(|x| u64::from(x.count_ones())).sum()
}

/// Number of trailing zero bits, or `None` if every limb is zero.
#[must_use]
pub fn trailing_zeros(a: &[Limb]) -> Option<u64> {
    let i = a.iter().position(|&x| x != 0)?;
    Some(limbs_to_bits(i) + u64::from(a[i].trailing_zeros()))
}

/// Number of significant bits of a normalised slice (0 for zero).
#[must_use]
pub fn bit_length(a: &[Limb]) -> u64 {
    debug_assert!(is_normalized(a));
    match a.last() {
        None => 0,
        Some(top) => limbs_to_bits(a.len()) - u64::from(top.leading_zeros()),
    }
}

/// `n * 64` as a bit count. Slice lengths are bounded by `isize::MAX` bytes, so
/// this cannot overflow.
pub(crate) fn limbs_to_bits(n: usize) -> u64 {
    u64::try_from(n).expect("usize fits u64") * u64::from(Limb::BITS)
}

#[cfg(test)]
mod tests;
