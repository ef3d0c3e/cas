//! `Float`: the midpoint type of balls (CAS-DOC-5 §3).
//!
//! A value is `±0.m × 2^exp`: the mantissa `m` is a fraction in `[1/2, 1)`,
//! stored top-aligned in limbs (the top limb's top bit set, the lowest limb
//! non-zero), and the exponent is an [`Integer`], so it never overflows. Zero
//! is the empty mantissa with exponent 0 and no sign. Each value has exactly
//! one representation, so `Eq` and `Hash` are structural and agree with the
//! value. Only finite values exist (§3.3).
//!
//! Exact values stay short: `3` is one limb, whatever precision it is later
//! rounded to, and up to 4 limbs (256 bits) the mantissa is stored inline.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

use crate::kernel;
use crate::natural::limb_vec::LimbVec;
use crate::{Integer, Limb, Natural, Rational};

mod arith;
mod convert;
mod divsqrt;
mod round;

#[cfg(test)]
mod tests;

pub use divsqrt::APPROX_ULPS;
pub use round::Round;

/// Bits per limb, as the `u64` that bit counts use.
const LIMB_BITS: u64 = Limb::BITS as u64;

/// Limbs of the stack buffers [`with_limbs`] hands out before falling back to
/// the heap.
const STACK_LIMBS: usize = 16;

/// Runs `f` on `n` zeroed limbs: on the stack up to [`STACK_LIMBS`], else on
/// the heap, so operations on small floats never allocate a temporary.
fn with_limbs<R>(n: usize, f: impl FnOnce(&mut [Limb]) -> R) -> R {
    if n <= STACK_LIMBS {
        let mut buf = [0; STACK_LIMBS];
        f(&mut buf[..n])
    } else {
        f(&mut vec![0; n])
    }
}

/// An arbitrary-precision binary floating-point number.
///
/// Operations come in two forms: **exact** ([`Float::add`], [`Float::mul`],
/// negation, [`Float::mul_2exp`], comparisons), whose result grows as needed,
/// and **rounded** ([`Float::add_round`] and the others), which return the
/// value correctly rounded to a precision in bits and a [`Round`] mode, with
/// the sign of `rounded - exact` (MPFR's ternary value) as an [`Ordering`].
/// Rounded results equal MPFR's for the same precision and mode.
#[derive(Clone)]
pub struct Float {
    mant: LimbVec,
    exp: Integer,
    neg: bool,
}

impl Float {
    /// Zero.
    pub const ZERO: Float = Float {
        mant: LimbVec::new(),
        exp: Integer::ZERO,
        neg: false,
    };

    /// `±limbs × 2^e`, exactly, for any little-endian limbs (zero limbs at
    /// either end allowed).
    pub(crate) fn from_raw(neg: bool, limbs: &[Limb], e: Integer) -> Float {
        let n = kernel::normalized_len(limbs);
        let Some(low) = limbs[..n].iter().position(|&x| x != 0) else {
            return Float::ZERO;
        };
        let lz = limbs[n - 1].leading_zeros();
        let mut mant = LimbVec::from_slice(&limbs[low..n]);
        let out = kernel::shl_assign(mant.as_mut_slice(), lz);
        debug_assert_eq!(out, 0);
        // The shift can clear the lowest limb (its bits all moved up).
        if mant.as_slice()[0] == 0 {
            mant = LimbVec::from_slice(&mant.as_slice()[1..]);
        }
        // The value has 64 n - lz bits above 2^e.
        let exp = e + Integer::from(LIMB_BITS * n as u64 - u64::from(lz));
        Float::from_parts(neg, mant, exp)
    }

    /// `±buf × 2^e` for any limbs (zero limbs at either end allowed), exact or
    /// rounded, with set bits below `buf` if `sticky` (then `buf` must keep
    /// `prec + 1` bits from its top set bit down). `buf` is used as scratch.
    pub(crate) fn build(
        neg: bool,
        buf: &mut [Limb],
        e: &Integer,
        sticky: bool,
        rnd: Option<(u64, Round)>,
    ) -> (Float, Ordering) {
        let n = kernel::normalized_len(buf);
        if n == 0 {
            debug_assert!(!sticky);
            return (Float::ZERO, Ordering::Equal);
        }
        let m = &mut buf[..n];
        let lz = m[n - 1].leading_zeros();
        let out = kernel::shl_assign(m, lz);
        debug_assert_eq!(out, 0);
        let exp = e + &Integer::from(LIMB_BITS * n as u64 - u64::from(lz));
        match rnd {
            Some((prec, mode)) => Float::finish_round(neg, m, sticky, exp, prec, mode),
            None => (Float::from_aligned(neg, m, exp), Ordering::Equal),
        }
    }

    /// `±m × 2^exp` for a top-aligned `m` (top bit set) whose low limbs may be
    /// zero; they are dropped.
    pub(crate) fn from_aligned(neg: bool, m: &[Limb], exp: Integer) -> Float {
        let low = m.iter().position(|&x| x != 0).unwrap_or(m.len());
        if low == m.len() {
            return Float::ZERO;
        }
        Float::from_parts(neg, LimbVec::from_slice(&m[low..]), exp)
    }

    /// `±0.m × 2^e` for a one-limb mantissa with its top bit set.
    pub(crate) fn from_limb(neg: bool, m: Limb, e: i64) -> Float {
        Float::from_parts(neg, LimbVec::from_pair(m, 0), Integer::from(e))
    }

    /// Wraps an already-normalised mantissa.
    fn from_parts(neg: bool, mant: LimbVec, exp: Integer) -> Float {
        let x = Float { mant, exp, neg };
        x.debug_check();
        x
    }

    /// `m × 2^e`, exactly.
    #[must_use]
    pub fn from_integer_exp(m: &Integer, e: &Integer) -> Float {
        let (neg, mag) = m.clone().into_parts();
        Float::from_raw(neg, mag.limbs(), e.clone())
    }

    /// `(m, e)` with `self = m × 2^e`, where `m` is the stored mantissa as an
    /// integer (its trailing zero bits kept) and `e = exp - 64 × limbs`. Zero
    /// gives `(0, 0)`.
    #[must_use]
    pub fn to_integer_exp(&self) -> (Integer, Integer) {
        if self.is_zero() {
            return (Integer::ZERO, Integer::ZERO);
        }
        let m = Integer::from_parts(self.neg, Natural::from_limbs(self.mant()));
        (m, &self.exp - &Integer::from(self.mant_bits_stored()))
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.mant.len() == 0
    }

    #[must_use]
    pub fn is_negative(&self) -> bool {
        self.neg
    }

    /// `Less`, `Equal` or `Greater` as the value is negative, zero or positive.
    #[must_use]
    pub fn sign(&self) -> Ordering {
        if self.is_zero() {
            Ordering::Equal
        } else if self.neg {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }

    /// The exponent `e` with `2^(e - 1) <= |self| < 2^e`; 0 for zero.
    #[must_use]
    pub fn exponent(&self) -> &Integer {
        &self.exp
    }

    /// The number of significant bits of the mantissa: the precision at which
    /// the value is exact (0 for zero).
    #[must_use]
    pub fn bit_length(&self) -> u64 {
        match self.mant().first() {
            None => 0,
            Some(&low) => self.mant_bits_stored() - u64::from(low.trailing_zeros()),
        }
    }

    /// The mantissa limbs, top-aligned (empty for zero).
    pub(crate) fn mant(&self) -> &[Limb] {
        self.mant.as_slice()
    }

    /// `64 ×` the number of mantissa limbs.
    fn mant_bits_stored(&self) -> u64 {
        LIMB_BITS * self.mant.len() as u64
    }

    /// `-self`.
    #[must_use]
    pub fn neg(&self) -> Float {
        let mut x = self.clone();
        x.neg = !x.neg && !x.is_zero();
        x
    }

    /// `|self|`.
    #[must_use]
    pub fn abs(&self) -> Float {
        let mut x = self.clone();
        x.neg = false;
        x
    }

    /// `self × 2^k`, exactly.
    #[must_use]
    pub fn mul_2exp(&self, k: &Integer) -> Float {
        if self.is_zero() {
            return Float::ZERO;
        }
        let mut x = self.clone();
        x.exp += k;
        x
    }

    /// Compares `|self|` and `|other|`.
    #[must_use]
    pub fn cmp_abs(&self, other: &Float) -> Ordering {
        match (self.is_zero(), other.is_zero()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }
        self.exp.cmp(&other.exp).then_with(|| {
            // Top-aligned: compare from the top, the shorter one padded with
            // zero limbs.
            let (a, b) = (self.mant(), other.mant());
            for (x, y) in a.iter().rev().zip(b.iter().rev()) {
                match x.cmp(y) {
                    Ordering::Equal => {}
                    o => return o,
                }
            }
            a.len().cmp(&b.len())
        })
    }

    #[inline]
    fn debug_check(&self) {
        debug_assert!(self.is_canonical(), "non-canonical Float: {self:?}");
    }

    /// Whether the representation invariant holds (for tests).
    pub(crate) fn is_canonical(&self) -> bool {
        let m = self.mant();
        match (m.first(), m.last()) {
            (None, _) => self.exp.is_zero() && !self.neg,
            (Some(&low), Some(&top)) => low != 0 && top >> (Limb::BITS - 1) == 1,
            _ => unreachable!(),
        }
    }
}

impl Default for Float {
    fn default() -> Self {
        Float::ZERO
    }
}

impl PartialEq for Float {
    fn eq(&self, other: &Float) -> bool {
        self.neg == other.neg && self.exp == other.exp && self.mant() == other.mant()
    }
}

impl Eq for Float {}

impl Hash for Float {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.neg.hash(state);
        self.exp.hash(state);
        self.mant().hash(state);
    }
}

impl Ord for Float {
    fn cmp(&self, other: &Float) -> Ordering {
        match self.sign().cmp(&other.sign()) {
            Ordering::Equal => {}
            o => return o,
        }
        let abs = self.cmp_abs(other);
        if self.neg { abs.reverse() } else { abs }
    }
}

impl PartialOrd for Float {
    fn partial_cmp(&self, other: &Float) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl std::ops::Neg for Float {
    type Output = Float;

    fn neg(mut self) -> Float {
        self.neg = !self.neg && !self.is_zero();
        self
    }
}

impl std::ops::Neg for &Float {
    type Output = Float;

    fn neg(self) -> Float {
        Float::neg(self)
    }
}

/// `-0x0.c000000000000000p+2`: the sign, the mantissa limbs in hex from the
/// top (limb boundaries marked by `_`), and the exponent in decimal.
impl fmt::Debug for Float {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        if self.neg {
            f.write_str("-")?;
        }
        f.write_str("0x0.")?;
        for (i, limb) in self.mant().iter().rev().enumerate() {
            if i > 0 {
                f.write_str("_")?;
            }
            write!(f, "{limb:016x}")?;
        }
        write!(f, "p{}", self.exp)
    }
}

impl From<&Integer> for Float {
    fn from(n: &Integer) -> Float {
        Float::from_integer_exp(n, &Integer::ZERO)
    }
}

impl From<Integer> for Float {
    fn from(n: Integer) -> Float {
        Float::from(&n)
    }
}

impl From<i64> for Float {
    fn from(n: i64) -> Float {
        Float::from(&Integer::from(n))
    }
}

impl From<&Float> for Rational {
    /// The exact value.
    ///
    /// # Panics
    ///
    /// If the exponent is so large that the integers don't fit in memory.
    fn from(x: &Float) -> Rational {
        x.to_rational()
    }
}

impl Float {
    /// The exact value as a fraction.
    ///
    /// # Panics
    ///
    /// If the exponent is so large that the integers don't fit in memory.
    #[must_use]
    pub fn to_rational(&self) -> Rational {
        let (m, e) = self.to_integer_exp();
        let shift = |e: &Integer| -> u64 {
            u64::try_from(&e.unsigned_abs()).expect("Float::to_rational: exponent too large")
        };
        if e.is_negative() {
            let den = Integer::from(Natural::power_of_two(shift(&e)));
            Rational::new(m, den)
        } else {
            let (neg, mag) = m.into_parts();
            Rational::from(Integer::from_parts(neg, &mag << shift(&e)))
        }
    }
}
