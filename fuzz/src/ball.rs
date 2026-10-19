//! Inputs and oracles of the ball function targets (CAS-113).
//!
//! - Floats next to the points where functions change behaviour (0, ±1,
//!   small integers and half-integers, multiples of π/2 rounded to a few
//!   bits), and generic ones from tiny to huge, some beyond MPFR's
//!   exponent range.
//! - Radii from 0 through about an ulp of the midpoint to enormous.
//! - Precisions across the lazy tables' and bit-burst thresholds.
//! - MPFR brackets, compared with a ball exactly ([`meets`]): a ball holding
//!   the true value must meet the bracket that holds it too.

use std::sync::Once;

use arbitrary::{Result, Unstructured};
use cas_num::ball::{Ball, Float, Mag, Round as OurRound};
use cas_num::{Integer, Natural};
use rug::Float as MpFloat;
use rug::float::{Constant, Round};

use crate::{gmp_int, limb};

/// How far from 1 a float's magnitude may go.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Range {
    /// Exponents within `±2^20`, and rarely far beyond MPFR's range.
    Wide,
    /// Magnitudes below `2^20` (and exact 0 or above `2^-300`), where MPC's
    /// functions neither overflow nor underflow.
    Narrow,
}

/// A precision: mostly up to 256 bits, sometimes up to the lazy tables'
/// limit (4600 bits), rarely past the bit-burst threshold (`10^4`).
pub fn precision(u: &mut Unstructured<'_>) -> Result<u64> {
    Ok(match u.int_in_range(0u8..=15)? {
        0 => u.int_in_range(4_097..=12_000)?,
        1 | 2 => u.int_in_range(257..=4_096)?,
        _ => u.int_in_range(1..=256)?,
    })
}

/// A non-zero mantissa of one limb, sometimes two to `max`.
fn mantissa(u: &mut Unstructured<'_>, max: usize) -> Result<Natural> {
    let n = match u.int_in_range(0u8..=7)? {
        0 => u.int_in_range(2..=max)?,
        1 => 2,
        _ => 1,
    };
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(limb(u)?);
    }
    let m = Natural::from_limbs(&v);
    Ok(if m.is_zero() { Natural::from(1u64) } else { m })
}

/// The binary exponent `t` of a value (`2^(t-1) <= |x| < 2^t`).
fn top(u: &mut Unstructured<'_>, range: Range) -> Result<Integer> {
    let t = match (range, u.int_in_range(0u8..=31)?) {
        (Range::Wide, 0) => {
            // Beyond MPFR's exponents (about ±2^30), next to i64's edges,
            // where the fixed-size paths hand over to the general ones, and
            // past them.
            let k = u.int_in_range(31u64..=70)?;
            let d = Integer::from(u.int_in_range(-700i64..=700)?);
            return Ok(Integer::from_parts(u.arbitrary()?, Natural::from(1u64) << k) + d);
        }
        (Range::Wide, 1 | 2) => u.int_in_range(-(1i64 << 20)..=1 << 20)?,
        (Range::Wide, 3..=8) => u.int_in_range(-300..=64)?,
        (Range::Narrow, 0..=8) => u.int_in_range(-300..=20)?,
        _ => u.int_in_range(-8..=8)?,
    };
    Ok(Integer::from(t))
}

/// A float next to a point where functions change behaviour: 0, a small
/// integer or half-integer, or a multiple of π/2 rounded to a few bits
/// (the hardest trigonometric reductions), as is or moved by a few units of
/// `2^-k`.
fn special(u: &mut Unstructured<'_>, range: Range) -> Result<Float> {
    let base = match u.int_in_range(0u8..=3)? {
        0 => Float::ZERO,
        1 => Float::from(u.int_in_range(-8i64..=8)?),
        2 => Float::from(2 * u.int_in_range(-64i64..=63)? + 1).mul_2exp(&Integer::from(-1i64)),
        _ => {
            let k = match (range, u.int_in_range(0u8..=3)?) {
                (Range::Wide, 0) => u.int_in_range(-(1i64 << 40)..=1 << 40)?,
                (_, 1) => u.int_in_range(-(1i64 << 18)..=1 << 18)?,
                _ => u.int_in_range(-12i64..=12)?,
            };
            let bits = 64 - k.unsigned_abs().leading_zeros() + u.int_in_range(2u32..=300)?;
            let x = MpFloat::with_val(bits, MpFloat::with_val(bits, Constant::Pi) * k) >> 1u32;
            from_mp(&x).expect("finite")
        }
    };
    if u.ratio(1u8, 3)? {
        return Ok(base);
    }
    let d = Float::from(u.int_in_range(1i64..=4)? * if u.arbitrary()? { -1 } else { 1 });
    let k = u.int_in_range(1i64..=400)?;
    Ok(base.add(&d.mul_2exp(&Integer::from(-k))))
}

/// A float: zero, a [`special`] point, or a mantissa of up to 8 limbs at an
/// exponent from `range`.
pub fn float(u: &mut Unstructured<'_>, range: Range) -> Result<Float> {
    match u.int_in_range(0u8..=7)? {
        0 => Ok(Float::ZERO),
        1 | 2 => special(u, range),
        _ => {
            let m = mantissa(u, 8)?;
            let bits = Integer::from(m.bit_length());
            let negative: bool = u.arbitrary()?;
            let t = top(u, range)?;
            Ok(Float::from_integer_exp(
                &Integer::from_parts(negative, m),
                &(&t - &bits),
            ))
        }
    }
}

/// A radius for `mid` at `prec` bits: zero, about an ulp of `mid` or a
/// sizeable part of it, absolute, or (in the wide range) enormous.
pub fn radius(u: &mut Unstructured<'_>, mid: &Float, prec: u64, range: Range) -> Result<Mag> {
    let t = mid.exponent();
    let e = match u.int_in_range(0u8..=7)? {
        0..=2 => return Ok(Mag::ZERO),
        3 | 4 => t - &Integer::from(prec) - &Integer::from(u.int_in_range(-8i64..=80)?),
        5 => t - &Integer::from(u.int_in_range(-4i64..=64)?),
        6 => Integer::from(u.int_in_range(-300i64..=8)?),
        _ if range == Range::Wide => Integer::from(u.int_in_range(8i64..=1 << 20)?),
        _ => Integer::from(u.int_in_range(8i64..=20)?),
    };
    let m = Integer::from(u.int_in_range(1u64..=(1 << 30) - 1)?);
    Ok(Mag::from_float(&Float::from_integer_exp(
        &m,
        &(&e - &Integer::from(30u64)),
    )))
}

/// Points of the ball `b`: its midpoint, both ends and `m + t r` for a
/// fuzzer-chosen `t` in `[-1, 1]` (the midpoint alone for an exact ball,
/// none for an infinite radius).
///
/// Each is rounded toward the midpoint to at most `max(2 prec + 256,
/// bits(m) + 2)` bits, which keeps it in the ball (the midpoint fits that
/// many) and bounds MPFR's work: across a huge exponent gap an end has as
/// many bits as the gap, and MPFR's `atanh` spent minutes on one of a
/// million bits (CAS-113).
pub fn points(u: &mut Unstructured<'_>, b: &Ball, prec: u64) -> Result<Vec<Float>> {
    let m = b.mid();
    let Some(r) = b.rad().to_float() else {
        return Ok(Vec::new());
    };
    if r.is_zero() {
        return Ok(vec![m.clone()]);
    }
    let t = Float::from(u.int_in_range(-(1i64 << 16)..=1 << 16)?).mul_2exp(&Integer::from(-16i64));
    let cap = (2 * prec + 256).max(m.bit_length() + 2);
    let inward = |x: Float| {
        if x.bit_length() <= cap {
            return x;
        }
        let toward_m = if x < *m {
            OurRound::Ceil
        } else {
            OurRound::Floor
        };
        x.round(cap, toward_m).0
    };
    let ends = [
        m.checked_sub(&r),
        m.checked_add(&r),
        m.checked_add(&r.mul(&t)),
    ];
    Ok(std::iter::once(m.clone())
        .chain(ends.into_iter().flatten().map(inward))
        .collect())
}

/// MPFR's exact value of `x`, when its exponent is in MPFR's range.
pub fn to_mp(x: &Float) -> Option<MpFloat> {
    let (m, e) = x.to_integer_exp();
    let m = gmp_int(&m);
    let bits = m.significant_bits().max(1);
    let e = i32::try_from(e.to_i64()?).ok()?;
    let top = i32::try_from(bits).ok()?.checked_add(e)?;
    // The range in force (about ±2^30), not the largest one MPFR allows:
    // past it, the shift below would overflow to an infinity.
    let (lo, hi) = (rug::float::exp_min(), rug::float::exp_max());
    if m != 0 && !(lo..=hi).contains(&top) {
        return None;
    }
    let mut f = MpFloat::with_val(bits, m);
    f <<= e;
    Some(f)
}

/// Our float of a finite MPFR value; `None` for NaN and the infinities.
pub fn from_mp(x: &MpFloat) -> Option<Float> {
    if x.is_zero() {
        return Some(Float::ZERO);
    }
    let (m, e) = x.to_integer_exp()?;
    let digits = m.to_digits::<u64>(rug::integer::Order::Lsf);
    let m = Integer::from_parts(m < 0, Natural::from_limbs(&digits));
    Some(Float::from_integer_exp(&m, &Integer::from(i64::from(e))))
}

/// What MPFR says about a real value.
#[derive(Clone, Debug)]
pub enum Oracle {
    /// `[lo, hi]`, rounded toward −∞ and +∞, holds the value.
    Bracket(Float, Float),
    /// The point is outside the domain.
    Nan,
    /// A pole, or a value beyond MPFR's exponent range.
    Infinite,
}

/// MPFR's bracket from `f` at `prec` bits, rounded toward −∞ and +∞.
pub fn bracket(prec: u64, f: impl Fn(u32, Round) -> MpFloat) -> Oracle {
    let p = u32::try_from(prec).expect("precision below 2^32");
    let (lo, hi) = (f(p, Round::Down), f(p, Round::Up));
    if lo.is_nan() || hi.is_nan() {
        return Oracle::Nan;
    }
    match (from_mp(&lo), from_mp(&hi)) {
        (Some(lo), Some(hi)) => {
            assert!(lo <= hi, "oracle bracket reversed: {lo:?} > {hi:?}");
            Oracle::Bracket(lo, hi)
        }
        _ => Oracle::Infinite,
    }
}

/// Whether `b` meets `[lo, hi]`, decided exactly at any exponent: one of
/// the bracket's ends lies in the ball, or the ball lies inside the bracket.
pub fn meets(b: &Ball, lo: &Float, hi: &Float) -> bool {
    b.contains(lo) || b.contains(hi) || (lo <= b.mid() && b.mid() <= hi)
}

/// Fills the caches once, beyond every precision an input asks for, so
/// that results do not depend on earlier inputs: π and log 2 at `2^17`
/// bits (the reductions of arguments up to `2^65536`, past which they stop),
/// and the elementary functions' tables at their largest size.
pub fn warm_up() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let _ = Ball::const_pi(1 << 17);
        let _ = Ball::const_log2(1 << 17);
        let x = Ball::exact(Float::from(3i64).mul_2exp(&Integer::from(-3i64)));
        for f in [Ball::exp, Ball::log, Ball::sin, Ball::atan] {
            let _ = f(&x, 4600);
        }
    });
}
