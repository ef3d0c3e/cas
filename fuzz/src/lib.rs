//! Input decoding shared by the fuzz targets: numbers built from fuzzer
//! bytes, as cas-num values paired with the same values in GMP (via `rug`).
//!
//! Limbs are drawn from a mix of boundary patterns and raw bytes, and
//! lengths are weighted towards the algorithm thresholds, where the
//! interesting cases live.

use arbitrary::{Result, Unstructured};
use cas_num::{Integer, Natural, Rational};

pub mod ball;

/// Lengths near which cas-num switches algorithms (Karatsuba 16, divide and
/// conquer division 40, Karatsuba squaring 44, Toom-3 144, REDC 144).
const THRESHOLDS: [usize; 5] = [16, 40, 44, 144, 192];

/// One limb: a boundary value or raw bytes.
pub fn limb(u: &mut Unstructured<'_>) -> Result<u64> {
    Ok(match u.int_in_range(0u8..=7)? {
        0 => 0,
        1 => u64::MAX,
        2 => 1 << 63,
        3 => 1,
        4 => (1 << u.int_in_range(0u32..=63)?) - 1,
        _ => u.arbitrary()?,
    })
}

/// A length in `0..=max`: mostly short, sometimes next to a threshold.
pub fn len(u: &mut Unstructured<'_>, max: usize) -> Result<usize> {
    let n = match u.int_in_range(0u8..=7)? {
        0 => {
            let t = *u.choose(&THRESHOLDS)?;
            t + u.int_in_range(0..=2)? - 1
        }
        1 => u.int_in_range(0..=max)?,
        _ => u.int_in_range(0..=4)?,
    };
    Ok(n.min(max))
}

/// A natural number of at most `max` limbs, with a uniform fill, all-ones
/// limbs, or mostly zero limbs.
pub fn natural(u: &mut Unstructured<'_>, max: usize) -> Result<(Natural, rug::Integer)> {
    let n = len(u, max)?;
    let fill = u.int_in_range(0u8..=3)?;
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(match fill {
            0 => u64::MAX,
            1 if u.ratio(3u8, 4)? => 0,
            _ => limb(u)?,
        });
    }
    let ours = Natural::from_limbs(&v);
    let gmp = rug::Integer::from_digits(&v, rug::integer::Order::Lsf);
    Ok((ours, gmp))
}

/// A signed integer of at most `max` limbs, often within `i64`.
pub fn integer(u: &mut Unstructured<'_>, max: usize) -> Result<(Integer, rug::Integer)> {
    if u.ratio(1u8, 3)? {
        let x: i64 = u.arbitrary()?;
        return Ok((Integer::from(x), rug::Integer::from(x)));
    }
    let (m, g) = natural(u, max)?;
    let negative: bool = u.arbitrary()?;
    let ours = Integer::from_parts(negative, m);
    Ok(if negative { (ours, -g) } else { (ours, g) })
}

/// A rational, reduced by both sides; `None` for a zero denominator.
pub fn rational(u: &mut Unstructured<'_>, max: usize) -> Result<Option<(Rational, rug::Rational)>> {
    let (n, ng) = integer(u, max)?;
    let (d, dg) = integer(u, max)?;
    if d.is_zero() {
        return Ok(None);
    }
    Ok(Some((Rational::new(n, d), rug::Rational::from((ng, dg)))))
}

/// GMP's value of a cas-num natural number, converted through limbs.
pub fn gmp(x: &Natural) -> rug::Integer {
    rug::Integer::from_digits(x.limbs(), rug::integer::Order::Lsf)
}

/// GMP's value of a cas-num integer.
pub fn gmp_int(x: &Integer) -> rug::Integer {
    let (negative, mag) = x.clone().into_parts();
    let m = gmp(&mag);
    if negative { -m } else { m }
}

/// GMP's value of a cas-num rational.
pub fn gmp_rat(x: &Rational) -> rug::Rational {
    rug::Rational::from((gmp_int(x.numer()), gmp_int(x.denom())))
}

/// Panics unless the two sides agree, printing the operation.
#[track_caller]
pub fn check<T: PartialEq + std::fmt::Debug>(op: &str, ours: T, gmp: T) {
    assert!(ours == gmp, "{op}: ours {ours:?}, gmp {gmp:?}");
}

/// Largest exponent magnitude of [`dyadic`], so exact results stay small.
pub const DYADIC_MAX_EXP: i64 = 512;

/// A dyadic rational `m·2^e`: a mantissa of at most `max` limbs and
/// `|e| <= DYADIC_MAX_EXP`, the exact values ball midpoints and radii take.
pub fn dyadic(u: &mut Unstructured<'_>, max: usize) -> Result<(Rational, rug::Rational)> {
    let (m, mg) = integer(u, max)?;
    let e = u.int_in_range(-DYADIC_MAX_EXP..=DYADIC_MAX_EXP)?;
    let pow = Integer::from_parts(false, Natural::from(1u64) << e.unsigned_abs());
    let one = Integer::from(1);
    let ours = if e >= 0 {
        Rational::new(m * pow, one)
    } else {
        Rational::new(m, pow)
    };
    let mut gmp = rug::Rational::from(mg);
    if e >= 0 {
        gmp <<= u32::try_from(e).unwrap();
    } else {
        gmp >>= u32::try_from(-e).unwrap();
    }
    Ok((ours, gmp))
}
