//! Decimal printing and parsing of balls (CAS-DOC-5 §6.1, §6.3).
//!
//! Printing follows Arb's `arb_printn`: the digits shown are correct to
//! within one unit in the last place, for every point of the ball. So
//! `1 ± 10^-30` prints as `1.000…` rather than nothing, which is what a
//! "only certain digits" rule would give at a digit boundary. When not even
//! one digit is known, the ball prints as `[± R]`, `R` bounding every point's
//! magnitude; [`DecimalStyle::Radius`] always shows the radius, as
//! `[digits ± R]` where `R` bounds the distance from the printed decimal.
//!
//! The digits come from one ball product `x · 10^s` at a few more bits than
//! the digits need: `D`, its midpoint rounded to an integer, is good when
//! the product's ball lies within 1 of it. `10^s` is powered as a ball, so
//! exponents far beyond what an exact power of 10 could hold cost only a few
//! squarings.
//!
//! Parsing takes `[+-]digits[.digits][(e|E)[+-]digits]` (a digit on either
//! side of the point) as the exact decimal it spells, rounded to nearest at
//! the requested precision with a radius of half an ulp (0 when exact).
#![allow(
    clippy::many_single_char_names,
    reason = "x, z, D, s, k, n, e as in the formulas of the comments"
)]

use std::fmt;

use super::{Ball, Float, Mag, Round};
use crate::error::{ParseErrorKind, ParseNumberError};
use crate::{Integer, Natural};

#[cfg(test)]
mod tests;

/// How [`Ball::to_decimal`] shows the radius.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum DecimalStyle {
    /// Only digits correct to ±1 in the last place; `[± R]` when there are
    /// none.
    #[default]
    Digits,
    /// The digits and a bound on their distance to every point,
    /// `[3.14159 ± 2.66e-6]`; exact decimals print without it.
    Radius,
}

/// `floor(log10(2) · 2^64)`.
const LOG10_2_FIXED: i128 = 0x4d10_4d42_7de7_fbcc;

/// Significant digits of a printed radius, rounded up.
const RADIUS_DIGITS: u64 = 3;

/// `10^n` as an exact integer.
fn pow10(n: u64) -> Integer {
    Integer::from(10u64).pow(n)
}

/// A ball containing `10^s`, exact while it fits `prec` bits.
fn ten_pow(s: i64, prec: u64) -> Ball {
    Ball::from_integer(&Integer::from(10u64)).pow(s, prec)
}

/// `x · 10^s` at `prec` bits (a division for `s < 0`, which keeps exact
/// quotients exact).
fn scale(x: &Ball, s: i64, prec: u64) -> Ball {
    if s >= 0 {
        x.mul(&ten_pow(s, prec), prec)
    } else {
        x.div(&ten_pow(-s, prec), prec)
    }
}

/// `floor(log10 |x|)` within one: for `|x|` in `[2^(e-1), 2^e)` it is
/// `floor((e - 1) log10 2)` or one more; `None` beyond `i64` exponents.
fn log10_estimate(x: &Float) -> Option<i64> {
    let e = x.exponent().to_i64()?;
    let t = (i128::from(e) - 1) * LOG10_2_FIXED;
    i64::try_from(t >> 64).ok()
}

/// `x` rounded to an integer: to nearest (ties away from zero), or up
/// (toward `+∞`).
pub(super) fn to_integer(x: &Float, up: bool) -> Integer {
    let (m, e) = x.to_integer_exp();
    if !e.is_negative() {
        let (neg, mag) = m.into_parts();
        let k = e
            .to_i64()
            .and_then(|k| u64::try_from(k).ok())
            .expect("bounded");
        return Integer::from_parts(neg, mag << k);
    }
    let k = (-e).to_i64().and_then(|k| u64::try_from(k).ok());
    let Some(k) = k.filter(|&k| k <= m.bit_length() + 1) else {
        // |x| < 1/2.
        return Integer::from(i64::from(up && !m.is_negative() && !m.is_zero()));
    };
    let d = Integer::from(Natural::from(1u64) << k);
    if up {
        m.div_rem_ceil(&d).0
    } else {
        let half = Integer::from(Natural::from(1u64) << (k - 1));
        let shifted = if m.is_negative() {
            &m - &half
        } else {
            &m + &half
        };
        shifted.div_rem_trunc(&d).0
    }
}

/// `k` significant digits of `x`: `(D, digits of |D|, e10, z)` with `|D|`
/// of exactly `k` digits, `x ≈ D · 10^(e10 - k + 1)`, and the ball `z = x ·
/// 10^(k - 1 - e10)` whose midpoint rounded is `D`. `up` rounds the upper
/// end up instead (for radii). `None` beyond `i64` exponents.
fn digits(x: &Ball, k: u64, up: bool) -> Option<(Integer, String, i64, Ball)> {
    let mut e10 = log10_estimate(x.mid())?;
    let prec = k * 3322 / 1000 + 64;
    let want = usize::try_from(k).ok()?;
    for _ in 0..4 {
        let s = i64::try_from(k).ok()?.checked_sub(1)?.checked_sub(e10)?;
        let z = scale(x, s, prec);
        let d = if up {
            to_integer(&z.mid().add(&z.rad().to_float()?), true)
        } else {
            to_integer(z.mid(), false)
        };
        // The string is needed anyway; its length is the digit count.
        let ds = d.unsigned_abs().to_string();
        match ds.len().cmp(&want) {
            std::cmp::Ordering::Greater => e10 += 1,
            std::cmp::Ordering::Less => e10 -= 1,
            std::cmp::Ordering::Equal => return Some((d, ds, e10, z)),
        }
    }
    None
}

/// `|z - D| + rad(z)`, the distance from the digits to every point, in
/// units of the last digit.
fn distance(z: &Ball, d: &Integer) -> Mag {
    Mag::from_float(&z.mid().sub(&Float::from(d))).add(z.rad())
}

/// The digits `ds` (with a minus sign if `neg`) with the first at decimal exponent
/// `e10`, in fixed notation from `1e-4` to the last digit shown (to `n`
/// digits for exact values, which may show exact zeros), scientific
/// otherwise. Exact values drop trailing zeros of the fraction.
fn format(neg: bool, mut ds: String, e10: i64, exact: bool, n: u64) -> String {
    if exact {
        let keep = ds.trim_end_matches('0').len().max(1);
        ds.truncate(keep);
    }
    let sign = if neg { "-" } else { "" };
    let len = i64::try_from(ds.len()).expect("fits");
    let width = if exact {
        i64::try_from(n).unwrap_or(i64::MAX).max(len)
    } else {
        len
    };
    if (-4..0).contains(&e10) {
        let zeros = usize::try_from(-e10 - 1).expect("small");
        return format!("{sign}0.{}{ds}", "0".repeat(zeros));
    }
    if (0..width).contains(&e10) {
        let int = usize::try_from(e10 + 1).expect("small");
        if int >= ds.len() {
            return format!("{sign}{ds}{}", "0".repeat(int - ds.len()));
        }
        return format!("{sign}{}.{}", &ds[..int], &ds[int..]);
    }
    let (first, rest) = ds.split_at(1);
    let dot = if rest.is_empty() { "" } else { "." };
    let esign = if e10 < 0 { '-' } else { '+' };
    format!("{sign}{first}{dot}{rest}e{esign}{}", e10.unsigned_abs())
}

/// `m` rounded up to `RADIUS_DIGITS` digits, `inf` for `+∞`.
fn format_upper(m: &Mag) -> String {
    let Some(x) = m.to_float() else {
        return "inf".to_owned();
    };
    if x.is_zero() {
        return "0".to_owned();
    }
    match digits(&Ball::exact(x.clone()), RADIUS_DIGITS, true) {
        Some((_, ds, e10, _)) => format(false, ds, e10, true, RADIUS_DIGITS),
        None => x.to_hex(),
    }
}

impl Ball {
    /// The ball in decimal with at most `n` significant digits, each correct
    /// to within one unit in the last place for every point of the ball (Arb
    /// semantics), as [`DecimalStyle`] says. Exact decimals drop trailing
    /// zeros. A ball whose exponent is beyond `i64` prints its midpoint and
    /// radius as hexadecimal floats.
    ///
    /// # Panics
    ///
    /// If `n` is 0.
    #[must_use]
    pub fn to_decimal(&self, n: u64, style: DecimalStyle) -> String {
        assert!(n >= 1, "Ball::to_decimal: no digits");
        if !self.is_finite() {
            return "[± inf]".to_owned();
        }
        if self.mid().is_zero() {
            return if self.is_exact() {
                "0".to_owned()
            } else {
                format!("[± {}]", format_upper(self.rad()))
            };
        }
        if self.mid().exponent().to_i64().is_none() {
            return format!("[{} ± {}]", self.mid().to_hex(), self.rad().to_hex());
        }
        // Digits known: about log10(|m| / r), one more to try first.
        let known = match self.rel_accuracy_bits() {
            i64::MAX => n,
            bits => u64::try_from((i128::from(bits) * LOG10_2_FIXED) >> 64)
                .map_or(0, |k| k.saturating_add(2))
                .min(n),
        };
        for k in (1..=known).rev() {
            let Some((d, ds, e10, z)) = digits(self, k, false) else {
                break;
            };
            let dist = distance(&z, &d);
            if dist > Mag::ONE {
                continue;
            }
            let exact = dist.is_zero();
            let text = format(d.is_negative(), ds, e10, exact, n);
            if style == DecimalStyle::Digits || exact {
                return text;
            }
            // The distance back in the ball's units: dist · 10^(e10 - k + 1).
            let unit = i64::try_from(k).map_or(i64::MIN, |k| e10 - k + 1);
            let r = scale(&Ball::exact(dist.to_float().expect("finite")), unit, 64);
            let r = Mag::from_float(r.mid()).add(r.rad());
            return format!("[{text} ± {}]", format_upper(&r));
        }
        let bound = Mag::from_float(self.mid()).add(self.rad());
        format!("[± {}]", format_upper(&bound))
    }

    /// The ball of a decimal number at `prec` bits: the exact decimal value
    /// rounded to nearest, with a radius of half an ulp, or 0 when the value
    /// is exact. Exponents whose exact power of 10 would exceed
    /// `max(2^20, 4 prec)` bits are scaled by ball arithmetic instead, with
    /// a radius a little above half an ulp.
    ///
    /// # Errors
    ///
    /// If `s` is empty or not `[+-]digits[.digits][(e|E)[+-]digits]` (a
    /// digit on at least one side of the point), or its exponent is beyond
    /// `i64`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    pub fn from_decimal(s: &str, prec: u64) -> Result<Ball, ParseNumberError> {
        assert!(prec >= 1, "Ball::from_decimal: precision 0");
        let invalid = || ParseNumberError::new(ParseErrorKind::InvalidDigit);
        if s.is_empty() {
            return Err(ParseNumberError::new(ParseErrorKind::Empty));
        }
        let (body, exp) = match s.find(['e', 'E']) {
            Some(i) => (&s[..i], Some(&s[i + 1..])),
            None => (s, None),
        };
        let (neg, body) = match body.as_bytes().first() {
            Some(b'-') => (true, &body[1..]),
            Some(b'+') => (false, &body[1..]),
            _ => (false, body),
        };
        let (int, frac) = body.split_once('.').unwrap_or((body, ""));
        let digit = |t: &str| t.bytes().all(|c| c.is_ascii_digit());
        if int.len() + frac.len() == 0 || !digit(int) || !digit(frac) {
            return Err(invalid());
        }
        let e = match exp {
            None => 0,
            Some(t) => {
                let unsigned = t.strip_prefix(['+', '-']).unwrap_or(t);
                if unsigned.is_empty() || !digit(unsigned) {
                    return Err(invalid());
                }
                Integer::from_str_radix(t.strip_prefix('+').unwrap_or(t), 10)?
                    .to_i64()
                    .ok_or_else(|| ParseNumberError::new(ParseErrorKind::ExponentRange))?
            }
        };
        let mut all = String::with_capacity(int.len() + frac.len());
        all.push_str(int);
        all.push_str(frac);
        let n = Integer::from_parts(neg, Natural::from_str_radix(&all, 10)?);
        let e = i64::try_from(frac.len())
            .ok()
            .and_then(|f| e.checked_sub(f))
            .ok_or_else(|| ParseNumberError::new(ParseErrorKind::ExponentRange))?;
        if n.is_zero() {
            return Ok(Ball::ZERO);
        }
        // value = n · 10^e.
        let cap = prec.saturating_mul(4).max(1 << 20);
        let fits = e.unsigned_abs().saturating_mul(3322) / 1000 <= cap;
        if !fits {
            let w = prec + 64;
            let x = Ball::from_integer(&n).round_to_accuracy(w);
            return Ok(scale(&x, e, w).round_to_accuracy(prec));
        }
        let (m, o) = if e >= 0 {
            Float::from(&(&n * &pow10(e.unsigned_abs()))).round(prec, Round::Nearest)
        } else {
            Float::from(&n).div_round(&Float::from(&pow10(e.unsigned_abs())), prec, Round::Nearest)
        };
        Ok(Ball::from_nearest(m, o, prec))
    }
}

/// `{}` prints with [`Ball::to_decimal`] in [`DecimalStyle::Digits`], `{:#}`
/// in [`DecimalStyle::Radius`]; the precision `{:.n}` sets the digits,
/// by default those of the midpoint's bit length, at least 17 (as for an
/// `f64`; an inexact ball shows only the digits it knows anyway).
impl fmt::Display for Ball {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = match f.precision() {
            Some(n) => u64::try_from(n).unwrap_or(u64::MAX).max(1),
            None => (self.mid().bit_length() * 30_103 / 100_000 + 1).max(17),
        };
        let style = if f.alternate() {
            DecimalStyle::Radius
        } else {
            DecimalStyle::Digits
        };
        f.write_str(&self.to_decimal(n, style))
    }
}
