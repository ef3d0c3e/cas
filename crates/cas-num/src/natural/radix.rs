//! Conversion to and from strings in radix 2 to 36.
//!
//! Power-of-two radices pack and unpack bits directly, in linear time. Other
//! radices work in chunks of `k` digits, where `radix^k` is the largest power
//! that fits a limb (the "big base"): the basecase divides by, or multiplies
//! by, the big base one limb pass at a time, in quadratic time. Large numbers
//! are split by divide and conquer on the powers `big_base^(2^i)`, so their
//! cost is that of a few multiplications or divisions of the full size.

use std::fmt;
use std::str::FromStr;

use super::Natural;
use crate::error::{ParseErrorKind, ParseNumberError};
use crate::{Limb, kernel};

/// Limbs from which [`Natural::to_string_radix`] splits by divide and conquer.
pub(crate) const TO_STRING_DC_THRESHOLD: usize = crate::tuning::TO_STRING_DC_THRESHOLD;

/// Digits from which parsing splits by divide and conquer.
pub(crate) const FROM_STRING_DC_THRESHOLD: usize = crate::tuning::FROM_STRING_DC_THRESHOLD;

const DIGITS: &[u8; 36] = b"0123456789abcdefghijklmnopqrstuvwxyz";

/// # Panics
///
/// If `radix` is not in `2..=36`.
pub(crate) fn check_radix(radix: u32) {
    assert!(
        (2..=36).contains(&radix),
        "radix must be in 2..=36, got {radix}"
    );
}

/// `(radix^k, k)` for the largest `k` with `radix^k` fitting a limb.
fn big_base(radix: u32) -> (Limb, usize) {
    let r = Limb::from(radix);
    let (mut base, mut k) = (r, 1);
    while let Some(next) = base.checked_mul(r) {
        base = next;
        k += 1;
    }
    (base, k)
}

#[cfg(test)]
mod tests;

impl Natural {
    /// The digits in `radix`, lowercase, without prefix or leading zeros
    /// (`"0"` for zero).
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    #[must_use]
    pub fn to_string_radix(&self, radix: u32) -> String {
        check_radix(radix);
        let mut digits = self.digits(radix);
        for d in &mut digits {
            *d = DIGITS[usize::from(*d)];
        }
        String::from_utf8(digits).expect("digits are ASCII")
    }

    /// Parses digits in `radix` (either case), with an optional leading `+`.
    /// Like the primitive integers' `from_str_radix`, nothing else is
    /// accepted: no whitespace, separators or prefix.
    ///
    /// # Errors
    ///
    /// [`ParseErrorKind::Empty`] for an empty string,
    /// [`ParseErrorKind::InvalidDigit`] for any other character, or for a
    /// sign with no digits.
    ///
    /// # Panics
    ///
    /// If `radix` is not in `2..=36`.
    pub fn from_str_radix(s: &str, radix: u32) -> Result<Natural, ParseNumberError> {
        check_radix(radix);
        let digits = match s.as_bytes() {
            [] => return Err(ParseNumberError::new(ParseErrorKind::Empty)),
            [b'+', rest @ ..] => rest,
            all => all,
        };
        Natural::from_digit_bytes(digits, radix)
    }

    /// Parses the unsigned digits `s` (no sign) in `radix`.
    pub(crate) fn from_digit_bytes(s: &[u8], radix: u32) -> Result<Natural, ParseNumberError> {
        let invalid = || ParseNumberError::new(ParseErrorKind::InvalidDigit);
        if s.is_empty() {
            return Err(invalid());
        }
        let mut values = Vec::with_capacity(s.len());
        for &c in s {
            let d = char::from(c).to_digit(radix).ok_or_else(invalid)?;
            values.push(u8::try_from(d).expect("a digit is below 36"));
        }
        let first = values.iter().position(|&d| d != 0).unwrap_or(values.len());
        Ok(from_digits(&values[first..], radix))
    }

    /// Digit values, most significant first, without leading zeros (`[0]`
    /// for zero).
    fn digits(&self, radix: u32) -> Vec<u8> {
        if self.is_zero() {
            return vec![0];
        }
        if radix.is_power_of_two() {
            return pow2_digits(self.limbs(), radix.trailing_zeros());
        }
        let (base, k) = big_base(radix);
        let digits = Digits { radix, base, k };
        let mut out = Vec::new();
        if self.limb_len() < TO_STRING_DC_THRESHOLD {
            digits.basecase(self.limbs(), None, &mut out);
        } else {
            let mut powers = vec![Natural::from(base)];
            while powers.last().expect("non-empty") <= self {
                let next = powers.last().expect("non-empty").square();
                powers.push(next);
            }
            digits.dc(self, &powers, powers.len() - 1, None, &mut out);
        }
        out
    }
}

/// Digit output in a radix that is not a power of two.
struct Digits {
    radix: u32,
    /// `radix^k`, the largest power that fits a limb.
    base: Limb,
    k: usize,
}

impl Digits {
    /// Appends the digits of `x` to `out`, left-padded with zeros to `pad`
    /// digits if given (`x` must fit), or without leading zeros otherwise.
    fn basecase(&self, x: &[Limb], pad: Option<usize>, out: &mut Vec<u8>) {
        if let [] | [_] = x {
            // A limb is below base^2: at most two chunks, with no division
            // pass and no allocation.
            let x = x.first().copied().unwrap_or(0);
            let chunks = [x % self.base, x / self.base];
            let n = if chunks[1] != 0 {
                2
            } else {
                usize::from(chunks[0] != 0)
            };
            return self.push_chunks(&chunks[..n], pad, out);
        }
        // Chunks of `k` digits, least significant first.
        let mut chunks = Vec::with_capacity(x.len() * 2);
        let mut a = x.to_vec();
        let mut t = vec![0; a.len()];
        while !a.is_empty() {
            let n = a.len();
            chunks.push(kernel::div_rem_1(&mut t[..n], &a, self.base));
            std::mem::swap(&mut a, &mut t);
            a.truncate(n);
            a.truncate(kernel::normalized_len(&a));
        }
        self.push_chunks(&chunks, pad, out);
    }

    /// Appends the digits of the number with chunks `chunks` (least
    /// significant first, the top one non-zero), padded as for
    /// [`Digits::basecase`].
    fn push_chunks(&self, chunks: &[Limb], pad: Option<usize>, out: &mut Vec<u8>) {
        let mut buf = [0; Limb::BITS as usize];
        let Some((&top, rest)) = chunks.split_last() else {
            out.resize(out.len() + pad.unwrap_or(0), 0);
            return;
        };
        let top = self.chunk_digits(top, &mut buf);
        let len = top.len() + self.k * rest.len();
        if let Some(pad) = pad {
            assert!(len <= pad, "value does not fit its padding");
            out.resize(out.len() + pad - len, 0);
        }
        out.extend_from_slice(top);
        for &c in rest.iter().rev() {
            let d = self.chunk_digits(c, &mut buf);
            out.resize(out.len() + self.k - d.len(), 0);
            out.extend_from_slice(d);
        }
    }

    /// The digits of a chunk `c < base` without leading zeros (none for
    /// zero), written to the end of `buf`.
    ///
    /// Radix 10 gets its own copy of the loop, where the division by the
    /// constant compiles to a multiplication; a division by a variable is a
    /// hardware `div` per digit, which dominated small conversions.
    fn chunk_digits<'b>(&self, c: Limb, buf: &'b mut [u8; Limb::BITS as usize]) -> &'b [u8] {
        match self.radix {
            10 => chunk_digits_in(c, 10, buf),
            r => chunk_digits_in(c, Limb::from(r), buf),
        }
    }

    /// Appends the digits of `x < powers[level]` (`powers[i] = base^(2^i)`),
    /// padded as for [`Digits::basecase`].
    fn dc(
        &self,
        x: &Natural,
        powers: &[Natural],
        level: usize,
        pad: Option<usize>,
        out: &mut Vec<u8>,
    ) {
        if level == 0 || x.limb_len() < TO_STRING_DC_THRESHOLD {
            return self.basecase(x.limbs(), pad, out);
        }
        // q, r < powers[level - 1], and r has exactly `low` digits.
        let (q, r) = x.div_rem(&powers[level - 1]);
        let low = self.k << (level - 1);
        if pad.is_none() && q.is_zero() {
            return self.dc(&r, powers, level - 1, None, out);
        }
        self.dc(&q, powers, level - 1, pad.map(|p| p - low), out);
        self.dc(&r, powers, level - 1, Some(low), out);
    }
}

#[allow(
    clippy::inline_always,
    reason = "inlined into each call, so that a constant radix divides by multiplication"
)]
#[inline(always)]
fn chunk_digits_in(mut c: Limb, r: Limb, buf: &mut [u8; Limb::BITS as usize]) -> &[u8] {
    let mut i = buf.len();
    while c != 0 {
        i -= 1;
        buf[i] = u8::try_from(c % r).expect("a digit is below 36");
        c /= r;
    }
    &buf[i..]
}

/// Digits of non-zero `x` in radix `2^bits`.
fn pow2_digits(x: &[Limb], bits: u32) -> Vec<u8> {
    let width = u64::from(bits);
    let n = kernel::bit_length(x).div_ceil(width);
    let mask = (1 << bits) - 1;
    (0..n)
        .rev()
        .map(|i| {
            let pos = i * width;
            let limb = usize::try_from(pos / 64).expect("bit index fits usize");
            let shift = pos % 64;
            let mut v = x[limb] >> shift;
            if shift + width > 64 {
                if let Some(next) = x.get(limb + 1) {
                    v |= next << (64 - shift);
                }
            }
            u8::try_from(v & mask).expect("masked to a digit")
        })
        .collect()
}

/// The value of digit values `d` (most significant first) in `radix`.
fn from_digits(d: &[u8], radix: u32) -> Natural {
    if d.is_empty() {
        return Natural::ZERO;
    }
    if radix.is_power_of_two() {
        return pow2_from_digits(d, radix.trailing_zeros());
    }
    let (base, k) = big_base(radix);
    if d.len() < FROM_STRING_DC_THRESHOLD {
        return parse_basecase(d, radix, k);
    }
    let mut powers = vec![Natural::from(base)];
    while k << powers.len() < d.len() {
        let next = powers.last().expect("non-empty").square();
        powers.push(next);
    }
    parse_dc(d, radix, k, &powers)
}

/// Digits split at `k 2^i` from the right, for the largest `i` leaving a
/// non-empty high part; `powers[i] = base^(2^i)`.
fn parse_dc(d: &[u8], radix: u32, k: usize, powers: &[Natural]) -> Natural {
    if d.len() < FROM_STRING_DC_THRESHOLD {
        return parse_basecase(d, radix, k);
    }
    let level = (0..powers.len())
        .rev()
        .find(|&i| k << i < d.len())
        .expect("a single chunk is below the threshold");
    let (high, low) = d.split_at(d.len() - (k << level));
    &parse_dc(high, radix, k, powers) * &powers[level] + parse_dc(low, radix, k, powers)
}

/// Quadratic parsing, `k` digits (one limb multiplication) at a time.
fn parse_basecase(d: &[u8], radix: u32, k: usize) -> Natural {
    let r = Limb::from(radix);
    let mut acc: Vec<Limb> = Vec::with_capacity(d.len() / k + 1);
    let first = match d.len() % k {
        0 => k,
        m => m,
    };
    let (head, rest) = d.split_at(first.min(d.len()));
    for chunk in std::iter::once(head).chain(rest.chunks(k)) {
        let (mut value, mut scale) = (0, 1);
        for &digit in chunk {
            value = value * r + Limb::from(digit);
            scale *= r;
        }
        let high = kernel::mul_1_assign(&mut acc, scale);
        let carry = kernel::add_1_assign(&mut acc, value);
        // acc scale + value < (acc + 1) scale, so this does not overflow.
        let top = high + carry;
        if top != 0 {
            acc.push(top);
        }
    }
    Natural::from_vec(acc)
}

/// The value of digit values `d` in radix `2^bits`.
fn pow2_from_digits(d: &[u8], bits: u32) -> Natural {
    let width = u64::from(bits);
    let total = u64::try_from(d.len()).expect("length fits u64") * width;
    let mut limbs = vec![0; usize::try_from(total.div_ceil(64)).expect("fits usize")];
    for (i, &digit) in d.iter().rev().enumerate() {
        let pos = u64::try_from(i).expect("index fits u64") * width;
        let limb = usize::try_from(pos / 64).expect("bit index fits usize");
        let shift = pos % 64;
        let v = Limb::from(digit);
        limbs[limb] |= v << shift;
        if shift + width > 64 {
            limbs[limb + 1] |= v >> (64 - shift);
        }
    }
    Natural::from_vec(limbs)
}

impl fmt::Display for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(true, "", &self.to_string_radix(10))
    }
}

impl fmt::Debug for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

impl fmt::Binary for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(true, "0b", &self.to_string_radix(2))
    }
}

impl fmt::Octal for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad_integral(true, "0o", &self.to_string_radix(8))
    }
}

impl FromStr for Natural {
    type Err = ParseNumberError;

    /// Decimal, as [`Natural::from_str_radix`] with radix 10.
    fn from_str(s: &str) -> Result<Natural, ParseNumberError> {
        Natural::from_str_radix(s, 10)
    }
}
