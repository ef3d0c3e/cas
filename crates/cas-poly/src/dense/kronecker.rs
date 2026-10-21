//! Kronecker substitution \[Sch82\]: a product of integer polynomials as one
//! product of integers. With every coefficient of the product below
//! `2^(w-1)` in magnitude, `a(2^w) b(2^w) = c(2^w)`, and the coefficients of
//! `c` are the digits of that integer in base `2^w`, read with signs: each
//! in `[-2^(w-1), 2^(w-1))`. So the polynomial product gets cas-num's integer
//! products, Toom and NTT, as in FLINT's `fmpz_poly_mul_KS`.

use cas_num::modular::ModWord;
use cas_num::{Integer, Natural};

/// Bits of the largest magnitude among `coeffs`.
pub(crate) fn max_bits(coeffs: &[Integer]) -> u64 {
    coeffs.iter().map(Integer::bit_length).max().unwrap_or(0)
}

/// The most coefficients in the shorter factor for which schoolbook beats
/// Kronecker substitution over ℤ, when a product of two coefficients has
/// `bits` bits. Measured in instructions (CAS-135): with products that fit
/// a word, cas-num's small integers make schoolbook cheap up to about 13
/// coefficients; past a word it loses from 3 or 4.
pub(crate) fn schoolbook_up_to(bits: u64) -> usize {
    match bits {
        ..=64 => 12,
        65..=96 => 4,
        _ => 2,
    }
}

/// As [`schoolbook_up_to`], modulo a word `m`: the digits of the integer
/// product grow with the modulus, from about 17 coefficients for a 17-bit
/// prime to 56 for a 64-bit one.
pub(crate) fn word_schoolbook_up_to(m: &ModWord) -> usize {
    match 64 - m.modulus().leading_zeros() {
        ..=20 => 16,
        21..=32 => 20,
        33..=62 => 36,
        _ => 56,
    }
}

/// The width of a digit: the product's coefficients are sums of at most
/// `min(n, m)` products of an `a`-bit and a `b`-bit magnitude, so below
/// `2^(a + b + bits(min(n, m)))`, and one more bit holds the sign.
fn width(a_bits: u64, b_bits: u64, terms: usize) -> u64 {
    a_bits + b_bits + u64::from(usize::BITS - terms.leading_zeros()) + 1
}

/// `dst |= src << offset`, in bits; `dst` long enough.
fn or_shifted(dst: &mut [u64], src: &[u64], offset: u64) {
    let limb = usize::try_from(offset / 64).expect("fits");
    let shift = offset % 64;
    if shift == 0 {
        for (d, s) in dst[limb..].iter_mut().zip(src) {
            *d |= s;
        }
    } else {
        for (i, &s) in src.iter().enumerate() {
            dst[limb + i] |= s << shift;
            dst[limb + i + 1] |= s >> (64 - shift);
        }
    }
}

/// The limbs of `n · 2^(w·count)` bits, room for packing `count` digits.
fn room(count: usize, w: u64) -> Vec<u64> {
    let bits = w * u64::try_from(count).expect("fits");
    vec![0; usize::try_from(bits.div_ceil(64)).expect("fits") + 2]
}

/// `p(2^w)` for integer coefficients, each below `2^w` in magnitude: the
/// positive and the negative coefficients packed apart, without carries,
/// then one subtraction.
fn pack(coeffs: &[Integer], w: u64) -> Integer {
    let mut positive = room(coeffs.len(), w);
    let mut negative: Option<Vec<u64>> = None;
    for (i, c) in coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let at = w * u64::try_from(i).expect("fits");
        let magnitude = c.unsigned_abs();
        if c.is_negative() {
            let dst = negative.get_or_insert_with(|| room(coeffs.len(), w));
            or_shifted(dst, magnitude.limbs(), at);
        } else {
            or_shifted(&mut positive, magnitude.limbs(), at);
        }
    }
    let p = Integer::from(Natural::from_vec(positive));
    match negative {
        Some(n) => p - Integer::from(Natural::from_vec(n)),
        None => p,
    }
}

/// Bits `start .. start + w` of the number of `limbs`, as a natural.
fn digit(limbs: &[u64], start: u64, w: u64) -> Natural {
    let first = usize::try_from(start / 64).expect("fits");
    let shift = start % 64;
    let len = usize::try_from(w.div_ceil(64)).expect("fits");
    let at = |i: usize| limbs.get(i).copied().unwrap_or(0);
    let mut out: Vec<u64> = (0..len)
        .map(|i| {
            let lo = at(first + i) >> shift;
            if shift == 0 {
                lo
            } else {
                lo | at(first + i + 1) << (64 - shift)
            }
        })
        .collect();
    if w % 64 != 0 {
        out[len - 1] &= (1u64 << (w % 64)) - 1;
    }
    Natural::from_vec(out)
}

/// The `count` coefficients whose value at `2^w` is `c`, each in
/// `[-2^(w-1), 2^(w-1))`: the digits of `|c|` in base `2^w`, a digit from
/// `2^(w-1)` up taken as negative with a carry into the next, then the sign
/// of `c` on all.
fn unpack(c: &Integer, w: u64, count: usize) -> Vec<Integer> {
    let negative = c.is_negative();
    let magnitude = c.unsigned_abs();
    let limbs = magnitude.limbs();
    let base = Natural::power_of_two(w);
    let mut carry = false;
    let mut out = Vec::with_capacity(count);
    for k in 0..count {
        let mut d = digit(limbs, w * u64::try_from(k).expect("fits"), w);
        if carry {
            d += 1;
        }
        // From 2^(w-1) up, the digit is d - 2^w, and 1 carries.
        carry = d.bit_length() >= w;
        let value = if carry {
            Integer::from_parts(true, &base - &d)
        } else {
            Integer::from(d)
        };
        out.push(if negative { -value } else { value });
    }
    debug_assert!(!carry, "a coefficient past the last");
    out
}

/// The coefficients of `a b`, `a.len() + b.len() - 1` of them.
pub(crate) fn mul(a: &[Integer], b: &[Integer]) -> Vec<Integer> {
    mul_bits(a, b, max_bits(a), max_bits(b))
}

/// [`mul`] with the [`max_bits`] of `a` and `b` known.
pub(crate) fn mul_bits(a: &[Integer], b: &[Integer], a_bits: u64, b_bits: u64) -> Vec<Integer> {
    let w = width(a_bits, b_bits, a.len().min(b.len()));
    let product = pack(a, w) * pack(b, w);
    unpack(&product, w, a.len() + b.len() - 1)
}

/// The coefficients of `a²`: one integer square.
#[cfg(test)]
pub(crate) fn square(a: &[Integer]) -> Vec<Integer> {
    square_bits(a, max_bits(a))
}

/// [`square`] with the [`max_bits`] of `a` known.
pub(crate) fn square_bits(a: &[Integer], bits: u64) -> Vec<Integer> {
    let w = width(bits, bits, a.len());
    let packed = pack(a, w);
    let product = &packed * &packed;
    unpack(&product, w, 2 * a.len() - 1)
}

/// The coefficients of `a b` modulo `m`, for residues in `0..m`: their
/// integer product's coefficients are non-negative and below
/// `min(n, m) (m - 1)²`, so the digits need no sign, and each is reduced.
pub(crate) fn mul_mod(a: &[u64], b: &[u64], m: &ModWord) -> Vec<u64> {
    let bits = u64::from(64 - m.modulus().leading_zeros());
    let terms = a.len().min(b.len());
    let w = 2 * bits + u64::from(usize::BITS - terms.leading_zeros());
    let pack = |p: &[u64]| {
        let mut limbs = room(p.len(), w);
        for (i, &c) in p.iter().enumerate() {
            or_shifted(&mut limbs, &[c], w * u64::try_from(i).expect("fits"));
        }
        Natural::from_vec(limbs)
    };
    let product = &pack(a) * &pack(b);
    let limbs = product.limbs();
    (0..a.len() + b.len() - 1)
        .map(|k| m.reduce_natural(&digit(limbs, w * u64::try_from(k).expect("fits"), w)))
        .collect()
}
