//! Square roots, nth roots and perfect powers.
//!
//! The square root is Zimmermann's Karatsuba square root (Brent and
//! Zimmermann, *Modern Computer Arithmetic*, Algorithm 1.12): the root of the
//! top half, recursively, then one division gives the low half of the root.
//! It costs a constant times a division of the same size. The nth root is
//! Newton's iteration from an overestimate, which is the root of the top
//! part, recursively, so a couple of iterations finish each level.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (m, s, r, q, u, k, x, n, p)"
)]

use super::Natural;
use super::bits::split_bit_index;
use crate::Limb;
use crate::modular::ModWord;
use crate::prime::is_prime_u64;

#[cfg(test)]
mod tests;

/// Odd primes whose product fits a limb; a square is a quadratic residue
/// modulo each.
const RESIDUE_PRIMES: [Limb; 14] = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];

const RESIDUE_PRODUCT: Limb = {
    let mut p = 1;
    let mut i = 0;
    while i < RESIDUE_PRIMES.len() {
        p *= RESIDUE_PRIMES[i];
        i += 1;
    }
    p
};

/// Bit `x` is set iff `x` is a square modulo `m`, for `m <= 64`.
const fn square_mask(m: Limb) -> u64 {
    let mut mask = 0;
    let mut x = 0;
    while x < m {
        mask |= 1 << (x * x % m);
        x += 1;
    }
    mask
}

const SQUARES_MOD_64: u64 = square_mask(64);

const SQUARE_MASKS: [u64; RESIDUE_PRIMES.len()] = {
    let mut masks = [0; RESIDUE_PRIMES.len()];
    let mut i = 0;
    while i < masks.len() {
        masks[i] = square_mask(RESIDUE_PRIMES[i]);
        i += 1;
    }
    masks
};

/// Primes `q ≡ 1 (mod p)` used to reject non-`p`th powers.
const FILTER_PRIMES: usize = 3;

impl Natural {
    /// `⌊√self⌋`.
    #[must_use]
    pub fn isqrt(&self) -> Natural {
        self.sqrt_rem().0
    }

    /// `(s, r)` with `s = ⌊√self⌋` and `r = self - s²`, so `0 <= r <= 2s`.
    #[must_use]
    pub fn sqrt_rem(&self) -> (Natural, Natural) {
        if let Some(x) = self.to_u128() {
            let s = x.isqrt();
            return (Natural::from(s), Natural::from(x - s * s));
        }
        // The algorithm needs a length of 4k - 1 or 4k bits; otherwise take
        // the root of 4 self, (2s + b)² + r' = 4 (s² + r).
        if matches!(self.bit_length() % 4, 1 | 2) {
            let (s, r) = sqrt_rem_normalized(&(self << 2));
            let r = if s.bit(0) {
                (r + &s + &s - 1) >> 2
            } else {
                r >> 2
            };
            (s >> 1, r)
        } else {
            sqrt_rem_normalized(self)
        }
    }

    /// `⌊self^(1/n)⌋`.
    ///
    /// # Panics
    ///
    /// If `n` is zero.
    #[must_use]
    pub fn nth_root(&self, n: u32) -> Natural {
        match n {
            0 => panic!("zeroth root"),
            1 => self.clone(),
            2 => self.isqrt(),
            _ => root(self, n),
        }
    }

    /// `(s, r)` with `s = ⌊self^(1/n)⌋` and `r = self - s^n`.
    ///
    /// # Panics
    ///
    /// If `n` is zero.
    #[must_use]
    pub fn root_rem(&self, n: u32) -> (Natural, Natural) {
        if n == 2 {
            return self.sqrt_rem();
        }
        let s = self.nth_root(n);
        let r = self - s.pow(u64::from(n));
        (s, r)
    }

    /// Whether `self` is the square of an integer (0 and 1 are).
    ///
    /// Residues modulo 64 and modulo the odd primes up to 47 reject all but
    /// about one in 6000 non-squares before the root is computed.
    #[must_use]
    pub fn is_perfect_square(&self) -> bool {
        let low = self.limbs().first().copied().unwrap_or(0);
        if (SQUARES_MOD_64 >> (low % 64)) & 1 == 0 {
            return false;
        }
        let r = self.rem_limb(RESIDUE_PRODUCT);
        let residue = |(&p, &mask): (&Limb, &u64)| (mask >> (r % p)) & 1 == 1;
        RESIDUE_PRIMES.iter().zip(&SQUARE_MASKS).all(residue) && self.sqrt_rem().1.is_zero()
    }

    /// Whether `self = a^k` for some integers `a` and `k >= 2`. As in GMP, 0
    /// and 1 are perfect powers.
    #[must_use]
    pub fn is_perfect_power(&self) -> bool {
        self.limb_len() == 0 || self.is_one() || self.perfect_power().is_some()
    }

    /// `(a, k)` with `self = a^k` and `k >= 2` as large as possible (so `a`
    /// is not itself a perfect power), or `None` if there is none. 0 and 1
    /// give `None`: they have no largest exponent.
    ///
    /// Tries the prime exponents in increasing order, taking each root as
    /// long as it is exact. A prime `p` is first checked against the
    /// trailing zeros (which `p` must divide) and the `p`th power residues
    /// modulo a few primes `q ≡ 1 (mod p)`; a non-power passes those with
    /// probability about `p^-3`, so roots are rarely computed in vain.
    #[must_use]
    pub fn perfect_power(&self) -> Option<(Natural, u64)> {
        if self.limb_len() == 0 || self.is_one() {
            return None;
        }
        let (mut base, mut exp) = (self.clone(), 1);
        let mut p = 2;
        // `base = a^p` with `a >= 2` needs at least p + 1 bits.
        while p < base.bit_length() {
            while let Some(a) = exact_root(&base, p) {
                base = a;
                exp *= p;
            }
            p = next_prime(p);
        }
        (exp > 1).then_some((base, exp))
    }
}

/// `sqrt_rem` for a length of `4k - 1` or `4k` bits, above 128.
///
/// With `β = 2^k`, `m = a₃β³ + a₂β² + a₁β + a₀` and `a₃ >= β/4`: the root of
/// the top half `(s', r') = sqrt_rem(a₃β + a₂)` and the division
/// `(q, u) = divrem(r'β + a₁, 2s')` give `s = s'β + q` and
/// `r = uβ + a₀ - q²`, which is negative at most once, by a single step.
fn sqrt_rem_normalized(m: &Natural) -> (Natural, Natural) {
    let k = m.bit_length().div_ceil(4);
    debug_assert!(k > 32 && m.bit_length() >= 4 * k - 1);
    let (s1, r1) = (m >> (2 * k)).sqrt_rem();
    let a1 = low_bits(&(m >> k), k);
    let a0 = low_bits(m, k);
    let (q, u) = ((r1 << k) + a1).div_rem(&(&s1 << 1));
    let s = (s1 << k) + &q;
    let t = (u << k) + a0;
    let q2 = q.square();
    if let Some(r) = t.checked_sub(&q2) {
        return (s, r);
    }
    // r + 2s - 1 with the root s - 1.
    let r = (t + &s + &s - 1)
        .checked_sub(&q2)
        .expect("one correction suffices");
    (s - 1, r)
}

/// `x mod 2^k`.
fn low_bits(x: &Natural, k: u64) -> Natural {
    let (limb, bit) = split_bit_index(k);
    let a = x.limbs();
    if limb >= a.len() {
        return x.clone();
    }
    let mut v = a[..=limb].to_vec();
    v[limb] &= (1 << bit) - 1;
    Natural::from_vec(v)
}

/// `⌊x^(1/n)⌋` for `n >= 3`.
///
/// Newton's iteration `y = ((n - 1) r + ⌊x / r^(n-1)⌋) / n` decreases
/// strictly from any `r > ⌊x^(1/n)⌋` and does not go below it, so it stops
/// exactly at the root. The start is `(⌊(x / 2^(nh))^(1/n)⌋ + 1) 2^h`, an
/// overestimate with about half the bits of the root correct.
fn root(x: &Natural, n: u32) -> Natural {
    if let Ok(v) = u64::try_from(x) {
        return Natural::from(root_u64(v, n));
    }
    // The root has at most this many bits; x >= 1, so one bit means 1.
    let bits = x.bit_length().div_ceil(u64::from(n));
    if bits == 1 {
        return Natural::one();
    }
    let h = bits / 2;
    let mut r = (root(&(x >> (h * u64::from(n))), n) + 1) << h;
    let n_1 = u64::from(n - 1);
    loop {
        let y = (&r * n_1 + x / r.pow(n_1)) / u64::from(n);
        if y >= r {
            return r;
        }
        r = y;
    }
}

/// `⌊v^(1/n)⌋` for `n >= 2`: a floating-point estimate, corrected.
fn root_u64(v: u64, n: u32) -> u64 {
    if v < 2 || n >= u64::BITS {
        return v.min(1);
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "an estimate, corrected below"
    )]
    let mut r = (v as f64).powf(1.0 / f64::from(n)) as u64;
    let fits = |r: u64| r.checked_pow(n).is_some_and(|p| p <= v);
    while !fits(r) {
        r -= 1;
    }
    while fits(r + 1) {
        r += 1;
    }
    r
}

/// `a` with `x = a^p`, for `x >= 1` and a prime `p`, if there is one.
fn exact_root(x: &Natural, p: u64) -> Option<Natural> {
    if p >= x.bit_length() {
        // a >= 2 needs p + 1 bits, and x = 1 is excluded by the callers.
        return None;
    }
    let t = x.trailing_zeros().expect("x >= 1");
    if t % p != 0 {
        return None;
    }
    let odd = x >> t;
    let r = if odd.is_one() {
        odd
    } else {
        if p == 2 {
            if !odd.is_perfect_square() {
                return None;
            }
        } else if !power_residues(&odd, p) {
            return None;
        }
        let n = u32::try_from(p).ok()?;
        let (r, rem) = odd.root_rem(n);
        if !rem.is_zero() {
            return None;
        }
        r
    };
    Some(r << (t / p))
}

/// Whether `x` is a `p`th power residue modulo each of the first few primes
/// `q ≡ 1 (mod p)` (a necessary condition for a `p`th power), for an odd
/// prime `p`.
fn power_residues(x: &Natural, p: u64) -> bool {
    let mut q = 1;
    for _ in 0..FILTER_PRIMES {
        // q - 1 must be even, and q stays in range for any bit length.
        loop {
            q += 2 * p;
            if is_prime_u64(q) {
                break;
            }
        }
        let a = x.rem_limb(q);
        if a != 0 && ModWord::new(q).pow(a, (q - 1) / p) != 1 {
            return false;
        }
    }
    true
}

/// The smallest prime above `p`.
fn next_prime(p: u64) -> u64 {
    let mut q = p + 1;
    while !is_prime_u64(q) {
        q += 1;
    }
    q
}
