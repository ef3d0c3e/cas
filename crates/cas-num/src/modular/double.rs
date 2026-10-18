//! Powers modulo a two-limb modulus: in Montgomery form on `u128` for an odd
//! one, and by the Chinese remainder theorem for an even one.
//!
//! [`super::ModRing`] serves any modulus, but its set-up (`R mod n` and
//! `R^2 mod n` by division) and heap buffers dominate a power modulo a
//! 65- to 128-bit number. This works on `u128` values with `R = 2^128`:
//! a product is four limb products, and REDC one more product by
//! `-n^(-1) mod 2^128` and the high half of one by `n`.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, e, m, n, o, s, t, x)"
)]

use super::ModWord;
use super::word::sliding_window;
use crate::Natural;

/// Montgomery arithmetic modulo an odd `n` in `(2^64, 2^128)`. Elements in
/// Montgomery form are `x R mod n`; addition, subtraction, negation and
/// halving are the same in either form.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ModDouble {
    n: u128,
    /// `-n^(-1) mod 2^128`.
    ninv: u128,
    /// `R mod n`, the element 1.
    one: u128,
    /// `R^2 mod n`, for conversion into Montgomery form.
    r2: u128,
}

/// `(low, high)` halves of the 256-bit `a b`.
#[inline]
#[expect(clippy::cast_possible_truncation, reason = "splitting double limbs")]
fn mul_wide(a: u128, b: u128) -> (u128, u128) {
    let (a0, a1) = (a as u64, (a >> 64) as u64);
    let (b0, b1) = (b as u64, (b >> 64) as u64);
    let p00 = u128::from(a0) * u128::from(b0);
    let p01 = u128::from(a0) * u128::from(b1);
    let p10 = u128::from(a1) * u128::from(b0);
    let p11 = u128::from(a1) * u128::from(b1);
    // The middle column: three terms below 2^64, so no overflow.
    let mid = (p00 >> 64) + u128::from(p01 as u64) + u128::from(p10 as u64);
    let low = (p00 & u128::from(u64::MAX)) | (mid << 64);
    let high = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    (low, high)
}

impl ModDouble {
    /// The arithmetic modulo `n`, for an odd `n` above `2^64`.
    pub(crate) fn new(n: u128) -> Self {
        debug_assert!(n & 1 == 1 && n >> 64 != 0);
        // n n = 1 mod 8: 3 correct bits, doubling with each Newton step.
        let mut inv = n;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u128.wrapping_sub(n.wrapping_mul(inv)));
        }
        debug_assert_eq!(n.wrapping_mul(inv), 1);
        // R mod n = (R - n) mod n, and R^2 = R mod n doubled 128 times.
        let one = 0u128.wrapping_sub(n) % n;
        let mut d = ModDouble {
            n,
            ninv: inv.wrapping_neg(),
            one,
            r2: one,
        };
        for _ in 0..128 {
            d.r2 = d.double(d.r2);
        }
        d
    }

    /// The Montgomery form of `a < n`.
    pub(crate) fn to_mont(self, a: u128) -> u128 {
        self.mul(a, self.r2)
    }

    /// The element 1, in Montgomery form.
    pub(crate) fn one(self) -> u128 {
        self.one
    }

    /// `x + y mod n`.
    pub(crate) fn add(self, x: u128, y: u128) -> u128 {
        let (s, over) = x.overflowing_add(y);
        if over || s >= self.n {
            s.wrapping_sub(self.n)
        } else {
            s
        }
    }

    /// `x - y mod n`.
    pub(crate) fn sub(self, x: u128, y: u128) -> u128 {
        let (d, under) = x.overflowing_sub(y);
        if under { d.wrapping_add(self.n) } else { d }
    }

    /// `-x mod n`.
    pub(crate) fn neg(self, x: u128) -> u128 {
        if x == 0 { 0 } else { self.n - x }
    }

    /// `x / 2 mod n`: `(x + n) / 2` for odd `x`, taken as `x / 2 + n / 2 + 1`
    /// (both odd) so that it cannot overflow.
    pub(crate) fn half(self, x: u128) -> u128 {
        if x & 1 == 0 {
            x >> 1
        } else {
            (x >> 1) + (self.n >> 1) + 1
        }
    }

    /// `x^e` in Montgomery form, with `x^0 = 1`.
    pub(crate) fn pow_mont(self, x: u128, e: &Natural) -> u128 {
        let bits = e.bit_length();
        if bits == 0 {
            return self.one;
        }
        sliding_window(bits, &|i| e.bit(i), x, |p, q| self.mul(p, q))
    }

    /// `x y R^(-1) mod n` for `x, y < n`: with `t = x y` and
    /// `m = t ninv mod R`, `t + m n` is divisible by `R`; its low halves
    /// cancel, carrying exactly when `t`'s is non-zero, and the high halves add
    /// up to less than `2n`.
    #[inline]
    pub(crate) fn mul(self, x: u128, y: u128) -> u128 {
        let (tl, th) = mul_wide(x, y);
        let m = tl.wrapping_mul(self.ninv);
        let (_, mh) = mul_wide(m, self.n);
        let (u, over) = th.overflowing_add(mh + u128::from(tl != 0));
        if over || u >= self.n {
            u.wrapping_sub(self.n)
        } else {
            u
        }
    }

    /// `2x mod n` for `x < n`.
    #[inline]
    fn double(self, x: u128) -> u128 {
        let (s, over) = x.overflowing_add(x);
        if over || s >= self.n {
            s.wrapping_sub(self.n)
        } else {
            s
        }
    }

    /// `a^e mod n` for `a < n`, with `a^0 = 1`: sliding windows in
    /// Montgomery form, and REDC of the result by 1 to leave it.
    pub(crate) fn pow(self, a: u128, e: &Natural) -> u128 {
        debug_assert!(a < self.n);
        self.mul(self.pow_mont(self.to_mont(a), e), 1)
    }
}

/// `a^e mod m` for an even `m` in `(2^64, 2^128)` and `a < m`, with
/// `a^0 = 1`: for `m = 2^s o`, `o` odd, the power modulo `2^s` by wrapping
/// arithmetic, the power modulo `o` by [`ModDouble`] or [`ModWord`], and the
/// two combined, `x = x_o + o ((x_2 - x_o) o^(-1) mod 2^s)`, below `m`.
pub(crate) fn pow_mod_even(a: u128, e: &Natural, m: u128) -> u128 {
    debug_assert!(m & 1 == 0 && m >> 64 != 0 && a < m);
    let s = m.trailing_zeros();
    let o = m >> s;
    let mask = (1u128 << s) - 1;
    let bits = e.bit_length();
    let mut x2 = 1u128;
    for i in (0..bits).rev() {
        x2 = x2.wrapping_mul(x2);
        if e.bit(i) {
            x2 = x2.wrapping_mul(a);
        }
    }
    let x2 = x2 & mask;
    let xo = if o == 1 {
        0
    } else if let Ok(word) = u64::try_from(o) {
        let w = ModWord::new(word);
        #[expect(clippy::cast_possible_truncation, reason = "reduced below o")]
        let r = (a % o) as u64;
        u128::from(w.pow_natural(r, e))
    } else {
        ModDouble::new(o).pow(a % o, e)
    };
    // o^(-1) mod 2^128 by Newton's iteration; its low s bits invert o mod 2^s.
    let mut inv = o;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u128.wrapping_sub(o.wrapping_mul(inv)));
    }
    let t = x2.wrapping_sub(xo).wrapping_mul(inv) & mask;
    xo + o * t
}

#[cfg(test)]
mod tests;
