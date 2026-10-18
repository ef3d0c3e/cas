//! Arithmetic modulo a single limb.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, m, p, q, t, u, x)"
)]

use super::ring::{limb_inverse, window_bits};
use crate::kernel::div::{div_2by1, reciprocal_2by1};
use crate::{Integer, Limb, Natural};

#[cfg(test)]
mod tests;

/// Arithmetic in ℤ/nℤ for a modulus `n` that fits a limb.
///
/// Elements are plain `u64` values in `0..n`, so they can be stored in
/// vectors and compared directly. Products are reduced by one division of a
/// double limb with a precomputed reciprocal (Möller–Granlund), which works
/// for every modulus, odd or even. Powers modulo an odd `n` run in Montgomery
/// form inside [`ModWord::pow`], with sliding windows.
///
/// The arithmetic methods require reduced operands (in `0..n`), which is
/// checked in debug builds only; otherwise the result is unspecified. The
/// `reduce*` methods bring any value into range.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct ModWord {
    n: Limb,
    /// `n << shift` has its top bit set.
    shift: u32,
    /// `reciprocal_2by1(n << shift)`.
    v: Limb,
    /// `-n^(-1) mod B` for an odd `n`, else 0.
    ninv: Limb,
}

impl ModWord {
    /// Arithmetic modulo `n`.
    ///
    /// # Panics
    ///
    /// If `n` is zero.
    #[must_use]
    pub fn new(n: Limb) -> ModWord {
        assert!(n != 0, "zero modulus");
        let shift = n.leading_zeros();
        ModWord {
            n,
            shift,
            v: reciprocal_2by1(n << shift),
            ninv: if n & 1 == 1 {
                limb_inverse(n).wrapping_neg()
            } else {
                0
            },
        }
    }

    /// The modulus `n`.
    #[must_use]
    pub const fn modulus(&self) -> Limb {
        self.n
    }

    /// `(hi B + lo) mod n` for `hi < n`.
    #[inline]
    fn rem_2by1(&self, hi: Limb, lo: Limb) -> Limb {
        debug_assert!(hi < self.n);
        let s = self.shift;
        // hi < n, so the shifted top limb stays below n << s.
        let (u1, u0) = if s == 0 {
            (hi, lo)
        } else {
            ((hi << s) | (lo >> (Limb::BITS - s)), lo << s)
        };
        div_2by1(u1, u0, self.n << s, self.v).1 >> s
    }

    /// `a mod n`.
    #[must_use]
    #[inline]
    pub fn reduce(&self, a: Limb) -> Limb {
        self.rem_2by1(0, a)
    }

    /// `a mod n` for a double limb.
    #[must_use]
    #[expect(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    pub fn reduce_u128(&self, a: u128) -> Limb {
        let hi = self.reduce((a >> Limb::BITS) as Limb);
        self.rem_2by1(hi, a as Limb)
    }

    /// `a mod n`.
    #[must_use]
    pub fn reduce_natural(&self, a: &Natural) -> Limb {
        a.rem_limb(self.n)
    }

    /// `a mod n`, in `0..n` also for negative `a`.
    #[must_use]
    pub fn reduce_integer(&self, a: &Integer) -> Limb {
        let r = match a.to_i64() {
            Some(x) => self.reduce(x.unsigned_abs()),
            None => self.reduce_natural(&a.unsigned_abs()),
        };
        if a.is_negative() { self.neg(r) } else { r }
    }

    /// `a + b mod n`.
    #[must_use]
    #[inline]
    pub fn add(&self, a: Limb, b: Limb) -> Limb {
        debug_assert!(a < self.n && b < self.n);
        let (s, carry) = a.overflowing_add(b);
        if carry || s >= self.n {
            s.wrapping_sub(self.n)
        } else {
            s
        }
    }

    /// `a - b mod n`.
    #[must_use]
    #[inline]
    pub fn sub(&self, a: Limb, b: Limb) -> Limb {
        debug_assert!(a < self.n && b < self.n);
        if a >= b {
            a - b
        } else {
            a.wrapping_sub(b).wrapping_add(self.n)
        }
    }

    /// `-a mod n`.
    #[must_use]
    #[inline]
    pub fn neg(&self, a: Limb) -> Limb {
        debug_assert!(a < self.n);
        if a == 0 { 0 } else { self.n - a }
    }

    /// `a b mod n`.
    #[must_use]
    #[inline]
    #[expect(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    pub fn mul(&self, a: Limb, b: Limb) -> Limb {
        debug_assert!(a < self.n && b < self.n);
        // a, b < n, so the high limb of the product is below n.
        let t = u128::from(a) * u128::from(b);
        self.rem_2by1((t >> Limb::BITS) as Limb, t as Limb)
    }

    /// `a^e mod n`, with `a^0 = 1 mod n`.
    #[must_use]
    pub fn pow(&self, a: Limb, e: u64) -> Limb {
        self.pow_bits(a, u64::from(u64::BITS - e.leading_zeros()), |i| {
            (e >> i) & 1 == 1
        })
    }

    /// `a^e mod n` for a multi-precision exponent.
    pub(crate) fn pow_natural(&self, a: Limb, e: &Natural) -> Limb {
        self.pow_bits(a, e.bit_length(), |i| e.bit(i))
    }

    /// `a^e` over the `bits` exponent bits given by `bit`: in Montgomery
    /// form with sliding windows for an odd modulus (from `a R` to `a^e R`,
    /// `R = B`, then back), else by left-to-right binary exponentiation.
    fn pow_bits(&self, a: Limb, bits: u64, bit: impl Fn(u64) -> bool) -> Limb {
        debug_assert!(a < self.n);
        if self.ninv == 0 || bits == 0 {
            let mut r = self.reduce(1);
            for i in (0..bits).rev() {
                r = self.mul(r, r);
                if bit(i) {
                    r = self.mul(r, a);
                }
            }
            return r;
        }
        let x = self.rem_2by1(a, 0);
        let r = sliding_window(bits, &bit, x, |p, q| self.mont_mul(p, q));
        self.mont_mul(r, 1)
    }

    /// `p q R^(-1) mod n` for an odd `n` and `p, q < n`: REDC of the double
    /// limb `t = p q`, with `m = t ninv mod B` making `t + m n` divisible by
    /// `B`. The low limbs cancel, carrying exactly when `t`'s is non-zero; the
    /// high limbs add up to less than `2n`.
    #[inline]
    #[expect(clippy::cast_possible_truncation, reason = "splitting double limbs")]
    fn mont_mul(&self, p: Limb, q: Limb) -> Limb {
        let t = u128::from(p) * u128::from(q);
        let m = (t as Limb).wrapping_mul(self.ninv);
        let mn = u128::from(m) * u128::from(self.n);
        let carry = Limb::from(t as Limb != 0);
        let (u, over) =
            ((t >> Limb::BITS) as Limb).overflowing_add((mn >> Limb::BITS) as Limb + carry);
        if over || u >= self.n {
            u.wrapping_sub(self.n)
        } else {
            u
        }
    }

    /// The inverse of `a` modulo `n`, or `None` if they are not coprime.
    /// Modulo 1 the inverse of 0 is 0.
    #[must_use]
    pub fn inv(&self, a: Limb) -> Option<Limb> {
        debug_assert!(a < self.n);
        // Extended Euclid on (n, a), tracking the coefficient of a; the
        // coefficients are bounded by n in absolute value.
        let (mut r0, mut r1) = (self.n, a);
        let (mut t0, mut t1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (t0, t1) = (t1, t0 - i128::from(q) * t1);
        }
        if r0 != 1 {
            return None;
        }
        #[expect(clippy::cast_possible_truncation, reason = "reduced into 0..n")]
        let t = t0.rem_euclid(i128::from(self.n)) as Limb;
        Some(t)
    }
}

/// `x^e` for the `bits`-bit exponent `e` given by `bit` (top bit set), by
/// left-to-right sliding windows of up to [`window_bits`] bits ending in a
/// one, with `mul` the product: the odd powers `x, x^3, …` are tabled, and
/// each window costs one product after its squarings.
pub(super) fn sliding_window<T: Copy + Default>(
    bits: u64,
    bit: &impl Fn(u64) -> bool,
    x: T,
    mul: impl Fn(T, T) -> T,
) -> T {
    // Width 8 already suits exponents of 4609 to 11520 bits.
    let w = window_bits(bits).min(8);
    // table[i] = x^(2i + 1).
    let mut table = [T::default(); 1 << 7];
    table[0] = x;
    let count = 1usize << (w - 1);
    if count > 1 {
        let x2 = mul(x, x);
        for i in 1..count {
            table[i] = mul(table[i - 1], x2);
        }
    }
    // The window of bits i - 1 down to l, at most w long and ending in a one.
    let window = |i: u64| {
        let mut l = i.saturating_sub(w);
        while !bit(l) {
            l += 1;
        }
        let v = (l..i)
            .rev()
            .fold(0usize, |v, j| (v << 1) | usize::from(bit(j)));
        (l, v)
    };
    let (mut i, v) = window(bits);
    let mut r = table[v >> 1];
    while i > 0 {
        if !bit(i - 1) {
            r = mul(r, r);
            i -= 1;
            continue;
        }
        let (l, v) = window(i);
        for _ in l..i {
            r = mul(r, r);
        }
        r = mul(r, table[v >> 1]);
        i = l;
    }
    r
}
