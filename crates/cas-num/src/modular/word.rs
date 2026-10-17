//! Arithmetic modulo a single limb.

use crate::kernel::div::{div_2by1, reciprocal_2by1};
use crate::{Integer, Limb, Natural};

#[cfg(test)]
mod tests;

/// Arithmetic in ℤ/nℤ for a modulus `n` that fits a limb.
///
/// Elements are plain `u64` values in `0..n`, so they can be stored in
/// vectors and compared directly. Products are reduced by one division of a
/// double limb with a precomputed reciprocal (Möller–Granlund), which works
/// for every modulus, odd or even.
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

    /// Left-to-right binary exponentiation over the `bits` exponent bits
    /// given by `bit`.
    fn pow_bits(&self, a: Limb, bits: u64, bit: impl Fn(u64) -> bool) -> Limb {
        debug_assert!(a < self.n);
        let mut r = self.reduce(1);
        for i in (0..bits).rev() {
            r = self.mul(r, r);
            if bit(i) {
                r = self.mul(r, a);
            }
        }
        r
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
