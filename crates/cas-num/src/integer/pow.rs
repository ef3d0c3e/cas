//! Powers, modular powers, roots and perfect powers.

use super::{Integer, Repr};

#[cfg(test)]
mod tests;

const NEGATIVE_SQRT: &str = "argument of integer square root cannot be negative";
const NEGATIVE_EVEN_ROOT: &str = "even root of a negative integer";

impl Integer {
    /// `self^e`, with `0^0 = 1`.
    ///
    /// # Panics
    ///
    /// If the bit length of the result overflows `u64`. (Like every
    /// allocation, running out of memory aborts.)
    #[must_use]
    pub fn pow(&self, e: u64) -> Integer {
        if let Repr::Small(x) = self.0 {
            if let Some(p) = u32::try_from(e).ok().and_then(|e| x.checked_pow(e)) {
                return Integer::from(p);
            }
        }
        let (negative, mag) = self.sign_mag();
        Integer::from_parts(negative && e % 2 == 1, mag.pow(e))
    }

    /// `self^e mod |m|`, in `0..|m|`. A negative exponent raises the inverse
    /// of `self`; `None` if `m` is zero, or if `e` is negative and `self` has
    /// no inverse modulo `m`.
    #[must_use]
    pub fn pow_mod(&self, e: &Integer, m: &Integer) -> Option<Integer> {
        if m.is_zero() {
            return None;
        }
        let base = if e.is_negative() {
            self.mod_inverse(m)?
        } else {
            self.div_rem_euclid(m).1
        };
        // In 0..|m|, so the sign is positive.
        let (_, base) = base.into_parts();
        let (_, m) = m.sign_mag();
        Some(Integer::from(base.pow_mod(&e.sign_mag().1, &m)))
    }

    /// `⌊√self⌋`.
    ///
    /// # Panics
    ///
    /// If `self` is negative.
    #[must_use]
    pub fn isqrt(&self) -> Integer {
        assert!(!self.is_negative(), "{NEGATIVE_SQRT}");
        Integer::from(self.sign_mag().1.isqrt())
    }

    /// `(s, r)` with `s = ⌊√self⌋` and `r = self - s²`.
    ///
    /// # Panics
    ///
    /// If `self` is negative.
    #[must_use]
    pub fn sqrt_rem(&self) -> (Integer, Integer) {
        assert!(!self.is_negative(), "{NEGATIVE_SQRT}");
        let (s, r) = self.sign_mag().1.sqrt_rem();
        (Integer::from(s), Integer::from(r))
    }

    /// The `n`th root, truncated towards zero (as GMP's `mpz_root`).
    ///
    /// # Panics
    ///
    /// If `n` is zero, or if `n` is even and `self` is negative.
    #[must_use]
    pub fn nth_root(&self, n: u32) -> Integer {
        let (negative, mag) = self.sign_mag();
        assert!(!negative || n % 2 == 1, "{NEGATIVE_EVEN_ROOT}");
        Integer::from_parts(negative, mag.nth_root(n))
    }

    /// `(s, r)` with `s` the `n`th root truncated towards zero and
    /// `r = self - s^n`, which has the sign of `self`.
    ///
    /// # Panics
    ///
    /// If `n` is zero, or if `n` is even and `self` is negative.
    #[must_use]
    pub fn root_rem(&self, n: u32) -> (Integer, Integer) {
        let (negative, mag) = self.sign_mag();
        assert!(!negative || n % 2 == 1, "{NEGATIVE_EVEN_ROOT}");
        let (s, r) = mag.root_rem(n);
        (
            Integer::from_parts(negative, s),
            Integer::from_parts(negative, r),
        )
    }

    /// Whether `self` is the square of an integer; never for a negative.
    #[must_use]
    pub fn is_perfect_square(&self) -> bool {
        !self.is_negative() && self.sign_mag().1.is_perfect_square()
    }

    /// Whether `self = a^k` for some integers `a` and `k >= 2` (as GMP's
    /// `mpz_perfect_power_p`): 0, 1 and -1 are, and a negative value only
    /// with an odd `k`.
    #[must_use]
    pub fn is_perfect_power(&self) -> bool {
        let (negative, mag) = self.sign_mag();
        if !negative || mag.is_one() {
            return mag.is_perfect_power();
        }
        // |self| = a^k with k largest is an odd power iff k has an odd factor.
        mag.perfect_power()
            .is_some_and(|(_, k)| k >> k.trailing_zeros() > 1)
    }
}
