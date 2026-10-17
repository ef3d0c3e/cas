//! Powers and modular powers.

use super::Natural;

#[cfg(test)]
mod tests;

impl Natural {
    /// `self^e`, with `0^0 = 1`.
    ///
    /// Factors of two are split off and restored with one shift; the odd part
    /// is raised by left-to-right binary exponentiation.
    ///
    /// # Panics
    ///
    /// If the bit length of the result overflows `u64`. (Like every
    /// allocation, running out of memory aborts.)
    #[must_use]
    pub fn pow(&self, e: u64) -> Natural {
        if e == 0 {
            return Natural::one();
        }
        let Some(t) = self.trailing_zeros() else {
            return Natural::ZERO;
        };
        let shift = t.checked_mul(e).expect("power too large");
        let odd = self >> t;
        if odd.is_one() {
            return Natural::power_of_two(shift);
        }
        // The result has at least this many bits.
        let bits = (odd.bit_length() - 1)
            .checked_mul(e)
            .and_then(|b| b.checked_add(shift))
            .and_then(|b| b.checked_add(1));
        assert!(bits.is_some(), "power too large");
        let mut r = odd.clone();
        for i in (0..u64::BITS - 1 - e.leading_zeros()).rev() {
            r = r.square();
            if (e >> i) & 1 == 1 {
                r = match *odd.limbs() {
                    [b] => r * b,
                    _ => &r * &odd,
                };
            }
        }
        r << shift
    }

    /// `self^e mod m`, in `0..m`.
    ///
    /// Left-to-right binary exponentiation with a full division after every
    /// product; the modular arithmetic ticket (CAS-14) replaces the reduction.
    ///
    /// # Panics
    ///
    /// If `m` is zero.
    #[must_use]
    pub fn pow_mod(&self, e: &Natural, m: &Natural) -> Natural {
        assert!(!m.is_zero(), "attempt to divide by zero");
        if m.is_one() {
            return Natural::ZERO;
        }
        let b = self % m;
        let mut r = Natural::one();
        for i in (0..e.bit_length()).rev() {
            r = r.square() % m;
            if e.bit(i) {
                r = &r * &b % m;
            }
        }
        r
    }
}
