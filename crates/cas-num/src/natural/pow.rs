//! Powers and modular powers.

use super::Natural;
use crate::modular::{ModDouble, ModFixed, ModRing, ModWord, pow_mod_even};

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
    /// A single-limb modulus uses [`ModWord`]; a two-limb one `u128`
    /// arithmetic (Montgomery form when odd, else split by the Chinese
    /// remainder theorem); an odd one of three to eight limbs Montgomery
    /// arithmetic on stack arrays; any other [`ModRing`] (Montgomery
    /// multiplication for odd moduli, sliding-window exponentiation).
    ///
    /// # Panics
    ///
    /// If `m` is zero.
    #[must_use]
    pub fn pow_mod(&self, e: &Natural, m: &Natural) -> Natural {
        assert!(!m.is_zero(), "attempt to divide by zero");
        if let [n] = *m.limbs() {
            let w = ModWord::new(n);
            return Natural::from(w.pow_natural(w.reduce_natural(self), e));
        }
        if let [lo, hi] = *m.limbs() {
            let n = u128::from(lo) | u128::from(hi) << 64;
            let a = self % m;
            let a = a
                .limbs()
                .iter()
                .rev()
                .fold(0u128, |v, &l| v << 64 | u128::from(l));
            return Natural::from(if lo & 1 == 1 {
                ModDouble::new(n).pow(a, e)
            } else {
                pow_mod_even(a, e, n)
            });
        }
        if m.bit(0) {
            match m.limb_len() {
                3 => return fixed_pow::<3>(self, e, m),
                4 => return fixed_pow::<4>(self, e, m),
                5 => return fixed_pow::<5>(self, e, m),
                6 => return fixed_pow::<6>(self, e, m),
                7 => return fixed_pow::<7>(self, e, m),
                8 => return fixed_pow::<8>(self, e, m),
                _ => {}
            }
        }
        let ring = ModRing::new(m.clone());
        ring.to_natural(&ring.pow(&ring.elem(self), e))
    }
}

/// `a^e mod m` for an odd `K`-limb `m`, on stack arrays.
fn fixed_pow<const K: usize>(a: &Natural, e: &Natural, m: &Natural) -> Natural {
    let ring = ModFixed::<K>::new(m);
    ring.value(ring.pow(ring.elem(a), e))
}
