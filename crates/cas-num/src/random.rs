//! Pseudo-random numbers (CAS-DOC-33 §1, CAS-370): the generator
//! xoshiro256\*\* \[BV18\], seeded through `SplitMix64` \[SLF14\] from an
//! integer of any size, and the natural numbers it draws uniformly.
//!
//! - **The state**: 256 bits, never all zero, of period `2^256 − 1`.
//! - **A seed**: `SplitMix64`'s output function absorbs its sign, its count
//!   of limbs and its limbs one word at a time, and `SplitMix64` from the
//!   word it ends with gives the four words of the state, which are never
//!   all zero: its outputs are a bijection of a counter, 0 only at 0.
//! - **Natural numbers**: `k` bits are whole words, the top one shifted
//!   down; a number below `n` is drawn by rejection, the bits of `n` until
//!   they fall below it, at most 2 tries on average, exact and unbiased.
//!
//! \[BV18\] D. Blackman, S. Vigna, *Scrambled linear pseudorandom number
//! generators*, 2018. \[SLF14\] G. L. Steele, D. Lea, C. H. Flood, *Fast
//! splittable pseudorandom number generators*, OOPSLA 2014.

use crate::{Integer, Limb, Natural};

/// `SplitMix64`'s increment, the golden ratio's fraction.
const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;

/// `SplitMix64`'s output function, a bijection of the words.
fn mix(z: u64) -> u64 {
    let z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    let z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// The generator xoshiro256\*\*.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Generator {
    s: [u64; 4],
}

impl Generator {
    /// The generator of `seed`: one sequence for each integer.
    #[must_use]
    pub fn seeded(seed: &Integer) -> Generator {
        let magnitude = seed.unsigned_abs();
        let limbs = magnitude.limbs();
        let mut h = 0u64;
        let count = u64::try_from(limbs.len()).unwrap_or(u64::MAX);
        let words = [u64::from(seed.is_negative()), count];
        for &w in words.iter().chain(limbs) {
            h = mix(h.wrapping_add(GOLDEN) ^ w);
        }
        let mut s = [0u64; 4];
        for word in &mut s {
            h = h.wrapping_add(GOLDEN);
            *word = mix(h);
        }
        debug_assert!(s != [0; 4], "SplitMix64 gave four zeros");
        Generator { s }
    }

    /// The generator at the state `s`, which is not all zero.
    #[cfg(test)]
    pub(crate) fn from_state(s: [u64; 4]) -> Generator {
        assert!(s != [0; 4], "Generator::from_state: the zero state");
        Generator { s }
    }

    /// The next 64 bits.
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// A natural number uniform below `2^k`.
    pub fn bits(&mut self, k: u64) -> Natural {
        let words = k.div_ceil(Limb::BITS.into());
        let mut limbs: Vec<Limb> = (0..words).map(|_| self.next_u64()).collect();
        let spare = words * u64::from(Limb::BITS) - k;
        if let Some(top) = limbs.last_mut() {
            *top >>= spare;
        }
        Natural::from_vec(limbs)
    }

    /// A natural number uniform below `n`.
    ///
    /// # Panics
    ///
    /// If `n` is 0.
    pub fn below(&mut self, n: &Natural) -> Natural {
        assert!(!n.is_zero(), "Generator::below: no number below 0");
        let k = n.bit_length();
        loop {
            let r = self.bits(k);
            if r < *n {
                return r;
            }
        }
    }
}

#[cfg(test)]
mod tests;
