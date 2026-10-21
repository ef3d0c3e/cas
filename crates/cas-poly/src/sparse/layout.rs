//! Packed exponent vectors \[MP07\]: the exponents of a term in fields of
//! machine words, so that comparing two monomials compares a few words, and
//! multiplying them adds words.

/// How the exponents of a term are packed: `nvars` fields of `bits` bits,
/// `per_word` of them in each 64-bit word, the first variable in the top
/// field of the first word. Comparing the words in order then compares the
/// exponent vectors lexicographically.
///
/// The top bit of each field is a guard \[MP07\], clear in a packed vector:
/// a sum of two vectors that sets one has overflowed that field, and a
/// difference computed with the guards set that clears one has gone
/// negative, so the subtrahend does not divide.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) struct Layout {
    nvars: usize,
    bits: u32,
    per_word: usize,
    words: usize,
    /// The guards of a word full of fields, and of the last word.
    guards: u64,
    last_guards: u64,
}

impl Layout {
    /// The narrowest fields: 8 bits, exponents up to 127.
    pub(crate) const MIN_BITS: u32 = 8;

    /// Fields of `bits` bits, guard included.
    ///
    /// # Panics
    ///
    /// If `bits` is not in `2..=64`.
    pub(crate) fn new(nvars: usize, bits: u32) -> Layout {
        assert!((2..=64).contains(&bits), "fields of {bits} bits");
        let per_word = (64 / bits) as usize;
        let words = nvars.div_ceil(per_word);
        let guards_of = |fields: usize| {
            (0..fields).fold(0u64, |g, f| {
                let shift = 64 - bits * (u32::try_from(f).expect("fits") + 1);
                g | 1 << (shift + bits - 1)
            })
        };
        let in_last = nvars - words.saturating_sub(1) * per_word;
        Layout {
            nvars,
            bits,
            per_word,
            words,
            guards: guards_of(per_word),
            last_guards: guards_of(in_last),
        }
    }

    /// The narrowest layout, from [`Layout::MIN_BITS`], whose fields hold
    /// `max`.
    ///
    /// # Panics
    ///
    /// If `max` is `2^63` or more.
    pub(crate) fn holding(nvars: usize, max: u64) -> Layout {
        let bits = (64 - max.leading_zeros() + 1).max(Layout::MIN_BITS);
        assert!(bits <= 64, "exponents below 2^63");
        Layout::new(nvars, bits)
    }

    pub(crate) fn nvars(&self) -> usize {
        self.nvars
    }

    pub(crate) fn bits(&self) -> u32 {
        self.bits
    }

    /// The words of a packed vector.
    #[inline]
    pub(crate) fn words(&self) -> usize {
        self.words
    }

    /// The largest exponent a field holds.
    pub(crate) fn max_exponent(&self) -> u64 {
        (1u64 << (self.bits - 1)) - 1
    }

    fn mask(&self) -> u64 {
        if self.bits == 64 {
            u64::MAX
        } else {
            (1u64 << self.bits) - 1
        }
    }

    /// The word and the shift of `var`'s field.
    fn position(&self, var: usize) -> (usize, u32) {
        let field = u32::try_from(var % self.per_word).expect("fits");
        (var / self.per_word, 64 - self.bits * (field + 1))
    }

    fn guards_of(&self, word: usize) -> u64 {
        if word + 1 == self.words {
            self.last_guards
        } else {
            self.guards
        }
    }

    /// The exponent of `var` in `mono`.
    pub(crate) fn get(&self, mono: &[u64], var: usize) -> u64 {
        let (w, shift) = self.position(var);
        (mono[w] >> shift) & self.mask()
    }

    /// Sets the exponent of `var` in `mono` to `e`, at most
    /// [`Layout::max_exponent`].
    pub(crate) fn set(&self, mono: &mut [u64], var: usize, e: u64) {
        debug_assert!(e <= self.max_exponent());
        let (w, shift) = self.position(var);
        mono[w] = (mono[w] & !(self.mask() << shift)) | e << shift;
    }

    /// `exps` packed into `out`, each at most [`Layout::max_exponent`].
    pub(crate) fn pack(&self, exps: &[u64], out: &mut [u64]) {
        debug_assert_eq!(exps.len(), self.nvars);
        out.fill(0);
        for (var, &e) in exps.iter().enumerate() {
            debug_assert!(e <= self.max_exponent(), "{e} in {} bits", self.bits);
            let (w, shift) = self.position(var);
            out[w] |= e << shift;
        }
    }

    /// The exponents of `mono`.
    pub(crate) fn unpack(&self, mono: &[u64]) -> Vec<u64> {
        (0..self.nvars).map(|var| self.get(mono, var)).collect()
    }

    /// The largest exponent in `mono`.
    pub(crate) fn max_in(&self, mono: &[u64]) -> u64 {
        (0..self.nvars)
            .map(|var| self.get(mono, var))
            .max()
            .unwrap_or(0)
    }

    /// `out = a + b`: the product of the monomials. Whether a field
    /// overflowed, which sets its guard.
    #[inline]
    pub(crate) fn add(&self, a: &[u64], b: &[u64], out: &mut [u64]) -> bool {
        let mut overflow = false;
        for (w, ((o, x), y)) in out.iter_mut().zip(a).zip(b).enumerate() {
            // Fields below 2^(bits-1) sum below 2^bits: no carry leaves one.
            *o = x + y;
            overflow |= *o & self.guards_of(w) != 0;
        }
        overflow
    }

    /// `out = b - a` when `a` divides `b`, every field of `a` at most `b`'s:
    /// with the guards of `b` set, a field of `a` above `b`'s clears its
    /// guard and borrows from nothing else \[MP07\].
    #[inline]
    pub(crate) fn divides(&self, a: &[u64], b: &[u64], out: &mut [u64]) -> bool {
        for (w, ((o, x), y)) in out.iter_mut().zip(a).zip(b).enumerate() {
            let g = self.guards_of(w);
            let d = (y | g) - x;
            if d & g != g {
                return false;
            }
            *o = d & !g;
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packing_round_trips_and_orders() {
        for bits in [8, 9, 13, 16, 21, 32, 33, 64] {
            for nvars in [0, 1, 2, 3, 7, 8, 9, 17] {
                let layout = Layout::new(nvars, bits);
                let top = layout.max_exponent();
                let exps: Vec<u64> = (0..nvars).map(|v| (v as u64 * 7919) % (top + 1)).collect();
                let mut packed = vec![0; layout.words()];
                layout.pack(&exps, &mut packed);
                assert_eq!(layout.unpack(&packed), exps, "{bits} bits, {nvars} vars");
                // Lexicographic order is word order.
                if nvars >= 2 {
                    let mut bigger = exps.clone();
                    bigger[nvars - 1] = (bigger[nvars - 1] + 1).min(top);
                    bigger[0] = (bigger[0] + 1) % (top + 1);
                    let mut other = vec![0; layout.words()];
                    layout.pack(&bigger, &mut other);
                    assert_eq!(packed.cmp(&other), exps.cmp(&bigger), "{bits} bits");
                }
            }
        }
    }

    #[test]
    fn guards_catch_overflow_and_non_divisors() {
        let layout = Layout::new(3, 8);
        let pack = |e: &[u64]| {
            let mut p = vec![0; layout.words()];
            layout.pack(e, &mut p);
            p
        };
        let mut out = vec![0; layout.words()];
        // 127 is the largest exponent in 8 bits: 127 + 0 fits, 127 + 1 does
        // not, in any field.
        assert!(!layout.add(&pack(&[127, 3, 0]), &pack(&[0, 124, 5]), &mut out));
        assert_eq!(layout.unpack(&out), [127, 127, 5]);
        for var in 0..3 {
            let mut e = [0; 3];
            e[var] = 127;
            let mut one = [0; 3];
            one[var] = 1;
            assert!(layout.add(&pack(&e), &pack(&one), &mut out), "field {var}");
        }
        // x y^2 divides x^3 y^2 z, not x^3 y z.
        assert!(layout.divides(&pack(&[1, 2, 0]), &pack(&[3, 2, 1]), &mut out));
        assert_eq!(layout.unpack(&out), [2, 0, 1]);
        assert!(!layout.divides(&pack(&[1, 2, 0]), &pack(&[3, 1, 1]), &mut out));
        assert!(!layout.divides(&pack(&[0, 0, 127]), &pack(&[127, 127, 126]), &mut out));
        assert_eq!(Layout::holding(2, 127).bits(), 8);
        assert_eq!(Layout::holding(2, 128).bits(), 9);
        assert_eq!(Layout::holding(2, (1 << 62) - 1).bits(), 63);
        assert_eq!(Layout::holding(2, 1 << 62).bits(), 64);
    }
}
