//! Pseudo-random numbers for the randomised algorithms: the splitting of
//! factors over finite fields, the evaluation points of Wang's algorithm and
//! of Zippel's sparse gcd. Their seeds are fixed, so that results repeat.

/// The `SplitMix64` generator \[SLF14\].
pub(crate) struct Random {
    state: u64,
}

impl Random {
    pub(crate) fn new(seed: u64) -> Random {
        Random { state: seed }
    }

    pub(crate) fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}
