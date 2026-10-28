use proptest::prelude::*;

use super::Generator;
use crate::{Integer, Natural};

/// The reference implementation's first outputs from the state `{1, 2, 3,
/// 4}`.
#[test]
fn reference_outputs() {
    let mut g = Generator::from_state([1, 2, 3, 4]);
    let outputs: Vec<u64> = (0..6).map(|_| g.next_u64()).collect();
    assert_eq!(
        outputs,
        [
            0x2d00,
            0,
            0x5a00_7080,
            0x10e0_0000_0000_9d80,
            0x10e0_b61c_e100_9d80,
            0x0870_021c_e143_ad00,
        ]
    );
}

/// One sequence for each seed: the same seed repeats it, and the seeds 0,
/// 1, −1, `2^64` and `−2^64` start five different ones.
#[test]
fn seeds() {
    let first = |seed: &Integer| Generator::seeded(seed).next_u64();
    let big = Integer::from(Natural::power_of_two(64));
    let seeds = [
        Integer::ZERO,
        Integer::ONE,
        -Integer::ONE,
        big.clone(),
        -big,
    ];
    for seed in &seeds {
        assert_eq!(Generator::seeded(seed), Generator::seeded(seed), "{seed}");
    }
    for (i, a) in seeds.iter().enumerate() {
        for b in &seeds[i + 1..] {
            assert_ne!(first(a), first(b), "{a} and {b}");
        }
    }
}

/// `bits(0)` is 0; the counts of `below(6)` over 60000 draws are within 5
/// standard deviations of 10000 each.
#[test]
fn uniform() {
    let mut g = Generator::seeded(&Integer::from(7));
    assert!(g.bits(0).is_zero());
    let six = Natural::from(6u64);
    let mut counts = [0u32; 6];
    for _ in 0..60_000 {
        let r = g.below(&six);
        let i = usize::try_from(&r).expect("small");
        counts[i] += 1;
    }
    // The standard deviation of a count: √(60000 · 1/6 · 5/6) ≈ 91.
    for c in counts {
        assert!((9_544..=10_456).contains(&c), "{counts:?}");
    }
}

proptest! {
    /// `bits(k)` is below `2^k`, and `below(n)` below `n`.
    #[test]
    fn ranges(seed in any::<i64>(), k in 0u64..300, limbs in prop::collection::vec(any::<u64>(), 1..5)) {
        let mut g = Generator::seeded(&Integer::from(seed));
        prop_assert!(g.bits(k) < Natural::power_of_two(k));
        let n = Natural::from_limbs(&limbs);
        if !n.is_zero() {
            for _ in 0..8 {
                prop_assert!(g.below(&n) < n);
            }
        }
    }
}
