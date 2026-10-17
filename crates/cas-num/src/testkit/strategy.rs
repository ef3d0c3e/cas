//! Proptest strategies producing numbers with adversarial shapes.
//!
//! Uniformly random limbs almost never hit the inputs that break arithmetic
//! (carries rippling through all-ones limbs, values just past a limb boundary,
//! operands of very different sizes), so the generators below weight those
//! shapes heavily. Every generated value is normalised.

use proptest::collection::vec;
use proptest::prelude::*;

use super::model::{Mag, Signed};

/// Default size bound, in limbs, for the unparameterised strategies.
pub(crate) const DEFAULT_MAX_LIMBS: usize = 40;

/// Limb values that sit on carry/borrow and sign boundaries.
const SPECIAL_LIMBS: [u64; 11] = [
    0,
    1,
    2,
    u64::MAX,
    u64::MAX - 1,
    1 << 63,
    (1 << 63) - 1,
    (1 << 63) + 1,
    1 << 32,
    (1 << 32) - 1,
    0x5555_5555_5555_5555,
];

/// A single limb: special values, single-bit and low-mask values, or uniform.
pub(crate) fn limb() -> impl Strategy<Value = u64> + Clone {
    prop_oneof![
        3 => proptest::sample::select(&SPECIAL_LIMBS[..]),
        1 => (0..64u32).prop_map(|k| 1u64 << k),
        1 => (0..64u32).prop_map(|k| (1u64 << k) - 1),
        4 => any::<u64>(),
    ]
}

/// A non-zero limb, for the most significant position.
pub(crate) fn top_limb() -> impl Strategy<Value = u64> + Clone {
    limb().prop_map(|x| if x == 0 { 1 } else { x })
}

/// A length in `1..=max`, weighted towards 1–4 limbs, where most special
/// cases (single limb, first carry into a new limb, inline storage) live.
pub(crate) fn len_up_to(max: usize) -> BoxedStrategy<usize> {
    assert!(max >= 1);
    let small = 1..=max.min(4);
    if max <= 4 {
        small.boxed()
    } else {
        prop_oneof![3 => small, 1 => 5..=max].boxed()
    }
}

/// Lengths `t - 1`, `t` and `t + 1` for each algorithm threshold `t`
/// (never below 1), so both sides of every crossover get exercised.
pub(crate) fn len_around(thresholds: &[usize]) -> BoxedStrategy<usize> {
    let mut lens: Vec<usize> = thresholds
        .iter()
        .flat_map(|&t| [t.saturating_sub(1), t, t + 1])
        .filter(|&n| n >= 1)
        .collect();
    lens.sort_unstable();
    lens.dedup();
    proptest::sample::select(lens).boxed()
}

/// `2^k`.
fn pow2(k: usize) -> Mag {
    let mut limbs = vec![0; k / 64 + 1];
    limbs[k / 64] = 1 << (k % 64);
    Mag::from_limbs(limbs)
}

/// `2^k - 1`: `k` one bits.
fn ones(k: usize) -> Mag {
    let mut limbs = vec![u64::MAX; k / 64];
    if k % 64 != 0 {
        limbs.push((1 << (k % 64)) - 1);
    }
    Mag::from_limbs(limbs)
}

/// `2^k + 1` (`k >= 1`).
fn pow2_plus_1(k: usize) -> Mag {
    let mut limbs = pow2(k).limbs().to_vec();
    limbs[0] |= 1;
    Mag::from_limbs(limbs)
}

/// A natural number with exactly `len` limbs (`len >= 1`), filled densely,
/// sparsely (mostly zero limbs), or with all-ones limbs.
pub(crate) fn natural_with_len(len: impl Strategy<Value = usize>) -> impl Strategy<Value = Mag> {
    len.prop_flat_map(|len| {
        assert!(len >= 1);
        let low = len - 1;
        let dense = vec(limb(), low).boxed();
        let sparse = vec(prop_oneof![6 => Just(0u64), 1 => limb()], low).boxed();
        let all_ones = Just(vec![u64::MAX; low]).boxed();
        let fill = prop_oneof![4 => dense, 2 => sparse, 1 => all_ones];
        (fill, top_limb()).prop_map(|(mut limbs, top)| {
            limbs.push(top);
            Mag::from_limbs(limbs)
        })
    })
}

/// A natural number of at most `max_limbs` limbs, mixing zero, one, single
/// limbs, powers of two, `2^k ± 1` around limb boundaries, all-ones values
/// and dense/sparse multi-limb values.
pub(crate) fn natural_up_to(max_limbs: usize) -> BoxedStrategy<Mag> {
    assert!(max_limbs >= 1);
    let max_bits = 64 * max_limbs;
    // Bit counts within one bit of a limb boundary: 64j - 1, 64j, 64j + 1.
    let boundary_bits =
        (1..=max_limbs, 0..3usize).prop_map(move |(j, o)| (64 * j + o - 1).min(max_bits - 1));
    prop_oneof![
        1 => Just(Mag::default()),
        1 => Just(Mag::from_u64(1)),
        3 => top_limb().prop_map(Mag::from_u64),
        1 => (0..max_bits).prop_map(pow2),
        2 => boundary_bits.clone().prop_map(pow2_plus_1),
        2 => boundary_bits.prop_map(ones),
        1 => (1..=max_bits).prop_map(ones),
        6 => natural_with_len(len_up_to(max_limbs)),
    ]
    .boxed()
}

/// [`natural_up_to`] with [`DEFAULT_MAX_LIMBS`].
pub(crate) fn natural() -> BoxedStrategy<Mag> {
    natural_up_to(DEFAULT_MAX_LIMBS)
}

/// A signed integer of at most `max_limbs` limbs.
pub(crate) fn signed_up_to(max_limbs: usize) -> BoxedStrategy<Signed> {
    (any::<bool>(), natural_up_to(max_limbs))
        .prop_map(|(neg, mag)| Signed::new(neg, mag))
        .boxed()
}

/// [`signed_up_to`] with [`DEFAULT_MAX_LIMBS`].
pub(crate) fn signed() -> BoxedStrategy<Signed> {
    signed_up_to(DEFAULT_MAX_LIMBS)
}

/// A pair of naturals of very different sizes: one of `max_limbs / 2 ..= max_limbs`
/// limbs, one of 1–3 limbs, in either order.
pub(crate) fn unbalanced_pair(max_limbs: usize) -> BoxedStrategy<(Mag, Mag)> {
    assert!(max_limbs >= 8);
    let long = natural_with_len(max_limbs / 2..=max_limbs);
    let short = natural_with_len(1..=3usize);
    (long, short, any::<bool>())
        .prop_map(|(long, short, swap)| if swap { (short, long) } else { (long, short) })
        .boxed()
}
