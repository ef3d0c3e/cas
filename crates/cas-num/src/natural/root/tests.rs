//! `Natural` roots and perfect powers checked against GMP, with the root
//! bounds checked directly.

use proptest::prelude::*;
use rug::Integer;
use rug::ops::Pow;

use crate::Natural;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

/// `s^n + d` for `d` in `-1..=1`, or one of the extremes `s^n + 2s` and
/// `(s + 1)^n` of a square root's remainder: the values where a root or its
/// correction step is off by one.
fn near_power(max_limbs: usize, n: u32) -> BoxedStrategy<Mag> {
    (strategy::natural_up_to(max_limbs), 0..5u8)
        .prop_map(move |(s, shape)| {
            let s = s.to_rug();
            let p = Integer::from((&s).pow(n));
            let x = match shape {
                0 => p - 1u32,
                1 => p,
                2 => p + 1u32,
                3 => p + Integer::from(&s * 2u32),
                _ => (s + 1u32).pow(n),
            };
            Mag::from_rug(&x.max(Integer::new()))
        })
        .boxed()
}

/// `s² + 64 j P`, with `P` the product of the residue primes: congruent to a
/// square modulo everything `is_perfect_square` filters by, but rarely a
/// square, so only the root itself can reject it.
fn filter_passing() -> BoxedStrategy<Mag> {
    (strategy::natural_up_to(DEFAULT_MAX_LIMBS / 2), 1..=1000u32)
        .prop_map(|(s, j)| {
            let p = Integer::from(super::RESIDUE_PRODUCT) * 64u32 * j;
            Mag::from_rug(&(s.to_rug().square() + p))
        })
        .boxed()
}

fn square_input() -> BoxedStrategy<Mag> {
    prop_oneof![
        3 => strategy::natural(),
        3 => near_power(DEFAULT_MAX_LIMBS / 2, 2),
        1 => filter_passing(),
    ]
    .boxed()
}

/// A root degree: small, around the 64-bit word, or large.
fn degree() -> BoxedStrategy<u32> {
    prop_oneof![
        6 => 1..=12u32,
        2 => 13..=200u32,
        1 => prop::sample::select(vec![1000, 4095, 65_536, u32::MAX]),
    ]
    .boxed()
}

/// A degree with a value near its powers; bases shrink as the degree grows.
fn root_input() -> BoxedStrategy<(Mag, u32)> {
    let near = (2..=40u32).prop_flat_map(|n| {
        let limbs = (60 / n as usize).max(1);
        (near_power(limbs, n), Just(n))
    });
    prop_oneof![
        2 => (strategy::natural(), degree()),
        2 => near,
    ]
    .boxed()
}

/// `a^k` with a small base and exponent, near misses, and powers with a
/// factor of two split off, which exercises the trailing-zero rule.
fn power_input() -> BoxedStrategy<Mag> {
    let power =
        (strategy::natural_up_to(4), 2..=24u32, -1i32..=1, 0..3u32).prop_map(|(a, k, d, twos)| {
            let a = a.to_rug() << twos;
            let x = a.pow(k) + d;
            Mag::from_rug(&x.max(Integer::new()))
        });
    let twos = (0..2000u32, -1i32..=1)
        .prop_map(|(t, d)| Mag::from_rug(&((Integer::from(1) << t) + d).max(Integer::new())));
    let mixed = (2..1000u32, 1..=6u32, 0..=6u32)
        .prop_map(|(m, k, t)| Mag::from_rug(&(Integer::from(m).pow(k) << (k * t))));
    prop_oneof![
        4 => power,
        2 => strategy::natural(),
        1 => twos,
        1 => mixed,
    ]
    .boxed()
}

differential! {
    fn sqrt_rem(a in square_input()) {
        ours: {
            let x = a.to_natural();
            let (s, r) = x.sqrt_rem();
            prop_assert!(s.square() <= x && x < (&s + 1).square(), "root bounds");
            prop_assert_eq!(&x.isqrt(), &s);
            (s, r)
        },
        gmp: a.to_rug().sqrt_rem(Integer::new()),
    }

    fn root_rem((a, n) in root_input()) {
        ours: {
            let x = a.to_natural();
            let (s, r) = x.root_rem(n);
            prop_assert_eq!(&x.nth_root(n), &s);
            (s, r)
        },
        gmp: a.to_rug().root_rem(Integer::new(), n),
    }

    fn is_perfect_square(a in square_input()) {
        ours: a.to_natural().is_perfect_square(),
        gmp: a.to_rug().is_perfect_square(),
    }

    fn is_perfect_power(a in power_input()) {
        ours: a.to_natural().is_perfect_power(),
        gmp: a.to_rug().is_perfect_power(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// The root bounds on large values (several recursion levels).
    #[test]
    fn sqrt_rem_large(a in strategy::natural_up_to(400)) {
        let x = a.to_natural();
        let (s, r) = x.sqrt_rem();
        prop_assert!(s.square() + &r == x && r <= &s + &s, "root bounds");
    }

    #[test]
    fn root_rem_large(a in strategy::natural_up_to(400), n in 3..=9u32) {
        let x = a.to_natural();
        let (s, r) = x.root_rem(n);
        prop_assert!(r < (&s + 1).pow(u64::from(n)) - s.pow(u64::from(n)), "root bounds");
    }
}

proptest! {
    /// The largest exponent: `a^k` is the value, `a` is not a perfect
    /// power, and there is one exactly when GMP finds a perfect power.
    #[test]
    fn perfect_power(a in power_input()) {
        let x = a.to_natural();
        let found = x.perfect_power();
        let gmp = a.to_rug() > 1 && a.to_rug().is_perfect_power();
        prop_assert_eq!(found.is_some(), gmp);
        if let Some((base, k)) = found {
            prop_assert!(k >= 2 && base.pow(k) == x, "not a power");
            prop_assert!(!base.is_perfect_power(), "exponent not largest");
        }
    }
}

#[test]
fn perfect_power_examples() {
    let pp = |x: u64| {
        Natural::from(x)
            .perfect_power()
            .map(|(a, k)| (u64::try_from(a).unwrap(), k))
    };
    assert_eq!(pp(0), None);
    assert_eq!(pp(1), None);
    assert_eq!(pp(2), None);
    assert_eq!(pp(4), Some((2, 2)));
    assert_eq!(pp(64), Some((2, 6)));
    assert_eq!(pp(72), None);
    assert_eq!(pp(729), Some((3, 6)));
    assert_eq!(pp(1 << 63), Some((2, 63)));
    assert_eq!(pp(3u64.pow(40)), Some((3, 40)));
    assert_eq!(pp(1_000_000), Some((10, 6)));
    assert!(Natural::ZERO.is_perfect_power() && Natural::one().is_perfect_power());
    let big = Natural::from(10u64).pow(1000);
    assert_eq!(big.perfect_power(), Some((Natural::from(10u64), 1000)));
    let prime_power = Natural::from(1_000_003u64).pow(97);
    assert_eq!(
        prime_power.perfect_power(),
        Some((Natural::from(1_000_003u64), 97))
    );
    assert_eq!((prime_power + 1).perfect_power(), None);
}

#[test]
fn root_edge_cases() {
    let two_64 = Natural::power_of_two(64);
    assert_eq!(Natural::ZERO.nth_root(7), Natural::ZERO);
    assert_eq!(two_64.nth_root(64), Natural::from(2u64));
    assert_eq!((&two_64 - 1).nth_root(64), Natural::one());
    assert_eq!(two_64.nth_root(65), Natural::one());
    assert_eq!(two_64.nth_root(u32::MAX), Natural::one());
    assert_eq!(Natural::from(u64::MAX).nth_root(u32::MAX), Natural::one());
    assert_eq!(
        Natural::from(u64::MAX).isqrt(),
        Natural::from(u64::from(u32::MAX))
    );
    let (s, r) = Natural::power_of_two(1001).sqrt_rem();
    assert_eq!(s.square() + r, Natural::power_of_two(1001));
}

#[test]
#[should_panic(expected = "zeroth root")]
fn zeroth_root_panics() {
    let _ = Natural::one().nth_root(0);
}

/// A degree for the kth root kernel, and its root's largest length in limbs:
/// enough for several Newton steps and the guard bits at small degrees,
/// a single limb (its first step) at large ones.
fn kth_degree() -> BoxedStrategy<(u32, usize)> {
    prop::sample::select(vec![(3, 40), (4, 30), (5, 24), (7, 17), (11, 10), (64, 1)]).boxed()
}

/// `(x, k)` around the kth power of a root `s` with up to the degree's
/// limbs (random, all ones or a power of two): `s^k + d` for `d` in
/// `-1..=1`, `(s + 1)^k - 1` (the largest remainder), or `s^k` plus a random
/// remainder.
fn kth_root_input() -> BoxedStrategy<(Mag, u32)> {
    kth_degree()
        .prop_flat_map(|(k, limbs)| {
            let bits = 1..=64 * u32::try_from(limbs).unwrap();
            let root = prop_oneof![
                3 => strategy::natural_up_to(limbs),
                1 => bits.clone().prop_map(|b| Mag::from_rug(&((Integer::from(1) << b) - 1u32))),
                1 => bits.prop_map(|b| Mag::from_rug(&(Integer::from(1) << (b - 1)))),
            ];
            (
                root,
                0..5u8,
                strategy::natural_up_to(limbs * k as usize),
                Just(k),
            )
        })
        .prop_map(|(s, shape, extra, k)| {
            let s = s.to_rug();
            let p = Integer::from((&s).pow(k));
            let next = Integer::from(&s + 1u32).pow(k);
            let x = match shape {
                0 => p - 1u32,
                1 => p,
                2 => p + 1u32,
                3 => next - 1u32,
                _ => {
                    let gap = Integer::from(&next - &p);
                    p + extra.to_rug() % gap
                }
            };
            (Mag::from_rug(&x.max(Integer::new())), k)
        })
        .boxed()
}

differential! {
    /// Every number of Newton steps, with and without the guard bits.
    fn kth_root_near_powers((a, k) in kth_root_input()) {
        ours: {
            let x = a.to_natural();
            let (s, r) = x.root_rem(k);
            prop_assert_eq!(&x.nth_root(k), &s);
            (s, r)
        },
        gmp: a.to_rug().root_rem(Integer::new(), k),
    }
}

/// Runs `f` and returns how often it reached `counter`.
fn count(counter: &'static std::thread::LocalKey<std::cell::Cell<u64>>, f: impl FnOnce()) -> u64 {
    counter.with(|c| c.set(0));
    f();
    counter.with(std::cell::Cell::get)
}

/// `nth_root` and `root_rem` of `x` against GMP.
fn check_root(x: &Integer, k: u32) {
    let (s, r) = Mag::from_rug(x).to_natural().root_rem(k);
    let (gs, gr) = x.clone().root_rem(Integer::new(), k);
    assert_eq!(
        (s.to_rug(), r.to_rug()),
        (gs.clone(), gr),
        "root_rem, k = {k}"
    );
    assert_eq!(
        Mag::from_rug(x).to_natural().nth_root(k).to_rug(),
        gs,
        "nth_root, k = {k}"
    );
}

/// Inputs that reach the kth root's rare branches: `(s + 1)^k - 1`, whose
/// remainder is the largest, makes the last step's estimate too large, and
/// all-ones values make it reach past its bits; exact powers make the
/// guard bits ambiguous, and with roots of about 46 bits the floating-point
/// estimate of the first step too small.
#[test]
fn kth_root_rare_branches() {
    use super::nth::{FIRST_INCREMENT, GUARD_CHECK, STEP_CLAMP, STEP_DECREMENT};
    let mut hits = [0; 4];
    for k in [3, 5, 7] {
        for i in 1..40u32 {
            let s = Integer::from(3).pow(40 + 7 * i) + i;
            hits[0] += count(&STEP_DECREMENT, || {
                check_root(&(Integer::from(&s + 1u32).pow(k) - 1u32), k);
            });
            let ones = (Integer::from(1) << (k * (64 + 13 * i))) - 1u32;
            hits[1] += count(&STEP_CLAMP, || check_root(&ones, k));
            hits[2] += count(&GUARD_CHECK, || check_root(&Integer::from((&s).pow(k)), k));
            let small = Integer::from(1u32 << 15) * (1u32 << 30) + Integer::from(3).pow(10 + i);
            hits[3] += count(&FIRST_INCREMENT, || check_root(&small.pow(k), k));
        }
    }
    assert!(hits.iter().all(|&h| h > 0), "{hits:?}");
}

/// Where a step's division would just pass a transform length, the step
/// takes fewer bits; the roots stay exact.
#[test]
fn kth_root_fits_transforms() {
    use super::nth::STEP_FIT;
    use crate::kernel::mul::{MUL_WRAP_NTT_THRESHOLD, wrap_len};
    use crate::testkit::TESTED_SIZE_MAX;
    if MUL_WRAP_NTT_THRESHOLD > TESTED_SIZE_MAX {
        return;
    }
    // With c = 8 for cube roots, the last step to a root of 128 f - 56
    // bits appends 64 f - 32, so its division works on f + 2 limbs, two
    // past the transform length f (at least 64, for two limbs of slack).
    let fit = u32::try_from(wrap_len(MUL_WRAP_NTT_THRESHOLD.max(64))).unwrap();
    let bits = 128 * fit - 56;
    let s = (Integer::from(1) << (bits - 1)) + Integer::from(3).pow(bits / 4);
    let fits = count(&STEP_FIT, || {
        check_root(&(Integer::from(&s + 1u32).pow(3) - 1u32), 3);
        check_root(&Integer::from((&s).pow(3)), 3);
    });
    assert!(fits > 0);
}

/// Exact powers (and ±1) whose root ends in a tiny last digit `e` below a
/// top part with random low bits. The last step's dividend `D` then exceeds
/// `e k s'^(k-1)` by a small amount, and its low bits, the input's next
/// bits, are what bring it there.
#[test]
fn kth_root_low_digits() {
    for k in [3u32, 5, 7] {
        let c = u32::BITS - (k - 1).leading_zeros() + 6;
        for n in [100u32, 129, 200, 250] {
            // The last step appends b bits (see `step_bits`).
            let b = (n - c) / 2;
            let half = Integer::from(1) << (n - b - 1);
            let top = &half + Integer::from(3).pow(2 * n) % &half;
            for e in 1..4u32 {
                let s = (Integer::from(&top) << b) + e;
                let p = Integer::from((&s).pow(k));
                check_root(&p, k);
                check_root(&(Integer::from(&p) - 1u32), k);
                check_root(&(p + 1u32), k);
            }
        }
    }
}

/// Cube roots of `L/2` to `L` limbs, `L` the shortest transform length from
/// which products modulo `B^m - 1` take the cyclic transform: their steps'
/// last products go modulo `B^m - 1`, at lengths just above a transform
/// length (the limbs above `m` from a short product) and below one
/// (directly). `(s + 1)^3 - 1` makes the last step's estimate too large, so
/// that the wrapped remainder comes out negative.
#[test]
fn kth_root_wraps_last_product() {
    use super::nth::{STEP_DECREMENT, WRAP_DIRECT, WRAP_LOW};
    use crate::kernel::mul::{MUL_WRAP_NTT_THRESHOLD, wrap_len};
    use crate::testkit::TESTED_SIZE_MAX;
    if 2 * MUL_WRAP_NTT_THRESHOLD > TESTED_SIZE_MAX {
        return;
    }
    let l = wrap_len(MUL_WRAP_NTT_THRESHOLD.max(8));
    let mut hits = [0; 3];
    for limbs in [l / 2, l / 2 + 1, 5 * l / 8, 3 * l / 4, l] {
        let bits = 64 * u32::try_from(limbs).unwrap();
        let s = (Integer::from(1) << (bits - 1)) + Integer::from(3).pow(bits / 3);
        for x in [
            Integer::from(&s + 1u32).pow(3) - 1u32,
            Integer::from((&s).pow(3)),
            (Integer::from(1) << (3 * bits)) - 1u32,
        ] {
            hits[0] += count(&WRAP_LOW, || {
                hits[1] += count(&WRAP_DIRECT, || {
                    hits[2] += count(&STEP_DECREMENT, || check_root(&x, 3));
                });
            });
        }
    }
    assert!(hits.iter().all(|&h| h > 0), "{hits:?}");
}
