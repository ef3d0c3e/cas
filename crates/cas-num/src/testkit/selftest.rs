//! Tests of the harness itself: conversions against an independent
//! construction, and coverage of the adversarial shapes.

use proptest::prelude::*;
use proptest::strategy::ValueTree;
use proptest::test_runner::TestRunner;
use rug::Integer;

use super::strategy::{self, DEFAULT_MAX_LIMBS};
use super::{Mag, Signed, ToRug, check_same, differential};

/// `sum limbs[i] * 2^(64 i)` built with GMP shifts and additions only, as an
/// oracle for the limb-slice conversion that does not share its code path.
fn horner(limbs: &[u64]) -> Integer {
    limbs
        .iter()
        .rev()
        .fold(Integer::new(), |acc, &l| (acc << 64u32) + l)
}

differential! {
    /// The trivial op: u64 into our limb form agrees with GMP.
    fn from_u64(a in any::<u64>()) {
        ours: Mag::from_u64(a),
        gmp: Integer::from(a),
    }

    /// ... and comes back out unchanged.
    fn from_u64_round_trip(a in any::<u64>()) {
        ours: Mag::from_rug(&Integer::from(a)).to_u64(),
        gmp: Some(Integer::from(a)),
    }

    fn from_i64(a in any::<i64>()) {
        ours: Signed::from_i64(a),
        gmp: Integer::from(a),
    }

    fn mag_to_rug_matches_horner(a in strategy::natural()) {
        ours: a.clone(),
        gmp: horner(a.limbs()),
    }

    fn signed_to_rug_matches_horner(a in strategy::signed()) {
        ours: a.clone(),
        gmp: if a.is_negative() { -horner(a.mag().limbs()) } else { horner(a.mag().limbs()) },
    }
}

proptest! {
    #[test]
    fn mag_round_trip(a in strategy::natural()) {
        prop_assert_eq!(Mag::from_rug(&a.to_rug()), a);
    }

    #[test]
    fn signed_round_trip(a in strategy::signed()) {
        prop_assert_eq!(Signed::from_rug(&a.to_rug()), a);
    }

    #[test]
    fn generated_values_are_normalised(a in strategy::signed()) {
        prop_assert_ne!(a.mag().limbs().last(), Some(&0));
        prop_assert!(!(a.is_negative() && a.mag().is_zero()));
        prop_assert!(a.mag().limbs().len() <= DEFAULT_MAX_LIMBS);
    }

    #[test]
    fn natural_with_len_is_exact(
        (len, a) in (1..64usize).prop_flat_map(|n| (Just(n), strategy::natural_with_len(Just(n))))
    ) {
        prop_assert_eq!(a.limbs().len(), len);
    }

    #[test]
    fn len_around_straddles(t in 2..1000usize) {
        let lens = strategy::len_around(&[t]);
        let mut runner = TestRunner::deterministic();
        for _ in 0..32 {
            let n = lens.new_tree(&mut runner).unwrap().current();
            prop_assert!((t - 1..=t + 1).contains(&n));
        }
    }

    #[test]
    fn unbalanced_pair_is_unbalanced((a, b) in strategy::unbalanced_pair(32)) {
        let (la, lb) = (a.limbs().len(), b.limbs().len());
        prop_assert!(la.max(lb) >= 16 && la.min(lb) <= 3);
    }
}

#[test]
fn natural_covers_adversarial_shapes() {
    let strat = strategy::natural();
    let mut runner = TestRunner::deterministic();
    let (mut zero, mut one, mut single, mut long) = (0, 0, 0, 0);
    let (mut all_ones, mut pow2, mut pow2_plus_1, mut sparse) = (0, 0, 0, 0);
    for _ in 0..10_000 {
        let a = strat.new_tree(&mut runner).unwrap().current();
        let limbs = a.limbs();
        let (low, top) = match limbs.split_last() {
            None => {
                zero += 1;
                continue;
            }
            Some((top, low)) => (low, *top),
        };
        if limbs == [1] {
            one += 1;
        }
        if limbs.len() == 1 {
            single += 1;
            continue;
        }
        long += usize::from(limbs.len() > 4);
        let low_all = |x| low.iter().all(|&l| l == x);
        all_ones += usize::from(low_all(u64::MAX) && top & top.wrapping_add(1) == 0);
        pow2 += usize::from(low_all(0) && top.is_power_of_two());
        pow2_plus_1 +=
            usize::from(low[0] == 1 && low[1..].iter().all(|&l| l == 0) && top.is_power_of_two());
        sparse += usize::from(low.iter().filter(|&&l| l == 0).count() * 2 > low.len());
    }
    for (name, n) in [
        ("zero", zero),
        ("one", one),
        ("single limb", single),
        ("more than 4 limbs", long),
        ("multi-limb 2^k - 1", all_ones),
        ("multi-limb 2^k", pow2),
        ("multi-limb 2^k + 1", pow2_plus_1),
        ("sparse", sparse),
    ] {
        assert!(n >= 100, "only {n}/10000 samples of shape {name}");
    }
}

#[test]
fn check_same_reports_mismatch() {
    let err = check_same(Mag::from_u64(1), Integer::from(2)).unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("differs from GMP") && msg.contains("0x1"),
        "{msg}"
    );
}
