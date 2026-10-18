//! `mulhi` against the exact high half from GMP: the result must be a lower
//! bound within `mulhi_error(n)`, at every algorithm and on all-ones
//! operands, where the dropped columns are largest.

use proptest::prelude::*;
use rug::Integer;

use super::{
    Mulhi, dc_split, mulhi, mulhi_algorithm, mulhi_basecase, mulhi_error, mulhi_scratch_len,
};
use crate::Limb;
use crate::kernel::mul::mullo::{MULLO_DC_THRESHOLD, MULLO_MUL_THRESHOLD};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, TESTED_SIZE_MAX, ToRug, tested};

/// Filler for output and scratch buffers, so results can't rely on zeroed memory.
const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn run_mulhi(a: &[Limb], b: &[Limb]) -> Integer {
    let mut r = vec![GARBAGE; a.len()];
    let mut scratch = vec![GARBAGE; mulhi_scratch_len(a.len())];
    mulhi(&mut r, a, b, &mut scratch);
    Mag::from_limbs(r).to_rug()
}

/// `floor(a b / B^n)`.
fn exact_high(a: &[Limb], b: &[Limb]) -> Integer {
    let n = u32::try_from(a.len()).unwrap();
    (Mag::from_limbs(a).to_rug() * Mag::from_limbs(b).to_rug()) >> (64 * n)
}

/// `Ok(error)` if `r` is a lower bound of the exact high half within the
/// documented error.
fn check(a: &[Limb], b: &[Limb], r: &Integer) -> Result<usize, TestCaseError> {
    let exact = exact_high(a, b);
    let err = Integer::from(&exact - r);
    let bound = mulhi_error(a.len());
    prop_assert!(
        err >= 0,
        "n = {}: above the exact high half by {}",
        a.len(),
        -err
    );
    prop_assert!(
        err <= bound,
        "n = {}: error {} beyond the bound {}",
        a.len(),
        err,
        bound
    );
    Ok(err.to_usize().unwrap())
}

/// Lengths around both thresholds and where the split recurses into its own
/// range.
fn mulhi_lens() -> BoxedStrategy<usize> {
    len_around(&[
        1,
        2,
        3,
        MULLO_DC_THRESHOLD,
        2 * MULLO_DC_THRESHOLD,
        4 * MULLO_DC_THRESHOLD,
        MULLO_MUL_THRESHOLD,
    ])
    .boxed()
}

proptest! {
    #[test]
    fn mulhi_within_bound(
        (a, b) in prop_oneof![
            mulhi_lens(),
            1..=200usize,
        ].prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
    ) {
        let r = run_mulhi(&a, &b);
        check(&a, &b, &r)?;
    }

    #[test]
    fn mulhi_basecase_within_bound(
        (a, b) in (1..=60usize).prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
    ) {
        let mut r = vec![GARBAGE; a.len()];
        let mut t = vec![GARBAGE; a.len() + 2];
        mulhi_basecase(&mut r, &a, &b, &mut t);
        let exact = exact_high(&a, &b);
        let err = &exact - Mag::from_limbs(r).to_rug();
        prop_assert!(err >= 0 && err < a.len().max(1));
    }
}

/// All-ones operands maximise every dropped column, at each algorithm.
#[test]
fn mulhi_all_ones() {
    for n in [
        1,
        2,
        3,
        17,
        MULLO_DC_THRESHOLD - 1,
        MULLO_DC_THRESHOLD,
        MULLO_DC_THRESHOLD + 1,
        3 * MULLO_DC_THRESHOLD + 1,
        MULLO_MUL_THRESHOLD,
    ] {
        let n = tested(n);
        let ones = vec![Limb::MAX; n];
        let r = run_mulhi(&ones, &ones);
        check(&ones, &ones, &r).unwrap();
    }
}

/// One limb and the full-product range are exact.
#[test]
fn mulhi_exact_cases() {
    let a = [0x1234_5678_9abc_def0];
    let b = [Limb::MAX];
    assert_eq!(run_mulhi(&a, &b), exact_high(&a, &b));
    assert_eq!(mulhi_error(1), 0);
    if MULLO_MUL_THRESHOLD <= TESTED_SIZE_MAX {
        let n = MULLO_MUL_THRESHOLD;
        assert_eq!(mulhi_algorithm(n), Mulhi::Full);
        assert_eq!(mulhi_error(n), 0);
    }
}

/// The dispatcher's choices, the split's range, and that the bound stays
/// small: a few hundred ulps at most up to the full-product threshold.
#[test]
fn mulhi_algorithm_and_bound_shape() {
    assert_eq!(mulhi_algorithm(1), Mulhi::Full);
    assert_eq!(mulhi_algorithm(2), Mulhi::Full);
    assert_eq!(mulhi_error(2), 0);
    if MULLO_DC_THRESHOLD < MULLO_MUL_THRESHOLD && MULLO_DC_THRESHOLD > 3 {
        assert_eq!(mulhi_algorithm(3), Mulhi::Basecase);
        assert_eq!(mulhi_algorithm(MULLO_DC_THRESHOLD - 1), Mulhi::Basecase);
        assert_eq!(mulhi_algorithm(MULLO_DC_THRESHOLD), Mulhi::Dc);
    }
    for n in 2..=tested(MULLO_MUL_THRESHOLD) {
        let n1 = dc_split(n);
        assert!((1..=n / 2).contains(&n1), "n = {n}, n1 = {n1}");
        assert!(
            mulhi_error(n) <= n.max(1000),
            "n = {n}: bound {}",
            mulhi_error(n)
        );
    }
}
