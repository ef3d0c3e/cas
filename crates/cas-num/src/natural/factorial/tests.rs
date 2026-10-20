//! Factorials checked against GMP.

use proptest::prelude::*;
use rug::Integer;

use crate::Natural;
use crate::testkit::{ToRug, differential};

/// Every `n` up to a few thousand, weighted towards the end of the one-limb
/// results and the first levels of the odd-part recursion.
fn argument() -> BoxedStrategy<u32> {
    prop_oneof![
        2 => 0..=70u32,
        1 => 70..=5000u32,
    ]
    .boxed()
}

differential! {
    fn factorial(n in argument()) {
        ours: Natural::factorial(u64::from(n)),
        gmp: Integer::factorial(n),
    }
}

#[test]
fn small_factorials() {
    let expected = [1u64, 1, 2, 6, 24, 120, 720, 5040, 40320];
    for (n, f) in (0..).zip(expected) {
        assert_eq!(Natural::factorial(n), Natural::from(f), "{n}!");
    }
    assert_eq!(
        Natural::factorial(20),
        Natural::from(2_432_902_008_176_640_000u64)
    );
}

/// Sizes where the product trees and the last products run through the
/// subquadratic multiplications.
#[test]
fn large_factorials() {
    for n in [30_000u32, 100_000] {
        let ours = Natural::factorial(u64::from(n)).to_rug();
        assert!(ours == Integer::from(Integer::factorial(n)), "{n}!");
    }
}
