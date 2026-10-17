//! Kernel semantics checked against GMP: each kernel's output limbs plus its
//! carry/borrow/shifted-out limb must equal the exact GMP result.

use std::cmp::Ordering;

use proptest::prelude::*;
use rug::Integer;
use rug::ops::RemRounding;

use super::*;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

/// Value of raw (possibly unnormalised) limbs.
fn val(a: &[Limb]) -> Integer {
    Integer::from_digits(a, rug::integer::Order::Lsf)
}

/// `2^(64 n)`.
fn base_pow(n: usize) -> Integer {
    Integer::from(1) << (64 * u32::try_from(n).unwrap())
}

/// `(x mod 2^(64 n), x div 2^(64 n))` for `x >= 0`.
fn split_at_limbs(x: Integer, n: usize) -> (Integer, Integer) {
    let (q, r) = x.div_rem_euc(base_pow(n));
    (r, q)
}

const MAX: usize = DEFAULT_MAX_LIMBS;

/// Runs `$kernel` (out-of-place) and `$assign` (in-place) on the same inputs,
/// checks they agree, and yields `(result, returned limb)`.
macro_rules! both {
    ($kernel:ident, $assign:ident, $a:expr, $($b:expr),*) => {{
        let a: &[Limb] = &$a;
        let mut r = vec![0xdead_beef; a.len()];
        let c = $kernel(&mut r, a, $($b),*);
        let mut r2 = a.to_vec();
        let c2 = $assign(&mut r2, $($b),*);
        prop_assert_eq!(&r, &r2, "in-place variant differs");
        prop_assert_eq!(c, c2, "in-place variant returns a different limb");
        (Mag::from_limbs(r), c)
    }};
}

differential! {
    fn add_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: both!(add_n, add_n_assign, a, &b),
        gmp: split_at_limbs(val(&a) + val(&b), a.len()),
    }

    fn add_matches((a, b) in strategy::raw_pair_longer_first(MAX)) {
        ours: both!(add, add_assign, a, &b),
        gmp: split_at_limbs(val(&a) + val(&b), a.len()),
    }

    fn add_1_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs), b in strategy::limb()) {
        ours: both!(add_1, add_1_assign, a, b),
        gmp: split_at_limbs(val(&a) + b, a.len()),
    }

    fn sub_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: both!(sub_n, sub_n_assign, a, &b),
        gmp: (
            (val(&a) - val(&b)).rem_euc(base_pow(a.len())),
            Integer::from(val(&a) < val(&b)),
        ),
    }

    fn sub_matches((a, b) in strategy::raw_pair_longer_first(MAX)) {
        ours: both!(sub, sub_assign, a, &b),
        gmp: (
            (val(&a) - val(&b)).rem_euc(base_pow(a.len())),
            Integer::from(val(&a) < val(&b)),
        ),
    }

    fn sub_1_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs), b in strategy::limb()) {
        ours: both!(sub_1, sub_1_assign, a, b),
        gmp: (
            (val(&a) - b).rem_euc(base_pow(a.len())),
            Integer::from(val(&a) < b),
        ),
    }

    fn shl_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs), s in 0..64u32) {
        ours: both!(shl, shl_assign, a, s),
        gmp: split_at_limbs(val(&a) << s, a.len()),
    }

    fn shr_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs), s in 0..64u32) {
        ours: both!(shr, shr_assign, a, s),
        gmp: (
            val(&a) >> s,
            // Bits shifted out, aligned to the top of a limb.
            (val(&a) & ((Integer::from(1) << s) - 1u32)) << (64 - s),
        ),
    }

    fn and_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: { let (r, ()) = both!(and_n, and_n_assign, a, &b); r },
        gmp: val(&a) & val(&b),
    }

    fn or_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: { let (r, ()) = both!(or_n, or_n_assign, a, &b); r },
        gmp: val(&a) | val(&b),
    }

    fn xor_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: { let (r, ()) = both!(xor_n, xor_n_assign, a, &b); r },
        gmp: val(&a) ^ val(&b),
    }

    fn cmp_n_matches((a, b) in strategy::raw_pair_same_len(MAX)) {
        ours: cmp_n(&a, &b),
        gmp: val(&a).cmp(&val(&b)),
    }

    fn cmp_matches(a in strategy::natural(), b in strategy::natural()) {
        ours: cmp(a.limbs(), b.limbs()),
        gmp: a.to_rug().cmp(&b.to_rug()),
    }

    fn normalized_len_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs)) {
        ours: normalized_len(&a),
        gmp: Mag::from_limbs(a.clone()).limbs().len(),
    }

    fn count_ones_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs)) {
        ours: count_ones(&a),
        gmp: val(&a).count_ones().unwrap(),
    }

    fn trailing_zeros_matches(a in (0..=MAX).prop_flat_map(strategy::raw_limbs)) {
        ours: trailing_zeros(&a),
        gmp: val(&a).find_one(0).map(Integer::from),
    }

    fn bit_length_matches(a in strategy::natural()) {
        ours: bit_length(a.limbs()),
        gmp: a.to_rug().significant_bits(),
    }
}

#[test]
fn add_1_on_empty_carries_the_limb() {
    assert_eq!(add_1(&mut [], &[], 7), 7);
    assert_eq!(add_1_assign(&mut [], 7), 7);
    assert_eq!(sub_1(&mut [], &[], 7), 1);
    assert_eq!(sub_1_assign(&mut [], 0), 0);
}

#[test]
fn carry_ripples_through_all_ones() {
    for n in 1..=8 {
        let mut a = vec![Limb::MAX; n];
        assert_eq!(add_1_assign(&mut a, 1), 1, "n = {n}");
        assert!(a.iter().all(|&x| x == 0));
        assert_eq!(sub_1_assign(&mut a, 1), 1, "n = {n}");
        assert!(a.iter().all(|&x| x == Limb::MAX));
        let b = a.clone();
        assert_eq!(add_n_assign(&mut a, &b), 1);
        assert_eq!(a[0], Limb::MAX - 1);
        assert!(a[1..].iter().all(|&x| x == Limb::MAX));
    }
}

#[test]
#[should_panic(expected = "length mismatch")]
fn length_mismatch_panics() {
    add_n(&mut [0; 2], &[1, 2], &[3]);
}

#[test]
fn cmp_orders_by_length_first() {
    assert_eq!(cmp(&[Limb::MAX], &[0, 1]), Ordering::Less);
    assert_eq!(cmp(&[], &[]), Ordering::Equal);
}

mod mul;
