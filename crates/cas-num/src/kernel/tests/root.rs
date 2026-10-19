//! The square root kernel checked against GMP, on random operands, on
//! squares and their neighbours (where the correction step runs), and at the
//! lengths where its division changes algorithm.

use std::cell::Cell;
use std::thread::LocalKey;

use proptest::prelude::*;
use rug::Integer;
use rug::ops::Pow;

use super::val;
use crate::Limb;
use crate::kernel::div::DC_DIV_THRESHOLD;
use crate::kernel::div::mu::MU_DIV_THRESHOLD;
use crate::kernel::root::{SQRT_CARRY, SQRT_CORRECTION, dc_sqrtrem, dc_sqrtrem_scratch_len};
use crate::kernel::{sqrtrem, sqrtrem_scratch_len};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, ToRug, differential, tested};

const GARBAGE: Limb = 0xdead_beef_dead_beef;

/// `sqrtrem` on `a`, with garbage in every output and scratch limb; checks
/// the returned length and that the remainder's limbs above it are zero.
fn run_sqrtrem(a: &[Limb]) -> (Mag, Mag) {
    let an = a.len();
    let mut s = vec![GARBAGE; an.div_ceil(2)];
    let mut r = vec![GARBAGE; an / 2 + 1];
    let mut scratch = vec![GARBAGE; sqrtrem_scratch_len(an)];
    let rn = sqrtrem(&mut s, &mut r, a, &mut scratch);
    assert!(rn == 0 || r[rn - 1] != 0, "remainder length not normalised");
    assert!(r[rn..].iter().all(|&x| x == 0), "limbs above the remainder");
    assert!(s.last() != Some(&0), "root's top limb is zero");
    (Mag::from_limbs(s), Mag::from_limbs(r))
}

fn gmp_sqrtrem(a: &[Limb]) -> (Integer, Integer) {
    val(a).sqrt_rem(Integer::new())
}

/// Root lengths: the shortest, then around the lengths at which the top
/// division (by about half the root) takes divide and conquer or Barrett's
/// method, or splits again.
fn root_lens() -> BoxedStrategy<usize> {
    prop_oneof![
        2 => 1..=12usize,
        1 => len_around(&[
            2 * DC_DIV_THRESHOLD,
            4 * DC_DIV_THRESHOLD,
            2 * MU_DIV_THRESHOLD,
            4 * MU_DIV_THRESHOLD,
        ]),
    ]
    .boxed()
}

/// An operand for a root of `n` limbs: random with any top limb (so either
/// length parity and any normalising shift), `s² + e` for `e` in `-1..=1`,
/// `s² + 2s` (the largest remainder) or `(s + 1)² - 1`, or all ones.
fn operand() -> BoxedStrategy<Vec<Limb>> {
    root_lens()
        .prop_flat_map(|n| {
            let an =
                prop_oneof![Just(2 * n), Just(2 * n - 1)].prop_filter("non-empty", |&an| an >= 1);
            let random = an
                .clone()
                .prop_flat_map(|an| strategy::natural_with_len(Just(an)))
                .prop_map(|a| a.limbs().to_vec());
            let near = (strategy::natural_with_len(Just(n)), 0..5u8).prop_map(|(s, shape)| {
                let s = s.to_rug();
                let p = Integer::from((&s).pow(2u32));
                let x = match shape {
                    0 => p - 1u32,
                    1 => p,
                    2 => p + 1u32,
                    3 => p + Integer::from(&s * 2u32),
                    _ => (s + 1u32).pow(2u32) - 1u32,
                };
                Mag::from_rug(&x.max(Integer::from(1))).limbs().to_vec()
            });
            let ones = an.prop_map(|an| vec![Limb::MAX; an]);
            prop_oneof![3 => random, 3 => near, 1 => ones]
        })
        .boxed()
}

differential! {
    fn sqrtrem_matches(a in operand()) {
        ours: run_sqrtrem(&a),
        gmp: gmp_sqrtrem(&a),
    }
}

/// `sqrtrem` on `a` against GMP.
fn check(a: &[Limb]) {
    let (s, r) = run_sqrtrem(a);
    assert_eq!(
        (s.to_rug(), r.to_rug()),
        gmp_sqrtrem(a),
        "{} limbs",
        a.len()
    );
}

fn count(counter: &'static LocalKey<Cell<u64>>, f: impl FnOnce()) -> u64 {
    counter.with(|c| c.set(0));
    f();
    counter.with(Cell::get)
}

/// `s² - 1` has the root `s - 1`, which the estimates overshoot, so the
/// correction step runs. When the top half's remainder is the largest,
/// `2s'`, as for all-ones operands, the low half's estimate is `B^l` and
/// carries into the top half.
#[test]
fn sqrtrem_corrects() {
    for n in [2, 3, 4, 5, 8, 2 * tested(DC_DIV_THRESHOLD) + 1] {
        let a = vec![Limb::MAX; 2 * n];
        let carries = count(&SQRT_CARRY, || check(&a));
        assert!(carries >= 1, "no carry at n = {n}");
    }
    let mut corrections = 0;
    for n in 2..40 {
        let s = Integer::from(3).pow(u32::try_from(40 * n).unwrap());
        let x = Integer::from(s.square_ref()) - 1u32;
        let a = Mag::from_rug(&x).limbs().to_vec();
        corrections += count(&SQRT_CORRECTION, || check(&a));
    }
    assert!(corrections >= 1);
}

/// The scratch length covers the recursion, whose divisions take Barrett's
/// method at large sizes: run on exactly that much.
#[test]
fn dc_sqrtrem_at_barrett_sizes() {
    let mu = tested(MU_DIV_THRESHOLD);
    for n in [2 * mu - 1, 2 * mu, 2 * mu + 3] {
        let mut a: Vec<Limb> = (0..2 * n as u64)
            .map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15))
            .collect();
        a[2 * n - 1] |= 1 << 62;
        let mut work = a.clone();
        let mut root = vec![GARBAGE; n];
        let mut scratch = vec![GARBAGE; dc_sqrtrem_scratch_len(n)];
        let top = dc_sqrtrem(&mut root, &mut work, &mut scratch);
        let mut rem = work[..n].to_vec();
        rem.push(Limb::from(top));
        assert_eq!((val(&root), val(&rem)), gmp_sqrtrem(&a), "n = {n}");
    }
}
