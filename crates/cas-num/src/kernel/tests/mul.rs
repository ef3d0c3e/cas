//! Multiplication kernels: checked against GMP through the size dispatcher, and
//! Karatsuba and Toom-3 called directly, checked against schoolbook (the
//! independent second oracle) at every small split shape.

use std::collections::HashSet;

use proptest::prelude::*;
use rug::Integer;

use super::{split_at_limbs, val};
use crate::Limb;
use crate::kernel::mul::mullo::{
    MULLO_DC_THRESHOLD, MULLO_MUL_THRESHOLD, Mullo, mullo, mullo_algorithm, mullo_basecase,
    mullo_scratch_len,
};
use crate::kernel::mul::split::{split_wrap_len, splits};
use crate::kernel::mul::{
    Algorithm, MUL_KARATSUBA_THRESHOLD, MUL_NTT_THRESHOLD, MUL_TOOM3_THRESHOLD,
    MUL_TOOM4_THRESHOLD, MUL_TOOM32_THRESHOLD, MUL_TOOM42_THRESHOLD, MUL_WRAP_NTT_THRESHOLD,
    MUL_WRAP_SPLIT_THRESHOLD, NTT_EFFICIENCY_EXPONENT, SQR_KARATSUBA_THRESHOLD, SQR_NTT_THRESHOLD,
    SQR_TOOM3_THRESHOLD, SQR_TOOM4_THRESHOLD, algorithm, karatsuba, karatsuba_scratch_len,
    mul_wrap, mul_wrap_scratch_len, sqr_algorithm, toom3, toom3_scratch_len, toom4,
    toom4_scratch_len, toom32, toom32_fits, toom32_scratch_len, toom42, toom42_fits,
    toom42_scratch_len, wrap_len, wraps_cyclically,
};
use crate::kernel::{
    addmul_1, mul, mul_1, mul_1_assign, mul_basecase, mul_scratch_len, sqr, sqr_basecase,
    sqr_scratch_len,
};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, TESTED_SIZE_MAX, ToRug, differential, tested};

/// Filler for output and scratch buffers, so results can't rely on zeroed memory.
const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn run_mul(a: &[Limb], b: &[Limb]) -> Mag {
    let mut r = vec![GARBAGE; a.len() + b.len()];
    let mut scratch = vec![GARBAGE; mul_scratch_len(a.len(), b.len())];
    mul(&mut r, a, b, &mut scratch);
    Mag::from_limbs(r)
}

fn run_sqr(a: &[Limb]) -> Mag {
    let mut r = vec![GARBAGE; 2 * a.len()];
    let mut scratch = vec![GARBAGE; sqr_scratch_len(a.len())];
    sqr(&mut r, a, &mut scratch);
    Mag::from_limbs(r)
}

fn basecase(a: &[Limb], b: &[Limb]) -> Vec<Limb> {
    let mut r = vec![GARBAGE; a.len() + b.len()];
    mul_basecase(&mut r, a, b);
    r
}

/// Signature of [`karatsuba`], [`toom3`] and [`toom4`].
type Recursive = fn(&mut [Limb], &[Limb], &[Limb], &mut [Limb], bool);

/// Calls `karatsuba`, `toom3` or `toom4` directly with an exactly sized
/// scratch.
fn run_direct(
    f: Recursive,
    scratch_len: fn(usize, usize, bool) -> usize,
    a: &[Limb],
    b: &[Limb],
    square: bool,
) -> Vec<Limb> {
    let mut r = vec![GARBAGE; a.len() + b.len()];
    let mut scratch = vec![GARBAGE; scratch_len(a.len(), b.len(), square)];
    f(&mut r, a, b, &mut scratch, square);
    r
}

/// Thresholds of the multiplication dispatcher, and the sizes where Toom-3
/// recurses into Karatsuba or Toom-3 again.
fn mul_thresholds() -> Vec<usize> {
    vec![
        1,
        2,
        MUL_KARATSUBA_THRESHOLD,
        2 * MUL_KARATSUBA_THRESHOLD,
        MUL_TOOM3_THRESHOLD,
        3 * MUL_KARATSUBA_THRESHOLD,
        3 * MUL_TOOM3_THRESHOLD,
        MUL_TOOM4_THRESHOLD,
        4 * MUL_TOOM4_THRESHOLD,
        MUL_NTT_THRESHOLD,
    ]
}

/// Operands `(a, b)` whose shorter length straddles a threshold and whose
/// longer length is equal, slightly longer, around the Karatsuba and Toom-3
/// balance limits (1.5x, 2x), the unbalanced dispatcher's (4/3, 7/4, 11/4)
/// or several times longer; in either order.
fn mul_operands() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    (len_around(&mul_thresholds()), 0..16usize, any::<bool>())
        .prop_flat_map(|(bn, shape, swap)| {
            let an = match shape {
                0 => bn,
                1 => bn + 1,
                2 => bn * 3 / 2 - 1,
                3 => bn * 3 / 2 + 1,
                4 => 2 * bn - 1,
                5 => 2 * bn,
                6 => 2 * bn + 1,
                7 => 3 * bn + 5,
                // The unbalanced dispatcher's boundaries: 4/3, 7/4, 11/4.
                8 => (4 * bn).div_ceil(3),
                9 => (4 * bn).div_ceil(3) - 1,
                10 => (7 * bn).div_ceil(4),
                11 => (7 * bn).div_ceil(4) - 1,
                12 => (11 * bn).div_ceil(4),
                13 => (11 * bn).div_ceil(4) - 1,
                // Pieces of 2 bn with a short and a long last one.
                14 => 6 * bn + 1,
                _ => 8 * bn - 1,
            }
            .max(bn);
            (strategy::raw_limbs(an), strategy::raw_limbs(bn), Just(swap))
        })
        .prop_map(|(a, b, swap)| if swap { (b, a) } else { (a, b) })
}

/// Lengths `(an, bn)` valid for Karatsuba: `an >= bn > ceil(an / 2)`.
fn karatsuba_lens(max: usize) -> impl Strategy<Value = (usize, usize)> {
    (2..=max).prop_flat_map(|an| (Just(an), an.div_ceil(2) + 1..=an))
}

/// Lengths `(an, bn)` valid for Toom-4: `an >= bn > 3 ceil(an / 4)`.
fn toom4_lens(max: usize) -> impl Strategy<Value = (usize, usize)> {
    (7..=max)
        .prop_filter("an has a valid b", |&an| an > 3 * an.div_ceil(4))
        .prop_flat_map(|an| (Just(an), 3 * an.div_ceil(4) + 1..=an))
}

/// Lengths `(an, bn)` valid for Toom-3: `an >= bn > 2 ceil(an / 3)`.
fn toom3_lens(max: usize) -> impl Strategy<Value = (usize, usize)> {
    (3..=max)
        .prop_filter("an = 4 has no valid b", |&an| an != 4)
        .prop_flat_map(|an| (Just(an), 2 * an.div_ceil(3) + 1..=an))
}

/// Lengths `(an, bn)` that Toom-32 or Toom-42 takes (`fits`), with `an`
/// between `bn` and `3 bn`.
fn unbalanced_lens(
    max: usize,
    fits: fn(usize, usize) -> bool,
) -> impl Strategy<Value = (usize, usize)> {
    (2..=max)
        .prop_flat_map(|bn| (bn..=3 * bn, Just(bn)))
        .prop_filter("the variant takes the shape", move |&(an, bn)| fits(an, bn))
}

/// Signature of [`toom32`] and [`toom42`].
type Unbalanced = fn(&mut [Limb], &[Limb], &[Limb], &mut [Limb]);

/// Calls `toom32` or `toom42` directly with an exactly sized scratch.
fn run_unbalanced(
    f: Unbalanced,
    scratch_len: fn(usize, usize) -> usize,
    a: &[Limb],
    b: &[Limb],
) -> Vec<Limb> {
    let mut r = vec![GARBAGE; a.len() + b.len()];
    let mut scratch = vec![GARBAGE; scratch_len(a.len(), b.len())];
    f(&mut r, a, b, &mut scratch);
    r
}

fn raw_pair(lens: (usize, usize)) -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    (strategy::raw_limbs(lens.0), strategy::raw_limbs(lens.1))
}

differential! {
    fn mul_1_matches(a in (0..=40usize).prop_flat_map(strategy::raw_limbs), b in strategy::limb()) {
        ours: both!(mul_1, mul_1_assign, a, b),
        gmp: split_at_limbs(val(&a) * b, a.len()),
    }

    fn addmul_1_matches(
        (r, a) in (0..=40usize).prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
        b in strategy::limb(),
    ) {
        ours: {
            let mut out = r.clone();
            let c = addmul_1(&mut out, &a, b);
            (Mag::from_limbs(out), c)
        },
        gmp: split_at_limbs(val(&r) + val(&a) * b, r.len()),
    }

    fn mul_basecase_matches(
        (a, b) in (0..=40usize, 0..=40usize)
            .prop_flat_map(|(n, m)| (strategy::raw_limbs(n), strategy::raw_limbs(m))),
    ) {
        ours: Mag::from_limbs(basecase(&a, &b)),
        gmp: val(&a) * val(&b),
    }

    fn sqr_basecase_matches(a in (0..=80usize).prop_flat_map(strategy::raw_limbs)) {
        ours: {
            let mut r = vec![GARBAGE; 2 * a.len()];
            sqr_basecase(&mut r, &a);
            Mag::from_limbs(r)
        },
        gmp: val(&a).square(),
    }

    fn mul_matches((a, b) in mul_operands()) {
        ours: run_mul(&a, &b),
        gmp: val(&a) * val(&b),
    }

    fn sqr_matches(
        a in len_around(&[
            1,
            SQR_KARATSUBA_THRESHOLD,
            SQR_TOOM3_THRESHOLD,
            3 * SQR_TOOM3_THRESHOLD,
            SQR_NTT_THRESHOLD,
        ])
        .prop_flat_map(strategy::raw_limbs),
    ) {
        ours: run_sqr(&a),
        gmp: val(&a).square(),
    }

    fn mul_unbalanced_matches((a, b) in strategy::unbalanced_pair(3 * tested(MUL_TOOM3_THRESHOLD))) {
        ours: run_mul(a.limbs(), b.limbs()),
        gmp: val(a.limbs()) * val(b.limbs()),
    }
}

proptest! {
    #[test]
    fn karatsuba_matches_basecase((a, b) in karatsuba_lens(80).prop_flat_map(raw_pair)) {
        let ours = run_direct(karatsuba, karatsuba_scratch_len, &a, &b, false);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    #[test]
    fn karatsuba_square_matches_basecase(a in (2..=80usize).prop_flat_map(strategy::raw_limbs)) {
        let ours = run_direct(karatsuba, karatsuba_scratch_len, &a, &a, true);
        prop_assert_eq!(ours, basecase(&a, &a));
    }

    #[test]
    fn toom3_matches_basecase((a, b) in toom3_lens(120).prop_flat_map(raw_pair)) {
        let ours = run_direct(toom3, toom3_scratch_len, &a, &b, false);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    #[test]
    fn toom3_square_matches_basecase(
        a in (3..=120usize).prop_filter("", |&n| n != 4).prop_flat_map(strategy::raw_limbs),
    ) {
        let ours = run_direct(toom3, toom3_scratch_len, &a, &a, true);
        prop_assert_eq!(ours, basecase(&a, &a));
    }

    #[test]
    fn toom4_matches_basecase((a, b) in toom4_lens(160).prop_flat_map(raw_pair)) {
        let ours = run_direct(toom4, toom4_scratch_len, &a, &b, false);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    #[test]
    fn toom4_square_matches_basecase(
        a in toom4_lens(160).prop_flat_map(|(an, _)| strategy::raw_limbs(an)),
    ) {
        let ours = run_direct(toom4, toom4_scratch_len, &a, &a, true);
        prop_assert_eq!(ours, basecase(&a, &a));
    }

    /// Toom-4's evaluations at -1 and -2 go negative and the interpolation's
    /// intermediates are largest with all-ones pieces; mixing all-ones and
    /// zero pieces reaches every sign pattern.
    #[test]
    fn toom4_extreme_parts(
        (an, bn) in toom4_lens(80),
        pa in [any::<bool>(); 4],
        pb in [any::<bool>(); 4],
        square in any::<bool>(),
    ) {
        let fill = |len: usize, parts: [bool; 4]| -> Vec<Limb> {
            let n = len.div_ceil(4);
            (0..len).map(|i| if parts[(i / n).min(3)] { Limb::MAX } else { 0 }).collect()
        };
        let a = fill(an, pa);
        let b = if square { a.clone() } else { fill(bn, pb) };
        let ours = run_direct(toom4, toom4_scratch_len, &a, &b, square);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    #[test]
    fn toom32_matches_basecase(
        (a, b) in unbalanced_lens(60, toom32_fits).prop_flat_map(raw_pair),
    ) {
        let ours = run_unbalanced(toom32, toom32_scratch_len, &a, &b);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    #[test]
    fn toom42_matches_basecase(
        (a, b) in unbalanced_lens(60, toom42_fits).prop_flat_map(raw_pair),
    ) {
        let ours = run_unbalanced(toom42, toom42_scratch_len, &a, &b);
        prop_assert_eq!(ours, basecase(&a, &b));
    }

    /// Toom-32's and Toom-42's evaluations at -1 and -2 go negative, and are
    /// largest with all-ones pieces; mixing all-ones and zero pieces of `a`
    /// (thirds or quarters) and `b` (halves) reaches every sign pattern.
    #[test]
    fn toom32_42_extreme_parts(
        (an, bn) in unbalanced_lens(40, |an, bn| toom32_fits(an, bn) || toom42_fits(an, bn)),
        pa in [any::<bool>(); 4],
        pb in [any::<bool>(); 2],
    ) {
        let fill = |len: usize, n: usize, parts: &[bool]| -> Vec<Limb> {
            (0..len).map(|i| if parts[(i / n).min(parts.len() - 1)] { Limb::MAX } else { 0 }).collect()
        };
        if toom32_fits(an, bn) {
            let n = an.div_ceil(3).max(bn.div_ceil(2));
            let (a, b) = (fill(an, n, &pa[..3]), fill(bn, n, &pb));
            prop_assert_eq!(run_unbalanced(toom32, toom32_scratch_len, &a, &b), basecase(&a, &b));
        }
        if toom42_fits(an, bn) {
            let n = an.div_ceil(4).max(bn.div_ceil(2));
            let (a, b) = (fill(an, n, &pa), fill(bn, n, &pb));
            prop_assert_eq!(run_unbalanced(toom42, toom42_scratch_len, &a, &b), basecase(&a, &b));
        }
    }

    /// Toom-3's evaluations go negative; all-ones and all-zero halves push
    /// them to their extremes.
    #[test]
    fn toom3_extreme_parts(
        (an, bn) in toom3_lens(60),
        pa in [any::<bool>(); 3],
        pb in [any::<bool>(); 3],
    ) {
        let fill = |len: usize, parts: [bool; 3]| -> Vec<Limb> {
            let n = len.div_ceil(3);
            (0..len).map(|i| if parts[(i / n).min(2)] { Limb::MAX } else { 0 }).collect()
        };
        let (a, b) = (fill(an, pa), fill(bn, pb));
        let ours = run_direct(toom3, toom3_scratch_len, &a, &b, false);
        prop_assert_eq!(ours, basecase(&a, &b));
    }
}

/// `a * b mod (B^m - 1)` from [`mul_wrap`], checked to be canonical.
fn run_wrap(m: usize, a: &[Limb], b: &[Limb]) -> Integer {
    let mut r = vec![GARBAGE; m];
    let mut scratch = vec![GARBAGE; mul_wrap_scratch_len(m, a.len(), b.len())];
    mul_wrap(&mut r, a, b, &mut scratch);
    let r = val(&r);
    assert!(r < wrap_modulus(m), "not canonical");
    r
}

fn wrap_modulus(m: usize) -> Integer {
    (Integer::from(1) << (64 * u32::try_from(m).unwrap())) - 1u32
}

/// Moduli `B^m - 1`: small ones of every length, powers of two and three
/// times powers of two around the cyclic transform's threshold, and even and
/// odd lengths around the split's threshold, some halving evenly three
/// times.
fn wrap_lens() -> BoxedStrategy<usize> {
    let t = tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two();
    let s = tested(MUL_WRAP_SPLIT_THRESHOLD).max(2);
    prop_oneof![
        1..=40usize,
        Just(t / 2),
        Just(t),
        Just(2 * t),
        Just(t + 1),
        Just(3 * t / 4),
        Just(3 * t / 2),
        Just(3 * t),
        s - 1..=s + 3,
        Just(2 * s + 2),
        Just(split_wrap_len(8 * s + 1)),
        Just(split_wrap_len(5 * s)),
    ]
    .boxed()
}

/// `split_wrap_len` rounds up so that three halvings (fewer near the
/// threshold) stay even, stays below twice its argument, and gives a
/// length that splits for two operands of that argument.
#[test]
fn split_wrap_len_halves_evenly() {
    let s = MUL_WRAP_SPLIT_THRESHOLD.max(2);
    for x in (1..=200).chain([s - 1, s, s + 1, 4 * s + 3, 8 * s + 5]) {
        let m = split_wrap_len(x);
        assert!(m >= x && (x < 2 || m < 2 * x), "x = {x}, m = {m}");
        if x >= 8 * s {
            assert_eq!(m % 8, 0, "x = {x}");
        }
        if x >= s && !wraps_cyclically(m) {
            assert!(splits(m, x, x), "x = {x}, m = {m}");
        }
    }
}

/// The split's edge values modulo `B^m - 1`, `m = 2h`, at split lengths:
/// operands whose residue modulo `B^h + 1` is `B^h` (that is, -1) or zero,
/// residues modulo `B^h - 1` that are zero or all ones, and products that
/// are multiples of either factor or of both.
#[test]
#[allow(
    clippy::many_single_char_names,
    reason = "m and h as in the split's notation"
)]
fn mul_wrap_split_edges() {
    let s = tested(MUL_WRAP_SPLIT_THRESHOLD).max(2);
    for m in [s, s + 1, 2 * s + 2, split_wrap_len(8 * s + 1)] {
        let m = m + m % 2;
        if m > TESTED_SIZE_MAX {
            continue;
        }
        let h = u32::try_from(m / 2).unwrap();
        let bh = Integer::from(1) << (64 * h);
        let modulus = wrap_modulus(m);
        let limbs = |v: &Integer| {
            let mut l = Mag::from_rug(&(v.clone() % &modulus)).limbs().to_vec();
            if l.is_empty() {
                l.push(0);
            }
            l
        };
        let values = [
            Integer::from(&bh - 1u32),
            bh.clone(),
            Integer::from(&bh + 1u32),
            Integer::from(&bh * 2u32),
            Integer::from(&bh * 3u32) - 1u32,
            Integer::from(&modulus - 1u32),
            Integer::from(&modulus - &bh),
            Integer::from(1),
            Integer::from(2),
            Integer::from(&modulus / 3u32),
        ];
        for x in &values {
            for y in &values {
                let (a, b) = (limbs(x), limbs(y));
                // Full-length copies too, so that the split, not the short
                // product, takes them.
                let mut af = a.clone();
                af.resize(m, 0);
                let mut bf = b.clone();
                bf.resize(m, 0);
                let want = val(&a) * val(&b) % &modulus;
                assert_eq!(run_wrap(m, &af, &bf), want, "m = {m}, {x} * {y}");
                assert_eq!(run_wrap(m, &a, &b), want, "m = {m}, {x} * {y}");
            }
        }
    }
}

/// `wrap_len` is the shortest power of two or three times one at least its
/// argument.
#[test]
fn wrap_len_is_shortest_transform_len() {
    for x in 1..=5000usize {
        let m = wrap_len(x);
        assert!(m >= x && crate::kernel::ntt::is_transform_len(m), "{x}");
        assert!(
            (x..m).all(|y| !crate::kernel::ntt::is_transform_len(y)),
            "{x}"
        );
    }
    assert!(wraps_cyclically(
        3 * MUL_WRAP_NTT_THRESHOLD.next_power_of_two()
    ));
}

proptest! {
    /// Operands shorter than the modulus, around it (one piece wraps) and up
    /// to three times longer (folded before a cyclic product), and all-ones
    /// ones, whose residue can come out as `B^m - 1`.
    #[test]
    fn mul_wrap_matches(
        (m, a, b) in wrap_lens()
            .prop_flat_map(|m| (Just(m), 1..=3 * m, prop_oneof![1..=m, m..=3 * m], any::<bool>()))
            .prop_flat_map(|(m, an, bn, ones)| {
                let limbs = move |n: usize| if ones {
                    Just(vec![Limb::MAX; n]).boxed()
                } else {
                    strategy::raw_limbs(n)
                };
                (Just(m), limbs(an), limbs(bn))
            }),
    ) {
        prop_assert_eq!(run_wrap(m, &a, &b), val(&a) * val(&b) % wrap_modulus(m));
    }
}

/// Products that are a multiple of `B^m - 1` come out as zero.
#[test]
fn mul_wrap_multiple_of_modulus() {
    for m in [
        1,
        2,
        3,
        8,
        tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two(),
    ] {
        // B^m - 1 = (B^(m/2) - 1)(B^(m/2) + 1) for even m; 3 divides it.
        let third = Mag::from_rug(&(wrap_modulus(m) / 3u32)).limbs().to_vec();
        assert_eq!(run_wrap(m, &third, &[3]), 0, "m = {m}");
        assert_eq!(run_wrap(m, &[6], &third), 0, "m = {m}");
        let ones = vec![Limb::MAX; m];
        assert_eq!(run_wrap(m, &ones, &ones), 0, "m = {m}");
    }
}

/// Each algorithm, including nested Toom-3 and chunking, on all-ones operands
/// (the worst case for carries). The dispatcher's choice is checked when the
/// configured thresholds are ordered as these shapes assume (as in the
/// default and `stress-small` configurations); shapes beyond test sizes are
/// skipped.
#[test]
#[allow(clippy::many_single_char_names, reason = "one letter per threshold")]
fn every_algorithm_on_all_ones() {
    let k = MUL_KARATSUBA_THRESHOLD;
    let t = MUL_TOOM3_THRESHOLD;
    let q = MUL_TOOM4_THRESHOLD;
    let f = MUL_NTT_THRESHOLD;
    let (p, h) = (MUL_TOOM32_THRESHOLD, MUL_TOOM42_THRESHOLD);
    // A length past the NTT threshold whose transform is tight (a power of
    // two just covering the product), so the NTT takes it whatever the
    // efficiency rule.
    let g = (f + 1).next_power_of_two();
    let cases = [
        (k - 1, k - 1, Algorithm::Basecase),
        (k, k, Algorithm::Karatsuba),
        (k + k / 4, k, Algorithm::Karatsuba),
        (k + k / 2, k, Algorithm::Karatsuba),
        (2 * k, k, Algorithm::Chunked),
        (p + p / 2, p, Algorithm::Toom32),
        (2 * h, h, Algorithm::Toom42),
        (2 * h + 1, h, Algorithm::Toom42),
        (3 * h, h, Algorithm::Chunked),
        (t, t, Algorithm::Toom3),
        (t + t / 2, t, Algorithm::Toom32),
        (q - 1, q - 1, Algorithm::Toom3),
        (q, q, Algorithm::Toom4),
        (4 * q - 1, 3 * q, Algorithm::Toom3),
        (4 * q, 3 * q, Algorithm::Toom32),
        (5 * t + 3, 2 * t, Algorithm::Toom42),
        (7 * t, 2 * t, Algorithm::Chunked),
        (f - 1, f - 1, Algorithm::Toom4),
        (g, g, Algorithm::Ntt),
        (2 * f + 1, f, Algorithm::Ntt),
    ];
    let max = 3 * TESTED_SIZE_MAX;
    let ordered = 2 * k <= t
        && k < p
        && k < h
        && 2 * p <= t
        && 2 * h <= t
        && t < q
        && 3 * q < f
        && 2 * g < max;
    let mut seen = HashSet::new();
    for (an, bn, expected) in cases {
        if an > max {
            continue;
        }
        if ordered {
            assert_eq!(algorithm(an, bn), expected, "({an}, {bn})");
        }
        seen.insert(algorithm(an, bn));
        let (a, b) = (vec![Limb::MAX; an], vec![Limb::MAX; bn]);
        assert_eq!(run_mul(&a, &b).to_rug(), val(&a) * val(&b), "({an}, {bn})");
    }
    if ordered {
        assert_eq!(seen.len(), 8);
    }
    // Around the threshold, where the efficiency rule picks either.
    for (an, bn) in [(f, f), (2 * f - 2, f), (3 * f / 2, f + 7)] {
        if an <= max {
            let (a, b) = (vec![Limb::MAX; an], vec![Limb::MAX; bn]);
            assert_eq!(run_mul(&a, &b).to_rug(), val(&a) * val(&b), "({an}, {bn})");
        }
    }

    for n in [
        SQR_KARATSUBA_THRESHOLD - 1,
        SQR_KARATSUBA_THRESHOLD,
        SQR_TOOM3_THRESHOLD,
        3 * SQR_TOOM3_THRESHOLD + 1,
        SQR_TOOM4_THRESHOLD.max(7),
        4 * SQR_TOOM4_THRESHOLD + 1,
        SQR_NTT_THRESHOLD,
    ] {
        if n > max {
            continue;
        }
        let a = vec![Limb::MAX; n];
        assert_eq!(run_sqr(&a).to_rug(), val(&a).square(), "{n}");
    }
    // Toom-4 squaring takes the lengths with n > 3 ceil(n / 4) (all from 13)
    // from its threshold on; the others stay with Toom-3.
    let (s3, s4, sf) = (SQR_TOOM3_THRESHOLD, SQR_TOOM4_THRESHOLD, SQR_NTT_THRESHOLD);
    let tight = (sf + 1).next_power_of_two();
    let first4 = (s4..s4 + 16)
        .find(|&n| n > 3 * n.div_ceil(4))
        .expect("13 qualifies, and every length from it");
    if SQR_KARATSUBA_THRESHOLD <= s3 && s3 < s4 && first4 < sf {
        assert_eq!(sqr_algorithm(s3), Algorithm::Toom3);
        assert_eq!(sqr_algorithm(s4 - 1), Algorithm::Toom3);
        assert_eq!(sqr_algorithm(first4), Algorithm::Toom4);
        let below = sf - 1;
        let fits = below > 3 * below.div_ceil(4);
        let expected = if fits {
            Algorithm::Toom4
        } else {
            Algorithm::Toom3
        };
        assert_eq!(sqr_algorithm(below), expected);
        assert_eq!(sqr_algorithm(tight), Algorithm::Ntt);
    }
}

/// The NTT's products by chunks of the longer operand, each by the shorter
/// one's prepared transform ([`ntt::chunked`]), against GMP: chunks that
/// divide the longer operand, a last chunk of one limb, and one shorter than
/// the shorter operand, with random and all-ones operands in either order.
/// Called directly, so every configuration runs them.
///
/// [`ntt::chunked`]: crate::kernel::ntt::chunked
#[test]
fn ntt_chunked_matches() {
    use crate::kernel::ntt;
    let mut state = 0x5555_aaaa_1234_5678_u64;
    let mut random = |n: usize| -> Vec<Limb> {
        (0..n)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                state.wrapping_mul(0x2545_f491_4f6c_dd1d)
            })
            .collect()
    };
    let shapes = [
        (640, 64),
        (641, 64),
        (700, 64),
        (8 * 2048, 2048),
        (8 * 2048 + 5, 2048),
    ];
    let mut chunked = 0;
    for (an, bn) in shapes {
        chunked += usize::from(ntt::chunked(an, bn).is_some());
        let (a, b) = (random(an), random(bn));
        for (a, b) in [(a, b), (vec![Limb::MAX; an], vec![Limb::MAX; bn])] {
            for (x, y) in [(&a, &b), (&b, &a)] {
                let mut r = vec![GARBAGE; an + bn];
                let mut scratch = vec![GARBAGE; ntt::mul_scratch_len(x.len(), y.len())];
                ntt::mul(&mut r, x, y, &mut scratch);
                assert_eq!(val(&r), val(&a) * val(&b), "({an}, {bn})");
            }
        }
    }
    assert_eq!(chunked, shapes.len());
}

/// Products of about 10^5 limbs through the dispatcher (the NTT, balanced
/// and unbalanced, and squaring), on pseudo-random and all-ones operands,
/// against GMP.
/// Skipped when the configuration leaves these sizes to quadratic algorithms.
#[test]
fn mul_at_1e5_limbs() {
    if MUL_NTT_THRESHOLD > 30_000 || SQR_NTT_THRESHOLD > 100_000 {
        return;
    }
    let mut state = 0x0123_4567_89ab_cdef_u64;
    let mut random = |n: usize| -> Vec<Limb> {
        (0..n)
            .map(|_| {
                // xorshift64*
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                state.wrapping_mul(0x2545_f491_4f6c_dd1d)
            })
            .collect()
    };
    let a = random(100_000);
    let b = random(99_999);
    let c = random(30_000);
    let ones = vec![Limb::MAX; 100_000];
    for (x, y) in [(&a, &b), (&a, &c), (&ones, &ones), (&ones, &a)] {
        assert_eq!(
            run_mul(x, y).to_rug(),
            val(x) * val(y),
            "({}, {})",
            x.len(),
            y.len()
        );
    }
    assert_eq!(algorithm(a.len(), b.len()), Algorithm::Ntt);
    // One transform of the whole product beats chunks at this size.
    assert_eq!(algorithm(a.len(), c.len()), Algorithm::Ntt);
    for x in [&a, &ones] {
        assert_eq!(run_sqr(x).to_rug(), val(x).square());
    }
}

#[test]
fn empty_operands() {
    assert_eq!(run_mul(&[], &[]), Mag::default());
    assert_eq!(run_mul(&[5, 6], &[]), Mag::default());
    assert_eq!(run_sqr(&[]), Mag::default());
    assert_eq!(mul_1(&mut [], &[], 7), 0);
}

#[test]
fn short_scratch_panics() {
    let n = tested(MUL_TOOM3_THRESHOLD.max(MUL_KARATSUBA_THRESHOLD));
    let need = mul_scratch_len(n, n);
    if need == 0 {
        // Every testable size is schoolbook, which needs no scratch.
        return;
    }
    let (a, b) = (vec![1; n], vec![1; n]);
    let mut r = vec![0; 2 * n];
    let panic = std::panic::catch_unwind(move || {
        mul(&mut r, &a, &b, &mut vec![0; need - 1]);
    })
    .expect_err("a short scratch must panic");
    let message = panic.downcast_ref::<&str>().copied().unwrap_or_default();
    assert!(message.contains("scratch too short"), "{message}");
}

/// Nanoseconds per call of `f`: best of 5 runs of at least 10 ms each.
pub(super) fn time(mut f: impl FnMut()) -> f64 {
    use std::time::{Duration, Instant};
    let mut iters = 1u32;
    loop {
        let start = Instant::now();
        for _ in 0..iters {
            f();
        }
        if start.elapsed() >= Duration::from_millis(10) {
            break;
        }
        iters *= 2;
    }
    (0..5)
        .map(|_| {
            let start = Instant::now();
            for _ in 0..iters {
                f();
            }
            start.elapsed().as_secs_f64() * 1e9 / f64::from(iters)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Threshold tuning and a comparison with GMP. Not a correctness test; run with
/// `cargo test --release -p cas-num tune_mul -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
fn tune_mul() {
    use std::hint::black_box;

    let operand = |n: usize| -> Vec<Limb> {
        (0..n as u64)
            .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0x5555)
            .collect()
    };
    let direct =
        |f: Recursive, scratch_len: fn(usize, usize, bool) -> usize, n: usize, square: bool| {
            let a = operand(n);
            let mut r = vec![0; 2 * n];
            let mut scratch = vec![0; scratch_len(n, n, square)];
            time(|| f(&mut r, black_box(&a), black_box(&a), &mut scratch, square))
        };

    println!("n  basecase  karatsuba | sqr_basecase  karatsuba_sqr   (ns)");
    for n in (8..=64).step_by(4) {
        let a = operand(n);
        let mut r = vec![0; 2 * n];
        let base = time(|| mul_basecase(&mut r, black_box(&a), black_box(&a)));
        let sqr_base = time(|| sqr_basecase(&mut r, black_box(&a)));
        let k = direct(karatsuba, karatsuba_scratch_len, n, false);
        let ks = direct(karatsuba, karatsuba_scratch_len, n, true);
        println!("{n} {base:.0} {k:.0} | {sqr_base:.0} {ks:.0}");
    }

    println!("n  karatsuba  toom3 | karatsuba_sqr  toom3_sqr   (ns)");
    for n in (48..=320).step_by(16) {
        let k = direct(karatsuba, karatsuba_scratch_len, n, false);
        let t = direct(toom3, toom3_scratch_len, n, false);
        let ks = direct(karatsuba, karatsuba_scratch_len, n, true);
        let ts = direct(toom3, toom3_scratch_len, n, true);
        println!("{n} {k:.0} {t:.0} | {ks:.0} {ts:.0}");
    }

    println!("n  ours_mul  gmp_mul  ratio | ours_sqr  gmp_sqr  ratio   (ns, with allocation)");
    for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        let (a, b) = (operand(n), operand(n + 1)[1..].to_vec());
        let (x, y) = (
            crate::Natural::from_limbs(&a),
            crate::Natural::from_limbs(&b),
        );
        let (gx, gy) = (val(&a), val(&b));
        let ours = time(|| drop(black_box(black_box(&x) * black_box(&y))));
        let gmp = time(|| drop(black_box(Integer::from(black_box(&gx) * black_box(&gy)))));
        let ours_sq = time(|| drop(black_box(black_box(&x).square())));
        let gmp_sq = time(|| drop(black_box(Integer::from(black_box(&gx).square_ref()))));
        println!(
            "{n} {ours:.0} {gmp:.0} {:.2} | {ours_sq:.0} {gmp_sq:.0} {:.2}",
            ours / gmp,
            ours_sq / gmp_sq
        );
    }
}

/// The Toom-3 / NTT crossover and a comparison with GMP at NTT sizes. Not a
/// correctness test; run with
/// `cargo test --release -p cas-num tune_ntt -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
#[allow(clippy::many_single_char_names, reason = "as in tune_mul")]
fn tune_ntt() {
    use std::hint::black_box;

    use crate::kernel::ntt;

    let operand = |n: usize, seed: u64| -> Vec<Limb> {
        (0..n as u64)
            .map(|i| (i ^ seed).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0x5555)
            .collect()
    };
    println!("n  toom3  ntt | toom3_sqr  ntt_sqr   (us)");
    for n in [
        256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2304, 2560, 2816, 3072, 3328, 3584, 3840,
        4096, 6144, 8192,
    ] {
        let (a, b) = (operand(n, 1), operand(n, 2));
        let mut r = vec![0; 2 * n];
        let mut s = vec![0; toom3_scratch_len(n, n, false).max(ntt::mul_scratch_len(n, n))];
        let t = time(|| toom3(&mut r, black_box(&a), black_box(&b), &mut s, false)) / 1e3;
        let f = time(|| ntt::mul(&mut r, black_box(&a), black_box(&b), &mut s)) / 1e3;
        let ts = time(|| toom3(&mut r, black_box(&a), black_box(&a), &mut s, true)) / 1e3;
        let fs = time(|| ntt::sqr(&mut r, black_box(&a), &mut s)) / 1e3;
        println!("{n} {t:.1} {f:.1} | {ts:.1} {fs:.1}");
    }

    println!("n  ours_mul  gmp_mul  ratio | ours_sqr  gmp_sqr  ratio   (us, with allocation)");
    for n in [4096, 8192, 16384, 32768, 65536, 100_000, 131_072, 262_144] {
        let (a, b) = (operand(n, 3), operand(n, 4));
        let (x, y) = (
            crate::Natural::from_limbs(&a),
            crate::Natural::from_limbs(&b),
        );
        let (gx, gy) = (val(&a), val(&b));
        let ours = time(|| drop(black_box(black_box(&x) * black_box(&y)))) / 1e3;
        let gmp = time(|| drop(black_box(Integer::from(black_box(&gx) * black_box(&gy))))) / 1e3;
        let ours_sq = time(|| drop(black_box(black_box(&x).square()))) / 1e3;
        let gmp_sq = time(|| drop(black_box(Integer::from(black_box(&gx).square_ref())))) / 1e3;
        println!(
            "{n} {ours:.0} {gmp:.0} {:.2} | {ours_sq:.0} {gmp_sq:.0} {:.2}",
            ours / gmp,
            ours_sq / gmp_sq
        );
    }
}

/// Schoolbook products and squares of all-ones operands, the largest carries
/// through both double-limb sums of a two-row pass, at every small shape: odd
/// and even row counts, the single last row, and squares whose triangle pairs
/// rows down to an empty multiplicand. Squares go up to the length where
/// [`sqr`] leaves schoolbook (at least 24, at most 128 for the stress
/// tunings): every unrolled length, and at each length the largest square,
/// whose doubling pass carries into the top limb.
#[test]
fn basecase_all_ones() {
    for an in 1..=12 {
        let a = vec![Limb::MAX; an];
        for bn in 0..=12 {
            let b = vec![Limb::MAX; bn];
            assert_eq!(val(&basecase(&a, &b)), val(&a) * val(&b), "({an}, {bn})");
        }
    }
    for n in 0..=SQR_KARATSUBA_THRESHOLD.clamp(24, 128) {
        let a = vec![Limb::MAX; n];
        let mut r = vec![GARBAGE; 2 * n];
        sqr_basecase(&mut r, &a);
        assert_eq!(val(&r), val(&a).square(), "{n}");
    }
}

/// From the NTT threshold on, balanced products take the NTT exactly where
/// the efficiency rule says it has caught up with Toom: both outcomes occur
/// between the threshold and twice it, the tight sizes on the NTT side.
#[test]
#[allow(clippy::cast_precision_loss, reason = "the rule's estimate")]
fn ntt_dispatch_follows_efficiency() {
    let f = MUL_NTT_THRESHOLD;
    if f > 1 << 20 {
        return;
    }
    let x = NTT_EFFICIENCY_EXPONENT as f64 / 10.0;
    let (mut ntt, mut toom) = (0, 0);
    for n in f..2 * f {
        let predicted = n as f64 * crate::kernel::ntt::efficiency(n, n).powf(x) >= f as f64;
        let chosen = algorithm(n, n) == Algorithm::Ntt;
        assert_eq!(chosen, predicted, "n = {n}");
        if chosen {
            ntt += 1;
        } else {
            toom += 1;
        }
    }
    if NTT_EFFICIENCY_EXPONENT > 0 && f >= 16 {
        assert!(ntt > 0 && toom > 0, "ntt {ntt}, toom {toom}");
    }
}

/// From both Toom thresholds on, every shape between `4/3 bn` and `11/4 bn`
/// goes to Toom-32, Toom-42 or the NTT: none falls through to the fallbacks
/// for tiny shapes. From `11/4 bn` on, the pieces are taken.
#[test]
fn unbalanced_shapes_have_a_toom() {
    let from = MUL_TOOM32_THRESHOLD.max(MUL_TOOM42_THRESHOLD).max(8);
    if from > TESTED_SIZE_MAX {
        return;
    }
    for bn in from..from + 300 {
        for an in (4 * bn).div_ceil(3)..=4 * bn {
            let chosen = algorithm(an, bn);
            let expected: &[Algorithm] = if 4 * an < 11 * bn {
                &[Algorithm::Toom32, Algorithm::Toom42, Algorithm::Ntt]
            } else {
                &[Algorithm::Chunked, Algorithm::Ntt]
            };
            assert!(expected.contains(&chosen), "({an}, {bn}): {chosen:?}");
        }
    }
}

/// `a b mod B^n` from [`mullo`] with an exactly sized, garbage-filled scratch.
fn run_mullo(a: &[Limb], b: &[Limb]) -> Vec<Limb> {
    let mut r = vec![GARBAGE; a.len()];
    let mut scratch = vec![GARBAGE; mullo_scratch_len(a.len())];
    mullo(&mut r, a, b, &mut scratch);
    r
}

/// The low `n` limbs of `a b`, by schoolbook.
fn low_half(a: &[Limb], b: &[Limb]) -> Vec<Limb> {
    let mut full = basecase(a, b);
    full.truncate(a.len());
    full
}

/// Lengths around both short-product thresholds, and where the divide and
/// conquer split recurses into its own range.
fn mullo_lens() -> BoxedStrategy<usize> {
    len_around(&[
        1,
        2,
        MULLO_DC_THRESHOLD,
        2 * MULLO_DC_THRESHOLD,
        4 * MULLO_DC_THRESHOLD,
        MULLO_MUL_THRESHOLD,
    ])
    .boxed()
}

proptest! {
    #[test]
    fn mullo_matches(
        (a, b) in mullo_lens().prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
    ) {
        prop_assert_eq!(run_mullo(&a, &b), low_half(&a, &b));
    }

    #[test]
    fn mullo_basecase_matches(
        (a, b) in (0..=40usize).prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
    ) {
        let mut r = vec![GARBAGE; a.len()];
        mullo_basecase(&mut r, &a, &b);
        prop_assert_eq!(r, low_half(&a, &b));
    }
}

/// All-ones operands (every carry out of `B^n` set) at each algorithm, and the
/// dispatcher's choice.
#[test]
fn mullo_all_ones() {
    for n in [
        1,
        2,
        3,
        MULLO_DC_THRESHOLD,
        MULLO_DC_THRESHOLD + 1,
        3 * MULLO_DC_THRESHOLD + 1,
        MULLO_MUL_THRESHOLD,
    ] {
        let n = tested(n);
        let ones = vec![Limb::MAX; n];
        assert_eq!(run_mullo(&ones, &ones), low_half(&ones, &ones), "n = {n}");
    }
    if MULLO_DC_THRESHOLD < MULLO_MUL_THRESHOLD && MULLO_MUL_THRESHOLD <= TESTED_SIZE_MAX {
        assert_eq!(mullo_algorithm(MULLO_DC_THRESHOLD - 1), Mullo::Basecase);
        assert_eq!(mullo_algorithm(MULLO_DC_THRESHOLD), Mullo::Dc);
        assert_eq!(mullo_algorithm(MULLO_MUL_THRESHOLD), Mullo::Full);
    }
}
