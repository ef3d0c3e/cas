//! Multiplication kernels: checked against GMP through the size dispatcher, and
//! Karatsuba and Toom-3 called directly, checked against schoolbook (the
//! independent second oracle) at every small split shape.

use std::collections::HashSet;

use proptest::prelude::*;
use rug::Integer;

use super::{split_at_limbs, val};
use crate::Limb;
use crate::kernel::mul::{
    Algorithm, MUL_KARATSUBA_THRESHOLD, MUL_NTT_THRESHOLD, MUL_TOOM3_THRESHOLD,
    MUL_TOOM4_THRESHOLD, MUL_WRAP_NTT_THRESHOLD, SQR_KARATSUBA_THRESHOLD, SQR_NTT_THRESHOLD,
    SQR_TOOM3_THRESHOLD, SQR_TOOM4_THRESHOLD, algorithm, karatsuba, karatsuba_scratch_len,
    mul_wrap, mul_wrap_scratch_len, sqr_algorithm, toom3, toom3_scratch_len, toom4,
    toom4_scratch_len,
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
/// balance limits (1.5x, 2x) or several times longer; in either order.
fn mul_operands() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    (len_around(&mul_thresholds()), 0..8usize, any::<bool>())
        .prop_flat_map(|(bn, shape, swap)| {
            let an = match shape {
                0 => bn,
                1 => bn + 1,
                2 => bn * 3 / 2 - 1,
                3 => bn * 3 / 2 + 1,
                4 => 2 * bn - 1,
                5 => 2 * bn,
                6 => 2 * bn + 1,
                _ => 3 * bn + 5,
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

    fn sqr_basecase_matches(a in (0..=40usize).prop_flat_map(strategy::raw_limbs)) {
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

/// Moduli `B^m - 1`: small ones of every length, and powers of two around the
/// cyclic transform's threshold.
fn wrap_lens() -> BoxedStrategy<usize> {
    let t = tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two();
    prop_oneof![1..=40usize, Just(t / 2), Just(t), Just(2 * t), Just(t + 1)].boxed()
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
    let cases = [
        (k - 1, k - 1, Algorithm::Basecase),
        (k, k, Algorithm::Karatsuba),
        (2 * k, k, Algorithm::Chunked),
        (2 * k + 1, k, Algorithm::Chunked),
        (t, t, Algorithm::Toom3),
        (t + t / 2, t, Algorithm::Karatsuba),
        (q - 1, q - 1, Algorithm::Toom3),
        (q, q, Algorithm::Toom4),
        (4 * q, 3 * q, Algorithm::Toom3),
        (5 * t + 3, 2 * t, Algorithm::Chunked),
        (f - 1, f - 1, Algorithm::Toom4),
        (f, f, Algorithm::Ntt),
        (2 * f - 2, f, Algorithm::Ntt),
        (2 * f + 1, f, Algorithm::Chunked),
    ];
    let max = 3 * TESTED_SIZE_MAX;
    let ordered = 2 * k <= t && t < q && 3 * q < f && 2 * f < max;
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
        assert_eq!(seen.len(), 6);
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
        assert_eq!(sqr_algorithm(sf), Algorithm::Ntt);
    }
}

/// Products of about 10^5 limbs through the dispatcher (the NTT, chunked
/// NTT, and squaring), on pseudo-random and all-ones operands, against GMP.
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
    assert_eq!(algorithm(a.len(), c.len()), Algorithm::Chunked);
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
/// rows down to an empty multiplicand.
#[test]
fn basecase_all_ones() {
    for an in 1..=12 {
        let a = vec![Limb::MAX; an];
        for bn in 0..=12 {
            let b = vec![Limb::MAX; bn];
            assert_eq!(val(&basecase(&a, &b)), val(&a) * val(&b), "({an}, {bn})");
        }
    }
    for n in 0..=24 {
        let a = vec![Limb::MAX; n];
        let mut r = vec![GARBAGE; 2 * n];
        sqr_basecase(&mut r, &a);
        assert_eq!(val(&r), val(&a).square(), "{n}");
    }
}
