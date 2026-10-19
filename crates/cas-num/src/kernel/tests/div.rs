//! Division kernels checked against GMP, including inputs built to reach the
//! rare quotient-correction branches (checked with the branch counters).
#![allow(
    clippy::many_single_char_names,
    reason = "textbook names, as in the kernels"
)]
#![allow(
    clippy::cast_possible_truncation,
    reason = "the small-word tests truncate random values on purpose"
)]

use std::cell::Cell;
use std::thread::LocalKey;

use proptest::prelude::*;
use rug::Complete;
use rug::Integer;
use rug::ops::{DivRounding, RemRounding};

use super::{base_pow, val};
use crate::Limb;
use crate::kernel::div::bdiv::{
    BDIV_MU_THRESHOLD, BDIV_Q_DC_THRESHOLD, BDIV_QR_DC_THRESHOLD, BINV_NEWTON_THRESHOLD,
    DIVEXACT_DIV_THRESHOLD, bdiv_prepared_divisor, bdiv_prepared_inverse, binvert,
    binvert_scratch_len, binvert_shared,
};
use crate::kernel::div::mu::{
    INV_NEWTON_THRESHOLD, MU_DIV_THRESHOLD, block_inverse, block_len, folds_blocks, invert,
    invert_scratch_len, mu_div_qr, mu_div_qr_scratch_len, mu_divappr_q, mu_divappr_q_scratch_len,
    prepares, rem_wrap_len,
};
use crate::kernel::div::{
    DC_DIV_THRESHOLD, bz_div_qr, bz_div_qr_scratch_len, bz_divappr_q, div_2by1, div_3by2,
    divappr_q, divappr_q_scratch_len, hits, reciprocal_2by1, reciprocal_3by2, sb_div_qr,
};
use crate::kernel::{
    bdiv_q, bdiv_q_scratch_len, div_rem, div_rem_1, div_rem_scratch_len, divexact, divexact_1,
    divexact_scratch_len, rem_1, shl, submul_1,
};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, ToRug, differential, tested};

const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn normalised_limb() -> impl Strategy<Value = Limb> {
    strategy::limb().prop_map(|x| x | (1 << 63))
}

fn run_div(a: &[Limb], d: &[Limb]) -> (Mag, Mag) {
    let (an, dn) = (a.len(), d.len());
    let mut q = vec![GARBAGE; an - dn + 1];
    let mut r = vec![GARBAGE; dn];
    let mut scratch = vec![GARBAGE; div_rem_scratch_len(an, dn)];
    div_rem(&mut q, &mut r, a, d, &mut scratch);
    (Mag::from_limbs(q), Mag::from_limbs(r))
}

fn gmp_div(a: &[Limb], d: &[Limb]) -> (Integer, Integer) {
    val(a).div_rem(val(d))
}

/// Limbs of a non-negative GMP value.
fn limbs_of(x: &Integer) -> Vec<Limb> {
    Mag::from_rug(x).limbs().to_vec()
}

/// Divisor lengths straddling the divide-and-conquer threshold and the sizes
/// where its recursion splits again.
fn divisor_lens() -> BoxedStrategy<usize> {
    len_around(&[
        1,
        2,
        3,
        DC_DIV_THRESHOLD,
        2 * DC_DIV_THRESHOLD,
        4 * DC_DIV_THRESHOLD,
        tested(MU_DIV_THRESHOLD),
    ])
}

/// Quotient lengths: none, one limb, around the divisor length (one block,
/// a partial block) and several blocks.
fn quotient_len(dn: usize) -> impl Strategy<Value = usize> {
    prop_oneof![
        Just(0),
        Just(1),
        Just(2),
        Just(dn.saturating_sub(1)),
        Just(dn),
        Just(dn + 1),
        Just(2 * dn + 3),
        0..=dn,
    ]
}

/// `(a, d)` for the kernel: `d` with a non-zero top limb (normalised or not),
/// `a` either raw limbs or `q d + e` for a small offset `e`, which puts the
/// quotient estimates right at their correction boundary.
fn div_operands() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    divisor_lens()
        .prop_flat_map(|dn| (Just(dn), quotient_len(dn)))
        .prop_flat_map(|(dn, qn)| {
            let d = strategy::natural_with_len(Just(dn));
            let raw = strategy::raw_limbs(dn + qn).prop_map(|a| (a, None)).boxed();
            let shaped = (
                strategy::natural_with_len(Just(qn.max(1))),
                prop_oneof![Just(-1i64), Just(0), Just(1), Just(-2), Just(i64::MAX)],
            )
                .prop_map(|(q, e)| (Vec::new(), Some((q, e))))
                .boxed();
            (d, prop_oneof![raw, shaped])
        })
        .prop_map(|(d, (raw, shaped))| {
            let d = d.limbs().to_vec();
            let a = match shaped {
                None => raw,
                Some((q, e)) => {
                    // e = i64::MAX stands for d - 1, the largest remainder.
                    let e = if e == i64::MAX {
                        val(&d) - 1u32
                    } else {
                        Integer::from(e)
                    };
                    let a = q.to_rug() * val(&d) + e;
                    let mut a = limbs_of(&a.max(Integer::new()));
                    a.resize(a.len().max(d.len()), 0);
                    a
                }
            };
            (a, d)
        })
}

/// `(u, d)` for `div_3by2`: `d` normalised and `(u2, u1) < (d1, d0)`.
fn div_3by2_operands() -> impl Strategy<Value = ([Limb; 3], [Limb; 2])> {
    (
        [strategy::limb(), strategy::limb(), strategy::limb()],
        (strategy::limb(), normalised_limb()),
        any::<bool>(),
    )
        .prop_map(|(u, (d0, d1), near)| {
            let d = [d0, d1];
            // Either reduce (u2, u1) mod d, or take it just below d.
            let top = if near {
                val(&d) - 1u32
            } else {
                val(&u[1..]) % val(&d)
            };
            let mut top = limbs_of(&top);
            top.resize(2, 0);
            ([u[0], top[0], top[1]], d)
        })
}

differential! {
    fn reciprocal_2by1_matches(d in normalised_limb()) {
        ours: reciprocal_2by1(d),
        gmp: (Integer::from(u128::MAX) / d) - base_pow(1),
    }

    fn reciprocal_3by2_matches(d1 in normalised_limb(), d0 in strategy::limb()) {
        ours: reciprocal_3by2(d1, d0),
        gmp: (base_pow(3) - 1u32) / val(&[d0, d1]) - base_pow(1),
    }

    fn div_2by1_matches(d in normalised_limb(), u1 in strategy::limb(), u0 in strategy::limb()) {
        ours: {
            let u1 = u1 % d;
            div_2by1(u1, u0, d, reciprocal_2by1(d))
        },
        gmp: val(&[u0, u1 % d]).div_rem(Integer::from(d)),
    }

    fn div_3by2_matches((u, d) in div_3by2_operands()) {
        ours: div_3by2(u[2], u[1], u[0], d[1], d[0], reciprocal_3by2(d[1], d[0])),
        gmp: val(&u).div_rem(val(&d)),
    }

    fn div_rem_1_matches(a in (0..=40usize).prop_flat_map(strategy::raw_limbs), d in strategy::top_limb()) {
        ours: {
            let mut q = vec![GARBAGE; a.len()];
            let r = div_rem_1(&mut q, &a, d);
            prop_assert_eq!(r, rem_1(&a, d), "rem_1 differs");
            (Mag::from_limbs(q), r)
        },
        gmp: val(&a).div_rem(Integer::from(d)),
    }

    fn divexact_1_matches(a in (0..=40usize).prop_flat_map(strategy::raw_limbs), d in strategy::top_limb()) {
        ours: {
            let mut q = vec![GARBAGE; a.len()];
            divexact_1(&mut q, &a, d).then(|| Mag::from_limbs(q))
        },
        gmp: val(&a).is_divisible(&Integer::from(d)).then(|| val(&a) / d),
    }

    /// Multiples of `d`, where the quotient must come out right.
    fn divexact_1_exact(q in (1..=40usize).prop_flat_map(strategy::raw_limbs), d in strategy::top_limb()) {
        ours: {
            let a = limbs_of(&(val(&q) * d));
            let mut out = vec![GARBAGE; a.len()];
            divexact_1(&mut out, &a, d).then(|| Mag::from_limbs(out))
        },
        gmp: Some(val(&q)),
    }

    fn submul_1_matches(
        (r, a) in (0..=40usize).prop_flat_map(|n| (strategy::raw_limbs(n), strategy::raw_limbs(n))),
        b in strategy::limb(),
    ) {
        ours: {
            let mut out = r.clone();
            let c = submul_1(&mut out, &a, b);
            (Mag::from_limbs(out), c)
        },
        gmp: {
            let x = val(&r) - val(&a) * b;
            let low = x.clone().rem_euc(base_pow(r.len()));
            let borrow = (Integer::from(&low - &x)).div_euc(base_pow(r.len()));
            (low, borrow)
        },
    }

    /// Algorithm D on its own, with its full output contract: quotient, top
    /// quotient limb, remainder in the low limbs and zeros above.
    fn sb_div_qr_matches(
        (n, d) in (2..=12usize, 0..=12usize).prop_flat_map(|(m, qn)| (
            strategy::raw_limbs(m + qn),
            (strategy::raw_limbs(m - 1), normalised_limb()),
        )),
    ) {
        ours: {
            let (mut d, top) = d.clone();
            d.push(top);
            let mut work = n.clone();
            let mut q = vec![GARBAGE; n.len() - d.len()];
            let qh = sb_div_qr(&mut q, &mut work, &d, reciprocal_3by2(top, d[d.len() - 2]));
            prop_assert!(work[d.len()..].iter().all(|&x| x == 0), "limbs above the remainder");
            q.push(qh);
            (Mag::from_limbs(q), Mag::from_limbs(&work[..d.len()]))
        },
        gmp: {
            let (mut d, top) = d.clone();
            d.push(top);
            gmp_div(&n, &d)
        },
    }

    fn div_rem_matches((a, d) in div_operands()) {
        ours: run_div(&a, &d),
        gmp: gmp_div(&a, &d),
    }
}

proptest! {
    /// Division undoes multiplication, with a remainder below the divisor.
    #[test]
    fn div_rem_identity((a, d) in div_operands()) {
        let (q, r) = run_div(&a, &d);
        prop_assert!(r.to_rug() < val(&d));
        prop_assert_eq!(q.to_rug() * val(&d) + r.to_rug(), val(&a));
    }
}

fn count(counter: &'static LocalKey<Cell<u64>>, f: impl FnOnce()) -> u64 {
    counter.with(|c| c.set(0));
    f();
    counter.with(Cell::get)
}

/// Checks `a / d` against GMP through the schoolbook kernel directly, whatever the
/// configured divide-and-conquer threshold; `d` normalised, `d.len() >= 2`.
fn check_schoolbook(a: &[Limb], d: &[Limb]) -> Mag {
    let dn = d.len();
    let mut n = a.to_vec();
    let mut q = vec![GARBAGE; a.len() - dn];
    let v = reciprocal_3by2(d[dn - 1], d[dn - 2]);
    let qh = sb_div_qr(&mut q, &mut n, d, v);
    q.push(qh);
    let (q, r) = (Mag::from_limbs(q), Mag::from_limbs(n[..dn].to_vec()));
    assert_eq!((q.to_rug(), r.to_rug()), gmp_div(a, d));
    q
}

/// Algorithm D's add-back step: the 3-by-2 estimate is one too large. This is
/// the 64-bit version of the Hacker's Delight example for `divmnu`.
#[test]
fn add_back_is_reached() {
    let a = [0, 0, 1 << 63, (1 << 63) - 1];
    let d = [1, 0, 1 << 63];
    let mut q = Mag::default();
    assert_eq!(count(&hits::ADD_BACK, || q = check_schoolbook(&a, &d)), 1);
    assert_eq!(q, Mag::from_u64(u64::MAX - 1));

    // In general a = q d - 1 with low divisor limbs makes the last estimate overshoot.
    let d = [5, 7, u64::MAX, 1 << 63];
    let a = limbs_of(&(Integer::from(0x1234_5678_9abc_u64) * val(&d) - 1u32));
    assert!(count(&hits::ADD_BACK, || drop(check_schoolbook(&a, &d))) >= 1);
}

/// Top window limbs equal to the divisor's: the estimate is capped at B - 1.
#[test]
fn capped_estimate_is_reached() {
    let d = [3, u64::MAX, 1 << 63];
    // a = B d - 1, so the quotient is B - 1.
    let a = limbs_of(&(val(&d) * base_pow(1) - 1u32));
    let mut q = Mag::default();
    assert_eq!(count(&hits::CAPPED, || q = check_schoolbook(&a, &d)), 1);
    assert_eq!(q, Mag::from_u64(u64::MAX));
}

/// Checks `a / d` against GMP through divide and conquer directly, whatever the
/// configured Barrett threshold; `d` normalised, `d.len() >= 2`.
fn check_bz(a: &[Limb], d: &[Limb]) {
    let dn = d.len();
    let mut n = a.to_vec();
    let mut q = vec![GARBAGE; a.len() - dn];
    let mut scratch = vec![GARBAGE; bz_div_qr_scratch_len(a.len(), dn)];
    let v = reciprocal_3by2(d[dn - 1], d[dn - 2]);
    let qh = bz_div_qr(&mut q, &mut n, d, v, &mut scratch);
    q.push(qh);
    let (q, r) = (Mag::from_limbs(q), Mag::from_limbs(n[..dn].to_vec()));
    assert_eq!((q.to_rug(), r.to_rug()), gmp_div(a, d));
}

/// The divide-and-conquer block estimates overshoot for `a = q d - 1`.
#[test]
fn dc_correction_is_reached() {
    if DC_DIV_THRESHOLD > crate::testkit::TESTED_SIZE_MAX {
        return;
    }
    // Divide and conquer needs two quotient limbs, and its blocks only
    // correct when they have at least two limbs (dn >= 4).
    let dn0 = DC_DIV_THRESHOLD.max(4);
    for dn in [dn0, 2 * dn0 + 1] {
        for qn in [(dn / 2).max(2), dn, 2 * dn + 1] {
            let d: Vec<Limb> = (0..dn as u64)
                .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15) | (1 << 63))
                .collect();
            let q: Vec<Limb> = (0..qn as u64).map(|i| !i).collect();
            let a = limbs_of(&(val(&q) * val(&d) - 1u32));
            assert!(
                count(&hits::DC_CORRECTION, || check_bz(&a, &d)) >= 1,
                "dn = {dn}, qn = {qn}"
            );
        }
    }
}

#[test]
fn single_limb_edge_cases() {
    assert_eq!(div_rem_1(&mut [], &[], 7), 0);
    assert_eq!(rem_1(&[], 7), 0);
    assert!(divexact_1(&mut [], &[], 7));
    let mut q = [0; 2];
    assert_eq!(div_rem_1(&mut q, &[u64::MAX, u64::MAX], 1), 0);
    assert_eq!(q, [u64::MAX, u64::MAX]);
    // Even divisor, odd dividend: rejected by the low bits alone.
    assert!(!divexact_1(&mut q, &[3, 0], 2));
}

#[test]
#[should_panic(expected = "zero top limb")]
fn unnormalised_divisor_panics() {
    let mut q = [0; 2];
    let mut r = [0; 2];
    div_rem(&mut q, &mut r, &[1, 2, 3], &[1, 0], &mut []);
}

#[test]
#[should_panic(expected = "division by zero")]
fn division_by_zero_limb_panics() {
    let _ = rem_1(&[1], 0);
}

/// xorshift64*, for reproducible bulk sampling outside proptest.
fn rng(seed: u64) -> impl FnMut() -> u64 {
    let mut s = seed;
    move || {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        s.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
}

/// The Möller–Granlund primitives on 8-bit words, exhaustively where
/// feasible: every boundary equality (probability ~2^-64 for 64-bit limbs)
/// occurs here.
#[test]
fn word8_primitives_exhaustive() {
    use crate::kernel::div::word8::{div_2by1, div_3by2, reciprocal_2by1, reciprocal_3by2};
    for d in 0x80..=u8::MAX {
        let v = reciprocal_2by1(d);
        assert_eq!(u16::from(v), (u16::MAX / u16::from(d)) - 0x100, "d = {d}");
        for u1 in 0..d {
            for u0 in 0..=u8::MAX {
                let u = u16::from(u1) << 8 | u16::from(u0);
                let expected = (u / u16::from(d), u % u16::from(d));
                let got = div_2by1(u1, u0, d, v);
                assert_eq!((u16::from(got.0), u16::from(got.1)), expected, "{u} / {d}");
            }
        }
    }
    let mut next = rng(1);
    for d1 in 0x80..=u8::MAX {
        for d0 in 0..=u8::MAX {
            let d = u32::from(d1) << 8 | u32::from(d0);
            let v = reciprocal_3by2(d1, d0);
            assert_eq!(u32::from(v), (0xff_ffff / d) - 0x100, "d = {d}");
            // (u2, u1) at both ends of its range, then random below d.
            let tops: Vec<u32> = [0, d - 1]
                .into_iter()
                .chain((0..62).map(|_| next() as u32 % d))
                .collect();
            for top in tops {
                let u0 = next() as u8;
                let u = top << 8 | u32::from(u0);
                let (q, r) = div_3by2((top >> 8) as u8, top as u8, u0, d1, d0, v);
                assert_eq!((u32::from(q), u32::from(r)), (u / d, u % d), "{u} / {d}");
            }
        }
    }
}

/// The same primitives on 16-bit words, sampled.
#[test]
fn word16_primitives_sampled() {
    use crate::kernel::div::word16::{div_2by1, div_3by2, reciprocal_2by1, reciprocal_3by2};
    let mut next = rng(2);
    for _ in 0..1 << 20 {
        let (d1, d0) = (next() as u16 | 0x8000, next() as u16);
        let v = reciprocal_2by1(d1);
        assert_eq!(u32::from(v), (u32::MAX / u32::from(d1)) - 0x1_0000);
        let (u1, u0) = (next() as u16 % d1, next() as u16);
        let u = u32::from(u1) << 16 | u32::from(u0);
        let (q, r) = div_2by1(u1, u0, d1, v);
        assert_eq!(
            (u32::from(q), u32::from(r)),
            (u / u32::from(d1), u % u32::from(d1))
        );

        let d = u64::from(d1) << 16 | u64::from(d0);
        let v = reciprocal_3by2(d1, d0);
        assert_eq!(u64::from(v), ((1 << 48) - 1) / d - 0x1_0000);
        let top = if next() % 8 == 0 { d - 1 } else { next() % d };
        let u0 = next() as u16;
        let u = top << 16 | u64::from(u0);
        let (q, r) = div_3by2((top >> 16) as u16, top as u16, u0, d1, d0, v);
        assert_eq!((u64::from(q), u64::from(r)), (u / d, u % d));
    }
}

/// [`dc_div_qr_n_scratch_len`]'s pass over the depths against the recursion
/// it stands for: lengths whose halves fall on both sides of
/// `DC_DIV_THRESHOLD` and around powers of two, with products of every
/// algorithm.
///
/// [`dc_div_qr_n_scratch_len`]: crate::kernel::div::dc_div_qr_n_scratch_len
#[test]
fn dc_div_scratch_matches_recursion() {
    use crate::kernel::div::dc_div_qr_n_scratch_len;
    use crate::kernel::mul::mul_scratch_len;

    fn block(k: usize, dn: usize) -> usize {
        if k == 1 {
            return 0;
        }
        let inner = if k < DC_DIV_THRESHOLD { 0 } else { dc(k) };
        let correction = if k == dn {
            0
        } else {
            dn + mul_scratch_len(k, dn - k)
        };
        inner.max(correction)
    }
    fn dc(k: usize) -> usize {
        block(k - k / 2, k).max(block(k / 2, k))
    }
    let lens = (2..=400).chain([
        511, 512, 513, 1000, 1023, 1025, 2047, 2049, 3001, 4097, 8191,
    ]);
    for k in lens.filter(|&k| k >= DC_DIV_THRESHOLD) {
        assert_eq!(dc_div_qr_n_scratch_len(k), dc(k), "k = {k}");
    }
}

/// Threshold tuning and a comparison with GMP. Not a correctness test; run with
/// `cargo test --release -p cas-num tune_div -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
fn tune_div() {
    use std::hint::black_box;

    use super::mul::time;
    use crate::kernel::div::{dc_div_qr_n, dc_div_qr_n_scratch_len};

    let mut next = rng(3);
    let mut operand = |n: usize| -> Vec<Limb> { (0..n).map(|_| next()).collect() };

    println!("n  schoolbook  divide_and_conquer   (ns, 2n by n limbs)");
    for n in (16..=160).step_by(8) {
        let mut d = operand(n);
        d[n - 1] |= 1 << 63;
        let a = operand(2 * n);
        let v = reciprocal_3by2(d[n - 1], d[n - 2]);
        let mut w = a.clone();
        let mut q = vec![0; n];
        let sb = time(|| {
            w.copy_from_slice(&a);
            sb_div_qr(&mut q, &mut w, black_box(&d), v);
        });
        let mut scratch = vec![0; dc_div_qr_n_scratch_len(n)];
        let dc = time(|| {
            w.copy_from_slice(&a);
            dc_div_qr_n(&mut q, &mut w, black_box(&d), v, &mut scratch);
        });
        println!("{n} {sb:.0} {dc:.0}");
    }

    println!("dn  qn  ours  gmp  ratio   (ns, with allocation)");
    for (dn, qn) in [
        (1, 1),
        (1, 16),
        (1, 256),
        (2, 2),
        (8, 8),
        (32, 32),
        (128, 128),
        (512, 512),
        (2048, 2048),
        (100, 1000),
        (1000, 100),
    ] {
        let mut d = operand(dn);
        d[dn - 1] |= 1;
        let a = operand(dn + qn);
        let (x, y) = (
            crate::Natural::from_limbs(&a),
            crate::Natural::from_limbs(&d),
        );
        let (gx, gy) = (val(&a), val(&d));
        let ours = time(|| drop(black_box(black_box(&x).div_rem(black_box(&y)))));
        let gmp = time(|| {
            drop(black_box(
                black_box(&gx).div_rem_ref(black_box(&gy)).complete(),
            ));
        });
        println!("{dn} {qn} {ours:.0} {gmp:.0} {:.2}", ours / gmp);
    }
}

/// `X = B^n + I` from [`invert`] on the normalised `a`, checked against
/// `A X < B^(2n) < A (X + 2)`.
fn check_invert(a: &[Limb]) -> Integer {
    let n = a.len();
    let mut x = vec![GARBAGE; n];
    let mut scratch = vec![GARBAGE; invert_scratch_len(n)];
    invert(&mut x, a, &mut scratch);
    let x = val(&x) + base_pow(n);
    let (a, b2n) = (val(a), base_pow(2 * n));
    assert!(Integer::from(&a * &x) < b2n, "A X >= B^2n for n = {n}");
    assert!(
        (x.clone() + 2u32) * a > b2n,
        "A (X + 2) <= B^2n for n = {n}"
    );
    x
}

/// Reciprocal lengths: every small one, and around the Newton threshold and
/// its first recursion, `2 (t - 1)`.
fn invert_lens() -> BoxedStrategy<usize> {
    let t = tested(INV_NEWTON_THRESHOLD);
    prop_oneof![len_around(&[1, 2, 3, t, 2 * t - 2, 4 * t]), 1..=40usize].boxed()
}

/// Normalised values with the shapes that stress a reciprocal: random, the
/// smallest (`B^n / 2`, reciprocal `2 B^n - 1`), the largest, and runs of
/// ones or zeros below a random top.
fn invert_operand(n: usize) -> BoxedStrategy<Vec<Limb>> {
    let top = 1 << 63;
    prop_oneof![
        strategy::raw_limbs(n).prop_map(move |mut a| {
            a[n - 1] |= top;
            a
        }),
        Just({
            let mut a = vec![0; n];
            a[n - 1] = top;
            a
        }),
        Just(vec![Limb::MAX; n]),
        (normalised_limb(), any::<bool>(), 0..n).prop_map(move |(t, ones, k)| {
            let mut a = vec![if ones { Limb::MAX } else { 0 }; n];
            a[n - 1] = t;
            a[k.min(n - 1)] ^= 1;
            a
        }),
    ]
    .boxed()
}

proptest! {
    #[test]
    fn invert_bound(a in invert_lens().prop_flat_map(invert_operand)) {
        check_invert(&a);
    }
}

/// Checks `a / d` against GMP through Barrett division directly, whatever the configured
/// thresholds; `d` normalised, `d.len() >= 2`, `a.len() > d.len()`.
fn check_mu(a: &[Limb], d: &[Limb]) {
    let dn = d.len();
    let mut n = a.to_vec();
    let mut q = vec![GARBAGE; a.len() - dn];
    let mut scratch = vec![GARBAGE; mu_div_qr_scratch_len(a.len(), dn)];
    let qh = mu_div_qr(&mut q, &mut n, d, &mut scratch);
    q.push(qh);
    let (q, r) = (Mag::from_limbs(q), Mag::from_limbs(n[..dn].to_vec()));
    assert_eq!(
        (q.to_rug(), r.to_rug()),
        gmp_div(a, d),
        "dn = {dn}, an = {}",
        a.len()
    );
}

/// `(a, d)` for Barrett division: `d` normalised, of a length around the
/// thresholds of the reciprocal it inverts, with `a` as in [`div_operands`].
fn mu_operands() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    let t = tested(INV_NEWTON_THRESHOLD);
    (len_around(&[2, 3, 4, 8, t, 2 * t, 3 * t]), any::<bool>())
        .prop_flat_map(|(dn, small)| {
            let dn = dn.max(2);
            let qn = prop_oneof![
                Just(1),
                Just(dn / 3),
                Just(dn / 3 + 1),
                Just(dn),
                Just(dn + 1),
                Just(2 * dn + 3),
                1..=2 * dn
            ];
            (
                invert_operand(dn),
                qn.prop_map(|q: usize| q.max(1)),
                Just(small),
            )
        })
        .prop_flat_map(|(d, qn, small)| {
            let dn = d.len();
            let raw = strategy::raw_limbs(dn + qn).prop_map(|a| (a, None)).boxed();
            let shaped = (
                strategy::natural_with_len(Just(qn)),
                prop_oneof![Just(0i64), Just(1), Just(i64::MAX)],
            )
                .prop_map(|(q, e)| (Vec::new(), Some((q, e))))
                .boxed();
            (Just(d), if small { raw } else { shaped })
        })
        .prop_map(|(d, (raw, shaped))| {
            let a = match shaped {
                None => raw,
                Some((q, e)) => {
                    let e = if e == i64::MAX {
                        val(&d) - 1u32
                    } else {
                        Integer::from(e)
                    };
                    let mut a = limbs_of(&(q.to_rug() * val(&d) + e));
                    a.resize(a.len().max(d.len() + 1), 0);
                    a
                }
            };
            (a, d)
        })
}

proptest! {
    #[test]
    fn mu_div_qr_matches((a, d) in mu_operands()) {
        check_mu(&a, &d);
    }
}

/// The reciprocal's decrease step: with `A`'s low limbs all ones, `A Xh`
/// reaches `B^(n + h)` for most tops.
#[test]
fn invert_decrement_is_reached() {
    if INV_NEWTON_THRESHOLD > crate::testkit::TESTED_SIZE_MAX {
        return;
    }
    let n = tested(INV_NEWTON_THRESHOLD).max(3);
    let mut hits = 0;
    for t in 0..20u64 {
        let mut a = vec![Limb::MAX; n];
        a[n - 1] = (1 << 63) | t.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 1;
        hits += count(&hits::INV_DECREMENT, || drop(check_invert(&a)));
    }
    assert!(hits >= 1);
}

/// The largest remainders `a = q d + d - 1` make Barrett's estimates fall
/// short, so the correction runs.
#[test]
fn mu_correction_is_reached() {
    if INV_NEWTON_THRESHOLD > crate::testkit::TESTED_SIZE_MAX {
        return;
    }
    let t = tested(INV_NEWTON_THRESHOLD);
    for dn in [2, 3, 8, t, 2 * t + 1] {
        let mut corrections = 0;
        for qn in [1, dn / 3 + 1, dn, 2 * dn + 3] {
            let d: Vec<Limb> = (0..dn as u64)
                .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 1)
                .chain([1 << 63])
                .skip(1)
                .collect();
            let q: Vec<Limb> = (0..qn as u64).map(|i| !i).collect();
            let a = limbs_of(&(val(&q) * val(&d) + val(&d) - 1u32));
            corrections += count(&hits::MU_CORRECTION, || check_mu(&a, &d));
        }
        // An estimate can come out exact, but not for every shape.
        assert!(corrections >= 1, "dn = {dn}");
    }
}

/// A divisor just above `B^dn / 2`, of a power-of-two length, leaves
/// remainders of about `2D = B^dn` before the correction: the residue modulo
/// `B^dn - 1` alone is then ambiguous, and the low limb decides.
#[test]
fn mu_wrapped_is_reached() {
    let t = tested(INV_NEWTON_THRESHOLD).next_power_of_two();
    for dn in [2, 4, 8, t, 2 * t] {
        let mut d = vec![0; dn];
        d[dn - 1] = 1 << 63;
        d[0] = 0x9e37_79b9;
        let q: Vec<Limb> = (0..dn as u64)
            .map(|i| i.wrapping_mul(0xbf58_476d_1ce4_e5b9))
            .collect();
        let a = limbs_of(&(val(&q) * val(&d) + val(&d) - 1u32));
        assert!(
            count(&hits::MU_WRAPPED, || check_mu(&a, &d)) >= 1,
            "dn = {dn}"
        );
    }
}

/// Divisions through [`div_rem`] on both sides of the configured Barrett
/// threshold, beyond the property tests' sizes: random and largest-remainder
/// dividends, quotients of one block, two halves and several blocks, and
/// dividends that the normalising shift carries past their length.
/// Skipped when the configuration leaves these sizes to divide and conquer.
#[test]
fn div_rem_around_mu_threshold() {
    let t = MU_DIV_THRESHOLD;
    if t > 20_000 {
        return;
    }
    let lens = |dn: usize| [t - 1, t, dn, 2 * dn + 3];
    for dn in [t, t + 1, (2 * t).next_power_of_two()] {
        for qn in lens(dn) {
            let d: Vec<Limb> = (0..dn as u64)
                .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ i << 7)
                .collect();
            let q: Vec<Limb> = (0..qn as u64)
                .map(|i| !i.wrapping_mul(0xbf58_476d_1ce4_e5b9))
                .collect();
            for a in [
                limbs_of(&(val(&q) * val(&d) + val(&d) - 1u32)),
                limbs_of(&(val(&q) * val(&d))),
            ] {
                let (q, r) = run_div(&a, &d);
                assert_eq!(
                    (q.to_rug(), r.to_rug()),
                    gmp_div(&a, &d),
                    "dn = {dn}, qn = {qn}"
                );
            }
            // A divisor shifted by 3 to normalise and a dividend with its top
            // bit set: the shift carries a limb out, whose quotient limb takes
            // a schoolbook step of its own.
            let mut d = d;
            d[dn - 1] = d[dn - 1] >> 3 | 1 << 60;
            let mut a: Vec<Limb> = (0..(dn + qn) as u64)
                .map(|i| i.wrapping_mul(0xd6e8_feb8_6659_fd93) ^ i >> 3)
                .collect();
            a[dn + qn - 1] |= 1 << 63;
            let (q, r) = run_div(&a, &d);
            assert_eq!(
                (q.to_rug(), r.to_rug()),
                gmp_div(&a, &d),
                "carried: dn = {dn}, qn = {qn}"
            );
        }
    }
}

/// The block reciprocal never exceeds `B^(dn + k) / D`, on divisors where the
/// truncated top limbs alone would give one above it: top limbs
/// `Dt = floor(B^(2k + 1) / X0)`, whose reciprocal is `X0 B` to within a few
/// units, and a tail of ones.
#[test]
fn block_inverse_is_below_reciprocal() {
    let t = tested(INV_NEWTON_THRESHOLD);
    for k in [1, 2, 3, 8, t, 2 * t + 1] {
        for s in [1, 3] {
            for seed in 1..=4u64 {
                let x0: Vec<Limb> = (0..k as u64)
                    .map(|i| (i + seed).wrapping_mul(0x9e37_79b9_7f4a_7c15))
                    .chain([1])
                    .collect();
                let dt = base_pow(2 * k) * (Integer::from(1) << 64) / val(&x0);
                let d = dt * base_pow(s) + base_pow(s) - 1u32;
                let d = limbs_of(&d);
                let dn = d.len();
                assert_eq!(dn, k + 1 + s);
                let mut inv = vec![GARBAGE; k];
                let mut scratch = vec![GARBAGE; 2 * k + 2 + invert_scratch_len(k + 1)];
                block_inverse(&mut inv, &d, &mut scratch);
                let x = val(&inv) + base_pow(k);
                let (d, top) = (val(&d), base_pow(dn + k));
                assert!(Integer::from(&x * &d) <= top, "k = {k}, s = {s}: X above");
                assert!((x + 2u32) * d > top, "k = {k}, s = {s}: X too small");
            }
        }
    }
}

/// Exact multiples of divisors of `B^m - 1` times a power-of-two length `m`
/// (all ones, and `2 (B^m - 1) / 3`, by multiples of 3): the window and the
/// product are both multiples of `B^m - 1`, so the residue modulo `B^m - 1`
/// is zero, which the block may hold as `B^m - 1`.
#[test]
fn mu_exact_multiple_of_wrap_modulus() {
    let t = tested(INV_NEWTON_THRESHOLD).next_power_of_two();
    for dn in [2, 4, 8, t] {
        let ones = vec![Limb::MAX; dn];
        let two_thirds = limbs_of(&(val(&ones) / 3u32 * 2u32));
        for d in [ones, two_thirds] {
            for qn in [1, 2, dn / 2 + 1, dn, 2 * dn + 1] {
                for q in [vec![Limb::MAX; qn], (1..=qn as u64).collect()] {
                    let mut a = limbs_of(&(val(&q) * 3u32 * val(&d)));
                    a.resize(a.len().max(dn + 1), 0);
                    check_mu(&a, &d);
                }
            }
        }
    }
}

/// `(quotient, divides)` from [`divexact`], with garbage-filled buffers.
fn run_divexact(a: &[Limb], d: &[Limb]) -> (Integer, bool) {
    let mut q = vec![GARBAGE; a.len() - d.len() + 1];
    let mut scratch = vec![GARBAGE; divexact_scratch_len(a.len(), d.len())];
    let exact = divexact(&mut q, a, d, &mut scratch);
    (val(&q), exact)
}

/// Divisor lengths around the divide-and-conquer and block thresholds (and
/// one recursion level deeper) and where the inverse takes Newton steps, with
/// every quotient shape.
fn bdiv_divisor_lens() -> BoxedStrategy<usize> {
    len_around(&[
        1,
        2,
        3,
        tested(BDIV_Q_DC_THRESHOLD),
        2 * tested(BDIV_Q_DC_THRESHOLD) + 1,
        tested(BDIV_QR_DC_THRESHOLD),
        2 * tested(BDIV_QR_DC_THRESHOLD) + 1,
        tested(DIVEXACT_DIV_THRESHOLD),
        tested(BDIV_MU_THRESHOLD),
        2 * tested(BINV_NEWTON_THRESHOLD),
        tested(BINV_NEWTON_THRESHOLD),
    ])
}

/// `(a, d)` with `d` odd and normalised: `a = q d + e` for `e` in `{0, 1,
/// -1, d - 1, d}` (a multiple again), padded with zero limbs or not, or raw.
fn bdiv_operands() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    bdiv_divisor_lens()
        .prop_flat_map(|dn| (Just(dn), quotient_len(dn), 0..7u8, any::<bool>()))
        .prop_flat_map(|(dn, qn, shape, pad)| {
            let d = strategy::natural_with_len(Just(dn));
            let q = strategy::natural_with_len(Just(qn.max(1)));
            let raw = strategy::raw_limbs(dn + qn);
            (d, q, raw, Just(shape), Just(pad))
        })
        .prop_map(|(d, q, raw, shape, pad)| {
            let mut d = d.limbs().to_vec();
            d[0] |= 1;
            let dv = val(&d);
            let a = match shape {
                0 => return (raw.clone(), d),
                1 => q.to_rug() * &dv + 1u32,
                2 => q.to_rug() * &dv - 1u32,
                3 => q.to_rug() * &dv + (dv.clone() - 1u32),
                4 => (q.to_rug() + 1u32) * &dv,
                _ => q.to_rug() * &dv,
            };
            let mut a = limbs_of(&a.max(Integer::new()));
            a.resize(a.len().max(d.len()) + usize::from(pad), 0);
            (a, d)
        })
}

proptest! {
    #[test]
    fn divexact_matches((a, d) in bdiv_operands()) {
        let (av, dv) = (val(&a), val(&d));
        let (q, exact) = run_divexact(&a, &d);
        prop_assert_eq!(exact, av.is_divisible(&dv));
        if exact {
            prop_assert_eq!(q, av / dv);
        }
    }

    #[test]
    fn binvert_matches(
        d in len_around(&[1, 2, tested(BINV_NEWTON_THRESHOLD), 2 * tested(BINV_NEWTON_THRESHOLD) + 1])
            .prop_flat_map(strategy::raw_limbs)
            .prop_filter("non-empty", |d| !d.is_empty()),
    ) {
        let mut d = d;
        d[0] |= 1;
        let n = d.len();
        let mut x = vec![GARBAGE; n];
        let mut scratch = vec![GARBAGE; binvert_scratch_len(n)];
        binvert(&mut x, &d, &mut scratch);
        let modulus = base_pow(n);
        prop_assert_eq!(val(&x) * val(&d) % modulus, 1);
    }
}

/// Quotient lengths for the all-ones tests: the 2-adic quotient has one limb
/// more, so these give it fewer limbs than, as many as and more than the
/// divisor, and two and three times as many (a first block of every kind).
fn extreme_quotient_lens(dn: usize) -> Vec<usize> {
    let mut v = vec![
        1,
        2,
        dn / 2 + 1,
        dn.max(2) - 1,
        dn,
        dn + 1,
        2 * dn - 1,
        2 * dn,
        3 * dn - 1,
        3 * dn + 1,
    ];
    v.sort_unstable();
    v.dedup();
    v
}

/// All-ones divisors and quotients (the largest carries), one bit off a
/// multiple, zero, the divisor itself and `d B^2`.
#[test]
fn divexact_extremes() {
    let (qr, t) = (tested(BDIV_QR_DC_THRESHOLD), tested(BDIV_MU_THRESHOLD));
    let mut dns = vec![1, 2, 3, qr, qr + 1, 2 * qr + 1, t, t + 1];
    dns.sort_unstable();
    dns.dedup();
    for dn in dns {
        let d = vec![Limb::MAX; dn];
        for qn in extreme_quotient_lens(dn) {
            let q = vec![Limb::MAX; qn];
            let a = limbs_of(&(val(&q) * val(&d)));
            assert_eq!(
                run_divexact(&a, &d),
                (val(&q), true),
                "dn = {dn}, qn = {qn}"
            );
            let mut off = a.clone();
            off[0] ^= 2;
            assert!(!run_divexact(&off, &d).1, "dn = {dn}, qn = {qn}: off by 2");
        }
        assert_eq!(run_divexact(&vec![0; dn + 2], &d), (Integer::new(), true));
        assert_eq!(run_divexact(&d, &d), (Integer::from(1), true));
        let mut shifted = vec![0; 2];
        shifted.extend_from_slice(&d);
        assert_eq!(run_divexact(&shifted, &d), (base_pow(2), true));
    }
}

/// Dividends whose Hensel remainder is zero with a carry out of the top, so
/// that only the carry tells them from multiples: `A = B^(qn + dn) - Q' D`
/// for an all-ones `D` and `Q' = B^qn - B^(qn - 2) - 1`. It has one limb
/// fewer than `qn + dn`, and modulo `D = B^dn - 1` it is `B^qn`, not zero.
/// Through schoolbook and divide and conquer, with the carry left by every
/// kind of last block.
#[test]
fn divexact_remainder_carry() {
    let qr = tested(BDIV_QR_DC_THRESHOLD);
    let mut dns = vec![2, 3, qr, qr + 1, 2 * qr + 1];
    dns.sort_unstable();
    dns.dedup();
    for dn in dns {
        let d = vec![Limb::MAX; dn];
        for qn in extreme_quotient_lens(dn).into_iter().filter(|&qn| qn >= 2) {
            let q_neg = base_pow(qn) - base_pow(qn - 2) - 1u32;
            let a = limbs_of(&(base_pow(qn + dn) - q_neg * val(&d)));
            assert_eq!(a.len(), qn + dn - 1, "dn = {dn}, qn = {qn}");
            assert!(!run_divexact(&a, &d).1, "dn = {dn}, qn = {qn}");
        }
    }
}

/// The 2-adic quotient from [`bdiv_q`] to `qn` limbs, with garbage-filled
/// buffers.
fn run_bdiv_q(a: &[Limb], d: &[Limb], qn: usize) -> Vec<Limb> {
    let mut q = vec![GARBAGE; qn];
    let mut scratch = vec![GARBAGE; bdiv_q_scratch_len(qn, d.len())];
    bdiv_q(&mut q, a, d, &mut scratch);
    q
}

proptest! {
    /// Any dividend: `q d = a mod B^qn` for quotients shorter than, as long
    /// as and longer than the divisor; multiples: the exact quotient.
    #[test]
    fn bdiv_q_matches((a, d) in bdiv_operands(), cut in 0..3usize) {
        let qn = a.len() - d.len() + 1;
        let qn = match cut { 0 => qn, 1 => qn.div_ceil(2), _ => (qn + d.len()).min(a.len()) };
        let q = run_bdiv_q(&a, &d, qn);
        let modulus = base_pow(qn);
        prop_assert_eq!(val(&q) * val(&d) % &modulus, val(&a) % &modulus);
        let (av, dv) = (val(&a), val(&d));
        if cut == 0 && av.is_divisible(&dv) {
            prop_assert_eq!(val(&q), av / dv);
        }
    }
}

/// All-ones divisors and quotients through every block shape of [`bdiv_q`],
/// by divide and conquer (odd and even halves, first blocks of every kind)
/// and by blocks.
#[test]
fn bdiv_q_extremes() {
    let (dc, qr, t) = (
        tested(BDIV_Q_DC_THRESHOLD),
        tested(BDIV_QR_DC_THRESHOLD),
        tested(BDIV_MU_THRESHOLD),
    );
    let mut dns = vec![
        1,
        2,
        dc,
        dc + 1,
        2 * dc + 1,
        qr + 1,
        2 * qr,
        t,
        t + 1,
        2 * t + 1,
    ];
    dns.sort_unstable();
    dns.dedup();
    for dn in dns {
        let d = vec![Limb::MAX; dn];
        for qn in extreme_quotient_lens(dn) {
            let q = vec![Limb::MAX; qn];
            let a = limbs_of(&(val(&q) * val(&d)));
            let got = run_bdiv_q(&a, &d, a.len() - dn + 1);
            assert_eq!(val(&got), val(&q), "dn = {dn}, qn = {qn}");
        }
    }
}

/// The 2-adic quotient where the block updates and the inverse's Newton
/// steps take products modulo `B^m - 1` (transform lengths from the wrap
/// threshold): random operands, and all-ones ones, whose products are zero
/// modulo `B^m - 1` (the residue that needs the borrow into the wrapped
/// limbs), against GMP. Skipped when the configuration leaves these sizes to
/// other algorithms.
#[test]
fn bdiv_q_wrapped_at_large_sizes() {
    use crate::kernel::mul::{MUL_WRAP_NTT_THRESHOLD, wrap_len};
    let big = (2 * MUL_WRAP_NTT_THRESHOLD + 1).max(BDIV_MU_THRESHOLD);
    if big > 12_000 {
        return;
    }
    let mut state = 0x0fed_cba9_8765_4321_u64;
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
    // An odd length (the inverse of about half of it wraps too), a power of
    // two and a radix-3 length.
    for dn in [big, big.next_power_of_two(), 3 * wrap_len(big.div_ceil(3))] {
        for qn in [dn, 3 * dn] {
            let mut d = random(dn);
            d[0] |= 1;
            d[dn - 1] |= 1 << 63;
            let q = random(qn);
            for (d, q) in [(d, q), (vec![Limb::MAX; dn], vec![Limb::MAX; qn])] {
                let a = limbs_of(&(val(&q) * val(&d)));
                let got = run_bdiv_q(&a, &d, a.len() - dn + 1);
                assert_eq!(val(&got), val(&q), "dn = {dn}, qn = {qn}");
            }
        }
    }
}

/// Newton's step and the 2-adic quotient's blocks where they share prepared
/// transforms ([`binvert_shared`], [`bdiv_prepared_inverse`] and
/// [`bdiv_prepared_divisor`]), against GMP:
/// around the shortest lengths where the model takes them (one limb more, a
/// power-of-two and a radix-3 length, longer and shorter quotients, which it
/// may take or not), with random and all-ones divisors. At a power-of-two
/// length the all-ones divisor's product by the inverse is a multiple of
/// `B^m - 1`, whose residue may come in either form. Skipped when the
/// configuration takes neither below 20000 limbs.
#[test]
fn bdiv_shared_transforms_at_large_sizes() {
    let mut state = 0x7766_5544_3322_1100_u64;
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
    if let Some(n) = (BINV_NEWTON_THRESHOLD..20_000).find(|&n| binvert_shared(n).is_some()) {
        let three = 3 * (n / 3).next_power_of_two();
        for n in [n, n + 1, n.next_power_of_two(), three] {
            let mut d = random(n);
            d[0] |= 1;
            for d in [d, vec![Limb::MAX; n]] {
                let mut x = vec![GARBAGE; n];
                let mut scratch = vec![GARBAGE; binvert_scratch_len(n)];
                binvert(&mut x, &d, &mut scratch);
                assert_eq!(val(&x) * val(&d) % base_pow(n), 1, "n = {n}");
            }
        }
    }
    // Top bits set: the dividend has qn + dn limbs, and the 2-adic quotient
    // qn + 1, the length the model is asked about.
    let take = |dn: usize| bdiv_prepared_inverse(dn + 1, dn).is_some();
    if let Some(dn) = (BDIV_MU_THRESHOLD..20_000).find(|&dn| take(dn)) {
        // Two halves, a short last block of a long quotient, and halves of a
        // quotient shorter than the divisor.
        let (mut inverse, mut divisor) = (0, 0);
        for (dn, qn) in [(dn, dn), (dn + 1, 3 * dn), (2 * dn, dn + 1)] {
            let mut d = random(dn);
            d[0] |= 1;
            d[dn - 1] |= 1 << 63;
            let mut q = random(qn);
            q[qn - 1] |= 1 << 63;
            for (d, q) in [(d, q), (vec![Limb::MAX; dn], vec![Limb::MAX; qn])] {
                let a = limbs_of(&(val(&q) * val(&d)));
                inverse += usize::from(bdiv_prepared_inverse(qn + 1, dn).is_some());
                divisor += usize::from(bdiv_prepared_divisor(qn + 1, dn).is_some());
                let got = run_bdiv_q(&a, &d, qn + 1);
                assert_eq!(val(&got), val(&q), "dn = {dn}, qn = {qn}");
            }
        }
        assert!(inverse >= 2 && divisor >= 2, "dn = {dn}");
    }
}

/// Barrett division where the blocks reuse the divisor's and the
/// reciprocal's transforms (divisors past the NTT threshold, two blocks or
/// more), and moduli `B^m - 1` of `3 2^k` limbs, against GMP; and the plan
/// prepares them where expected. Skipped when the configuration leaves
/// these sizes to other algorithms.
#[test]
fn mu_div_prepared_at_large_sizes() {
    use crate::kernel::mul::{MUL_NTT_THRESHOLD, MUL_WRAP_NTT_THRESHOLD, wrap_len};
    let dn_big = (2 * MUL_NTT_THRESHOLD + 2).max(MU_DIV_THRESHOLD);
    if dn_big > 20_000 || MUL_WRAP_NTT_THRESHOLD > 8192 {
        return;
    }
    let mut state = 0x0123_4567_89ab_cdef_u64;
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
    // Even and odd divisors, 3 2^k + 1 just above a radix-3 length and
    // 3 2^k - 1 below one, with quotients that split into equal blocks or
    // leave a short last one, whichever count the plan takes.
    let three = 3 * (dn_big / 3).next_power_of_two();
    for dn in [dn_big, dn_big + 2, three - 1, three + 1] {
        for qn in [dn, 2 * dn + 5, dn / 2 + 1] {
            let mut d = random(dn);
            d[dn - 1] |= 1 << 63;
            let a = random(dn + qn);
            let mut q = vec![GARBAGE; qn];
            let mut n = a.clone();
            let mut scratch = vec![GARBAGE; mu_div_qr_scratch_len(dn + qn, dn)];
            let qh = mu_div_qr(&mut q, &mut n, &d, &mut scratch);
            let (gq, gr) = val(&a).div_rem(val(&d));
            let ours_q = val(&q) + (Integer::from(qh) << (64 * u32::try_from(qn).unwrap()));
            assert_eq!(ours_q, gq, "q: dn = {dn}, qn = {qn}");
            assert_eq!(val(&n[..dn]), gr, "r: dn = {dn}, qn = {qn}");
        }
    }
    // Blocks of k = a power of two past the NTT threshold, so the
    // reciprocal's products take the NTT and both are prepared: three for
    // qn = 3k - 1 and dn = 5k/4, the last one limb short (it takes the
    // reciprocal's top limbs, so not the prepared transform).
    let k = MUL_NTT_THRESHOLD.next_power_of_two();
    let (qn, dn) = (3 * k - 1, k + k / 4);
    if wrap_len(dn) >= MUL_WRAP_NTT_THRESHOLD && 2 * k <= 20_000 {
        assert_eq!((block_len(qn, dn), prepares(qn, dn)), (k, (true, true)));
        let mut d = random(dn);
        d[dn - 1] |= 1 << 63;
        let a = random(dn + qn);
        let mut q = vec![GARBAGE; qn];
        let mut n = a.clone();
        let mut scratch = vec![GARBAGE; mu_div_qr_scratch_len(dn + qn, dn)];
        let qh = mu_div_qr(&mut q, &mut n, &d, &mut scratch);
        let (gq, gr) = val(&a).div_rem(val(&d));
        let ours_q = val(&q) + (Integer::from(qh) << (64 * u32::try_from(qn).unwrap()));
        assert_eq!(ours_q, gq, "q: short last block");
        assert_eq!(val(&n[..dn]), gr, "r: short last block");
    }
    // Two blocks of a power-of-two length past the NTT threshold prepare
    // both; one block neither.
    let dn = 2 * MUL_NTT_THRESHOLD.next_power_of_two();
    if wrap_len(dn) >= MUL_WRAP_NTT_THRESHOLD {
        assert_eq!(
            (block_len(dn, dn), prepares(dn, dn)),
            (dn / 2, (true, true))
        );
        assert!(prepares(2 * dn + 5, dn).0);
    }
    assert_eq!(
        (block_len(dn / 4, dn), prepares(dn / 4, dn)),
        (dn / 4, (false, false))
    );
}

/// `divappr_q` on `(a, d)`, with garbage in its output and scratch: checks
/// that the quotient is the true one or one more.
fn check_divappr(a: &[Limb], d: &[Limb]) {
    let (an, dn) = (a.len(), d.len());
    let mut q = vec![GARBAGE; an - dn + 1];
    let mut scratch = vec![GARBAGE; divappr_q_scratch_len(an, dn)];
    divappr_q(&mut q, a, d, &mut scratch);
    let exact = gmp_div(a, d).0;
    let got = val(&q);
    assert!(
        got == exact || got == Integer::from(&exact + 1u32),
        "an = {an}, dn = {dn}: {} above",
        got - exact
    );
}

proptest! {
    #[test]
    fn divappr_q_within_one((a, d) in div_operands()) {
        prop_assume!(d.last() != Some(&0));
        check_divappr(&a, &d);
    }

    /// Barrett's estimate without the last remainder: at most 5 below.
    #[test]
    fn mu_divappr_q_bound((a, d) in mu_operands()) {
        let dn = d.len();
        let mut n = a.clone();
        let mut q = vec![GARBAGE; a.len() - dn];
        let mut scratch = vec![GARBAGE; mu_divappr_q_scratch_len(a.len(), dn)];
        let qh = mu_divappr_q(&mut q, &mut n, &d, &mut scratch);
        q.push(qh);
        let below = (gmp_div(&a, &d).0 - val(&q)).to_i64().unwrap_or(i64::MIN);
        prop_assert!((0..=5).contains(&below), "{} below", below);
    }

    /// Divide and conquer without the last block's correction: never below,
    /// and at most 2 per level of the recursion above.
    #[test]
    fn bz_divappr_q_bound((a, d) in div_operands()) {
        let dn = d.len();
        prop_assume!(dn >= 2 && a.len() > dn);
        let shift = d[dn - 1].leading_zeros();
        let mut ds = vec![0; dn];
        shl(&mut ds, &d, shift);
        let mut n = vec![0; a.len() + 1];
        n[a.len()] = shl(&mut n[..a.len()], &a, shift);
        let mut q = vec![GARBAGE; n.len() - dn];
        let mut scratch = vec![GARBAGE; bz_div_qr_scratch_len(n.len(), dn)];
        let v = reciprocal_3by2(ds[dn - 1], ds[dn - 2]);
        let qh = bz_divappr_q(&mut q, &mut n.clone(), &ds, v, &mut scratch);
        q.push(qh);
        let above = val(&q) - gmp_div(&a, &d).0;
        let bound = 2 * (usize::BITS - dn.leading_zeros()) + 2;
        prop_assert!(above >= 0 && above <= bound, "{} above, bound {}", above, bound);
    }
}

/// The largest quotient a length allows: `B^(qn + 1) - 1`, from an all-ones
/// dividend over `B^(dn - 1)`. The estimate may reach `B^(qn + 1)`, which
/// the quotient's limbs cannot hold; it is clamped to the exact quotient.
#[test]
fn divappr_q_largest_quotient() {
    let t = tested(DC_DIV_THRESHOLD).max(2);
    for (an, dn) in [
        (4, 2),
        (9, 3),
        (3 * t, t + 1),
        (4 * t + 5, 2 * t),
        (2 * t, t - 1),
    ] {
        let a = vec![Limb::MAX; an];
        let mut d = vec![0; dn];
        d[dn - 1] = 1;
        check_divappr(&a, &d);
        d[0] = 1;
        check_divappr(&a, &d);
    }
}

/// The guarded, approximate path at Barrett's sizes, with all-ones and
/// random operands and exact multiples.
#[test]
fn divappr_q_at_barrett_sizes() {
    let t = tested(MU_DIV_THRESHOLD);
    let mut next = rng(114);
    for (qn, dn) in [(t, t), (t + 1, t + 2), (2 * t, t), (t - 2, 3 * t)] {
        let d: Vec<Limb> = (0..dn).map(|_| next()).collect();
        let a: Vec<Limb> = (0..qn + dn).map(|_| next()).collect();
        check_divappr(&a, &d);
        check_divappr(&vec![Limb::MAX; qn + dn], &vec![Limb::MAX; dn]);
        let q: Vec<Limb> = (0..qn).map(|_| next()).collect();
        let mut exact = limbs_of(&(val(&q) * val(&d)));
        exact.resize(qn + dn, 0);
        check_divappr(&exact, &d);
    }
}

/// Plans with the divisor's transform prepared at a modulus shorter than
/// their blocks, which are folded first: the first such divisor and
/// quotient lengths past a transform length, for exact and approximate
/// quotients.
#[test]
fn mu_folds_blocks() {
    use crate::kernel::mul::MUL_WRAP_NTT_THRESHOLD;
    let low = MUL_WRAP_NTT_THRESHOLD.max(MU_DIV_THRESHOLD).max(4);
    // Divisors just past the transform lengths, quotients a few times longer.
    let lens: Vec<usize> = (2..12)
        .flat_map(|j| [1 << j, 3 << (j - 1)])
        .filter(|&len| len >= low && 2 * len <= crate::testkit::TESTED_SIZE_MAX)
        .flat_map(|len| [len + 1, len + 2, len + len / 16, len + len / 8])
        .collect();
    let mut next = rng(1142);
    for exact in [true, false] {
        // Plans as the cost model draws them; a configuration may have none.
        let Some((qn, dn)) = lens
            .iter()
            .flat_map(|&dn| [2, 3, 4, 6].map(|b| (b * dn, dn)))
            .find(|&(qn, dn)| folds_blocks(qn, dn, exact))
        else {
            continue;
        };
        let mut d: Vec<Limb> = (0..dn).map(|_| next()).collect();
        d[dn - 1] |= 1 << 63;
        let a: Vec<Limb> = (0..qn + dn).map(|_| next()).collect();
        if exact {
            check_mu(&a, &d);
        } else {
            let mut n = a.clone();
            let mut q = vec![GARBAGE; qn];
            let mut scratch = vec![GARBAGE; mu_divappr_q_scratch_len(a.len(), dn)];
            let qh = mu_divappr_q(&mut q, &mut n, &d, &mut scratch);
            q.push(qh);
            let below = (gmp_div(&a, &d).0 - val(&q)).to_i64().unwrap_or(i64::MIN);
            assert!((0..=5).contains(&below), "{below} below");
        }
    }
}

/// Divisors just past a transform length take their remainder modulo the
/// length below, with the low limbs from a short product: one limb past it
/// (two low limbs) and an eighth past it, with all-ones operands, random
/// ones and exact multiples, through the exact and approximate quotients.
#[test]
fn mu_wrap_below_divisor() {
    use crate::kernel::mul::MUL_WRAP_NTT_THRESHOLD;
    let fit = MUL_WRAP_NTT_THRESHOLD
        .max(MU_DIV_THRESHOLD)
        .max(16)
        .next_power_of_two();
    if fit > crate::testkit::TESTED_SIZE_MAX {
        return;
    }
    let mut next = rng(1141);
    for dn in [fit + 1, fit + fit / 8] {
        assert!(rem_wrap_len(dn) < dn, "dn = {dn}");

        for qn in [dn, dn / 2, 3 * dn] {
            let mut d: Vec<Limb> = (0..dn).map(|_| next()).collect();
            d[dn - 1] |= 1 << 63;
            let a: Vec<Limb> = (0..qn + dn).map(|_| next()).collect();
            check_mu(&a, &d);
            check_divappr(&a, &d);
            check_mu(&vec![Limb::MAX; qn + dn], &vec![Limb::MAX; dn]);
            let q: Vec<Limb> = (0..qn).map(|_| next()).collect();
            let mut exact = limbs_of(&(val(&q) * val(&d)));
            exact.resize(qn + dn, 0);
            check_mu(&exact, &d);
            check_divappr(&exact, &d);
        }
    }
}

/// An all-ones quotient with the largest remainder, `d - 1`, has an all-ones
/// guard quotient, so the approximate last block, which may be above, carries
/// into the limbs above it: through `divappr_q` at divide and conquer's
/// sizes, and through `bz_divappr_q` with blocks of all of `d`.
#[test]
fn divappr_carry_is_reached() {
    if DC_DIV_THRESHOLD > crate::testkit::TESTED_SIZE_MAX {
        return;
    }
    let t = DC_DIV_THRESHOLD.max(2);
    let mut next = rng(1143);
    let mut carries = 0;
    for dn in [t, t + 1, 2 * t + 1, 4 * t] {
        for qn in [dn - 2, dn, 3 * dn] {
            for _ in 0..8 {
                let mut d: Vec<Limb> = (0..dn).map(|_| next()).collect();
                d[dn - 1] |= 1 << 63;
                let q = vec![Limb::MAX; qn];
                let mut a = limbs_of(&(val(&q) * val(&d) + val(&d) - 1u32));
                a.resize(qn + dn, 0);
                carries += count(&hits::APPR_CARRY, || {
                    check_divappr(&a, &d);
                    let mut n = a.clone();
                    let mut q = vec![GARBAGE; qn];
                    let mut scratch = vec![GARBAGE; bz_div_qr_scratch_len(n.len(), dn)];
                    let v = reciprocal_3by2(d[dn - 1], d[dn - 2]);
                    let qh = bz_divappr_q(&mut q, &mut n, &d, v, &mut scratch);
                    q.push(qh);
                    let above = (val(&q) - gmp_div(&a, &d).0).to_i64().unwrap_or(-1);
                    assert!((0..=64).contains(&above), "{above} above");
                });
            }
        }
    }
    assert!(carries > 0);
}
