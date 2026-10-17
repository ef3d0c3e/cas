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
use crate::kernel::div::{
    DC_DIV_THRESHOLD, div_2by1, div_3by2, hits, reciprocal_2by1, reciprocal_3by2, sb_div_qr,
};
use crate::kernel::{div_rem, div_rem_1, div_rem_scratch_len, divexact_1, rem_1, submul_1};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, ToRug, differential};

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

fn check_div(a: &[Limb], d: &[Limb]) {
    let (q, r) = run_div(a, d);
    assert_eq!((q.to_rug(), r.to_rug()), gmp_div(a, d));
}

/// [`check_div`] through the schoolbook kernel directly, whatever the
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
                count(&hits::DC_CORRECTION, || check_div(&a, &d)) >= 1,
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
