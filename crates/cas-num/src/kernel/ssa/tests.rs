//! Schönhage–Strassen checked against GMP: the Fermat-ring helpers on edge
//! residues, the modular product at several splits (with and without
//! recursion), and full products at every shape.

use proptest::prelude::*;
use rug::Integer;
use rug::ops::RemRounding;

use super::*;
use crate::testkit::strategy;

/// Filler for output and scratch buffers, so results can't rely on zeroed memory.
const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn val(a: &[Limb]) -> Integer {
    Integer::from_digits(a, rug::integer::Order::Lsf)
}

/// `2^(64 t) + 1`.
fn fermat(t: usize) -> Integer {
    (Integer::from(1) << (64 * u32::try_from(t).unwrap())) + 1
}

/// A normalised residue of `t + 1` limbs from `x mod F`.
fn residue(x: &Integer, t: usize) -> Vec<Limb> {
    let r = x % fermat(t);
    let mut v = r.to_digits::<Limb>(rug::integer::Order::Lsf);
    v.resize(t + 1, 0);
    v
}

fn is_normalised(x: &[Limb]) -> bool {
    let t = x.len() - 1;
    x[t] == 0 || (x[t] == 1 && x[..t].iter().all(|&l| l == 0))
}

/// Residues of `t + 1` limbs biased to the edges: 0, 1, 2^n' - 1, 2^n' = -1,
/// and values with all-ones or all-zero limbs.
fn edge_residue(t: usize) -> impl Strategy<Value = Integer> {
    let f = fermat(t);
    prop_oneof![
        Just(Integer::ZERO),
        Just(Integer::from(1)),
        Just(Integer::from(&f - 2)),
        Just(Integer::from(&f - 1)),
        strategy::raw_limbs(t).prop_map(|v| val(&v)),
        strategy::raw_limbs(t).prop_map(|v| val(&v) % Integer::from(1u32 << 20)),
    ]
}

fn residue_pair() -> impl Strategy<Value = (usize, Integer, Integer)> {
    (1..=6usize).prop_flat_map(|t| (Just(t), edge_residue(t), edge_residue(t)))
}

proptest! {
    #[test]
    fn ring_helpers((t, x, y) in residue_pair(), s in any::<usize>()) {
        let f = fermat(t);
        let (rx, ry) = (residue(&x, t), residue(&y, t));
        let check = |got: &[Limb], want: Integer| -> Result<(), TestCaseError> {
            prop_assert!(is_normalised(got));
            prop_assert_eq!(val(got), want.rem_euc(&f));
            Ok(())
        };
        let mut z = rx.clone();
        add_mod(&mut z, &ry);
        check(&z, Integer::from(&x + &y))?;
        let mut z = rx.clone();
        sub_mod(&mut z, &ry);
        check(&z, Integer::from(&x - &y))?;
        let (mut z, mut d) = (rx.clone(), vec![GARBAGE; t + 1]);
        add_sub(&mut z, &ry, &mut d);
        check(&z, Integer::from(&x + &y))?;
        check(&d, Integer::from(&x - &y))?;
        let mut z = rx.clone();
        neg_mod(&mut z);
        check(&z, Integer::from(-&x))?;
        let s = s % (2 * 64 * t);
        let mut z = vec![GARBAGE; t + 1];
        let mut tmp = vec![GARBAGE; t + 1];
        shift_mod(&mut z, &rx, s, &mut tmp);
        check(&z, Integer::from(&x << u32::try_from(s).unwrap()))?;
        // Any small top limb.
        let mut z = rx.clone();
        z[t] += 2;
        norm(&mut z);
        check(&z, &x + (Integer::from(&f - 1) << 1))?;
    }

    /// Pointwise products, plain and recursive, including the residue
    /// `2^n' = -1` on either side.
    #[test]
    fn mul_mod_matches(
        t in prop_oneof![1..=8usize, Just(FERMAT_THRESHOLD), Just(2 * FERMAT_THRESHOLD)],
        pick in (0..4usize, 0..4usize),
        seed in any::<u64>(),
        square in any::<bool>(),
    ) {
        let f = fermat(t);
        let make = |pick: usize, salt: u64| -> Integer {
            match pick {
                0 => Integer::from(&f - 1),
                1 => Integer::from(&f - 2),
                _ => {
                    let v: Vec<Limb> = (0..t as u64)
                        .map(|i| (i ^ seed ^ salt).wrapping_mul(0x9e37_79b9_7f4a_7c15))
                        .collect();
                    val(&v)
                }
            }
        };
        let (x, y) = (make(pick.0, 1), make(pick.1, 2));
        let (rx, ry) = (residue(&x, t), residue(&y, t));
        let b = if square { None } else { Some(&ry[..]) };
        let mut r = vec![GARBAGE; t + 1];
        let mut scratch = vec![GARBAGE; pointwise_scratch_len(t, square)];
        mul_mod(&mut r, &rx, b, &mut scratch);
        prop_assert!(is_normalised(&r));
        let want = if square { Integer::from(&x * &x) } else { Integer::from(&x * &y) };
        prop_assert_eq!(val(&r), want % &f);
    }

    /// The final fold, with the high part positive (and above or below the
    /// low part) or negative.
    #[test]
    fn fold_matches(
        n in 1..=6usize,
        lo in strategy::raw_limbs(6),
        hi in any::<u64>(),
        hi_len in 1..=3usize,
        negative in any::<bool>(),
    ) {
        let hi_len = hi_len.min(n + 1);
        // |hi| < 2^(64 (hi_len - 1)) + 2^63 keeps the sign limb meaningful.
        let mut h = vec![0; hi_len];
        h[0] = hi;
        if hi_len == 1 {
            h[0] >>= 1;
        }
        let mut hi_val = val(&h);
        if negative {
            hi_val = -hi_val;
            for l in &mut h {
                *l = !*l;
            }
            crate::kernel::add_1_assign(&mut h, 1);
        }
        let mut acc = lo[..n].to_vec();
        acc.extend_from_slice(&h);
        let mut out = vec![GARBAGE; n + 1];
        fold(&mut out, &mut acc);
        prop_assert!(is_normalised(&out));
        prop_assert_eq!(val(&out), (val(&lo[..n]) - hi_val).rem_euc(&fermat(n)));
    }

    /// The modular product at small splits; the pieces cover all of `n`.
    #[test]
    fn fermat_mul_matches(
        (k, m) in (4..=6u32, 1..=6usize),
        seed in any::<u64>(),
        square in any::<bool>(),
        extreme in 0..3usize,
    ) {
        let n = m << k;
        let make = |salt: u64| -> Vec<Limb> {
            match extreme {
                0 => vec![Limb::MAX; n],
                _ => (0..n as u64)
                    .map(|i| (i ^ seed ^ salt).wrapping_mul(0x9e37_79b9_7f4a_7c15))
                    .collect(),
            }
        };
        let (a, b) = (make(1), make(2));
        let b = if square { None } else { Some(&b[..]) };
        let mut out = vec![GARBAGE; n + 1];
        let mut scratch = vec![GARBAGE; fermat_scratch_len(n, k, square)];
        fermat_mul(&mut out, &a, b, k, &mut scratch);
        prop_assert!(is_normalised(&out));
        let want = val(&a) * val(b.unwrap_or(&a)) % fermat(n);
        prop_assert_eq!(val(&out), want);
    }

    #[test]
    fn mul_matches(
        (a, b) in (1..=700usize, 1..=700usize)
            .prop_flat_map(|(n, k)| (strategy::raw_limbs(n), strategy::raw_limbs(k))),
    ) {
        let mut r = vec![GARBAGE; a.len() + b.len()];
        let mut scratch = vec![GARBAGE; mul_scratch_len(a.len(), b.len())];
        mul(&mut r, &a, &b, &mut scratch);
        prop_assert_eq!(val(&r), val(&a) * val(&b));
    }

    #[test]
    fn sqr_matches(a in (1..=700usize).prop_flat_map(strategy::raw_limbs)) {
        let mut r = vec![GARBAGE; 2 * a.len()];
        let mut scratch = vec![GARBAGE; sqr_scratch_len(a.len())];
        sqr(&mut r, &a, &mut scratch);
        prop_assert_eq!(val(&r), val(&a).square());
    }
}

/// Pointwise products large enough to recurse, on all-ones and pseudo-random
/// operands, as products and squares.
#[test]
fn recursive_pointwise() {
    let k = 4;
    let m = FERMAT_THRESHOLD / 2;
    let n = m << k;
    assert!(inner_k(residue_len(m, k)).is_some());
    let random: Vec<Limb> = (0..n as u64)
        .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0x5555)
        .collect();
    let ones = vec![Limb::MAX; n];
    for (a, b) in [(&ones, &ones), (&random, &ones), (&random, &random)] {
        let mut out = vec![GARBAGE; n + 1];
        let mut scratch = vec![GARBAGE; fermat_scratch_len(n, k, false)];
        fermat_mul(&mut out, a, Some(b), k, &mut scratch);
        assert_eq!(val(&out), val(a) * val(b) % fermat(n));
        let mut scratch = vec![GARBAGE; fermat_scratch_len(n, k, true)];
        fermat_mul(&mut out, a, None, k, &mut scratch);
        assert_eq!(val(&out), val(a).square() % fermat(n));
    }
}

/// All-ones operands give the largest coefficients; the shapes cover every
/// row of the split table.
#[test]
fn all_ones() {
    for (an, bn) in [
        (1, 1),
        (2, 1),
        (17, 16),
        (300, 299),
        (1500, 1100),
        (5000, 4000),
    ] {
        let (a, b) = (vec![Limb::MAX; an], vec![Limb::MAX; bn]);
        let mut r = vec![GARBAGE; an + bn];
        let mut scratch = vec![GARBAGE; mul_scratch_len(an, bn)];
        mul(&mut r, &a, &b, &mut scratch);
        assert_eq!(val(&r), val(&a) * val(&b), "({an}, {bn})");
    }
}

/// The ticket's scale: 10^5-limb products and squares against GMP.
#[test]
fn at_1e5_limbs() {
    let a: Vec<Limb> = (0..100_000u64)
        .map(|i| i.wrapping_mul(0x2545_f491_4f6c_dd1d) ^ (i << 7))
        .collect();
    let b = vec![Limb::MAX; 99_999];
    let mut r = vec![GARBAGE; a.len() + b.len()];
    mul(
        &mut r,
        &a,
        &b,
        &mut vec![GARBAGE; mul_scratch_len(a.len(), b.len())],
    );
    assert_eq!(val(&r), val(&a) * val(&b));
    let mut r = vec![GARBAGE; 2 * a.len()];
    sqr(&mut r, &a, &mut vec![GARBAGE; sqr_scratch_len(a.len())]);
    assert_eq!(val(&r), val(&a).square());
}

/// Split tuning and a comparison with the NTT and GMP. Not a correctness
/// test; run with
/// `cargo test --release -p cas-num tune_ssa -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
fn tune_ssa() {
    use std::hint::black_box;

    use crate::kernel::tests::mul::time;

    let operand = |n: usize, seed: u64| -> Vec<Limb> {
        (0..n as u64)
            .map(|i| (i ^ seed).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0x5555)
            .collect()
    };
    let fermat_time = |n: usize, k: u32| {
        let (a, b) = (operand(n, 1), operand(n, 2));
        let mut out = vec![0; n + 1];
        let mut scratch = vec![0; fermat_scratch_len(n, k, false)];
        time(|| {
            fermat_mul(
                &mut out,
                black_box(&a),
                Some(black_box(&b)),
                k,
                &mut scratch,
            );
        }) / 1e3
    };

    println!("modular products: t  plain  k=4..8 (us)");
    for t in [256usize, 384, 512, 768, 1024, 1536, 2048, 3072] {
        let (a, b) = (operand(t, 1), operand(t, 2));
        let mut p = vec![0; 2 * t];
        let mut s = vec![0; crate::kernel::mul_scratch_len(t, t)];
        let plain = time(|| crate::kernel::mul(&mut p, black_box(&a), black_box(&b), &mut s)) / 1e3;
        let ks: Vec<String> = (4..=8)
            .map(|k| format!("{:.0}", fermat_time(t.next_multiple_of(1 << k), k)))
            .collect();
        println!("{t} {plain:.0} {}", ks.join(" "));
    }

    println!("products: L  k=4..12 (us)");
    for l in [2048usize, 4096, 8192, 16384, 32768, 65536, 131_072, 262_144] {
        let ks: Vec<String> = (4..=12)
            .map(|k| {
                let n = l.div_ceil(1 << k) << k;
                format!("{:.0}", fermat_time(n, k))
            })
            .collect();
        println!("{l} {}", ks.join(" "));
    }

    println!("n  ssa  ntt  gmp | ssa_sqr  ntt_sqr  gmp_sqr   (us, no allocation for ours)");
    for n in [
        1024usize, 2048, 3072, 4096, 8192, 16384, 32768, 50_000, 65536, 100_000, 131_072,
    ] {
        let (a, b) = (operand(n, 3), operand(n, 4));
        let mut r = vec![0; 2 * n];
        let mut s = vec![0; mul_scratch_len(n, n).max(crate::kernel::ntt::mul_scratch_len(n, n))];
        let ssa = time(|| mul(&mut r, black_box(&a), black_box(&b), &mut s)) / 1e3;
        let ntt =
            time(|| crate::kernel::ntt::mul(&mut r, black_box(&a), black_box(&b), &mut s)) / 1e3;
        let (ga, gb) = (val(&a), val(&b));
        let gmp = time(|| drop(black_box(Integer::from(black_box(&ga) * black_box(&gb))))) / 1e3;
        let ssa_sq = time(|| sqr(&mut r, black_box(&a), &mut s)) / 1e3;
        let ntt_sq = time(|| crate::kernel::ntt::sqr(&mut r, black_box(&a), &mut s)) / 1e3;
        let gmp_sq = time(|| drop(black_box(Integer::from(black_box(&ga).square_ref())))) / 1e3;
        println!("{n} {ssa:.0} {ntt:.0} {gmp:.0} | {ssa_sq:.0} {ntt_sq:.0} {gmp_sq:.0}");
    }
}
