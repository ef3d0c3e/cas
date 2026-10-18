//! The constants against MPFR brackets at `prec + 64` bits, the tail bounds
//! against the true tails, and the cache's rounding and growth rules.
#![allow(
    clippy::many_single_char_names,
    reason = "the series' t, d, q and the precisions w, n of the formulas"
)]

use std::sync::RwLock;
use std::time::Instant;

use rug::float::{Constant, Round as MpRound};
use rug::{Float as MpFloat, Rational as MpRat};

use super::bsplit::{self, Series};
use super::{
    Atanh, Ball, Cache, Chudnovsky, E, atanh_tail, cached, chudnovsky_tail, e_tail, e_uncached,
    empty, log2_uncached, pi_uncached,
};
use crate::Rational;
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses};

type Oracle = fn(u32, MpRound) -> MpFloat;
type Ours = fn(u64) -> Ball;

fn mp_pi(p: u32, r: MpRound) -> MpFloat {
    MpFloat::with_val_round(p, Constant::Pi, r).0
}

fn mp_log2(p: u32, r: MpRound) -> MpFloat {
    MpFloat::with_val_round(p, Constant::Log2, r).0
}

fn mp_e(p: u32, r: MpRound) -> MpFloat {
    let mut x = MpFloat::with_val(p, 1);
    x.exp_round(r);
    x
}

const CONSTANTS: [(&str, Ours, Oracle); 3] = [
    ("pi", pi_uncached, mp_pi),
    ("log2", log2_uncached, mp_log2),
    ("e", e_uncached, mp_e),
];

/// Fails unless `b` meets MPFR's bracket at `prec + 64` bits and knows at
/// least `min_bits` bits.
fn check(name: &str, b: &Ball, prec: u64, oracle: Oracle, min_bits: i64) {
    let mid = b.mid().to_rational().to_rug();
    let rad = mag_value(b.rad());
    let br = bracket(u32::try_from(prec).unwrap(), oracle);
    if let Err(e) = check_encloses(&mid, &rad, &br) {
        panic!("{name} at {prec} bits: {e}");
    }
    let bits = b.rel_accuracy_bits();
    assert!(
        bits >= min_bits,
        "{name} at {prec} bits: only {bits} bits known"
    );
    assert!(
        b.mid().bit_length() <= prec + super::GUARD,
        "{name}: midpoint too long"
    );
}

#[test]
fn every_small_precision() {
    for (name, f, oracle) in CONSTANTS {
        for prec in 1..=300 {
            let b = f(prec);
            #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
            check(name, &b, prec, oracle, (prec + super::GUARD) as i64 - 6);
        }
    }
}

/// Precisions around powers of two and ten, up to `max`.
fn sweep(max: u64) -> Vec<u64> {
    let mut v: Vec<u64> = (6..24).map(|k| 1 << k).collect();
    v.extend([1000, 3333, 10_000, 31_623, 100_000, 316_228, 1_000_000]);
    v.extend(v.clone().into_iter().map(|p| p + 1));
    v.retain(|&p| p <= max);
    v.sort_unstable();
    v
}

#[test]
fn large_precisions() {
    // 10^6 bits takes seconds in release, minutes in debug.
    let max = if cfg!(debug_assertions) {
        20_000
    } else {
        1_000_000
    };
    for (name, f, oracle) in CONSTANTS {
        for prec in sweep(max) {
            #[allow(clippy::cast_possible_wrap, reason = "precisions below 2^63")]
            check(
                name,
                &f(prec),
                prec,
                oracle,
                (prec + super::GUARD) as i64 - 6,
            );
        }
    }
}

/// The partial sum of `n` terms.
fn partial<S: Series>(s: &S, n: u64) -> MpRat {
    let (t, d, _) = bsplit::sum(s, n);
    Rational::new(t, d).to_rug()
}

fn mag_value(m: &super::Mag) -> MpRat {
    m.to_float().expect("finite radius").to_rational().to_rug()
}

/// `|truth - partial(n)| <= tail(n)`, with `truth` from MPFR at 8000 bits
/// (its error, below 2^-7990, is far inside the slack of every bound here).
fn check_tail(name: &str, truth: &MpFloat, partial: &MpRat, tail: &MpRat) {
    let truth = truth.to_rational().expect("finite");
    let slack = MpRat::from((1, rug::Integer::from(1) << 7990u32));
    let err = MpRat::from(&truth - partial).abs();
    assert!(
        MpRat::from(&err + &slack) <= *tail,
        "{name}: tail {err} above the bound {tail}"
    );
}

#[test]
fn tail_bounds() {
    let p = 8000;
    // S = 426880 sqrt(10005) / pi.
    let s_true =
        MpFloat::with_val(p, 10_005).sqrt() * 426_880u32 / MpFloat::with_val(p, Constant::Pi);
    for n in 1..=40 {
        let name = format!("Chudnovsky, {n} terms");
        check_tail(
            &name,
            &s_true,
            &partial(&Chudnovsky, n),
            &mag_value(&chudnovsky_tail(n)),
        );
    }
    for m in [2, 3, 26, 4801, 8749] {
        let truth = (MpFloat::with_val(p, 1) / u32::try_from(m).unwrap()).atanh();
        for n in 1..=40 {
            let name = format!("atanh(1/{m}), {n} terms");
            let tail = mag_value(&atanh_tail(m, n));
            check_tail(&name, &truth, &partial(&Atanh(m), n), &tail);
        }
    }
    let truth = MpFloat::with_val(p, 1).exp();
    for n in 1..=200 {
        let (t, d, q) = bsplit::sum(&E, n);
        let name = format!("e, {n} terms");
        let sum = Rational::new(t, d).to_rug();
        check_tail(&name, &truth, &sum, &mag_value(&e_tail(n, &q)));
    }
}

#[test]
fn cache_rounds_and_grows() {
    // A private cache, so that other tests' requests cannot interleave.
    let cache: Cache = RwLock::new(empty());
    let prec_of = |c: &Cache| c.read().unwrap().prec;
    let b = cached(&cache, 1000, pi_uncached);
    check("cached pi", &b, 1000, mp_pi, 999);
    assert_eq!(prec_of(&cache), 1000);
    // Every lower request rounds the cached ball and still contains π.
    for prec in 1..=1000 {
        let b = cached(&cache, prec, pi_uncached);
        #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
        check("cached pi", &b, prec, mp_pi, prec as i64 - 1);
        assert!(b.mid().bit_length() <= prec, "not rounded to {prec} bits");
    }
    assert_eq!(prec_of(&cache), 1000, "a lower request recomputed");
    // Above it: at least twice the cached precision.
    let _ = cached(&cache, 1001, pi_uncached);
    assert_eq!(prec_of(&cache), 2000);
    let _ = cached(&cache, 5000, pi_uncached);
    assert_eq!(prec_of(&cache), 5000);
    check(
        "cached pi",
        &cached(&cache, 4999, pi_uncached),
        4999,
        mp_pi,
        4998,
    );
}

#[test]
fn public_entries() {
    // Threads racing on the global caches at mixed precisions.
    let handles: Vec<_> = (0..8u64)
        .map(|i| {
            std::thread::spawn(move || {
                for prec in [64, 300 + 50 * i, 2000 - 100 * i, 128] {
                    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
                    let min = prec as i64 - 1;
                    check("pi", &Ball::const_pi(prec), prec, mp_pi, min);
                    check("log2", &Ball::const_log2(prec), prec, mp_log2, min);
                    check("e", &Ball::const_e(prec), prec, mp_e, min);
                }
            })
        })
        .collect();
    for h in handles {
        h.join().unwrap();
    }
    super::clear_constant_cache();
    check("pi", &Ball::const_pi(100), 100, mp_pi, 99);
}

#[test]
#[should_panic(expected = "precision 0")]
fn precision_zero() {
    let _ = Ball::const_pi(0);
}

/// Where the time goes in π: the splitting tree level by level (each level
/// is the total minus the time of its subtrees) and the final steps. Run
/// with `--ignored --nocapture` in release.
#[test]
#[ignore = "timing, not a test"]
fn profile() {
    let time = |f: &mut dyn FnMut()| {
        let t = Instant::now();
        f();
        t.elapsed().as_secs_f64()
    };
    for prec in [10_000u64, 100_000, 1_000_000, 10_000_000] {
        let w = prec + super::GUARD;
        let n = w / super::CHUD_BITS_PER_TERM + 2;
        let mut sums = None;
        let total = time(&mut || sums = Some(bsplit::sum(&Chudnovsky, n)));
        let (t, d, _) = sums.unwrap();
        println!("pi at {prec} bits: {n} terms, splitting {total:.4} s");
        // Subtrees of 2^j terms: time(levels above 2^j) = total - Σ subtrees.
        let mut j = 0;
        while (1u64 << j) < n {
            j += 1;
            let size = 1u64 << j;
            let chunks = time(&mut || {
                let mut lo = 0;
                while lo < n {
                    let hi = (lo + size).min(n);
                    std::hint::black_box(bsplit::split(&Chudnovsky, lo, hi, true));
                    lo = hi;
                }
            });
            println!(
                "  subtrees of {size:>8} terms: {chunks:.4} s, levels above: {:.4} s",
                total - chunks
            );
        }
        let q = time(&mut || {
            std::hint::black_box(bsplit::quotient(&t, &d, w));
        });
        let rest = time(&mut || {
            std::hint::black_box(pi_uncached(prec));
        });
        println!("  T/Q division {q:.4} s; whole pi {rest:.4} s");
    }
}
