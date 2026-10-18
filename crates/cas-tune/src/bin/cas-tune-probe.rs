//! Times one cas-num operation at several sizes, as built with the tuning
//! configuration under test. Run by `cas-tune`, which rebuilds this binary
//! for every candidate configuration.
//!
//! `cas-tune-probe <op> <min-us> <size>...` prints `<size> <nanoseconds>` per
//! size: the best of five runs, each of at least `min-us` microseconds.
//! Sizes are limbs, or decimal digits for `from_str`.
#![allow(
    clippy::many_single_char_names,
    reason = "operands named as in the kernels (a, b, d, q, r, s)"
)]

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::{Duration, Instant};

use cas_num::modular::ModRing;
use cas_num::{Natural, kernel};

/// Deterministic pseudo-random limbs (`SplitMix64`), top limb non-zero.
fn limbs(len: usize, seed: u64) -> Vec<u64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ len as u64;
    let mut v: Vec<u64> = (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        })
        .collect();
    if let Some(top) = v.last_mut() {
        *top |= 1 << 63;
    }
    v
}

/// Nanoseconds per call of `f`: best of 5 runs of at least `min` each.
fn time(min: Duration, mut f: impl FnMut()) -> f64 {
    let mut iters = 1u32;
    loop {
        let start = Instant::now();
        for _ in 0..iters {
            f();
        }
        if start.elapsed() >= min {
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

fn measure(op: &str, n: usize, min: Duration) -> f64 {
    match op {
        "mul" => {
            let (a, b) = (limbs(n, 1), limbs(n, 2));
            let mut r = vec![0; 2 * n];
            let mut s = vec![0; kernel::mul_scratch_len(n, n)];
            time(min, || {
                kernel::mul(&mut r, black_box(&a), black_box(&b), &mut s);
            })
        }
        "sqr" => {
            let a = limbs(n, 3);
            let mut r = vec![0; 2 * n];
            let mut s = vec![0; kernel::sqr_scratch_len(n)];
            time(min, || kernel::sqr(&mut r, black_box(&a), &mut s))
        }
        "mullo" => {
            let (a, b) = (limbs(n, 1), limbs(n, 2));
            let mut r = vec![0; n];
            let mut s = vec![0; kernel::mullo_scratch_len(n)];
            time(min, || {
                kernel::mullo(&mut r, black_box(&a), black_box(&b), &mut s);
            })
        }
        "divexact" | "bdiv_q" => {
            // A 2n-limb multiple of an odd n-limb divisor.
            let mut d = limbs(n, 5);
            d[0] |= 1;
            let q = limbs(n, 4);
            let mut a = vec![0; 2 * n];
            let mut s = vec![0; kernel::mul_scratch_len(n, n)];
            kernel::mul(&mut a, &q, &d, &mut s);
            let mut quot = vec![0; n + 1];
            if op == "bdiv_q" {
                let mut s = vec![0; kernel::bdiv_q_scratch_len(n + 1, n)];
                time(min, || {
                    kernel::bdiv_q(&mut quot, black_box(&a), black_box(&d), &mut s);
                })
            } else {
                let mut s = vec![0; kernel::divexact_scratch_len(2 * n, n)];
                time(min, || {
                    assert!(kernel::divexact(
                        &mut quot,
                        black_box(&a),
                        black_box(&d),
                        &mut s
                    ));
                })
            }
        }
        "div" => {
            // A 2n-limb dividend by a normalised n-limb divisor.
            let (a, d) = (limbs(2 * n, 4), limbs(n, 5));
            let (mut q, mut r) = (vec![0; n + 1], vec![0; n]);
            let mut s = vec![0; kernel::div_rem_scratch_len(2 * n, n)];
            time(min, || {
                kernel::div_rem(&mut q, &mut r, black_box(&a), black_box(&d), &mut s);
            })
        }
        "redc" => {
            // One Montgomery product modulo an odd n-limb modulus.
            let mut m = limbs(n, 6);
            m[0] |= 1;
            let ring = ModRing::new(Natural::from_limbs(&m));
            let x = ring.elem(&Natural::from_limbs(&limbs(n, 7)));
            let y = ring.elem(&Natural::from_limbs(&limbs(n, 8)));
            time(min, || {
                drop(black_box(ring.mul(black_box(&x), black_box(&y))));
            })
        }
        "powm" => {
            // x^e modulo an odd n-limb modulus, with an n-limb exponent.
            let mut m = limbs(n, 6);
            m[0] |= 1;
            let ring = ModRing::new(Natural::from_limbs(&m));
            let x = ring.elem(&Natural::from_limbs(&limbs(n, 7)));
            let e = Natural::from_limbs(&limbs(n, 9));
            time(min, || {
                drop(black_box(ring.pow(black_box(&x), black_box(&e))));
            })
        }
        "gcd" => {
            let (a, b) = (
                Natural::from_limbs(&limbs(n, 11)),
                Natural::from_limbs(&limbs(n, 12)),
            );
            time(min, || drop(black_box(black_box(&a).gcd(black_box(&b)))))
        }
        "to_string" => {
            let x = Natural::from_limbs(&limbs(n, 9));
            time(min, || drop(black_box(black_box(&x).to_string())))
        }
        "from_str" => {
            let mut digits = String::from("9");
            for l in limbs(n.div_ceil(19), 10) {
                write!(digits, "{:019}", l % 10_000_000_000_000_000_000).expect("String");
            }
            digits.truncate(n.max(1));
            time(min, || {
                drop(black_box(black_box(&digits).parse::<Natural>()));
            })
        }
        _ => panic!("unknown op `{op}`"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [op, min_us, sizes @ ..] = args.as_slice() else {
        panic!("usage: cas-tune-probe <op> <min-us> <size>...");
    };
    let min = Duration::from_micros(min_us.parse().expect("min-us"));
    for size in sizes {
        let n: usize = size.parse().expect("size");
        println!("{n} {:.1}", measure(op, n, min));
    }
}
