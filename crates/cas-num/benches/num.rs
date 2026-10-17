//! Benchmarks of cas-num against GMP (through `rug`), op by op.
//!
//! Every benchmark has an `ours` and a `gmp` variant on the same inputs, so
//! `cargo bench -p cas-num` prints the two times side by side. Sizes sweep
//! powers of two up to `CAS_BENCH_MAX_LIMBS` limbs (default 4096; the
//! ticket's full sweep is `CAS_BENCH_MAX_LIMBS=100000`), which places the
//! algorithm thresholds on real data. Pass a filter to run part of it, e.g.
//! `cargo bench -p cas-num -- mul/`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;
use std::time::Duration;

use cas_num::{Integer, Natural};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

/// Deterministic pseudo-random limbs (`SplitMix64`), so runs compare.
fn limbs(len: usize, seed: u64) -> Vec<u64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ len as u64;
    (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        })
        .collect()
}

/// A value of exactly `len` limbs (top limb non-zero), ours and GMP's.
fn operand(len: usize, seed: u64) -> (Natural, rug::Integer) {
    let mut v = limbs(len, seed);
    if let Some(top) = v.last_mut() {
        *top |= 1 << 63;
    }
    let ours = Natural::from_limbs(&v);
    let gmp = rug::Integer::from_digits(&v, rug::integer::Order::Lsf);
    (ours, gmp)
}

/// Powers of two from 1 to the configured maximum, plus the maximum.
fn sizes(cap: usize) -> Vec<usize> {
    let max = std::env::var("CAS_BENCH_MAX_LIMBS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4096usize)
        .min(cap);
    let mut v: Vec<usize> = (0..usize::BITS)
        .map(|k| 1 << k)
        .take_while(|&n| n <= max)
        .collect();
    if v.last() != Some(&max) {
        v.push(max);
    }
    v
}

/// Runs `ours` and `gmp` for every size, on inputs from `make`.
fn sweep<I, O, G>(
    c: &mut Criterion,
    name: &str,
    cap: usize,
    make: impl Fn(usize) -> I,
    ours: O,
    gmp: G,
) where
    O: Fn(&I),
    G: Fn(&I),
{
    let mut group = c.benchmark_group(name);
    for n in sizes(cap) {
        let input = make(n);
        group.bench_with_input(BenchmarkId::new("ours", n), &input, |b, i| {
            b.iter(|| ours(i));
        });
        group.bench_with_input(BenchmarkId::new("gmp", n), &input, |b, i| {
            b.iter(|| gmp(i));
        });
    }
    group.finish();
}

fn pair(n: usize) -> ((Natural, rug::Integer), (Natural, rug::Integer)) {
    (operand(n, 1), operand(n, 2))
}

fn arithmetic(c: &mut Criterion) {
    sweep(
        c,
        "add",
        usize::MAX,
        pair,
        |((a, _), (b, _))| {
            black_box(a + b);
        },
        |((_, a), (_, b))| {
            black_box(rug::Integer::from(a + b));
        },
    );
    sweep(
        c,
        "mul",
        usize::MAX,
        pair,
        |((a, _), (b, _))| {
            black_box(a * b);
        },
        |((_, a), (_, b))| {
            black_box(rug::Integer::from(a * b));
        },
    );
    sweep(
        c,
        "sqr",
        usize::MAX,
        |n| operand(n, 3),
        |(a, _)| {
            black_box(a.square());
        },
        |(_, a)| {
            black_box(rug::Integer::from(a.square_ref()));
        },
    );
    // A 2n-limb dividend by an n-limb divisor.
    sweep(
        c,
        "div_rem",
        usize::MAX,
        |n| (operand(2 * n, 4), operand(n, 5)),
        |((a, _), (d, _))| {
            black_box(a.div_rem(d));
        },
        |((_, a), (_, d))| {
            black_box(<(rug::Integer, rug::Integer)>::from(a.div_rem_ref(d)));
        },
    );
}

fn number_theory(c: &mut Criterion) {
    sweep(
        c,
        "gcd",
        usize::MAX,
        pair,
        |((a, _), (b, _))| {
            black_box(a.gcd(b));
        },
        |((_, a), (_, b))| {
            black_box(rug::Integer::from(a.gcd_ref(b)));
        },
    );
    // Odd modulus and full-size exponent; quadratic per step, so capped.
    sweep(
        c,
        "pow_mod",
        64,
        |n| {
            let (mut m, mut mg) = operand(n, 6);
            if !m.bit(0) {
                m += 1u64;
                mg += 1;
            }
            (operand(n, 7), operand(n, 8), (m, mg))
        },
        |((a, _), (e, _), (m, _))| {
            black_box(a.pow_mod(e, m));
        },
        |((_, a), (_, e), (_, m))| {
            black_box(rug::Integer::from(a.pow_mod_ref(e, m).unwrap()));
        },
    );
    // A prime of n limbs: the full Baillie-PSW test runs.
    sweep(
        c,
        "primality",
        64,
        |n| {
            let p = operand(n, 9).0.next_prime();
            let pg = rug::Integer::from_digits(p.limbs(), rug::integer::Order::Lsf);
            (p, pg)
        },
        |(p, _)| {
            black_box(p.primality());
        },
        |(_, p)| {
            black_box(p.is_probably_prime(0));
        },
    );
}

fn radix(c: &mut Criterion) {
    sweep(
        c,
        "to_string",
        usize::MAX,
        |n| operand(n, 10),
        |(a, _)| {
            black_box(a.to_string());
        },
        |(_, a)| {
            black_box(a.to_string());
        },
    );
    sweep(
        c,
        "from_str",
        usize::MAX,
        |n| operand(n, 11).1.to_string(),
        |s| {
            black_box(s.parse::<Natural>().unwrap());
        },
        |s| {
            black_box(s.parse::<rug::Integer>().unwrap());
        },
    );
}

/// Values below 2^63, the dominant case in a CAS: a batch of 1024 sums and
/// products, against GMP and against plain `i128`.
fn small_integers(c: &mut Criterion) {
    const N: usize = 1024;
    let raw = limbs(2 * N, 12);
    let vals: Vec<i64> = raw
        .iter()
        .map(|&x| i64::from_ne_bytes(x.to_ne_bytes()) >> 1)
        .collect();
    let ours: Vec<Integer> = vals.iter().map(|&x| Integer::from(x)).collect();
    let gmp: Vec<rug::Integer> = vals.iter().map(|&x| rug::Integer::from(x)).collect();
    let mut group = c.benchmark_group("small");
    for (op, f_ours, f_gmp, f_prim) in [
        (
            "add",
            (|a: &Integer, b: &Integer| a + b) as fn(&Integer, &Integer) -> Integer,
            (|a: &rug::Integer, b: &rug::Integer| rug::Integer::from(a + b))
                as fn(&rug::Integer, &rug::Integer) -> rug::Integer,
            (|a: i64, b: i64| i128::from(a) + i128::from(b)) as fn(i64, i64) -> i128,
        ),
        (
            "mul",
            |a, b| a * b,
            |a, b| rug::Integer::from(a * b),
            |a, b| i128::from(a) * i128::from(b),
        ),
    ] {
        group.bench_function(BenchmarkId::new(op, "ours"), |b| {
            b.iter(|| {
                for p in ours.chunks_exact(2) {
                    black_box(f_ours(&p[0], &p[1]));
                }
            });
        });
        group.bench_function(BenchmarkId::new(op, "gmp"), |b| {
            b.iter(|| {
                for p in gmp.chunks_exact(2) {
                    black_box(f_gmp(&p[0], &p[1]));
                }
            });
        });
        group.bench_function(BenchmarkId::new(op, "i128"), |b| {
            b.iter(|| {
                for p in vals.chunks_exact(2) {
                    black_box(f_prim(black_box(p[0]), black_box(p[1])));
                }
            });
        });
    }
    group.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(20)
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(1))
}

criterion_group! {
    name = benches;
    config = config();
    targets = arithmetic, number_theory, radix, small_integers
}
criterion_main!(benches);
