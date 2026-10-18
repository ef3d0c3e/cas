//! Benchmarks of binary floats and ball arithmetic against MPFR (through
//! `rug::Float`), CAS-DOC-5 §11.
//!
//! Every group times the MPFR baseline as `op/mpfr/bits`; the `Float`, `Mag`
//! and `Ball` tickets add `op/ours/bits` next to it on the same inputs.
//! Precisions sweep 64, 128, 256 … bits up to `CAS_BENCH_MAX_BITS` (default
//! 4096); constants and functions run at decimal digit counts up to
//! `CAS_BENCH_MAX_DIGITS` (default 10^4; the tickets' sweeps go to 10^6).
//! Pass a filter to run part of it: `cargo bench -p cas-num --bench ball -- '^mul/'`.
//!
//! Targets (CAS-DOC-5 §11): `Float` add/mul within 1.3× of MPFR at 1–4 limbs,
//! ball arithmetic within 2× of MPFR's plain floats, constants and elementary
//! functions within 2× after the fast-path tickets.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use rug::Float;
use rug::float::{Constant, Round};
use rug::ops::AssignRound;

fn env_max(var: &str, default: u32) -> u32 {
    std::env::var(var)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// Powers of two from 64 bits to `CAS_BENCH_MAX_BITS`.
fn precisions() -> Vec<u32> {
    let max = env_max("CAS_BENCH_MAX_BITS", 4096);
    (6..u32::BITS)
        .map(|k| 1 << k)
        .take_while(|&p| p <= max)
        .collect()
}

/// Powers of ten from 100 digits to `CAS_BENCH_MAX_DIGITS`.
fn digit_counts() -> Vec<u32> {
    let max = env_max("CAS_BENCH_MAX_DIGITS", 10_000);
    (2..10)
        .map(|k| 10u32.pow(k))
        .take_while(|&d| d <= max)
        .collect()
}

/// Bits for `digits` decimal digits: `⌈digits·log₂ 10⌉ + 1`, with log₂ 10
/// rounded up to six decimals.
fn bits(digits: u32) -> u32 {
    let b = (u64::from(digits) * 3_321_929).div_ceil(1_000_000) + 1;
    u32::try_from(b).expect("digit count too large")
}

/// Deterministic operands with all `prec` bits significant: `sqrt(seed)`
/// rounded to `prec` bits, scaled to land near 1.
fn operand(prec: u32, seed: u32) -> Float {
    let x = Float::with_val(prec, Float::with_val(prec, seed).sqrt_ref());
    x / seed.isqrt()
}

/// One group per binary operation at every precision. MPFR rounds toward
/// zero, as ball midpoints do (§5.2).
fn arithmetic(c: &mut Criterion) {
    type Op = fn(&mut Float, &Float, &Float);
    let ops: [(&str, Op); 4] = [
        ("add", |r, a, b| _ = r.assign_round(a + b, Round::Zero)),
        ("mul", |r, a, b| _ = r.assign_round(a * b, Round::Zero)),
        ("div", |r, a, b| _ = r.assign_round(a / b, Round::Zero)),
        ("sqrt", |r, a, _| {
            _ = r.assign_round(a.sqrt_ref(), Round::Zero);
        }),
    ];
    for (name, op) in ops {
        let mut group = c.benchmark_group(name);
        for prec in precisions() {
            let (a, b) = (operand(prec, 2), operand(prec, 3));
            let mut r = Float::new(prec);
            group.bench_function(BenchmarkId::new("mpfr", prec), |bench| {
                bench.iter(|| op(&mut r, black_box(&a), black_box(&b)));
            });
        }
        group.finish();
    }
}

/// π, log 2 and e from scratch. MPFR caches constants, so each iteration
/// clears the cache first; the cost of clearing is negligible next to the
/// computation from 100 digits up.
fn constants(c: &mut Criterion) {
    let consts = [("pi", Constant::Pi), ("log2", Constant::Log2)];
    for (name, k) in consts {
        let mut group = c.benchmark_group(name);
        for d in digit_counts() {
            group.bench_function(BenchmarkId::new("mpfr", d), |bench| {
                bench.iter(|| {
                    rug::float::free_cache(rug::float::FreeCache::All);
                    black_box(Float::with_val(bits(d), k))
                });
            });
        }
        group.finish();
    }
    let mut group = c.benchmark_group("e");
    for d in digit_counts() {
        group.bench_function(BenchmarkId::new("mpfr", d), |bench| {
            bench.iter(|| black_box(Float::with_val(bits(d), 1).exp()));
        });
    }
    group.finish();
}

/// Elementary functions at an argument with all bits significant.
fn functions(c: &mut Criterion) {
    type Fun = fn(&Float) -> Float;
    let funs: [(&str, Fun); 4] = [
        ("exp", |x| x.clone().exp()),
        ("log", |x| x.clone().ln()),
        ("sin", |x| x.clone().sin()),
        ("atan", |x| x.clone().atan()),
    ];
    for (name, f) in funs {
        let mut group = c.benchmark_group(name);
        for d in digit_counts() {
            let x = operand(bits(d), 2);
            group.bench_function(BenchmarkId::new("mpfr", d), |bench| {
                bench.iter(|| f(black_box(&x)));
            });
        }
        group.finish();
    }
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
    targets = arithmetic, constants, functions
}
criterion_main!(benches);
