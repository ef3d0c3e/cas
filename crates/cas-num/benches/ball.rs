//! Benchmarks of binary floats and ball arithmetic against MPFR (through
//! `rug::Float`), CAS-DOC-5 §11.
//!
//! Every group times the MPFR baseline as `op/mpfr/bits` next to ours,
//! `op/ours/bits` (correctly rounded, toward zero) and, for division and
//! square root, `op/approx/bits` (the midpoint entries), and `op/ball/bits`
//! (balls with a radius of about 2^-20 relative), on the same inputs.
//! Precisions sweep 64, 128, 256 … bits up to `CAS_BENCH_MAX_BITS` (default
//! 4096); constants and functions run at decimal digit counts up to
//! `CAS_BENCH_MAX_DIGITS` (default 10^4; the tickets' sweeps go to 10^6).
//! The `mulhi` group times the short products of the midpoint kernels
//! against the full product, on limb counts.
//! Pass a filter to run part of it: `cargo bench -p cas-num --bench ball -- '^mul/'`.
//!
//! Targets (CAS-DOC-5 §11): `Float` add/mul within 1.3× of MPFR at 1–4 limbs,
//! ball arithmetic within 2× of MPFR's plain floats, constants and elementary
//! functions within 2× after the fast-path tickets.
#![allow(missing_docs, reason = "benchmark binary")]
#![allow(clippy::many_single_char_names, reason = "operands a, b, x, y, r")]

use std::hint::black_box;
use std::time::Duration;

use cas_num::{ball, kernel};
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

/// The short products the midpoints use against the full product, at
/// 1 … `CAS_BENCH_MAX_LIMBS` (default 4096) limbs: `mulhi/{mulhi,mullo,mul}/n`.
fn short_products(c: &mut Criterion) {
    let max = env_max("CAS_BENCH_MAX_LIMBS", 4096) as usize;
    let mut group = c.benchmark_group("mulhi");
    for n in (0..usize::BITS)
        .map(|k| 1usize << k)
        .take_while(|&n| n <= max)
    {
        let (a, b) = (limbs(n, 1), limbs(n, 2));
        let mut r = vec![0; 2 * n];
        let mut scratch = vec![0; kernel::mulhi_scratch_len(n)];
        group.bench_function(BenchmarkId::new("mulhi", n), |bench| {
            bench.iter(|| kernel::mulhi(&mut r[..n], black_box(&a), black_box(&b), &mut scratch));
        });
        let mut scratch = vec![0; kernel::mullo_scratch_len(n)];
        group.bench_function(BenchmarkId::new("mullo", n), |bench| {
            bench.iter(|| kernel::mullo(&mut r[..n], black_box(&a), black_box(&b), &mut scratch));
        });
        let mut scratch = vec![0; kernel::mul_scratch_len(n, n)];
        group.bench_function(BenchmarkId::new("mul", n), |bench| {
            bench.iter(|| kernel::mul(&mut r, black_box(&a), black_box(&b), &mut scratch));
        });
    }
    group.finish();
}

/// Radius operations: `mag/{add,mul,div,sub_lower,exp,log1p}`, on values
/// with `i64` exponents (target: a few ns, no allocation).
fn mag_ops(c: &mut Criterion) {
    use ball::Mag;
    let mag = |x: f64| Mag::from_float(&ball::Float::from_f64(x).expect("finite"));
    let (a, b) = (mag(1.234_567_8e-20), mag(9.876_543_2e-22));
    let (small, mid) = (mag(3.0e-5), mag(2.75));
    let mut group = c.benchmark_group("mag");
    group.bench_function("add", |bench| {
        bench.iter(|| black_box(&a).add(black_box(&b)));
    });
    group.bench_function("mul", |bench| {
        bench.iter(|| black_box(&a).mul(black_box(&b)));
    });
    group.bench_function("div", |bench| {
        bench.iter(|| black_box(&a).div(black_box(&b)));
    });
    group.bench_function("sub_lower", |bench| {
        bench.iter(|| black_box(&a).sub_lower(black_box(&b)));
    });
    group.bench_function("exp/small", |bench| bench.iter(|| black_box(&small).exp()));
    group.bench_function("exp/mid", |bench| bench.iter(|| black_box(&mid).exp()));
    group.bench_function("log1p/mid", |bench| bench.iter(|| black_box(&mid).log1p()));
    group.finish();
}

/// One group per binary operation at every precision. MPFR rounds toward
/// zero, as ball midpoints do (§5.2).
fn arithmetic(c: &mut Criterion) {
    type Op = fn(&mut Float, &Float, &Float);
    type OurOp = fn(&ball::Float, &ball::Float, u64) -> ball::Float;
    // (group, MPFR, ours correctly rounded, ours approximate)
    let ops: [(&str, Op, OurOp, Option<OurOp>); 4] = [
        (
            "add",
            |r, a, b| _ = r.assign_round(a + b, Round::Zero),
            |a, b, p| a.add_round(b, p, ball::Round::Down).0,
            None,
        ),
        (
            "mul",
            |r, a, b| _ = r.assign_round(a * b, Round::Zero),
            |a, b, p| a.mul_round(b, p, ball::Round::Down).0,
            None,
        ),
        (
            "div",
            |r, a, b| _ = r.assign_round(a / b, Round::Zero),
            |a, b, p| a.div_round(b, p, ball::Round::Down).0,
            Some(|a, b, p| a.div_approx(b, p).expect("non-zero")),
        ),
        (
            "sqrt",
            |r, a, _| {
                _ = r.assign_round(a.sqrt_ref(), Round::Zero);
            },
            |a, _, p| a.sqrt_round(p, ball::Round::Down).0,
            Some(|a, _, p| a.sqrt_approx(p).expect("positive")),
        ),
    ];
    for (name, op, ours, approx) in ops {
        let mut group = c.benchmark_group(name);
        for prec in precisions() {
            let (a, b) = (operand(prec, 2), operand(prec, 3));
            let (x, y) = (to_ours(&a), to_ours(&b));
            group.bench_function(BenchmarkId::new("ours", prec), |bench| {
                bench.iter(|| ours(black_box(&x), black_box(&y), u64::from(prec)));
            });
            if let Some(ball_op) = balls_op(name) {
                // Balls with a radius of about 2^-20 relative, as after a few
                // operations.
                let ball = |x: &ball::Float| {
                    ball::Ball::new(
                        x.clone(),
                        ball::Mag::from_float(x).mul_2exp(&cas_num::Integer::from(-20)),
                    )
                };
                let (bx, by) = (ball(&x), ball(&y));
                group.bench_function(BenchmarkId::new("ball", prec), |bench| {
                    bench.iter(|| ball_op(black_box(&bx), black_box(&by), u64::from(prec)));
                });
            }
            if let Some(approx) = approx {
                group.bench_function(BenchmarkId::new("approx", prec), |bench| {
                    bench.iter(|| approx(black_box(&x), black_box(&y), u64::from(prec)));
                });
            }
            let mut r = Float::new(prec);
            group.bench_function(BenchmarkId::new("mpfr", prec), |bench| {
                bench.iter(|| op(&mut r, black_box(&a), black_box(&b)));
            });
        }
        group.finish();
    }
}

/// The ball operation of a group (target: within 2× of MPFR's plain float).
fn balls_op(name: &str) -> Option<fn(&ball::Ball, &ball::Ball, u64) -> ball::Ball> {
    Some(match name {
        "add" => |a, b, p| a.add(b, p),
        "mul" => |a, b, p| a.mul(b, p),
        "div" => |a, b, p| a.div(b, p),
        "sqrt" => |a, _, p| a.sqrt(p),
        _ => return None,
    })
}

/// Machine balls against 53-bit balls on the same values (radius about
/// 2^-20 relative), the fast path's target being 2–5 ns per operation.
fn machine(c: &mut Criterion) {
    type Machine = fn(&ball::MachineBall) -> Option<ball::MachineBall>;
    type Libm = fn(f64) -> f64;
    type Fallback = fn(&ball::Ball, u64) -> ball::Ball;
    let a = ball::MachineBall::new(std::f64::consts::PI, 3e-6).expect("regular");
    let b = ball::MachineBall::new(std::f64::consts::E, 2.6e-6).expect("regular");
    let (ba, bb) = (a.to_ball(), b.to_ball());
    let mut group = c.benchmark_group("machine");
    group.bench_function(BenchmarkId::new("add", "machine"), |bench| {
        bench.iter(|| black_box(&a).add(black_box(&b)));
    });
    group.bench_function(BenchmarkId::new("mul", "machine"), |bench| {
        bench.iter(|| black_box(&a).mul(black_box(&b)));
    });
    group.bench_function(BenchmarkId::new("div", "machine"), |bench| {
        bench.iter(|| black_box(&a).div(black_box(&b)));
    });
    group.bench_function(BenchmarkId::new("sqrt", "machine"), |bench| {
        bench.iter(|| black_box(&a).sqrt());
    });
    // Elementary functions: the double-double kernels, libm (no proven
    // bound, for scale) and the 64-bit Ball fallback they replace.
    let x = ball::MachineBall::new(0.7, 1e-12).expect("regular");
    let funs: [(&str, Machine, Libm, Fallback); 4] = [
        ("exp", ball::MachineBall::exp, f64::exp, ball::Ball::exp),
        ("log", ball::MachineBall::log, f64::ln, ball::Ball::log),
        ("sin", ball::MachineBall::sin, f64::sin, ball::Ball::sin),
        ("atan", ball::MachineBall::atan, f64::atan, ball::Ball::atan),
    ];
    for (name, ours, libm, fallback) in funs {
        group.bench_function(BenchmarkId::new(name, "machine"), |bench| {
            bench.iter(|| ours(black_box(&x)));
        });
        group.bench_function(BenchmarkId::new(name, "libm"), |bench| {
            bench.iter(|| libm(black_box(0.7)));
        });
        group.bench_function(BenchmarkId::new(name, "via_ball"), |bench| {
            bench.iter(|| black_box(&x).via_ball(fallback));
        });
    }
    for name in ["add", "mul", "div", "sqrt"] {
        let ball_op = balls_op(name).expect("ball op");
        group.bench_function(BenchmarkId::new(name, "ball"), |bench| {
            bench.iter(|| ball_op(black_box(&ba), black_box(&bb), 53));
        });
    }
    group.finish();
}

/// The same value as our float, through its integer mantissa and exponent.
fn to_ours(x: &Float) -> ball::Float {
    let (m, e) = x.to_integer_exp().expect("finite");
    let mag = cas_num::Natural::from_limbs(&m.as_abs().to_digits::<u64>(rug::integer::Order::Lsf));
    let m = cas_num::Integer::from_parts(m < 0, mag);
    ball::Float::from_integer_exp(&m, &cas_num::Integer::from(i64::from(e)))
}

/// π, log 2 and e from scratch. Both sides cache constants, so each
/// iteration clears the caches first; the cost of clearing is negligible
/// next to the computation from 100 digits up. MPFR has no `e` constant:
/// its row is `exp(1)`.
fn constants(c: &mut Criterion) {
    type Ours = fn(u64) -> ball::Ball;
    let consts: [(&str, Constant, Ours); 2] = [
        ("pi", Constant::Pi, ball::Ball::const_pi),
        ("log2", Constant::Log2, ball::Ball::const_log2),
    ];
    for (name, k, ours) in consts {
        let mut group = c.benchmark_group(name);
        for d in digit_counts() {
            group.bench_function(BenchmarkId::new("ours", d), |bench| {
                bench.iter(|| {
                    ball::clear_constant_cache();
                    black_box(ours(u64::from(bits(d))))
                });
            });
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
        group.bench_function(BenchmarkId::new("ours", d), |bench| {
            bench.iter(|| {
                ball::clear_constant_cache();
                black_box(ball::Ball::const_e(u64::from(bits(d))))
            });
        });
        group.bench_function(BenchmarkId::new("mpfr", d), |bench| {
            bench.iter(|| black_box(Float::with_val(bits(d), 1).exp()));
        });
    }
    group.finish();
}

/// π printed to `d` digits and a `d`-digit decimal parsed, against MPFR's
/// `to_string_radix` and `Float::parse` at the same precision.
fn decimal(c: &mut Criterion) {
    let mut group = c.benchmark_group("decimal");
    for d in digit_counts() {
        let p = bits(d);
        let ours = ball::Ball::const_pi(u64::from(p) + 16);
        let theirs = Float::with_val(p + 16, Constant::Pi);
        let n = usize::try_from(d).expect("small");
        group.bench_function(BenchmarkId::new("print/ours", d), |bench| {
            bench.iter(|| black_box(&ours).to_decimal(u64::from(d), ball::DecimalStyle::Digits));
        });
        group.bench_function(BenchmarkId::new("print/mpfr", d), |bench| {
            bench.iter(|| black_box(&theirs).to_string_radix(10, Some(n)));
        });
        let text = theirs.to_string_radix(10, Some(n));
        group.bench_function(BenchmarkId::new("parse/ours", d), |bench| {
            bench.iter(|| ball::Ball::from_decimal(black_box(&text), u64::from(p)));
        });
        group.bench_function(BenchmarkId::new("parse/mpfr", d), |bench| {
            bench.iter(|| Float::with_val(p, Float::parse(black_box(&text)).expect("decimal")));
        });
    }
    group.finish();
}

/// Elementary functions at an argument with all bits significant.
fn functions(c: &mut Criterion) {
    type Fun = fn(&Float) -> Float;
    type Ours = Option<fn(&ball::Ball, u64) -> ball::Ball>;
    let funs: [(&str, Fun, Ours); 4] = [
        ("exp", |x| x.clone().exp(), Some(ball::Ball::exp)),
        ("log", |x| x.clone().ln(), Some(ball::Ball::log)),
        ("sin", |x| x.clone().sin(), Some(ball::Ball::sin)),
        ("atan", |x| x.clone().atan(), Some(ball::Ball::atan)),
    ];
    for (name, f, ours) in funs {
        let mut group = c.benchmark_group(name);
        for d in digit_counts() {
            let x = operand(bits(d), 2);
            if let Some(g) = ours {
                // Constants are cached on both sides after the first call.
                let b = ball::Ball::exact(to_ours(&x));
                group.bench_function(BenchmarkId::new("ours", d), |bench| {
                    bench.iter(|| g(black_box(&b), u64::from(bits(d))));
                });
            }
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
    targets = short_products, arithmetic, mag_ops, machine, constants, decimal, functions
}
criterion_main!(benches);
