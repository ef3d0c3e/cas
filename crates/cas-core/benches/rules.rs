//! Rule dispatch and loading (CAS-56, CAS-DOC-4 §10).
//!
//! - `exact/n`: `f[k]` against n exact definitions `f[i] = i`: one hash
//!   lookup, flat across n.
//! - `filtered/n`: `f[n - 1, y]` against n rules `f[i, x_] := x`, where the
//!   first-argument filter skips all but the last.
//! - `matched/n`: `f[y, n - 1]` against n rules `f[x_, i] := x`, which the
//!   filter cannot tell apart: the matcher tries each.
//! - `load/…/n`: n definitions into a new context, exact or patterned on one
//!   head, or on n heads; `shared` keeps every state it passes through, as a
//!   session with a history would. Linear in n.
//! - `up/…`: a sum of n terms evaluated with no up-rule, and with one up-rule
//!   on a symbol that is not in it.
//!
//! Pass a filter to run part of it: `cargo bench -p cas-core --bench rules
//! -- '^load/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_core::{Context, Expr, Symbol};
use criterion::{BatchSize, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap()
}

fn call(f: Symbol, args: &[Expr]) -> Expr {
    Expr::apply(&f.into(), args)
}

/// `Pattern[x, Blank[]]`.
fn blank() -> Expr {
    read("Pattern[x, Blank[]]")
}

fn dispatch(c: &mut Criterion) {
    let f = Symbol::new("f");
    let (x, y) = (Expr::from(Symbol::new("x")), Expr::from(Symbol::new("y")));
    let mut group = c.benchmark_group("exact");
    for n in [10, 1_000, 100_000] {
        let mut ctx = Context::new();
        for i in 0..n {
            ctx.define(&call(f, &[i.into()]), &i.into()).unwrap();
        }
        let subject = call(f, &[(n / 2).into()]);
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &subject,
            |bencher, subject| {
                bencher.iter(|| ctx.eval(black_box(subject)).unwrap());
            },
        );
    }
    group.finish();
    let mut group = c.benchmark_group("filtered");
    for n in [10, 100] {
        let mut ctx = Context::new();
        for i in 0..n {
            ctx.define_delayed(&call(f, &[i.into(), blank()]), &x)
                .unwrap();
        }
        let subject = call(f, &[(n - 1).into(), y.clone()]);
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &subject,
            |bencher, subject| {
                bencher.iter(|| ctx.eval(black_box(subject)).unwrap());
            },
        );
    }
    group.finish();
    let mut group = c.benchmark_group("matched");
    for n in [10, 100] {
        let mut ctx = Context::new();
        for i in 0..n {
            ctx.define_delayed(&call(f, &[blank(), i.into()]), &x)
                .unwrap();
        }
        let subject = call(f, &[y.clone(), (n - 1).into()]);
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &subject,
            |bencher, subject| {
                bencher.iter(|| ctx.eval(black_box(subject)).unwrap());
            },
        );
    }
    group.finish();
}

/// What `load` defines: each left-hand side with its right-hand side.
type Definitions = Vec<(Expr, Expr)>;

/// The `i`-th definition of a kind.
type Definition = fn(i64) -> (Expr, Expr);

fn load(c: &mut Criterion) {
    let sizes = [1_000, 10_000, 100_000];
    let kinds: [(&str, Definition); 3] = [
        ("exact", |i| (call(Symbol::new("f"), &[i.into()]), i.into())),
        ("patterned", |i| {
            (
                call(Symbol::new("f"), &[i.into(), blank()]),
                Symbol::new("x").into(),
            )
        }),
        ("heads", |i| {
            (
                call(Symbol::new(&format!("h{i}")), &[blank()]),
                Symbol::new("x").into(),
            )
        }),
    ];
    let mut group = c.benchmark_group("load");
    group.sample_size(10);
    for (kind, definition) in kinds {
        for n in sizes {
            let definitions: Definitions = (0..n).map(definition).collect();
            group.throughput(Throughput::Elements(n.unsigned_abs()));
            group.bench_with_input(BenchmarkId::new(kind, n), &definitions, |b, d| {
                b.iter_batched(
                    Context::new,
                    |mut ctx| {
                        for (lhs, rhs) in d {
                            ctx.define_delayed(lhs, rhs).unwrap();
                        }
                        ctx
                    },
                    BatchSize::LargeInput,
                );
            });
        }
    }
    // Each state kept: every definition copies what it changes.
    for n in sizes {
        let definitions: Definitions = (0..n).map(kinds[2].1).collect();
        group.throughput(Throughput::Elements(n.unsigned_abs()));
        group.bench_with_input(BenchmarkId::new("shared", n), &definitions, |b, d| {
            b.iter_batched(
                || (Context::new(), Vec::with_capacity(d.len())),
                |(mut ctx, mut history)| {
                    for (lhs, rhs) in d {
                        history.push(ctx.clone());
                        ctx.define_delayed(lhs, rhs).unwrap();
                    }
                    (ctx, history)
                },
                BatchSize::LargeInput,
            );
        });
    }
    group.finish();
}

/// `Plus[f[n - 1], …, f[0]]`.
fn sum(n: i64) -> Expr {
    let f = Symbol::new("f");
    let terms: Vec<Expr> = (0..n).rev().map(|i| call(f, &[i.into()])).collect();
    call(Symbol::PLUS, &terms)
}

fn up(c: &mut Criterion) {
    let n = 10_000;
    let e = sum(n);
    let mut group = c.benchmark_group("up");
    group.throughput(Throughput::Elements(n.unsigned_abs()));
    let plain = Context::new();
    let mut with_up = Context::new();
    let phi = Symbol::new("phi");
    with_up
        .define_up(phi, &read("Power[phi, 2]"), &read("Plus[phi, 1]"))
        .unwrap();
    for (name, ctx) in [("none", plain), ("one", with_up)] {
        group.bench_with_input(BenchmarkId::new(name, n), &e, |b, e| {
            // A new state each time, so that no stamp answers: any change
            // makes one.
            b.iter_batched(
                || {
                    let mut fresh = ctx.clone();
                    fresh.set_limits(ctx.limits());
                    fresh
                },
                |fresh| fresh.eval(black_box(e)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, dispatch, load, up);
criterion_main!(benches);
