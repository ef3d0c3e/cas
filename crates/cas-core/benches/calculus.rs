//! `D` as system rules against the size of what it differentiates (CAS-59,
//! CAS-DOC-4 §10). Each iteration evaluates in a new context, with no
//! stamps to reuse.
//!
//! - `d/polynomial`: `D[1 + 2 x + 3 x^2 + …, x]` with `n` terms: linear, the
//!   sum rule mapping `D` over all the terms at once.
//! - `d/nested`: `D[Sin[Sin[…Sin[x]…]], x]`, `n` deep, the chain rule at each
//!   level: a product of `n` factors, on a thread of its own with a large
//!   stack (`Limits::deep`), whose start the smallest size mostly measures.
//!   About cubic: each level sorts its factors again, and the canonical
//!   order compares `Cos[Sin[…]]` factors down their nesting.
//! - `d/product`: `D[(x + 1) (x + 2) … (x + n), x]`: the product rule splits
//!   off one factor at a time, so it tends to quadratic; also on a thread of
//!   its own.
//!
//! Pass a filter to run part of it: `cargo bench -p cas-core --bench
//! calculus -- '^d/nested/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_core::{Context, Expr, Limits, Symbol};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

fn x() -> Expr {
    Symbol::new("x").into()
}

fn apply(head: Symbol, args: &[Expr]) -> Expr {
    Expr::apply(&head.into(), args)
}

/// `D[f, x]`.
fn d(f: Expr) -> Expr {
    apply(Symbol::D, &[f, x()])
}

/// `1 + 2 x + … + n x^(n - 1)`.
fn polynomial(n: i64) -> Expr {
    let terms: Vec<Expr> = (0..n)
        .map(|i| {
            let power = apply(Symbol::POWER, &[x(), i.into()]);
            apply(Symbol::TIMES, &[(i + 1).into(), power])
        })
        .collect();
    apply(Symbol::PLUS, &terms)
}

/// `Sin[Sin[…Sin[x]…]]`, `n` deep.
fn nested(n: usize) -> Expr {
    (0..n).fold(x(), |e, _| apply(Symbol::SIN, &[e]))
}

/// `(x + 1) (x + 2) … (x + n)`.
fn product(n: i64) -> Expr {
    let factors: Vec<Expr> = (1..=n)
        .map(|i| apply(Symbol::PLUS, &[x(), i.into()]))
        .collect();
    apply(Symbol::TIMES, &factors)
}

fn deep() -> Context {
    let mut ctx = Context::new();
    ctx.set_limits(Limits::deep());
    ctx
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("d");
    for n in [10, 100, 1000, 10_000] {
        let e = d(polynomial(n));
        group.throughput(Throughput::Elements(n.unsigned_abs()));
        group.bench_with_input(BenchmarkId::new("polynomial", n), &e, |b, e| {
            b.iter(|| Context::new().eval(black_box(e)).unwrap());
        });
    }
    for n in [10, 30, 100, 300] {
        let e = d(nested(n));
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("nested", n), &e, |b, e| {
            b.iter(|| deep().eval(black_box(e)).unwrap());
        });
    }
    for n in [10, 30, 100] {
        let e = d(product(n));
        group.throughput(Throughput::Elements(n.unsigned_abs()));
        group.bench_with_input(BenchmarkId::new("product", n), &e, |b, e| {
            b.iter(|| deep().eval(black_box(e)).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
