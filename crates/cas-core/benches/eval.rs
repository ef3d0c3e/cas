//! Evaluation time against the size of the expression (CAS-53, CAS-DOC-4
//! §10). The expression is a sum of distinct terms `f[i]`, out of order.
//!
//! - `eval/terms` evaluates it in a new state, with no stamps to reuse:
//!   linear in the terms, plus the sort.
//! - `reeval/terms` evaluates the result again in the same state: the root's
//!   fixpoint stamp answers at once, so the time is flat across sizes.
//!
//! Pass a filter to run part of it: `cargo bench -p cas-core --bench eval --
//! '^reeval/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_core::{Context, Expr, Symbol};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

/// `Plus[f[n - 1], …, f[1], f[0]]`.
fn sum(n: i64) -> Expr {
    let f = Expr::from(Symbol::new("f"));
    let terms: Vec<Expr> = (0..n).rev().map(|i| Expr::apply(&f, &[i.into()])).collect();
    Expr::apply(&Symbol::PLUS.into(), &terms)
}

fn bench(c: &mut Criterion) {
    let sizes = [1_000, 10_000, 100_000];
    let mut group = c.benchmark_group("eval");
    for n in sizes {
        let e = sum(n);
        group.throughput(Throughput::Elements(n.unsigned_abs()));
        group.bench_with_input(BenchmarkId::new("terms", n), &e, |b, e| {
            b.iter(|| Context::new().eval(black_box(e)).unwrap());
        });
    }
    group.finish();
    let mut group = c.benchmark_group("reeval");
    for n in sizes {
        let ctx = Context::new();
        let r = ctx.eval(&sum(n)).unwrap();
        group.bench_with_input(BenchmarkId::new("terms", n), &r, |b, r| {
            b.iter(|| ctx.eval(black_box(r)).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
