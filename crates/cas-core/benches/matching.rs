//! Matching under `Flat` and `Orderless` heads (CAS-57, CAS-DOC-4 §3.3).
//!
//! - `split/n`: `x_ + y_` against a sum of n terms: `x` takes the first
//!   term and `y`, the last node, all the others. Linear.
//! - `literal/n`: `f[0] + r___` against the sum: the literal is removed,
//!   then `r` takes the rest. Linear.
//! - `head/n`: `g[x_] + r___` against the sum and one `g[1]`: only terms
//!   with head `g` are tried. Linear.
//! - `worst/n`: three sequences under an `Orderless` head and a condition
//!   that never holds, so every way of splitting n arguments is tried
//!   (about 3^n), until the budget of 10^6 steps stops it with
//!   `MatchLimit` (from n = 12), always at the same step: every iteration
//!   must give what the first gave.
//!
//! Pass a filter to run part of it: `cargo bench -p cas-core --bench
//! matching -- '^worst/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_core::{Attributes, Context, EvalError, Expr, Pattern, Symbol};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap()
}

/// `Plus[f[0], …, f[n - 1]]`, with `extra` terms, evaluated.
fn sum(ctx: &Context, n: i64, extra: &[Expr]) -> Expr {
    let f = Expr::from(Symbol::new("f"));
    let mut terms: Vec<Expr> = (0..n).map(|i| Expr::apply(&f, &[i.into()])).collect();
    terms.extend_from_slice(extra);
    ctx.eval(&Expr::apply(&Symbol::PLUS.into(), &terms))
        .unwrap()
}

fn linear(c: &mut Criterion) {
    let ctx = Context::new();
    let cases = [
        (
            "split",
            "Plus[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
            false,
        ),
        (
            "literal",
            "Plus[f[0], Pattern[r, BlankNullSequence[]]]",
            false,
        ),
        (
            "head",
            "Plus[g[Pattern[x, Blank[]]], Pattern[r, BlankNullSequence[]]]",
            true,
        ),
    ];
    for (name, pattern, with_g) in cases {
        let p = Pattern::with_context(&read(pattern), &ctx).unwrap();
        let mut group = c.benchmark_group(name);
        for n in [100, 10_000] {
            let extra = if with_g {
                vec![read("g[1]")]
            } else {
                Vec::new()
            };
            let e = sum(&ctx, n, &extra);
            group.throughput(Throughput::Elements(n.unsigned_abs()));
            group.bench_with_input(BenchmarkId::from_parameter(n), &e, |b, e| {
                b.iter(|| assert!(ctx.matches(black_box(e), &p).unwrap()));
            });
        }
        group.finish();
    }
}

fn worst(c: &mut Criterion) {
    let mut ctx = Context::new();
    let o = Symbol::new("o");
    ctx.set_attributes(o, Attributes::ORDERLESS).unwrap();
    let p = Pattern::with_context(
        &read(
            "Condition[o[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]], Pattern[z, BlankSequence[]]], False]",
        ),
        &ctx,
    )
    .unwrap();
    let mut group = c.benchmark_group("worst");
    group.sample_size(10);
    for n in [6, 8, 10, 12] {
        let args: Vec<Expr> = (0..n).map(Expr::from).collect();
        let e = Expr::apply(&o.into(), &args);
        let expected = ctx.matches(&e, &p);
        assert!(matches!(expected, Ok(false) | Err(EvalError::MatchLimit)));
        group.bench_with_input(BenchmarkId::from_parameter(n), &e, |b, e| {
            b.iter(|| assert_eq!(ctx.matches(black_box(e), &p), expected));
        });
    }
    group.finish();
}

criterion_group!(benches, linear, worst);
criterion_main!(benches);
