//! The polynomial builtins through evaluation (CAS-141, CAS-140): the
//! conversion to cas-poly, the algorithm, the way back and the evaluation
//! of the result, each call in a new context.
//!
//! - `expand_power/K`: `Expand[(1 + x + y + z)^K]`.
//! - `expand_fateman/K`: `Expand[f (f + 1)]` for `f = (1 + x + y + z + t)^K`,
//!   Fateman's product through expressions.
//! - `together/N`: `Together` of `1/(x + 1) + … + 1/(x + N)`, one
//!   denominator of degree `N`.
//! - `together_cancel/N`: `Together` of the sum for `k = 1, …, N` of
//!   `(x^2 - k^2 y^2)/(x + k y)`, each term cancelling to `x - k y`.
//!
//! The full Fateman product, `K = 20`, is the ignored `poly_workload` test,
//! counted in instructions under `perf stat`. Pass a filter to run part of
//! it: `cargo bench -p cas-core --bench poly -- '^together'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_core::{Context, Expr, Symbol};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap()
}

/// `Together` of `1/(x + 1) + … + 1/(x + n)`.
fn fractions(n: usize) -> Expr {
    let terms: Vec<String> = (1..=n)
        .map(|k| format!("Power[Plus[x, {k}], -1]"))
        .collect();
    read(&format!("Together[Plus[{}]]", terms.join(", ")))
}

/// `Together` of the sum of `(x^2 - k^2 y^2)/(x + k y)`.
fn cancelling(n: usize) -> Expr {
    let terms: Vec<String> = (1..=n)
        .map(|k| {
            format!(
                "Times[Plus[Power[x, 2], Times[-{}, Power[y, 2]]], Power[Plus[x, Times[{k}, y]], -1]]",
                k * k
            )
        })
        .collect();
    read(&format!("Together[Plus[{}]]", terms.join(", ")))
}

/// `Expand[f (f + 1)]` for `f = (1 + x + y + z + t)^k`.
fn fateman(k: usize) -> Expr {
    read(&format!(
        "Expand[Times[Power[Plus[1, x, y, z, t], {k}], Plus[Power[Plus[1, x, y, z, t], {k}], 1]]]"
    ))
}

fn bench(c: &mut Criterion) {
    let eval = |e: &Expr| Context::new().eval(black_box(e)).unwrap();
    let mut group = c.benchmark_group("expand_power");
    for k in [5, 10, 20] {
        let e = read(&format!("Expand[Power[Plus[1, x, y, z], {k}]]"));
        group.bench_with_input(BenchmarkId::from_parameter(k), &e, |b, e| {
            b.iter(|| eval(e));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("expand_fateman");
    group.sample_size(10);
    for k in [5, 10] {
        let e = fateman(k);
        group.bench_with_input(BenchmarkId::from_parameter(k), &e, |b, e| {
            b.iter(|| eval(e));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("together");
    for n in [5, 10, 20] {
        let e = fractions(n);
        group.bench_with_input(BenchmarkId::from_parameter(n), &e, |b, e| {
            b.iter(|| eval(e));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("together_cancel");
    for n in [5, 10, 20] {
        let e = cancelling(n);
        assert_eq!(eval(&e).head(), Some(&Expr::from(Symbol::PLUS)));
        group.bench_with_input(BenchmarkId::from_parameter(n), &e, |b, e| {
            b.iter(|| eval(e));
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
