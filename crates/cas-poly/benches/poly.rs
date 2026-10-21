//! Polynomial benchmarks (CAS-141, CAS-DOC-8 §13), wall time by criterion,
//! at sizes that run in seconds. The full sizes are the ignored workload
//! tests (`products_workload`, `sparse_workload`, `gcd_workload`), counted
//! in instructions under `perf stat`; their numbers are in
//! `doc/algorithms`.
//!
//! - `dense_mul/BITS/LENGTH`: the product of two dense polynomials over ℤ
//!   of `LENGTH` coefficients of `BITS` bits, of either sign, as
//!   `DenseRing` chooses it: Kronecker substitution from three
//!   coefficients at these sizes. The 100- and 200-bit products are
//!   CAS-70's workload.
//! - `fateman/K`: `f (f + 1)` for `f = (1 + x + y + z + t)^K`, sparse.
//! - `pearce/K`: Monagan and Pearce's `(1 + x + y + 2z^2 + 3t^3 + 5u^5)^K
//!   (1 + u + t + 2z^2 + 3y^3 + 5x^5)^K`, sparse.
//! - `gcd_dense/DEGREE`: the gcd of `f g` and `f h` for pseudo-random `f`,
//!   `g`, `h` of that degree with 64-bit coefficients.
//! - `gcd_sparse/N`: the same for `f = (1 + x + y + z)^N + x`, and `g`, `h`
//!   likewise with `y` and `z`.
//!
//! The factors are built outside the timed loop. Pass a filter to run part
//! of it: `cargo bench -p cas-poly --bench poly -- '^fateman/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use cas_num::algebra::{Integers, Ring};
use cas_num::{Integer, Natural};
use cas_poly::{Dense, DenseRing, Sparse, SparseRing};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

/// Pseudo-random words, splitmix64: the same factors in every run.
struct Words(u64);

impl Words {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// An integer of `bits` bits, of either sign.
    fn integer(&mut self, bits: u64) -> Integer {
        let limbs: Vec<u64> = (0..bits.div_ceil(64)).map(|_| self.next()).collect();
        let n = Natural::from_vec(limbs) >> (bits.div_ceil(64) * 64 - bits);
        Integer::from_parts(self.next() & 1 == 1, n)
    }

    /// A dense polynomial of `length` coefficients of `bits` bits.
    fn dense(&mut self, ring: &DenseRing<Integers>, length: usize, bits: u64) -> Dense<Integer> {
        ring.poly((0..length).map(|_| self.integer(bits)).collect())
    }
}

/// The sum of the variables of `ring`, and 1.
fn one_plus_variables(ring: &SparseRing<Integers>) -> Sparse<Integer> {
    (0..ring.nvars()).fold(ring.one(), |s, v| ring.add(&s, &ring.var(v)))
}

/// `c x_v^e` in the five variables of `ring`.
fn monomial(ring: &SparseRing<Integers>, c: i64, v: usize, e: u64) -> Sparse<Integer> {
    let mut exps = vec![0; ring.nvars()];
    exps[v] = e;
    ring.monomial(Integer::from(c), &exps)
}

fn dense_products(c: &mut Criterion) {
    let ring = DenseRing::new(Integers);
    let mut words = Words(141);
    for bits in [100, 200] {
        let mut group = c.benchmark_group(format!("dense_mul/{bits}"));
        for length in [16, 128, 1024, 8192] {
            let (a, b) = (
                words.dense(&ring, length, bits),
                words.dense(&ring, length, bits),
            );
            group.throughput(Throughput::Elements(length as u64));
            group.bench_with_input(
                BenchmarkId::from_parameter(length),
                &(a, b),
                |bench, (a, b)| {
                    bench.iter(|| ring.mul(black_box(a), black_box(b)));
                },
            );
        }
        group.finish();
    }
}

fn sparse_products(c: &mut Criterion) {
    let mut group = c.benchmark_group("fateman");
    group.sample_size(10);
    let ring = SparseRing::new(Integers, 4);
    for k in [5, 10] {
        let f = ring.pow(&one_plus_variables(&ring), k);
        let g = ring.add(&f, &ring.one());
        group.bench_with_input(BenchmarkId::from_parameter(k), &(f, g), |bench, (f, g)| {
            bench.iter(|| ring.mul(black_box(f), black_box(g)));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("pearce");
    group.sample_size(10);
    let ring = SparseRing::new(Integers, 5);
    let sum = |parts: &[Sparse<Integer>]| parts.iter().fold(ring.one(), |s, p| ring.add(&s, p));
    let m = |c, v, e| monomial(&ring, c, v, e);
    let f1 = sum(&[m(1, 0, 1), m(1, 1, 1), m(2, 2, 2), m(3, 3, 3), m(5, 4, 5)]);
    let g1 = sum(&[m(1, 4, 1), m(1, 3, 1), m(2, 2, 2), m(3, 1, 3), m(5, 0, 5)]);
    for k in [4, 8] {
        let (f, g) = (ring.pow(&f1, k), ring.pow(&g1, k));
        group.bench_with_input(BenchmarkId::from_parameter(k), &(f, g), |bench, (f, g)| {
            bench.iter(|| ring.mul(black_box(f), black_box(g)));
        });
    }
    group.finish();
}

#[allow(
    clippy::many_single_char_names,
    reason = "gcd(f g, f h), as written in mathematics"
)]
fn gcds(c: &mut Criterion) {
    let mut group = c.benchmark_group("gcd_dense");
    group.sample_size(20);
    let ring = DenseRing::new(Integers);
    let mut words = Words(137);
    for degree in [100, 400] {
        let mut random = || words.dense(&ring, degree + 1, 64);
        let (f, g, h) = (random(), random(), random());
        let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &(a, b),
            |bench, (a, b)| {
                bench.iter(|| ring.gcd(black_box(a), black_box(b)));
            },
        );
    }
    group.finish();
    let mut group = c.benchmark_group("gcd_sparse");
    group.sample_size(20);
    let ring = SparseRing::new(Integers, 3);
    for n in [5, 8] {
        let base = ring.pow(&one_plus_variables(&ring), n);
        let (f, g, h) = (
            ring.add(&base, &ring.var(0)),
            ring.add(&base, &ring.var(1)),
            ring.add(&base, &ring.var(2)),
        );
        let (a, b) = (ring.mul(&f, &g), ring.mul(&f, &h));
        group.bench_with_input(BenchmarkId::from_parameter(n), &(a, b), |bench, (a, b)| {
            bench.iter(|| ring.gcd(black_box(a), black_box(b)));
        });
    }
    group.finish();
}

criterion_group!(benches, dense_products, sparse_products, gcds);
criterion_main!(benches);
