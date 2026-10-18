//! Benchmarks of cas-num against GMP (through `rug`), op by op.
//!
//! Every benchmark has an `ours` and a `gmp` variant on the same inputs, so
//! `cargo bench -p cas-num` prints the two times side by side. Sizes sweep
//! powers of two up to `CAS_BENCH_MAX_LIMBS` limbs (default 4096; the
//! ticket's full sweep is `CAS_BENCH_MAX_LIMBS=100000`), which places the
//! algorithm thresholds on real data. Pass a filter to run part of it, e.g.
//! `cargo bench -p cas-num -- mul/`. The `num_mimalloc` bench runs the same
//! suite with mimalloc as the global allocator.
#![allow(missing_docs, reason = "benchmark binary")]

include!("common/suite.rs");
