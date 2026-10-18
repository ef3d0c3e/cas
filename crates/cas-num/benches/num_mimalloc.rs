//! The benchmarks of `num.rs` with mimalloc as the global allocator, which
//! cas-num's binaries are meant to use (CAS-DOC-2 §4.2). GMP's own buffers
//! come from C `malloc` in both: `rug` does not route them through Rust's
//! allocator, so only our side changes.
#![allow(missing_docs, reason = "benchmark binary")]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

include!("common/suite.rs");
