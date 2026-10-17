//! Arbitrary-precision numerics for the CAS engine.
//!
//! Correctness is the primary goal: every operation is differentially tested
//! against GMP (test-only dependency). See the "Design decisions" page in Lific.

#[cfg(test)]
mod testkit;
