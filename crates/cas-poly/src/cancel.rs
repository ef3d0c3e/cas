//! Stopping long computations (CAS-DOC-8 §11.4).

use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};

/// What long loops (products, and later Hensel lifting and recombination)
/// poll to stop early. cas-core implements it over the session's interrupt
/// flag, so that an interrupted session stops inside them.
pub trait Cancel {
    /// Whether the computation should stop now.
    fn cancelled(&self) -> bool;
}

/// The default: never cancelled.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub struct Never;

impl Cancel for Never {
    #[inline]
    fn cancelled(&self) -> bool {
        false
    }
}

/// A flag set from another thread, as cas-core's interrupt is.
impl Cancel for AtomicBool {
    #[inline]
    fn cancelled(&self) -> bool {
        self.load(Ordering::Relaxed)
    }
}

/// A computation stopped because its [`Cancel`] said so.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Cancelled;

impl fmt::Display for Cancelled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("cancelled")
    }
}

impl std::error::Error for Cancelled {}
