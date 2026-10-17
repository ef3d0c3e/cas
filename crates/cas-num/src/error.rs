//! Error types.

use std::fmt;

/// A conversion failed because the value is out of range for the target type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct OutOfRangeError(pub(crate) ());

impl fmt::Display for OutOfRangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("value out of range for the target type")
    }
}

impl std::error::Error for OutOfRangeError {}
