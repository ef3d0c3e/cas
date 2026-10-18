//! Unsigned arbitrary-precision integers.

mod bits;
mod convert;
mod div;
mod gcd;
pub(crate) mod limb_vec;
mod mul;
mod ops;
mod pow;
mod radix;
mod root;

#[cfg(test)]
mod tests;

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

use crate::{Limb, kernel};
use limb_vec::LimbVec;
pub(crate) use radix::check_radix;

/// Scratch limbs that kernels called on short operands take from the stack
/// instead of the heap.
const STACK_SCRATCH: usize = 32;

/// Runs `f` on `len` zeroed limbs, on the stack up to [`STACK_SCRATCH`].
fn with_scratch<R>(len: usize, f: impl FnOnce(&mut [Limb]) -> R) -> R {
    if len <= STACK_SCRATCH {
        f(&mut [0; STACK_SCRATCH][..len])
    } else {
        f(&mut vec![0; len])
    }
}

/// An arbitrary-precision natural number (non-negative integer).
///
/// Stored as little-endian 64-bit limbs, inline up to 256 bits and on the heap
/// beyond. The representation is normalised (no zero most-significant limb,
/// zero has no limbs), so equality, ordering and hashing are by value.
///
/// Arithmetic follows `std` conventions: operators are implemented for every
/// owned/borrowed combination and reuse owned buffers where possible;
/// subtraction below zero panics, and [`Natural::checked_sub`] returns `None`
/// instead.
#[derive(Clone, Default)]
pub struct Natural {
    limbs: LimbVec,
}

impl Natural {
    /// Zero.
    pub const ZERO: Natural = Natural {
        limbs: LimbVec::new(),
    };

    /// Zero.
    #[must_use]
    pub const fn new() -> Self {
        Self::ZERO
    }

    /// One.
    #[must_use]
    pub fn one() -> Self {
        Self::from(1u64)
    }

    /// Builds a value from little-endian limbs; zero high limbs are allowed.
    #[must_use]
    pub fn from_limbs(limbs: &[Limb]) -> Self {
        let n = kernel::normalized_len(limbs);
        Self::from_limb_vec(LimbVec::from_slice(&limbs[..n]))
    }

    /// Builds a value from little-endian limbs, reusing the vector's buffer
    /// when the value does not fit inline; zero high limbs are allowed.
    #[must_use]
    pub fn from_vec(mut limbs: Vec<Limb>) -> Self {
        limbs.truncate(kernel::normalized_len(&limbs));
        Self::from_limb_vec(LimbVec::from_vec(limbs))
    }

    /// `2^k`.
    ///
    /// # Panics
    ///
    /// If the result would not fit in memory.
    #[must_use]
    pub fn power_of_two(k: u64) -> Self {
        let (limb, bit) = bits::split_bit_index(k);
        let mut limbs = LimbVec::zeroed(limb + 1);
        limbs.as_mut_slice()[limb] = 1 << bit;
        Self::from_limb_vec(limbs)
    }

    /// Little-endian limbs, without zero high limbs (empty for zero).
    #[must_use]
    pub fn limbs(&self) -> &[Limb] {
        self.limbs.as_slice()
    }

    /// Number of limbs (0 for zero).
    #[must_use]
    pub fn limb_len(&self) -> usize {
        self.limbs.len()
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.limb_len() == 0
    }

    #[must_use]
    pub fn is_one(&self) -> bool {
        self.limbs() == [1]
    }

    /// Releases spare capacity; values of up to 256 bits move back inline.
    pub fn shrink_to_fit(&mut self) {
        self.limbs.shrink_to_fit();
    }

    /// Wraps already-normalised storage.
    fn from_limb_vec(limbs: LimbVec) -> Self {
        let n = Natural { limbs };
        n.debug_check();
        n
    }

    /// Drops zero high limbs after an operation that may have produced them.
    fn normalize(&mut self) {
        let n = kernel::normalized_len(self.limbs());
        self.limbs.truncate(n);
    }

    #[inline]
    fn debug_check(&self) {
        debug_assert!(
            kernel::is_normalized(self.limbs()),
            "unnormalised Natural, limbs {:x?}",
            self.limbs()
        );
    }
}

impl PartialEq for Natural {
    fn eq(&self, other: &Self) -> bool {
        self.limbs() == other.limbs()
    }
}

impl Eq for Natural {}

impl Hash for Natural {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.limbs().hash(state);
    }
}

impl Ord for Natural {
    fn cmp(&self, other: &Self) -> Ordering {
        kernel::cmp(self.limbs(), other.limbs())
    }
}

impl PartialOrd for Natural {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::LowerHex for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_hex(self, f, false)
    }
}

impl fmt::UpperHex for Natural {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_hex(self, f, true)
    }
}

fn fmt_hex(n: &Natural, f: &mut fmt::Formatter<'_>, upper: bool) -> fmt::Result {
    use fmt::Write as _;
    let mut s = String::with_capacity(16 * n.limb_len().max(1));
    match n.limbs().split_last() {
        None => s.push('0'),
        Some((top, rest)) if upper => {
            write!(s, "{top:X}")?;
            rest.iter().rev().try_for_each(|l| write!(s, "{l:016X}"))?;
        }
        Some((top, rest)) => {
            write!(s, "{top:x}")?;
            rest.iter().rev().try_for_each(|l| write!(s, "{l:016x}"))?;
        }
    }
    f.pad_integral(true, "0x", &s)
}
