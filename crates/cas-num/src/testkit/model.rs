//! Neutral limb-vector values exchanged between our types and GMP.

use std::fmt;

use rug::Integer;
use rug::integer::Order;

use super::ToRug;

/// A natural number as little-endian 64-bit limbs.
///
/// Always normalised: no zero high limb, and zero is the empty vector, so
/// derived equality is value equality.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub(crate) struct Mag(Vec<u64>);

impl Mag {
    /// Builds a value from little-endian limbs, dropping zero high limbs.
    pub(crate) fn from_limbs(limbs: impl Into<Vec<u64>>) -> Self {
        let mut limbs = limbs.into();
        while limbs.last() == Some(&0) {
            limbs.pop();
        }
        Mag(limbs)
    }

    pub(crate) fn from_u64(x: u64) -> Self {
        Self::from_limbs(vec![x])
    }

    /// Little-endian limbs, normalised (empty for zero).
    pub(crate) fn limbs(&self) -> &[u64] {
        &self.0
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The value if it fits a single limb.
    pub(crate) fn to_u64(&self) -> Option<u64> {
        match self.0.as_slice() {
            [] => Some(0),
            [x] => Some(*x),
            _ => None,
        }
    }

    /// Converts a non-negative GMP integer; panics on a negative one.
    pub(crate) fn from_rug(x: &Integer) -> Self {
        assert!(*x >= 0, "Mag::from_rug on negative value {x}");
        Self::from_limbs(x.to_digits::<u64>(Order::Lsf))
    }
}

impl ToRug for Mag {
    type Rug = Integer;

    fn to_rug(&self) -> Integer {
        Integer::from_digits(&self.0, Order::Lsf)
    }
}

/// Hex, most significant limb first, limbs separated by `_`: failure output
/// shows limb boundaries, which is where most bugs live.
impl fmt::Debug for Mag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Some((top, rest)) = self.0.split_last() else {
            return f.write_str("0x0");
        };
        write!(f, "0x{top:x}")?;
        for limb in rest.iter().rev() {
            write!(f, "_{limb:016x}")?;
        }
        write!(f, " ({} limbs)", self.0.len())
    }
}

/// A signed integer: sign and magnitude. Zero is never negative.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub(crate) struct Signed {
    negative: bool,
    mag: Mag,
}

impl Signed {
    pub(crate) fn new(negative: bool, mag: Mag) -> Self {
        Signed {
            negative: negative && !mag.is_zero(),
            mag,
        }
    }

    pub(crate) fn from_i64(x: i64) -> Self {
        Self::new(x < 0, Mag::from_u64(x.unsigned_abs()))
    }

    pub(crate) fn is_negative(&self) -> bool {
        self.negative
    }

    pub(crate) fn mag(&self) -> &Mag {
        &self.mag
    }

    pub(crate) fn from_rug(x: &Integer) -> Self {
        Self::new(*x < 0, Mag::from_limbs(x.to_digits::<u64>(Order::Lsf)))
    }
}

impl ToRug for Signed {
    type Rug = Integer;

    fn to_rug(&self) -> Integer {
        let mag = self.mag.to_rug();
        if self.negative { -mag } else { mag }
    }
}

impl fmt::Debug for Signed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.negative {
            f.write_str("-")?;
        }
        self.mag.fmt(f)
    }
}
