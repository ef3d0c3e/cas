//! Conversions between [`Natural`] and primitive integers.

use super::Natural;
use super::limb_vec::LimbVec;
use crate::{Limb, OutOfRangeError};

macro_rules! from_unsigned {
    ($($t:ty),*) => {$(
        impl From<$t> for Natural {
            fn from(x: $t) -> Natural {
                Natural::from_limb_vec(LimbVec::from_pair(Limb::try_from(x).expect("fits a limb"), 0))
            }
        }
    )*};
}

from_unsigned!(u8, u16, u32, u64, usize);

impl From<u128> for Natural {
    #[expect(clippy::cast_possible_truncation, reason = "splitting into limbs")]
    fn from(x: u128) -> Natural {
        let lo = x as Limb;
        let hi = (x >> Limb::BITS) as Limb;
        Natural::from_limb_vec(LimbVec::from_pair(lo, hi))
    }
}

macro_rules! try_from_signed {
    ($($t:ty => $u:ty),*) => {$(
        impl TryFrom<$t> for Natural {
            type Error = OutOfRangeError;

            /// Fails on negative values.
            fn try_from(x: $t) -> Result<Natural, OutOfRangeError> {
                <$u>::try_from(x).map(Natural::from).map_err(|_| OutOfRangeError(()))
            }
        }
    )*};
}

try_from_signed!(i8 => u8, i16 => u16, i32 => u32, i64 => u64, i128 => u128, isize => usize);

impl Natural {
    /// The value as a `u128`, if it fits.
    pub(super) fn to_u128(&self) -> Option<u128> {
        match *self.limbs() {
            [] => Some(0),
            [lo] => Some(u128::from(lo)),
            [lo, hi] => Some(u128::from(hi) << Limb::BITS | u128::from(lo)),
            _ => None,
        }
    }
}

macro_rules! try_into_primitive {
    ($($t:ty),*) => {$(
        impl TryFrom<&Natural> for $t {
            type Error = OutOfRangeError;

            fn try_from(x: &Natural) -> Result<$t, OutOfRangeError> {
                x.to_u128().and_then(|v| <$t>::try_from(v).ok()).ok_or(OutOfRangeError(()))
            }
        }

        impl TryFrom<Natural> for $t {
            type Error = OutOfRangeError;

            fn try_from(x: Natural) -> Result<$t, OutOfRangeError> {
                <$t>::try_from(&x)
            }
        }
    )*};
}

try_into_primitive!(
    u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize
);
