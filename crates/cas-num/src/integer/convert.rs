//! Conversions between [`Integer`], [`Natural`] and primitive integers.

use super::{Integer, Repr};
use crate::{Natural, OutOfRangeError};

macro_rules! from_small {
    ($($t:ty),*) => {$(
        impl From<$t> for Integer {
            fn from(x: $t) -> Integer {
                Integer(Repr::Small(i64::from(x)))
            }
        }
    )*};
}

from_small!(i8, i16, i32, i64, u8, u16, u32);

impl From<isize> for Integer {
    fn from(x: isize) -> Integer {
        Integer::from(i64::try_from(x).expect("isize fits i64"))
    }
}

impl From<u64> for Integer {
    fn from(x: u64) -> Integer {
        Integer::from(u128::from(x))
    }
}

impl From<usize> for Integer {
    fn from(x: usize) -> Integer {
        Integer::from(u64::try_from(x).expect("usize fits u64"))
    }
}

impl From<u128> for Integer {
    fn from(x: u128) -> Integer {
        match i64::try_from(x) {
            Ok(x) => Integer(Repr::Small(x)),
            Err(_) => Integer::from_parts(false, Natural::from(x)),
        }
    }
}

impl From<i128> for Integer {
    fn from(x: i128) -> Integer {
        match i64::try_from(x) {
            Ok(x) => Integer(Repr::Small(x)),
            Err(_) => Integer::from_parts(x < 0, Natural::from(x.unsigned_abs())),
        }
    }
}

impl From<Natural> for Integer {
    fn from(x: Natural) -> Integer {
        Integer::from_parts(false, x)
    }
}

impl From<&Natural> for Integer {
    fn from(x: &Natural) -> Integer {
        match u64::try_from(x) {
            Ok(x) => Integer::from(x),
            Err(_) => Integer::from_parts(false, x.clone()),
        }
    }
}

impl TryFrom<Integer> for Natural {
    type Error = OutOfRangeError;

    /// Fails on negative values.
    fn try_from(x: Integer) -> Result<Natural, OutOfRangeError> {
        match x.into_parts() {
            (false, mag) => Ok(mag),
            (true, _) => Err(OutOfRangeError(())),
        }
    }
}

impl TryFrom<&Integer> for Natural {
    type Error = OutOfRangeError;

    /// Fails on negative values.
    fn try_from(x: &Integer) -> Result<Natural, OutOfRangeError> {
        if x.is_negative() {
            Err(OutOfRangeError(()))
        } else {
            Ok(x.unsigned_abs())
        }
    }
}

impl Integer {
    /// The value as an `i128`, if it fits.
    fn to_i128(&self) -> Option<i128> {
        match &self.0 {
            Repr::Small(x) => Some(i128::from(*x)),
            Repr::Big(b) => {
                let m = u128::try_from(&b.mag).ok()?;
                if b.negative {
                    0i128.checked_sub_unsigned(m)
                } else {
                    i128::try_from(m).ok()
                }
            }
        }
    }
}

macro_rules! try_into_primitive {
    ($($t:ty),*) => {$(
        impl TryFrom<&Integer> for $t {
            type Error = OutOfRangeError;

            fn try_from(x: &Integer) -> Result<$t, OutOfRangeError> {
                x.to_i128().and_then(|v| <$t>::try_from(v).ok()).ok_or(OutOfRangeError(()))
            }
        }

        impl TryFrom<Integer> for $t {
            type Error = OutOfRangeError;

            fn try_from(x: Integer) -> Result<$t, OutOfRangeError> {
                <$t>::try_from(&x)
            }
        }
    )*};
}

try_into_primitive!(u8, u16, u32, u64, usize, i8, i16, i32, i64, i128, isize);

// u128 values above i128::MAX do not go through to_i128.
impl TryFrom<&Integer> for u128 {
    type Error = OutOfRangeError;

    fn try_from(x: &Integer) -> Result<u128, OutOfRangeError> {
        if x.is_negative() {
            return Err(OutOfRangeError(()));
        }
        u128::try_from(&*x.sign_mag().1)
    }
}

impl TryFrom<Integer> for u128 {
    type Error = OutOfRangeError;

    fn try_from(x: Integer) -> Result<u128, OutOfRangeError> {
        u128::try_from(&x)
    }
}
