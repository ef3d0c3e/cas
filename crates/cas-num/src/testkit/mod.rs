//! Differential test harness: every operation is checked against GMP (via `rug`).
//!
//! Test-only (`#[cfg(test)]`); `rug` is a dev-dependency and never reaches a
//! runtime build. See the crate README for usage.
//!
//! - [`model`]: [`Mag`] and [`Signed`], plain little-endian limb vectors used as
//!   the neutral exchange format between our types and GMP.
//! - [`strategy`]: proptest strategies with adversarial shapes.
//! - [`ToRug`] / [`check_same`] / [`differential!`]: comparison plumbing.
//!
//! Conversions go through limb slices, never through strings, so the
//! conversion itself is not one of the things under test.

pub(crate) mod model;
pub(crate) mod strategy;

mod selftest;

use std::cmp::Ordering;
use std::fmt::Debug;

use proptest::test_runner::TestCaseError;
use rug::Integer;

pub(crate) use model::{Mag, Signed};

use crate::Natural;

/// Maps a value of ours to the equivalent GMP value, for comparison.
///
/// Implemented by every cas-num type and by the primitive types operations
/// return, so `ours` and `gmp` sides of a differential test can be compared
/// whatever their shape (a number, a `bool`, a quotient/remainder pair…).
pub(crate) trait ToRug {
    type Rug: PartialEq + Debug;

    fn to_rug(&self) -> Self::Rug;
}

macro_rules! to_rug_via_integer {
    ($($t:ty),*) => {$(
        impl ToRug for $t {
            type Rug = Integer;

            fn to_rug(&self) -> Integer {
                Integer::from(*self)
            }
        }
    )*};
}

to_rug_via_integer!(
    u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize
);

macro_rules! to_rug_identity {
    ($($t:ty),*) => {$(
        impl ToRug for $t {
            type Rug = $t;

            fn to_rug(&self) -> $t {
                *self
            }
        }
    )*};
}

to_rug_identity!(bool, Ordering);

impl<T: ToRug + ?Sized> ToRug for &T {
    type Rug = T::Rug;

    fn to_rug(&self) -> T::Rug {
        (**self).to_rug()
    }
}

impl<T: ToRug> ToRug for Option<T> {
    type Rug = Option<T::Rug>;

    fn to_rug(&self) -> Self::Rug {
        self.as_ref().map(ToRug::to_rug)
    }
}

impl<A: ToRug, B: ToRug> ToRug for (A, B) {
    type Rug = (A::Rug, B::Rug);

    fn to_rug(&self) -> Self::Rug {
        (self.0.to_rug(), self.1.to_rug())
    }
}

impl ToRug for Natural {
    type Rug = Integer;

    fn to_rug(&self) -> Integer {
        Integer::from_digits(self.limbs(), rug::integer::Order::Lsf)
    }
}

impl Mag {
    /// The same value as a [`Natural`].
    pub(crate) fn to_natural(&self) -> Natural {
        Natural::from_limbs(self.limbs())
    }
}

/// Fails the current proptest case unless `ours` and `gmp` denote the same value.
///
/// `gmp` may be anything convertible into the GMP-side type, which covers
/// rug's incomplete-computation values (`Integer::from(&a + &b)` is not needed,
/// `&a + &b` is enough when the target is `Integer`).
pub(crate) fn check_same<O, G>(ours: O, gmp: G) -> Result<(), TestCaseError>
where
    O: ToRug + Debug,
    G: Into<O::Rug>,
{
    let expected = gmp.into();
    let actual = ours.to_rug();
    if actual == expected {
        Ok(())
    } else {
        Err(TestCaseError::fail(format!(
            "differs from GMP\n  ours:     {ours:?}\n  as GMP:   {actual:?}\n  expected: {expected:?}"
        )))
    }
}

/// Declares property tests asserting that an operation of ours agrees with GMP.
///
/// ```ignore
/// differential! {
///     /// Doc comments and attributes are forwarded.
///     fn add(a in strategy::natural(), b in strategy::natural()) {
///         ours: a.to_ours() + b.to_ours(),
///         gmp: a.to_rug() + b.to_rug(),
///     }
/// }
/// ```
///
/// Each entry becomes a `#[test]` running proptest over the given strategies.
/// `ours` must implement [`ToRug`] + `Debug`; `gmp` must convert into its
/// `ToRug::Rug`. A leading `#![proptest_config(..)]` is forwarded to proptest.
macro_rules! differential {
    (
        $(#![$config:meta])?
        $(
            $(#[$meta:meta])*
            fn $name:ident($($arg:pat in $strat:expr),+ $(,)?) {
                ours: $ours:expr,
                gmp: $gmp:expr $(,)?
            }
        )*
    ) => {
        ::proptest::proptest! {
            $(#![$config])?
            $(
                $(#[$meta])*
                #[test]
                fn $name($($arg in $strat),+) {
                    $crate::testkit::check_same($ours, $gmp)?;
                }
            )*
        }
    };
}

pub(crate) use differential;
