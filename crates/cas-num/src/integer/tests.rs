//! `Integer` checked against GMP. Every result goes through `ToRug`, which
//! also asserts the canonical representation (inline iff it fits an `i64`).

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use proptest::prelude::*;
use rug::Integer as Gmp;
use rug::ops::DivRounding;

use super::{Integer, Rounding};
use crate::Natural;
use crate::testkit::strategy::{self, signed_near_i64_bounds};
use crate::testkit::{Signed, ToRug, all_forms, differential};

fn nonzero() -> BoxedStrategy<Signed> {
    strategy::integer()
        .prop_filter("non-zero", |d| !d.mag().is_zero())
        .boxed()
}

fn gmp_div(a: &Signed, d: &Signed, rounding: Rounding) -> (Gmp, Gmp) {
    let (a, d) = (a.to_rug(), d.to_rug());
    match rounding {
        Rounding::Trunc => a.div_rem(d),
        Rounding::Floor => a.div_rem_floor(d),
        Rounding::Ceil => a.div_rem_ceil(d),
        Rounding::Euclid => a.div_rem_euc(d),
    }
}

fn rounding() -> impl Strategy<Value = Rounding> {
    prop_oneof![
        Just(Rounding::Trunc),
        Just(Rounding::Floor),
        Just(Rounding::Ceil),
        Just(Rounding::Euclid),
    ]
}

fn hash_of(x: &Integer) -> u64 {
    let mut h = DefaultHasher::new();
    x.hash(&mut h);
    h.finish()
}

differential! {
    fn from_parts_round_trip(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            let (negative, mag) = x.clone().into_parts();
            prop_assert_eq!(negative, a.is_negative());
            prop_assert_eq!(Integer::from_parts(negative, mag.clone()), x.clone());
            prop_assert_eq!(x.unsigned_abs(), mag);
            x
        },
        gmp: a.to_rug(),
    }

    fn from_i64(x in any::<i64>()) {
        ours: Integer::from(x),
        gmp: Gmp::from(x),
    }

    fn from_u64(x in any::<u64>()) {
        ours: Integer::from(x),
        gmp: Gmp::from(x),
    }

    fn from_i128(x in any::<i128>()) {
        ours: Integer::from(x),
        gmp: Gmp::from(x),
    }

    fn from_u128(x in any::<u128>()) {
        ours: Integer::from(x),
        gmp: Gmp::from(x),
    }

    fn from_natural(a in strategy::natural()) {
        ours: {
            let x = Integer::from(&a.to_natural());
            prop_assert_eq!(&Integer::from(a.to_natural()), &x);
            x
        },
        gmp: a.to_rug(),
    }

    fn to_natural(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            let n = Natural::try_from(&x).ok();
            prop_assert_eq!(&Natural::try_from(x).ok(), &n);
            n
        },
        gmp: (!a.is_negative()).then(|| a.to_rug()),
    }

    fn to_i64(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            prop_assert_eq!(x.to_i64(), i64::try_from(&x).ok());
            i64::try_from(x).ok()
        },
        gmp: a.to_rug().to_i64().map(Gmp::from),
    }

    fn to_u64(a in strategy::integer()) {
        ours: u64::try_from(&a.to_integer()).ok(),
        gmp: a.to_rug().to_u64().map(Gmp::from),
    }

    fn to_i128(a in strategy::integer()) {
        ours: i128::try_from(&a.to_integer()).ok(),
        gmp: a.to_rug().to_i128().map(Gmp::from),
    }

    fn to_u128(a in strategy::integer()) {
        ours: u128::try_from(&a.to_integer()).ok(),
        gmp: a.to_rug().to_u128().map(Gmp::from),
    }

    fn to_i32(a in strategy::integer()) {
        ours: i32::try_from(&a.to_integer()).ok(),
        gmp: a.to_rug().to_i32().map(Gmp::from),
    }

    fn cmp(a in strategy::integer(), b in strategy::integer()) {
        ours: a.to_integer().cmp(&b.to_integer()),
        gmp: a.to_rug().cmp(&b.to_rug()),
    }

    fn sign(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            prop_assert_eq!(x.is_negative(), a.to_rug() < 0);
            prop_assert_eq!(x.is_positive(), a.to_rug() > 0);
            prop_assert_eq!(x.is_zero(), a.to_rug() == 0);
            prop_assert_eq!(x.is_one(), a.to_rug() == 1);
            x.sign()
        },
        gmp: a.to_rug().cmp0(),
    }

    fn neg(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            let n = -&x;
            prop_assert_eq!(&(-x), &n);
            n
        },
        gmp: -a.to_rug(),
    }

    fn abs(a in strategy::integer()) {
        ours: a.to_integer().abs(),
        gmp: a.to_rug().abs(),
    }

    fn add(a in strategy::integer(), b in strategy::integer()) {
        ours: all_forms!(a.to_integer(), +, +=, b.to_integer()),
        gmp: a.to_rug() + b.to_rug(),
    }

    fn sub(a in strategy::integer(), b in strategy::integer()) {
        ours: all_forms!(a.to_integer(), -, -=, b.to_integer()),
        gmp: a.to_rug() - b.to_rug(),
    }

    fn mul(a in strategy::integer(), b in strategy::integer()) {
        ours: all_forms!(a.to_integer(), *, *=, b.to_integer()),
        gmp: a.to_rug() * b.to_rug(),
    }

    fn div(a in strategy::integer(), d in nonzero()) {
        ours: all_forms!(a.to_integer(), /, /=, d.to_integer()),
        gmp: a.to_rug().div_trunc(d.to_rug()),
    }

    fn rem(a in strategy::integer(), d in nonzero()) {
        ours: all_forms!(a.to_integer(), %, %=, d.to_integer()),
        gmp: a.to_rug().div_rem(d.to_rug()).1,
    }

    fn div_rem(a in strategy::integer(), d in nonzero(), rounding in rounding()) {
        ours: {
            let (x, y) = (a.to_integer(), d.to_integer());
            let r = x.checked_div_rem(&y, rounding).unwrap();
            let named = match rounding {
                Rounding::Trunc => x.div_rem_trunc(&y),
                Rounding::Floor => x.div_rem_floor(&y),
                Rounding::Ceil => x.div_rem_ceil(&y),
                Rounding::Euclid => x.div_rem_euclid(&y),
            };
            prop_assert_eq!(&named, &r);
            r
        },
        gmp: gmp_div(&a, &d, rounding),
    }

    fn checked_div_rem_by_zero(a in strategy::integer(), rounding in rounding()) {
        ours: a.to_integer().checked_div_rem(&Integer::ZERO, rounding),
        gmp: None::<(Gmp, Gmp)>,
    }

    fn div_exact(q in strategy::integer(), d in nonzero()) {
        ours: (q.to_integer() * d.to_integer()).div_exact(&d.to_integer()),
        gmp: q.to_rug(),
    }

    fn is_divisible_by(a in strategy::integer(), d in strategy::integer()) {
        ours: a.to_integer().is_divisible_by(&d.to_integer()),
        gmp: a.to_rug().is_divisible(&d.to_rug()),
    }

    fn is_divisible_by_shaped(q in strategy::integer(), d in nonzero(), plus_one in any::<bool>()) {
        ours: (q.to_integer() * d.to_integer() + Integer::from(u8::from(plus_one)))
            .is_divisible_by(&d.to_integer()),
        gmp: (q.to_rug() * d.to_rug() + u32::from(plus_one)).is_divisible(&d.to_rug()),
    }
}

proptest! {
    #[test]
    fn div_rem_identity(a in strategy::integer(), d in nonzero(), rounding in rounding()) {
        let (x, y) = (a.to_integer(), d.to_integer());
        let (q, r) = x.checked_div_rem(&y, rounding).unwrap();
        prop_assert_eq!(&q * &y + &r, x);
        prop_assert!(r.unsigned_abs() < y.unsigned_abs());
    }

    /// Values computed through the boxed path equal, and hash like, the same
    /// values built directly.
    #[test]
    fn eq_and_hash_are_by_value(a in signed_near_i64_bounds(), b in signed_near_i64_bounds()) {
        let (x, y) = (a.to_integer(), b.to_integer());
        let big = Integer::from(u128::MAX);
        let roundabout = &(&(&x + &big) - &big) + &y;
        let direct = &x + &y;
        prop_assert_eq!(&roundabout, &direct);
        prop_assert_eq!(hash_of(&roundabout), hash_of(&direct));
    }
}

/// Every pair of values around the promotion edges, for every operation:
/// exhaustive where proptest only samples.
#[test]
fn promotion_edges_exhaustive() {
    let mut values = Vec::new();
    for base in [0i128, 1 << 31, 1 << 32, 1 << 62, 1 << 63, 1 << 64] {
        for offset in -2..=2 {
            values.push(base + offset);
            values.push(-(base + offset));
        }
    }
    values.extend([i128::from(i64::MIN), i128::from(i64::MAX), 3, -3, 7, -7]);
    values.sort_unstable();
    values.dedup();

    for &a in &values {
        let (x, ga) = (Integer::from(a), Gmp::from(a));
        assert_eq!((-&x).to_rug(), -ga.clone(), "-{a}");
        for &b in &values {
            let (y, gb) = (Integer::from(b), Gmp::from(b));
            assert_eq!((&x + &y).to_rug(), Gmp::from(&ga + &gb), "{a} + {b}");
            assert_eq!((&x - &y).to_rug(), Gmp::from(&ga - &gb), "{a} - {b}");
            assert_eq!((&x * &y).to_rug(), Gmp::from(&ga * &gb), "{a} * {b}");
            assert_eq!(x.cmp(&y), ga.cmp(&gb), "{a} cmp {b}");
            if b == 0 {
                continue;
            }
            let (sa, sb) = (Signed::from_rug(&ga), Signed::from_rug(&gb));
            for rounding in [
                Rounding::Trunc,
                Rounding::Floor,
                Rounding::Ceil,
                Rounding::Euclid,
            ] {
                assert_eq!(
                    x.checked_div_rem(&y, rounding).to_rug(),
                    Some(gmp_div(&sa, &sb, rounding)),
                    "{a} / {b} ({rounding:?})"
                );
            }
        }
    }
}

#[test]
fn i64_min_over_minus_one_promotes() {
    let min = Integer::from(i64::MIN);
    let q = &min / &Integer::from(-1);
    assert_eq!(q, Integer::from(1u64 << 63));
    assert!(q.to_i64().is_none());
    assert_eq!(&min % &Integer::from(-1), Integer::ZERO);
    assert_eq!(min.div_exact(&Integer::from(-1)), q);
    assert!(min.is_divisible_by(&Integer::from(-1)));
    // And back: -(2^63) demotes to the inline form.
    assert_eq!((-q).to_i64(), Some(i64::MIN));
}

#[test]
fn layout_is_sixteen_bytes() {
    assert_eq!(size_of::<Integer>(), 16);
}

#[test]
fn hex_formatting() {
    assert_eq!(format!("{:?}", Integer::from(-255)), "-255");
    assert_eq!(format!("{:x}", Integer::from(u128::MAX)), "f".repeat(32));
    assert_eq!(
        format!("{:#X}", Integer::from(-(1i128 << 64))),
        "-0x10000000000000000"
    );
    assert_eq!(format!("{:x}", Integer::ZERO), "0");
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn div_by_zero_panics() {
    let _ = Integer::ONE / Integer::ZERO;
}

/// Debug builds check `div_exact`'s contract; release builds return an
/// unspecified value (see the `Natural` tests).
#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "division is not exact")]
fn div_exact_panics_when_inexact() {
    let _ = Integer::from(-7).div_exact(&Integer::from(2));
}

#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "division is not exact")]
fn div_exact_big_panics_when_inexact() {
    let _ = Integer::from(i128::MIN + 1).div_exact(&Integer::from(1i128 << 64));
}

#[test]
fn zero_divides_only_zero() {
    assert!(Integer::ZERO.is_divisible_by(&Integer::ZERO));
    assert!(!Integer::from(u128::MAX).is_divisible_by(&Integer::ZERO));
}
