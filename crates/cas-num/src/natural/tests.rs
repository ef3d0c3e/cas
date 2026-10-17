//! `Natural` checked against GMP; every operator is exercised in all its
//! owned/borrowed and in-place forms, which must agree with each other.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use proptest::prelude::*;
use rug::Integer;

use super::Natural;
use super::limb_vec::INLINE;
use crate::Limb;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

/// Evaluates `x op y` in every owned/borrowed form plus both in-place forms,
/// fails unless they all agree, and yields the result.
macro_rules! all_forms {
    ($x:expr, $op:tt, $op_assign:tt, $y:expr) => {{
        let (x, y) = ($x, $y);
        let r = &x $op &y;
        let forms = [
            x.clone() $op &y,
            &x $op y.clone(),
            x.clone() $op y.clone(),
            { let mut z = x.clone(); z $op_assign &y; z },
            { let mut z = x.clone(); z $op_assign y.clone(); z },
        ];
        for (i, f) in forms.iter().enumerate() {
            prop_assert_eq!(f, &r, "form {} differs from &x {} &y", i, stringify!($op));
        }
        r
    }};
}

/// Same for an operation with a primitive right operand (`Limb`, shift amount).
macro_rules! all_forms_prim {
    ($x:expr, $op:tt, $op_assign:tt, $y:expr) => {{
        let (x, y) = ($x, $y);
        let r = &x $op y;
        let forms = [x.clone() $op y, { let mut z = x.clone(); z $op_assign y; z }];
        for (i, f) in forms.iter().enumerate() {
            prop_assert_eq!(f, &r, "form {} differs from &x {} y", i, stringify!($op));
        }
        r
    }};
}

fn hash_of(x: &Natural) -> u64 {
    let mut h = DefaultHasher::new();
    x.hash(&mut h);
    h.finish()
}

/// Lengths on both sides of the inline/heap boundary.
fn around_inline() -> BoxedStrategy<Mag> {
    strategy::natural_with_len(strategy::len_around(&[INLINE, INLINE + 1])).boxed()
}

/// Shift amounts up to a few limbs past the largest generated value, weighted
/// towards multiples of 64 and their neighbours.
fn shift() -> impl Strategy<Value = u64> {
    let max = 64 * (DEFAULT_MAX_LIMBS as u64 + 2);
    prop_oneof![
        2 => (0..=DEFAULT_MAX_LIMBS as u64 + 2, -1i64..=1)
            .prop_map(|(j, d)| (64 * j).saturating_add_signed(d)),
        3 => 0..max,
    ]
}

/// Two values of the same length; in half the cases `b` copies the top
/// `shared` limbs of `a`.
fn same_len_sharing_top() -> impl Strategy<Value = (Mag, Mag)> {
    (1..=DEFAULT_MAX_LIMBS)
        .prop_flat_map(|n| {
            let m = || strategy::natural_with_len(Just(n));
            (m(), m(), 0..=n, any::<bool>())
        })
        .prop_map(|(a, b, shared, share)| {
            if !share {
                return (a, b);
            }
            let mut limbs = b.limbs().to_vec();
            let n = limbs.len();
            limbs[n - shared..].copy_from_slice(&a.limbs()[n - shared..]);
            (a, Mag::from_limbs(limbs))
        })
}

differential! {
    fn from_limbs_normalises(a in strategy::natural(), junk in 0..4usize) {
        ours: {
            let mut limbs = a.limbs().to_vec();
            limbs.resize(limbs.len() + junk, 0);
            let x = Natural::from_limbs(&limbs);
            prop_assert_eq!(&x, &Natural::from_vec(limbs));
            x
        },
        gmp: a.to_rug(),
    }

    fn from_u64(x in any::<u64>()) {
        ours: Natural::from(x),
        gmp: Integer::from(x),
    }

    fn from_u128(x in any::<u128>()) {
        ours: Natural::from(x),
        gmp: Integer::from(x),
    }

    fn try_from_i128(x in any::<i128>()) {
        ours: Natural::try_from(x).ok(),
        gmp: (x >= 0).then(|| Integer::from(x)),
    }

    fn to_u64(a in strategy::natural_up_to(3)) {
        ours: u64::try_from(&a.to_natural()).ok(),
        gmp: a.to_rug().to_u64().map(Integer::from),
    }

    fn to_u128(a in strategy::natural_up_to(3)) {
        ours: u128::try_from(a.to_natural()).ok(),
        gmp: a.to_rug().to_u128().map(Integer::from),
    }

    fn to_i64(a in strategy::natural_up_to(3)) {
        ours: i64::try_from(&a.to_natural()).ok(),
        gmp: a.to_rug().to_i64().map(Integer::from),
    }

    fn to_i128(a in strategy::natural_up_to(3)) {
        ours: i128::try_from(&a.to_natural()).ok(),
        gmp: a.to_rug().to_i128().map(Integer::from),
    }

    fn power_of_two(k in 0..64 * DEFAULT_MAX_LIMBS as u64) {
        ours: Natural::power_of_two(k),
        gmp: Integer::from(1) << u32::try_from(k).unwrap(),
    }

    fn cmp(a in strategy::natural(), b in strategy::natural()) {
        ours: a.to_natural().cmp(&b.to_natural()),
        gmp: a.to_rug().cmp(&b.to_rug()),
    }

    fn add(a in strategy::natural(), b in strategy::natural()) {
        ours: all_forms!(a.to_natural(), +, +=, b.to_natural()),
        gmp: a.to_rug() + b.to_rug(),
    }

    fn add_unbalanced((a, b) in strategy::unbalanced_pair(DEFAULT_MAX_LIMBS)) {
        ours: all_forms!(a.to_natural(), +, +=, b.to_natural()),
        gmp: a.to_rug() + b.to_rug(),
    }

    fn add_around_inline(a in around_inline(), b in around_inline()) {
        ours: all_forms!(a.to_natural(), +, +=, b.to_natural()),
        gmp: a.to_rug() + b.to_rug(),
    }

    fn add_limb(a in strategy::natural(), b in strategy::limb()) {
        ours: all_forms_prim!(a.to_natural(), +, +=, b),
        gmp: a.to_rug() + b,
    }

    fn checked_sub(a in strategy::natural(), b in strategy::natural()) {
        ours: {
            let (x, y) = (a.to_natural(), b.to_natural());
            let r = x.checked_sub(&y);
            let mut z = x.clone();
            let done = z.checked_sub_assign(&y);
            prop_assert_eq!(done, r.is_some());
            prop_assert_eq!(&z, r.as_ref().unwrap_or(&x), "failed checked_sub_assign must not modify");
            r
        },
        gmp: {
            let d = a.to_rug() - b.to_rug();
            (d >= 0).then_some(d)
        },
    }

    /// Subtraction of a smaller value: the full operator surface.
    fn sub((a, b) in (strategy::natural(), strategy::natural())
        .prop_map(|(a, b)| if a.to_rug() >= b.to_rug() { (a, b) } else { (b, a) }))
    {
        ours: all_forms!(a.to_natural(), -, -=, b.to_natural()),
        gmp: a.to_rug() - b.to_rug(),
    }

    fn sub_unbalanced((a, b) in strategy::unbalanced_pair(DEFAULT_MAX_LIMBS)) {
        ours: {
            let (x, y) = (a.to_natural(), b.to_natural());
            let (x, y) = if x >= y { (x, y) } else { (y, x) };
            all_forms!(x, -, -=, y)
        },
        gmp: (a.to_rug() - b.to_rug()).abs(),
    }

    fn sub_around_inline(a in around_inline(), b in around_inline()) {
        ours: a.to_natural().checked_sub(&b.to_natural()),
        gmp: {
            let d = a.to_rug() - b.to_rug();
            (d >= 0).then_some(d)
        },
    }

    fn sub_limb(a in strategy::natural(), b in strategy::limb()) {
        ours: {
            let x = a.to_natural();
            let r = x.checked_sub_limb(b);
            if r.is_some() {
                prop_assert_eq!(&r, &Some(all_forms_prim!(x, -, -=, b)));
            }
            r
        },
        gmp: {
            let d = a.to_rug() - b;
            (d >= 0).then_some(d)
        },
    }

    fn shl(a in strategy::natural(), k in shift()) {
        ours: all_forms_prim!(a.to_natural(), <<, <<=, k),
        gmp: a.to_rug() << u32::try_from(k).unwrap(),
    }

    fn shr(a in strategy::natural(), k in shift()) {
        ours: all_forms_prim!(a.to_natural(), >>, >>=, k),
        gmp: a.to_rug() >> u32::try_from(k).unwrap(),
    }

    fn shr_huge(a in strategy::natural(), k in u64::from(u32::MAX)..) {
        ours: all_forms_prim!(a.to_natural(), >>, >>=, k),
        gmp: Integer::new(),
    }

    fn and(a in strategy::natural(), b in strategy::natural()) {
        ours: all_forms!(a.to_natural(), &, &=, b.to_natural()),
        gmp: a.to_rug() & b.to_rug(),
    }

    fn or(a in strategy::natural(), b in strategy::natural()) {
        ours: all_forms!(a.to_natural(), |, |=, b.to_natural()),
        gmp: a.to_rug() | b.to_rug(),
    }

    fn xor(a in strategy::natural(), b in strategy::natural()) {
        ours: all_forms!(a.to_natural(), ^, ^=, b.to_natural()),
        gmp: a.to_rug() ^ b.to_rug(),
    }

    /// Equal-length operands, often sharing high limbs so xor cancels them.
    fn xor_same_len((a, b) in same_len_sharing_top()) {
        ours: all_forms!(a.to_natural(), ^, ^=, b.to_natural()),
        gmp: a.to_rug() ^ b.to_rug(),
    }

    fn bit_length(a in strategy::natural()) {
        ours: a.to_natural().bit_length(),
        gmp: a.to_rug().significant_bits(),
    }

    fn trailing_zeros(a in strategy::natural()) {
        ours: a.to_natural().trailing_zeros(),
        gmp: a.to_rug().find_one(0).map(Integer::from),
    }

    fn count_ones(a in strategy::natural()) {
        ours: a.to_natural().count_ones(),
        gmp: a.to_rug().count_ones().unwrap(),
    }

    fn bit(a in strategy::natural(), k in shift()) {
        ours: a.to_natural().bit(k),
        gmp: a.to_rug().get_bit(u32::try_from(k).unwrap()),
    }

    fn is_power_of_two(a in strategy::natural()) {
        ours: a.to_natural().is_power_of_two(),
        gmp: a.to_rug().is_power_of_two(),
    }

    fn mul(a in strategy::natural(), b in strategy::natural()) {
        ours: all_forms!(a.to_natural(), *, *=, b.to_natural()),
        gmp: a.to_rug() * b.to_rug(),
    }

    fn mul_unbalanced((a, b) in strategy::unbalanced_pair(DEFAULT_MAX_LIMBS)) {
        ours: all_forms!(a.to_natural(), *, *=, b.to_natural()),
        gmp: a.to_rug() * b.to_rug(),
    }

    fn mul_around_inline(a in around_inline(), b in around_inline()) {
        ours: all_forms!(a.to_natural(), *, *=, b.to_natural()),
        gmp: a.to_rug() * b.to_rug(),
    }

    /// Sizes past every algorithm threshold (the kernel tests cover the
    /// crossovers themselves).
    fn mul_large(a in natural_beyond_thresholds(), b in natural_beyond_thresholds()) {
        ours: a.to_natural() * b.to_natural(),
        gmp: a.to_rug() * b.to_rug(),
    }

    fn mul_limb(a in strategy::natural(), b in strategy::limb()) {
        ours: all_forms_prim!(a.to_natural(), *, *=, b),
        gmp: a.to_rug() * b,
    }

    fn square(a in prop_oneof![strategy::natural(), natural_beyond_thresholds()]) {
        ours: {
            let x = a.to_natural();
            let s = x.square();
            // `&x * &x` takes the squaring path too.
            prop_assert_eq!(&(&x * &x), &s);
            s
        },
        gmp: a.to_rug().square(),
    }
}

/// A natural of up to three times the largest multiplication threshold.
fn natural_beyond_thresholds() -> BoxedStrategy<Mag> {
    let max =
        3 * crate::kernel::mul::MUL_TOOM3_THRESHOLD.max(crate::kernel::mul::SQR_TOOM3_THRESHOLD);
    strategy::natural_with_len(1..=max).boxed()
}

proptest! {
    #[test]
    fn eq_and_hash_ignore_representation(a in strategy::natural_up_to(INLINE)) {
        // Push the value through heap storage and back down: it stays on the
        // heap after shrinking, and must still equal and hash like the inline one.
        let inline = a.to_natural();
        let big = Natural::power_of_two(64 * 8);
        let mut heap = &inline + &big;
        heap -= &big;
        prop_assert!(!heap.limbs.is_inline());
        prop_assert_eq!(&heap, &inline);
        prop_assert_eq!(hash_of(&heap), hash_of(&inline));
        prop_assert_eq!(heap.cmp(&inline), std::cmp::Ordering::Equal);
        heap.shrink_to_fit();
        prop_assert!(heap.limbs.is_inline());
        prop_assert_eq!(&heap.clone(), &inline);
    }

    #[test]
    fn clone_moves_small_values_inline(a in strategy::natural_up_to(INLINE)) {
        let big = Natural::power_of_two(64 * 8);
        let heap = &a.to_natural() + &big - &big;
        prop_assert!(heap.clone().limbs.is_inline());
    }

    #[test]
    fn add_sub_round_trip(a in strategy::natural(), b in strategy::natural()) {
        let (x, y) = (a.to_natural(), b.to_natural());
        prop_assert_eq!(&(&x + &y) - &y, x);
    }

    #[test]
    fn mul_is_associative_and_distributive(
        a in strategy::natural(),
        b in strategy::natural(),
        c in strategy::natural(),
    ) {
        let (x, y, z) = (a.to_natural(), b.to_natural(), c.to_natural());
        prop_assert_eq!(&(&x * &y) * &z, &x * &(&y * &z));
        prop_assert_eq!(&x * &(&y + &z), &(&x * &y) + &(&x * &z));
    }

    #[test]
    fn shl_shr_round_trip(a in strategy::natural(), k in shift()) {
        let x = a.to_natural();
        prop_assert_eq!(&(&x << k) >> k, x);
    }
}

#[test]
fn carries_across_limb_boundaries() {
    for n in 1..=8u64 {
        let ones = Natural::power_of_two(64 * n) - 1;
        assert_eq!(ones.limb_len() as u64, n);
        assert!(ones.limbs().iter().all(|&x| x == Limb::MAX));
        // Carry ripples out of every limb into a new one.
        assert_eq!(&ones + 1, Natural::power_of_two(64 * n));
        assert_eq!(&ones + &Natural::one(), Natural::power_of_two(64 * n));
        assert_eq!(&ones + &ones, Natural::power_of_two(64 * n + 1) - 2);
        // Borrow ripples through every limb and the top limb vanishes.
        assert_eq!(Natural::power_of_two(64 * n) - &ones, Natural::one());
        assert_eq!((Natural::power_of_two(64 * n) - 1).limb_len() as u64, n);
        // Shifts across the boundary.
        assert_eq!(&ones << 1, Natural::power_of_two(64 * n + 1) - 2);
        assert_eq!((&ones << 64).limb_len() as u64, n + 1);
        assert_eq!(&ones >> (64 * n - 1), Natural::one());
    }
}

#[test]
fn zero_and_one() {
    assert!(Natural::ZERO.is_zero() && Natural::new().limbs().is_empty());
    assert!(Natural::one().is_one());
    assert_eq!(Natural::default(), Natural::ZERO);
    assert_eq!(Natural::ZERO.bit_length(), 0);
    assert_eq!(Natural::ZERO.trailing_zeros(), None);
    assert_eq!(&Natural::ZERO << 1000, Natural::ZERO);
    assert_eq!(Natural::ZERO - 0, Natural::ZERO);
    assert_eq!(Natural::ZERO + 0, Natural::ZERO);
}

#[test]
#[should_panic(expected = "attempt to subtract with overflow")]
fn sub_underflow_panics() {
    let _ = Natural::from(1u64) - &Natural::from(2u64);
}

#[test]
#[should_panic(expected = "attempt to subtract with overflow")]
fn sub_limb_underflow_panics() {
    let _ = Natural::ZERO - 1;
}

#[test]
fn hex_formatting() {
    let x = Natural::from(0x1_0000_0000_0000_00ab_u128);
    assert_eq!(format!("{x:x}"), "100000000000000ab");
    assert_eq!(format!("{x:#X}"), "0x100000000000000AB");
    assert_eq!(format!("{:?}", Natural::ZERO), "0x0");
    assert_eq!(format!("{:>6x}", Natural::from(255u64)), "    ff");
}

#[test]
fn natural_is_small() {
    assert_eq!(size_of::<Natural>(), 32);
}
