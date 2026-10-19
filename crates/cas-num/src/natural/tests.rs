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
use crate::testkit::{Mag, ToRug, all_forms, all_forms_prim, differential};

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

    fn div(a in strategy::natural(), d in divisor()) {
        ours: all_forms!(a.to_natural(), /, /=, d.to_natural()),
        gmp: a.to_rug() / d.to_rug(),
    }

    fn rem(a in strategy::natural(), d in divisor()) {
        ours: all_forms!(a.to_natural(), %, %=, d.to_natural()),
        gmp: a.to_rug() % d.to_rug(),
    }

    fn div_rem_large(a in natural_beyond_thresholds(), d in natural_beyond_thresholds()) {
        ours: a.to_natural().div_rem(&d.to_natural()),
        gmp: a.to_rug().div_rem(d.to_rug()),
    }

    fn checked_div_rem(a in strategy::natural(), d in strategy::natural()) {
        ours: a.to_natural().checked_div_rem(&d.to_natural()),
        gmp: (!d.is_zero()).then(|| a.to_rug().div_rem(d.to_rug())),
    }

    fn div_limb(a in strategy::natural(), d in strategy::top_limb()) {
        ours: all_forms_prim!(a.to_natural(), /, /=, d),
        gmp: a.to_rug() / d,
    }

    fn rem_limb(a in strategy::natural(), d in strategy::top_limb()) {
        ours: {
            let x = a.to_natural();
            let r = &x % d;
            prop_assert_eq!(x.clone() % d, r);
            prop_assert_eq!(x.div_rem_limb(d).1, r);
            let mut y = x;
            y %= d;
            prop_assert_eq!(y, Natural::from(r));
            r
        },
        gmp: a.to_rug() % d,
    }

    fn div_exact(q in strategy::natural(), d in divisor()) {
        ours: (q.to_natural() * d.to_natural()).div_exact(&d.to_natural()),
        gmp: q.to_rug(),
    }

    /// Even divisors too: both sides lose the divisor's factors of two.
    fn div_exact_even(q in strategy::natural(), d in divisor(), twos in 0..130u64) {
        ours: {
            let d = d.to_natural() << twos;
            (q.to_natural() * &d).div_exact(&d)
        },
        gmp: q.to_rug(),
    }

    fn div_exact_limb(q in strategy::natural(), d in strategy::top_limb()) {
        ours: (q.to_natural() * d).div_exact_limb(d),
        gmp: q.to_rug(),
    }

    fn is_divisible_by(a in strategy::natural(), d in strategy::natural()) {
        ours: a.to_natural().is_divisible_by(&d.to_natural()),
        gmp: a.to_rug().is_divisible(&d.to_rug()),
    }

    /// Multiples, and multiples plus one: both answers occur often.
    fn is_divisible_by_shaped(q in strategy::natural(), d in divisor(), plus_one in any::<bool>()) {
        ours: {
            let a = q.to_natural() * d.to_natural() + Limb::from(plus_one);
            let d = d.to_natural();
            if let [x] = d.limbs() {
                prop_assert_eq!(a.is_divisible_by_limb(*x), a.is_divisible_by(&d));
            }
            a.is_divisible_by(&d)
        },
        gmp: (q.to_rug() * d.to_rug() + u32::from(plus_one)).is_divisible(&d.to_rug()),
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

/// A non-zero natural, sometimes with few limbs relative to `natural()`.
fn divisor() -> BoxedStrategy<Mag> {
    prop_oneof![
        strategy::natural().prop_filter("non-zero", |d| !d.is_zero()),
        strategy::natural_with_len(1..=3usize),
    ]
    .boxed()
}

/// A natural of up to three times the largest multiplication threshold.
fn natural_beyond_thresholds() -> BoxedStrategy<Mag> {
    let max = 3 * crate::testkit::tested(
        crate::kernel::mul::MUL_TOOM3_THRESHOLD.max(crate::kernel::mul::SQR_TOOM3_THRESHOLD),
    );
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
    fn div_rem_identity(a in strategy::natural(), d in divisor()) {
        let (x, y) = (a.to_natural(), d.to_natural());
        let (q, r) = x.div_rem(&y);
        prop_assert!(r < y);
        prop_assert_eq!(&q * &y + &r, x);
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
    assert_eq!(format!("{:?}", Natural::ZERO), "0");
    assert_eq!(format!("{:>6x}", Natural::from(255u64)), "    ff");
}

/// Four inline limbs and a length (CAS-40: 4 limbs beat 3 on 2-4 limb values
/// for 8 more bytes); a big `Integer` stays one pointer and a tag.
#[test]
fn natural_is_small() {
    assert_eq!(size_of::<Natural>(), 40);
    assert_eq!(size_of::<crate::Integer>(), 16);
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn div_by_zero_panics() {
    let _ = Natural::one() / Natural::ZERO;
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn rem_by_zero_limb_panics() {
    let _ = Natural::one() % 0;
}

/// `div_exact` where Hensel division works in blocks, against GMP: quotients
/// as long as `an - dn + 1` limbs (the dividend's top limbs reach the
/// divisor) and one limb shorter, longer and shorter than the divisor.
/// Skipped when the configuration starts the blocks above 20000 limbs.
#[test]
fn div_exact_in_blocks() {
    let t = crate::kernel::div::bdiv::BDIV_MU_THRESHOLD;
    if t > 20_000 {
        return;
    }
    let mut state = 0x1357_9bdf_0246_8ace_u64;
    let mut random = |n: usize, top: Limb| -> Vec<Limb> {
        let mut v: Vec<Limb> = (0..n)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                state.wrapping_mul(0x2545_f491_4f6c_dd1d)
            })
            .collect();
        v[n - 1] = top;
        v
    };
    for (qn, dn) in [(t + 7, t + 7), (2 * t, t + 1), (t / 2 + 1, t)] {
        for (q_top, d_top) in [(Limb::MAX, 1), (1, Limb::MAX)] {
            let (q, d) = (random(qn, q_top), random(dn, d_top));
            let (q, d) = (Natural::from_limbs(&q), Natural::from_limbs(&d));
            let a = &q * &d;
            assert_eq!(a.div_exact(&d).to_rug(), q.to_rug(), "qn = {qn}, dn = {dn}");
        }
    }
}

/// Debug builds check `div_exact`'s contract.
#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "division is not exact")]
fn div_exact_panics_when_inexact() {
    let _ = Natural::from(u128::MAX).div_exact(&Natural::from(1u128 << 64));
}

#[test]
#[cfg(debug_assertions)]
#[should_panic(expected = "division is not exact")]
fn div_exact_limb_panics_when_inexact() {
    let _ = Natural::from(7u64).div_exact_limb(2);
}

/// Release builds trust it, as GMP's `mpz_divexact` does: an inexact
/// division returns an unspecified value instead of panicking, for divisors
/// with more trailing zeros than the dividend, longer than it, or odd.
#[test]
#[cfg(not(debug_assertions))]
fn div_exact_inexact_returns() {
    let a = Natural::from(u128::MAX);
    for d in [
        Natural::from(1u128 << 64),
        Natural::from(u128::MAX) * 3u64,
        Natural::from(3u128 << 64 | 5),
    ] {
        let _ = a.div_exact(&d);
    }
    let _ = Natural::from(7u64).div_exact_limb(2);
}

#[test]
fn zero_divides_only_zero() {
    assert!(Natural::ZERO.is_divisible_by(&Natural::ZERO));
    assert!(!Natural::one().is_divisible_by(&Natural::ZERO));
    assert!(Natural::ZERO.is_divisible_by_limb(0));
    assert!(!Natural::one().is_divisible_by_limb(0));
}
