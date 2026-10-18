//! MPFR oracle for `Float` and ball arithmetic (CAS-DOC-5 §11).
//!
//! - [`Dyadic`]: `±m·2^e` as a neutral exchange value, the binary-float
//!   counterpart of [`Signed`]. The exponent is a [`Signed`] too, so values
//!   beyond MPFR's exponent range are representable.
//! - [`bracket`]: encloses a true value between two MPFR results computed
//!   `GUARD` bits above the precision under test, rounded toward −∞ and +∞.
//! - [`check_encloses`]: a ball, given as an exact midpoint and radius,
//!   against a bracket.
//! - [`dyadic`] and friends: proptest strategies with adversarial shapes.
//!
//! Rounding names differ: rug's `Round::Down`/`Round::Up` are toward −∞/+∞
//! (our `Floor`/`Ceil`), and rug's `Round::Zero`/`Round::AwayZero` are our
//! `Down`/`Up`.
//!
//! MPFR exponents are `i32` in rug (by default about `±2^30`). Tests of
//! exponents near and beyond `±2^63` therefore compare against MPFR on the
//! operands with a common power of two removed: `x·2^E op y·2^E` against
//! `(x op y)·2^E`, which every operation of interest commutes with.

use std::cmp::Ordering;
use std::fmt;

use proptest::prelude::*;
use proptest::test_runner::TestCaseError;
use rug::float::Round;
use rug::{Float, Integer, Rational};

use super::model::{Mag, Signed};
use super::{ToRug, strategy};

/// Extra bits the oracle computes above the precision under test.
pub(crate) const GUARD: u32 = 64;

/// A binary float `±mant·2^exp` with an odd mantissa, or zero.
///
/// Normalised as stated (zero has a zero exponent), so derived equality is
/// value equality.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub(crate) struct Dyadic {
    mant: Signed,
    exp: Signed,
}

impl Dyadic {
    /// `mant·2^exp`, normalised (trailing zero bits move into the exponent).
    pub(crate) fn new(mant: &Integer, exp: &Integer) -> Self {
        if *mant == 0 {
            return Dyadic::default();
        }
        let tz = mant.find_one(0).expect("nonzero");
        Dyadic {
            mant: Signed::from_rug(&(mant.clone() >> tz)),
            exp: Signed::from_rug(&(exp.clone() + tz)),
        }
    }

    pub(crate) fn mant(&self) -> &Signed {
        &self.mant
    }

    pub(crate) fn exp(&self) -> &Signed {
        &self.exp
    }

    /// The mantissa's magnitude, for tests that compare limbs.
    pub(crate) fn mant_mag(&self) -> &Mag {
        self.mant.mag()
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.mant.mag().is_zero()
    }

    /// The exact value of a finite MPFR float; `None` for NaN and infinities.
    pub(crate) fn from_rug(x: &Float) -> Option<Self> {
        let (m, e) = x.to_integer_exp()?;
        Some(Self::new(&m, &Integer::from(e)))
    }

    /// The same value times `2^k`.
    pub(crate) fn mul_2exp(&self, k: &Integer) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        Dyadic {
            mant: self.mant.clone(),
            exp: Signed::from_rug(&(self.exp.to_rug() + k)),
        }
    }

    /// Bit length of the mantissa (0 for zero).
    pub(crate) fn mant_bits(&self) -> u32 {
        self.mant.to_rug().significant_bits()
    }

    /// The exact MPFR float, when the exponent is inside MPFR's range.
    pub(crate) fn to_rug(&self) -> Option<Float> {
        let prec = self.mant_bits().max(1);
        let mut f = Float::with_val(prec, self.mant.to_rug());
        let e = self.exp.to_rug().to_i32()?;
        let (lo, hi) = rug::float::allowed_exp_range();
        // `f` has exponent `prec`; the result has exponent `prec + e`.
        let top = i32::try_from(prec).ok()?.checked_add(e)?;
        if !self.is_zero() && !(lo..=hi).contains(&top) {
            return None;
        }
        f <<= e;
        Some(f)
    }

    /// The exact rational value. Panics past a few million bits: callers keep
    /// exponents small when they need this.
    pub(crate) fn to_rational(&self) -> Rational {
        let e = self
            .exp
            .to_rug()
            .to_i32()
            .expect("exponent too large for an exact rational");
        assert!(
            e.unsigned_abs() < 1 << 24,
            "exponent too large for an exact rational: {e}"
        );
        let m = Rational::from(self.mant.to_rug());
        if e >= 0 {
            m << e
        } else {
            m >> e.unsigned_abs()
        }
    }
}

/// `±0x…p±e`: the mantissa in limb-grouped hex (see [`Mag`]'s `Debug`).
impl fmt::Debug for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} * 2^{}", self.mant, self.exp.to_rug())
    }
}

/// An exact interval `[lo, hi]` known to contain a true value.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Bracket {
    pub(crate) lo: Rational,
    pub(crate) hi: Rational,
}

impl Bracket {
    pub(crate) fn width(&self) -> Rational {
        Rational::from(&self.hi - &self.lo)
    }

    pub(crate) fn contains(&self, x: &Rational) -> bool {
        self.lo <= *x && *x <= self.hi
    }
}

/// Encloses `f`'s true value: `f(prec + GUARD, Down)` and `f(prec + GUARD, Up)`.
///
/// `f` computes the value in MPFR at the given precision and rounding, as one
/// correctly rounded operation or as a chain of directed ones (a chain must
/// itself round every step in the given direction for the result to bound
/// the value). Panics if MPFR returns a non-finite value.
pub(crate) fn bracket(prec: u32, f: impl Fn(u32, Round) -> Float) -> Bracket {
    let p = prec + GUARD;
    let exact = |r: Round| {
        let x = f(p, r);
        x.to_rational()
            .unwrap_or_else(|| panic!("oracle returned non-finite {x}"))
    };
    let (lo, hi) = (exact(Round::Down), exact(Round::Up));
    assert!(lo <= hi, "oracle bracket reversed: [{lo}, {hi}]");
    Bracket { lo, hi }
}

/// Fails unless the ball `mid ± rad` is consistent with the true value in `b`.
///
/// When `rad ≥ hi − lo`, the ball must contain all of `[lo, hi]`: the oracle is
/// then sharper than the ball, so the true value could be anywhere in it.
/// Otherwise the ball must at least meet `[lo, hi]` (the true value lies in
/// both).
pub(crate) fn check_encloses(
    mid: &Rational,
    rad: &Rational,
    b: &Bracket,
) -> Result<(), TestCaseError> {
    assert!(*rad >= 0, "negative radius {rad}");
    let (ball_lo, ball_hi) = (Rational::from(mid - rad), Rational::from(mid + rad));
    let ok = if *rad >= b.width() {
        ball_lo <= b.lo && b.hi <= ball_hi
    } else {
        ball_lo <= b.hi && b.lo <= ball_hi
    };
    if ok {
        Ok(())
    } else {
        Err(TestCaseError::fail(format!(
            "ball does not enclose the oracle bracket\n  ball:   {mid} ± {rad}\n  \
             as [lo, hi]: [{ball_lo}, {ball_hi}]\n  oracle: [{}, {}]",
            b.lo, b.hi
        )))
    }
}

/// Fails unless `actual` (ours, as a `Dyadic`) equals MPFR's `expected` and
/// the inexact flags agree in sign (`Ordering` as MPFR's ternary value).
pub(crate) fn check_rounded(
    actual: (&Dyadic, Ordering),
    expected: (&Float, Ordering),
) -> Result<(), TestCaseError> {
    let want = Dyadic::from_rug(expected.0).expect("oracle returned non-finite value");
    if *actual.0 == want && actual.1 == expected.1 {
        Ok(())
    } else {
        Err(TestCaseError::fail(format!(
            "differs from MPFR\n  ours:     {:?} ({:?})\n  expected: {:?} ({:?})",
            actual.0, actual.1, want, expected.1
        )))
    }
}

/// A float with a mantissa of at most `max_limbs` limbs and an exponent
/// in `exp_range`, with adversarial mantissa shapes (all-ones, powers of two,
/// `2^k ± 1`) from [`strategy::signed_up_to`].
pub(crate) fn dyadic_with(
    max_limbs: usize,
    exp_range: std::ops::RangeInclusive<i64>,
) -> BoxedStrategy<Dyadic> {
    (strategy::signed_up_to(max_limbs), exp_range)
        .prop_map(|(m, e)| Dyadic::new(&m.to_rug(), &Integer::from(e)))
        .boxed()
}

/// [`dyadic_with`] at up to 4 limbs (the inline `LimbVec` range) with
/// exponents mostly near 0, sometimes at the edge of MPFR's range.
pub(crate) fn dyadic() -> BoxedStrategy<Dyadic> {
    prop_oneof![
        6 => dyadic_with(4, -300..=300),
        1 => dyadic_with(12, -4000..=4000),
        1 => dyadic_with(2, -(1 << 29)..=(1 << 29)),
    ]
    .boxed()
}

/// Two floats whose exponents differ by up to `max_gap` bits, so that one
/// operand falls partly or entirely below the other's rounding position.
pub(crate) fn far_apart_pair(max_gap: i64) -> BoxedStrategy<(Dyadic, Dyadic)> {
    (
        dyadic_with(4, -64..=64),
        dyadic_with(4, -64..=64),
        0..=max_gap,
    )
        .prop_map(|(a, b, gap)| (a, b.mul_2exp(&Integer::from(-gap))))
        .boxed()
}

/// An exponent shift near or beyond `±2^63`, for the scaling tests described
/// in the module docs.
pub(crate) fn huge_shift() -> BoxedStrategy<Integer> {
    let edge = |k: u32| Integer::from(1) << k;
    prop_oneof![
        (-3i64..=3).prop_map(move |d| edge(63) + d),
        (-3i64..=3).prop_map(move |d| -edge(63) + d),
        (64u32..200, any::<bool>()).prop_map(move |(k, neg)| if neg { -edge(k) } else { edge(k) }),
    ]
    .boxed()
}

#[cfg(test)]
mod tests {
    use super::*;

    proptest! {
        #[test]
        fn rug_round_trip(a in dyadic()) {
            let f = a.to_rug().expect("strategy stays in MPFR's range");
            prop_assert_eq!(Dyadic::from_rug(&f), Some(a));
        }

        #[test]
        fn rational_matches_rug(a in dyadic_with(4, -300..=300)) {
            let f = a.to_rug().unwrap();
            prop_assert_eq!(a.to_rational(), f.to_rational().unwrap());
        }

        #[test]
        fn mant_mag_is_abs_mantissa(a in dyadic()) {
            prop_assert_eq!(a.mant_mag(), &Mag::from_rug(&a.mant().to_rug().abs()));
        }

        /// The second operand sits `0..=gap` bits lower, and the shift is exact.
        #[test]
        fn far_apart_shifts_exactly((a, b) in far_apart_pair(1 << 20)) {
            prop_assert!(a.to_rug().is_some());
            if !b.is_zero() {
                prop_assert!(b.exp().to_rug() <= 64);
                prop_assert!(b.exp().to_rug() >= -64 - (1 << 20) - 256);
            }
        }

        #[test]
        fn huge_shift_leaves_i64_or_borders_it(k in huge_shift()) {
            let edge = Integer::from(1) << 63u32;
            prop_assert!(k.clone().abs() >= edge - 3u32);
        }

        #[test]
        fn mantissa_is_odd(a in dyadic()) {
            prop_assert!(a.is_zero() || a.mant().to_rug().is_odd());
            prop_assert!(!a.is_zero() || a.exp().to_rug() == 0);
        }

        /// The bracket of `a / b` contains the exact quotient and is at most
        /// 2 ulps wide at `prec + GUARD` bits.
        #[test]
        fn bracket_contains_quotient(
            a in dyadic_with(4, -300..=300),
            b in dyadic_with(4, -300..=300),
            prec in 1u32..300,
        ) {
            prop_assume!(!b.is_zero());
            let (fa, fb) = (a.to_rug().unwrap(), b.to_rug().unwrap());
            let br = bracket(prec, |p, r| Float::with_val_round(p, &fa / &fb, r).0);
            let q = a.to_rational() / b.to_rational();
            prop_assert!(br.contains(&q));
            if q != 0 {
                let ulp = q.clone().abs() >> (prec + GUARD - 1);
                prop_assert!(br.width() <= ulp);
            }
        }
    }

    /// `sqrt(2)`: brackets nest across precisions and square around 2.
    #[test]
    fn bracket_of_sqrt2() {
        let mut prev: Option<Bracket> = None;
        for prec in [1, 53, 64, 200, 1000] {
            let b = bracket(prec, |p, r| {
                Float::with_val_round(p, Float::with_val(2, 2).sqrt_ref(), r).0
            });
            assert!(
                Rational::from(b.lo.square_ref()) <= 2 && Rational::from(b.hi.square_ref()) >= 2
            );
            if let Some(p) = &prev {
                assert!(p.lo <= b.lo && b.hi <= p.hi, "brackets must nest");
            }
            prev = Some(b);
        }
    }

    #[test]
    fn encloses_accepts_and_rejects() {
        let b = Bracket {
            lo: Rational::from((99, 100)),
            hi: Rational::from((101, 100)),
        };
        let one = Rational::from(1);
        // Wide ball containing the bracket.
        assert!(check_encloses(&one, &Rational::from((1, 10)), &b).is_ok());
        // Wide ball missing part of it.
        assert!(check_encloses(&Rational::from((11, 10)), &Rational::from((1, 10)), &b).is_err());
        // Narrow ball inside the bracket: consistent.
        assert!(check_encloses(&one, &Rational::from((1, 1000)), &b).is_ok());
        // Narrow ball outside the bracket.
        assert!(check_encloses(&Rational::from(2), &Rational::from((1, 1000)), &b).is_err());
        // Exact ball at an endpoint of a zero-width bracket.
        let exact = Bracket {
            lo: one.clone(),
            hi: one.clone(),
        };
        assert!(check_encloses(&one, &Rational::new(), &exact).is_ok());
        assert!(check_encloses(&Rational::from(2), &Rational::new(), &exact).is_err());
    }

    #[test]
    fn rounded_compares_value_and_ternary() {
        // A borrowed operand keeps the division incomplete, so it rounds once
        // at 10 bits (an owned one would divide at its own 2 bits first).
        let one = Float::with_val(2, 1);
        let (f, o) = Float::with_val_round(10, &one / 3u32, Round::Nearest);
        assert_ne!(o, Ordering::Equal);
        let d = Dyadic::from_rug(&f).unwrap();
        assert!(check_rounded((&d, o), (&f, o)).is_ok());
        assert!(check_rounded((&d, o.reverse()), (&f, o)).is_err());
        assert!(check_rounded((&Dyadic::default(), o), (&f, o)).is_err());
    }

    #[test]
    fn huge_shift_round_trip() {
        let a = Dyadic::new(&Integer::from(3), &Integer::from(-1));
        let k = Integer::from(1) << 70u32;
        let big = a.mul_2exp(&k);
        assert_eq!(big.to_rug(), None, "beyond MPFR's range");
        assert_eq!(big.mul_2exp(&Integer::from(-&k)), a);
    }
}
