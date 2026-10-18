//! `Mag` against exact rational arithmetic (and MPFR brackets for exp, expm1,
//! log1p): every upper bound must be `>=` the exact value and within a factor
//! `1 + 2^-28` of it, every lower bound `<=` and within `1 - 2^-28`.
#![allow(
    clippy::many_single_char_names,
    clippy::needless_pass_by_value,
    reason = "operands a, b, x and their values p, q; results passed as built"
)]

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use super::{HIGH, LOW, Mag};
use crate::testkit::float::{bracket, huge_shift};
use crate::testkit::{Signed, ToRug};
use crate::{Integer, Rational};

/// `2^-28`, the tightness of one operation.
fn tight() -> Rational {
    Rational::new(Integer::ONE, Integer::from(1i64 << 28))
}

fn value(m: &Mag) -> Option<Rational> {
    m.to_float().map(|f| f.to_rational())
}

/// Finite values: edge and random mantissas, exponents mostly small.
fn finite(exps: std::ops::RangeInclusive<i64>) -> BoxedStrategy<Mag> {
    let man = prop_oneof![
        1 => Just(LOW),
        1 => Just(HIGH - 1),
        1 => Just(LOW + 1),
        5 => LOW..HIGH,
    ];
    (man, exps)
        .prop_map(|(m, e)| Mag::new(m, Integer::from(e)))
        .boxed()
}

/// Finite values plus zero and `+∞`.
fn mag() -> BoxedStrategy<Mag> {
    prop_oneof![
        1 => Just(Mag::ZERO),
        1 => Just(Mag::INF),
        12 => finite(-300..=300),
    ]
    .boxed()
}

/// Fails unless `got` is an upper bound of `exact` within `1 + slack × 2^-28`.
fn check_upper(what: &str, got: &Mag, exact: &Rational, slack: i64) -> Result<(), TestCaseError> {
    let v = value(got)
        .ok_or_else(|| TestCaseError::fail(format!("{what}: +inf for finite {exact}")))?;
    prop_assert!(v >= *exact, "{}: {:?} below the exact {}", what, got, exact);
    let limit = exact * &(&Rational::ONE + &(&tight() * &Rational::from(Integer::from(slack))));
    prop_assert!(
        v <= limit,
        "{}: {:?} too far above the exact {}",
        what,
        got,
        exact
    );
    Ok(())
}

/// Fails unless `got` is a lower bound of `exact` within `1 - slack × 2^-28`.
fn check_lower(what: &str, got: &Mag, exact: &Rational, slack: i64) -> Result<(), TestCaseError> {
    let v =
        value(got).ok_or_else(|| TestCaseError::fail(format!("{what}: +inf as a lower bound")))?;
    prop_assert!(v <= *exact, "{}: {:?} above the exact {}", what, got, exact);
    let limit = exact * &(&Rational::ONE - &(&tight() * &Rational::from(Integer::from(slack))));
    prop_assert!(
        v >= limit,
        "{}: {:?} too far below the exact {}",
        what,
        got,
        exact
    );
    Ok(())
}

/// MPFR's value of `x` at 200 bits rounded `rnd`, as a rational.
fn mpfr_bracket(x: &Mag, f: impl Fn(&MpFloat, u32, MpRound) -> MpFloat) -> (Rational, Rational) {
    let xf = x.to_float().unwrap();
    let xr = xf.to_rational().to_rug();
    let mp = MpFloat::with_val(400, &xr);
    let b = bracket(200, |p, r| f(&mp, p, r));
    let ours = |q: &rug::Rational| {
        let (n, d) = q.clone().into_numer_denom();
        Rational::new(
            Signed::from_rug(&n).to_integer(),
            Signed::from_rug(&d).to_integer(),
        )
    };
    (ours(&b.lo), ours(&b.hi))
}

proptest! {
    /// Zero and `+∞` mixed with finite values: never a panic, and the
    /// results follow the extended rules (0 × ∞ is 0 for an exact zero).
    #[test]
    fn specials(a in mag(), b in mag()) {
        for r in [a.add(&b), a.add_lower(&b), a.mul(&b), a.mul_lower(&b), a.div(&b), a.div_lower(&b), a.sub(&b), a.sub_lower(&b)] {
            prop_assert!(r.is_zero() || r.is_inf() || (LOW..HIGH).contains(&r.man));
        }
        prop_assert!(a.add(&b) >= a.add_lower(&b) && a.mul(&b) >= a.mul_lower(&b));
        prop_assert!(a.div(&b) >= a.div_lower(&b) && a.sub(&b) >= a.sub_lower(&b));
        prop_assert_eq!(a.add(&b), b.add(&a));
        prop_assert_eq!(a.mul(&b), b.mul(&a));
    }

    #[test]
    fn arithmetic_bounds(a in finite(-300..=300), b in finite(-300..=300)) {
        let (p, q) = (value(&a).unwrap(), value(&b).unwrap());
        check_upper("add", &a.add(&b), &(&p + &q), 1)?;
        check_lower("add_lower", &a.add_lower(&b), &(&p + &q), 1)?;
        check_upper("mul", &a.mul(&b), &(&p * &q), 1)?;
        check_lower("mul_lower", &a.mul_lower(&b), &(&p * &q), 1)?;
        check_upper("div", &a.div(&b), &(&p / &q), 1)?;
        check_lower("div_lower", &a.div_lower(&b), &(&p / &q), 1)?;
        let diff = &p - &q;
        if diff.is_positive() {
            check_upper("sub", &a.sub(&b), &diff, 1)?;
            check_lower("sub_lower", &a.sub_lower(&b), &diff, 1)?;
        } else {
            prop_assert!(a.sub_lower(&b).is_zero());
            prop_assert!(value(&a.sub(&b)).unwrap() <= Rational::ZERO + &(&q * &tight()), "sub of a smaller value");
        }
        prop_assert_eq!(a.cmp(&b), p.cmp(&q));
    }

    /// Operands close in size, where sums carry and differences cancel.
    #[test]
    fn close_operands(m in LOW..HIGH, d in 0u32..2000, e in -10i64..10, gap in 0i64..70) {
        let a = Mag::new(m, Integer::from(e));
        let n = (m + d).min(HIGH - 1);
        let b = Mag::new(n, Integer::from(e - gap));
        let (p, q) = (value(&a).unwrap(), value(&b).unwrap());
        check_upper("add", &a.add(&b), &(&p + &q), 1)?;
        let diff = &p - &q;
        if diff.is_positive() {
            check_upper("sub", &a.sub(&b), &diff, 1)?;
            check_lower("sub_lower", &a.sub_lower(&b), &diff, 1)?;
        }
    }

    #[test]
    fn pow_bound(a in finite(-40..=40), n in 0u64..50) {
        let p = value(&a).unwrap();
        let exact = p.pow(i64::try_from(n).unwrap());
        // (1 + 2^-29)^(2n) <= 1 + (n + 1) 2^-28 for these n.
        let slack = i64::try_from(n).unwrap() + 1;
        check_upper("pow", &a.pow(n), &exact, slack)?;
    }

    #[test]
    fn conversions(m in crate::testkit::strategy::signed_up_to(3), e in -200i64..200) {
        let x = crate::ball::Float::from_integer_exp(&m.to_integer(), &Integer::from(e));
        let exact = x.to_rational().abs();
        if exact.is_zero() {
            prop_assert!(Mag::from_float(&x).is_zero() && Mag::from_float_lower(&x).is_zero());
        } else {
            check_upper("from_float", &Mag::from_float(&x), &exact, 1)?;
            check_lower("from_float_lower", &Mag::from_float_lower(&x), &exact, 1)?;
        }
        let n = m.to_integer();
        let nv = Rational::from(n.abs());
        if !n.is_zero() {
            check_upper("from_integer", &Mag::from_integer(&n), &nv, 1)?;
            check_lower("from_integer_lower", &Mag::from_integer_lower(&n), &nv, 1)?;
        }
    }

    #[test]
    fn exp_bounds(a in finite(-60..=28)) {
        let (lo, hi) = mpfr_bracket(&a, |x, p, r| MpFloat::with_val_round(p, x.exp_ref(), r).0);
        let got = a.exp();
        let v = value(&got).unwrap();
        prop_assert!(v >= hi, "exp {:?} below e^x", got);
        prop_assert!(v <= &lo * &(&Rational::ONE + &tight()), "exp {:?} too far above", got);
    }

    #[test]
    fn expm1_bounds(a in finite(-60..=28)) {
        let (lo, hi) = mpfr_bracket(&a, |x, p, r| MpFloat::with_val_round(p, x.exp_m1_ref(), r).0);
        let got = a.expm1();
        let v = value(&got).unwrap();
        prop_assert!(v >= hi, "expm1 {:?} below e^x - 1", got);
        prop_assert!(v <= &lo * &(&Rational::ONE + &tight()), "expm1 {:?} too far above", got);
    }

    #[test]
    fn log1p_bounds(a in prop_oneof![finite(-60..=70), finite(70..=1000)]) {
        let (lo, hi) = mpfr_bracket(&a, |x, p, r| MpFloat::with_val_round(p, x.ln_1p_ref(), r).0);
        let got = a.log1p();
        let v = value(&got).unwrap();
        prop_assert!(v >= hi, "log1p {:?} below log(1 + x)", got);
        prop_assert!(v <= &lo * &(&Rational::ONE + &tight()), "log1p {:?} too far above", got);
    }

    /// Exponents near and beyond ±2^63: scaling commutes with every
    /// operation (and cancels in products against the inverse scale).
    #[test]
    fn huge_exponents(a in finite(-100..=100), b in finite(-100..=100), k in huge_shift()) {
        let k = Signed::from_rug(&k).to_integer();
        let (ak, bk) = (a.mul_2exp(&k), b.mul_2exp(&k));
        prop_assert_eq!(ak.add(&bk), a.add(&b).mul_2exp(&k));
        prop_assert_eq!(ak.add_lower(&bk), a.add_lower(&b).mul_2exp(&k));
        prop_assert_eq!(ak.sub(&bk), a.sub(&b).mul_2exp(&k));
        prop_assert_eq!(ak.sub_lower(&bk), a.sub_lower(&b).mul_2exp(&k));
        prop_assert_eq!(ak.div(&bk), a.div(&b));
        prop_assert_eq!(ak.mul(&b.mul_2exp(&-&k)), a.mul(&b));
        prop_assert_eq!(ak.cmp(&bk), a.cmp(&b));
        let f = a.to_float().unwrap().mul_2exp(&k);
        prop_assert_eq!(Mag::from_float(&f), Mag::from_float(&a.to_float().unwrap()).mul_2exp(&k));
    }
}

#[test]
fn zero_and_infinity() {
    let one = Mag::ONE;
    let (z, inf) = (Mag::ZERO, Mag::INF);
    assert_eq!(value(&one), Some(Rational::ONE));
    assert!(z < one && one < inf && z < inf);
    assert_eq!(one.add(&inf), inf);
    assert_eq!(one.mul(&inf), inf);
    assert_eq!(z.mul(&inf), z);
    assert_eq!(one.div(&z), inf);
    assert_eq!(one.div(&inf), z);
    assert_eq!(inf.div(&inf), inf);
    assert_eq!(inf.div_lower(&inf), z);
    assert_eq!(inf.sub(&inf), inf);
    assert_eq!(inf.sub_lower(&inf), z);
    assert_eq!(one.sub(&one), z);
    assert_eq!(one.sub(&inf), z);
    assert_eq!(z.exp(), one);
    assert_eq!(inf.exp(), inf);
    assert_eq!(z.expm1(), z);
    assert_eq!(z.log1p(), z);
    assert_eq!(inf.log1p(), inf);
    assert_eq!(one.pow(0), one);
    assert_eq!(z.pow(0), one);
    assert_eq!(inf.pow(3), inf);
    assert_eq!(Mag::from_float(&crate::ball::Float::ZERO), z);
    assert_eq!(inf.to_float(), None);
    assert_eq!(
        value(&Mag::power_of_two(&Integer::from(-5))),
        Some(Rational::new(Integer::ONE, Integer::from(32)))
    );
}

/// exp at the edges of its paths: 2^-31 (the tiny bound), 2^62 (saturation),
/// and huge exponents.
#[test]
fn exp_saturation() {
    let tiny = Mag::new(HIGH - 1, Integer::from(-31));
    let next = Mag::ONE.next_up();
    assert_eq!(tiny.exp(), next);
    assert_eq!(Mag::new(LOW, Integer::from(63)).exp(), Mag::INF);
    assert_eq!(
        Mag::new(LOW, Integer::from(crate::Natural::power_of_two(70))).exp(),
        Mag::INF
    );
    assert_eq!(
        Mag::new(LOW, -(Integer::from(crate::Natural::power_of_two(70)))).exp(),
        next
    );
    // e^(2^61) is finite, with an exponent near 2^61 / ln 2.
    let big = Mag::new(LOW, Integer::from(62)).exp();
    assert!(big.is_finite());
    let want = 3_326_628_274_461_080_622i64; // floor(2^61 / ln 2)
    let e = big.exponent().to_i64().unwrap();
    assert!((want..=want + 2).contains(&e), "exponent {e}");
    // log1p with an exponent beyond i64: about exp × ln 2.
    let huge = Mag::new(LOW, Integer::from(crate::Natural::power_of_two(70)));
    let l = huge.log1p();
    assert!(l.is_finite() && l.exponent() == &Integer::from(70));
}
