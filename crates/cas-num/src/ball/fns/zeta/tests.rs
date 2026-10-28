//! Tests of `Ball::zeta` (CAS-244) against MPFR's zeta: values with
//! `prec - 3` bits, the exact values, the pole, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses, dyadic_with};

fn mpfr(x: &MpFloat, p: u32, r: MpRound) -> MpFloat {
    MpFloat::with_val_round(p, x.zeta_ref(), r).0
}

/// `x` as an MPFR float, exactly.
fn to_mp(x: &Float) -> MpFloat {
    MpFloat::with_val(10_000, &x.to_rational().to_rug())
}

/// `ζ(x)` for an exact `x` at `prec` bits encloses MPFR's bracket and has
/// `prec - 3` bits.
fn check(x: &Float, prec: u64) -> Result<(), String> {
    let b = Ball::exact(x.clone()).zeta(prec);
    if !b.is_finite() {
        return Err(format!("ζ({x:?}) at {prec}: indeterminate"));
    }
    let xr = to_mp(x);
    let br = bracket(u32::try_from(prec).unwrap(), |p, r| mpfr(&xr, p, r));
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br).map_err(|e| format!("ζ({x:?}) at {prec}: {e}"))?;
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let want = prec as i64 - 3;
    if !b.mid().is_zero() && b.rel_accuracy_bits() < want {
        return Err(format!(
            "ζ({x:?}) at {prec}: {} bits: {b:?}",
            b.rel_accuracy_bits()
        ));
    }
    Ok(())
}

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

/// Borwein's algorithm above 1/2, the functional equation below.
#[test]
fn values() {
    let mut points: Vec<Float> = [
        2.0, 3.0, 4.0, 0.5, 0.75, 1.5, 2.5, 10.0, 100.5, -0.5, -1.0, -3.0, -1.5, -20.5, -100.25,
        0.25, 0.4999,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    // Next to the pole and to 0.
    points.push(Float::from(1i64).add(&Float::from(1i64).mul_2exp(&Integer::from(-50i64))));
    points.push(Float::from(1i64).sub(&Float::from(1i64).mul_2exp(&Integer::from(-30i64))));
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-40i64)));
    points.push(Float::from(-1i64).mul_2exp(&Integer::from(-40i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for prec in [24, 53, 200, 1000] {
            if let Err(e) = check(x, prec) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Past 3000 bits Borwein's tail bound lost `0.04 n` bits: `ζ(-63.5)` at
/// 8489 bits had 8425 (`ball_fns` fuzz, CAS-248).
#[test]
fn high_precision() {
    let mut wrong = Vec::new();
    for x in [2.5, 0.75, -63.5] {
        if let Err(e) = check(&float(x), 4000) {
            wrong.push(e);
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn exact_values() {
    let at = |x: i64| Ball::exact(Float::from(x)).zeta(64);
    assert_eq!(
        at(0),
        Ball::exact(Float::from(-1i64).mul_2exp(&-Integer::ONE))
    );
    for n in [-2i64, -4, -100] {
        assert_eq!(at(n), Ball::ZERO, "{n}");
    }
    // Above the precision: 1 to within 2^-σ.
    let big = Ball::exact(Float::from(5000i64)).zeta(100);
    assert!(big.contains(&Float::from(1i64)) && big.rel_accuracy_bits() >= 97);
    // Huge arguments, without rationals (CAS-249): a trivial zero of 2^40
    // bits, and a radius a rational holds past 2^47.
    let huge = Float::from(1i64).mul_2exp(&Integer::from(1i64 << 40));
    assert_eq!(Ball::exact(huge.neg()).zeta(64), Ball::ZERO);
    check(&Float::from(1i64).mul_2exp(&Integer::from(47i64)), 2).unwrap();
    // Tiny arguments: 1 − s rounded, not exact over 2^68 bits (CAS-252),
    // nor a million (a timeout, CAS-254); ζ next to ζ(0) = −1/2.
    let small = Float::from(-1i64).mul_2exp(&Integer::from(-1_048_559i64));
    check(&small, 200).unwrap();
    let e = -(&Integer::from(2).pow(68) + &Integer::from(496));
    let tiny = Ball::exact(Float::from(-16383i64).mul_2exp(&e)).zeta(49);
    let half = Float::from(-1i64).mul_2exp(&-Integer::ONE);
    assert!(
        tiny.is_finite() && tiny.rel_accuracy_bits() >= 46,
        "{tiny:?}"
    );
    assert!(tiny.mid().sub(&half).abs() <= Float::from(1i64).mul_2exp(&Integer::from(-40i64)));
}

/// Exact arguments of more than 64 bits next to 1/2: the ends rounded to
/// 64 bits made `1/2 − 2^−66` straddle 1/2, indeterminate (CAS-280).
#[test]
fn next_to_one_half() {
    let bit = |e: i64| Float::from(1i64).mul_2exp(&Integer::from(e));
    let half = bit(-1);
    let mut wrong = Vec::new();
    for x in [
        half.sub(&bit(-66)),
        half.add(&bit(-66)),
        half.sub(&bit(-200)),
    ] {
        for prec in [2, 53] {
            if let Err(e) = check(&x, prec) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn the_pole() {
    assert_eq!(Ball::exact(Float::from(1i64)).zeta(64), Ball::INDETERMINATE);
    let near = Ball::new(float(1.0), Mag::power_of_two(&Integer::from(-20i64)));
    assert_eq!(near.zeta(64), Ball::INDETERMINATE);
    // Straddling 1/2.
    let half = Ball::new(float(0.5), Mag::power_of_two(&Integer::from(-20i64)));
    assert_eq!(half.zeta(64), Ball::INDETERMINATE);
}

/// A ball's value holds the values at its ends.
#[test]
fn radii() {
    for (m, e) in [
        (2.5, -10i64),
        (0.75, -12),
        (-1.25, -8),
        (-7.5, -20),
        (30.0, -30),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let b = Ball::new(float(m), r.clone()).zeta(100);
        assert!(b.is_finite(), "ζ({m} ± 2^{e})");
        let rf = r.to_float().unwrap();
        for end in [float(m).sub(&rf), float(m).add(&rf)] {
            let xr = to_mp(&end);
            let br = bracket(100, |p, rd| mpfr(&xr, p, rd));
            let mid = b.mid().to_rational().to_rug();
            let rad = b.rad().to_float().unwrap().to_rational().to_rug();
            check_encloses(&mid, &rad, &br)
                .unwrap_or_else(|err| panic!("ζ({m} ± 2^{e}) at {end:?}: {err}"));
        }
    }
}

proptest! {
    #[test]
    fn against_mpfr(
        x in dyadic_with(2, 0..=0),
        top in -30i64..=12,
        prec in prop_oneof![2u64..120, 120u64..500],
    ) {
        // Exponents up to 12: |ζ| within MPFR's range at the negative ones.
        let x = x.to_ours();
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let one = Float::from(1i64);
        let half = one.mul_2exp(&-Integer::ONE);
        if x == one || x == half {
            return Ok(());
        }
        check(&x, prec).map_err(TestCaseError::fail)?;
    }
}
