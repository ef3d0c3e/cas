//! Tests of `Ball::gamma` (CAS-242) against MPFR's gamma: values with
//! `prec - 3` bits, the poles, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses, dyadic_with};

fn mpfr(x: &MpFloat, p: u32, r: MpRound) -> MpFloat {
    MpFloat::with_val_round(p, x.gamma_ref(), r).0
}

/// `x` as an MPFR float, exactly.
fn to_mp(x: &Float) -> MpFloat {
    MpFloat::with_val(10_000, &x.to_rational().to_rug())
}

/// `Γ(x)` for an exact `x` at `prec` bits encloses MPFR's bracket and has
/// `prec - 3` bits.
fn check(x: &Float, prec: u64) -> Result<(), String> {
    let b = Ball::exact(x.clone()).gamma(prec);
    if !b.is_finite() {
        return Err(format!("Γ({x:?}) at {prec}: indeterminate"));
    }
    let xr = to_mp(x);
    let br = bracket(u32::try_from(prec).unwrap(), |p, r| mpfr(&xr, p, r));
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br).map_err(|e| format!("Γ({x:?}) at {prec}: {e}"))?;
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let want = prec as i64 - 3;
    if !b.mid().is_zero() && b.rel_accuracy_bits() < want {
        return Err(format!(
            "Γ({x:?}) at {prec}: {} bits: {b:?}",
            b.rel_accuracy_bits()
        ));
    }
    Ok(())
}

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

/// Brent's series for small arguments, Stirling's for large ones, the
/// recurrence below 1/2 and above −64, the reflection below.
#[test]
fn values() {
    let mut points: Vec<Float> = [
        0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 0.75, 0.3, 0.001, 10.5, 40.0, 100.25, 171.5, 1000.5,
        123_456.75, -0.5, -1.5, -2.25, -20.5, -63.75, -100.5, -1000.25,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-50i64)));
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-1000i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for prec in [24, 53, 200, 1000, 3000] {
            if let Err(e) = check(x, prec) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // Exactly: Γ(5) = 24.
    assert!(
        Ball::exact(Float::from(5i64))
            .gamma(64)
            .contains_integer(&Integer::from(24))
    );
}

/// The poles, and balls that meet them, are indeterminate.
#[test]
fn poles() {
    for n in [0i64, -1, -2, -10, -1000] {
        assert_eq!(
            Ball::exact(Float::from(n)).gamma(64),
            Ball::INDETERMINATE,
            "{n}"
        );
    }
    let near = |m: f64, e: i64| Ball::new(float(m), Mag::power_of_two(&Integer::from(e)));
    assert_eq!(near(-1.0, -4).gamma(64), Ball::INDETERMINATE);
    assert_eq!(near(0.25, -1).gamma(64), Ball::INDETERMINATE);
    assert_eq!(near(-3.25, 0).gamma(64), Ball::INDETERMINATE);
    // Close to a pole but not on it: finite.
    assert!(near(-1.5, -4).gamma(64).is_finite());
    // Closer than 64 bits tell: decided exactly (−8 + 2^-164 at 228 bits
    // was indeterminate, CAS-250), and a radius that reaches the pole.
    let next = Float::from(-8i64).add(&Float::from(1i64).mul_2exp(&Integer::from(-164i64)));
    check(&next, 228).unwrap();
    let reach = Ball::new(next.clone(), Mag::power_of_two(&Integer::from(-163i64)));
    assert_eq!(reach.gamma(64), Ball::INDETERMINATE);
    let short = Ball::new(next, Mag::power_of_two(&Integer::from(-170i64)));
    assert!(short.gamma(64).is_finite());
    // A tiny negative argument, beyond MPFR's exponents: Γ(−ε) ≈ −1/ε,
    // shifted by two steps without a rational (CAS-253).
    let tiny = Float::from(-1i64).mul_2exp(&-Integer::from(1i64 << 40));
    let g = Ball::exact(tiny).gamma(64);
    assert!(g.is_finite() && g.rel_accuracy_bits() >= 61, "{g:?}");
    let scaled = g.mid().mul_2exp(&-Integer::from(1i64 << 40));
    let near = Float::from(1i64).mul_2exp(&Integer::from(-60i64));
    assert!(scaled.add(&Float::from(1i64)).abs() <= near, "{g:?}");
}

/// A ball's value holds the values at its ends.
#[test]
fn radii() {
    for (m, e) in [
        (1.5, -10i64),
        (0.75, -6),
        (7.25, -20),
        (-2.5, -8),
        (300.5, -30),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let b = Ball::new(float(m), r.clone()).gamma(100);
        assert!(b.is_finite(), "Γ({m} ± 2^{e})");
        let rf = r.to_float().unwrap();
        for end in [float(m).sub(&rf), float(m).add(&rf)] {
            let xr = to_mp(&end);
            let br = bracket(100, |p, rd| mpfr(&xr, p, rd));
            let mid = b.mid().to_rational().to_rug();
            let rad = b.rad().to_float().unwrap().to_rational().to_rug();
            check_encloses(&mid, &rad, &br)
                .unwrap_or_else(|err| panic!("Γ({m} ± 2^{e}) at {end:?}: {err}"));
        }
    }
}

proptest! {
    #[test]
    fn against_mpfr(x in dyadic_with(2, -40..=12), prec in prop_oneof![2u64..120, 120u64..600, Just(1500u64)]) {
        let x = x.to_ours();
        // Below 2^12, where MPFR's gamma has a value (two limbs of mantissa
        // reach 2^140), and away from the poles.
        if x.is_zero()
            || *x.exponent() > Integer::from(12)
            || (x.is_negative() && x.to_rational().is_integer())
        {
            return Ok(());
        }
        check(&x, prec).map_err(TestCaseError::fail)?;
    }
}

/// The Bernoulli numbers from the tangent numbers.
#[test]
fn bernoulli_numbers() {
    use crate::Rational;
    use crate::ball::bernoulli;
    let q = |n: i64, d: i64| Some(Rational::new(Integer::from(n), Integer::from(d)));
    assert_eq!(bernoulli(0), q(1, 1));
    assert_eq!(bernoulli(1), q(-1, 2));
    assert_eq!(bernoulli(2), q(1, 6));
    assert_eq!(bernoulli(3), q(0, 1));
    assert_eq!(bernoulli(4), q(-1, 30));
    assert_eq!(bernoulli(6), q(1, 42));
    assert_eq!(bernoulli(10), q(5, 66));
    assert_eq!(bernoulli(12), q(-691, 2730));
    assert_eq!(bernoulli(20), q(-174_611, 330));
    assert_eq!(bernoulli(100_000), None);
}
