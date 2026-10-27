//! Tests of `Ball::erf` and `Ball::erfc` (CAS-246) against MPFR's: values
//! with `prec - 3` bits, small to large arguments, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses, dyadic_with};

type Oracle = fn(&MpFloat, u32, MpRound) -> MpFloat;
type Ours = fn(&Ball, u64) -> Ball;

const FUNS: [(&str, Ours, Oracle); 2] = [
    ("erf", Ball::erf, |x, p, r| {
        MpFloat::with_val_round(p, x.erf_ref(), r).0
    }),
    ("erfc", Ball::erfc, |x, p, r| {
        MpFloat::with_val_round(p, x.erfc_ref(), r).0
    }),
];

/// `x` as an MPFR float, exactly.
fn to_mp(x: &Float) -> MpFloat {
    MpFloat::with_val(10_000, &x.to_rational().to_rug())
}

/// `f(x)` for an exact `x` at `prec` bits encloses MPFR's bracket and has
/// `prec - 3` bits.
fn check(name: &str, f: Ours, oracle: Oracle, x: &Float, prec: u64) -> Result<(), String> {
    let b = f(&Ball::exact(x.clone()), prec);
    if !b.is_finite() {
        return Err(format!("{name}({x:?}) at {prec}: indeterminate"));
    }
    let xr = to_mp(x);
    let br = bracket(u32::try_from(prec).unwrap(), |p, r| oracle(&xr, p, r));
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br).map_err(|e| format!("{name}({x:?}) at {prec}: {e}"))?;
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let want = prec as i64 - 3;
    if !b.mid().is_zero() && b.rel_accuracy_bits() < want {
        return Err(format!(
            "{name}({x:?}) at {prec}: {} bits: {b:?}",
            b.rel_accuracy_bits()
        ));
    }
    Ok(())
}

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

/// The series, `1 − erf` with its extra bits, and the asymptotic expansion.
#[test]
fn values() {
    let mut points: Vec<Float> = [
        0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 27.5, 60.0, 100.0, -0.5, -3.0, -10.0, 0.001,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-60i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for (name, f, oracle) in FUNS {
            for prec in [24, 53, 200, 1000, 3000] {
                if let Err(e) = check(name, f, oracle, x, prec) {
                    wrong.push(e);
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // At 0.
    assert_eq!(Ball::exact(Float::ZERO).erf(64), Ball::ZERO);
    assert_eq!(
        Ball::exact(Float::ZERO).erfc(64),
        Ball::exact(Float::from(1i64))
    );
}

/// A ball's value holds the values at its ends.
#[test]
fn radii() {
    for (m, e) in [
        (0.5, -10i64),
        (-1.25, -12),
        (3.0, -20),
        (8.0, -30),
        (0.0, -6),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let rf = r.to_float().unwrap();
        for (name, f, oracle) in FUNS {
            let b = f(&Ball::new(float(m), r.clone()), 100);
            assert!(b.is_finite(), "{name}({m} ± 2^{e})");
            for end in [float(m).sub(&rf), float(m).add(&rf)] {
                let xr = to_mp(&end);
                let br = bracket(100, |p, rd| oracle(&xr, p, rd));
                let mid = b.mid().to_rational().to_rug();
                let rad = b.rad().to_float().unwrap().to_rational().to_rug();
                check_encloses(&mid, &rad, &br)
                    .unwrap_or_else(|err| panic!("{name}({m} ± 2^{e}) at {end:?}: {err}"));
            }
        }
    }
}

proptest! {
    #[test]
    fn against_mpfr(
        x in dyadic_with(2, 0..=0),
        top in -40i64..=6,
        prec in prop_oneof![2u64..120, 120u64..600],
    ) {
        // Below 2^6: erfc's value stays in MPFR's exponent range.
        let x = x.to_ours();
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        for (name, f, oracle) in FUNS {
            check(name, f, oracle, &x, prec).map_err(TestCaseError::fail)?;
        }
    }
}
