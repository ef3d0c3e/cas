//! Tests of `Ball::ei`, `Ball::si` and `Ball::ci` (CAS-350): `ei` against
//! MPFR's `eint`, `si` and `ci` against their series summed in MPFR far
//! above the precision; values with `prec - 3` bits next to the zeros too,
//! poles and domains, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Constant;
use rug::ops::NegAssign;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{Bracket, bracket, check_encloses, dyadic_with};

/// `x` as an MPFR float, exactly.
fn to_mp(x: &Float) -> MpFloat {
    MpFloat::with_val(10_000, &x.to_rational().to_rug())
}

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

/// MPFR's bracket of `Ei(x)` at `prec` bits.
fn ei_bracket(x: &Float, prec: u64) -> Bracket {
    let xr = to_mp(x);
    bracket(u32::try_from(prec).unwrap(), |p, r| {
        MpFloat::with_val_round(p, xr.eint_ref(), r).0
    })
}

/// `Si(x)` or `Ci(x)` by the series in MPFR, at enough bits for an
/// absolute error far below `2^-(prec + 64)` of the result's size: a
/// reference, bracketed by that much.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    reason = "small test sizes"
)]
fn series_bracket(x: &Float, prec: u64, cosine: bool) -> Bracket {
    let size = x.to_f64().abs();
    // The terms reach e^|x|, the result may be as small as 1/|x|.
    let p = (prec as f64 + 2.0 * size * std::f64::consts::LOG2_E + 256.0) as u32;
    let xr = MpFloat::with_val(p, to_mp(x));
    let square = MpFloat::with_val(p, &xr * &xr);
    let (mut power, mut sum, mut j) = if cosine {
        (MpFloat::with_val(p, 1), MpFloat::new(p), 0u32)
    } else {
        (xr.clone(), xr.clone(), 1u32)
    };
    let tiny = MpFloat::with_val(p, MpFloat::i_exp(1, -(p as i32)));
    loop {
        j += 2;
        power *= &square;
        power /= (j - 1) * j;
        power.neg_assign();
        let term = MpFloat::with_val(p, &power / j);
        sum += &term;
        if f64::from(j) > 2.0 * size && term.clone().abs() < tiny {
            break;
        }
    }
    if cosine {
        sum += MpFloat::with_val(p, Constant::Euler);
        sum += MpFloat::with_val(p, xr.ln_ref());
    }
    let slack = MpFloat::with_val(p, MpFloat::i_exp(1, -(prec as i32) - 64));
    let lo = MpFloat::with_val(p, &sum - &slack);
    let hi = MpFloat::with_val(p, &sum + &slack);
    Bracket {
        lo: lo.to_rational().unwrap(),
        hi: hi.to_rational().unwrap(),
    }
}

type Ours = fn(&Ball, u64) -> Ball;

/// `f(x)` for an exact `x` at `prec` bits meets the bracket and has
/// `prec - 3` bits.
fn check(name: &str, f: Ours, x: &Float, prec: u64, oracle: &Bracket) -> Result<(), String> {
    let b = f(&Ball::exact(x.clone()), prec);
    if !b.is_finite() {
        return Err(format!("{name}({x:?}) at {prec}: indeterminate"));
    }
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, oracle).map_err(|e| format!("{name}({x:?}) at {prec}: {e}"))?;
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

/// The series, the asymptotic expansions on both sides of 0, and the zero
/// of `Ei` at 0.3725….
#[test]
fn ei_values() {
    let mut points: Vec<Float> = [
        0.5,
        1.0,
        2.0,
        0.372_507_410_781_366_6,
        5.0,
        10.0,
        30.0,
        50.0,
        100.0,
        300.0,
        1000.0,
        1e5,
        -0.5,
        -1.0,
        -5.0,
        -30.0,
        -100.0,
        -1000.0,
        -1e5,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-60i64)));
    points.push(Float::from(-1i64).mul_2exp(&Integer::from(-60i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for prec in [24, 53, 200, 1000, 3000] {
            if let Err(e) = check("ei", Ball::ei, x, prec, &ei_bracket(x, prec)) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The series, and the auxiliary functions' expansions from about 60 at
/// low precisions on; the first zero of `Ci`, 0.6165….
#[test]
fn si_ci_values() {
    let mut points: Vec<Float> = [
        0.5,
        1.0,
        2.0,
        0.616_505_485_620_716_2,
        3.384_180_422_551_186,
        5.0,
        10.0,
        30.0,
        80.0,
        200.0,
        700.0,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-60i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for prec in [24, 53, 200, 1000] {
            for (name, f, cosine) in [("si", Ball::si as Ours, false), ("ci", Ball::ci, true)] {
                let oracle = series_bracket(x, prec, cosine);
                if let Err(e) = check(name, f, x, prec, &oracle) {
                    wrong.push(e);
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Far out, against giac's values to 50 digits.
#[test]
fn si_ci_far() {
    let x = float(1e5);
    let si = Ball::exact(x.clone()).si(150);
    let ci = Ball::exact(x).ci(150);
    let near = |b: &Ball, digits: &str| {
        let want = MpFloat::with_val(400, MpFloat::parse(digits).unwrap());
        let mid = MpFloat::with_val(400, b.mid().to_rational().to_rug());
        let error = MpFloat::with_val(400, &mid - &want).abs();
        // giac's 50 digits, and the ball's own radius.
        let radius = MpFloat::with_val(400, b.rad().to_float().unwrap().to_rational().to_rug());
        let allowed = MpFloat::with_val(400, &radius + MpFloat::with_val(400, &want * 1e-48).abs());
        assert!(error <= allowed, "{b:?} against {digits}");
        assert!(b.rel_accuracy_bits() >= 147, "{b:?}");
    };
    near(&si, "1.5708063203993941228391714019773038489126828018428");
    near(
        &ci,
        "0.35758791572935135690333184729128548684762053177217e-6",
    );
}

/// `Si` is odd and 0 at 0; `Ei` and `Ci` have their pole at 0, and `Ci` is
/// for positive arguments.
#[test]
fn domains() {
    assert_eq!(Ball::exact(Float::ZERO).si(64), Ball::ZERO);
    for x in [0.5, 3.0, 150.0] {
        let plus = Ball::exact(float(x)).si(100);
        let minus = Ball::exact(float(-x)).si(100);
        assert_eq!(minus, plus.neg(), "si(-{x})");
    }
    let around = Ball::new(Float::ZERO, Mag::power_of_two(&Integer::from(-10)));
    let near = Ball::new(float(0.001), Mag::power_of_two(&Integer::from(-8)));
    for b in [Ball::exact(Float::ZERO), around.clone(), near] {
        assert_eq!(b.ei(64), Ball::INDETERMINATE, "ei({b:?})");
        assert_eq!(b.ci(64), Ball::INDETERMINATE, "ci({b:?})");
    }
    assert!(around.si(64).is_finite());
    assert_eq!(Ball::exact(float(-1.0)).ci(64), Ball::INDETERMINATE);
    // Past exp's range Ei is indeterminate above 0, about 0 below.
    let huge = Float::from(1i64).mul_2exp(&Integer::from(1i64 << 20));
    assert_eq!(Ball::exact(huge.clone()).ei(64), Ball::INDETERMINATE);
    let below = Ball::exact(huge.neg()).ei(64);
    assert!(
        below.is_finite() && below.contains(&Float::ZERO),
        "{below:?}"
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
        (120.0, -40),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let rf = r.to_float().unwrap();
        let ball = Ball::new(float(m), r);
        for end in [float(m).sub(&rf), float(m).add(&rf)] {
            let mut cases = vec![("ei", ball.ei(100), ei_bracket(&end, 100))];
            cases.push(("si", ball.si(100), series_bracket(&end, 100, false)));
            if m > 0.0 {
                cases.push(("ci", ball.ci(100), series_bracket(&end, 100, true)));
            }
            for (name, b, oracle) in cases {
                assert!(b.is_finite(), "{name}({m} ± 2^{e})");
                let mid = b.mid().to_rational().to_rug();
                let rad = b.rad().to_float().unwrap().to_rational().to_rug();
                check_encloses(&mid, &rad, &oracle)
                    .unwrap_or_else(|err| panic!("{name}({m} ± 2^{e}) at {end:?}: {err}"));
            }
        }
    }
}

proptest! {
    #[test]
    fn ei_against_mpfr(
        x in dyadic_with(2, 0..=0),
        top in -40i64..=10,
        prec in prop_oneof![2u64..120, 120u64..600],
    ) {
        let x = x.to_ours();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        check("ei", Ball::ei, &x, prec, &ei_bracket(&x, prec)).map_err(TestCaseError::fail)?;
    }

    #[test]
    fn si_ci_against_series(
        x in dyadic_with(1, 0..=0),
        top in -20i64..=7,
        prec in 2u64..300,
    ) {
        let x = x.to_ours().abs();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        check("si", Ball::si, &x, prec, &series_bracket(&x, prec, false))
            .map_err(TestCaseError::fail)?;
        check("ci", Ball::ci, &x, prec, &series_bracket(&x, prec, true))
            .map_err(TestCaseError::fail)?;
    }
}
