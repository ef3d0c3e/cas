//! Tests of `Ball::digamma` and `Ball::polygamma` (CAS-351): `ψ` against
//! MPFR's digamma, `ψ^(m)` against its exact values at integers and
//! half-integers through MPFR's zeta, the recurrence and the reflection,
//! the poles, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;

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

/// `b` meets the bracket and, unless it is 0, has `prec - 3` bits.
fn check(what: &str, b: &Ball, prec: u64, oracle: &Bracket) -> Result<(), String> {
    if !b.is_finite() {
        return Err(format!("{what} at {prec}: indeterminate"));
    }
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, oracle).map_err(|e| format!("{what} at {prec}: {e}"))?;
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let want = prec as i64 - 3;
    if !b.mid().is_zero() && b.rel_accuracy_bits() < want {
        return Err(format!(
            "{what} at {prec}: {} bits: {b:?}",
            b.rel_accuracy_bits()
        ));
    }
    Ok(())
}

/// MPFR's bracket of `ψ(x)`.
fn digamma_bracket(x: &Float, prec: u64) -> Bracket {
    let xr = to_mp(x);
    bracket(u32::try_from(prec).unwrap(), |p, r| {
        MpFloat::with_val_round(p, xr.digamma_ref(), r).0
    })
}

/// `ψ^(m)(x)` for `x = n/2`, an integer or a half-integer, from `ψ^(m)(1)
/// = (−1)^(m+1) m! ζ(m + 1)` or `ψ^(m)(1/2) = (−1)^(m+1) m! (2^(m+1) − 1)
/// ζ(m + 1)` and the recurrence `ψ^(m)(x + 1) = ψ^(m)(x) + (−1)^m m!/x^(m+1)`,
/// in MPFR at `prec + 256` bits, bracketed by `2^-(prec + 64)` of it.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "small test sizes"
)]
fn exact_bracket(n: i64, m: u32, prec: u64) -> Bracket {
    let p = u32::try_from(prec).unwrap() + 256;
    let factorial = MpFloat::with_val(p, MpFloat::factorial(m));
    let zeta = MpFloat::with_val(p, MpFloat::with_val(p, m + 1).zeta());
    let half = n % 2 != 0;
    // The start, 1 or 1/2.
    let mut value = MpFloat::with_val(p, &factorial * &zeta);
    if half {
        let scale = MpFloat::with_val(p, MpFloat::i_exp(1, m as i32 + 1)) - 1u32;
        value *= scale;
    }
    if m % 2 == 0 {
        value = -value;
    }
    let sign = if m % 2 == 0 { 1 } else { -1 };
    // x runs from the start to n/2 by steps of 1.
    let start = if half { 1 } else { 2 };
    let step = |x2: i64| MpFloat::with_val(p, MpFloat::with_val(p, x2) / 2u32);
    let power = |x: &MpFloat| -> MpFloat {
        let mut y = MpFloat::with_val(p, 1);
        for _ in 0..=m {
            y *= x;
        }
        y
    };
    let mut x2 = start;
    while x2 < n {
        // ψ(x + 1) = ψ(x) + (−1)^m m!/x^(m+1).
        let term = MpFloat::with_val(p, &factorial / power(&step(x2)));
        value += term * sign;
        x2 += 2;
    }
    while x2 > n {
        // ψ(x − 1) = ψ(x) − (−1)^m m!/(x − 1)^(m+1).
        x2 -= 2;
        let term = MpFloat::with_val(p, &factorial / power(&step(x2)));
        value -= term * sign;
    }
    let unit = MpFloat::with_val(p, MpFloat::i_exp(1, -(prec as i32) - 64));
    let slack = MpFloat::with_val(p, &value * &unit).abs();
    Bracket {
        lo: MpFloat::with_val(p, &value - &slack).to_rational().unwrap(),
        hi: MpFloat::with_val(p, &value + &slack).to_rational().unwrap(),
    }
}

/// The shift, the expansion at once for large arguments, the reflection
/// below −64, and the zero at 1.4616….
#[test]
fn digamma_values() {
    let mut points: Vec<Float> = [
        0.5,
        1.0,
        1.461_632_144_968_362_3,
        2.0,
        3.5,
        10.0,
        100.25,
        1e5,
        1e30,
        -0.5,
        -1.5,
        -2.25,
        -20.5,
        -63.75,
        -100.5,
        -1000.25,
        -0.504_083_008_264_455_4,
    ]
    .iter()
    .map(|&x| float(x))
    .collect();
    points.push(Float::from(1i64).mul_2exp(&Integer::from(-50i64)));
    let mut wrong = Vec::new();
    for x in &points {
        for prec in [24, 53, 200, 1000, 3000] {
            let b = Ball::exact(x.clone()).digamma(prec);
            if let Err(e) = check(&format!("ψ({x:?})"), &b, prec, &digamma_bracket(x, prec)) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Orders 1 to 8 and 40 at integers and half-integers on both sides of 0,
/// the reflection's among them.
#[test]
fn polygamma_values() {
    let mut wrong = Vec::new();
    for m in [1u32, 2, 3, 4, 7, 8, 40] {
        for n in [2i64, 1, 3, 10, 41, -1, -3, -7, -201, 400] {
            let x = Float::from(n).mul_2exp(&-Integer::ONE);
            for prec in [24, 53, 200, 1000] {
                let b = Ball::exact(x.clone()).polygamma(u64::from(m), prec);
                let what = format!("ψ^({m})({n}/2)");
                if let Err(e) = check(&what, &b, prec, &exact_bracket(n, m, prec)) {
                    wrong.push(e);
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The poles, balls meeting them, and orders past the bound.
#[test]
fn poles() {
    for x in [0.0, -1.0, -5.0, -1000.0] {
        for m in [0, 1, 4] {
            assert_eq!(
                Ball::exact(float(x)).polygamma(m, 64),
                Ball::INDETERMINATE,
                "ψ^({m})({x})"
            );
        }
    }
    let around = Ball::new(float(-3.0), Mag::power_of_two(&Integer::from(-20)));
    assert_eq!(around.digamma(64), Ball::INDETERMINATE);
    let wide = Ball::new(float(-2.5), Mag::ONE);
    assert_eq!(wide.digamma(64), Ball::INDETERMINATE);
    assert_eq!(
        Ball::exact(float(1.5)).polygamma(257, 64),
        Ball::INDETERMINATE
    );
    assert!(Ball::exact(float(1.5)).polygamma(256, 64).is_finite());
}

/// A ball's value holds the values at its ends.
#[test]
fn radii() {
    for (m, e) in [
        (0.5, -10i64),
        (1.25, -12),
        (-1.25, -12),
        (30.0, -20),
        (-80.5, -30),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let rf = r.to_float().unwrap();
        let b = Ball::new(float(m), r).digamma(100);
        assert!(b.is_finite(), "ψ({m} ± 2^{e})");
        for end in [float(m).sub(&rf), float(m).add(&rf)] {
            let mid = b.mid().to_rational().to_rug();
            let rad = b.rad().to_float().unwrap().to_rational().to_rug();
            check_encloses(&mid, &rad, &digamma_bracket(&end, 100))
                .unwrap_or_else(|err| panic!("ψ({m} ± 2^{e}) at {end:?}: {err}"));
        }
    }
}

proptest! {
    #[test]
    fn digamma_against_mpfr(
        x in dyadic_with(2, 0..=0),
        top in -30i64..=12,
        prec in prop_oneof![2u64..120, 120u64..600],
    ) {
        let x = x.to_ours();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let b = Ball::exact(x.clone());
        // MPFR's digamma is NaN at the poles, where ours is indeterminate.
        if b.digamma(64) == Ball::INDETERMINATE {
            prop_assert!(x.is_negative());
            return Ok(());
        }
        check(&format!("ψ({x:?})"), &b.digamma(prec), prec, &digamma_bracket(&x, prec))
            .map_err(TestCaseError::fail)?;
    }

    /// `ψ^(m)(x + 1) − ψ^(m)(x)` holds `(−1)^m m!/x^(m+1)`.
    #[test]
    fn recurrence(
        x in dyadic_with(1, 0..=0),
        top in -6i64..=8,
        m in 0u64..12,
        prec in 2u64..300,
    ) {
        let x = x.to_ours();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let b = Ball::exact(x);
        let one = Ball::from_integer(&Integer::ONE);
        let (now, next) = (b.polygamma(m, prec), b.add(&one, prec + 64).polygamma(m, prec));
        if !now.is_finite() || !next.is_finite() {
            return Ok(());
        }
        let q = prec + 64;
        let i = i64::try_from(m + 1).unwrap();
        let factorial = Ball::from_integer(&Integer::from(crate::Natural::factorial(m)));
        let mut step = factorial.div(&b.pow(i, q), q);
        if m % 2 == 1 {
            step = step.neg();
        }
        let difference = next.sub(&now, q);
        prop_assert!(difference.overlaps(&step), "ψ^({})({:?}): {:?} against {:?}", m, b, difference, step);
    }
}
