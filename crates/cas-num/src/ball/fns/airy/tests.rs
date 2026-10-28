//! Tests of Airy's functions (CAS-380): `Ai` against MPFR's `ai` on both
//! sides of the switch to Bessel functions; `Ai′`, `Bi` and `Bi′` through
//! the Wronskian `Ai Bi′ − Ai′ Bi = 1/π`; the values at 0; radii; the
//! exponent bound.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::ops::Pow;

use crate::Integer;
use crate::ball::{Airy, Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{Bracket, bracket, check_encloses, dyadic_with};

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

fn exact(x: f64) -> Ball {
    Ball::exact(float(x))
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

/// `Ai` against MPFR's, from the series to the expansions.
#[test]
fn against_mpfr() {
    let mut wrong = Vec::new();
    for prec in [53u64, 128, 300] {
        let p = u32::try_from(prec).unwrap();
        for x in [-200.0, -20.0, -5.5, -1.0, 0.0, 0.5, 2.0, 7.25, 30.0, 100.0] {
            let xf = float(x);
            let xr = MpFloat::with_val(10_000, &xf.to_rational().to_rug());
            let ai = exact(x).airy(Airy::Ai, prec);
            let oracle = bracket(p, |bits, r| MpFloat::with_val_round(bits, xr.ai_ref(), r).0);
            if let Err(e) = check(&format!("Ai({x})"), &ai, prec, &oracle) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `Ai(0) = 3^(−2/3)/Γ(2/3)`, `Ai′(0) = −3^(−1/3)/Γ(1/3)`, `Bi(0) = √3
/// Ai(0)` and `Bi′(0) = −√3 Ai′(0)`, in MPFR at 2000 bits.
#[test]
fn at_zero() {
    let p = 2000;
    let gamma = |a: u32| MpFloat::with_val(p, rug::Rational::from((a, 3))).gamma();
    let three = MpFloat::with_val(p, 3);
    let c1 = three
        .clone()
        .pow(MpFloat::with_val(p, rug::Rational::from((-2, 3))))
        / gamma(2);
    let c2 = three
        .clone()
        .pow(MpFloat::with_val(p, rug::Rational::from((-1, 3))))
        / gamma(1);
    let root3 = three.sqrt();
    let values = [c1.clone(), -c2.clone(), root3.clone() * &c1, root3 * &c2];
    for (which, value) in [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime]
        .into_iter()
        .zip(values)
    {
        let found = Ball::ZERO.airy(which, 200);
        let exact = value.to_rational().unwrap();
        let slack = exact.clone().abs() >> 1900u32;
        let oracle = Bracket {
            lo: rug::Rational::from(&exact - &slack),
            hi: rug::Rational::from(&exact + &slack),
        };
        check(&format!("{which:?}(0)"), &found, 200, &oracle).unwrap();
    }
}

/// `Ai Bi′ − Ai′ Bi = 1/π` (DLMF 9.2.7), on both sides of the switch.
#[test]
fn wronskian() {
    for prec in [53u64, 128, 300] {
        let q = prec + 64;
        for x in [-150.5, -12.0, -0.75, 0.0, 1.5, 9.0, 25.0, 80.0] {
            let b = exact(x);
            let [ai, aip, bi, bip] =
                [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime].map(|which| b.airy(which, prec));
            let w = ai.mul(&bip, q).sub(&aip.mul(&bi, q), q);
            let expected = Ball::from_integer(&Integer::ONE).div(&Ball::const_pi(q), q);
            assert!(w.overlaps(&expected), "x = {x} at {prec}: {w:?}");
            #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
            let want = prec as i64 - 6;
            assert!(w.rel_accuracy_bits() >= want, "x = {x} at {prec}: {w:?}");
        }
    }
}

/// A ball holds the values at its ends; past the exponent bound,
/// indeterminate.
/// The fuzz target `bessel_fns`'s timeout (CAS-420): the ball of midpoint
/// 0 and radius `758134843·2^−10`, about `7.4·10⁵`, at 49 bits. The four
/// functions are finite, their balls hold the values at points inside, and
/// Airy's Wronskian holds over the ball.
#[test]
fn wide_ball_meeting_zero() {
    let radius = Mag::from_float(&Float::from_integer_exp(
        &Integer::from(758_134_843),
        &Integer::from(-10),
    ));
    let x = Ball::new(Float::ZERO, radius);
    let prec = 49;
    let values = [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime].map(|w| x.airy(w, prec));
    for (w, y) in [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime]
        .iter()
        .zip(&values)
    {
        assert!(y.is_finite(), "{w:?}: {y:?}");
        for p in [-700_000.0, -1.5, 0.0, 2.0, 700_000.0] {
            let inside = exact(p).airy(*w, 64);
            assert!(
                y.overlaps(&inside),
                "{w:?} at {p}: {inside:?} outside {y:?}"
            );
        }
    }
    let q = 128;
    let [ai, aip, bi, bip] = values;
    let w = ai.mul(&bip, q).sub(&aip.mul(&bi, q), q);
    let expected = Ball::from_integer(&Integer::ONE).div(&Ball::const_pi(q), q);
    assert!(w.overlaps(&expected), "the Wronskian is {w:?}");
}

#[test]
fn radii_and_bounds() {
    let r = 2f64.powi(-20);
    for x in [-9.5, 0.25, 3.0, 40.0] {
        let ball = Ball::new(float(x), Mag::from_float(&float(r)));
        for which in [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime] {
            let whole = ball.airy(which, 64);
            for end in [x - r, x + r] {
                let at = exact(end).airy(which, 64);
                assert!(
                    whole.overlaps(&at),
                    "{which:?}({x} ± 2^-20): {whole:?} misses {at:?}"
                );
            }
        }
    }
    let far = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(5000)));
    assert!(!far.airy(Airy::Ai, 64).is_finite());
}

proptest! {
    /// `Ai Bi′ − Ai′ Bi = 1/π` at random points.
    #[test]
    fn wronskian_anywhere(
        x in dyadic_with(1, 0..=0),
        top in -6i64..=8,
        negative in any::<bool>(),
        prec in 2u64..160,
    ) {
        let x = x.to_ours().abs();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let b = Ball::exact(if negative { x.neg() } else { x });
        let q = prec + 64;
        let [ai, aip, bi, bip] = [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime]
            .map(|which| b.airy(which, prec));
        let w = ai.mul(&bip, q).sub(&aip.mul(&bi, q), q);
        let expected = Ball::from_integer(&Integer::ONE).div(&Ball::const_pi(q), q);
        prop_assert!(w.overlaps(&expected), "x = {:?}: {:?}", b, w);
    }
}
