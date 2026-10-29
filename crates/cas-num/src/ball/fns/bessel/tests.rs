//! Tests of the Bessel functions of balls (CAS-378, CAS-379): `J_n` and
//! `Y_n` against MPFR's `jn` and `yn` across both methods; `I_n` and `K_n`
//! against 50-digit values of their integrals and by their Wronskian; the
//! recurrences, the symmetries, the values at 0, the domains, radii, and
//! huge arguments.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::ops::Pow;

use super::MAX_ORDER;
use crate::Integer;
use crate::ball::{Ball, Float, Mag};
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

/// The bracket of a value given to 50 digits: within `10^-48` of it.
fn decimal(text: &str) -> Bracket {
    let (mantissa, exponent) = text.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    let (whole, fraction) = mantissa.split_once('.').unwrap();
    let digits: rug::Integer = format!("{whole}{fraction}").parse().unwrap();
    let scale = exponent - i32::try_from(fraction.len()).unwrap();
    let ten = rug::Rational::from(10);
    let value = rug::Rational::from(digits) * ten.clone().pow(scale);
    let slack = value.clone().abs() * ten.pow(-48);
    Bracket {
        lo: rug::Rational::from(&value - &slack),
        hi: value + slack,
    }
}

/// `J_n` and `Y_n` against MPFR, on both sides of the threshold of the
/// expansions, for orders below and above the argument.
#[test]
fn against_mpfr() {
    let mut wrong = Vec::new();
    for prec in [53u64, 128, 300] {
        let p = u32::try_from(prec).unwrap();
        for x in [0.5, 1.0, 3.75, 10.0, 30.5, 52.5, 100.0, 1000.0, -7.25] {
            let xf = float(x);
            let xr = MpFloat::with_val(10_000, &xf.to_rational().to_rug());
            for n in [0i32, 1, 2, 5, 17, 100] {
                let j = exact(x).bessel_j(i64::from(n), prec);
                let oracle = bracket(p, |bits, r| {
                    MpFloat::with_val_round(bits, xr.jn_ref(n), r).0
                });
                if let Err(e) = check(&format!("J_{n}({x})"), &j, prec, &oracle) {
                    wrong.push(e);
                }
                if x > 0.0 {
                    let y = exact(x).bessel_y(i64::from(n), prec);
                    let oracle = bracket(p, |bits, r| {
                        MpFloat::with_val_round(bits, xr.yn_ref(n), r).0
                    });
                    if let Err(e) = check(&format!("Y_{n}({x})"), &y, prec, &oracle) {
                        wrong.push(e);
                    }
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `I_n` and `K_n` against `(1/π) ∫_0^π e^(x cos θ) cos(nθ) dθ` and
/// `∫_0^∞ e^(−x cosh t) cosh(nt) dt` (DLMF 10.32.3, 10.32.9), summed by the
/// trapezoidal rule to 88 digits in Python's `Decimal`.
#[test]
fn against_integrals() {
    let mut wrong = Vec::new();
    for prec in [53u64, 128] {
        for (n, x, value) in [
            (
                0,
                2.0,
                "2.27958530233606726743720444081153335328584110278546e+0",
            ),
            (
                3,
                10.0,
                "1.75838071661085323810039006274328871584602427073211e+3",
            ),
            (
                1,
                0.5,
                "2.57894305390896316362479659523209634187743149640795e-1",
            ),
            (
                7,
                20.0,
                "1.25628736861788495660474905039140225747224419751979e+7",
            ),
            (
                0,
                60.0,
                "5.89407705560980116827881744033390473797898302031621e+24",
            ),
            (
                1,
                60.0,
                "5.84475158839046828133517287252083704489777978998393e+24",
            ),
        ] {
            let i = exact(x).bessel_i(n, prec);
            if let Err(e) = check(&format!("I_{n}({x})"), &i, prec, &decimal(value)) {
                wrong.push(e);
            }
        }
        for (n, x, value) in [
            (
                0,
                2.0,
                "1.13893872749533435652719574932481832998326624388809e-1",
            ),
            (
                1,
                0.5,
                "1.65644112000330089369644540317409151153410075946408e+0",
            ),
            (
                5,
                3.0,
                "9.37773602386808030572692957110100586631159423906341e-1",
            ),
            (
                0,
                0.125,
                "2.20786908674497006000761214484982058638524650976281e+0",
            ),
            (
                2,
                10.0,
                "2.15098170069327687306645644239671272492068461808732e-5",
            ),
            (
                0,
                60.0,
                "1.41389784055910780909564629824170120316175653486573e-27",
            ),
            (
                1,
                200.0,
                "1.22874237347298581204459101041646009828967124256215e-88",
            ),
        ] {
            let k = exact(x).bessel_k(n, prec);
            if let Err(e) = check(&format!("K_{n}({x})"), &k, prec, &decimal(value)) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The values at 0, negative orders and arguments, and the domains.
#[test]
fn symmetries_and_domains() {
    let zero = Ball::ZERO;
    assert_eq!(zero.bessel_j(0, 64), Ball::from_integer(&Integer::ONE));
    assert_eq!(zero.bessel_j(3, 64), Ball::ZERO);
    assert_eq!(zero.bessel_i(0, 64), Ball::from_integer(&Integer::ONE));
    assert_eq!(zero.bessel_i(-2, 64), Ball::ZERO);
    for f in [Ball::bessel_y, Ball::bessel_k] {
        assert!(!f(&zero, 0, 64).is_finite());
        assert!(!f(&exact(-1.5), 1, 64).is_finite());
        assert!(
            !f(
                &Ball::new(float(0.5), crate::ball::Mag::from_float(&float(1.0))),
                1,
                64
            )
            .is_finite()
        );
    }
    let x = exact(2.5);
    let minus = exact(-2.5);
    assert_eq!(x.bessel_j(-3, 100), x.bessel_j(3, 100).neg());
    assert_eq!(x.bessel_j(-2, 100), x.bessel_j(2, 100));
    assert_eq!(x.bessel_y(-3, 100), x.bessel_y(3, 100).neg());
    assert_eq!(x.bessel_i(-3, 100), x.bessel_i(3, 100));
    assert_eq!(x.bessel_k(-3, 100), x.bessel_k(3, 100));
    assert_eq!(minus.bessel_j(3, 100), x.bessel_j(3, 100).neg());
    assert_eq!(minus.bessel_i(2, 100), x.bessel_i(2, 100));
    let past = i64::try_from(MAX_ORDER).unwrap() + 1;
    assert!(!x.bessel_j(past, 64).is_finite());
}

/// A ball holds the values at its ends.
#[test]
fn radii() {
    let r = 2f64.powi(-20);
    for x in [0.75, 9.5, 120.0] {
        let ball = Ball::new(float(x), crate::ball::Mag::from_float(&float(r)));
        for n in [0i64, 1, 4] {
            for (name, f) in [
                ("J", Ball::bessel_j as fn(&Ball, i64, u64) -> Ball),
                ("Y", Ball::bessel_y),
                ("I", Ball::bessel_i),
                ("K", Ball::bessel_k),
            ] {
                let whole = f(&ball, n, 64);
                for end in [x - r, x + r] {
                    let at = f(&exact(end), n, 64);
                    assert!(
                        whole.overlaps(&at),
                        "{name}_{n}({x} ± 2^-20): {whole:?} misses {at:?}"
                    );
                }
            }
        }
    }
}

/// Huge and tiny arguments stay finite with the right sizes.
#[test]
fn far_arguments() {
    let huge = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(100)));
    let j = huge.bessel_j(0, 64);
    assert!(j.is_finite() && j.rel_accuracy_bits() >= 60, "{j:?}");
    assert!(*j.mid().exponent() <= Integer::from(-48), "{j:?}");
    let k = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(40))).bessel_k(0, 64);
    assert!(k.is_finite() && !k.mid().is_negative(), "{k:?}");
    let tiny = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(-200)));
    let y = tiny.bessel_y(2, 64);
    assert!(y.is_finite() && y.mid().is_negative(), "{y:?}");
}

/// The fuzz target `bessel_fns`'s timeout (CAS-381): the order −1024 at
/// `x = −1647099/4`, `π·2^17` rounded to 21 bits, where MPFR's `jn` takes
/// its series for minutes. Here Hankel's `J_0` and `J_1` start a forward
/// recurrence: finite with `prec − 3` bits, `J_(−1024) = J_1024`, under the
/// envelope `√(2/(π x))`, and the recurrence across 1022 to 1024 holds.
#[test]
fn largest_order_far_out() {
    let x = exact(-411_774.75);
    let prec = 89;
    let j = x.bessel_j(-1024, prec);
    assert!(j.is_finite() && j.rel_accuracy_bits() >= 86, "{j:?}");
    assert!(j.overlaps(&x.bessel_j(1024, prec)), "{j:?}");
    assert!(j.mid().to_f64().abs() < 1.0 / 400.0, "{j:?}");
    let q = 128;
    let (a, b, c) = (
        x.bessel_j(1022, q),
        x.bessel_j(1023, q),
        x.bessel_j(1024, q),
    );
    let two_n = Ball::from_integer(&Integer::from(2046));
    let (left, right) = (a.add(&c, q), two_n.mul(&b, q).div(&x, q));
    assert!(left.overlaps(&right), "{left:?} against {right:?}");
}

/// The fuzz target `bessel_fns`'s timeout (CAS-430): `J_(−12)` of the
/// ball `0 ± 320748.1142578125`, whose ends MPFR's `jn` took to its series,
/// as its expansion reads a negative order unsigned. At either end
/// `J_(−12) = J_12`, against Python's Hankel expansion at 90 digits
/// (giac's −0.0013942607073), and the ball's value holds it.
#[test]
fn negative_order_far_out() {
    let end = float(320_748.114_257_812_5);
    let reference: crate::Rational = "-139426070730140819687383090870494387896608120448887/\
         100000000000000000000000000000000000000000000000000000"
        .parse()
        .unwrap();
    let reference = Ball::from_rational(&reference, 256);
    for x in [end.clone(), end.neg()] {
        let j = Ball::exact(x).bessel_j(-12, 128);
        assert!(j.is_finite() && j.rel_accuracy_bits() >= 125, "{j:?}");
        assert!(j.overlaps(&reference), "{j:?} against {reference:?}");
    }
    let wide = Ball::new(Float::ZERO, Mag::from_float(&end));
    let j = wide.bessel_j(-12, 232);
    assert!(j.is_finite() && j.overlaps(&reference), "{j:?}");
}

proptest! {
    /// The recurrences: `C_(n−1) + C_(n+1) = (2n/x) C_n` for `J` and `Y`,
    /// `I_(n−1) − I_(n+1) = (2n/x) I_n` and `K_(n+1) − K_(n−1) = (2n/x)
    /// K_n`.
    #[test]
    fn recurrences(
        x in dyadic_with(1, 0..=0),
        top in -4i64..=8,
        n in 1i64..30,
        prec in 2u64..200,
    ) {
        let x = x.to_ours().abs();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let b = Ball::exact(x);
        let q = prec + 64;
        let scale = Ball::from_integer(&Integer::from(2 * n)).div(&b, q);
        for (name, f, sign) in [
            ("J", Ball::bessel_j as fn(&Ball, i64, u64) -> Ball, 1i64),
            ("Y", Ball::bessel_y, 1),
            ("I", Ball::bessel_i, -1),
            ("K", Ball::bessel_k, 0),
        ] {
            let (before, now, after) = (f(&b, n - 1, prec), f(&b, n, prec), f(&b, n + 1, prec));
            if !(before.is_finite() && now.is_finite() && after.is_finite()) {
                continue;
            }
            let (left, right) = match sign {
                1 => (before.add(&after, q), scale.mul(&now, q)),
                -1 => (before.sub(&after, q), scale.mul(&now, q)),
                _ => (after.sub(&before, q), scale.mul(&now, q)),
            };
            prop_assert!(left.overlaps(&right), "{}_{}({:?}): {:?} against {:?}", name, n, b, left, right);
        }
    }

    /// `I_n K_(n+1) + I_(n+1) K_n = 1/x` (DLMF 10.28.2).
    #[test]
    fn wronskian(
        x in dyadic_with(1, 0..=0),
        top in -4i64..=9,
        n in 0i64..20,
        prec in 2u64..200,
    ) {
        let x = x.to_ours().abs();
        if x.is_zero() {
            return Ok(());
        }
        let x = x.mul_2exp(&(&Integer::from(top) - x.exponent()));
        let b = Ball::exact(x);
        let q = prec + 64;
        let w = b.bessel_i(n, prec).mul(&b.bessel_k(n + 1, prec), q)
            .add(&b.bessel_i(n + 1, prec).mul(&b.bessel_k(n, prec), q), q);
        let expected = Ball::from_integer(&Integer::ONE).div(&b, q);
        prop_assert!(w.overlaps(&expected), "ν = {}, x = {:?}: {:?}", n, b, w);
    }
}
