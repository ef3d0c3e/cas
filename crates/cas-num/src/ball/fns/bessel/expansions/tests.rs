//! Tests of the asymptotic expansions (CAS-377): half-integer orders, where
//! the expansions end and are the closed forms; Hankel's `J_0`, `J_1`,
//! `Y_0` and `Y_1` against MPFR; `K_ν` and `I_n` against 50-digit values of
//! their integrals; the Wronskian of `I_ν` and `K_ν`.

use rug::Float as MpFloat;
use rug::ops::Pow;

use super::{Order, hankel, i_expansion, k_expansion, threshold};
use crate::ball::{Ball, Float};
use crate::testkit::ToRug;
use crate::testkit::float::{Bracket, bracket, check_encloses};
use crate::{Integer, Rational};

fn order(p: i64, q: i64) -> Order {
    Order::new(Rational::new(Integer::from(p), Integer::from(q)))
}

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

fn exact(x: f64) -> Ball {
    Ball::exact(float(x))
}

/// `b` meets the bracket and, unless it is 0, has `prec - 8` bits.
fn check(what: &str, b: &Ball, prec: u64, oracle: &Bracket) -> Result<(), String> {
    if !b.is_finite() {
        return Err(format!("{what} at {prec}: indeterminate"));
    }
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, oracle).map_err(|e| format!("{what} at {prec}: {e}"))?;
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let want = prec as i64 - 8;
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

/// Where the expansions end: `ν = 1/2` and `3/2` give the closed forms,
/// `√(2/(πx)) sin x` and the like.
#[test]
fn half_integer_orders() {
    let q = 256;
    for x in [60.0, 123.25, 1e6] {
        let xf = float(x);
        let xb = Ball::exact(xf.clone());
        let pi = Ball::const_pi(q);
        let one = Ball::from_integer(&Integer::ONE);
        let (sin, cos) = xb.sin_cos(q);
        let root = Ball::from_integer(&Integer::from(2))
            .div(&pi.mul(&xb, q), q)
            .sqrt(q);
        let decay = pi
            .div(&xb.mul_2exp(&Integer::ONE), q)
            .sqrt(q)
            .mul(&xb.neg().exp(q), q);
        let (sinh, cosh) = {
            let (up, down) = (xb.exp(q), xb.neg().exp(q));
            (
                up.sub(&down, q).mul_2exp(&-Integer::ONE),
                up.add(&down, q).mul_2exp(&-Integer::ONE),
            )
        };
        let inverse = one.div(&xb, q);
        let cases = [
            ((1, 2), root.mul(&sin, q), root.mul(&cos, q).neg()),
            (
                (3, 2),
                root.mul(&sin.mul(&inverse, q).sub(&cos, q), q),
                root.mul(&cos.mul(&inverse, q).add(&sin, q), q).neg(),
            ),
        ];
        for ((p, d), j, y) in cases {
            let (found_j, found_y) = hankel(&order(p, d), &Ball::exact(xf.clone()), 128).unwrap();
            assert!(
                found_j.overlaps(&j),
                "J_{p}/{d}({x}): {found_j:?} against {j:?}"
            );
            assert!(
                found_y.overlaps(&y),
                "Y_{p}/{d}({x}): {found_y:?} against {y:?}"
            );
        }
        let k_half = k_expansion(&order(1, 2), &Ball::exact(xf.clone()), 128).unwrap();
        assert!(k_half.overlaps(&decay), "K_1/2({x}): {k_half:?}");
        let k_three = k_expansion(&order(3, 2), &Ball::exact(xf.clone()), 128).unwrap();
        let expected = decay.mul(&one.add(&inverse, q), q);
        assert!(k_three.overlaps(&expected), "K_3/2({x}): {k_three:?}");
        if x < 1e5 {
            let i_half = i_expansion(&order(1, 2), &Ball::exact(xf.clone()), 128).unwrap();
            assert!(
                i_half.overlaps(&root.mul(&sinh, q)),
                "I_1/2({x}): {i_half:?}"
            );
            let i_three = i_expansion(&order(3, 2), &Ball::exact(xf.clone()), 128).unwrap();
            let expected = root.mul(&cosh.sub(&sinh.mul(&inverse, q), q), q);
            assert!(i_three.overlaps(&expected), "I_3/2({x}): {i_three:?}");
        }
    }
}

/// `J_0`, `J_1`, `Y_0` and `Y_1` past the threshold, against MPFR's.
#[test]
fn hankel_against_mpfr() {
    let mut wrong = Vec::new();
    for prec in [53u64, 128, 300] {
        let start = threshold(prec).ceil();
        for x in [
            start,
            start + 0.5,
            2.0 * start + 0.25,
            1000.0,
            12345.5,
            1_073_741_824.0,
        ] {
            let xf = float(x);
            let xr = MpFloat::with_val(10_000, &xf.to_rational().to_rug());
            for n in [0i32, 1] {
                let Some((j, y)) = hankel(&order(i64::from(n), 1), &Ball::exact(xf.clone()), prec)
                else {
                    wrong.push(format!("n = {n}, x = {x} at {prec}: no expansion"));
                    continue;
                };
                let p = u32::try_from(prec).unwrap();
                let j_oracle = bracket(p, |bits, r| {
                    MpFloat::with_val_round(bits, xr.jn_ref(n), r).0
                });
                let y_oracle = bracket(p, |bits, r| {
                    MpFloat::with_val_round(bits, xr.yn_ref(n), r).0
                });
                if let Err(e) = check(&format!("J_{n}({x})"), &j, prec, &j_oracle) {
                    wrong.push(e);
                }
                if let Err(e) = check(&format!("Y_{n}({x})"), &y, prec, &y_oracle) {
                    wrong.push(e);
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `K_ν` at 60 and 200 and `I_n` at 60, against `∫_0^∞ e^(−x cosh t)
/// cosh(νt) dt` and `(1/π) ∫_0^π e^(x cos θ) cos(nθ) dθ` (DLMF 10.32.9,
/// 10.32.3), summed by the trapezoidal rule to 88 digits in Python's
/// `Decimal`.
#[test]
fn against_integrals() {
    let k_values = [
        (
            (0, 1),
            60.0,
            "1.41389784055910780909564629824170120316175653486573e-27",
        ),
        (
            (0, 1),
            200.0,
            "1.22568197977653345166005408750744722858043232343307e-88",
        ),
        (
            (1, 1),
            60.0,
            "1.42563202651710432321438861042921740915123025597281e-27",
        ),
        (
            (1, 1),
            200.0,
            "1.22874237347298581204459101041646009828967124256215e-88",
        ),
        (
            (1, 3),
            60.0,
            "1.41519688056236612385696865176189271968750311017083e-27",
        ),
        (
            (1, 3),
            200.0,
            "1.22602164737766005300228073145356630524635083230093e-88",
        ),
        (
            (2, 3),
            60.0,
            "1.41910112741720786592017985995118658378458333907760e-27",
        ),
        (
            (2, 3),
            200.0,
            "1.22704121413571054384926893048938000591278851527090e-88",
        ),
        (
            (4, 3),
            60.0,
            "1.43482553720123415618525728941520761400333337363505e-27",
        ),
        (
            (4, 3),
            200.0,
            "1.23112795296030274402594319959422522693027635137857e-88",
        ),
    ];
    let mut wrong = Vec::new();
    for ((p, q), x, value) in k_values {
        let k = k_expansion(&order(p, q), &exact(x), 128).unwrap();
        if let Err(e) = check(&format!("K_{p}/{q}({x})"), &k, 128, &decimal(value)) {
            wrong.push(e);
        }
    }
    for (n, value) in [
        (
            0,
            "5.89407705560980116827881744033390473797898302031621e+24",
        ),
        (
            1,
            "5.84475158839046828133517287252083704489777978998393e+24",
        ),
    ] {
        let i = i_expansion(&order(n, 1), &exact(60.0), 128).unwrap();
        if let Err(e) = check(&format!("I_{n}(60)"), &i, 128, &decimal(value)) {
            wrong.push(e);
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `I_ν K_(ν+1) + I_(ν+1) K_ν = 1/x` (DLMF 10.28.2), for integer and
/// non-integer orders.
#[test]
fn wronskian() {
    let q = 200;
    for (p, d) in [(0, 1), (1, 3), (2, 3), (1, 1), (5, 2), (7, 1)] {
        let next = (p + d, d);
        for x in [80.0, 150.5, 4000.0] {
            let xf = float(x);
            let (i0, i1) = (
                i_expansion(&order(p, d), &Ball::exact(xf.clone()), 128).unwrap(),
                i_expansion(&order(next.0, next.1), &Ball::exact(xf.clone()), 128).unwrap(),
            );
            let (k0, k1) = (
                k_expansion(&order(p, d), &Ball::exact(xf.clone()), 128).unwrap(),
                k_expansion(&order(next.0, next.1), &Ball::exact(xf.clone()), 128).unwrap(),
            );
            let w = i0.mul(&k1, q).add(&i1.mul(&k0, q), q);
            let expected = Ball::from_integer(&Integer::ONE).div(&Ball::exact(xf.clone()), q);
            assert!(w.overlaps(&expected), "ν = {p}/{d}, x = {x}: {w:?}");
            assert!(w.rel_accuracy_bits() >= 110, "ν = {p}/{d}, x = {x}: {w:?}");
        }
    }
}

/// Below the threshold the terms grow before they reach the precision.
#[test]
fn small_arguments() {
    assert!(hankel(&order(0, 1), &exact(5.0), 128).is_none());
    assert!(k_expansion(&order(0, 1), &exact(5.0), 128).is_none());
    assert!(i_expansion(&order(0, 1), &exact(5.0), 128).is_none());
}
