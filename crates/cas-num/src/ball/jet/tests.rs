//! Tests of jets (CAS-389): the coefficients of the elementary functions of
//! `x₀ + s` against their closed forms, the inverse functions at 0 against
//! their series, and identities.

use super::Jet;
use crate::ball::Ball;
use crate::{Integer, Rational};

const PREC: u64 = 128;
const N: usize = 9;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn ball(n: i64, d: i64) -> Ball {
    Ball::from_rational(&q(n, d), PREC)
}

fn factorial(k: usize) -> Ball {
    (1..=k).fold(Ball::from_integer(&Integer::ONE), |acc, j| {
        acc.mul(&Ball::from_integer(&Integer::from(j as u64)), PREC)
    })
}

/// Each coefficient of `found` meets `expected`'s, and has 100 bits.
fn agrees(what: &str, found: &Jet, expected: &[Ball]) {
    assert_eq!(found.len(), expected.len(), "{what}: lengths");
    for (k, (f, e)) in found.coeffs().iter().zip(expected).enumerate() {
        assert!(f.overlaps(e), "{what}: coefficient {k} is {f:?}, not {e:?}");
        assert!(
            f.mid().is_zero() || f.rel_accuracy_bits() >= 100 || f.rad().is_zero(),
            "{what}: coefficient {k} is {f:?}"
        );
    }
}

/// `e^(x₀ + s)`, `log(x₀ + s)`, `sin`, `cos`, `√` and powers at `x₀ = 3/4`.
#[test]
fn closed_forms() {
    let x0 = ball(3, 4);
    let u = Jet::variable(x0.clone(), N);
    let e = x0.exp(PREC);
    agrees(
        "exp",
        &u.exp(PREC),
        &(0..N)
            .map(|k| e.div(&factorial(k), PREC))
            .collect::<Vec<_>>(),
    );
    let mut logs = vec![x0.log(PREC)];
    for k in 1..N {
        // (−1)^(k+1)/(k x₀^k).
        let kk = i64::try_from(k).unwrap();
        let c = ball(if k % 2 == 1 { 1 } else { -1 }, kk).div(&x0.pow(kk, PREC), PREC);
        logs.push(c);
    }
    agrees("log", &u.log(PREC), &logs);
    let (s0, c0) = x0.sin_cos(PREC);
    let (sin, cos) = u.sin_cos(PREC);
    let trig = |k: usize, start: &Ball, other: &Ball, sign: i64| {
        // The k-th derivative of sin is sin(x + kπ/2).
        let value = match k % 4 {
            0 => start.clone(),
            1 => other.mul(&Ball::from_integer(&Integer::from(sign)), PREC),
            2 => start.neg(),
            _ => other.mul(&Ball::from_integer(&Integer::from(-sign)), PREC),
        };
        value.div(&factorial(k), PREC)
    };
    agrees(
        "sin",
        &sin,
        &(0..N).map(|k| trig(k, &s0, &c0, 1)).collect::<Vec<_>>(),
    );
    agrees(
        "cos",
        &cos,
        &(0..N).map(|k| trig(k, &c0, &s0, -1)).collect::<Vec<_>>(),
    );
    // (x₀ + s)^a: C(a, k) x₀^(a − k), for a = 1/2 by sqrt and pow, and
    // a = 5/3 by pow.
    for (name, a, jet) in [
        ("sqrt", q(1, 2), u.sqrt(PREC)),
        ("pow 1/2", q(1, 2), u.pow(&ball(1, 2), PREC)),
        ("pow 5/3", q(5, 3), u.pow(&ball(5, 3), PREC)),
    ] {
        let ab = Ball::from_rational(&a, PREC);
        let lead = x0.log(PREC).mul(&ab, PREC).exp(PREC);
        let mut expected = Vec::with_capacity(N);
        let mut binomial = Ball::from_integer(&Integer::ONE);
        for k in 0..N {
            if k > 0 {
                let kk = i64::try_from(k).unwrap();
                let top = ab.sub(&Ball::from_integer(&Integer::from(kk - 1)), PREC);
                binomial = binomial
                    .mul(&top, PREC)
                    .div(&Ball::from_integer(&Integer::from(kk)), PREC);
            }
            let kk = i64::try_from(k).unwrap();
            expected.push(binomial.mul(&lead, PREC).div(&x0.pow(kk, PREC), PREC));
        }
        agrees(name, &jet, &expected);
    }
}

/// The inverse functions at 0: their Maclaurin series.
#[test]
fn inverse_functions() {
    let s = Jet::variable(Ball::ZERO, 8);
    let zero = Ball::ZERO;
    let series = |c: [Rational; 4]| -> Vec<Ball> {
        // c₁ s + c₃ s³ + c₅ s⁵ + c₇ s⁷.
        let mut out = vec![zero.clone(); 8];
        for (i, ci) in c.iter().enumerate() {
            out[2 * i + 1] = Ball::from_rational(ci, PREC);
        }
        out
    };
    agrees(
        "atan",
        &s.atan(PREC),
        &series([q(1, 1), q(-1, 3), q(1, 5), q(-1, 7)]),
    );
    agrees(
        "atanh",
        &s.atanh(PREC),
        &series([q(1, 1), q(1, 3), q(1, 5), q(1, 7)]),
    );
    agrees(
        "asin",
        &s.asin(PREC),
        &series([q(1, 1), q(1, 6), q(3, 40), q(5, 112)]),
    );
    agrees(
        "asinh",
        &s.asinh(PREC),
        &series([q(1, 1), q(-1, 6), q(3, 40), q(-5, 112)]),
    );
    // tan s = s + s³/3 + 2s⁵/15 + 17 s⁷/315, tanh with alternating signs.
    agrees(
        "tan",
        &s.tan(PREC),
        &series([q(1, 1), q(1, 3), q(2, 15), q(17, 315)]),
    );
    agrees(
        "tanh",
        &s.tanh(PREC),
        &series([q(1, 1), q(-1, 3), q(2, 15), q(-17, 315)]),
    );
    // e^(sin s) = 1 + s + s²/2 − s⁴/8 − s⁵/15 − s⁶/240 + s⁷/90.
    let expected: Vec<Ball> = [
        q(1, 1),
        q(1, 1),
        q(1, 2),
        q(0, 1),
        q(-1, 8),
        q(-1, 15),
        q(-1, 240),
        q(1, 90),
    ]
    .iter()
    .map(|c| Ball::from_rational(c, PREC))
    .collect();
    agrees("exp sin", &s.sin_cos(PREC).0.exp(PREC), &expected);
}

/// Identities at `u = 2/5 + 3s − s² + s⁴/7`: `exp(log u) = u`, `sin² + cos²
/// = 1`, `cosh² − sinh² = 1`, `(√u)² = u`, `u³` by `powi` and products,
/// `u^(−2) u² = 1`, and `asin + acos = π/2` at `u/2`.
#[test]
fn identities() {
    let u = Jet::new(vec![
        ball(2, 5),
        ball(3, 1),
        ball(-1, 1),
        Ball::ZERO,
        ball(1, 7),
        Ball::ZERO,
    ]);
    let one = Jet::constant(Ball::from_integer(&Integer::ONE), u.len());
    let same = |what: &str, a: &Jet, b: &Jet| agrees(what, a, b.coeffs());
    same("exp log", &u.log(PREC).exp(PREC), &u);
    let (s, c) = u.sin_cos(PREC);
    same(
        "sin² + cos²",
        &s.mul(&s, PREC).add(&c.mul(&c, PREC), PREC),
        &one,
    );
    let (sh, ch) = u.sinh_cosh(PREC);
    same(
        "cosh² − sinh²",
        &ch.mul(&ch, PREC).sub(&sh.mul(&sh, PREC), PREC),
        &one,
    );
    let root = u.sqrt(PREC);
    same("√u²", &root.mul(&root, PREC), &u);
    same("u³", &u.powi(3, PREC), &u.mul(&u, PREC).mul(&u, PREC));
    same(
        "u^−2 u²",
        &u.powi(-2, PREC).mul(&u.powi(2, PREC), PREC),
        &one,
    );
    let half = u.scale(&ball(1, 2), PREC);
    let half_pi = Jet::constant(Ball::const_pi(PREC).mul_2exp(&-Integer::ONE), u.len());
    same(
        "asin + acos",
        &half.asin(PREC).add(&half.acos(PREC), PREC),
        &half_pi,
    );
    // acosh(cosh u) = u for u₀ > 0, and the derivative and integral.
    same("acosh cosh", &ch.acosh(PREC), &u);
    let back = u.derivative(PREC).integral(u.coeffs()[0].clone(), PREC);
    same("∫ u′", &back, &u);
}
