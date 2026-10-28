//! Jets of balls (CAS-DOC-35 §5, CAS-390), on random jets `u` with their
//! constant in each function's domain:
//! - the constant of `f(u)` is the ball function at `u₀`;
//! - `f(u)′ = f′(u) u′`, coefficient by coefficient, the identity the
//!   recurrences come from, `f′` written with other jets: `exp(u)`, `1/u`,
//!   `cos u`, `1/(1 + u²)`, …;
//! - two precisions agree.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::ball::{Ball, Jet};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = run(&mut Unstructured::new(data));
});

fn ball(u: &mut Unstructured<'_>, low: i64, high: i64) -> Result<Ball> {
    let n = u.int_in_range(low..=high)?;
    let d = u.int_in_range(1..=16i64)?;
    Ok(Ball::from_rational(
        &Rational::new(Integer::from(n), Integer::from(d)),
        256,
    ))
}

/// The coefficients of two jets meet.
fn meet(what: &str, a: &Jet, b: &Jet) {
    for (k, (x, y)) in a.coeffs().iter().zip(b.coeffs()).enumerate() {
        assert!(
            x.overlaps(y),
            "{what}: coefficient {k}: {x:?} against {y:?}"
        );
    }
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let prec = u.int_in_range(16u64..=400)?;
    let n = u.int_in_range(1usize..=12)?;
    let which = u.int_in_range(0u8..=14)?;
    // The constant in the domain: positive, inside (−1, 1), above 1, or
    // anywhere.
    let u0 = match which {
        1 | 4 | 5 => ball(u, 1, 64)?,
        8 | 11 => ball(u, -15, 15)?.mul(
            &Ball::from_rational(&Rational::new(Integer::ONE, Integer::from(16)), 256),
            256,
        ),
        12 => ball(u, 17, 64)?,
        _ => ball(u, -64, 64)?,
    };
    let mut coeffs = vec![u0];
    for _ in 1..n {
        coeffs.push(ball(u, -16, 16)?);
    }
    let v = Jet::new(coeffs);
    let one = Jet::constant(Ball::from_integer(&Integer::ONE), n);
    let q = prec + 64;
    let (f, derivative, value): (fn(&Jet, u64) -> Jet, Jet, Ball) = match which {
        0 => (Jet::exp, v.exp(q), v.coeffs()[0].exp(q)),
        1 => (Jet::log, one.div(&v, q), v.coeffs()[0].log(q)),
        2 => (
            |j: &Jet, p| j.sin_cos(p).0,
            v.sin_cos(q).1,
            v.coeffs()[0].sin(q),
        ),
        3 => (
            |j: &Jet, p| j.sin_cos(p).1,
            v.sin_cos(q).0.neg(),
            v.coeffs()[0].cos(q),
        ),
        4 => (
            Jet::sqrt,
            v.sqrt(q)
                .scale(&Ball::from_integer(&Integer::from(2)), q)
                .powi(-1, q),
            v.coeffs()[0].sqrt(q),
        ),
        5 => {
            let a = Ball::from_rational(&Rational::new(Integer::from(5), Integer::from(3)), 256);
            (
                |j: &Jet, p| {
                    j.pow(
                        &Ball::from_rational(
                            &Rational::new(Integer::from(5), Integer::from(3)),
                            256,
                        ),
                        p,
                    )
                },
                v.pow(&a.sub(&Ball::from_integer(&Integer::ONE), q), q)
                    .scale(&a, q),
                v.coeffs()[0].log(q).mul(&a, q).exp(q),
            )
        }
        6 => (
            |j: &Jet, p| j.powi(3, p),
            v.powi(2, q)
                .scale(&Ball::from_integer(&Integer::from(3)), q),
            v.coeffs()[0].pow(3, q),
        ),
        7 => (
            Jet::atan,
            one.div(&v.mul(&v, q).add(&one, q), q),
            v.coeffs()[0].atan(q),
        ),
        8 => (
            Jet::asin,
            one.div(&one.sub(&v.mul(&v, q), q).sqrt(q), q),
            v.coeffs()[0].asin(q),
        ),
        9 => (
            |j: &Jet, p| j.sinh_cosh(p).0,
            v.sinh_cosh(q).1,
            v.coeffs()[0].sinh(q),
        ),
        10 => (
            Jet::asinh,
            one.div(&v.mul(&v, q).add(&one, q).sqrt(q), q),
            v.coeffs()[0].asinh(q),
        ),
        11 => (
            Jet::atanh,
            one.div(&one.sub(&v.mul(&v, q), q), q),
            v.coeffs()[0].atanh(q),
        ),
        12 => (
            Jet::acosh,
            one.div(&v.mul(&v, q).sub(&one, q).sqrt(q), q),
            v.coeffs()[0].acosh(q),
        ),
        13 => {
            let t = v.tan(q);
            (Jet::tan, t.mul(&t, q).add(&one, q), v.coeffs()[0].tan(q))
        }
        _ => {
            let t = v.tanh(q);
            (Jet::tanh, one.sub(&t.mul(&t, q), q), v.coeffs()[0].tanh(q))
        }
    };
    let y = f(&v, prec);
    assert!(
        y.coeffs()[0].overlaps(&value),
        "f #{which}: {:?} against {value:?}",
        y.coeffs()[0]
    );
    // f(u)′ = f′(u) u′.
    let left = y.derivative(q);
    let right = derivative.mul(&v.derivative(q), q);
    meet(&format!("f #{which}: f(u)′"), &left, &right);
    meet(
        &format!("f #{which}: two precisions"),
        &y,
        &f(&v, prec + 64),
    );
    Ok(())
}
