//! Tests of the integrands on complex balls (CAS-323): values, the poles
//! and cuts kept away in the analytic mode and covered in the real one,
//! switches continued from a segment, and what does not compile.

#![allow(
    clippy::needless_pass_by_value,
    reason = "the helpers take the rationals the cases build"
)]

use cas_num::ball::{Ball, ComplexBall, Float, Mag};
use cas_num::{Integer, Rational};

use super::{Mode, Node};
use crate::{Context, Expr};

const PREC: u64 = 128;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// The ball `re + i im ± r` in each part, `r = 2^-k`.
fn disk(re: Rational, im: Rational, k: i64) -> ComplexBall {
    let radius = Mag::from_integer(&Integer::from(1)).mul_2exp(&Integer::from(-k));
    ComplexBall::new(
        Ball::from_rational(&re, PREC).add_error(&radius),
        Ball::from_rational(&im, PREC).add_error(&radius),
    )
}

fn point(re: Rational) -> ComplexBall {
    ComplexBall::from_real(Ball::from_rational(&re, PREC))
}

fn compile(f: &str) -> Option<Node> {
    let ctx = Context::new();
    let f = ctx.eval(&read(f)).unwrap();
    Node::compile(&ctx, &f, &Expr::symbol("x"), PREC)
}

fn at(f: &str, z: &ComplexBall, mode: Mode) -> Option<ComplexBall> {
    compile(f).expect("compiles").eval(z, mode, PREC)
}

#[test]
fn values() {
    let v = at("Plus[1, Power[x, 2]]", &point(q(2, 1)), Mode::Analytic).unwrap();
    assert!(v.contains_complex(&cas_num::Complex::from(q(5, 1))));
    // e^(iπ) = −1.
    let pi = ComplexBall::from_real(Ball::const_pi(PREC));
    let v = at("Exp[Times[I, x]]", &pi, Mode::Analytic).unwrap();
    assert!(v.contains_complex(&cas_num::Complex::from(q(-1, 1))));
    let v = at(
        "Times[Power[x, 1/2], Log[x]]",
        &point(q(4, 1)),
        Mode::Analytic,
    )
    .unwrap();
    assert!(v.is_finite() && v.is_real());
}

/// Poles and cuts: kept away in the analytic mode, covered in the real
/// one where the value stays finite.
#[test]
fn poles_and_cuts() {
    let near = |re: Rational, im: Rational| disk(re, im, 3);
    for (f, z, analytic, real) in [
        ("Log[x]", near(q(2, 1), q(0, 1)), true, true),
        ("Log[x]", near(q(-1, 1), q(0, 1)), false, true),
        ("Log[x]", near(q(-1, 1), q(1, 1)), true, true),
        ("Power[x, 1/2]", near(q(-1, 1), q(0, 1)), false, true),
        ("Power[x, -1]", near(q(0, 1), q(0, 1)), false, false),
        (
            "Tan[x]",
            ComplexBall::from_real(Ball::const_pi(PREC).mul_2exp(&Integer::from(-1))),
            false,
            false,
        ),
        ("ArcTan[x]", near(q(0, 1), q(2, 1)), false, true),
        ("ArcTan[x]", near(q(1, 2), q(0, 1)), true, true),
        ("ArcSin[x]", near(q(2, 1), q(0, 1)), false, true),
        ("ArcSin[x]", near(q(1, 2), q(1, 1)), true, true),
        ("ArcCosh[x]", near(q(1, 2), q(0, 1)), false, true),
        ("Power[2, x]", near(q(1, 1), q(1, 1)), true, true),
        ("Power[x, x]", near(q(-1, 1), q(0, 1)), false, true),
    ] {
        assert_eq!(
            at(f, &z, Mode::Analytic).is_some(),
            analytic,
            "{f} analytic"
        );
        assert_eq!(at(f, &z, Mode::Real).is_some(), real, "{f} real");
    }
}

/// Switches: nothing in the analytic mode, their values in the real one,
/// continued from a segment where their argument's sign is known.
#[test]
fn switches() {
    let segment = |lo: i64, hi: i64| {
        let m = Ball::from_rational(&q(lo + hi, 2), PREC);
        let r = Mag::from_integer(&Integer::from(hi - lo)).mul_2exp(&Integer::from(-1));
        ComplexBall::from_real(m.add_error(&r))
    };
    let abs = compile("Abs[Plus[x, -1/2]]").unwrap();
    assert!(abs.has_switch());
    assert!(abs.eval(&point(q(1, 1)), Mode::Analytic, PREC).is_none());
    let real = abs.eval(&segment(-1, 2), Mode::Real, PREC).unwrap();
    assert!(real.re().contains_rational(&q(0, 1)) && real.re().contains_rational(&q(3, 2)));
    // Over [1, 2], x − 1/2 > 0: |x − 1/2| is x − 1/2; over [−2, 0] it is
    // 1/2 − x; over [0, 1] it is not continued.
    let right = abs.resolved(&segment(1, 2), PREC).unwrap();
    assert!(!right.has_switch());
    let v = right.eval(&point(q(3, 1)), Mode::Analytic, PREC).unwrap();
    assert!(v.re().contains_rational(&q(5, 2)));
    let left = abs.resolved(&segment(-2, 0), PREC).unwrap();
    let v = left.eval(&point(q(-3, 1)), Mode::Analytic, PREC).unwrap();
    assert!(v.re().contains_rational(&q(7, 2)));
    assert!(abs.resolved(&segment(0, 1), PREC).is_none());
    let sign = compile("Sign[x]").unwrap();
    let v = sign.eval(&segment(-1, 1), Mode::Real, PREC).unwrap();
    assert!(v.re().contains_rational(&q(-1, 1)) && v.re().contains_rational(&q(1, 1)));
    let step = compile("HeavisideTheta[x]")
        .unwrap()
        .resolved(&segment(1, 3), PREC)
        .unwrap();
    let v = step.eval(&point(q(-5, 1)), Mode::Analytic, PREC).unwrap();
    assert!(v.re().contains_rational(&q(1, 1)));
    let _ = Float::ZERO;
}

/// What does not compile: an unknown head, a parameter.
#[test]
fn what_does_not_compile() {
    for f in ["f[x]", "Times[a, x]", "Gamma[x]", "Power[x, 100000000]"] {
        assert!(compile(f).is_none(), "{f}");
    }
}
