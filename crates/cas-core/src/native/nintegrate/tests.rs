//! Tests of `NIntegrate` (CAS-DOC-26 §4, CAS-324): the ball contains the
//! integral, known to 44 digits or exactly, and is accurate; switches,
//! an end where the integrand is not smooth, complex values, the forms of
//! the call, and what stays.

use cas_num::ball::{Ball, ComplexBall, Mag};
use cas_num::{Integer, Rational};

use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The ball of a decimal `0.d…`, its last digit's half unit around it.
fn reference(digits: &str) -> Ball {
    let (sign, digits) = digits.strip_prefix('-').map_or((1, digits), |d| (-1, d));
    let fraction = digits.strip_prefix("0.").expect("a decimal below 1");
    let numer: Integer = fraction.parse().expect("digits");
    let places = u64::try_from(fraction.len()).unwrap();
    let scale = Integer::from(10).pow(places);
    let q = Rational::new(numer * Integer::from(sign), scale);
    let half_unit = Mag::from_integer(&Integer::from(1))
        .div(&Mag::from_integer(&Integer::from(10).pow(places)));
    Ball::from_rational(&q, 400).add_error(&half_unit)
}

/// The ball `NIntegrate` gives for `input`, without a warning.
fn integrated(ctx: &Context, input: &str) -> ComplexBall {
    let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
    let found = found.unwrap();
    assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    match found.view() {
        View::Ball(b) => ComplexBall::from_real(b.clone()),
        View::ComplexBall(z) => z.clone(),
        _ => panic!("{input}: {found:?}"),
    }
}

/// Integrals known to 44 digits: the ball holds them, to 48 bits at least
/// at the default precision.
#[test]
fn known_values() {
    let ctx = Context::new();
    for (input, value) in [
        (
            "NIntegrate[Sin[Sin[x]], x, 0, 1]",
            "0.43060610312069060491237735524846578643360804",
        ),
        (
            "NIntegrate[Times[Exp[Times[-1, Power[x, 2]]], Cos[Times[3, x]]], x, 0, 2]",
            "0.090009235156271950217905593517629542216415123068",
        ),
        (
            "NIntegrate[Exp[Times[-1, Power[x, 2]]], x, 0, 1]",
            "0.746824132812427025399467436131853005354499687",
        ),
        (
            "NIntegrate[Power[Plus[1, Power[x, 4]], -1], x, 0, 1]",
            "0.866972987339911037573995163882870713652175368",
        ),
        (
            "NIntegrate[Times[Log[Plus[1, x]], Power[Plus[1, Power[x, 2]], -1]], x, 0, 1]",
            "0.272198261287950266312586112279701743417322963",
        ),
    ] {
        let found = integrated(&ctx, input);
        assert!(found.is_real(), "{input}: {found:?}");
        assert!(found.re().overlaps(&reference(value)), "{input}: {found:?}");
        assert!(found.re().rel_accuracy_bits() >= 48, "{input}: {found:?}");
    }
}

/// Exact values: switches split, a square root at an end, a complex
/// integrand, the list form, reversed bounds, real numbers as bounds.
#[test]
fn exact_values() {
    let ctx = Context::new();
    let q = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
    for (input, re, im) in [
        ("NIntegrate[Abs[Sin[x]], x, 0, Pi]", q(2, 1), q(0, 1)),
        ("NIntegrate[Sign[Plus[x, -1/3]], x, 0, 1]", q(1, 3), q(0, 1)),
        (
            "NIntegrate[HeavisideTheta[Plus[x, -1/4]], x, 0, 1]",
            q(3, 4),
            q(0, 1),
        ),
        ("NIntegrate[Power[x, 1/2], x, 0, 1]", q(2, 3), q(0, 1)),
        ("NIntegrate[Exp[Times[I, x]], x, 0, Pi]", q(0, 1), q(2, 1)),
        ("NIntegrate[Power[x, 2], List[x, 0, 1]]", q(1, 3), q(0, 1)),
        ("NIntegrate[Power[x, 2], x, 1, 0]", q(-1, 3), q(0, 1)),
        ("NIntegrate[Power[x, -1], x, 1, E]", q(1, 1), q(0, 1)),
    ] {
        let found = integrated(&ctx, input);
        assert!(
            found.re().contains_rational(&re) && found.im().contains_rational(&im),
            "{input}: {found:?}"
        );
        assert!(found.rel_accuracy_bits() >= 40, "{input}: {found:?}");
    }
}

/// What stays: an unknown function, a parameter, an infinite bound,
/// without a warning; a pole inside, past the budget, with one.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "NIntegrate[f[x], x, 0, 1]",
        "NIntegrate[Times[a, x], x, 0, 1]",
        "NIntegrate[x, x, 0, Infinity]",
        "NIntegrate[x, Plus[1, x], 0, 1]",
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
    let e = read("NIntegrate[Power[x, -1], x, -1, 1]");
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head());
    assert!(!report.warnings.is_empty());
}
