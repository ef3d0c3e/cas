//! Tests of the Bessel and Airy builtins (CAS-382): exact values, the
//! symmetries of orders and arguments, balls against cas-num's, what stays,
//! and the derivatives.

use cas_num::ball::{Airy, Ball, Float};

use crate::{Context, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The values at 0, and the signs of negative orders and arguments.
#[test]
fn exact_values() {
    check(&[
        ("BesselJ[0, 0]", "1"),
        ("BesselJ[3, 0]", "0"),
        ("BesselI[0, 0]", "1"),
        ("BesselI[2, 0]", "0"),
        ("BesselJ[-3, x]", "Times[-1, BesselJ[3, x]]"),
        ("BesselJ[-2, x]", "BesselJ[2, x]"),
        ("BesselY[-1, x]", "Times[-1, BesselY[1, x]]"),
        ("BesselI[-3, x]", "BesselI[3, x]"),
        ("BesselK[-2, x]", "BesselK[2, x]"),
        ("BesselJ[3, -2]", "Times[-1, BesselJ[3, 2]]"),
        ("BesselJ[2, -2]", "BesselJ[2, 2]"),
        ("BesselI[1, -1/2]", "Times[-1, BesselI[1, 1/2]]"),
        ("AiryAi[0]", "Times[Power[3, -2/3], Power[Gamma[2/3], -1]]"),
        ("AiryBi[0]", "Times[Power[3, -1/6], Power[Gamma[2/3], -1]]"),
        (
            "AiryAiPrime[0]",
            "Times[-1, Power[3, -1/3], Power[Gamma[1/3], -1]]",
        ),
        (
            "AiryBiPrime[0]",
            "Times[Power[3, 1/6], Power[Gamma[1/3], -1]]",
        ),
    ]);
}

/// What stays: exact arguments, symbols, other orders, `Y` and `K` at 0
/// and below.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "BesselJ[0, 1]",
        "BesselJ[n, x]",
        "BesselJ[1/2, x]",
        "BesselY[0, 0]",
        "BesselK[1, 0]",
        "BesselY[0, -2]",
        "BesselK[0, -1/2]",
        "AiryAi[x]",
        "AiryBi[2]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
    // A negative ball for Y stays a call.
    let found = ctx.eval(&read("N[BesselY[0, -5/2]]")).unwrap();
    assert_eq!(
        found.head(),
        Some(&Expr::from(Symbol::BESSEL_Y)),
        "{found:?}"
    );
}

/// `N` of each function holds cas-num's value of it.
#[test]
fn balls() {
    let ctx = Context::new();
    let x = Ball::exact(Float::from(3i64).mul_2exp(&-cas_num::Integer::ONE));
    let cases: [(&str, Ball); 8] = [
        ("BesselJ[2, 3/2]", x.bessel_j(2, 53)),
        ("BesselY[1, 3/2]", x.bessel_y(1, 53)),
        ("BesselI[0, 3/2]", x.bessel_i(0, 53)),
        ("BesselK[3, 3/2]", x.bessel_k(3, 53)),
        ("AiryAi[3/2]", x.airy(Airy::Ai, 53)),
        ("AiryBi[3/2]", x.airy(Airy::Bi, 53)),
        ("AiryAiPrime[3/2]", x.airy(Airy::AiPrime, 53)),
        ("AiryBiPrime[3/2]", x.airy(Airy::BiPrime, 53)),
    ];
    for (input, expected) in cases {
        let found = ctx.eval(&Expr::call(Symbol::N, [read(input)])).unwrap();
        let View::Ball(b) = found.view() else {
            panic!("N[{input}]: {found:?}");
        };
        assert!(
            b.overlaps(&expected),
            "N[{input}]: {b:?} against {expected:?}"
        );
        assert!(b.rel_accuracy_bits() >= 50, "N[{input}]: {b:?}");
    }
    // N[AiryAi[0]] through its closed form.
    let found = ctx.eval(&read("N[AiryAi[0]]")).unwrap();
    let View::Ball(b) = found.view() else {
        panic!("N[AiryAi[0]]: {found:?}");
    };
    assert!(b.overlaps(&Ball::ZERO.airy(Airy::Ai, 53)), "{b:?}");
}

/// The derivatives: `(C_(n−1) ∓ C_(n+1))/2` and Airy's equation.
#[test]
fn derivatives() {
    check(&[
        ("D[BesselJ[0, x], x]", "Times[-1, BesselJ[1, x]]"),
        (
            "D[BesselJ[n, x], x]",
            "Times[1/2, Plus[BesselJ[Plus[-1, n], x], Times[-1, BesselJ[Plus[1, n], x]]]]",
        ),
        ("D[BesselI[0, x], x]", "BesselI[1, x]"),
        ("D[BesselK[0, x], x]", "Times[-1, BesselK[1, x]]"),
        (
            "D[BesselY[1, Times[2, x]], x]",
            "Plus[BesselY[0, Times[2, x]], Times[-1, BesselY[2, Times[2, x]]]]",
        ),
        ("D[AiryAi[x], x]", "AiryAiPrime[x]"),
        ("D[AiryAiPrime[x], x]", "Times[x, AiryAi[x]]"),
        (
            "D[AiryBi[Times[3, x]], x]",
            "Times[3, AiryBiPrime[Times[3, x]]]",
        ),
        ("D[AiryBiPrime[x], x]", "Times[x, AiryBi[x]]"),
    ]);
}
