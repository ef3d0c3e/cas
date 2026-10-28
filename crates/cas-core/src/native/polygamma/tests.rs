//! Tests of `PolyGamma` and `EulerGamma` (CAS-355): giac's exact values at
//! integers, half-integers and thirds, the poles, the orders, balls against
//! giac's 50 digits, the derivatives of `PolyGamma`, `Gamma` and `Beta`,
//! and what stays.

use cas_num::ball::Ball;

use crate::expr::View;
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same once put together.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        let difference = ctx
            .eval(&Expr::call(
                Symbol::TOGETHER,
                [Expr::call(
                    Symbol::EXPAND,
                    [found.clone() - expected.clone()],
                )],
            ))
            .unwrap();
        if found != expected && difference != Expr::from(0) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `Psi(1)`, `Psi(3)`, `Psi(1/2)`, `Psi(5/2)`, `Psi(-1/2)`,
/// `Psi(2,1)`, `Psi(1,3)` and `Psi(1/2,1)`.
#[test]
fn exact_values() {
    check(&[
        ("PolyGamma[1]", "Times[-1, EulerGamma]"),
        ("PolyGamma[3]", "Plus[3/2, Times[-1, EulerGamma]]"),
        (
            "PolyGamma[1/2]",
            "Plus[Times[-1, EulerGamma], Times[-2, Log[2]]]",
        ),
        (
            "PolyGamma[5/2]",
            "Plus[8/3, Times[-1, EulerGamma], Times[-2, Log[2]]]",
        ),
        (
            "PolyGamma[-1/2]",
            "Plus[2, Times[-1, EulerGamma], Times[-2, Log[2]]]",
        ),
        ("PolyGamma[2, 1]", "Plus[-1, Times[1/6, Power[Pi, 2]]]"),
        ("PolyGamma[1, 3]", "Times[1/15, Power[Pi, 4]]"),
        ("PolyGamma[1/2, 1]", "Times[1/2, Power[Pi, 2]]"),
        ("PolyGamma[1, 2]", "Times[-2, Zeta[3]]"),
        ("PolyGamma[x, 0]", "PolyGamma[x]"),
        ("PolyGamma[Infinity]", "Infinity"),
        ("PolyGamma[Infinity, 2]", "0"),
        // The poles.
        ("PolyGamma[0]", "ComplexInfinity"),
        ("PolyGamma[-3, 2]", "ComplexInfinity"),
        // What stays.
        ("PolyGamma[n]", "PolyGamma[n]"),
        ("PolyGamma[Plus[1, x]]", "PolyGamma[Plus[1, x]]"),
        ("PolyGamma[1/3, 1]", "PolyGamma[1/3, 1]"),
        ("PolyGamma[1/13]", "PolyGamma[1/13]"),
        ("PolyGamma[2, 1/2]", "PolyGamma[2, 1/2]"),
    ]);
}

/// `e` as a real ball.
fn ball_of(e: &Expr) -> Ball {
    match e.view() {
        View::Ball(b) => b.clone(),
        _ => panic!("{e:?} is not a real ball"),
    }
}

/// Gauss's theorem at thirds, quarters and fifths, and balls, against
/// giac's 50 digits; `EulerGamma`'s value and sign.
#[test]
fn numbers() {
    let mut ctx = Context::new();
    ctx.set_working_precision(150);
    for (input, digits) in [
        (
            "PolyGamma[1/3]",
            "-3.1320337800208063229964190742872688541554282967204",
        ),
        (
            "PolyGamma[2/3]",
            "-1.3182344157865884724023408166451131218713620486277",
        ),
        (
            "PolyGamma[1/4]",
            "-4.2274535333762654080895301460966835773672444387082",
        ),
        (
            "PolyGamma[3/5]",
            "-1.5406192138931904147606639488062319415105342546896",
        ),
        (
            "PolyGamma[7/3]",
            "0.61796621997919367700358092571273114584457170327958",
        ),
        (
            "PolyGamma[5/2, 2]",
            "-0.23620405164172740300374166856770727811721550017431",
        ),
        (
            "PolyGamma[N[-5/4]]",
            "3.7141391202135278303731132371828193068299249606669",
        ),
        (
            "EulerGamma",
            "0.57721566490153286060651209008240243104215933593992",
        ),
    ] {
        let found = ctx.eval(&read(&format!("N[{input}]"))).unwrap();
        let b = ball_of(&found);
        let want = Ball::from_decimal(digits, 200).unwrap();
        assert!(b.overlaps(&want), "{input}: {b:?} against {digits}");
        assert!(b.rel_accuracy_bits() >= 140, "{input}: {b:?}");
    }
    check(&[
        ("Greater[EulerGamma, 1/2]", "True"),
        ("Less[EulerGamma, Times[1/5, Pi]]", "True"),
    ]);
    // A ball meeting a pole.
    let found = ctx.eval(&read("PolyGamma[Ball[-2, 1/1024]]")).unwrap();
    assert_eq!(found, read("Indeterminate"));
}

/// giac's `diff(Psi(x),x)`, `diff(Psi(x,2),x)`, `diff(Gamma(x),x)` and
/// `diff(Beta(x,y),x)`.
#[test]
fn derivatives() {
    check(&[
        ("D[PolyGamma[x], x]", "PolyGamma[x, 1]"),
        ("D[PolyGamma[x, 2], x]", "PolyGamma[x, 3]"),
        (
            "D[PolyGamma[Power[x, 2]], x]",
            "Times[2, x, PolyGamma[Power[x, 2], 1]]",
        ),
        ("D[Gamma[x], x]", "Times[Gamma[x], PolyGamma[x]]"),
        (
            "D[Beta[x, y], x]",
            "Times[Beta[x, y], Plus[PolyGamma[x], Times[-1, PolyGamma[Plus[x, y]]]]]",
        ),
    ]);
}
