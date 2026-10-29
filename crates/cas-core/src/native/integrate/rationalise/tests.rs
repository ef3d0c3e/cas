//! Tests of the substitutions that rationalise integrands (CAS-DOC-40): on
//! giac's examples, the antiderivative whole, its derivative the integrand
//! at three points; giac's forms where they are this CAS's; what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a` and `b`, numbers, are equal balls once `N` gives them.
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let n = |e: &Expr| ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]));
    let (Ok(a), Ok(b)) = (n(a), n(b)) else {
        return false;
    };
    let gap = Expr::call(Symbol::ABS, [a - b]);
    let tolerance = read("1/1000000000000000000000");
    ctx.eval(&Expr::call(Symbol::LESS, [gap, tolerance]))
        .is_ok_and(|t| t.as_symbol() == Some(Symbol::TRUE))
}

/// Each integrand's antiderivative: whole, and its derivative the integrand
/// at `x` = 3/7, 5/2 and −4/3.
fn check_derivatives(integrands: &[&str]) {
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let mut wrong = Vec::new();
    for f in integrands {
        let f = read(f);
        let found = ctx
            .eval(&Expr::call(Symbol::INTEGRATE, [f.clone(), x.clone()]))
            .unwrap();
        if found.contains(&Expr::from(Symbol::INTEGRATE)) {
            wrong.push(format!("{f:?}: {found:?}"));
            continue;
        }
        let derivative = ctx
            .eval(&Expr::call(Symbol::D, [found.clone(), x.clone()]))
            .unwrap();
        for (p, q) in [(3, 7), (5, 2), (-4, 3)] {
            let point = [(x.clone(), Expr::from(p) / Expr::from(q))];
            let (d, g) = (derivative.substitute(&point), f.substitute(&point));
            if !same_value(&ctx, &d, &g) {
                wrong.push(format!("{f:?}: {found:?} at {p}/{q}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Rational functions of exponentials and hyperbolic functions (CAS-433).
#[test]
fn exponentials() {
    check_derivatives(&[
        // giac's x-ln(exp(x)+1), atan(exp(x)), atan(exp(x)), exp(x)-ln(exp(x)+1).
        "Power[Plus[1, Exp[x]], -1]",
        "Times[Exp[x], Power[Plus[1, Exp[Times[2, x]]], -1]]",
        "Power[Plus[Exp[x], Exp[Times[-1, x]]], -1]",
        "Times[Exp[Times[2, x]], Power[Plus[1, Exp[x]], -1]]",
        // giac's 2*atan(exp(x)) and 2/(exp(2*x)+1)+x.
        "Power[Cosh[x], -1]",
        "Power[Tanh[x], 2]",
        "Power[Sinh[x], 2]",
        // Slopes 1/2 and 1: u = e^(x/2).
        "Times[Power[E, Times[1/2, x]], Power[Plus[1, Power[E, x]], -1]]",
        // A constant e^b, a factor of the rational function.
        "Times[Power[E, Plus[1, x]], Power[Plus[1, Power[E, Times[2, x]]], -1]]",
    ]);
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "Integrate[Power[Plus[Exp[x], Exp[Times[-1, x]]], -1], x]",
        ))
        .unwrap();
    assert_eq!(found, ctx.eval(&read("ArcTan[Exp[x]]")).unwrap());
}

/// Rational functions of sines and cosines (CAS-434): Bioche's three rules
/// and Weierstrass's substitution.
#[test]
fn sines_and_cosines() {
    check_derivatives(&[
        // u = cos x: 1/sin, giac's ln|tan(x/2)|; sin/(1 + cos), giac's
        // -ln(cos(x)+1); sin³/cos², giac's cos(x)+1/cos(x); cos²/sin; tan³.
        "Power[Sin[x], -1]",
        "Times[Sin[x], Power[Plus[1, Cos[x]], -1]]",
        "Times[Power[Sin[x], 3], Power[Cos[x], -2]]",
        "Times[Power[Cos[x], 2], Power[Sin[x], -1]]",
        "Power[Tan[x], 3]",
        // u = sin x: 1/cos, sec.
        "Power[Cos[x], -1]",
        "Sec[x]",
        // u = tan x: tan², 1/(1 + sin²), 1/cos(2x)².
        "Power[Tan[x], 2]",
        "Power[Plus[1, Power[Sin[x], 2]], -1]",
        "Power[Cos[Times[2, x]], -2]",
        // Weierstrass: 1/(2 + cos), 1/(1 + sin), 1/(3 + 5 sin), at 3x + 1.
        "Power[Plus[2, Cos[x]], -1]",
        "Power[Plus[1, Sin[x]], -1]",
        "Power[Plus[3, Times[5, Sin[x]]], -1]",
        "Power[Plus[2, Cos[Plus[1, Times[3, x]]]], -1]",
    ]);
    let ctx = Context::new();
    for (input, expected) in [
        // giac's tan(x)-x, ArcTan[Tan[x]] back as x.
        (
            "Integrate[Power[Tan[x], 2], x]",
            "Plus[Tan[x], Times[-1, x]]",
        ),
        // giac's cos(x)+1/cos(x).
        (
            "Integrate[Times[Power[Sin[x], 3], Power[Cos[x], -2]], x]",
            "Plus[Cos[x], Power[Cos[x], -1]]",
        ),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
}

/// What the substitutions leave: `x` outside the exponentials, a slope
/// that is not rational, an exponent that is not linear; two arguments of
/// sines, `x` outside them.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Integrate[Times[x, Power[Plus[1, Exp[x]], -1]], x]",
        "Integrate[Power[Plus[1, Exp[Times[Pi, x]], Exp[x]], -1], x]",
        "Integrate[Power[Plus[1, Exp[Power[x, 2]]], -1], x]",
        "Integrate[Power[Plus[2, Sin[x], Cos[Times[2, x]]], -1], x]",
        "Integrate[Times[x, Power[Plus[2, Cos[x]], -1]], x]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert!(
            found.contains(&Expr::from(Symbol::INTEGRATE)),
            "{input}: {found:?}"
        );
    }
}
