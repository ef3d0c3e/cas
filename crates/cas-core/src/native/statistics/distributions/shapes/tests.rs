//! Tests of the gamma, beta, Weibull and Cauchy distributions (CAS-470):
//! giac's densities, cumulative functions and quantiles, the supports, and
//! what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
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

/// The value of `input` is within `10⁻¹¹` of `expected`, giac's decimal.
fn check_near(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let gap = Expr::call(
            Symbol::N,
            [Expr::call(
                Symbol::ABS,
                [Expr::call(
                    Symbol::PLUS,
                    [found.clone(), Expr::from(-1) * read(expected)],
                )],
            )],
        );
        let small = Expr::call(Symbol::LESS, [gap, read("Power[10, -11]")]);
        if ctx.eval(&small).unwrap().as_symbol() != Some(Symbol::TRUE) {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn densities() {
    check(&[
        // giac's x*exp(-3*x)*9 and 12*x*(1-x)^2.
        ("GammaPDF[2, 3, x]", "Times[9, x, Exp[Times[-3, x]]]"),
        (
            "BetaPDF[2, 3, x]",
            "Times[12, x, Power[Plus[1, Times[-1, x]], 2]]",
        ),
        // giac's 2/3*x/3*exp(-(x/3)^2) and 1/pi/(x^2+1).
        (
            "WeibullPDF[2, 3, x]",
            "Times[2/9, x, Exp[Times[-1, Power[Times[1/3, x], 2]]]]",
        ),
        ("CauchyPDF[x]", "Power[Times[Pi, Plus[1, Power[x, 2]]], -1]"),
        (
            "CauchyPDF[1, 2, x]",
            "Times[2, Power[Times[Pi, Plus[Power[Plus[x, -1], 2], 4]], -1]]",
        ),
        // Outside the supports, 0, where giac's formulas go on.
        ("GammaPDF[2, 3, -1]", "0"),
        ("BetaPDF[2, 3, 2]", "0"),
        ("WeibullPDF[2, 3, 1, 1/2]", "0"),
        // Invalid parameters stay.
        ("GammaPDF[-1, 3, x]", "GammaPDF[-1, 3, x]"),
        ("CauchyPDF[0, 0, x]", "CauchyPDF[0, 0, x]"),
    ]);
}

#[test]
fn cumulative_functions() {
    check(&[
        // giac's igamma(2,3*x,1), here in closed form: 1 − e^(−3x)(1 + 3x).
        (
            "GammaCDF[2, 3, x]",
            "Plus[1, Times[-1, Exp[Times[-3, x]], Plus[1, Times[3, x]]]]",
        ),
        ("GammaCDF[2, 3, -1]", "0"),
        // giac's 0.5248 = 328/625, and 0.6352 = 397/625 on [1/10, 1/2].
        ("BetaCDF[2, 3, 2/5]", "328/625"),
        ("BetaCDF[2, 3, 1/10, 1/2]", "397/625"),
        ("BetaCDF[2, 3, 2]", "1"),
        // giac's 1-exp(-(x/3)^2), 1-exp(-((x-1)/3)^2), 1/2+atan(x)/pi.
        (
            "WeibullCDF[2, 3, x]",
            "Plus[1, Times[-1, Exp[Times[-1, Power[Times[1/3, x], 2]]]]]",
        ),
        (
            "WeibullCDF[2, 3, 1, x]",
            "Plus[1, Times[-1, Exp[Times[-1, Power[Times[1/3, Plus[x, -1]], 2]]]]]",
        ),
        ("CauchyCDF[x]", "Plus[1/2, Times[Power[Pi, -1], ArcTan[x]]]"),
        ("CauchyCDF[-1, 1]", "1/2"),
    ]);
}

#[test]
fn quantiles() {
    check(&[
        // giac's 3.0, and 3 √(log 2) for giac's 2.49766383347.
        ("CauchyICDF[1, 2, 3/4]", "3"),
        ("CauchyICDF[1/2]", "0"),
        ("WeibullICDF[2, 3, 1/2]", "Times[3, Power[Log[2], 1/2]]"),
        (
            "WeibullICDF[2, 3, 1, 1/2]",
            "Plus[1, Times[3, Power[Log[2], 1/2]]]",
        ),
        ("GammaICDF[2, 3, 0]", "0"),
        ("BetaICDF[2, 3, 1]", "1"),
        ("GammaICDF[a, 3, 1/2]", "GammaICDF[a, 3, 1/2]"),
    ]);
    // giac's 0.559448996672 and 0.385727568132.
    check_near(&[
        ("GammaICDF[2, 3, 1/2]", "0.559448996672"),
        ("BetaICDF[2, 3, 1/2]", "0.385727568132"),
    ]);
}

#[test]
fn variates() {
    let ctx = Context::new();
    for input in [
        "RandomGamma[2, 3]",
        "RandomBeta[2, 3]",
        "RandomWeibull[2, 3]",
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert!(found.is_number(), "{input}: {found:?}");
        let positive = ctx
            .eval(&Expr::call(Symbol::GREATER, [found.clone(), Expr::from(0)]))
            .unwrap();
        assert_eq!(
            positive.as_symbol(),
            Some(Symbol::TRUE),
            "{input}: {found:?}"
        );
    }
}
