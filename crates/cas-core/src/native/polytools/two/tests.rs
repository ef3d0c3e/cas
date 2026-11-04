//! Tests of the polynomial tools II (CAS-482): giac's answers.

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

#[test]
fn roots_and_coefficients() {
    check(&[
        // giac's [-i,i], and exact rational roots with multiplicities.
        (
            "PolynomialRoots[List[1, 0, 1]]",
            "List[Complex[0, -1], Complex[0, 1]]",
        ),
        (
            "PolynomialRoots[Times[Power[Plus[x, -1], 2], Plus[Times[2, x], 1]]]",
            "List[-1/2, 1, 1]",
        ),
        // giac's poly1[1,-6,11,-6].
        ("FromRoots[List[1, 2, 3]]", "List[1, -6, 11, -6]"),
        (
            "FromRoots[List[Complex[0, -1], Complex[0, 1]]]",
            "List[1, 0, 1]",
        ),
        // giac's 11 twice and 13.
        ("HornerEval[List[1, 2, 3], 2]", "11"),
        ("HornerEval[Plus[Power[x, 3], Times[2, x], 1], 2]", "13"),
        (
            "HornerEval[List[1, 2, 3], a]",
            "Plus[3, Times[2, a], Power[a, 2]]",
        ),
        // giac's [1,6,14,13] and x^3+6*x^2+14*x+13.
        ("TaylorShift[List[1, 0, 2, 1], 2]", "List[1, 6, 14, 13]"),
        (
            "TaylorShift[Plus[Power[x, 3], Times[2, x], 1], 2]",
            "Plus[13, Times[14, x], Times[6, Power[x, 2]], Power[x, 3]]",
        ),
        // giac's [0,2,3,3], and a pole.
        (
            "RootMultiplicities[Times[Power[x, 2], Power[Plus[x, -3], 3]]]",
            "List[0, 2, 3, 3]",
        ),
        (
            "RootMultiplicities[Times[Power[Plus[x, -1], 2], Power[Plus[x, 2], -1]]]",
            "List[-2, -1, 1, 2]",
        ),
        // The names e2r, r2e and ldegree are CoefficientList's,
        // FromCoefficients's and Valuation's.
    ]);
    // The other roots are balls, real ones first: x³ − 1.
    let ctx = Context::new();
    let roots = ctx
        .eval(&read("PolynomialRoots[Plus[Power[x, 3], -1]]"))
        .unwrap();
    assert_eq!(roots.args().len(), 3, "{roots:?}");
    let last = &roots.args()[2];
    let near_one = ctx
        .eval(&Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [last.clone() - Expr::from(1)]),
                read("Power[10, -10]"),
            ],
        ))
        .unwrap();
    assert_eq!(near_one.as_symbol(), Some(Symbol::TRUE), "{roots:?}");
}
