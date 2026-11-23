//! Tests of `HeavisideTheta` and `DiracDelta`, and of their calculus
//! (CAS-273).

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
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
fn values() {
    check(
        &Context::new(),
        &[
            ("HeavisideTheta[0]", "1/2"),
            ("HeavisideTheta[2]", "1"),
            ("HeavisideTheta[-1]", "0"),
            ("HeavisideTheta[Pi]", "1"),
            ("HeavisideTheta[Plus[1, Times[-1, Power[2, 1/2]]]]", "0"),
            ("HeavisideTheta[x]", "HeavisideTheta[x]"),
            (
                "HeavisideTheta[Complex[1, 1]]",
                "HeavisideTheta[Complex[1, 1]]",
            ),
            ("DiracDelta[1]", "0"),
            ("DiracDelta[-1/2]", "0"),
            ("DiracDelta[0]", "ComplexInfinity"),
            ("DiracDelta[x]", "DiracDelta[x]"),
        ],
    );
    let positive = Context::new().assuming(read("Greater[x, 0]"));
    check(
        &positive,
        &[
            ("HeavisideTheta[x]", "1"),
            ("HeavisideTheta[Times[-1, x]]", "0"),
            ("DiracDelta[x]", "0"),
        ],
    );
}

#[test]
fn calculus() {
    check(
        &Context::new(),
        &[
            ("D[HeavisideTheta[x], x]", "DiracDelta[x]"),
            (
                "D[HeavisideTheta[Plus[Power[x, 2], -1]], x]",
                "Times[2, x, DiracDelta[Plus[Power[x, 2], -1]]]",
            ),
            (
                "Integrate[HeavisideTheta[Plus[x, -1]], x]",
                "Times[Plus[x, -1], HeavisideTheta[Plus[x, -1]]]",
            ),
            (
                "Integrate[DiracDelta[Times[2, x]], x]",
                "Times[1/2, HeavisideTheta[Times[2, x]]]",
            ),
            // Sifting.
            (
                "Integrate[Times[DiracDelta[Plus[x, -1]], Sin[x]], x, 0, 3]",
                "Sin[1]",
            ),
            ("Integrate[DiracDelta[x], x, -1, 1]", "1"),
            ("Integrate[DiracDelta[Plus[x, -5]], x, 0, 3]", "0"),
            (
                "Integrate[DiracDelta[Plus[Times[2, x], -1]], x, 0, 3]",
                "1/2",
            ),
            (
                "Integrate[Times[x, DiracDelta[Plus[x, -2]]], x, 3, 0]",
                "-2",
            ),
            (
                "Integrate[Times[Power[E, x], DiracDelta[Plus[x, -1]]], x, 0, Infinity]",
                "E",
            ),
        ],
    );
}

/// `LinearStep` (CAS-570), giac's answers but the downward quadratic.
#[test]
fn linear_step() {
    let h = |u: &str| format!("HeavisideTheta[{u}]");
    check(
        &Context::new(),
        &[
            (
                "LinearStep[Times[HeavisideTheta[Plus[x, -2]], HeavisideTheta[Plus[3, Times[-1, x]]]], x]",
                "Plus[HeavisideTheta[Plus[x, -2]], Times[-1, HeavisideTheta[Plus[x, -3]]]]",
            ),
            (
                "LinearStep[Piecewise[Less[x, -1], 0, Less[x, 1], Plus[x, 1], Plus[2, Times[-1, Power[Plus[x, -1], 2]]]], x]",
                "Plus[Times[HeavisideTheta[Plus[x, 1]], Plus[x, 1]], Times[HeavisideTheta[Plus[x, -1]], Plus[x, Times[-1, Power[x, 2]]]]]",
            ),
            // The variable defaults to x.
            (
                "LinearStep[Abs[x]]",
                "Plus[Times[2, x, HeavisideTheta[x]], Times[-1, x]]",
            ),
            (
                "LinearStep[Abs[Plus[Abs[x], Times[-1, Abs[Plus[x, -1]]]]]]",
                &format!(
                    "Plus[1, Times[-2, x, {}], Times[{}, Plus[2, Times[-2, x]]], Times[{}, Plus[-2, Times[4, x]]]]",
                    h("x"),
                    h("Plus[x, -1]"),
                    h("Plus[x, -1/2]")
                ),
            ),
            (
                "LinearStep[Times[Sign[Plus[t, -1]], t], t]",
                "Plus[Times[2, t, HeavisideTheta[Plus[t, -1]]], Times[-1, t]]",
            ),
            (
                "LinearStep[Times[HeavisideTheta[Plus[Times[2, x], -1]], HeavisideTheta[x]], x]",
                "HeavisideTheta[Plus[x, -1/2]]",
            ),
            (
                "LinearStep[Power[HeavisideTheta[x], 2], x]",
                "HeavisideTheta[x]",
            ),
            (
                "LinearStep[HeavisideTheta[Plus[Power[x, 2], -1]], x]",
                "Plus[1, HeavisideTheta[Plus[x, -1]], Times[-1, HeavisideTheta[Plus[x, 1]]]]",
            ),
            // giac's box has the opposite sign here.
            (
                "LinearStep[HeavisideTheta[Plus[1, Times[-1, Power[x, 2]]]], x]",
                "Plus[HeavisideTheta[Plus[x, 1]], Times[-1, HeavisideTheta[Plus[x, -1]]]]",
            ),
            ("LinearStep[f, 2]", "LinearStep[f, 2]"),
        ],
    );
}

/// `SimplifyDirac` (CAS-572), giac's answers.
#[test]
fn simplify_dirac() {
    check(
        &Context::new(),
        &[
            (
                "SimplifyDirac[DiracDelta[Times[a, x]], x]",
                "Times[Power[Abs[a], -1], DiracDelta[x]]",
            ),
            (
                "SimplifyDirac[Times[DiracDelta[x, 1], x]]",
                "Times[-1, DiracDelta[x]]",
            ),
            (
                "SimplifyDirac[Times[Derivative[1][DiracDelta][x], x]]",
                "Times[-1, DiracDelta[x]]",
            ),
            (
                "SimplifyDirac[DiracDelta[Plus[Power[x, 2], Times[-1, x], -6]]]",
                "Plus[Times[1/5, DiracDelta[Plus[x, 2]]], Times[1/5, DiracDelta[Plus[x, -3]]]]",
            ),
            ("SimplifyDirac[DiracDelta[Times[-1, x]]]", "DiracDelta[x]"),
            ("SimplifyDirac[Times[DiracDelta[x], Sin[x]]]", "0"),
            ("SimplifyDirac[DiracDelta[Plus[Power[x, 2], 1]]]", "0"),
            (
                "SimplifyDirac[Times[DiracDelta[Plus[x, -1], 2], Power[x, 3]]]",
                "Plus[Times[6, DiracDelta[Plus[x, -1]]], Times[-6, DiracDelta[Plus[x, -1], 1]], DiracDelta[Plus[x, -1], 2]]",
            ),
            (
                "SimplifyDirac[Times[DiracDelta[Plus[Times[2, t], -2], 1], Power[t, 2]], t]",
                "Plus[Times[-1/2, DiracDelta[Plus[t, -1]]], Times[1/4, DiracDelta[Plus[t, -1], 1]]]",
            ),
            (
                "SimplifyDirac[Times[DiracDelta[x], Power[x, -1]]]",
                "Times[-1, DiracDelta[x, 1]]",
            ),
            (
                "SimplifyDirac[Plus[Times[DiracDelta[x], Cos[x]], Power[Sin[x], 2], Power[Cos[x], 2]]]",
                "Plus[1, DiracDelta[x]]",
            ),
            ("SimplifyDirac[Power[x, 2]]", "Power[x, 2]"),
        ],
    );
}
