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
