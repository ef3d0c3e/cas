//! Tests of implicit differentiation (CAS-412): giac's examples, higher
//! orders, and what stays.

use crate::{Context, Expr};

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
fn derivatives() {
    let circle = "Equation[Plus[Power[x, 2], Power[y, 2]], 1]";
    let second = format!("ImplicitD[{circle}, y, x, x]");
    let third = format!("ImplicitD[{circle}, y, List[x, 3]]");
    check(&[
        (
            &format!("ImplicitD[{circle}, y, x]"),
            "Times[-1, x, Power[y, -1]]",
        ),
        // giac's other order: dx/dy.
        (
            &format!("ImplicitD[{circle}, x, y]"),
            "Times[-1, y, Power[x, -1]]",
        ),
        (
            &second,
            "Together[Times[Plus[Times[-1, Power[x, 2]], Times[-1, Power[y, 2]]], Power[y, -3]]]",
        ),
        (
            &third,
            "Together[Times[Plus[Times[-3, Power[x, 3]], Times[-3, x, Power[y, 2]]], Power[y, -5]]]",
        ),
        (
            "ImplicitD[Equation[Plus[Times[Power[x, 2], y], Power[y, 2]], 1], y, x]",
            "Together[Times[-2, x, y, Power[Plus[Power[x, 2], Times[2, y]], -1]]]",
        ),
        // The other symbols are constants.
        (
            "ImplicitD[Equation[R, Times[P, V, Power[T, -1]]], P, T]",
            "Times[P, Power[T, -1]]",
        ),
        // An expression is read as itself = 0, and `==` as `=`.
        (
            "ImplicitD[Plus[Power[x, 2], Power[y, 2], -1], y, x]",
            "Times[-1, x, Power[y, -1]]",
        ),
        (
            "ImplicitD[Equal[Plus[Power[x, 2], Power[y, 2]], 1], y, x]",
            "Times[-1, x, Power[y, -1]]",
        ),
        // What stays: one variable twice, a number for a variable, two
        // variables of differentiation.
        ("ImplicitD[Plus[x, y], y, y]", "ImplicitD[Plus[x, y], y, y]"),
        ("ImplicitD[Plus[x, y], 1, x]", "ImplicitD[Plus[x, y], 1, x]"),
        (
            "ImplicitD[Plus[x, y, z], y, x, z]",
            "ImplicitD[Plus[x, y, z], y, x, z]",
        ),
    ]);
}
