//! Tests of `Simplify` (CAS-229): the complexity, the search through the
//! rational candidates, and what is simplified piece by piece.

use super::complexity;
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

/// §2's measure, on expressions as read.
#[test]
fn complexities() {
    for (text, expected) in [
        ("x", 1),
        ("0", 1),
        ("7", 1),
        ("-7", 2),
        ("1234", 4),
        ("-1234", 5),
        ("1000000000000000000000000000000", 31),
        ("1/2", 3),
        ("-12/7", 5),
        ("Plus[1, x]", 3),
        ("Exp[Plus[x, y]]", 4),
        ("Times[Exp[x], Exp[y]]", 5),
        ("Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]", 14),
    ] {
        assert_eq!(complexity(&read(text)), expected, "{text}");
    }
}

/// Rational normalization, inside functions, and piece by piece.
#[test]
fn rational_simplification() {
    check(&[
        (
            "Simplify[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]",
            "Plus[1, x]",
        ),
        (
            "Simplify[Times[x, Power[Plus[x, Power[x, 2]], -1]]]",
            "Power[Plus[1, x], -1]",
        ),
        (
            "Simplify[Times[Plus[-1, Power[x, 3]], Power[Plus[1, x, Power[x, 2]], -1]]]",
            "Plus[-1, x]",
        ),
        (
            "Simplify[Plus[Power[Plus[x, 1], 2], Times[-1, Power[x, 2]], Times[-2, x]]]",
            "1",
        ),
        (
            "Simplify[Power[Plus[Power[x, -1], Power[y, -1]], -1]]",
            "Times[x, y, Power[Plus[x, y], -1]]",
        ),
        (
            "Simplify[Times[Plus[Power[a, 2], Times[-1, Power[b, 2]]], Power[Plus[a, b], -1]]]",
            "Plus[a, Times[-1, b]]",
        ),
        // The arguments of functions first.
        (
            "Simplify[Sin[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]]",
            "Sin[Plus[1, x]]",
        ),
        // Lists and equations piece by piece.
        (
            "Simplify[List[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]], y]]",
            "List[Plus[1, x], y]",
        ),
        (
            "Simplify[Equation[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]], y]]",
            "Equation[Plus[1, x], y]",
        ),
        // Already the simplest: itself.
        ("Simplify[x]", "x"),
        ("Simplify[Plus[1, x]]", "Plus[1, x]"),
        ("Simplify[Times[2, x]]", "Times[2, x]"),
        ("Simplify[Power[Plus[1, x], 10]]", "Power[Plus[1, x], 10]"),
        // A held argument stays as it is.
        (
            "Simplify[Hold[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]]",
            "Hold[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]",
        ),
    ]);
}

/// A simplification is never more complex than what it simplifies.
#[test]
fn never_more_complex() {
    let ctx = Context::new();
    for text in [
        "Power[Plus[1, x], 10]",
        "Times[Plus[1, x], Plus[2, x], Plus[3, x]]",
        "Plus[Power[x, 5], Times[-1, Power[y, 5]]]",
        "Times[Power[Plus[1, x], -1], Plus[1, Power[x, 2]]]",
        "Plus[Sin[x], Times[2, Cos[y]], Power[x, 1/2]]",
        "Times[Plus[Power[x, 2], 1], Power[Plus[Power[x, 4], -1], -1]]",
    ] {
        let e = ctx.eval(&read(text)).unwrap();
        let simplified = ctx
            .eval(&Expr::call(crate::Symbol::SIMPLIFY, [e.clone()]))
            .unwrap();
        assert!(
            complexity(&simplified) <= complexity(&e),
            "{text}: {simplified:?}"
        );
    }
}
