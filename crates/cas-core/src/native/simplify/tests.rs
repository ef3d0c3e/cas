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

/// Exponentials and logarithms (CAS-230).
#[test]
fn exponentials_and_logarithms() {
    check(&[
        ("Simplify[Exp[Times[2, Log[x]]]]", "Power[x, 2]"),
        ("Simplify[Exp[Plus[x, Log[2]]]]", "Times[2, Exp[x]]"),
        (
            "Simplify[Times[Exp[x], Power[Exp[Plus[-1, x]], -1]]]",
            "Exp[1]",
        ),
        ("Simplify[Times[Exp[x], Exp[y]]]", "Exp[Plus[x, y]]"),
        ("Simplify[Plus[Log[2], Log[3]]]", "Log[6]"),
        ("Simplify[Plus[Log[4], Times[-2, Log[2]]]]", "0"),
        ("Simplify[Times[Log[8], Power[Log[2], -1]]]", "3"),
        // A fractional power of x out of an exponent only for a positive
        // x: in the Real domain an odd root of a negative is real.
        (
            "Simplify[Exp[Times[1/3, Log[x]]]]",
            "Exp[Times[1/3, Log[x]]]",
        ),
        // Log[x^2] is 2 Log[Abs[x]] for a real x of unknown sign.
        (
            "Simplify[Plus[Log[Power[x, 2]], Times[-2, Log[x]]]]",
            "Plus[Log[Power[x, 2]], Times[-2, Log[x]]]",
        ),
    ]);
    // With x positive, Log[x^2] is 2 Log[x], and x^(1/3) comes out.
    let ctx = Context::new().assuming(read("Greater[x, 0]"));
    for (input, expected) in [
        ("Simplify[Plus[Log[Power[x, 2]], Times[-2, Log[x]]]]", "0"),
        ("Simplify[Exp[Times[1/3, Log[x]]]]", "Power[x, 1/3]"),
        ("Simplify[Plus[Log[x], Log[y]]]", "Plus[Log[x], Log[y]]"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
}

/// Radicals (CAS-232): denominators rationalized, square roots denested.
#[test]
fn radicals() {
    check(&[
        (
            "Simplify[Power[Plus[1, Power[2, 1/2]], -1]]",
            "Plus[-1, Power[2, 1/2]]",
        ),
        (
            "Simplify[Power[Plus[Power[3, 1/2], Times[-1, Power[2, 1/2]]], -1]]",
            "Plus[Power[2, 1/2], Power[3, 1/2]]",
        ),
        (
            "Simplify[Power[Plus[3, Times[2, Power[2, 1/2]]], 1/2]]",
            "Plus[1, Power[2, 1/2]]",
        ),
        (
            "Simplify[Power[Plus[5, Times[2, Power[6, 1/2]]], 1/2]]",
            "Plus[Power[2, 1/2], Power[3, 1/2]]",
        ),
        (
            "Simplify[Power[Plus[3, Times[-2, Power[2, 1/2]]], 1/2]]",
            "Plus[-1, Power[2, 1/2]]",
        ),
        (
            "Simplify[Times[Plus[1, Power[x, 1/2]], Plus[-1, Power[x, 1/2]]]]",
            "Plus[-1, x]",
        ),
        // A root that does not denest stays.
        (
            "Simplify[Power[Plus[1, Power[2, 1/2]], 1/2]]",
            "Power[Plus[1, Power[2, 1/2]], 1/2]",
        ),
    ]);
}

/// Trigonometric simplification (CAS-234).
#[test]
fn trigonometric() {
    check(&[
        ("Simplify[Plus[Power[Sin[x], 2], Power[Cos[x], 2]]]", "1"),
        (
            "Simplify[Times[Sin[Times[2, x]], Power[Sin[x], -1]]]",
            "Times[2, Cos[x]]",
        ),
        ("Simplify[Times[Tan[x], Cos[x]]]", "Sin[x]"),
        ("Simplify[Times[Sin[x], Power[Cos[x], -1]]]", "Tan[x]"),
        (
            "Simplify[Times[Plus[1, Power[Tan[x], 2]], Power[Cos[x], 2]]]",
            "1",
        ),
        (
            "Simplify[Plus[Sin[Plus[x, y]], Times[-1, Sin[x], Cos[y]], Times[-1, Cos[x], Sin[y]]]]",
            "0",
        ),
        (
            "Simplify[Plus[Power[Sin[x], -2], Times[-1, Power[Tan[x], -2]]]]",
            "1",
        ),
        (
            "Simplify[Plus[Cos[Times[2, x]], Times[2, Power[Sin[x], 2]]]]",
            "1",
        ),
        (
            "Simplify[Times[Power[Sin[x], 2], Power[Plus[1, Times[-1, Cos[x]]], -1]]]",
            "Plus[1, Cos[x]]",
        ),
        (
            "Simplify[Plus[Power[Sin[x], 3], Times[Sin[x], Power[Cos[x], 2]]]]",
            "Sin[x]",
        ),
        (
            "Simplify[Times[Cos[x], Power[Plus[1, Times[-1, Power[Sin[x], 2]]], -1]]]",
            "Power[Cos[x], -1]",
        ),
        (
            "Simplify[Plus[Times[Sin[x], Tan[x]], Cos[x]]]",
            "Power[Cos[x], -1]",
        ),
        (
            "Simplify[Times[Tan[Times[2, x]], Plus[1, Times[-1, Power[Tan[x], 2]]]]]",
            "Times[2, Tan[x]]",
        ),
        (
            "Simplify[Plus[Times[Sin[x], Cos[y]], Times[Cos[x], Sin[y]]]]",
            "Sin[Plus[x, y]]",
        ),
    ]);
}

/// Factorials and binomials (CAS-235).
#[test]
fn factorials() {
    check(&[
        (
            "Simplify[Times[Factorial[n], Power[Factorial[Plus[-1, n]], -1]]]",
            "n",
        ),
        (
            "Simplify[Times[Factorial[Plus[1, n]], Power[Factorial[n], -1]]]",
            "Plus[1, n]",
        ),
        (
            "Simplify[Times[Factorial[n], Power[Factorial[Plus[1, n]], -1]]]",
            "Power[Plus[1, n], -1]",
        ),
        (
            "Simplify[Plus[Factorial[n], Times[-1, n, Factorial[Plus[-1, n]]]]]",
            "0",
        ),
        (
            "Simplify[Times[Plus[2, n], Power[Factorial[Plus[2, n]], -1]]]",
            "Power[Factorial[Plus[1, n]], -1]",
        ),
        (
            "Simplify[Times[Binomial[n, k], Power[Binomial[n, Plus[-1, k]], -1]]]",
            "Times[Power[k, -1], Plus[1, Times[-1, k], n]]",
        ),
        // A binomial alone stays one.
        ("Simplify[Binomial[n, k]]", "Binomial[n, k]"),
    ]);
}
