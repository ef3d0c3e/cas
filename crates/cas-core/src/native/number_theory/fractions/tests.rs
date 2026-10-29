//! Tests of continued fractions (CAS-424): rationals, periods of quadratic
//! irrationals, numeric quotients, and the way back.

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
fn quotients() {
    check(&[
        ("ContinuedFraction[355/113]", "List[3, 7, 16]"),
        ("ContinuedFraction[-7/3]", "List[-3, 1, 2]"),
        ("ContinuedFraction[5]", "List[5]"),
        ("ContinuedFraction[355/113, 2]", "List[3, 7]"),
        ("ContinuedFraction[Power[2, 1/2]]", "List[1, List[2]]"),
        (
            "ContinuedFraction[Power[7, 1/2]]",
            "List[2, List[1, 1, 1, 4]]",
        ),
        (
            "ContinuedFraction[Plus[1/2, Times[1/2, Power[5, 1/2]]]]",
            "List[List[1]]",
        ),
        (
            "ContinuedFraction[Plus[1, Times[-1, Power[2, 1/2]]]]",
            "List[-1, 1, 1, List[2]]",
        ),
        // (1 + √5)/2 as it evaluates, a product (the oracle's finding).
        (
            "ContinuedFraction[Times[1/2, Plus[1, Power[5, 1/2]]]]",
            "List[List[1]]",
        ),
        (
            "ContinuedFraction[Times[-2/3, Plus[1, Power[3, 1/2]]]]",
            "ContinuedFraction[Plus[-2/3, Times[-2/3, Power[3, 1/2]]]]",
        ),
        // (1 + √2)/3, its form scaled so that Q divides N − P².
        (
            "ContinuedFraction[Plus[1/3, Times[1/3, Power[2, 1/2]]]]",
            "List[0, 1, List[4, 8]]",
        ),
        (
            "ContinuedFraction[Power[13, 1/2]]",
            "List[3, List[1, 1, 1, 1, 6]]",
        ),
        // A period of 532572 quotients, by Python's: past 2^16, it stays.
        (
            "ContinuedFraction[Power[1000000000039, 1/2]]",
            "ContinuedFraction[Power[1000000000039, 1/2]]",
        ),
        ("ContinuedFraction[Pi, 6]", "List[3, 7, 15, 1, 292, 1]"),
        ("ContinuedFraction[E, 8]", "List[2, 1, 2, 1, 1, 4, 1, 1]"),
        ("ContinuedFraction[Power[2, 1/2], 4]", "List[1, 2, 2, 2]"),
        // What stays: a number not rational or quadratic without n.
        ("ContinuedFraction[Pi]", "ContinuedFraction[Pi]"),
    ]);
}

#[test]
fn back() {
    check(&[
        ("FromContinuedFraction[List[3, 7, 16]]", "355/113"),
        ("FromContinuedFraction[List[1, List[2]]]", "Power[2, 1/2]"),
        // giac's form, its period started late.
        (
            "FromContinuedFraction[List[1, 2, List[2]]]",
            "Power[2, 1/2]",
        ),
        (
            "FromContinuedFraction[List[2, List[1, 1, 1, 4]]]",
            "Power[7, 1/2]",
        ),
        (
            "FromContinuedFraction[List[List[1]]]",
            "Plus[1/2, Times[1/2, Power[5, 1/2]]]",
        ),
        (
            "FromContinuedFraction[List[-1, 1, 1, List[2]]]",
            "Plus[1, Times[-1, Power[2, 1/2]]]",
        ),
        (
            "FromContinuedFraction[ContinuedFraction[Plus[1/3, Times[1/3, Power[2, 1/2]]]]]",
            "Plus[1/3, Times[1/3, Power[2, 1/2]]]",
        ),
    ]);
}
