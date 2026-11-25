//! Tests of `IntegerPart`, `FractionalPart` and `Truncate` (CAS-622).

use std::cmp::Ordering;

use crate::{Context, Expr, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
#[test]
fn exact() {
    let ctx = Context::new();
    let cases = [
        ("IntegerPart[7/2]", "3"),
        ("IntegerPart[-7/2]", "-3"),
        ("IntegerPart[Power[2, 1/2]]", "1"),
        ("IntegerPart[List[7/2, -7/2]]", "List[3, -3]"),
        ("List[IntegerPart[7/2, Power[2, 1/2]]]", "List[3, 1]"),
        ("IntegerPart[x]", "IntegerPart[x]"),
        ("IntegerPart[]", "IntegerPart[]"),
        ("FractionalPart[1/2]", "1/2"),
        ("FractionalPart[-1/2]", "-1/2"),
        ("FractionalPart[7/3]", "1/3"),
        ("FractionalPart[Power[2, 1/2]]", "Plus[-1, Power[2, 1/2]]"),
        ("FractionalPart[x]", "FractionalPart[x]"),
        ("Truncate[7/2]", "3"),
        ("Truncate[-7/2]", "-3"),
        ("Truncate[-0x1.3p+2]", "-4"),
        ("Truncate[Power[2, 1/2], 3]", "707/500"),
        ("Truncate[Times[-1, Power[2, 1/2]], 3]", "-707/500"),
        ("Truncate[Pi, 4]", "6283/2000"),
        ("Truncate[1/3, 2]", "33/100"),
        ("Truncate[2/3, 1]", "3/5"),
        ("Truncate[1/3, 0]", "0"),
        ("Truncate[123, -1]", "100"),
        (
            "Truncate[List[Power[2, 1/2], 7/2], 2]",
            "List[141/100, 7/2]",
        ),
        (
            "Truncate[Plus[Power[2, 1/2], Times[I, Power[5, 1/2]]], 4]",
            "Plus[7071/5000, Times[I, 559/250]]",
        ),
        ("FractionalPart[Complex[7/2, -1/3]]", "Complex[1/2, -1/3]"),
        ("Truncate[x]", "Truncate[x]"),
        ("Truncate[x, 2]", "Truncate[x, 2]"),
        ("Truncate[1/3, y]", "Truncate[1/3, y]"),
        // `10^b` past the limits stays.
        ("Truncate[1/3, 1000000000]", "Truncate[1/3, 1000000000]"),
        (
            "Truncate[1/3, -1000000000000]",
            "Truncate[1/3, -1000000000000]",
        ),
        ("Truncate[0, -2]", "Truncate[0, -2]"),
    ];
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

/// A decimal gives a decimal, compared by its midpoint (C0028).
#[test]
fn decimals() {
    let ctx = Context::new();
    let cases = [
        ("IntegerPart[0x1.1p+2]", 4.0),
        ("IntegerPart[-0x1.3p+2]", -4.0),
        ("FractionalPart[0x1.4p+0]", 0.25),
        ("FractionalPart[-0x1.4p+1]", -0.5),
        ("Truncate[0x1.1p+2, 2]", 4.25),
        ("Truncate[0x1.81cp+13, -2]", 12000.0),
    ];
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let ok = match found.view() {
            View::Ball(b) => b.mid().to_f64().total_cmp(&expected) == Ordering::Equal,
            _ => false,
        };
        if !ok {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
