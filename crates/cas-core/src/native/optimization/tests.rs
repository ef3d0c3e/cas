//! Tests of the optimization builtins (CAS-415 …): giac's examples, higher
//! derivatives, saddles, constraints, and what stays.

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
fn extrema() {
    check(&[
        (
            "Extrema[Plus[Power[x, 3], Times[-3, x]], x]",
            "List[List[1], List[-1]]",
        ),
        (
            "Extrema[Plus[Power[x, 2], Power[y, 2], Times[-2, x]], List[x, y]]",
            "List[List[List[1, 0]], List[]]",
        ),
        // A saddle at the origin and a minimum at (1, 1).
        (
            "Extrema[Plus[Power[x, 3], Power[y, 3], Times[-3, x, y]], List[x, y]]",
            "List[List[List[1, 1]], List[]]",
        ),
        // Higher derivatives: x⁴ has a minimum, x³ nothing.
        ("Extrema[Power[x, 4], x]", "List[List[0], List[]]"),
        ("Extrema[Power[x, 3], x]", "List[List[], List[]]"),
        (
            "Extrema[Times[-1, Power[x, 6]], x]",
            "List[List[], List[0]]",
        ),
        // A range keeps the points within it.
        (
            "Extrema[Plus[Power[x, 3], Times[-3, x]], Equation[x, Range[0, Infinity]]]",
            "List[List[1], List[]]",
        ),
        // Lagrange's multipliers: giac's examples.
        (
            "Extrema[Times[x, y, z], Equation[Plus[x, y, z], 1], List[x, y, z]]",
            "List[List[], List[List[1/3, 1/3, 1/3]]]",
        ),
        (
            "Extrema[Plus[Power[x, 2], Power[y, 2]], Equation[Times[x, y], 1], \
             List[Equation[x, Range[0, Infinity]], Equation[y, Range[0, Infinity]]]]",
            "List[List[List[1, 1]], List[]]",
        ),
        (
            "Extrema[Times[x, y], Equation[Plus[x, y], 1], List[x, y]]",
            "List[List[], List[List[1/2, 1/2]]]",
        ),
        // What stays: infinitely many critical points.
        (
            "Extrema[Times[Power[x, 2], Power[y, 2]], List[x, y]]",
            "Extrema[Times[Power[x, 2], Power[y, 2]], List[x, y]]",
        ),
        ("Extrema[x, 1]", "Extrema[x, 1]"),
    ]);
}
