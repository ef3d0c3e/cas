//! Tests of the global optima (CAS-416): giac's examples, infinite ends,
//! points, and what stays.

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
fn one_variable() {
    check(&[
        ("Minimize[Plus[Power[x, 2], Times[-2, x], 3], x]", "2"),
        (
            "Minimize[Plus[Power[x, 3], Times[-3, x]], Equation[x, Range[-3, 3]]]",
            "-18",
        ),
        (
            "Maximize[Plus[Power[x, 3], Times[-3, x]], Equation[x, Range[-3, 3]]]",
            "18",
        ),
        (
            "Minimize[Plus[Power[x, 4], Times[-1, Power[x, 2]]], Equation[x, Range[-3, 3]], \
             coordinates]",
            "List[-1/4, List[List[Times[-1/2, Power[2, 1/2]]], List[Times[1/2, Power[2, 1/2]]]]]",
        ),
        // giac's is 0, its critical value; the infimum is -∞.
        ("Minimize[Power[x, 3], x]", "Times[-1, Infinity]"),
        ("Maximize[Power[E, Times[-1, Power[x, 2]]], x]", "1"),
        // An infimum not attained, with a warning.
        ("Minimize[Power[E, Times[-1, Power[x, 2]]], x]", "0"),
        ("FMin[Plus[Power[x, 2], Times[-2, x], 1], x]", "1"),
        ("FMax[Plus[Times[-1, Power[x, 2]], Times[4, x]], x]", "2"),
        (
            "FMin[Plus[Power[x, 2], Times[-2, x], 1], Equation[x, Range[2, 3]]]",
            "2",
        ),
    ]);
}

#[test]
fn several_variables() {
    check(&[
        (
            "Minimize[Plus[Power[x, 2], Times[-3, x], Power[y, 2], Times[3, y], 3], \
             List[Equation[x, Range[2, 4]], Equation[y, Range[-4, -2]]], point]",
            "List[-1, List[List[2, -2]]]",
        ),
        (
            "Minimize[Plus[Times[2, Power[x, 2]], Times[-1, Power[y, 2]], Times[6, y]], \
             LessEqual[Plus[Power[x, 2], Power[y, 2]], 16], List[x, y]]",
            "-40",
        ),
        (
            "Maximize[Plus[Times[2, Power[x, 2]], Times[-1, Power[y, 2]], Times[6, y]], \
             LessEqual[Plus[Power[x, 2], Power[y, 2]], 16], List[x, y]]",
            "35",
        ),
        // Coercive: the line is unbounded, x² + y²'s level sets are not.
        (
            "Minimize[Plus[Power[x, 2], Power[y, 2]], Equation[Plus[x, y], 1], List[x, y]]",
            "1/2",
        ),
        // What stays: a region not shown bounded.
        (
            "Maximize[Times[x, y], Equation[Plus[x, y], 1], List[x, y]]",
            "Maximize[Times[x, y], Equation[Plus[x, y], 1], List[x, y]]",
        ),
    ]);
}
