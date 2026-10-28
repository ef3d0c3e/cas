//! Tests of the Diophantine equations (CAS-425): giac's examples, linear
//! systems, products, sums of two squares and Pell's equations.

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
fn pairs_and_squares() {
    check(&[
        ("IntegerSolveLinear[21, 28, 14]", "List[-2, 2]"),
        ("IntegerSolveLinear[4, 6, 3]", "IntegerSolveLinear[4, 6, 3]"),
        ("PrimeSumOfSquares[13]", "List[3, 2]"),
        ("PrimeSumOfSquares[17]", "List[4, 1]"),
        ("PrimeSumOfSquares[2]", "List[1, 1]"),
        // Primes past 10⁹, their squares found by Python's enumeration.
        ("PrimeSumOfSquares[1000000009]", "List[31400, 3747]"),
        ("PrimeSumOfSquares[1000000021]", "List[31614, 745]"),
        ("PrimeSumOfSquares[7]", "PrimeSumOfSquares[7]"),
        ("PrimeSumOfSquares[21]", "PrimeSumOfSquares[21]"),
    ]);
}

#[test]
fn linear() {
    check(&[
        // giac's [x=-1+5*_Z0,y=2-3*_Z0].
        (
            "IntegerSolve[Equation[Plus[Times[3, x], Times[5, y]], 7]]",
            "List[Equation[x, Plus[-1, Times[5, n_0]]], Equation[y, Plus[2, Times[-3, n_0]]]]",
        ),
        (
            "IntegerSolve[List[Equation[Plus[x, y], 5], Equation[Plus[x, Times[-1, y]], 1]]]",
            "List[Equation[x, 3], Equation[y, 2]]",
        ),
        // No solution: 2x + 4y = 3.
        (
            "IntegerSolve[Equation[Plus[Times[2, x], Times[4, y]], 3]]",
            "List[]",
        ),
        (
            "IntegerSolve[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2]]]",
            "List[]",
        ),
    ]);
    // 6x + 9y + 15z = 3 has two parameters; each instance satisfies it.
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "IntegerSolve[Equation[Plus[Times[6, x], Times[9, y], Times[15, z]], 3]]",
        ))
        .unwrap();
    for (a, b) in [(0, 0), (1, -2), (-3, 5), (7, 4)] {
        let value = |k: usize| {
            format!("ReplaceAll[Part[{found:?}, {k}, 2], List[Rule[n_0, {a}], Rule[n_1, {b}]]]")
        };
        let sum = ctx
            .eval(&read(&format!(
                "Plus[Times[6, {}], Times[9, {}], Times[15, {}]]",
                value(1),
                value(2),
                value(3)
            )))
            .unwrap();
        assert_eq!(sum, Expr::from(3), "{found:?} at ({a}, {b})");
    }
}

#[test]
fn nonlinear() {
    check(&[
        (
            "IntegerSolve[Equation[Times[x, y], 6]]",
            "List[List[Equation[x, 1], Equation[y, 6]], List[Equation[x, 2], Equation[y, 3]], \
             List[Equation[x, 3], Equation[y, 2]], List[Equation[x, 6], Equation[y, 1]], \
             List[Equation[x, -1], Equation[y, -6]], List[Equation[x, -2], Equation[y, -3]], \
             List[Equation[x, -3], Equation[y, -2]], List[Equation[x, -6], Equation[y, -1]]]",
        ),
        (
            "Length[IntegerSolve[Equation[Plus[Power[x, 2], Power[y, 2]], 25]]]",
            "12",
        ),
        (
            "Part[IntegerSolve[Equation[Plus[Power[x, 2], Power[y, 2]], 25]], 1]",
            "List[Equation[x, -5], Equation[y, 0]]",
        ),
        (
            "IntegerSolve[Equation[Plus[Power[x, 2], Power[y, 2]], 3]]",
            "List[]",
        ),
    ]);
    // Pell: x² − 2y² = 1 from (3, 2), x² − 13 y² = 1 from (649, 180): at
    // n_0 = 1 the fundamental solution, at 2 another.
    let ctx = Context::new();
    for (d, first) in [(2, "List[3, 2]"), (13, "List[649, 180]"), (3, "List[2, 1]")] {
        let found = ctx
            .eval(&read(&format!(
                "IntegerSolve[Equation[Plus[Power[x, 2], Times[-{d}, Power[y, 2]]], 1]]"
            )))
            .unwrap();
        let at = |n: i64| {
            ctx.eval(&read(&format!(
                "Expand[ReplaceAll[List[Part[{found:?}, 1, 2], Part[{found:?}, 2, 2]], Rule[n_0, {n}]]]"
            )))
            .unwrap()
        };
        assert_eq!(at(1), read(first), "D = {d}: {found:?}");
        let second = at(2);
        let check = ctx
            .eval(&read(&format!(
                "Plus[Power[Part[{second:?}, 1], 2], Times[-{d}, Power[Part[{second:?}, 2], 2]]]"
            )))
            .unwrap();
        assert_eq!(check, Expr::from(1), "D = {d}: {second:?}");
    }
}
