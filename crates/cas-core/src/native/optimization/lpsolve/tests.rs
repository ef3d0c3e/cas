//! Tests of the linear programs (CAS-417, CAS-418): giac's examples,
//! unbounded and infeasible programs, integers and binaries.

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
fn programs() {
    check(&[
        (
            "LPSolve[Plus[Times[2, x], y, Times[-1, z], 4], List[LessEqual[x, 1], \
             GreaterEqual[y, 2], Equation[Plus[x, Times[3, y], Times[-1, z]], 2], \
             LessEqual[Plus[Times[2, x], Times[-1, y], z], 8], \
             LessEqual[Plus[Times[-1, x], y], 5]]]",
            "List[-4, List[Equation[x, 0], Equation[y, 5], Equation[z, 13]]]",
        ),
        (
            "LPSolve[Plus[Times[-7, x], Times[2, y]], List[LessEqual[Plus[Times[4, x], \
             Times[-12, y]], 20], LessEqual[Plus[Times[-1, x], Times[3, y]], 3]], \
             Equation[x, Range[-5, 5]], Equation[y, Range[0, Infinity]], Maximize]",
            "List[21, List[Equation[x, -3], Equation[y, 0]]]",
        ),
        (
            "LPSolve[Plus[Times[-1, x], Times[-1, y]], List[LessEqual[y, Plus[Times[3, x], 1/2]], \
             LessEqual[y, Plus[Times[-5, x], 2]]], Equation[Assume, lp_nonnegative]]",
            "List[-5/4, List[Equation[x, 3/16], Equation[y, 17/16]]]",
        ),
        (
            "LPSolve[Plus[Times[-6, x], Times[4, y], z], List[LessEqual[Plus[Times[5, x], \
             Times[-10, y]], 20], Equation[Plus[Times[2, z], Times[-3, y]], 6], \
             LessEqual[Plus[Times[-1, x], Times[3, y]], 3]], Equation[x, Range[1, 20]], \
             Equation[y, Range[0, Infinity]]]",
            "List[-133/2, List[Equation[x, 18], Equation[y, 7], Equation[z, 27/2]]]",
        ),
        // Unbounded and infeasible, as giac's.
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[Plus[x, y], 2]], \
             Equation[Assume, lp_nonnegative], Maximize]",
            "List[Infinity, List[]]",
        ),
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[Plus[x, y], 2], LessEqual[Plus[x, y], 1]]]",
            "List[]",
        ),
    ]);
}

#[test]
fn integers() {
    check(&[
        (
            "LPSolve[Plus[Times[-5, x], Times[-7, y]], List[LessEqual[Plus[Times[7, x], y], 35], \
             LessEqual[Plus[Times[-1, x], Times[3, y]], 6]], Equation[Assume, integer]]",
            "List[-41, List[Equation[x, 4], Equation[y, 3]]]",
        ),
        (
            "LPSolve[Plus[Times[8, x1], Times[11, x2], Times[6, x3], Times[4, x4]], \
             List[LessEqual[Plus[Times[5, x1], Times[7, x2], Times[4, x3], Times[3, x4]], 14]], \
             Equation[Assume, lp_binary], Maximize]",
            "List[21, List[Equation[x1, 0], Equation[x2, 1], Equation[x3, 1], Equation[x4, 1]]]",
        ),
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[x, 1/3], GreaterEqual[y, 1/2]], \
             Equation[Assume, integer]]",
            "List[2, List[Equation[x, 1], Equation[y, 1]]]",
        ),
        // What stays: a nonlinear objective.
        (
            "LPSolve[Times[x, y], List[LessEqual[x, 1]]]",
            "LPSolve[Times[x, y], List[LessEqual[x, 1]]]",
        ),
    ]);
}
