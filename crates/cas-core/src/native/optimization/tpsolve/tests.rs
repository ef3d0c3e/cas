//! Tests of the transportation problems (CAS-598): giac's answers, the
//! optimum against `LPSolve`, and the inputs that leave the call.

use crate::{Context, Expr, Symbol};

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

/// giac's answers to its probe rows.
#[test]
fn problems() {
    check(&[
        (
            "List[TPSolve[List[10, 20], List[15, 15], List[List[1, 2], List[3, 1]]]]",
            "List[40, List[List[10, 0], List[5, 15]]]",
        ),
        (
            "List[TPSolve[List[10, 20], List[10, 10], List[List[1, 2], List[3, 1]]]]",
            "List[20, List[List[10, 0], List[0, 10]]]",
        ),
        (
            "List[TPSolve[List[12, 17, 11], List[10, 10, 10, 10], List[List[500, 750, 300, 450], \
             List[650, 800, 400, 600], List[400, 700, 500, 550]]]]",
            "List[20200, List[List[0, 0, 2, 10], List[0, 9, 8, 0], List[10, 1, 0, 0]]]",
        ),
        (
            "List[TPSolve[List[95, 70, 165, 165], List[195, 150, 30, 45, 75], \
             List[List[15, M, 45, M, 0], List[12, 40, M, M, 0], List[0, 15, 25, 25, 0], \
             List[M, 0, M, 12, 0]]]]",
            "List[2820, List[List[20, 0, 0, 0, 75], List[70, 0, 0, 0, 0], \
             List[105, 0, 30, 30, 0], List[0, 150, 0, 15, 0]]]",
        ),
        (
            "List[TPSolve[List[1, 1, 1, 1], List[1, 1, 1, 1], List[List[10, 12, 9, 11], \
             List[5, 10, 7, 8], List[12, 14, 13, 11], List[8, 15, 11, 9]]]]",
            "List[37, List[List[0, 0, 1, 0], List[1, 0, 0, 0], List[0, 1, 0, 0], \
             List[0, 0, 0, 1]]]",
        ),
    ]);
}

/// The cost is the optimum of the linear program.
#[test]
fn optimum() {
    check(&[
        (
            "First[List[TPSolve[List[10, 20], List[15, 15], List[List[1, 2], List[3, 1]]]]]",
            "First[LPSolve[Plus[Times[1, x11], Times[2, x12], Times[3, x21], Times[1, x22]], \
             List[Equation[Plus[x11, x12], 10], Equation[Plus[x21, x22], 20], \
             Equation[Plus[x11, x21], 15], Equation[Plus[x12, x22], 15]], \
             Equation[Assume, lp_nonnegative]]]",
        ),
        (
            "First[List[TPSolve[List[12, 17, 11], List[10, 10, 10, 10], \
             List[List[500, 750, 300, 450], List[650, 800, 400, 600], \
             List[400, 700, 500, 550]]]]]",
            "First[LPSolve[Plus[Times[500, x11], Times[750, x12], Times[300, x13], \
             Times[450, x14], Times[650, x21], Times[800, x22], Times[400, x23], \
             Times[600, x24], Times[400, x31], Times[700, x32], Times[500, x33], \
             Times[550, x34]], List[Equation[Plus[x11, x12, x13, x14], 12], \
             Equation[Plus[x21, x22, x23, x24], 17], Equation[Plus[x31, x32, x33, x34], 11], \
             Equation[Plus[x11, x21, x31], 10], Equation[Plus[x12, x22, x32], 10], \
             Equation[Plus[x13, x23, x33], 10], Equation[Plus[x14, x24, x34], 10]], \
             Equation[Assume, lp_nonnegative]]]",
        ),
    ]);
}

/// Inputs that leave the call: a decimal cost, a negative supply, `2 M`,
/// two symbols, a wrong shape, a wrong arity; `M` in the plan and a zero
/// supply (its basis does not span, as in giac's method) warn and leave.
#[test]
fn left() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "TPSolve[List[10, 20], List[15, 15], List[List[N[3/2], 2], List[3, 1]]]",
        "TPSolve[List[10, -20], List[15, 15], List[List[1, 2], List[3, 1]]]",
        "TPSolve[List[10, 20], List[15, 15], List[List[Times[2, M], 2], List[3, 1]]]",
        "TPSolve[List[10, 20], List[15, 15], List[List[M, 2], List[3, K]]]",
        "TPSolve[List[10, 20], List[15, 15], List[List[1, 2]]]",
        "TPSolve[List[], List[], List[]]",
        "TPSolve[List[10, 20], List[15, 15]]",
        "TPSolve[List[6], List[6], List[List[M]]]",
        "TPSolve[List[0, 10], List[5, 5], List[List[1, 2], List[3, 1]]]",
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        if found.head().and_then(Expr::as_symbol) != Some(Symbol::TP_SOLVE) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
