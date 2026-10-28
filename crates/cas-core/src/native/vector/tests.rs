//! Tests of the differential operators (CAS-408): giac's examples and what
//! stays.

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
fn operators() {
    check(&[
        (
            "Grad[Plus[Times[2, Power[x, 2], y], Times[-1, x, Power[z, 3]]], List[x, y, z]]",
            "List[Plus[Times[4, x, y], Times[-1, Power[z, 3]]], Times[2, Power[x, 2]], \
             Times[-3, x, Power[z, 2]]]",
        ),
        (
            "Divergence[List[Plus[Power[x, 2], y], Plus[x, y, z], Plus[Power[z, 3], Power[x, 2]]], \
             List[x, y, z]]",
            "Plus[1, Times[2, x], Times[3, Power[z, 2]]]",
        ),
        (
            "Divergence[List[Times[x, y], Times[y, z]], List[x, y]]",
            "Plus[y, z]",
        ),
        (
            "Curl[List[Times[2, x, y], Times[x, z], Times[y, z]], List[x, y, z]]",
            "List[Plus[z, Times[-1, x]], 0, Plus[z, Times[-2, x]]]",
        ),
        (
            "Curl[List[Times[y, z], Times[x, z], Times[x, y]], List[x, y, z]]",
            "List[0, 0, 0]",
        ),
        (
            "Laplacian[Plus[Times[Power[x, 2], y], Power[z, 3]], List[x, y, z]]",
            "Plus[Times[2, y], Times[6, z]]",
        ),
        (
            "Laplacian[3]",
            "List[List[2, -1, 0], List[-1, 2, -1], List[0, -1, 2]]",
        ),
        (
            "Hessian[Times[x, y, z], List[x, y, z]]",
            "List[List[0, z, y], List[z, 0, x], List[y, x, 0]]",
        ),
        (
            "Hessian[Plus[Times[Power[x, 2], y], z], List[x, y, z]]",
            "List[List[Times[2, y], Times[2, x], 0], List[Times[2, x], 0, 0], List[0, 0, 0]]",
        ),
        // What stays: no variables, a variable that is not a symbol, a
        // field of another length, a plane curl, a list to Grad.
        ("Grad[Times[x, y]]", "Grad[Times[x, y]]"),
        (
            "Grad[Times[x, y], List[1, x]]",
            "Grad[Times[x, y], List[1, x]]",
        ),
        (
            "Grad[Times[x, y], List[x, x]]",
            "Grad[Times[x, y], List[x, x]]",
        ),
        (
            "Divergence[List[x], List[x, y]]",
            "Divergence[List[x], List[x, y]]",
        ),
        (
            "Curl[List[x, y], List[x, y]]",
            "Curl[List[x, y], List[x, y]]",
        ),
        (
            "Grad[List[x, y], List[x, y]]",
            "Grad[List[x, y], List[x, y]]",
        ),
        ("Laplacian[0]", "Laplacian[0]"),
    ]);
}
