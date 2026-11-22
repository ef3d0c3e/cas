//! Tests of the polynomial tools III (CAS-540, CAS-543): giac's answers.

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

/// Inputs that leave the call: the call of its evaluated (canonically
/// ordered) arguments.
fn unevaluated(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        let head = e.head().expect("a call").clone();
        let args: Vec<Expr> = e.args().iter().map(|a| ctx.eval(a).unwrap()).collect();
        if found != Expr::apply(&head, &args) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn root_of() {
    check(&[
        // Polynomials in one variable become giac's lists.
        (
            "RootOf[Plus[x, 2], Plus[Power[x, 2], -2]]",
            "RootOf[List[1, 2], List[1, 0, -2]]",
        ),
        (
            "RootOf[Plus[y, 2], Plus[Power[y, 2], -2]]",
            "RootOf[List[1, 2], List[1, 0, -2]]",
        ),
        (
            "RootOf[Plus[Power[x, 4], x, 1]]",
            "RootOf[List[1, 0], List[1, 0, 0, 1, 1]]",
        ),
        (
            "RootOf[List[1, 0], Plus[Power[x, 2], -2]]",
            "RootOf[List[1, 0], List[1, 0, -2]]",
        ),
    ]);
    unevaluated(&[
        "RootOf[List[1, 0], List[1, 0, 10, 0, 1]]",
        "RootOf[List[], List[1, 0, -2]]",
        // Not monic: giac errs.
        "RootOf[List[1, 0], List[2, 0, -4]]",
        // Two variables.
        "RootOf[Plus[x, 2], Plus[Power[y, 2], -2]]",
    ]);
}

#[test]
fn factored_from_roots() {
    check(&[
        (
            "FactoredFromRoots[List[1, 2, -1, 1]]",
            "Times[Power[Plus[x, -1], 2], Plus[x, 1]]",
        ),
        (
            "FactoredFromRoots[List[1, 2, 0, 1, 3, -1]]",
            "Times[Power[Plus[x, -1], 2], x, Power[Plus[x, -3], -1]]",
        ),
        ("FactoredFromRoots[List[]]", "1"),
        ("FactoredFromRoots[List[1, 0]]", "1"),
        ("FactoredFromRoots[List[2, 3], y]", "Power[Plus[y, -2], 3]"),
    ]);
    unevaluated(&["FactoredFromRoots[List[1, 2, 3]]"]);
}

#[test]
fn homogenize() {
    check(&[
        (
            "Homogenize[Plus[Power[x, 2], Times[-1, y]], z]",
            "Plus[Power[x, 2], Times[-1, y, z]]",
        ),
        (
            "Homogenize[Plus[Power[x, 2], -1]]",
            "Plus[Power[x, 2], Times[-1, Power[t, 2]]]",
        ),
        (
            "Homogenize[Plus[Times[x, y], -1], y]",
            "Plus[Times[x, y], Times[-1, Power[y, 2]]]",
        ),
        (
            "Homogenize[Plus[Power[x, 3], Times[x, y], 1], z]",
            "Plus[Power[x, 3], Times[x, y, z], Power[z, 3]]",
        ),
        (
            "Homogenize[Plus[Power[x, 2], Power[y, 3]]]",
            "Plus[Times[Power[x, 2], t], Power[y, 3]]",
        ),
        (
            "Homogenize[Plus[Power[x, 2], -1], Power[z, 2]]",
            "Plus[Power[x, 2], Times[-1, Power[z, 4]]]",
        ),
        ("Homogenize[5]", "5"),
        (
            "Homogenize[Plus[Power[x, -1], 1]]",
            "Times[Plus[x, t], Power[x, -1]]",
        ),
    ]);
    unevaluated(&["Homogenize[List[Plus[Power[x, 2], -1], x], z]"]);
}

#[test]
fn total_degree() {
    check(&[
        (
            "TotalDegree[Plus[Times[Power[x, 2], y], y], List[x, y]]",
            "3",
        ),
        (
            "TotalDegree[Plus[Times[Power[x, 3], y], Times[x, y]], List[x, y]]",
            "4",
        ),
        (
            "TotalDegree[Plus[Times[Power[x, 2], y, z], y], List[x, y]]",
            "3",
        ),
        ("TotalDegree[Plus[Times[Power[x, 2], y], y]]", "2"),
        ("TotalDegree[Plus[Times[Power[x, 2], y], y], x]", "2"),
        ("TotalDegree[Plus[Times[Power[x, 2], y], y], List[z]]", "0"),
        ("TotalDegree[0, List[x, y]]", "0"),
        ("TotalDegree[7, List[x, y]]", "0"),
    ]);
    unevaluated(&["TotalDegree[Plus[Power[x, -1], y], List[x, y]]"]);
}

#[test]
fn factor_xn() {
    check(&[
        (
            "FactorXn[Plus[Power[x, 3], Power[x, 2]]]",
            "Times[Power[x, 3], Plus[1, Power[x, -1]]]",
        ),
        (
            "FactorXn[Plus[Power[x, 4], -1]]",
            "Times[Power[x, 4], Plus[1, Times[-1, Power[x, -4]]]]",
        ),
        (
            "FactorXn[Plus[Times[2, x], 1]]",
            "Times[x, Plus[2, Power[x, -1]]]",
        ),
        (
            "FactorXn[Plus[Times[Power[x, 2], y], y]]",
            "Times[Power[x, 2], Plus[y, Times[y, Power[x, -2]]]]",
        ),
        ("FactorXn[Plus[Power[y, 2], y]]", "Plus[Power[y, 2], y]"),
        ("FactorXn[5]", "5"),
        ("FactorXn[0]", "0"),
    ]);
    unevaluated(&["FactorXn[Plus[Power[x, -1], 1]]"]);
}

#[test]
fn irreducible_polynomial_q() {
    check(&[
        ("IrreduciblePolynomialQ[Plus[Power[x, 2], 1]]", "1"),
        ("IrreduciblePolynomialQ[Plus[Power[x, 2], -1]]", "0"),
        // Over ℚ, though Factor splits it with radicals.
        ("IrreduciblePolynomialQ[Plus[Power[x, 2], -2]]", "1"),
        ("IrreduciblePolynomialQ[Plus[Power[x, 4], 4]]", "0"),
        (
            "IrreduciblePolynomialQ[Plus[Times[4, Power[x, 2]], -4]]",
            "0",
        ),
        ("IrreduciblePolynomialQ[Power[x, 2]]", "0"),
        ("IrreduciblePolynomialQ[Plus[Times[2, x], 2]]", "1"),
        (
            "IrreduciblePolynomialQ[Plus[Times[Power[x, 2], Power[y, 2]], -1]]",
            "0",
        ),
        (
            "IrreduciblePolynomialQ[Plus[Power[x, 2], Power[y, 2]]]",
            "1",
        ),
        (
            "IrreduciblePolynomialQ[Plus[Times[Power[x, 2], y], 1]]",
            "1",
        ),
        ("IrreduciblePolynomialQ[Plus[Power[x, 2], 1], x]", "1"),
        // giac's 1: x y has two factors.
        ("IrreduciblePolynomialQ[Times[x, y]]", "0"),
        ("IrreduciblePolynomialQ[5]", "1"),
        ("IrreduciblePolynomialQ[0]", "1"),
    ]);
    unevaluated(&["IrreduciblePolynomialQ[Power[x, -1]]"]);
}
