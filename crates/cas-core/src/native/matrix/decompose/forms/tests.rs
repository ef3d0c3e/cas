//! Tests of the Hessenberg form and the constructions (CAS-476): giac's
//! answers, the similarity, and the resultant.

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
fn hessenberg() {
    let a = "List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 10]]";
    check(&[
        // giac's P and H.
        (
            &format!("HessenbergDecomposition[{a}]"),
            "List[List[List[1, 0, 0], List[0, 4/7, 1], List[0, 1, 0]], \
             List[List[1, 29/7, 2], List[7, 102/7, 8], List[0, 26/49, 3/7]]]",
        ),
        // P H P⁻¹ = A, with a symbol.
        (
            "Together[Times[Part[HessenbergDecomposition[List[List[a, 1, 2], List[1, 0, 1], \
             List[2, 1, 3]]], 1], Part[HessenbergDecomposition[List[List[a, 1, 2], List[1, 0, 1], \
             List[2, 1, 3]]], 2], Inverse[Part[HessenbergDecomposition[List[List[a, 1, 2], \
             List[1, 0, 1], List[2, 1, 3]]], 1]]]]",
            "List[List[a, 1, 2], List[1, 0, 1], List[2, 1, 3]]",
        ),
        // Already Hessenberg: P = I.
        (
            "HessenbergDecomposition[List[List[1, 2], List[3, 4]]]",
            "List[List[List[1, 0], List[0, 1]], List[List[1, 2], List[3, 4]]]",
        ),
    ]);
}

#[test]
fn constructions() {
    check(&[
        // giac's [[0,0,-5],[1,0,-2],[0,1,0]] and [[0,-3/2],[1,0]].
        (
            "CompanionMatrix[Plus[Power[x, 3], Times[2, x], 5], x]",
            "List[List[0, 0, -5], List[1, 0, -2], List[0, 1, 0]]",
        ),
        (
            "CompanionMatrix[Plus[Times[2, Power[x, 2]], 3], x]",
            "List[List[0, -3/2], List[1, 0]]",
        ),
        // giac's sylvester matrices; the determinant is the resultant.
        (
            "SylvesterMatrix[Plus[Power[x, 2], 1], Plus[x, 2], x]",
            "List[List[1, 0, 1], List[1, 2, 0], List[0, 1, 2]]",
        ),
        (
            "SylvesterMatrix[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]",
            "List[List[1, 0, 1, 1, 0], List[0, 1, 0, 1, 1], List[2, 0, 3, 0, 0], \
             List[0, 2, 0, 3, 0], List[0, 0, 2, 0, 3]]",
        ),
        (
            "Det[SylvesterMatrix[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]]",
            "Resultant[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]",
        ),
        // giac's [[5,12],[21,32]].
        (
            "HadamardProduct[List[List[1, 2], List[3, 4]], List[List[5, 6], List[7, 8]]]",
            "List[List[5, 12], List[21, 32]]",
        ),
        (
            "HadamardProduct[List[x, 2], List[y, 3]]",
            "List[Times[x, y], 6]",
        ),
    ]);
}
