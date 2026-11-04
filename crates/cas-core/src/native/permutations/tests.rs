//! Tests of the permutations (CAS-463): giac's examples, the cycle that
//! records the size, products of cycles and permutations, and what stays.

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
fn cycles() {
    check(&[
        // giac's [[0,1,3,5],[2,4]], [[3],[1,2]], [[2],[0,1]], [[0]] and
        // [[1,2]]: n − 1 fixed records the size.
        (
            "PermutationCycles[List[1, 3, 4, 5, 2, 0]]",
            "List[List[0, 1, 3, 5], List[2, 4]]",
        ),
        (
            "PermutationCycles[List[0, 2, 1, 3]]",
            "List[List[3], List[1, 2]]",
        ),
        (
            "PermutationCycles[List[1, 0, 2]]",
            "List[List[2], List[0, 1]]",
        ),
        ("PermutationCycles[List[0]]", "List[List[0]]"),
        ("PermutationCycles[List[0, 2, 1]]", "List[List[1, 2]]"),
        // giac's [0,3,4,5,2,1], [1,0,2] and [2,1,0].
        (
            "CyclesPermutation[List[List[1, 3, 5], List[2, 4]]]",
            "List[0, 3, 4, 5, 2, 1]",
        ),
        (
            "CyclesPermutation[List[List[2], List[0, 1]]]",
            "List[1, 0, 2]",
        ),
        ("CyclePermutation[List[0, 2]]", "List[2, 1, 0]"),
        // Cycles that meet stay.
        (
            "CyclesPermutation[List[List[0, 1], List[1, 2]]]",
            "CyclesPermutation[List[List[0, 1], List[1, 2]]]",
        ),
    ]);
}

#[test]
fn products_and_inverses() {
    check(&[
        // giac's [5,0,4,1,2,3], [5,3,1] and [1,2,0].
        (
            "InversePermutation[List[1, 3, 4, 5, 2, 0]]",
            "List[5, 0, 4, 1, 2, 3]",
        ),
        ("InverseCycle[List[1, 3, 5]]", "List[5, 3, 1]"),
        ("InverseCycle[List[0, 2, 1]]", "List[1, 2, 0]"),
        // giac's p1op2, p1oc2, c1op2 and c1oc2.
        (
            "PermutationProduct[List[1, 2, 0], List[2, 1, 0]]",
            "List[0, 2, 1]",
        ),
        (
            "PermutationCycleProduct[List[1, 2, 0], List[0, 2]]",
            "List[0, 2, 1]",
        ),
        (
            "CyclePermutationProduct[List[0, 2], List[1, 2, 0]]",
            "List[1, 0, 2]",
        ),
        ("CycleProduct[List[0, 1], List[1, 2]]", "List[1, 2, 0]"),
        // Not permutations.
        (
            "PermutationProduct[List[1, 1], List[0, 1]]",
            "PermutationProduct[List[1, 1], List[0, 1]]",
        ),
    ]);
}

#[test]
fn invariants() {
    check(&[
        // giac's 4, 2, 1, -1, 1 and 1, 0, 1.
        ("PermutationOrder[List[1, 3, 4, 5, 2, 0]]", "4"),
        ("PermutationOrder[List[1, 0, 2, 4, 3, 5]]", "2"),
        ("PermutationSignature[List[1, 3, 4, 5, 2, 0]]", "1"),
        ("PermutationSignature[List[0, 2, 1, 3]]", "-1"),
        ("PermutationSignature[List[0]]", "1"),
        ("PermutationQ[List[1, 3, 4, 5, 2, 0]]", "True"),
        ("PermutationQ[List[1, 1, 0]]", "False"),
        ("PermutationQ[List[]]", "True"),
        ("PermutationQ[x]", "False"),
        // giac's [[0,1,0],[0,0,1],[1,0,0]].
        (
            "PermutationMatrix[List[1, 2, 0]]",
            "List[List[0, 1, 0], List[0, 0, 1], List[1, 0, 0]]",
        ),
    ]);
}
