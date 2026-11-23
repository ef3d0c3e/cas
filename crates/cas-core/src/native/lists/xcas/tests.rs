//! Tests of Xcas's list commands (CAS-555 to CAS-557): giac's answers from
//! the phase 1.3 probes, and the calls its errors leave.

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

/// Each input in full form stays as it is written (giac's errors).
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let found = ctx.eval(&read(input)).unwrap();
        if found != read(input) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn positions() {
    check(&[
        ("FindPositions[1, List[3, x, 1, 2, 1, 3]]", "List[2, 4]"),
        ("FindPositions[5, List[1, 2]]", "List[]"),
        ("FindPositions[1, List[List[1, 2], 1]]", "List[1]"),
        ("FindPositions[List[1, 2, 1], 1]", "List[0, 2]"),
        ("FindPositions[1, List[1.0, 2, 1]]", "List[0, 2]"),
        ("FirstPosition[List[3, x, 1, 2, 1, 3], 1]", "2"),
        ("FirstPosition[List[0, 1, 3, 2, 4, 2, 5], 2]", "3"),
        ("FirstPosition[1, List[3, 1]]", "1"),
        ("FirstPosition[List[3, x, 1], x]", "1"),
    ]);
    stays(&[
        "FindPositions[x, y]",
        "FirstPosition[List[1, 2], 5]",
        "FirstPosition[List[1, 2, 1], 1, 1]",
    ]);
}

#[test]
fn counts() {
    check(&[
        ("CountEqual[1, List[-2, 1, 0, 1, 2, -3]]", "2"),
        ("CountLess[1, List[-2, -1, 0, 1, 2, 3]]", "3"),
        ("CountGreater[1, List[-2, -1, 0, 1, 2, 3]]", "2"),
        ("CountEqual[4, List[List[3, 4], List[4, 1]]]", "2"),
        (
            "CountEqual[4, List[List[3, 4], List[4, 1]], MatrixRow]",
            "List[1, 1]",
        ),
        (
            "CountEqual[1, List[List[1, 1], List[1, 2]], MatrixColumn]",
            "List[2, 1]",
        ),
        (
            "CountLess[2, List[List[1, 3], List[0, 5]], MatrixColumn]",
            "List[2, 0]",
        ),
        ("CountEqual[x, List[x, y, x]]", "2"),
        ("CountLess[x, List[1, 2]]", "0"),
        ("CountGreater[0, List[Sqrt[2], Times[-1, Pi], x]]", "1"),
        ("CountGreater[0, List[List[1, 2], List[3]]]", "3"),
        ("CountEqual[1, List[1.0, 1]]", "2"),
        ("CountEqual[1, 1]", "1"),
        ("CountEqual[1, List[]]", "0"),
    ]);
    stays(&["CountEqual[1]", "CountEqual[1, List[1, 2], MatrixRow]"]);
}

#[test]
fn differences_and_sums() {
    check(&[
        ("Differences[List[1, 4, 8, 9]]", "List[3, 4, 1]"),
        ("Differences[List[1, 8, 4, 9]]", "List[7, -4, 5]"),
        (
            "Differences[List[List[1, 2], List[4, 8]]]",
            "List[List[3, 6]]",
        ),
        (
            "Differences[List[x, y, z]]",
            "List[Plus[y, Times[-1, x]], Plus[z, Times[-1, y]]]",
        ),
        (
            "AccumulateHeadTail[List[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 2]",
            "List[3, 3, 4, 5, 6, 7, 17]",
        ),
        ("AccumulateHeadTail[List[1, 2, 3], 1, 1]", "List[1, 2, 3]"),
        (
            "AccumulateHeadTail[List[a, b, c, d], 2, 1]",
            "List[Plus[a, b], c, d]",
        ),
        (
            "Lengths[List[List[1, 2, 3], List[1, 2], List[1]]]",
            "List[3, 2, 1]",
        ),
        ("Lengths[List[]]", "List[]"),
    ]);
    stays(&[
        "Differences[List[1]]",
        "Differences[List[]]",
        "Differences[5]",
        "AccumulateHeadTail[List[1, 2, 3], 2, 2]",
        "AccumulateHeadTail[List[1, 2, 3, 4], 0, 0]",
        "AccumulateHeadTail[List[1, 2, 3], 1]",
        "AccumulateHeadTail[5, 1, 1]",
        "Lengths[List[1, List[1, 2]]]",
        "Lengths[5]",
    ]);
}

#[test]
fn builders() {
    check(&[
        ("MakeVector[1, 2, 3]", "List[1, 2, 3]"),
        ("MakeVector[1]", "List[1]"),
        ("MakeVector[]", "List[]"),
        ("MakeVector[List[1, 2]]", "List[1, 2]"),
        ("List[MakeSequence[List[1, 2, 3]]]", "List[1, 2, 3]"),
        ("MakeSequence[1]", "List[1]"),
        ("List[MakeSequence[1, 2]]", "List[1, 2]"),
        ("NewList[4]", "List[0, 0, 0, 0]"),
        ("NewList[0]", "List[]"),
        ("NewList[-1]", "List[0]"),
        (
            "Enumerate[List[a, b, c]]",
            "List[List[0, a], List[1, b], List[2, c]]",
        ),
        (
            "Enumerate[List[a, b, c], 1]",
            "List[List[1, a], List[2, b], List[3, c]]",
        ),
        (
            "Enumerate[List[a, b], x]",
            "List[List[x, a], List[Plus[1, x], b]]",
        ),
        ("Enumerate[List[]]", "List[]"),
        ("Copy[Plus[x, 1]]", "Plus[x, 1]"),
        ("Copy[List[List[1, 2], x]]", "List[List[1, 2], x]"),
    ]);
    stays(&[
        "NewList[x]",
        "NewList[3/2]",
        "NewList[2, 3]",
        "NewList[65537]",
        "Enumerate[x]",
    ]);
}
