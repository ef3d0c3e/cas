//! Tests of Xcas's list commands (CAS-555 to CAS-560): giac's answers from
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
        ("FindPositions[1, List[0x1p+0, 2, 1]]", "List[0, 2]"),
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
        ("CountEqual[1, List[0x1p+0, 1]]", "2"),
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

#[test]
fn maps_and_folds() {
    check(&[
        (
            "ApplyEach[Function[List[x], Power[x, 3]], List[1, 2, 3]]",
            "List[1, 8, 27]",
        ),
        (
            "ApplyEach[Function[List[x], Plus[x, 1]], List[List[1, 2, 3], List[1, 2, 3]], Matrix]",
            "List[List[2, 3, 4], List[2, 3, 4]]",
        ),
        ("ApplyEach[f, List[1, 2]]", "List[f[1], f[2]]"),
        ("ApplyEach[f, List[], Matrix]", "List[]"),
        ("FoldLeft[F, i, a, b, c]", "F[F[F[i, a], b], c]"),
        ("FoldRight[F, i, a, b, c]", "F[a, F[b, F[c, i]]]"),
        ("FoldLeft[F, i, a]", "F[i, a]"),
        (
            "FoldLeft[Function[List[x, y], Plus[x, y]], 0, 1, 2, 3]",
            "6",
        ),
        (
            "FoldRight[Function[List[x, y], Plus[x, Times[-1, y]]], 0, 1, 2, 3]",
            "2",
        ),
    ]);
    stays(&[
        "ApplyEach[f, x]",
        "ApplyEach[f, List[1, 2], Matrix]",
        "ApplyEach[f, List[List[1, 2], List[3]], Matrix]",
        "ApplyEach[f, List[1, 2], List[3, 4]]",
        "FoldLeft[F]",
        "FoldLeft[F, i]",
        "FoldRight[]",
    ]);
}

#[test]
fn products() {
    check(&[
        ("Multiply[List[1, 2], List[3, 4, 5]]", "11"),
        ("Multiply[List[4, 1], List[-4, 2]]", "-14"),
        (
            "Multiply[List[List[4, 1], List[-4, 1]], List[List[4, 1], List[-4, 1]]]",
            "List[List[12, 5], List[-20, -3]]",
        ),
        (
            "Multiply[List[1, 2], List[List[1, 0], List[0, 1]]]",
            "List[1, 2]",
        ),
        (
            "Multiply[List[List[1, 2], List[3, 4]], List[5, 6]]",
            "List[17, 39]",
        ),
        ("Multiply[2, List[1, 2]]", "List[2, 4]"),
        (
            "Multiply[List[1, 2], List[3, 4], List[5, 6]]",
            "List[55, 66]",
        ),
        ("Multiply[2, 3, 4]", "24"),
        ("Multiply[List[1, 2]]", "2"),
        ("Multiply[]", "1"),
        ("Multiply[x, y]", "Times[x, y]"),
    ]);
    stays(&["Multiply[List[List[1, 2], List[3, 4]], List[List[1, 2, 3]]]"]);
}

#[test]
fn insertions_and_substitutions() {
    check(&[
        ("InsertAt[List[0, 1, 2, 3], 2, 5]", "List[0, 1, 5, 2, 3]"),
        ("InsertAt[List[0, 1, 2, 3], 0, 5]", "List[5, 0, 1, 2, 3]"),
        ("InsertAt[List[0, 1, 2, 3], 4, 5]", "List[0, 1, 2, 3, 5]"),
        ("InsertAt[List[], 0, 5]", "List[5]"),
        (
            "SubstituteOperands[List[0, 1, 2, 3, 4], Equation[2, Plus[Sqrt[2], 1]]]",
            "List[0, 1, Plus[Sqrt[2], 1], 3, 4]",
        ),
        (
            "SubstituteOperands[List[List[1, 2], List[3, 4]], Equation[List[1, 1], 5]]",
            "List[List[1, 2], List[3, 5]]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[1, Null]]",
            "List[0, 2, 3]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[0, Null], Equation[1, Null]]",
            "List[2, 3]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[1, Null], Equation[1, 9]]",
            "List[0, 9, 2, 3]",
        ),
        (
            "SubstituteOperands[List[List[1, 2], List[3, 4]], Equation[List[1, 1], Null]]",
            "List[List[1, 2], List[3]]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[-1, 9]]",
            "List[0, 1, 2, 9]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], 1, 9]",
            "List[0, 9, 2, 3]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[x, 9], a]",
            "List[0, 1, 2, 3]",
        ),
        (
            "SubstituteOperands[List[0, 1, 2, 3], Equation[List[1, 0], 9]]",
            "List[0, 1, 2, 3]",
        ),
    ]);
    stays(&[
        "InsertAt[List[0, 1, 2, 3], 5, 5]",
        "InsertAt[List[0, 1, 2, 3], -1, 5]",
        "InsertAt[List[0, 1], x, 5]",
        "InsertAt[x, 0, 1]",
        "InsertAt[List[0, 1], 0]",
        "SubstituteOperands[List[0, 1, 2, 3], Equation[4, 9]]",
        "SubstituteOperands[List[0, 1, 2, 3], Equation[-5, 9]]",
        "SubstituteOperands[List[0, 1, 2, 3], Equation[List[], 9]]",
        "SubstituteOperands[List[List[1, 2], List[3, 4]], Equation[List[1, 5], 9]]",
        "SubstituteOperands[List[0, 1, 2, 3]]",
        "SubstituteOperands[5, Equation[1, 2]]",
    ]);
}

/// libstdc++'s heaps: giac's lists after `heapify` and `heappop` in place,
/// its `heappush` answers (CAS-560).
#[test]
fn heaps() {
    check(&[
        (
            "Heapify[List[3, 1, 4, 1, 5, 9, 2, 6]]",
            "List[9, 6, 4, 1, 5, 3, 2, 1]",
        ),
        (
            "Heapify[List[1, 2, 3, 4, 5, 6, 7]]",
            "List[7, 5, 6, 4, 2, 1, 3]",
        ),
        ("Heapify[List[]]", "List[]"),
        (
            "Heapify[List[1, 2, 3], Function[List[a, b], Greater[a, b]]]",
            "List[1, 2, 3]",
        ),
        ("HeapPop[List[5, 3, 4, 1]]", "List[4, 3, 1]"),
        (
            "HeapPop[List[9, 7, 5, 6, 1, 4, 2, 3, 1]]",
            "List[7, 6, 5, 3, 1, 4, 2, 1]",
        ),
        ("HeapPop[List[7]]", "List[]"),
        // libgiac crashes here; the port of pop_heap on a list that is not
        // a heap.
        (
            "HeapPop[List[9, 6, 5, 1, 1, 4, 2, 3]]",
            "List[6, 3, 5, 1, 1, 4, 2]",
        ),
        (
            "HeapPush[List[9, 6, 5, 1, 1, 4, 2, 3], 7]",
            "List[9, 7, 5, 6, 1, 4, 2, 3, 1]",
        ),
        (
            "HeapPush[List[9, 6, 5, 1, 1, 4, 2, 3], 10]",
            "List[10, 9, 5, 6, 1, 4, 2, 3, 1]",
        ),
        ("HeapPush[List[1], 7]", "List[7, 1]"),
        ("HeapPush[List[], 7]", "List[7]"),
        (
            "HeapPush[List[3, 1], 2, Function[List[a, b], Greater[a, b]]]",
            "List[2, 1, 3]",
        ),
        (
            "HeapPush[List[1, 3, 2], 0, Function[List[a, b], Greater[a, b]]]",
            "List[0, 1, 2, 3]",
        ),
        (
            "HeapPush[List[Power[2, -1], Power[3, -1]], Times[2, Power[3, -1]]]",
            "List[Times[2, Power[3, -1]], Power[3, -1], Power[2, -1]]",
        ),
    ]);
    stays(&[
        "HeapPush[List[5, 3], x]",
        "Heapify[List[x, y]]",
        "Heapify[5]",
        "HeapPush[5, 1]",
        "HeapPush[List[1]]",
        "HeapPop[5]",
        "HeapPop[List[]]",
    ]);
}
