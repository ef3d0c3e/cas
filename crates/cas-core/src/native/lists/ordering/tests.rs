//! Tests of Xcas's sorting commands (CAS-622): giac's answers from the
//! phase 1.6 T8 probes, and the calls its errors leave.

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
fn ascending_and_descending() {
    check(&[
        ("SortAscending[List[3, 4, 2]]", "List[2, 3, 4]"),
        ("SortAscending[3, 4, 2]", "List[2, 3, 4]"),
        ("SortAscending[List[2, 1, 2, 1]]", "List[1, 1, 2, 2]"),
        ("SortAscending[List[x, 1, y]]", "List[1, x, y]"),
        ("SortAscending[List[5]]", "List[5]"),
        ("SortAscending[List[]]", "List[]"),
        ("SortDescending[List[3, 4, 2]]", "List[4, 3, 2]"),
        ("SortDescending[3, 4, 2]", "List[4, 3, 2]"),
        ("SortDescending[List[x, 1, y]]", "List[y, x, 1]"),
        (
            "SortDescending[List[\"a\", \"c\", \"b\"]]",
            "List[\"c\", \"b\", \"a\"]",
        ),
        (
            "SortDescending[List[Power[2, -1], 0x1.3333333333333p-2, 1]]",
            "List[1, Power[2, -1], 0x1.3333333333333p-2]",
        ),
        // Ties reversed, as giac's.
        (
            "SortDescending[List[0x1.8p+0, Times[3, Power[2, -1]]]]",
            "List[Times[3, Power[2, -1]], 0x1.8p+0]",
        ),
    ]);
}

#[test]
fn matrices() {
    check(&[
        (
            "SortAscending[List[List[3, 4, 2], List[6, 4, 5]]]",
            "List[List[2, 3, 4], List[5, 6, 4]]",
        ),
        (
            "SortAscending[List[3, 4, 2], List[6, 4, 5]]",
            "List[List[2, 3, 4], List[5, 6, 4]]",
        ),
        (
            "SortAscending[List[List[2, 1, 2], List[3, 2, 1]]]",
            "List[List[1, 2, 2], List[2, 3, 1]]",
        ),
        (
            "SortDescending[List[List[3, 4, 2], List[6, 4, 5]]]",
            "List[List[4, 3, 2], List[4, 6, 5]]",
        ),
        (
            "SortDescending[List[List[2, 1, 2], List[3, 2, 1]]]",
            "List[List[2, 2, 1], List[1, 3, 2]]",
        ),
        (
            "SortAscending[List[List[x, 3, 1], List[1, 2, 3]]]",
            "List[List[1, 3, x], List[3, 2, 1]]",
        ),
        ("SortAscending[List[List[]]]", "List[List[]]"),
        // Rows of unequal lengths: a list of two lists.
        (
            "SortAscending[List[List[3, 1], List[6, 5, 4]]]",
            "List[List[3, 1], List[6, 5, 4]]",
        ),
    ]);
    stays(&["SortAscending[3]", "SortAscending[x]", "SortDescending[]"]);
}

#[test]
fn permutations() {
    check(&[
        (
            "SortPermutation[List[30, 25, 40, 10, 20]]",
            "List[3, 4, 1, 0, 2]",
        ),
        ("SortPermutation[List[3, 1, 2, 1]]", "List[1, 3, 2, 0]"),
        (
            "SortPermutation[List[\"b\", \"a\", \"c\"]]",
            "List[1, 0, 2]",
        ),
        ("SortPermutation[List[x, 1, y]]", "List[1, 0, 2]"),
        ("SortPermutation[List[5]]", "List[0]"),
        ("SortPermutation[List[]]", "List[]"),
        (
            "SortPermutation[List[\"ab\", \"cd\", \"ad\", \"bc\", \"de\"], List[3, 4, 1, 0, 2]]",
            "List[\"bc\", \"de\", \"cd\", \"ab\", \"ad\"]",
        ),
        (
            "SortPermutation[List[3, 1, 2], List[2, 0, 1]]",
            "List[2, 3, 1]",
        ),
        ("SortPermutation[List[5], List[0]]", "List[5]"),
    ]);
    stays(&[
        "SortPermutation[3]",
        "SortPermutation[List[3, 1, 2], 5]",
        "SortPermutation[List[3, 1, 2], List[2, 0]]",
        "SortPermutation[List[3, 1, 2], List[2, 0, 0]]",
        "SortPermutation[List[3, 1, 2], List[3, 0, 1]]",
        "SortPermutation[List[3, 1, 2], List[0x1p+0, 0, 2]]",
        "SortPermutation[List[1, 2], List[0, 1], List[1, 0]]",
    ]);
}
