//! Tests of the list operations (CAS-DOC-28): giac's answers to each, and
//! what stays.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same expression.
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

/// giac's `seq` and `makelist` (CAS-334).
#[test]
fn building() {
    check(&[
        ("Table[Power[k, 2], k, 1, 5]", "List[1, 4, 9, 16, 25]"),
        (
            "Table[Power[k, 2], Equation[k, Range[1, 5]]]",
            "List[1, 4, 9, 16, 25]",
        ),
        ("Table[k, k, 1, 10, 3]", "List[1, 4, 7, 10]"),
        ("Table[k, k, 5, 1, -1]", "List[5, 4, 3, 2, 1]"),
        ("Table[k, k, 1, 2, 1/2]", "List[1, 3/2, 2]"),
        ("Table[k, k, 3, 1]", "List[]"),
        (
            "Table[Times[a, k], k, 1, 3]",
            "List[a, Times[2, a], Times[3, a]]",
        ),
        (
            "MakeList[Function[List[x], Power[x, 2]], 1, 9, 2]",
            "List[1, 9, 25, 49, 81]",
        ),
    ]);
}

/// giac's `sort`, `revlist`, `rotate`, `concat`, `append`, `prepend`,
/// `head`, `tail`, `flatten` and `dim` (CAS-335).
#[test]
fn reordering_and_reshaping() {
    check(&[
        ("Sort[List[3, 1, 2]]", "List[1, 2, 3]"),
        (
            "Sort[List[3, 1/2, Power[2, 1/2]]]",
            "List[1/2, Power[2, 1/2], 3]",
        ),
        ("Sort[List[b, a, c]]", "List[a, b, c]"),
        ("Sort[List[x, 2, Pi, 1]]", "List[1, 2, Pi, x]"),
        (
            "Sort[List[3, 1, 2], Function[List[a, b], Greater[a, b]]]",
            "List[3, 2, 1]",
        ),
        (
            "Sort[List[List[2, 1], List[1, 3]]]",
            "List[List[1, 3], List[2, 1]]",
        ),
        ("Reverse[List[1, 2, 3]]", "List[3, 2, 1]"),
        ("Rotate[List[1, 2, 3]]", "List[3, 1, 2]"),
        ("Rotate[List[1, 2, 3, 4], 2]", "List[3, 4, 1, 2]"),
        ("Rotate[List[1, 2, 3, 4], -1]", "List[4, 1, 2, 3]"),
        ("Join[List[1], List[2], List[3]]", "List[1, 2, 3]"),
        ("Join[List[1, 2], 3]", "List[1, 2, 3]"),
        ("Append[List[1, 2], 3]", "List[1, 2, 3]"),
        ("Prepend[List[1, 2], 0]", "List[0, 1, 2]"),
        ("First[List[1, 2, 3]]", "1"),
        ("Rest[List[1, 2, 3]]", "List[2, 3]"),
        (
            "Flatten[List[List[1, 2], List[3, List[4, 5]]]]",
            "List[1, 2, 3, 4, 5]",
        ),
        ("Dimensions[List[1, 2, 3]]", "3"),
        (
            "Dimensions[List[List[1, 2, 3], List[4, 5, 6]]]",
            "List[2, 3]",
        ),
    ]);
}

/// giac's `select`, `remove`, `count`, `member`, `zip`, `cumsum`, and
/// `sum` and `product` of a list (CAS-336).
#[test]
fn filtering_and_combining() {
    check(&[
        (
            "Select[Function[List[x], Greater[x, 1]], List[1, 2, 3]]",
            "List[2, 3]",
        ),
        ("Select[EvenQ, List[1, 2, 3, 4]]", "List[2, 4]"),
        ("Reject[2, List[1, 2, 2, 3]]", "List[1, 3]"),
        (
            "Reject[Function[List[x], Greater[x, 2]], List[1, 2, 3, 4]]",
            "List[1, 2]",
        ),
        (
            "Count[Function[List[x], Greater[x, 1]], List[1, 2, 3]]",
            "2",
        ),
        ("Count[2, List[1, 2, 2, 3]]", "2"),
        ("MemberIndex[2, List[1, 2, 3]]", "2"),
        ("MemberIndex[5, List[1, 2, 3]]", "0"),
        (
            "Zip[Function[List[x, y], Plus[x, y]], List[1, 2, 3], List[4, 5, 6]]",
            "List[5, 7, 9]",
        ),
        (
            "Accumulate[List[a, b, c]]",
            "List[a, Plus[a, b], Plus[a, b, c]]",
        ),
        ("Sum[List[1, 2, 3, a]]", "Plus[6, a]"),
        ("Product[List[1, 2, 3, a]]", "Times[6, a]"),
    ]);
}

/// What stays: symbolic bounds, a step 0, an empty list's first or rest,
/// an undecided predicate or comparison, what is no list.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Table[k, k, 1, n]",
        "Table[k, k, 1, 2, 0]",
        "First[List[]]",
        "Rest[List[]]",
        "Select[Function[List[x], Greater[x, a]], List[1, 2]]",
        "Sort[List[3, 1], Function[List[a, b], f[a, b]]]",
        "Reverse[x]",
        "Dimensions[x]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
