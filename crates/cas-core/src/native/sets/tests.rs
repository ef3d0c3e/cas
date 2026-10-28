//! Tests of the set operations (CAS-365): giac's `union`, `intersect`,
//! `minus` and `is_included` on sets and lists, membership, and what stays.

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

/// giac's `set[1,2] union set[2,3]`, `[1,2] union [2,3]`,
/// `set[1,2] intersect set[2,3]`, `[1,2,3] minus [2]` and
/// `set[x,y] union set[y,z]`.
#[test]
fn operations() {
    check(&[
        (
            "Union[FiniteSet[1, 2], FiniteSet[2, 3]]",
            "FiniteSet[1, 2, 3]",
        ),
        ("Union[List[1, 2], List[2, 3]]", "FiniteSet[1, 2, 3]"),
        ("Union[FiniteSet[1, 2], List[5]]", "FiniteSet[1, 2, 5]"),
        (
            "Union[FiniteSet[1], FiniteSet[2], FiniteSet[3]]",
            "FiniteSet[1, 2, 3]",
        ),
        ("Union[FiniteSet[], FiniteSet[1]]", "FiniteSet[1]"),
        (
            "Intersection[FiniteSet[1, 2], FiniteSet[2, 3]]",
            "FiniteSet[2]",
        ),
        ("Intersection[List[1, 2], List[2, 3]]", "FiniteSet[2]"),
        (
            "Intersection[FiniteSet[1, 2, 3], FiniteSet[2, 3], FiniteSet[3]]",
            "FiniteSet[3]",
        ),
        (
            "Complement[FiniteSet[1, 2, 3], FiniteSet[2]]",
            "FiniteSet[1, 3]",
        ),
        ("Complement[List[1, 2, 3], List[2]]", "FiniteSet[1, 3]"),
        (
            "Complement[FiniteSet[1, 2, 3], FiniteSet[1], FiniteSet[3]]",
            "FiniteSet[2]",
        ),
        (
            "Union[FiniteSet[x, y], FiniteSet[y, z]]",
            "FiniteSet[x, y, z]",
        ),
        // A symbol among the operands leaves the call.
        ("Union[a, FiniteSet[1]]", "Union[a, FiniteSet[1]]"),
        ("Complement[FiniteSet[1], b]", "Complement[FiniteSet[1], b]"),
    ]);
}

/// giac's `is_included([1],[1,2])`, and membership with `in`.
#[test]
fn inclusion_and_membership() {
    check(&[
        ("IncludedQ[List[1], List[1, 2]]", "True"),
        ("IncludedQ[FiniteSet[1, 3], FiniteSet[1, 2]]", "False"),
        ("IncludedQ[FiniteSet[], FiniteSet[1]]", "True"),
        ("IncludedQ[FiniteSet[x], FiniteSet[x, y]]", "True"),
        (
            "IncludedQ[FiniteSet[x], FiniteSet[y]]",
            "IncludedQ[FiniteSet[x], FiniteSet[y]]",
        ),
        ("Element[2, FiniteSet[1, 2]]", "True"),
        ("Element[3, FiniteSet[1, 2]]", "False"),
        ("Element[3, List[1, 2]]", "False"),
        ("Element[x, FiniteSet[x, 1]]", "True"),
        ("Element[x, FiniteSet[y]]", "Element[x, FiniteSet[y]]"),
        // The number sets still.
        ("Element[2, Integers]", "True"),
    ]);
}
