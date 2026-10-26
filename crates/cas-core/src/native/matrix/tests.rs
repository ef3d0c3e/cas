//! Tests of the matrix basics (CAS-218): giac 1.9's `tran`, `idn`,
//! `trace`, `dot` and `cross`, and what stays.

use crate::trace::EvalOptions;
use crate::{Context, EvalError, Expr, Warning};

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

const M23: &str = "List[List[1, 2, 3], List[4, 5, 6]]";

#[test]
fn transposes() {
    check(&[
        (
            &format!("Transpose[{M23}]"),
            "List[List[1, 4], List[2, 5], List[3, 6]]",
        ),
        (
            "Transpose[List[1, 2, 3]]",
            "List[List[1], List[2], List[3]]",
        ),
        (
            "Transpose[List[List[a, b], List[c, d]]]",
            "List[List[a, c], List[b, d]]",
        ),
        ("Transpose[x]", "Transpose[x]"),
    ]);
}

#[test]
fn identities_and_traces() {
    check(&[
        (
            "IdentityMatrix[3]",
            "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]",
        ),
        ("IdentityMatrix[1]", "List[List[1]]"),
        ("IdentityMatrix[0]", "IdentityMatrix[0]"),
        ("IdentityMatrix[n]", "IdentityMatrix[n]"),
        ("Tr[List[List[1, 2], List[3, 4]]]", "5"),
        ("Tr[List[List[a, b], List[c, d]]]", "Plus[a, d]"),
        (&format!("Tr[{M23}]"), &format!("Tr[{M23}]")),
    ]);
}

#[test]
fn dot_products() {
    check(&[
        ("Dot[List[1, 2, 3], List[4, 5, 6]]", "32"),
        // The conjugate of the first, as giac's dot([i,0],[1,0]) is -i.
        ("Dot[List[I, 0], List[1, 0]]", "Times[-1, I]"),
        ("Dot[List[1, 0], List[I, 0]]", "I"),
        (
            "Dot[List[Plus[1, I], 2], List[3, Times[4, I]]]",
            "Plus[3, Times[5, I]]",
        ),
        // Real symbols are their own conjugates; of matrices, entrywise.
        (
            "Dot[List[a, b], List[c, d]]",
            "Plus[Times[a, c], Times[b, d]]",
        ),
        ("Dot[List[List[1, 2], List[3, 4]], IdentityMatrix[2]]", "5"),
    ]);
}

#[test]
fn cross_products() {
    check(&[
        ("Cross[List[1, 0, 0], List[0, 1, 0]]", "List[0, 0, 1]"),
        ("Cross[List[1, 2], List[3, 4]]", "List[0, 0, -2]"),
        (
            "Cross[List[a, b, c], List[x, y, z]]",
            "List[Plus[Times[b, z], Times[-1, c, y]], Plus[Times[c, x], Times[-1, a, z]], \
             Plus[Times[a, y], Times[-1, b, x]]]",
        ),
    ]);
}

/// Dimensions that do not fit are shape errors; a trace needs a square.
#[test]
fn shapes() {
    let ctx = Context::new();
    for input in [
        "Dot[List[1, 2], List[1, 2, 3]]",
        "Cross[List[1, 2, 3, 4], List[1, 2, 3, 4]]",
    ] {
        assert_eq!(ctx.eval(&read(input)), Err(EvalError::Shape), "{input}");
    }
    let (_, report) = ctx.eval_with(&read(&format!("Tr[{M23}]")), &EvalOptions::default());
    let messages: Vec<String> = report.warnings.iter().map(Warning::to_string).collect();
    assert!(
        messages.iter().any(|m| m.contains("a square matrix")),
        "{messages:?}"
    );
}
