//! Tests of the matrix constructors and norms (CAS-DOC-30 §4): giac's
//! answers.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

#[test]
fn constructors_and_norms() {
    let ctx = Context::new();
    for (input, expected) in [
        (
            "DiagonalMatrix[List[1, 2, 3]]",
            "List[List[1, 0, 0], List[0, 2, 0], List[0, 0, 3]]",
        ),
        ("DiagonalMatrix[List[List[1, 2], List[3, 4]]]", "List[1, 4]"),
        (
            "HilbertMatrix[3]",
            "List[List[1, 1/2, 1/3], List[1/2, 1/3, 1/4], List[1/3, 1/4, 1/5]]",
        ),
        (
            "VandermondeMatrix[List[1, 2, 3]]",
            "List[List[1, 1, 1], List[1, 2, 4], List[1, 3, 9]]",
        ),
        ("Norm[List[3, 4]]", "5"),
        ("Norm[List[1, I]]", "Power[2, 1/2]"),
        (
            "FrobeniusNorm[List[List[1, 2], List[3, 4]]]",
            "Power[30, 1/2]",
        ),
    ] {
        assert_eq!(
            ctx.eval(&read(input)).unwrap(),
            ctx.eval(&read(expected)).unwrap(),
            "{input}"
        );
    }
    for input in ["HilbertMatrix[0]", "Norm[x]", "DiagonalMatrix[List[]]"] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap().head(), e.head(), "{input}");
    }
}
