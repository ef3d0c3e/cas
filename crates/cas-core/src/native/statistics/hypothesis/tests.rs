//! Tests of the hypothesis tests (CAS-DOC-52): giac's answers where its
//! conventions are this CAS's, and the calls that stay.

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

/// Each input stays as it is.
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e, ctx.eval(&read(input)).unwrap());
        assert_eq!(e, read(input), "{input}: {e:?}");
    }
}

/// `ZTest` and `TTest` (CAS-503).
#[test]
fn location_tests() {
    check(&[
        // giac's normalt([10,30],.5,.02,'!=',0.1): 0.
        ("ZTest[List[10, 30], 1/2, 1/50, \"!=\", 1/10]", "0"),
        // giac's normalt([0.48,50],0.5,0.1,'<'): z = −√2, p = 0.0786.
        ("ZTest[List[12/25, 50], 1/2, 1/10, \"<\"]", "1"),
        ("ZTest[List[12/25, 50], 1/2, 1/10, \"<\", 1/10]", "0"),
        // z = 3 √3 / 2 ≈ 2.6: rejected both ways at 5 %, but not "<".
        ("ZTest[List[3, 2, 5], 1, 1, \">\"]", "0"),
        ("ZTest[List[3, 2, 5], 1, 1]", "0"),
        ("ZTest[List[3, 2, 5], 1, 1, \"<\"]", "1"),
        // giac's studentt([1,2,3,4,5],2,'>'): t = √2, 4 df, p = 0.115.
        ("TTest[List[1, 2, 3, 4, 5], 2, \">\"]", "1"),
        ("TTest[List[1, 2, 3, 4, 5], 0, \">\"]", "0"),
        // giac's studentt([0.48,20],0.5,0.1,'<').
        ("TTest[List[12/25, 20], 1/2, 1/10, \"<\"]", "1"),
        // Decimal data: balls.
        ("ZTest[List[0.48, 50], 0.5, 0.1, \"<\"]", "1"),
    ]);
    stays(&[
        // No σ for a Z test.
        "ZTest[List[1, 2, 3], 1]",
        // An alternative that is not one.
        "ZTest[List[1, 2, 3], 1, 1, \"=\"]",
        // A level outside (0, 1).
        "ZTest[List[1, 2, 3], 1, 1, \"<\", 2]",
        // A size that is not a positive integer; a negative σ.
        "ZTest[List[1, 1/2], 1, 1]",
        "ZTest[List[1, 2, 3], 1, -1]",
        // Symbols.
        "TTest[List[a, b, c], 1]",
    ]);
}
