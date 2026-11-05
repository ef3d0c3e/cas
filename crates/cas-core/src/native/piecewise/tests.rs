use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if found == ctx.eval(&read(expected)).unwrap() => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Evaluation (CAS-513), giac's answers.
#[test]
fn evaluation() {
    let sign = "Piecewise[Less[x, 0], -1, Less[x, 1], 0]";
    check(&[
        (sign, sign),
        (&format!("ReplaceAll[{sign}, Rule[x, 1/2]]"), "0"),
        (&format!("ReplaceAll[{sign}, Rule[x, -2]]"), "-1"),
        (&format!("ReplaceAll[{sign}, Rule[x, 3]]"), "Indeterminate"),
        ("Piecewise[True, 1, 2]", "1"),
        ("Piecewise[False, 1, 2]", "2"),
        ("Piecewise[1]", "1"),
        ("Piecewise[False, 1]", "Indeterminate"),
        // An undecided condition keeps its pair; True after it becomes the
        // default, and what follows goes.
        (
            "Piecewise[Less[x, 0], 1, False, 2, True, 3, Less[x, 5], 4]",
            "Piecewise[Less[x, 0], 1, 3]",
        ),
        (
            "Piecewise[Less[x, 0], 1, Less[x, 1], 2, 3]",
            "Piecewise[Less[x, 0], 1, Less[x, 1], 2, 3]",
        ),
        // Values evaluate; numbers go through.
        (
            "ReplaceAll[Piecewise[Less[x, 0], Times[-1, x], x], Rule[x, -1/3]]",
            "1/3",
        ),
        (
            "ReplaceAll[Piecewise[Less[Power[x, 2], 1], 0, 1], Rule[x, 1/2]]",
            "0",
        ),
        ("Piecewise[]", "Indeterminate"),
    ]);
}

/// Conditions decided under the assumptions.
#[test]
fn assumptions() {
    let ctx = Context::new().assuming(read("Greater[a, 0]"));
    let found = ctx.eval(&read("Piecewise[Less[a, 0], 1, 2]")).unwrap();
    assert_eq!(found, read("2"));
}
