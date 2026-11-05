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

/// Derivatives and antiderivatives branch by branch (CAS-514), giac's.
#[test]
fn calculus() {
    check(&[
        (
            "D[Piecewise[Less[x, 0], Power[x, 2], Less[y, 0], x, 1], x]",
            "Piecewise[Less[x, 0], Times[2, x], Less[y, 0], 1, 0]",
        ),
        (
            "D[Plus[x, Times[3, Piecewise[Less[x, 0], Sin[x], Power[x, 2]]]], x]",
            "Plus[1, Times[3, Piecewise[Less[x, 0], Cos[x], Times[2, x]]]]",
        ),
        ("D[Piecewise[Less[x, 0], 1, 2], y]", "0"),
        (
            "D[Piecewise[Less[x, 0], Times[-1, x], Power[x, 3]], List[x, 2]]",
            "Piecewise[Less[x, 0], 0, Times[6, x]]",
        ),
        (
            "Integrate[Piecewise[Less[x, 0], Times[-1, x], Power[x, 2]], x]",
            "Piecewise[Less[x, 0], Times[-1/2, Power[x, 2]], Times[1/3, Power[x, 3]]]",
        ),
        (
            "Integrate[Times[2, Piecewise[Less[x, 0], 1, x]], x]",
            "Times[2, Piecewise[Less[x, 0], x, Times[1/2, Power[x, 2]]]]",
        ),
        // A branch that does not integrate leaves the call.
        (
            "Integrate[Piecewise[Less[x, 0], Power[x, x], x], x]",
            "Integrate[Piecewise[Less[x, 0], Power[x, x], x], x]",
        ),
        // Not mapped by its head alone.
        (
            "PiecewiseMap[Function[Times[2, Slot[1]]], x]",
            "Times[2, x]",
        ),
    ]);
}
