//! Tests of the curves (CAS-410): giac's examples, circles and helices, and
//! what stays.

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

#[test]
fn curvatures() {
    check(&[
        ("Curvature[List[Cos[t], Sin[t]], t]", "1"),
        (
            "Curvature[List[t, Power[t, 2]], t]",
            "Times[2, Power[Plus[1, Times[4, Power[t, 2]]], -3/2]]",
        ),
        // Signed in the plane, as giac's.
        (
            "Curvature[List[t, Times[-1, Power[t, 2]]], t]",
            "Times[-2, Power[Plus[1, Times[4, Power[t, 2]]], -3/2]]",
        ),
        (
            "Curvature[Power[x, 2], x]",
            "Times[2, Power[Plus[1, Times[4, Power[x, 2]]], -3/2]]",
        ),
        (
            "Curvature[List[t, Power[t, 2]], t, 1]",
            "Times[2/25, Power[5, 1/2]]",
        ),
        ("Curvature[List[Cos[t], Sin[t], t], t, 0]", "1/2"),
        // A helix of radius a and pitch b: a/(a² + b²).
        (
            "Curvature[List[Times[2, Cos[t]], Times[2, Sin[t]], Times[3, t]], t, 0]",
            "2/13",
        ),
    ]);
}

#[test]
fn frames() {
    check(&[
        (
            "Frenet[List[t, Power[t, 2]], t, 1]",
            "List[Times[2/25, Power[5, 1/2]], List[-4, 7/2], \
             List[Power[5, -1/2], Times[2, Power[5, -1/2]]], \
             List[Times[-2, Power[5, -1/2]], Power[5, -1/2]]]",
        ),
        (
            "Frenet[List[Cos[t], Sin[t], t], t, 0]",
            "List[1/2, List[-1, 0, 0], 1/2, List[0, Power[2, -1/2], Power[2, -1/2]], \
             List[-1, 0, 0], List[0, Times[-1, Power[2, -1/2]], Power[2, -1/2]]]",
        ),
        // The parabola's evolute, through giac's centre at t = 1.
        (
            "Evolute[List[t, Power[t, 2]], t]",
            "List[Times[-4, Power[t, 3]], Together[Plus[1/2, Times[3, Power[t, 2]]]]]",
        ),
        (
            "Evolute[Power[x, 2], x]",
            "List[Times[-4, Power[x, 3]], Together[Plus[1/2, Times[3, Power[x, 2]]]]]",
        ),
        // A circle's evolute is its centre.
        (
            "Evolute[List[Plus[1, Times[3, Cos[t]]], Times[3, Sin[t]]], t]",
            "List[1, 0]",
        ),
        // What stays: a curve of four coordinates, a `t` that is no
        // symbol, a space curve's evolute.
        (
            "Curvature[List[t, t, t, t], t]",
            "Curvature[List[t, t, t, t], t]",
        ),
        ("Curvature[List[t, t], 1]", "Curvature[List[t, t], 1]"),
        ("Evolute[List[t, t, t], t]", "Evolute[List[t, t, t], t]"),
    ]);
}
