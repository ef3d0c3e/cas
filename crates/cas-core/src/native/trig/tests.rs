//! Tests of the trigonometric rewritings (CAS-233): their forms where giac's
//! are the plain ones, and their values against their arguments' at a
//! point.

use crate::{Context, Expr, Symbol, View};

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

/// Each rewriting `f[e]` has the value of `e` at `x = 2/9`, `y = 7/3`:
/// each part of the difference within 10⁻¹⁰ by `N`, as the exponential
/// forms of real values are complex balls, which `==` does not tell apart
/// from real ones, and the modulus of a complex ball about 0 is unbounded.
fn check_values(f: &str, inputs: &[&str]) {
    let ctx = Context::new();
    let point = [(read("x"), read("2/9")), (read("y"), read("7/3"))];
    let mut wrong = Vec::new();
    for input in inputs {
        let e = ctx.eval(&read(input)).unwrap();
        let rewritten = ctx.eval(&read(&format!("{f}[{input}]"))).unwrap();
        let difference = (e - rewritten.clone()).substitute(&point);
        let difference = ctx.eval(&Expr::call(Symbol::N, [difference])).unwrap();
        let parts = match difference.view() {
            View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
            _ => vec![difference],
        };
        let close = parts.into_iter().all(|p| {
            let gap = Expr::call(Symbol::ABS, [p]);
            let less = Expr::call(Symbol::LESS, [gap, read("1/10000000000")]);
            ctx.eval(&less).unwrap().as_symbol() == Some(Symbol::TRUE)
        });
        if !close {
            wrong.push(format!("{f}[{input}]: {rewritten:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn expansions() {
    check(&[
        ("TrigExpand[Sin[Times[2, x]]]", "Times[2, Cos[x], Sin[x]]"),
        (
            "TrigExpand[Cos[Times[3, x]]]",
            "Plus[Times[-3, Cos[x]], Times[4, Power[Cos[x], 3]]]",
        ),
        (
            "TrigExpand[Sin[Plus[x, y]]]",
            "Plus[Times[Cos[y], Sin[x]], Times[Cos[x], Sin[y]]]",
        ),
        (
            "TrigExpand[Cos[Plus[x, y]]]",
            "Plus[Times[Cos[x], Cos[y]], Times[-1, Sin[x], Sin[y]]]",
        ),
        // Neither a sum nor a multiple: as it is.
        ("TrigExpand[Sin[x]]", "Sin[x]"),
        ("TrigExpand[Sin[Times[1/2, x]]]", "Sin[Times[1/2, x]]"),
        ("TrigExpand[List[Sin[x], 1]]", "List[Sin[x], 1]"),
    ]);
    check_values(
        "TrigExpand",
        &[
            "Sin[Times[5, x]]",
            "Cos[Times[6, x]]",
            "Tan[Times[2, x]]",
            "Tan[Times[3, x]]",
            "Tan[Plus[x, y]]",
            "Cos[Plus[Times[2, x], y]]",
            "Sin[Plus[x, y, 1]]",
        ],
    );
}

#[test]
fn reductions() {
    check(&[
        (
            "TrigReduce[Times[Sin[x], Cos[x]]]",
            "Times[1/2, Sin[Times[2, x]]]",
        ),
        (
            "TrigReduce[Power[Sin[x], 2]]",
            "Plus[1/2, Times[-1/2, Cos[Times[2, x]]]]",
        ),
        (
            "TrigReduce[Power[Cos[x], 3]]",
            "Plus[Times[3/4, Cos[x]], Times[1/4, Cos[Times[3, x]]]]",
        ),
        ("TrigReduce[Plus[Power[Sin[x], 2], Power[Cos[x], 2]]]", "1"),
        // sin(y − x) and sin(x − y) meet.
        (
            "TrigReduce[Plus[Times[Sin[x], Cos[y]], Times[Cos[x], Sin[y]]]]",
            "Sin[Plus[x, y]]",
        ),
        // One factor: as it is.
        ("TrigReduce[Times[a, Sin[x]]]", "Times[a, Sin[x]]"),
    ]);
    check_values(
        "TrigReduce",
        &[
            "Times[Sin[x], Sin[y]]",
            "Times[Power[Sin[x], 2], Cos[x]]",
            "Times[Cos[x], Sin[y], Cos[y]]",
            "Power[Sin[x], 5]",
            "Times[3, Sin[x], Cos[Times[2, x]]]",
        ],
    );
}

#[test]
fn exponentials() {
    check(&[
        (
            "ExpToTrig[Exp[Times[I, x]]]",
            "Plus[Cos[x], Times[I, Sin[x]]]",
        ),
        ("ExpToTrig[TrigToExp[Cos[x]]]", "Cos[x]"),
        // A real exponent: as it is.
        ("ExpToTrig[Exp[x]]", "Exp[x]"),
    ]);
    check_values(
        "TrigToExp",
        &[
            "Sin[x]",
            "Cos[x]",
            "Tan[x]",
            "Plus[Sin[x], Cos[Times[2, y]]]",
        ],
    );
    // Evaluation leaves -(cos x - i sin x) in the round trip of Sin, as
    // giac's exp2trig leaves its forms: the same value.
    check_values(
        "ExpToTrig",
        &["Exp[Plus[x, Times[I, y]]]", "TrigToExp[Sin[x]]"],
    );
}
