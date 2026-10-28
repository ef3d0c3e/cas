//! What the tests of the rewriting forms share: forms compared after
//! evaluation, and values at a point.

use crate::{Context, Expr, Symbol, View};

pub(super) fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
pub(super) fn check(cases: &[(&str, &str)]) {
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
/// each part of the difference within 10⁻¹⁰ by `N`, as a complex ball's
/// modulus about 0 is unbounded.
pub(super) fn check_values(f: &str, inputs: &[&str]) {
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
