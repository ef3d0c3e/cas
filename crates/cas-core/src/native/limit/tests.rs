//! Tests of `Limit` (CAS-205): giac 1.9's limits by series, at points and
//! at infinity, from each side and both, with logarithms; what stays, and
//! the warning for a leading coefficient that cannot be told from 0.

use crate::trace::EvalOptions;
use crate::{Context, Expr, Warning};

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

/// giac's `limit` on what series decide.
#[test]
fn limits_by_series() {
    check(&[
        // sin(x)/x, (x^2-1)/(x-1), (1-cos(x))/x^2, (exp(x)-1-x)/x^2.
        ("Limit[Times[Sin[x], Power[x, -1]], x, 0]", "1"),
        ("Limit[Times[Sin[x], Power[x, -1]], Equation[x, 0]]", "1"),
        ("Limit[Times[Sin[x], Power[x, -1]], Rule[x, 0]]", "1"),
        (
            "Limit[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]], x, 1]",
            "2",
        ),
        (
            "Limit[Times[Plus[1, Times[-1, Cos[x]]], Power[x, -2]], x, 0]",
            "1/2",
        ),
        (
            "Limit[Times[Plus[-1, Exp[x], Times[-1, x]], Power[x, -2]], x, 0]",
            "1/2",
        ),
        // 1/x: infinity from both sides, +infinity and -infinity from one.
        ("Limit[Power[x, -1], x, 0]", "ComplexInfinity"),
        ("Limit[Power[x, -1], x, 0, 1]", "Infinity"),
        (
            "Limit[Power[x, -1], Equation[x, 0], -1]",
            "Times[-1, Infinity]",
        ),
        // At infinity: sqrt(x^2+x)-x, (2x^2+1)/(3x^2-x), atan(x) at -inf,
        // (1+1/x)^x, a*x of unknown sign.
        (
            "Limit[Plus[Power[Plus[x, Power[x, 2]], 1/2], Times[-1, x]], x, Infinity]",
            "1/2",
        ),
        (
            "Limit[Times[Plus[1, Times[2, Power[x, 2]]], Power[Plus[Times[-1, x], Times[3, Power[x, 2]]], -1]], x, Infinity]",
            "2/3",
        ),
        (
            "Limit[ArcTan[x], x, Times[-1, Infinity]]",
            "Times[-1/2, Pi]",
        ),
        (
            "Limit[Power[Plus[1, Power[x, -1]], x], x, Infinity]",
            "Exp[1]",
        ),
        ("Limit[Times[a, x], x, Infinity]", "ComplexInfinity"),
        // Logarithms: x*ln(x) at 0+, x^x at 0+, ln(x)/x, x^(1/x) and
        // ln(ln(x))/ln(x) at infinity.
        ("Limit[Times[x, Log[x]], x, 0, 1]", "0"),
        ("Limit[Power[x, x], x, 0, 1]", "1"),
        ("Limit[Times[Log[x], Power[x, -1]], x, Infinity]", "0"),
        ("Limit[Power[x, Power[x, -1]], x, Infinity]", "1"),
        (
            "Limit[Times[Log[Log[x]], Power[Log[x], -1]], x, Infinity]",
            "0",
        ),
        // tan(x) at pi/2: infinity, +infinity from the left.
        ("Limit[Tan[x], x, Times[1/2, Pi]]", "ComplexInfinity"),
        ("Limit[Tan[x], x, Times[1/2, Pi], -1]", "Infinity"),
        // x/abs(x) and sqrt(x) at 0+.
        ("Limit[Times[x, Power[Abs[x], -1]], x, 0, 1]", "1"),
        ("Limit[Power[x, 1/2], x, 0, 1]", "0"),
    ]);
}

/// Gruntz's algorithm where series have an essential singularity
/// (CAS-206): giac's exp(x)/x^10, x*exp(x) at -inf, exp(x)/exp(2x), and
/// the classic exp(x+exp(-x))-exp(x), whose class holds three
/// exponentials, one inside another.
#[test]
fn limits_by_gruntz() {
    check(&[
        (
            "Limit[Times[Exp[x], Power[x, -10]], x, Infinity]",
            "Infinity",
        ),
        (
            "Limit[Times[Power[x, 10], Exp[Times[-1, x]]], x, Infinity]",
            "0",
        ),
        ("Limit[Times[x, Exp[x]], x, Times[-1, Infinity]]", "0"),
        (
            "Limit[Times[Exp[x], Power[Exp[Times[2, x]], -1]], x, Infinity]",
            "0",
        ),
        ("Limit[Exp[Power[x, -1]], x, 0, 1]", "Infinity"),
        ("Limit[Exp[Power[x, -1]], x, 0, -1]", "0"),
        ("Limit[Exp[Times[-1, Power[x, -1]]], x, 0, 1]", "0"),
        (
            "Limit[Times[x, Power[Exp[Power[x, 1/2]], -1]], x, Infinity]",
            "0",
        ),
        (
            "Limit[Plus[Exp[Plus[x, Exp[Times[-1, x]]]], Times[-1, Exp[x]]], x, Infinity]",
            "1",
        ),
    ]);
}

/// What stays: giac's `bounded_function(1)`, sides that differ, a variable
/// that is no symbol.
#[test]
fn limits_that_stay() {
    let ctx = Context::new();
    for input in [
        "Limit[Exp[Power[x, -1]], x, 0]",
        "Limit[Sin[Power[x, -1]], x, 0]",
        "Limit[Times[x, Power[Abs[x], -1]], x, 0]",
        "Limit[f[x], x, 0]",
        "Limit[x, 2, 0]",
        "Limit[x, x, 0, 2]",
    ] {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e.head(), read(input).head(), "{input}: {e:?}");
    }
}

/// A leading coefficient that is 0 but not seen as such: `N` cannot tell
/// it from 0, and the call stays with a warning (CAS-DOC-10 §6).
#[test]
fn an_undecided_leading_coefficient_warns() {
    let ctx = Context::new();
    let input =
        read("Limit[Times[Plus[-1, Power[Cos[1], 2], Power[Sin[1], 2]], Power[x, -1]], x, 0, 1]");
    let (found, report) = ctx.eval_with(&input, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), input.head());
    let messages: Vec<String> = report.warnings.iter().map(Warning::to_string).collect();
    assert!(
        messages
            .iter()
            .any(|m| m.starts_with("the limit's leading coefficient cannot be told from 0")),
        "{messages:?}"
    );
}
