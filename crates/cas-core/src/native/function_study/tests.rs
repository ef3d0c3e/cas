//! Tests of the function study builtins (CAS-581, CAS-583): giac's
//! answers to `domain` and `convex`, and what stays.

use cas_num::Rational;

use super::at_least_alpha;
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

/// giac's forms: one condition, the whole line, excluded values alone or
/// joined to the conditions, a list of intervals.
#[test]
fn domains() {
    check(&[
        ("FunctionDomain[Log[Plus[x, 1]]]", "Greater[x, -1]"),
        ("FunctionDomain[Power[x, -1]]", "Unequal[x, 0]"),
        ("FunctionDomain[Sin[x]]", "x"),
        ("FunctionDomain[Plus[Power[x, 2], 1]]", "x"),
        (
            "FunctionDomain[Power[Plus[Power[x, 2], -1], 1/2]]",
            "List[LessEqual[x, -1], GreaterEqual[x, 1]]",
        ),
        (
            "FunctionDomain[Times[Log[x], Power[Plus[x, -2], -1]]]",
            "List[And[Greater[x, 0], Unequal[x, 2]]]",
        ),
        ("FunctionDomain[Power[x, 1/2]]", "GreaterEqual[x, 0]"),
        (
            "FunctionDomain[Times[Power[x, 1/2], Power[Plus[x, 1], -1]]]",
            "List[GreaterEqual[x, 0]]",
        ),
        ("FunctionDomain[Log[Plus[t, -1]], t]", "Greater[t, 1]"),
        ("FunctionDomain[Tan[x]]", "Unequal[x, Times[1/2, Pi]]"),
    ]);
}

/// Mode 1: the excluded values, giac's singular values.
#[test]
fn singular_values() {
    check(&[
        ("FunctionDomain[Power[x, 1/2], x, 1]", "List[0]"),
        ("FunctionDomain[Log[x], 1]", "List[0]"),
        ("FunctionDomain[Sin[x], x, 1]", "List[]"),
        ("FunctionDomain[ArcSin[x], x, 1]", "List[-1, 1]"),
        ("FunctionDomain[ArcTanh[x], x, 1]", "List[1, -1]"),
    ]);
}

/// What stays: a mode other than 0 and 1, a variable that is no symbol.
#[test]
fn domains_that_stay() {
    check(&[
        (
            "FunctionDomain[Log[x], x, 2]",
            "FunctionDomain[Log[x], x, 2]",
        ),
        (
            "FunctionDomain[Log[x], Plus[x, 1]]",
            "FunctionDomain[Log[x], Plus[x, 1]]",
        ),
    ]);
}

/// One variable by the sign of `f''`; several by the Hessian's `L D Lᵀ`.
#[test]
fn convexity() {
    check(&[
        ("Convex[Power[x, 2], x]", "True"),
        ("Convex[Times[-1, Power[x, 2]], x]", "False"),
        ("Convex[Power[x, 3], x]", "GreaterEqual[x, 0]"),
        ("Convex[Power[x, 4], List[x]]", "True"),
        ("Convex[Times[a, Power[x, 2]], x]", "GreaterEqual[a, 0]"),
        ("Convex[Times[x, y], List[x, y]]", "False"),
        ("Convex[Plus[Power[x, 2], Power[y, 2]], List[x, y]]", "True"),
        (
            "Convex[Plus[Power[x, 2], Times[-1, Power[y, 2]]], List[x, y]]",
            "False",
        ),
        ("Convex[Times[x, y, z], List[x, y, z]]", "False"),
        (
            "Convex[Plus[Times[a, Power[x, 2]], Times[b, Power[y, 2]]], List[x, y]]",
            "List[GreaterEqual[a, 0], GreaterEqual[b, 0]]",
        ),
        (
            "Convex[Plus[Power[x1, 2], Times[2, Power[x2, 2]], Times[a, Power[x3, 2]], \
             Times[-2, x1, x2], Times[2, x1, x3], Times[-6, x2, x3]], List[x1, x2, x3], Simplify]",
            "List[GreaterEqual[a, 5]]",
        ),
    ]);
}

/// What stays: variables that are no symbols (giac's `u(t)`), an unknown
/// option.
#[test]
fn convexity_that_stays() {
    check(&[
        (
            "Convex[Power[u[t], 2], u[t]]",
            "Convex[Power[u[t], 2], u[t]]",
        ),
        ("Convex[Power[x, 2], x, y]", "Convex[Power[x, 2], x, y]"),
    ]);
}

/// The pivot tests against `α = (1 + √17)/8 ≈ 0.6404`, exactly.
#[test]
fn alpha() {
    let q = |n: i64, d: i64| &Rational::from(n) / &Rational::from(d);
    assert!(at_least_alpha(&q(2, 1), &q(2, 1)));
    assert!(!at_least_alpha(&q(1, 1), &q(2, 1)));
    assert!(at_least_alpha(&q(6405, 10000), &q(1, 1)));
    assert!(!at_least_alpha(&q(6403, 10000), &q(1, 1)));
    assert!(at_least_alpha(&q(0, 1), &q(0, 1)));
}
