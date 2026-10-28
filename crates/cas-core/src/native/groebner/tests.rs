//! Tests of `GroebnerBasis` and `GroebnerReduce` (CAS-288): the answers'
//! forms, and what stays.

use crate::context::Context;
use crate::expr::Expr;
use crate::symbol::Symbol;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, answer)` in full form, both evaluated: the same expression.
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

/// Reduced bases, primitive, in increasing order of the leading monomials.
#[test]
fn bases() {
    let circle = "List[Plus[Power[x, 2], Power[y, 2], -1], Plus[x, Times[-1, y]]]";
    let clo = "List[Plus[Power[x, 3], Times[-2, x, y]], \
               Plus[Times[Power[x, 2], y], Times[-2, Power[y, 2]], x]]";
    check(&[
        // revlex by default: x − y below 2y² − 1.
        (
            &format!("GroebnerBasis[{circle}, List[x, y]]"),
            "List[Plus[x, Times[-1, y]], Plus[Times[2, Power[y, 2]], -1]]",
        ),
        (
            &format!("GroebnerBasis[{circle}, List[x, y], Lexicographic]"),
            "List[Plus[Times[2, Power[y, 2]], -1], Plus[x, Times[-1, y]]]",
        ),
        // Equations.
        (
            "GroebnerBasis[List[Equation[Plus[Power[x, 2], Power[y, 2]], 1], \
             Equation[x, y]], List[x, y]]",
            "List[Plus[x, Times[-1, y]], Plus[Times[2, Power[y, 2]], -1]]",
        ),
        (
            &format!("GroebnerBasis[{clo}, List[x, y]]"),
            "List[Plus[Times[2, Power[y, 2]], Times[-1, x]], Times[x, y], Power[x, 2]]",
        ),
        (
            &format!("GroebnerBasis[{clo}, List[x, y], Lexicographic]"),
            "List[Power[y, 3], Plus[x, Times[-2, Power[y, 2]]]]",
        ),
        (
            "GroebnerBasis[List[Plus[Power[x, 2], y], Plus[Times[x, y], 1]], List[x, y], \
             DegreeLexicographic]",
            "List[Plus[Power[y, 2], Times[-1, x]], Plus[Times[x, y], 1], Plus[Power[x, 2], y]]",
        ),
        // The whole ring, the zero ideal, a rational coefficient.
        (
            "GroebnerBasis[List[Plus[Power[x, 2], -1], Plus[x, -2]], List[x]]",
            "List[1]",
        ),
        ("GroebnerBasis[List[0], List[x]]", "List[]"),
        (
            "GroebnerBasis[List[Plus[x, -1/2], y], List[x, y]]",
            "List[y, Plus[Times[2, x], -1]]",
        ),
    ]);
}

/// Exact remainders, giac's.
#[test]
fn remainders() {
    let basis = "List[Plus[Times[2, Power[y, 2]], -1], Plus[x, Times[-1, y]]]";
    check(&[
        (
            &format!("GroebnerReduce[Plus[Power[x, 3], Power[y, 2]], {basis}, List[x, y]]"),
            "Plus[Times[1/2, y], 1/2]",
        ),
        (
            &format!("GroebnerReduce[Times[x, y], {basis}, List[x, y]]"),
            "1/2",
        ),
        (
            "GroebnerReduce[Power[x, 2], List[Plus[x, Times[-1, y]]], List[x, y], Lexicographic]",
            "Power[y, 2]",
        ),
        (
            "GroebnerReduce[Times[Power[x, 2], y], List[Plus[Times[x, y], -1]], List[x, y]]",
            "x",
        ),
    ]);
}

/// What is not a call of polynomials in the variables stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // A parameter.
        "GroebnerBasis[List[Plus[Times[x, y], -1], Plus[Power[x, 2], Times[-1, a]]], List[x, y]]",
        // A function of a variable.
        "GroebnerBasis[List[Plus[Sin[x], y]], List[x, y]]",
        // A negative power.
        "GroebnerBasis[List[Plus[Power[x, -1], y]], List[x, y]]",
        // No order, no variables.
        "GroebnerBasis[List[x], List[x], Foo]",
        "GroebnerBasis[List[x], List[x, x]]",
        "GroebnerBasis[List[x], List[]]",
        "GroebnerBasis[x, List[x]]",
    ] {
        let answer = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            answer.head(),
            Some(&Expr::from(Symbol::GROEBNER_BASIS)),
            "{input}: {answer:?}"
        );
    }
}
