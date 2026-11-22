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

/// Elimination ideals (CAS-292), as reduced bases under `revlex` in the
/// other symbols, in canonical order.
#[test]
fn eliminations() {
    check(&[
        // The cusp x = t², y = t³, one symbol or a list.
        (
            "Eliminate[List[Equation[x, Power[t, 2]], Equation[y, Power[t, 3]]], t]",
            "List[Plus[Power[x, 3], Times[-1, Power[y, 2]]]]",
        ),
        (
            "Eliminate[List[Equation[x, Power[t, 2]], Equation[y, Power[t, 3]]], List[t]]",
            "List[Plus[Power[x, 3], Times[-1, Power[y, 2]]]]",
        ),
        // giac's [2*x-4], primitive.
        (
            "Eliminate[List[Equation[Plus[x, y], 1], Equation[Plus[x, Times[-1, y]], 3]], \
             List[y]]",
            "List[Plus[x, -2]]",
        ),
        // Parameters are variables: b x² − a² y under revlex, a > b > x > y.
        (
            "Eliminate[List[Equation[x, Times[a, t]], Equation[y, Times[b, Power[t, 2]]]], t]",
            "List[Plus[Times[b, Power[x, 2]], Times[-1, Power[a, 2], y]]]",
        ),
        // Two relations: y + z and 2z² − 1, giac's [2*y^2-1, y+z].
        (
            "Eliminate[List[Plus[Power[x, 2], Power[y, 2], -1], Plus[x, Times[-1, y]], \
             Plus[x, z]], List[x]]",
            "List[Plus[y, z], Plus[Times[2, Power[z, 2]], -1]]",
        ),
        // The twisted cubic: giac's [y-x^2, z-x^3], reduced.
        (
            "Eliminate[List[Plus[x, Times[-1, t]], Plus[y, Times[-1, Power[t, 2]]], \
             Plus[z, Times[-1, Power[t, 3]]]], t]",
            "List[Plus[Power[y, 2], Times[-1, x, z]], Plus[Times[x, y], Times[-1, z]], \
             Plus[Power[x, 2], Times[-1, y]]]",
        ),
        (
            "Eliminate[List[Plus[Times[x, y], -1], Plus[Power[x, 2], -2]], x]",
            "List[Plus[Times[2, Power[y, 2]], -1]]",
        ),
        // No relation left, inconsistent equations, everything eliminated.
        ("Eliminate[Equation[x, Power[t, 2]], t]", "List[]"),
        (
            "Eliminate[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2]], List[y]]",
            "List[1]",
        ),
        (
            "Eliminate[List[Plus[x, Times[-1, y]], Plus[y, Times[-1, z]]], List[x, y, z]]",
            "List[]",
        ),
        // A variable that does not occur.
        (
            "Eliminate[List[Plus[Power[x, 2], -1]], List[y]]",
            "List[Plus[Power[x, 2], -1]]",
        ),
    ]);
    let ctx = Context::new();
    for input in [
        // Not polynomials.
        "Eliminate[List[Plus[x, Times[-1, Cos[t]]], Plus[y, Times[-1, Sin[t]]]], t]",
        // No variables to eliminate.
        "Eliminate[List[x], List[]]",
        "Eliminate[List[x], 1]",
        "Eliminate[List[x], List[t, t]]",
    ] {
        let answer = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            answer.head(),
            Some(&Expr::from(Symbol::ELIMINATE)),
            "{input}: {answer:?}"
        );
    }
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

/// `IdealMemberQ` (CAS-541): 1 or 0 by the reduced basis, a list for a
/// list; without variables, in the symbols of both.
#[test]
fn memberships() {
    let g = "List[Power[y, 2], Plus[Power[x, 2], Times[2, x, y]]]";
    check(&[
        (
            "IdealMemberQ[Plus[Power[x, 3], -1], List[Plus[x, -1]], List[x]]",
            "1",
        ),
        (
            "IdealMemberQ[Plus[Power[x, 3], -1], List[Plus[x, -2]], List[x]]",
            "0",
        ),
        (
            "IdealMemberQ[Plus[Power[x, 3], -1], List[Plus[x, -1]]]",
            "1",
        ),
        (
            "IdealMemberQ[List[Times[x, y], Plus[x, 1]], List[Power[x, 2], y], List[x, y]]",
            "List[1, 0]",
        ),
        (
            &format!("IdealMemberQ[Power[Plus[x, y], 2], {g}, List[x, y]]"),
            "1",
        ),
        (
            &format!("IdealMemberQ[Plus[x, y], {g}, List[x, y], Lexicographic]"),
            "0",
        ),
        // giac's 0: it reduces by the generators as given.
        (
            "IdealMemberQ[Power[x, 2], List[Plus[x, y], Plus[x, Times[-1, y]]], List[x, y]]",
            "1",
        ),
    ]);
}

/// `GroebnerVariables` (CAS-541): first appearance, stably sorted by
/// `∏ (deg + 1)`; `{}` for anything but a list.
#[test]
fn variables() {
    check(&[
        (
            "GroebnerVariables[List[Plus[Power[x, 2], Times[-1, y]], z]]",
            "List[y, z, x]",
        ),
        (
            "GroebnerVariables[List[Plus[Times[x, y], -1], Plus[Power[x, 3], Times[-1, Power[y, 2]]]]]",
            "List[y, x]",
        ),
        ("GroebnerVariables[List[3, x]]", "List[x]"),
        ("GroebnerVariables[List[]]", "List[]"),
        ("GroebnerVariables[Plus[Power[x, 2], y]]", "List[]"),
    ]);
}

/// `Syzygies` (CAS-541): `{S, V, A, X}`, the `Aᵢ` as the reader makes them.
#[test]
fn syzygies() {
    check(&[
        (
            "Syzygies[List[x, y]]",
            "List[List[Plus[Times[$A1, x], Times[-1, $A0, y]], Power[$A1, 2], Times[$A0, $A1], \
             Power[$A0, 2]], List[$A0, $A1, x, y], List[$A0, $A1], List[x, y]]",
        ),
        (
            "First[Syzygies[List[Plus[x, -1], Plus[Power[x, 2], -1]], List[x]]]",
            "List[Plus[$A0, Times[-1, $A1], Times[$A0, x]], Power[$A1, 2], Times[$A0, $A1], \
             Power[$A0, 2]]",
        ),
    ]);
}

/// `GroebnerBasis`'s `RationalUnivariateRepresentation` (CAS-541), with the
/// order in either place; the basis when there is none.
#[test]
fn representations() {
    let circle = "List[Plus[Power[x, 2], Power[y, 2], -1], Plus[x, Times[-1, y]]]";
    let rur = "RationalUnivariateRepresentation";
    let answer = format!("List[{rur}, y, Plus[1/2, Times[-1, Power[y, 2]]], Times[-2, y], -1, -1]");
    check(&[
        (
            &format!("GroebnerBasis[{circle}, List[x, y], {rur}]"),
            &answer,
        ),
        (
            &format!("GroebnerBasis[{circle}, List[x, y], Lexicographic, {rur}]"),
            &answer,
        ),
        (
            &format!("GroebnerBasis[{circle}, List[x, y], {rur}, Lexicographic]"),
            &answer,
        ),
        (
            &format!(
                "GroebnerBasis[List[Plus[Power[x, 2], -2], Plus[Power[y, 2], -3]], List[x, y], {rur}]"
            ),
            &format!(
                "List[{rur}, Plus[x, y], Plus[-1, Times[10, Power[x, 2]], Times[-1, Power[x, 4]]], \
                 Plus[Times[20, x], Times[-4, Power[x, 3]]], Plus[-8, Times[-8, Power[x, 2]]], \
                 Plus[12, Times[-12, Power[x, 2]]]]"
            ),
        ),
        (
            &format!("GroebnerBasis[List[Plus[Power[x, 3], -2]], List[x], {rur}]"),
            "List[Plus[Power[x, 3], -2]]",
        ),
        (
            &format!("GroebnerBasis[List[Plus[Power[x, 2], Power[y, 2], -1]], List[x, y], {rur}]"),
            "List[Plus[Power[x, 2], Power[y, 2], -1]]",
        ),
    ]);
}

/// What the builtins of CAS-541 cannot take stays.
#[test]
fn new_calls_stay() {
    let ctx = Context::new();
    let rur = "RationalUnivariateRepresentation";
    for (input, head) in [
        // Symbols only.
        (
            "GroebnerVariables[List[Plus[Sin[x], y]]]",
            Symbol::GROEBNER_VARIABLES,
        ),
        // A variable named as an Aᵢ; no generators.
        ("Syzygies[List[$A0, x]]", Symbol::SYZYGIES),
        ("Syzygies[List[]]", Symbol::SYZYGIES),
        // No symbols, a parameter.
        ("IdealMemberQ[1, List[2]]", Symbol::IDEAL_MEMBER_Q),
        (
            "IdealMemberQ[x, List[Plus[x, a]], List[x]]",
            Symbol::IDEAL_MEMBER_Q,
        ),
        // The option twice.
        (
            &format!("GroebnerBasis[List[x, y], List[x, y], {rur}, {rur}]"),
            Symbol::GROEBNER_BASIS,
        ),
    ] {
        let answer = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            answer.head(),
            Some(&Expr::from(head)),
            "{input}: {answer:?}"
        );
    }
}
