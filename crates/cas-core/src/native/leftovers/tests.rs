//! Tests of giac's leftovers (CAS-640 to CAS-643): giac's answers, as its
//! `caseval` prints them, and the calls that stay where giac errs.

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

/// Inputs that stay as they are, giac's errors and `undef`s.
fn stay(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let found = ctx.eval(&read(input)).unwrap();
        if found != read(input) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn roots_and_angles() {
    check(&[
        // giac's 2^(1/3), 8, 2*i, sqrt(x), undef.
        ("NthRoot[3, 2]", "Power[2, 1/3]"),
        ("NthRoot[1/3, 2]", "8"),
        ("NthRoot[2, -4]", "Complex[0, 2]"),
        ("NthRoot[2, x]", "Power[x, 1/2]"),
        ("NthRoot[n, x]", "Power[x, Power[n, -1]]"),
        ("NthRoot[0, 2]", "Indeterminate"),
        // giac's pi/4, pi, -pi/2, 0, pi/4-2*pi/2, pi/3.
        ("ArcTan2[1, 1]", "Times[1/4, Pi]"),
        ("ArcTan2[0, -1]", "Pi"),
        ("ArcTan2[-1, 0]", "Times[-1/2, Pi]"),
        ("ArcTan2[0, 0]", "0"),
        ("ArcTan2[-1, -1]", "Times[-3/4, Pi]"),
        ("ArcTan2[Power[3, 1/2], 1]", "Times[1/3, Pi]"),
        // giac's atan(y/x)+(1-sign(x))*sign(y)*pi/2.
        (
            "ArcTan2[y, x]",
            "Plus[ArcTan[Times[y, Power[x, -1]]], \
             Times[Plus[1, Times[-1, Sign[x]]], Sign[y], Pi, 1/2]]",
        ),
    ]);
}

#[test]
fn has_and_remainders() {
    check(&[
        // giac's 1, 0, 1, 1, 0 (names only), 0, 0.
        ("Has[Plus[x, y], x]", "1"),
        ("Has[Plus[x, y], n]", "0"),
        ("Has[f[x], f]", "1"),
        ("Has[List[x, 1], x]", "1"),
        ("Has[Plus[x, y], Plus[x, y]]", "0"),
        ("Has[Sin[x], Sin]", "0"),
        ("Has[3, x]", "0"),
        // giac's 1, 1/21, x-2*floor(x/2), undef.
        ("FloorMod[10, 3]", "1"),
        ("FloorMod[1/3, 1/7]", "1/21"),
        ("FloorMod[x, 2]", "Plus[x, Times[-2, Floor[Times[1/2, x]]]]"),
        ("FloorMod[7, 0]", "Indeterminate"),
        // giac's 0.5,3.0 and -0.5,-3.0.
        ("FractionalIntegerParts[7/2]", "Sequence[N[1/2], N[3]]"),
        ("FractionalIntegerParts[-7/2]", "Sequence[N[-1/2], N[-3]]"),
    ]);
    stay(&[
        "FractionalIntegerParts[x]",
        "FractionalIntegerParts[Complex[2, 1]]",
    ]);
}

#[test]
fn integers() {
    check(&[
        // giac's 15/128, 1/9, 2*p*(1-p), 3/8.
        (
            "MultinomialProbability[10, List[1/2, 1/2], List[3, 7]]",
            "15/128",
        ),
        (
            "MultinomialProbability[4, List[1/2, 1/3, 1/6], List[1, 2, 1]]",
            "1/9",
        ),
        (
            "MultinomialProbability[2, List[p, Plus[1, Times[-1, p]]], List[1, 1]]",
            "Times[2, p, Plus[1, Times[-1, p]]]",
        ),
        (
            "MultinomialProbability[3, List[1/2, 1/2], List[1, 2]]",
            "3/8",
        ),
        // giac's 6, 3, 1, 2, 333334.
        ("ModularOrder[Modular[3, 7]]", "6"),
        ("ModularOrder[Modular[2, 7]]", "3"),
        ("ModularOrder[Modular[1, 7]]", "1"),
        ("ModularOrder[Modular[6, 7]]", "2"),
        ("ModularOrder[Modular[3, 1000003]]", "333334"),
        // giac's [1,[[2,3],[3,1],[5,1]]], its sign -1, [0,[]], [1,[]].
        (
            "SignedFactorList[120]",
            "List[1, List[List[2, 3], List[3, 1], List[5, 1]]]",
        ),
        (
            "SignedFactorList[-120]",
            "List[-1, List[List[2, 3], List[3, 1], List[5, 1]]]",
        ),
        (
            "SignedFactorList[-12]",
            "List[-1, List[List[2, 2], List[3, 1]]]",
        ),
        ("SignedFactorList[0]", "List[0, List[]]"),
        ("SignedFactorList[1]", "List[1, List[]]"),
        ("SignedFactorList[-1]", "List[-1, List[]]"),
    ]);
    stay(&[
        "MultinomialProbability[3, List[1/2, 1/2], List[1, 1]]",
        "MultinomialProbability[2, List[1/2], List[1, 1]]",
        "ModularOrder[Modular[2, 8]]",
        "ModularOrder[3]",
        "SignedFactorList[x]",
        "SignedFactorList[5/2]",
    ]);
}

#[test]
fn rates_and_boxes() {
    check(&[
        // giac's (f(x+h)-(f(x)))/h and undef; with h = 1/10 giac's
        // decimals, here exact.
        (
            "AverageRateOfChange[f[x], x, h]",
            "Times[Plus[f[Plus[x, h]], Times[-1, f[x]]], Power[h, -1]]",
        ),
        ("AverageRateOfChange[Power[x, 2], x, 0]", "Indeterminate"),
        (
            "AverageRateOfChange[Power[x, 2], x, 1/10]",
            "Times[10, Plus[Power[Plus[x, 1/10], 2], Times[-1, Power[x, 2]]]]",
        ),
        // giac's x=((-2) .. 2),y=((-5) .. 5).
        (
            "BoxConstraints[List[x, y], List[List[-2, 2], List[-5, 5]]]",
            "Sequence[Equation[x, Range[-2, 2]], Equation[y, Range[-5, 5]]]",
        ),
    ]);
    stay(&[
        "BoxConstraints[List[x, y], List[List[0, 1]]]",
        "BoxConstraints[List[x], List[List[0]]]",
    ]);
}
