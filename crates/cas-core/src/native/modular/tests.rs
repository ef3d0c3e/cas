//! Tests of modular values and their arithmetic (CAS-450): giac's residues,
//! inverses and products, expressions reduced, giac's distributed
//! coefficients read back as one value, and what stays.

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
fn residues() {
    check(&[
        // giac's 3 % 7, -2 % 7 for 5, -3 % 7 for 1/2.
        ("Modular[3, 7]", "Modular[3, 7]"),
        ("Modular[5, 7]", "Modular[-2, 7]"),
        ("Modular[1/2, 7]", "Modular[-3, 7]"),
        ("Modular[-7, 7]", "Modular[0, 7]"),
        // n/2 itself is the top of (−n/2, n/2].
        ("Modular[-2, 4]", "Modular[2, 4]"),
        ("Modular[2, 6]", "Modular[2, 6]"),
        // Not invertible: giac's error.
        ("Modular[1/3, 6]", "Indeterminate"),
        // No modulus.
        ("Modular[3, 1]", "Modular[3, 1]"),
        ("Modular[3, n]", "Modular[3, n]"),
    ]);
}

#[test]
fn arithmetic() {
    check(&[
        // giac's (3 % 7)^-1 = -2 % 7 and (3 % 7)*5 = 1 % 7.
        ("Power[Modular[3, 7], -1]", "Modular[-2, 7]"),
        ("Times[Modular[3, 7], 5]", "Modular[1, 7]"),
        ("Plus[Modular[3, 7], Modular[6, 7]]", "Modular[2, 7]"),
        ("Plus[Modular[3, 7], 1/2]", "Modular[0, 7]"),
        // 3^1007 mod 101 = 66, -35 symmetric, by Python's pow.
        ("Power[Modular[3, 101], 1007]", "Modular[-35, 101]"),
        ("Power[Modular[3, 6], -1]", "Indeterminate"),
        // giac's "Mod are different": the sum stays.
        (
            "Plus[Modular[3, 7], Modular[2, 5]]",
            "Plus[Modular[2, 5], Modular[3, 7]]",
        ),
    ]);
}

/// Derivatives and substitutions of modular expressions (CAS-452).
#[test]
fn calculus() {
    check(&[
        // giac's diff((x^3+2*x) % 5,x), 3x^2 + 2.
        (
            "D[Modular[Plus[Power[x, 3], Times[2, x]], 5], x]",
            "Modular[Plus[2, Times[3, Power[x, 2]]], 5]",
        ),
        // The p-th power's derivative vanishes.
        ("D[Modular[Power[x, 5], 5], x]", "Modular[0, 5]"),
        (
            "ReplaceAll[Modular[Plus[Power[x, 2], 1], 7], Rule[x, 3]]",
            "Modular[3, 7]",
        ),
        (
            "ReplaceAll[Modular[Plus[Power[x, 2], 1], 7], Rule[x, Plus[y, 1]]]",
            "Modular[Plus[2, Times[2, y], Power[y, 2]], 7]",
        ),
    ]);
}

#[test]
fn expressions() {
    check(&[
        // giac's ((x+1)^5) % 5 and (x+y)^3 % 3: Frobenius.
        (
            "Modular[Power[Plus[x, 1], 5], 5]",
            "Modular[Plus[1, Power[x, 5]], 5]",
        ),
        (
            "Modular[Power[Plus[x, y], 3], 3]",
            "Modular[Plus[Power[x, 3], Power[y, 3]], 3]",
        ),
        // giac's printed (1 % 3)*x^2+(1 % 3)*x-1 % 3, read back as one value.
        (
            "Plus[Times[Modular[1, 3], Power[x, 2]], Times[Modular[1, 3], x], Modular[-1, 3]]",
            "Modular[Plus[-1, x, Power[x, 2]], 3]",
        ),
        // Coefficients reduced, vanishing terms dropped.
        (
            "Modular[Plus[Times[7, Power[x, 2]], Times[3, x], 10], 5]",
            "Modular[Plus[Times[2, Power[x, 2]], Times[-2, x]], 5]",
        ),
        (
            "Power[Modular[Plus[x, 1], 2], 2]",
            "Modular[Plus[1, Power[x, 2]], 2]",
        ),
    ]);
}
