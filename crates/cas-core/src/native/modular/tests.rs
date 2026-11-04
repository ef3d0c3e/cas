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

/// Over a prime field (CAS-453): giac's factorizations, gcds and divisions.
#[test]
fn prime_fields() {
    check(&[
        // giac's (x^2+x-1)(x^2-x-1) % 3.
        (
            "Factor[Modular[Plus[Power[x, 4], 1], 3]]",
            "Times[Plus[-1, x, Power[x, 2]], Plus[-1, Times[-1, x], Power[x, 2]]]",
        ),
        // x^8 − x over 𝔽₂: every irreducible of degree 1 and 3.
        (
            "Factor[Modular[Plus[Power[x, 8], Times[-1, x]], 2]]",
            "Times[x, Plus[1, x], Plus[1, x, Power[x, 3]], Plus[1, Power[x, 2], Power[x, 3]]]",
        ),
        // giac's ((1 % 2)*x+1 % 2)^2.
        (
            "Factor[Modular[Plus[Power[x, 2], 1], 2]]",
            "Power[Plus[1, x], 2]",
        ),
        // The leading coefficient apart, −1 for 2 modulo 3.
        (
            "Factor[Modular[Plus[Times[2, Power[x, 2]], 2], 3]]",
            "Times[-1, Plus[1, Power[x, 2]]]",
        ),
        // giac's (x+1) % 5, (-1 % 5)*x+1 % 5 and (1 % 5)*x.
        (
            "GCD[Modular[Plus[Power[x, 2], -1], 5], Modular[Plus[Power[x, 2], Times[3, x], 2], 5]]",
            "Modular[Plus[1, x], 5]",
        ),
        (
            "PolynomialRemainder[Modular[Plus[Power[x, 3], 1], 5], Modular[Plus[Power[x, 2], 1], 5], x]",
            "Modular[Plus[1, Times[-1, x]], 5]",
        ),
        (
            "PolynomialQuotient[Modular[Plus[Power[x, 3], 1], 5], Modular[Plus[Power[x, 2], 1], 5], x]",
            "Modular[x, 5]",
        ),
        // A composite modulus stays.
        (
            "Factor[Modular[Plus[Power[x, 2], 1], 4]]",
            "Factor[Modular[Plus[1, Power[x, 2]], 4]]",
        ),
    ]);
}

/// Matrices over ℤ/nℤ (CAS-457): giac's determinants, inverses and
/// characteristic polynomials.
#[test]
fn matrices() {
    let m = "Modular[List[List[1, 2], List[3, 4]], 7]";
    check(&[
        // giac's -2 % 7 and -2 % 6.
        (&format!("Det[{m}]"), "Modular[-2, 7]"),
        (
            "Det[Modular[List[List[1, 2], List[3, 4]], 6]]",
            "Modular[-2, 6]",
        ),
        // Mixed entries are taken in.
        (
            "Det[List[List[Modular[1, 7], 2], List[3, 4]]]",
            "Modular[-2, 7]",
        ),
        // giac's [[-2 % 7,1 % 7],[-2 % 7,3 % 7]], and the power -1 the same.
        (
            &format!("Inverse[{m}]"),
            "Modular[List[List[-2, 1], List[-2, 3]], 7]",
        ),
        (
            &format!("Power[{m}, -1]"),
            "Modular[List[List[-2, 1], List[-2, 3]], 7]",
        ),
        // giac's [[1 % 6,2 % 6],[3 % 6,-1 % 6]]: the determinant -1 is
        // invertible modulo 6.
        (
            "Inverse[Modular[List[List[1, 2], List[3, 5]], 6]]",
            "Modular[List[List[1, 2], List[3, -1]], 6]",
        ),
        // Singular modulo 7, and modulo 5 though not over ℚ (det −5): the
        // call stays.
        (
            "Inverse[Modular[List[List[1, 2], List[2, 4]], 7]]",
            "Inverse[List[List[Modular[1, 7], Modular[2, 7]], List[Modular[2, 7], Modular[-3, 7]]]]",
        ),
        (
            "Inverse[Modular[List[List[1, 2], List[3, 1]], 5]]",
            "Inverse[List[List[Modular[1, 5], Modular[2, 5]], List[Modular[-2, 5], Modular[1, 5]]]]",
        ),
        // giac's (1 % 7)*x^2+(2 % 7)*x-2 % 7, and the coefficients alone.
        (
            &format!("CharacteristicPolynomial[{m}, x]"),
            "Modular[Plus[Power[x, 2], Times[2, x], -2], 7]",
        ),
        (
            &format!("CharacteristicPolynomial[{m}]"),
            "Modular[List[1, 2, -2], 7]",
        ),
    ]);
}
