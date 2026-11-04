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

/// Reduced forms over 𝔽ₚ and modular linear systems (CAS-458): giac's
/// `rref`, `rank`, `ker`, `image` and `linsolve`.
#[test]
fn reduced_forms() {
    check(&[
        (
            "RowReduce[Modular[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 10]], 5]]",
            "Modular[IdentityMatrix[3], 5]",
        ),
        // giac's [[1 % 7,2 % 7],[0 % 7,0 % 7]]; rank 1 modulo 3.
        (
            "RowReduce[Modular[List[List[1, 2], List[2, 4]], 7]]",
            "Modular[List[List[1, 2], List[0, 0]], 7]",
        ),
        (
            "MatrixRank[Modular[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], 3]]",
            "1",
        ),
        // giac's [[2 % 5,-1,0],[-2 % 5,0,-1]] and [[1 % 5,0 % 5],[0 % 5,1 % 5]].
        (
            "NullSpace[Modular[List[List[1, 2, 3], List[2, 4, 6]], 5]]",
            "List[List[Modular[2, 5], -1, 0], List[Modular[-2, 5], 0, -1]]",
        ),
        (
            "ColumnSpace[Modular[List[List[1, 2, 3], List[2, 4, 7]], 5]]",
            "Modular[List[List[1, 0], List[0, 1]], 5]",
        ),
        // giac's [0 % 7,-3 % 7], one side modular or all of them.
        (
            "LinSolve[Modular[List[Equation[Plus[x, Times[2, y]], 1], \
             Equation[Plus[Times[3, x], Times[4, y]], 2]], 7], List[x, y]]",
            "Modular[List[0, -3], 7]",
        ),
        (
            "LinSolve[List[Equation[Plus[x, Times[2, y]], Modular[1, 7]], \
             Equation[Plus[Times[3, x], Times[4, y]], 2]], List[x, y]]",
            "Modular[List[0, -3], 7]",
        ),
        (
            "LinSolve[Modular[List[List[1, 2], List[3, 4]], 7], List[1, 2]]",
            "Modular[List[0, -3], 7]",
        ),
        // giac's [] and [(-1 % 5)*y+1 % 5,y].
        (
            "LinSolve[Modular[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2]], 5], \
             List[x, y]]",
            "List[]",
        ),
        (
            "LinSolve[Modular[List[Equation[Plus[x, y], 1], \
             Equation[Plus[Times[2, x], Times[2, y]], 2]], 5], List[x, y]]",
            "List[Modular[Plus[1, Times[-1, y]], 5], y]",
        ),
        // A composite modulus: the call stays.
        (
            "RowReduce[Modular[List[List[2, 4], List[1, 3]], 6]]",
            "RowReduce[List[List[Modular[2, 6], Modular[-2, 6]], List[Modular[1, 6], Modular[3, 6]]]]",
        ),
    ]);
}

/// Roots modulo a prime and eigenvalues over 𝔽ₚ (CAS-459): giac's `solve`,
/// whose plain residues are modular values here, and `eigenvals`.
#[test]
fn roots_and_eigenvalues() {
    check(&[
        // giac's [-3, 3], from an expression or an equation.
        (
            "Solve[Modular[Plus[Power[x, 2], -2], 7], x]",
            "Modular[List[-3, 3], 7]",
        ),
        (
            "Solve[Equation[Modular[Power[x, 2], 7], Modular[2, 7]], x]",
            "Modular[List[-3, 3], 7]",
        ),
        // giac's [], [-2, 1] (a double root once) and every element of 𝔽₅.
        ("Solve[Modular[Plus[Power[x, 2], 1], 3], x]", "List[]"),
        (
            "Solve[Modular[Times[Power[Plus[x, -1], 2], Plus[x, 2]], 5], x]",
            "Modular[List[-2, 1], 5]",
        ),
        (
            "Solve[Modular[Plus[Power[x, 5], Times[-1, x]], 5], x]",
            "Modular[List[-2, -1, 0, 1, 2], 5]",
        ),
        // giac's 1 % 7,3 % 7; -1 % 3,1 % 3; 1 % 5,1 % 5.
        (
            "Eigenvalues[Modular[List[List[2, 1], List[1, 2]], 7]]",
            "Modular[List[1, 3], 7]",
        ),
        (
            "Eigenvalues[Modular[List[List[0, 1], List[1, 0]], 3]]",
            "Modular[List[-1, 1], 3]",
        ),
        (
            "Eigenvalues[Modular[List[List[1, 1], List[0, 1]], 5]]",
            "Modular[List[1, 1], 5]",
        ),
        // x² + 1 does not split over 𝔽₃, where giac builds GF(3, 2).
        (
            "Eigenvalues[Modular[List[List[0, 1], List[-1, 0]], 3]]",
            "Eigenvalues[List[List[Modular[0, 3], Modular[1, 3]], List[Modular[-1, 3], Modular[0, 3]]]]",
        ),
    ]);
}

/// Chinese remainders of modular values and powers of polynomials modulo a
/// prime and a polynomial (CAS-460): giac's `ichinrem` and `powmod`.
#[test]
fn remainders_and_powers() {
    check(&[
        // giac's -7 % 15, 8 % 105 and -3 % 12.
        (
            "ChineseRemainder[Modular[2, 3], Modular[3, 5]]",
            "Modular[-7, 15]",
        ),
        (
            "ChineseRemainder[List[Modular[2, 3], Modular[3, 5], Modular[1, 7]]]",
            "Modular[8, 105]",
        ),
        (
            "ChineseRemainder[Modular[1, 4], Modular[3, 6]]",
            "Modular[-3, 12]",
        ),
        // giac's "No Integer Solution": the call stays.
        (
            "ChineseRemainder[Modular[2, 4], Modular[1, 6]]",
            "ChineseRemainder[Modular[2, 4], Modular[1, 6]]",
        ),
        // giac's -x, x^2, 3*x+1 and 1.
        (
            "PowerMod[Plus[x, 1], 10, 3, Plus[Power[x, 2], 1], x]",
            "Times[-1, x]",
        ),
        (
            "PowerMod[x, 100, 2, Plus[Power[x, 3], x, 1], x]",
            "Power[x, 2]",
        ),
        (
            "PowerMod[Plus[Times[3, x], 1], 7, 7, Plus[Power[x, 2], 3], x]",
            "Plus[1, Times[3, x]]",
        ),
        (
            "PowerMod[Plus[Power[x, 2], 2], 0, 5, Plus[Power[x, 3], 1]]",
            "1",
        ),
    ]);
}
