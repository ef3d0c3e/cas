//! Tests of division, pseudoprimes and cyclotomic polynomials (CAS-422):
//! giac's examples, Carmichael numbers, and Φ₁₀₅'s coefficients −2.

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
fn division_and_pseudoprimes() {
    check(&[
        ("QuotientRemainder[17, -5]", "List[-3, 2]"),
        ("QuotientRemainder[-17, 5]", "List[-4, 3]"),
        ("SymmetricMod[7, 5]", "2"),
        ("SymmetricMod[8, 5]", "-2"),
        ("SymmetricMod[-7, 5]", "-2"),
        ("SymmetricMod[5, 2]", "1"),
        ("SymmetricMod[9, -6]", "3"),
        ("PseudoprimeQ[13]", "2"),
        ("PseudoprimeQ[15]", "0"),
        // A Carmichael number, and Mersenne's 2^127 − 1 past 2^64.
        ("PseudoprimeQ[561]", "0"),
        ("PseudoprimeQ[Plus[Power[2, 127], -1]]", "1"),
        ("PseudoprimeQ[1]", "0"),
        // What stays.
        ("QuotientRemainder[1, 0]", "QuotientRemainder[1, 0]"),
        ("SymmetricMod[x, 5]", "SymmetricMod[x, 5]"),
    ]);
}

#[test]
fn cyclotomic_polynomials() {
    check(&[
        ("Cyclotomic[12]", "List[1, 0, -1, 0, 1]"),
        ("Cyclotomic[6, x]", "Plus[1, Times[-1, x], Power[x, 2]]"),
        ("Cyclotomic[1, x]", "Plus[-1, x]"),
        ("Cyclotomic[2, x]", "Plus[1, x]"),
        // Φ₁₀₅(2), by the Möbius product in Python's integers.
        (
            "ReplaceAll[Cyclotomic[105, x], Rule[x, 2]]",
            "473474689919911",
        ),
        ("Cyclotomic[0]", "Cyclotomic[0]"),
    ]);
    // Φ₁₀₅, of degree 48, has −2 at x⁷ and x⁴¹, the first coefficient
    // outside {−1, 0, 1}.
    let ctx = Context::new();
    let p = ctx.eval(&read("Cyclotomic[105]")).unwrap();
    let coefficients: Vec<i64> = p.args().iter().rev().map(|c| c.as_i64().unwrap()).collect();
    assert_eq!(coefficients.len(), 49);
    let twos: Vec<usize> = (0..49).filter(|&k| coefficients[k] == -2).collect();
    assert_eq!(twos, [7, 41]);
}
