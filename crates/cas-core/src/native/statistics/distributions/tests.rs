//! Tests of the probability distributions (CAS-DOC-23): masses,
//! cumulative probabilities and inverses, exact; what stays.

use crate::context::Context;
use crate::expr::Expr;

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

/// Each input evaluates to itself.
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}

/// The binomial distribution: `binomial(n, p, k)`, giac's exact answers
/// (CAS-301).
#[test]
fn binomial() {
    check(&[
        ("Binomial[10, 3/10, 2]", "466948881/2000000000"),
        ("Binomial[10, 3/10, 11]", "0"),
        ("Binomial[10, 3/10, -1]", "0"),
        ("Binomial[10, 3/10, 5/2]", "0"),
        // The coefficient, as before.
        ("Binomial[10, 2]", "45"),
        // Symbols: expressions, as giac's.
        (
            "Binomial[10, p, 2]",
            "Times[45, Power[p, 2], Power[Plus[1, Times[-1, p]], 8]]",
        ),
        (
            "Binomial[10, 1/3, k]",
            "Times[Binomial[10, k], Power[1/3, k], Power[2/3, Plus[10, Times[-1, k]]]]",
        ),
        // p = 1: all of the mass at n, 0⁰ taken as 1.
        ("Binomial[3, 1, 3]", "1"),
        ("Binomial[3, 1, 2]", "0"),
        ("BinomialCDF[10, 3/10, 2]", "478478483/1250000000"),
        ("BinomialCDF[10, 3/10, 2, 4]", "1400846643/2000000000"),
        ("BinomialCDF[10, 3/10, 4, 2]", "0"),
        ("BinomialCDF[10, 1/2, 10]", "1"),
        ("BinomialCDF[10, 1/2, -1]", "0"),
        ("BinomialCDF[10, 1/2, 5/2]", "BinomialCDF[10, 1/2, 2]"),
        ("BinomialCDF[10, 1/2, Infinity]", "1"),
        ("BinomialICDF[10, 3/10, 1/2]", "3"),
        ("BinomialICDF[10, 1/2, 0]", "0"),
        ("BinomialICDF[10, 1/2, 1]", "10"),
        ("BinomialICDF[10, 1/2, 999/1000]", "9"),
        ("BinomialICDF[10, 0, 1]", "0"),
    ]);
}

/// Poisson's distribution: exact, where giac's answers are decimals
/// (CAS-301).
#[test]
fn poisson() {
    check(&[
        ("PoissonPDF[2, 3]", "Times[4/3, Power[E, -2]]"),
        ("PoissonPDF[2, -1]", "0"),
        (
            "PoissonPDF[2, k]",
            "Times[Power[E, -2], Power[2, k], Power[Factorial[k], -1]]",
        ),
        ("PoissonCDF[2, 3]", "Times[19/3, Power[E, -2]]"),
        ("PoissonCDF[1, 1, 2]", "Times[3/2, Power[E, -1]]"),
        ("PoissonCDF[1, -1]", "0"),
        ("PoissonICDF[2, 1/2]", "2"),
        ("PoissonICDF[1/2, 9/10]", "1"),
        ("PoissonICDF[2, 0]", "0"),
        ("PoissonICDF[2, 1]", "Infinity"),
    ]);
}

/// The geometric distribution, from 1, and the negative binomial, the
/// failures before the n-th success (CAS-301).
#[test]
fn geometric_and_negative_binomial() {
    check(&[
        ("GeometricPDF[1/3, 2]", "2/9"),
        ("GeometricPDF[1/3, 0]", "0"),
        (
            "GeometricPDF[1/3, k]",
            "Times[1/3, Power[2/3, Plus[-1, k]]]",
        ),
        ("GeometricCDF[1/3, 2]", "5/9"),
        ("GeometricCDF[1/3, 2, 4]", "38/81"),
        ("GeometricCDF[1/3, 0]", "0"),
        ("GeometricICDF[1/3, 1/2]", "2"),
        // (1/2)^2 = 1 − 3/4 exactly.
        ("GeometricICDF[1/2, 3/4]", "2"),
        ("GeometricICDF[1/3, 0]", "1"),
        ("GeometricICDF[1/3, 1]", "Infinity"),
        ("GeometricICDF[1, 1]", "1"),
        ("NegativeBinomialPDF[3, 1/2, 2]", "3/16"),
        ("NegativeBinomialCDF[3, 1/2, 2]", "1/2"),
        ("NegativeBinomialCDF[3, 1/2, 1, 2]", "3/8"),
        ("NegativeBinomialICDF[3, 1/2, 1/2]", "2"),
        ("NegativeBinomialICDF[3, 1, 1]", "0"),
    ]);
}

/// Parameters that are numbers out of range, symbols where a number is
/// needed, sums past the bound (CAS-301).
#[test]
fn discrete_stays() {
    stays(&[
        "Binomial[10, 3/2, 2]",
        "Binomial[-1, 1/2, 2]",
        "PoissonPDF[-1, 2]",
        "PoissonPDF[0, 2]",
        "GeometricPDF[0, 2]",
        "NegativeBinomialPDF[0, 1/2, 2]",
        "BinomialICDF[10, p, 1/2]",
        "BinomialICDF[10, 1/2, 3/2]",
        "BinomialICDF[10, 1/2, q]",
        "PoissonCDF[2, x]",
        "BinomialCDF[10000, 1/2, 5000]",
        "PoissonICDF[10000, 1/2]",
    ]);
}
