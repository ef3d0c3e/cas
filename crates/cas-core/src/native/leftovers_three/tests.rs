//! Tests of giac's list, type, format, date and polynomial leftovers
//! (CAS-644, CAS-645, CAS-647): giac's answers, as its `caseval` prints
//! them, and the calls that stay where giac errs.

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

/// Inputs that stay as they are, giac's errors.
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
fn lists() {
    check(&[
        ("LastElement[1, 2, 3]", "3"),
        ("LastElement[List[1, 2, 3]]", "3"),
        ("LastElement[\"bonjour\"]", "\"r\""),
        ("LastElement[List[]]", "List[]"),
        ("LastElement[x]", "x"),
        ("ArithmeticRange[4]", "List[0, 1, 2, 3]"),
        ("ArithmeticRange[0]", "List[]"),
        ("ArithmeticRange[1, 11, 3]", "List[1, 4, 7, 10]"),
        ("ArithmeticRange[11, 7, -1]", "List[11, 10, 9, 8]"),
        ("ArithmeticRange[1, 1]", "List[]"),
        ("ShiftList[List[0, 1, 2, 3], 2]", "List[2, 3, 0, 0]"),
        ("ShiftList[List[0, 1, 2, 3]]", "List[0, 0, 1, 2]"),
        ("ShiftList[List[0, 1, 2, 3], -2]", "List[0, 0, 0, 1]"),
        ("ShiftList[List[0, 1, 2, 3], 5]", "List[0, 0, 0, 0]"),
        ("ShiftList[List[1, 2, 3], Rational[3, 2]]", "List[0, 1, 2]"),
        ("ShiftList[5, 1]", "2"),
        ("ShiftList[5]", "10"),
        ("ShiftList[-5, 1]", "-3"),
        ("ShiftList[\"abcd\", 1]", "\"bcd \""),
        ("ShiftList[\"abcd\"]", "\"   a\""),
        ("DeleteAt[List[0, 1, 2, 3], 2]", "List[0, 1, 3]"),
        ("DeleteAt[\"0123\", 2]", "\"013\""),
        ("DeleteAt[List[0, 1, 2, 3], 4]", "List[0, 1, 2, 3]"),
        ("DeleteAt[List[0, 1, 2, 3], -1]", "List[0, 1, 2, 3]"),
        ("DeleteAt[List[0, 1, 2, 3], List[0, 2]]", "List[1, 3]"),
        ("DeleteAt[List[0, 1, 2, 3], 1, 2]", "List[0, 3]"),
        ("DeleteAt[List[0, 1, 2, 3], Range[1, 2]]", "List[0, 3]"),
        ("DeleteAt[\"0123\", 5]", "\"\""),
        ("DeleteAt[\"0123\", \"13\"]", "\"02\""),
        ("ElementAt[List[10, 11, 12], 1]", "11"),
        ("ElementAt[List[10, 11, 12], -1]", "12"),
        ("ElementAt[List[List[1, 2], List[3, 4]], List[1, 0]]", "3"),
        ("ElementAt[List[List[1, 2], List[3, 4]], 1]", "List[3, 4]"),
        ("ElementAt[\"abc\", 1]", "\"b\""),
        ("Divide[3, 5]", "Rational[3, 5]"),
        ("Divide[x, y]", "Times[x, Power[y, -1]]"),
    ]);
    stay(&[
        "ArithmeticRange[1, 2, 0]",
        "ArithmeticRange[x]",
        "ArithmeticRange[5, 1]",
        "ElementAt[List[10, 11, 12], 3]",
        "ElementAt[List[List[1, 2], List[3, 4]], List[1, 5]]",
        "DeleteAt[List[0, 1, 2], x]",
    ]);
}

#[test]
fn types_formats_and_dates() {
    check(&[
        ("GetType[\"aze\"]", "\"STR\""),
        ("GetType[Complex[2, 1]]", "\"EXPR\""),
        ("GetType[Function[x, Cos[Times[2, x]]]]", "\"FUNC\""),
        ("GetType[x]", "\"VAR\""),
        ("GetType[Pi]", "\"VAR\""),
        ("GetType[List[1, 2]]", "\"LIST\""),
        ("GetType[List[List[1, 2], List[3, 4]]]", "\"MAT\""),
        ("GetType[List[List[1], List[2, 3]]]", "\"LIST\""),
        ("GetType[Plus[x, 1]]", "\"EXPR\""),
        ("GetType[3]", "\"NUM\""),
        ("GetType[Rational[1, 2]]", "\"NUM\""),
        ("GetType[True]", "\"NUM\""),
        ("GetType[Sin]", "\"OTHER\""),
        ("FormatNumber[3, \"f2\"]", "\"3.00\""),
        ("FormatNumber[Rational[1, 3], \"f4\"]", "\"0.3333\""),
        ("FormatNumber[Rational[-5, 2], \"f0\"]", "\"-2\""),
        ("FormatNumber[Rational[5, 2], \"f0\"]", "\"2\""),
        ("FormatNumber[123456789, \"s3\"]", "\"1.23e+08\""),
        ("FormatNumber[100, \"s3\"]", "\"100\""),
        ("FormatNumber[2, \"e2\"]", "\"2.00e+00\""),
        ("FormatNumber[Rational[1, 100000], \"s3\"]", "\"1e-05\""),
        ("FormatNumber[5, \"f\"]", "\"invalid format\""),
        ("FormatNumber[5, \"x3\"]", "\"5\""),
        ("DayOfWeek[21, 4, 2014]", "1"),
        ("DayOfWeek[15, 10, 1582]", "5"),
        ("DayOfWeek[1, 1, 2000]", "6"),
        ("DayOfWeek[29, 2, 2001]", "4"),
        ("DayOfWeek[1, 1, 1]", "1"),
        ("DayOfWeek[1, 1, -1]", "6"),
    ]);
    stay(&[
        "DayOfWeek[1, 13, 2000]",
        "DayOfWeek[32, 1, 2000]",
        "DayOfWeek[Rational[3, 2], 1, 2000]",
        "DayOfWeek[1, 1]",
    ]);
}

#[test]
fn polynomials() {
    check(&[
        ("RootsWithMultiplicity[3]", "List[]"),
        ("RootsWithMultiplicity[Power[Plus[x, -1], -1]]", "List[]"),
        (
            "RootsWithMultiplicity[Times[Power[Plus[x, -1], 2], Power[Plus[x, 2], -1]]]",
            "List[List[1, 2]]",
        ),
        (
            "TruncateDegree[Power[Plus[Power[x, 2], x], 2], 3]",
            "Plus[Times[2, Power[x, 3]], Power[x, 2]]",
        ),
        ("TruncateDegree[Plus[1, x, Power[x, 2]], 1]", "Plus[x, 1]"),
        (
            "TruncateDegree[Power[Plus[x, y], 3], 2]",
            "Plus[Times[3, Power[x, 2], y], Times[3, x, Power[y, 2]], Power[y, 3]]",
        ),
        ("TruncateDegree[Sin[x], 3]", "Sin[x]"),
        (
            "TruncateDegree[Plus[Times[x, Sin[x]], Power[x, 3]], 2]",
            "Times[x, Sin[x]]",
        ),
        ("TruncateDegree[Plus[Power[x, 7], Power[x, 6], 1]]", "1"),
        ("TruncateDegree[Plus[Power[x, 2], 1], -1]", "0"),
        (
            "TruncateDegree[Power[Plus[x, y], 2], List[y], 1]",
            "Plus[Times[2, x, y], Power[x, 2]]",
        ),
        ("PolynomialDivisors[6]", "List[1, 6]"),
        ("PolynomialDivisors[Power[x, 2]]", "List[1, x, Power[x, 2]]"),
        (
            "PolynomialDivisors[Plus[Power[x, 2], -1]]",
            "List[1, Plus[x, -1], Plus[x, 1], Times[Plus[x, -1], Plus[x, 1]]]",
        ),
        (
            "PolynomialDivisors[Plus[Times[2, Power[x, 2]], -2]]",
            "List[1, 2, Plus[x, -1], Times[2, Plus[x, -1]], Plus[x, 1], \
             Times[2, Plus[x, 1]], Times[Plus[x, -1], Plus[x, 1]], \
             Times[2, Plus[x, -1], Plus[x, 1]]]",
        ),
        (
            "PolynomialDivisors[List[6, Power[x, 2]]]",
            "List[List[1, 6], List[1, x, Power[x, 2]]]",
        ),
        ("PolynomialDivisors[0]", "List[1, 0]"),
    ]);
    stay(&[
        "TruncateDegree[Power[Plus[1, x], -1], 2]",
        "TruncateDegree[Plus[Power[x, 2], 1], Rational[3, 2]]",
    ]);
}
