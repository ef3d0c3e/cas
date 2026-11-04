//! Tests of the frequencies, classes and other descriptive statistics
//! (CAS-491, CAS-492): giac's answers, exact.

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
fn frequencies_and_classes() {
    check(&[
        // giac's [[1,0.25],[2,0.5],[3,0.25]], exact.
        (
            "Frequencies[List[3, 1, 2, 2]]",
            "List[List[1, 1/4], List[2, 1/2], List[3, 1/4]]",
        ),
        (
            "CumulatedFrequencies[List[3, 1, 2, 2]]",
            "List[List[1, 1/4], List[2, 3/4], List[3, 1]]",
        ),
        ("Frequencies[List[a, b, a]]", "Frequencies[List[a, b, a]]"),
        // giac's [[0.0 .. 3.0,2],[3.0 .. 6.0,3],[6.0 .. 9.0,2]].
        (
            "Classes[List[1, 2, 3, 4, 5, 6, 7], 0, 3]",
            "List[List[Range[0, 3], 2], List[Range[3, 6], 3], List[Range[6, 9], 2]]",
        ),
        (
            "Classes[List[1/2, 3/2, 5/2, 27/10], 0, 1]",
            "List[List[Range[0, 1], 1], List[Range[1, 2], 1], List[Range[2, 3], 2]]",
        ),
        // An empty class between kept; from the least value by 1.
        (
            "Classes[List[1, 7], 0, 3]",
            "List[List[Range[0, 3], 1], List[Range[3, 6], 0], List[Range[6, 9], 1]]",
        ),
        (
            "Classes[List[3, 1, 2]]",
            "List[List[Range[1, 2], 1], List[Range[2, 3], 1], List[Range[3, 4], 1]]",
        ),
        // A value below the start stays.
        ("Classes[List[1, 7], 2, 3]", "Classes[List[1, 7], 2, 3]"),
    ]);
}

#[test]
fn moments() {
    check(&[
        // giac's [5/3,15/(2*sqrt(57))], sqrt(30)/2 and the window means.
        (
            "CovarianceCorrelation[List[1, 2, 3], List[2, 4, 7]]",
            "List[5/3, Times[15, Power[Times[2, Power[57, 1/2]], -1]]]",
        ),
        (
            "CovarianceCorrelation[List[List[1, 2], List[2, 4], List[3, 7]]]",
            "List[5/3, Times[15, Power[Times[2, Power[57, 1/2]], -1]]]",
        ),
        (
            "RootMeanSquare[List[1, 2, 3, 4]]",
            "Times[1/2, Power[30, 1/2]]",
        ),
        ("RootMeanSquare[List[3, 4]]", "Times[5, Power[2, -1/2]]"),
        (
            "MovingAverage[List[1, 2, 3, 4, 5, 6], 2]",
            "List[3/2, 5/2, 7/2, 9/2, 11/2]",
        ),
        ("MovingAverage[List[1, 2, 3, 4, 5], 3]", "List[2, 3, 4]"),
        (
            "MovingAverage[List[1, 2], 3]",
            "MovingAverage[List[1, 2], 3]",
        ),
    ]);
}
