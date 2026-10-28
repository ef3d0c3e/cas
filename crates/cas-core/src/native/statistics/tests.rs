//! Tests of the statistics builtins (CAS-DOC-22): the data forms, the
//! answers, and what stays.

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

const ONE_TO_FOUR: &str = "List[1, 2, 3, 4]";
const WEIGHTS: &str = "List[1, 1, 1, 5]";
const MATRIX: &str = "List[List[1, 2], List[3, 4], List[5, 7]]";

/// Means of lists, with frequencies, of columns; symbols, radicals.
#[test]
fn means() {
    check(&[
        (&format!("Mean[{ONE_TO_FOUR}]"), "5/2"),
        (&format!("Mean[{ONE_TO_FOUR}, {WEIGHTS}]"), "13/4"),
        (&format!("Mean[{MATRIX}]"), "List[3, 13/3]"),
        // Row frequencies, which giac refuses.
        (
            "Mean[List[List[1, 2], List[3, 4]], List[1, 3]]",
            "List[5/2, 7/2]",
        ),
        ("Mean[List[a, b, c]]", "Times[1/3, Plus[a, b, c]]"),
        ("Mean[List[1/2, 1/3]]", "5/12"),
        (
            "Mean[List[Power[2, 1/2], 1]]",
            "Times[1/2, Plus[1, Power[2, 1/2]]]",
        ),
        // An expected value.
        (
            "Mean[List[x, y], List[p, Plus[1, Times[-1, p]]]]",
            "Plus[Times[p, x], Times[Plus[1, Times[-1, p]], y]]",
        ),
    ]);
}

/// Population variances and deviations, and the sample deviation (giac's
/// `stddevp`, with `Σ f − 1`).
#[test]
fn deviations() {
    check(&[
        (&format!("PopulationVariance[{ONE_TO_FOUR}]"), "5/4"),
        (
            &format!("PopulationVariance[{ONE_TO_FOUR}, {WEIGHTS}]"),
            "19/16",
        ),
        (&format!("PopulationVariance[{MATRIX}]"), "List[8/3, 38/9]"),
        ("PopulationVariance[List[1/2, 1/3]]", "1/144"),
        (
            "PopulationVariance[List[1, 2, 3], List[1/2, 1/3, 1/6]]",
            "5/9",
        ),
        // Expanded: a polynomial in the data.
        (
            "PopulationVariance[List[a, b]]",
            "Plus[Times[1/4, Power[a, 2]], Times[-1/2, a, b], Times[1/4, Power[b, 2]]]",
        ),
        (
            "PopulationVariance[List[Power[2, 1/2], 1]]",
            "Plus[3/4, Times[-1/2, Power[2, 1/2]]]",
        ),
        (
            &format!("PopulationStandardDeviation[{ONE_TO_FOUR}]"),
            "Times[1/2, Power[5, 1/2]]",
        ),
        (
            &format!("PopulationStandardDeviation[{MATRIX}]"),
            "List[Power[8/3, 1/2], Power[38/9, 1/2]]",
        ),
        ("PopulationStandardDeviation[List[5]]", "0"),
        (
            &format!("StandardDeviation[{ONE_TO_FOUR}]"),
            "Power[5/3, 1/2]",
        ),
        (
            &format!("StandardDeviation[{ONE_TO_FOUR}, {WEIGHTS}]"),
            "Power[19/14, 1/2]",
        ),
        (
            "StandardDeviation[List[1, 2, 3], List[2, 2, 2]]",
            "Power[4/5, 1/2]",
        ),
        (
            &format!("StandardDeviation[{MATRIX}]"),
            "List[2, Power[19/3, 1/2]]",
        ),
    ]);
}

/// Not data, or frequencies that leave nothing to divide by.
#[test]
fn what_stays() {
    stays(&[
        "Mean[List[]]",
        "Mean[x]",
        "Mean[List[1, 2], List[1]]",
        "Mean[List[1, 2], x]",
        "Mean[List[List[1, 2], List[3]]]",
        "Mean[List[1, List[2]]]",
        "Mean[List[List[]]]",
        "Mean[List[1, 2], List[1, -1]]",
        "Mean[List[1, 2], List[1, 2], List[3, 4]]",
        "PopulationVariance[List[1, 2], List[0, 0]]",
        "StandardDeviation[List[5]]",
        "StandardDeviation[List[1, 2], List[1/2, 1/2]]",
    ]);
}

/// Medians: the middle value, or the half-sum of the two middle values
/// where giac answers the lower one; with frequencies, by the cumulative
/// frequency (CAS-295).
#[test]
fn medians() {
    check(&[
        ("Median[List[3, 1, 2]]", "2"),
        ("Median[List[4, 1, 3, 2]]", "5/2"),
        ("Median[List[1, 2, 3, 5, 10, 4]]", "7/2"),
        ("Median[List[2, 1, 2, 1]]", "3/2"),
        (&format!("Median[{ONE_TO_FOUR}, {WEIGHTS}]"), "4"),
        // Half the total reached exactly at 2; a zero frequency left out.
        ("Median[List[1, 2, 3], List[0, 1, 1]]", "5/2"),
        ("Median[List[1, 2, 3], List[1, 0, 0]]", "1"),
        ("Median[List[1, 2], List[1/3, 2/3]]", "2"),
        ("Median[List[Power[2, 1/2], 1, Pi]]", "Power[2, 1/2]"),
        (&format!("Median[{MATRIX}]"), "List[3, 4]"),
    ]);
}

/// Quantiles, `x₍⌈np⌉₎`, the minimum and maximum at 0 and 1 where giac
/// answers `undef`; quartiles (CAS-295).
#[test]
fn quantiles() {
    let ten = "List[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]";
    check(&[
        (&format!("Quantile[{ten}, 1/4]"), "3"),
        (&format!("Quantile[{ten}, 0]"), "1"),
        (&format!("Quantile[{ten}, 1]"), "10"),
        (&format!("Quantile[{ten}, 95/100]"), "10"),
        (&format!("Quantile[{ten}, 1/10]"), "1"),
        (&format!("Quantile[{ten}, 11/100]"), "2"),
        (&format!("Quantile[{ten}, Decimal[25, -2]]"), "3"),
        ("Quantile[List[3, 1, 2], List[1, 2, 1], 1/4]", "1"),
        ("Quantile[List[1, 2, 3, 4, 5, 6, 7, 8], 1/3]", "3"),
        ("Quantile[List[Power[2, 1/2], 1, Pi], 1/3]", "1"),
        (&format!("Quartiles[{ten}]"), "List[1, 3, 11/2, 8, 10]"),
        (
            "Quartiles[List[1, 2, 3, 4, 5, 6, 7, 8]]",
            "List[1, 2, 9/2, 6, 8]",
        ),
        ("Quartiles[List[3, 1, 2]]", "List[1, 1, 2, 3, 3]"),
        (
            &format!("Quartiles[{MATRIX}]"),
            "List[List[1, 1, 3, 5, 5], List[2, 2, 4, 7, 7]]",
        ),
        (&format!("FirstQuartile[{ten}]"), "3"),
        (&format!("ThirdQuartile[{ten}]"), "8"),
        ("ThirdQuartile[List[1, 5, 2, 8, 4]]", "5"),
        (
            "FirstQuartile[List[1, 2, 3, 4, 5, 6, 7, 8], List[1, 1, 1, 1, 1, 1, 1, 1]]",
            "2",
        ),
    ]);
}

/// Values that comparisons do not order, frequencies that are not
/// non-negative numbers, fractions outside `[0, 1]` (CAS-295).
#[test]
fn unordered_stays() {
    stays(&[
        "Median[List[x, 1, 2]]",
        "Median[List[1, 2], List[1, -1]]",
        "Median[List[1, 2], List[0, 0]]",
        "Median[List[1, 2], List[a, 1]]",
        "Quantile[List[1, 2], 3/2]",
        "Quantile[List[1, 2], -1/2]",
        "Quantile[List[1, 2], p]",
        "Quantile[List[1, 2]]",
        "Quartiles[List[x, y]]",
        "FirstQuartile[List[]]",
    ]);
}

const XS: &str = "List[1, 2, 3, 4]";
const YS: &str = "List[2, 4, 5, 9]";
const PAIRS: &str = "List[List[1, 2], List[2, 4], List[3, 5], List[4, 9]]";
const PAIR_WEIGHTS: &str = "List[1, 1, 2, 1]";

/// Population covariances and correlations, of two lists or a matrix of
/// pairs, with frequencies or not (CAS-296).
#[test]
fn two_variables() {
    check(&[
        (&format!("PopulationCovariance[{XS}, {YS}]"), "11/4"),
        (
            &format!("PopulationCovariance[{XS}, {YS}, {PAIR_WEIGHTS}]"),
            "11/5",
        ),
        (&format!("PopulationCovariance[{PAIRS}]"), "11/4"),
        (
            &format!("PopulationCovariance[{PAIRS}, {PAIR_WEIGHTS}]"),
            "11/5",
        ),
        (&format!("PopulationCovariance[{XS}, {XS}]"), "5/4"),
        (
            "PopulationCovariance[List[1, 2, 3], List[a, b, c]]",
            "Plus[Times[-1/3, a], Times[1/3, c]]",
        ),
        (
            &format!("Correlation[{XS}, {YS}]"),
            "Times[11, Power[130, -1/2]]",
        ),
        (
            &format!("Correlation[{XS}, {YS}, {PAIR_WEIGHTS}]"),
            "Times[11/26, Power[5, 1/2]]",
        ),
        (
            &format!("Correlation[{PAIRS}]"),
            "Times[11, Power[130, -1/2]]",
        ),
        ("Correlation[List[1, 2, 3], List[2, 4, 6]]", "1"),
        ("Correlation[List[1, 2, 3], List[3, 2, 1]]", "-1"),
        // Constant data: 0/0, giac's `undef`.
        ("Correlation[List[1, 1, 1], List[2, 4, 6]]", "Indeterminate"),
    ]);
    stays(&[
        "PopulationCovariance[List[1, 2], List[1]]",
        "PopulationCovariance[List[List[1, 2, 3], List[4, 5, 6]]]",
        "PopulationCovariance[List[1, 2]]",
        "PopulationCovariance[List[1, 2], List[3, 4], List[1, -1]]",
        "Correlation[x, y]",
    ]);
}
