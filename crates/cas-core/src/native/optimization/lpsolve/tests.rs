//! Tests of the linear programs (CAS-417, CAS-418, CAS-600, CAS-601):
//! giac's examples and options, unbounded and infeasible programs,
//! integers, binaries and decimals.

use crate::trace::EvalOptions;
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
fn programs() {
    check(&[
        (
            "LPSolve[Plus[Times[2, x], y, Times[-1, z], 4], List[LessEqual[x, 1], \
             GreaterEqual[y, 2], Equation[Plus[x, Times[3, y], Times[-1, z]], 2], \
             LessEqual[Plus[Times[2, x], Times[-1, y], z], 8], \
             LessEqual[Plus[Times[-1, x], y], 5]]]",
            "List[-4, List[Equation[x, 0], Equation[y, 5], Equation[z, 13]]]",
        ),
        (
            "LPSolve[Plus[Times[-7, x], Times[2, y]], List[LessEqual[Plus[Times[4, x], \
             Times[-12, y]], 20], LessEqual[Plus[Times[-1, x], Times[3, y]], 3]], \
             Equation[x, Range[-5, 5]], Equation[y, Range[0, Infinity]], Maximize]",
            "List[21, List[Equation[x, -3], Equation[y, 0]]]",
        ),
        (
            "LPSolve[Plus[Times[-1, x], Times[-1, y]], List[LessEqual[y, Plus[Times[3, x], 1/2]], \
             LessEqual[y, Plus[Times[-5, x], 2]]], Equation[Assume, LPNonNegative]]",
            "List[-5/4, List[Equation[x, 3/16], Equation[y, 17/16]]]",
        ),
        (
            "LPSolve[Plus[Times[-6, x], Times[4, y], z], List[LessEqual[Plus[Times[5, x], \
             Times[-10, y]], 20], Equation[Plus[Times[2, z], Times[-3, y]], 6], \
             LessEqual[Plus[Times[-1, x], Times[3, y]], 3]], Equation[x, Range[1, 20]], \
             Equation[y, Range[0, Infinity]]]",
            "List[-133/2, List[Equation[x, 18], Equation[y, 7], Equation[z, 27/2]]]",
        ),
        // Unbounded and infeasible, as giac's.
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[Plus[x, y], 2]], \
             Equation[Assume, LPNonNegative], Maximize]",
            "List[Infinity, List[]]",
        ),
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[Plus[x, y], 2], LessEqual[Plus[x, y], 1]]]",
            "List[]",
        ),
    ]);
}

#[test]
fn integers() {
    check(&[
        (
            "LPSolve[Plus[Times[-5, x], Times[-7, y]], List[LessEqual[Plus[Times[7, x], y], 35], \
             LessEqual[Plus[Times[-1, x], Times[3, y]], 6]], Equation[Assume, integer]]",
            "List[-41, List[Equation[x, 4], Equation[y, 3]]]",
        ),
        (
            "LPSolve[Plus[Times[8, x1], Times[11, x2], Times[6, x3], Times[4, x4]], \
             List[LessEqual[Plus[Times[5, x1], Times[7, x2], Times[4, x3], Times[3, x4]], 14]], \
             Equation[Assume, LPBinary], Maximize]",
            "List[21, List[Equation[x1, 0], Equation[x2, 1], Equation[x3, 1], Equation[x4, 1]]]",
        ),
        (
            "LPSolve[Plus[x, y], List[GreaterEqual[x, 1/3], GreaterEqual[y, 1/2]], \
             Equation[Assume, integer]]",
            "List[2, List[Equation[x, 1], Equation[y, 1]]]",
        ),
        // What stays: a nonlinear objective.
        (
            "LPSolve[Times[x, y], List[LessEqual[x, 1]]]",
            "LPSolve[Times[x, y], List[LessEqual[x, 1]]]",
        ),
    ]);
}

/// `-7x + 2y` under `[4x - 12y ≤ 20, -x + 3y ≤ 3]` with the options `o`.
fn bounded(o: &str) -> String {
    format!(
        "LPSolve[Plus[Times[-7, x], Times[2, y]], List[LessEqual[Plus[Times[4, x], \
         Times[-12, y]], 20], LessEqual[Plus[Times[-1, x], Times[3, y]], 3]], {o}]"
    )
}

/// `-5x - 7y` under `[7x + y ≤ 35, -x + 3y ≤ 6]` with the options `o`.
fn lattice(o: &str) -> String {
    format!(
        "LPSolve[Plus[Times[-5, x], Times[-7, y]], List[LessEqual[Plus[Times[7, x], y], 35], \
         LessEqual[Plus[Times[-1, x], Times[3, y]], 6]], {o}]"
    )
}

/// The knapsack `8x1 + 11x2 + 6x3 + 4x4` under `5x1 + 7x2 + 4x3 + 3x4 ≤ 14`.
fn knapsack(o: &str) -> String {
    format!(
        "LPSolve[Plus[Times[8, x1], Times[11, x2], Times[6, x3], Times[4, x4]], \
         List[LessEqual[Plus[Times[5, x1], Times[7, x2], Times[4, x3], Times[3, x4]], 14]], {o}]"
    )
}

/// `(input, expected)` pairs of owned strings, as [`check`]'s.
fn check_owned(cases: &[(String, &str)]) {
    let cases: Vec<(&str, &str)> = cases.iter().map(|(i, e)| (i.as_str(), *e)).collect();
    check(&cases);
}

/// The calls that stay, unevaluated, with no warning.
fn stays(inputs: &[String]) {
    let ctx = Context::new();
    for input in inputs {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// giac's forms of maximize, the ranges and the assumptions (CAS-600),
/// giac's answers (probe1, probe6).
#[test]
fn options() {
    let max = "List[21, List[Equation[x, -3], Equation[y, 0]]]";
    let min = "List[-35, List[Equation[x, 5], Equation[y, 0]]]";
    let box_ = "Equation[x, Range[-5, 5]], Equation[y, Range[0, Infinity]]";
    let int = "List[-41, List[Equation[x, 4], Equation[y, 3]]]";
    check_owned(&[
        (bounded(&format!("{box_}, LPMaximize")), max),
        (bounded(&format!("{box_}, Equation[LPMaximize, True]")), max),
        (bounded(&format!("{box_}, Equation[Maximize, True]")), max),
        (
            bounded(&format!("{box_}, Equation[LPMaximize, False]")),
            min,
        ),
        (bounded(&format!("{box_}, Maximize, Minimize")), min),
        // Ranges intersect; a list of variables shares one.
        (bounded("Equation[List[x, y], Range[-5, 5]]"), min),
        (
            bounded(&format!("{box_}, Equation[x, Range[-2, 7]], Maximize")),
            "List[44/3, List[Equation[x, -2], Equation[y, 1/3]]]",
        ),
        (lattice("Equation[Assume, LPInteger]"), int),
        (lattice("Equation[LPAssume, LPNonNegInt]"), int),
        (lattice("Equation[LPIntegerVariables, List[x, y]]"), int),
        (
            lattice("Equation[LPIntegerVariables, List[x]]"),
            "List[-130/3, List[Equation[x, 4], Equation[y, 10/3]]]",
        ),
        (
            lattice("Equation[LPBinaryVariables, List[x]], Equation[x, Range[-3, 5]]"),
            "List[-64/3, List[Equation[x, 1], Equation[y, 7/3]]]",
        ),
        // The assumption after the range, as giac's.
        (
            lattice("Equation[x, Range[-3, 5]], Equation[Assume, LPBinary]"),
            "List[-12, List[Equation[x, 1], Equation[y, 1]]]",
        ),
        (
            knapsack(
                "Equation[LPBinaryVariables, List[x1, x2, x3, x4]], Equation[LPMaximize, True]",
            ),
            "List[21, List[Equation[x1, 0], Equation[x2, 1], Equation[x3, 1], Equation[x4, 1]]]",
        ),
        (
            knapsack(
                "Equation[LPBinaryVariables, List[x1, x2]], Equation[x3, Range[0, 1]], \
                 Equation[x4, Range[0, 1]], Equation[LPMaximize, True]",
            ),
            "List[22, List[Equation[x1, 1], Equation[x2, 1], Equation[x3, 1/2], Equation[x4, 0]]]",
        ),
    ]);
    let at = |q: &str| format!("List[{q}, List[Equation[x, {q}]]]");
    let (zero, minus_one) = (at("0"), at("-1"));
    check(&[
        (
            "LPSolve[x, List[GreaterEqual[x, -1]], Equation[LPNonNegative, True]]",
            zero.as_str(),
        ),
        (
            "LPSolve[x, List[GreaterEqual[x, -1]], Equation[LPNonNegative, False]]",
            minus_one.as_str(),
        ),
        (
            "LPSolve[x, List[GreaterEqual[x, -1]], Equation[LPAssume, LPNonNegative]]",
            zero.as_str(),
        ),
        (
            "LPSolve[x, List[GreaterEqual[x, -1/2]], Equation[Assume, LPBinary]]",
            zero.as_str(),
        ),
    ]);
    stays(&[
        bounded(&format!("{box_}, Equation[LPMaximize, 1]")),
        bounded(&format!("{box_}, Equation[Maximize, x]")),
        bounded("Equation[x, Range[-5, 5]], Equation[z, Range[0, 1]]"),
        bounded("Equation[List[x, 1], Range[-5, 5]]"),
        bounded("foo"),
        bounded("Equation[LPSimplex, True]"),
        lattice("Equation[Assume, foo]"),
        lattice("Equation[LPIntegerVariables, List[1]]"),
        lattice("Equation[LPBinaryVariables, List[x, z]]"),
        lattice("Equation[LPNonNegative, 1]"),
        // Branch and bound's keys, until T3.
        lattice("Equation[Assume, LPInteger], Equation[LPNodeSelect, LPDepthFirst]"),
        lattice("Equation[Assume, LPInteger], Equation[LPDepthLimit, 1]"),
    ]);
}

/// A word alone is ignored, as giac's, with one warning that names it:
/// `lp_integer` gives the relaxation.
#[test]
fn words_alone() {
    let ctx = Context::new();
    for (word, name) in [
        ("LPInteger", "lp_integer"),
        ("LPSimplex", "lp_simplex"),
        ("LPNonNegInt", "lp_nonnegint"),
    ] {
        let e = read(&lattice(word));
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        let expected = ctx
            .eval(&read("List[-47, List[Equation[x, 9/2], Equation[y, 7/2]]]"))
            .unwrap();
        assert_eq!(found.unwrap(), expected, "{word}");
        assert_eq!(report.warnings.len(), 1, "{word}: {:?}", report.warnings);
        let text = format!("{:?}", report.warnings[0]);
        assert!(text.contains(name), "{word}: {text}");
    }
}

/// giac's `lp_method`, `lp_verbose` and `lp_presolve` change nothing, and
/// decimals are made exact (CAS-601), giac's answers (probe1, probe5).
#[test]
fn methods() {
    let program = |o: &str| {
        format!(
            "LPSolve[Plus[Times[2, x], y, Times[-1, z], 4], List[LessEqual[x, 1], \
             GreaterEqual[y, 2], Equation[Plus[x, Times[3, y], Times[-1, z]], 2], \
             LessEqual[Plus[Times[2, x], Times[-1, y], z], 8], \
             LessEqual[Plus[Times[-1, x], y], 5]], {o}]"
        )
    };
    let answer = "List[-4, List[Equation[x, 0], Equation[y, 5], Equation[z, 13]]]";
    // 2.5x + y under [x + y ≥ 1.5, x - y ≤ 0.5].
    let decimals = |o: &str| {
        format!(
            "LPSolve[Plus[Times[0x1.4p+1, x], y], List[GreaterEqual[Plus[x, y], 0x1.8p+0], \
             LessEqual[Plus[x, Times[-1, y]], 0x1p-1]], {o}]"
        )
    };
    let decimal_answer = "List[3/2, List[Equation[x, 0], Equation[y, 3/2]]]";
    check_owned(&[
        (program("Equation[LPMethod, LPSimplex]"), answer),
        (program("Equation[LPMethod, LPInteriorPoint]"), answer),
        (program("Equation[LPMethod, Exact]"), answer),
        (program("Equation[LPMethod, float]"), answer),
        (program("Equation[LPPresolve, False]"), answer),
        (program("Equation[LPVerbose, True]"), answer),
        (program("LPVerbose, LPPresolve"), answer),
        (decimals("Equation[Assume, LPNonNegative]"), decimal_answer),
        (
            decimals("Equation[Assume, LPNonNegative], Equation[LPMethod, float]"),
            decimal_answer,
        ),
        (
            decimals("Equation[x, Range[0, 0x1.8p-1]], Equation[y, Range[0, Infinity]]"),
            decimal_answer,
        ),
    ]);
    stays(&[
        program("Equation[LPMethod, foo]"),
        program("Equation[LPMethod, True]"),
        program("Equation[LPVerbose, 1]"),
        program("Equation[LPPresolve, x]"),
        // A ball Exact leaves, past the bits of a number.
        "LPSolve[Times[0x1p+100000000, x], List[GreaterEqual[x, 1]]]".to_owned(),
        "LPSolve[x, List[GreaterEqual[x, 1]], Equation[x, Range[0, 0x1p+100000000]]]".to_owned(),
    ]);
    // No warning for the words that do nothing.
    let ctx = Context::new();
    let e = read(&program("LPVerbose, LPPresolve"));
    let (_, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);
}
