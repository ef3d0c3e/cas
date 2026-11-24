//! Tests of the optimization builtins (CAS-415 …): giac's examples, higher
//! derivatives, saddles, constraints, and what stays.

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
fn extrema() {
    check(&[
        (
            "Extrema[Plus[Power[x, 3], Times[-3, x]], x]",
            "List[List[1], List[-1]]",
        ),
        (
            "Extrema[Plus[Power[x, 2], Power[y, 2], Times[-2, x]], List[x, y]]",
            "List[List[List[1, 0]], List[]]",
        ),
        // A saddle at the origin and a minimum at (1, 1).
        (
            "Extrema[Plus[Power[x, 3], Power[y, 3], Times[-3, x, y]], List[x, y]]",
            "List[List[List[1, 1]], List[]]",
        ),
        // Higher derivatives: x⁴ has a minimum, x³ nothing.
        ("Extrema[Power[x, 4], x]", "List[List[0], List[]]"),
        ("Extrema[Power[x, 3], x]", "List[List[], List[]]"),
        // A periodic function's critical points repeat across the range
        // (CAS-431), giac's [[3/2*pi],[1/2*pi,5/2*pi]].
        (
            "Extrema[Sin[x], Equation[x, Range[0, 10]]]",
            "List[List[Times[3/2, Pi]], List[Times[1/2, Pi], Times[5/2, Pi]]]",
        ),
        (
            "Extrema[Times[-1, Power[x, 6]], x]",
            "List[List[], List[0]]",
        ),
        // A range keeps the points within it.
        (
            "Extrema[Plus[Power[x, 3], Times[-3, x]], Equation[x, Range[0, Infinity]]]",
            "List[List[1], List[]]",
        ),
        // Lagrange's multipliers: giac's examples.
        (
            "Extrema[Times[x, y, z], Equation[Plus[x, y, z], 1], List[x, y, z]]",
            "List[List[], List[List[1/3, 1/3, 1/3]]]",
        ),
        (
            "Extrema[Plus[Power[x, 2], Power[y, 2]], Equation[Times[x, y], 1], \
             List[Equation[x, Range[0, Infinity]], Equation[y, Range[0, Infinity]]]]",
            "List[List[List[1, 1]], List[]]",
        ),
        (
            "Extrema[Times[x, y], Equation[Plus[x, y], 1], List[x, y]]",
            "List[List[], List[List[1/2, 1/2]]]",
        ),
        // What stays: infinitely many critical points.
        (
            "Extrema[Times[Power[x, 2], Power[y, 2]], List[x, y]]",
            "Extrema[Times[Power[x, 2], Power[y, 2]], List[x, y]]",
        ),
        ("Extrema[x, 1]", "Extrema[x, 1]"),
    ]);
}

/// `NLPSolve`'s answer as doubles: the value and each variable's name and
/// value.
fn optimum(input: &str) -> (f64, Vec<(String, f64)>) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    let number = |e: &Expr| match ctx
        .eval(&Expr::call(crate::Symbol::N, [e.clone()]))
        .unwrap()
        .view()
    {
        crate::View::Ball(b) => b.mid().to_f64(),
        _ => panic!("{input}: {found:?}"),
    };
    let [value, point] = found.args() else {
        panic!("{input}: {found:?}");
    };
    let coordinates = point
        .args()
        .iter()
        .map(|eq| match eq.args() {
            [x, v] => (format!("{x:?}"), number(v)),
            _ => panic!("{input}: {found:?}"),
        })
        .collect();
    (number(value), coordinates)
}

/// A program in full form, its optimum's value and point.
type Program<'a> = (&'a str, f64, &'a [(&'a str, f64)]);

/// Nonlinear programs (CAS-519): local optima certified as
/// Karush–Kuhn–Tucker points, against their closed forms, under
/// inequalities (where giac's `nlpsolve` fails), equations and bounds.
#[test]
fn nonlinear_programs() {
    let third = std::f64::consts::PI / 3.0;
    let cases: [Program<'_>; 9] = [
        (
            "NLPSolve[Plus[Power[x, 2], Power[y, 2]], List[GreaterEqual[Plus[x, y], 1]]]",
            0.5,
            &[("x", 0.5), ("y", 0.5)],
        ),
        (
            "NLPSolve[Plus[Power[Plus[-1, x], 2], Power[Plus[-2, y], 2]], \
             List[LessEqual[Plus[x, y], 1]]]",
            2.0,
            &[("x", 0.0), ("y", 1.0)],
        ),
        (
            "NLPSolve[x1, List[LessEqual[Plus[Power[x1, 2], Power[x2, 2]], 1]]]",
            -1.0,
            &[("x1", -1.0), ("x2", 0.0)],
        ),
        // McCormick's function: x1 = 1/2 − π/3, x2 = −1/2 − π/3.
        (
            "NLPSolve[Plus[Sin[Plus[x1, x2]], Power[Plus[x1, Times[-1, x2]], 2], \
             Times[-3/2, x1], Times[5/2, x2], 1], Equation[x1, Range[-3/2, 4]], \
             Equation[x2, Range[-3, 3]]]",
            -(3.0_f64.sqrt()) / 2.0 - third,
            &[("x1", 0.5 - third), ("x2", -0.5 - third)],
        ),
        (
            "NLPSolve[Times[Sin[x], Cos[y]], List[Equation[x, Range[0, 3]], \
             Equation[y, Range[0, 3]]]]",
            3.0_f64.cos(),
            &[("x", std::f64::consts::FRAC_PI_2), ("y", 3.0)],
        ),
        (
            "NLPSolve[Times[x, y], List[Equation[Plus[x, y], 2]], Maximize]",
            1.0,
            &[("x", 1.0), ("y", 1.0)],
        ),
        (
            "NLPSolve[Power[Plus[-2, x], 2], List[Equation[x, Range[0, 1]]]]",
            1.0,
            &[("x", 1.0)],
        ),
        (
            "NLPSolve[Power[Plus[1, x], 2], Equation[NLPNonNegative, True]]",
            1.0,
            &[("x", 0.0)],
        ),
        // Hock and Schittkowski's problem 19, from its initial point.
        (
            "NLPSolve[Plus[Power[Plus[-10, x1], 3], Power[Plus[-20, x2], 3]], \
             List[GreaterEqual[Plus[Power[Plus[-5, x1], 2], Power[Plus[-5, x2], 2]], 100], \
             LessEqual[Plus[Power[Plus[-6, x1], 2], Power[Plus[-5, x2], 2]], 8281/100]], \
             Equation[NLPInitialPoint, List[Equation[x1, 201/10], Equation[x2, 146/25]]]]",
            -6_961.813_875_580_135,
            &[("x1", 14.095), ("x2", 0.842_960_789_215_480_2)],
        ),
    ];
    for (input, value, point) in cases {
        let (v, xs) = optimum(input);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * (1.0 + b.abs());
        assert!(close(v, value), "{input}: {v}, expected {value}");
        assert_eq!(xs.len(), point.len(), "{input}: {xs:?}");
        for ((name, a), (expected_name, b)) in xs.iter().zip(point) {
            assert_eq!(name, expected_name, "{input}");
            assert!(close(*a, *b), "{input}: {name} = {a}, expected {b}");
        }
    }
    // Unbounded, with a parameter: the call stays.
    for input in [
        "NLPSolve[Power[x, 3], List[]]",
        "NLPSolve[Times[a, Power[x, 2]], List[]]",
    ] {
        check(&[(input, input)]);
    }
    // A linear program is LPSolve's (CAS-604): infeasible, `{}`.
    check(&[(
        "NLPSolve[x, List[GreaterEqual[x, 1], LessEqual[x, 0]]]",
        "List[]",
    )]);
}

/// A program as [`Program`], its input built.
type OwnedProgram<'a> = (String, f64, &'a [(&'a str, f64)]);

/// `(x − 1)² + (y − 2)²` with `options`, in full form.
fn bowl(options: &str) -> String {
    format!("NLPSolve[Plus[Power[Plus[-1, x], 2], Power[Plus[-2, y], 2]], {options}]")
}

/// Every ball's radius in `e` at most `r` times its midpoint's magnitude,
/// or `r` below magnitude 1.
fn within(e: &Expr, r: f64) -> bool {
    match e.view() {
        crate::View::Ball(b) => {
            let radius = b.rad().to_float().map_or(f64::INFINITY, |f| f.to_f64());
            radius <= r * b.mid().to_f64().abs().max(1.0)
        }
        _ => e.args().iter().all(|a| within(a, r)),
    }
}

/// giac's options of `nlpsolve` (CAS-604, optimization.cc:6650-6832),
/// matched by their builtins: each method is ours, the words that change
/// nothing, a region's midpoint as the start, `nlp_precision` met.
#[test]
fn nonlinear_options() {
    let origin: &[(&str, f64)] = &[("x", 1.0), ("y", 2.0)];
    let box3 = "Equation[List[x, y], Range[0, 3]]";
    let mut cases: Vec<OwnedProgram<'_>> = vec![
        (bowl(box3), 0.0, origin),
        (
            bowl("Equation[x, Range[0, 3]], Equation[x, Range[-1, 2]], Equation[y, Range[0, 3]]"),
            0.0,
            origin,
        ),
        (
            bowl("List[LessEqual[Plus[x, y], 2]]"),
            0.5,
            &[("x", 0.5), ("y", 1.5)],
        ),
        (
            "NLPSolve[Plus[Power[Plus[1, x], 2], Power[Plus[-2, y], 2]], \
             Equation[Assume, NLPNonNegative]]"
                .to_owned(),
            1.0,
            &[("x", 0.0), ("y", 2.0)],
        ),
        (
            "NLPSolve[Times[-1, Plus[Power[Plus[-1, x], 2], Power[Plus[-2, y], 2]]], \
             Equation[List[x, y], Range[0, 3]], NLPMaximize]"
                .to_owned(),
            0.0,
            origin,
        ),
        (
            "NLPSolve[Times[-1, Plus[Power[Plus[-1, x], 2], Power[Plus[-2, y], 2]]], \
             Equation[List[x, y], Range[0, 3]], Equation[NLPMaximize, True]]"
                .to_owned(),
            0.0,
            origin,
        ),
        (
            bowl(
                "Equation[NLPInitialPoint, List[Equation[x, Range[0, 1]], Equation[y, Range[0, 1]]]]",
            ),
            0.0,
            origin,
        ),
        (
            bowl("Equation[NLPInitialPoint, List[Equation[x, 5], Equation[y, 5]]]"),
            0.0,
            origin,
        ),
    ];
    for option in [
        "Equation[NLPIterationLimit, 500]",
        "Equation[NLPTolerance, 1/1000]",
        "Equation[NLPPrecision, 0x1p-40]",
        "Equation[NLPPresolve, False]",
        "NLPPresolve",
        "NLPVerbose",
        "Equation[NLPVerbose, 2]",
        "Equation[NLPVerbose, True]",
        "Equation[NLPMethod, \"automatic\"]",
        "Equation[NLPMethod, \"ipt\"]",
        "Equation[NLPMethod, \"nelder-mead\"]",
        "Equation[NLPMethod, \"de\"]",
        "Equation[NLPMethod, \"cobyla\"]",
        "Equation[NLPMethod, List[\"nm\", Equation[\"size\", 1]]]",
    ] {
        cases.push((bowl(&format!("{box3}, {option}")), 0.0, origin));
    }
    let ctx = Context::new();
    for (input, value, point) in &cases {
        let (v, xs) = optimum(input);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * (1.0 + b.abs());
        assert!(close(v, *value), "{input}: {v}, expected {value}");
        assert_eq!(xs.len(), point.len(), "{input}: {xs:?}");
        for ((name, a), (expected_name, b)) in xs.iter().zip(*point) {
            assert_eq!(name, expected_name, "{input}");
            assert!(close(*a, *b), "{input}: {name} = {a}, expected {b}");
        }
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
        assert!(within(&found.unwrap(), 1e-12), "{input}");
    }
    // A precision past the working one: the certification raises its own.
    let fine = bowl(&format!("{box3}, Equation[NLPPrecision, 0x1p-120]"));
    let found = ctx.eval(&read(&fine)).unwrap();
    assert!(within(&found, 2f64.powi(-120)), "{fine}: {found:?}");
    assert!(!within(
        &ctx.eval(&read(&bowl(box3))).unwrap(),
        2f64.powi(-120)
    ));
}

/// What leaves `NLPSolve` (CAS-604): a warning that names the option or
/// method, or none for a form that is no option.
#[test]
fn nonlinear_options_that_stay() {
    let ctx = Context::new();
    let box3 = "Equation[List[x, y], Range[0, 3]]";
    for (option, named) in [
        ("Equation[NLPIterationLimit, 2]", Some("nlp_iterationlimit")),
        ("Equation[NLPMethod, \"foo\"]", Some("foo")),
        ("Equation[NLPMethod, cobyla]", Some("cobyla")),
        ("Equation[NLPMethod, List[nm]]", Some("nm")),
        ("Equation[NLPMethod, List[\"nm\", Equation[size, 1]]]", None),
        ("Equation[NLPIterationLimit, 0]", None),
        ("Equation[NLPIterationLimit, 3/2]", None),
        ("Equation[NLPPrecision, 0]", None),
        ("Equation[NLPPrecision, z]", None),
        ("Equation[NLPTolerance, -1]", None),
        ("Equation[NLPVerbose, 5]", None),
        ("Equation[NLPPresolve, 1]", None),
        ("Equation[NLPMaximize, 1]", None),
        ("Equation[Assume, LPBinary]", None),
        ("Equation[NLPNonNegative, 1]", None),
        // A builtin key that is no option is never a constraint.
        ("Equation[LPMaximize, True]", None),
        ("Equation[NLPInteger, 1]", None),
        ("NLPInteger", None),
        (
            "Equation[NLPInitialPoint, List[Equation[x, List[1, 2]], Equation[y, List[1, 2]]]]",
            None,
        ),
        (
            "Equation[NLPInitialPoint, List[Equation[x, Range[0, 1]], Equation[y, 1]]]",
            None,
        ),
        ("Equation[NLPInitialPoint, List[Equation[z, 1]]]", None),
        ("Equation[NLPIntegerVariables, List[z]]", None),
        ("Equation[NLPIntegerVariables, x]", None),
        // A nonlinear program with integer variables: not yet (T4.c2).
        ("Equation[Assume, integer]", None),
        ("Equation[Assume, NLPNonNegInt]", None),
        ("Equation[Assume, NLPBinary]", None),
        ("Equation[NLPBinaryVariables, List[x]]", None),
    ] {
        let input = bowl(&format!("{box3}, {option}"));
        let e = read(&input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        match named {
            Some(name) => {
                assert_eq!(report.warnings.len(), 1, "{input}: {:?}", report.warnings);
                let text = format!("{:?}", report.warnings[0]);
                assert!(text.contains(name), "{input}: {text}");
            }
            None => assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings),
        }
    }
}

/// A linear objective under linear constraints is `LPSolve`'s, exact
/// (CAS-604, optimization.cc:6080-6085), the point in sorted order;
/// giac's answers (T4 probes).
#[test]
fn linear_nonlinear_programs() {
    let rows = "List[LessEqual[Plus[x, Times[2, y]], 9/2], LessEqual[Plus[Times[3, x], y], 13/2]]";
    check(&[
        (
            "NLPSolve[Plus[x, y], List[LessEqual[Plus[x, Times[2, y]], 4], \
             LessEqual[Plus[Times[3, x], y], 6]], Equation[Assume, NLPNonNegative], Maximize]",
            "List[14/5, List[Equation[x, 8/5], Equation[y, 6/5]]]",
        ),
        (
            &format!(
                "NLPSolve[Plus[Times[2, x], y], {rows}, Equation[Assume, NLPNonNegInt], Maximize]"
            ),
            "List[4, List[Equation[x, 2], Equation[y, 0]]]",
        ),
        (
            &format!(
                "NLPSolve[Plus[Times[2, x], y], {rows}, Equation[List[x, y], Range[0, 10]], \
                 Equation[Assume, integer], Equation[NLPMaximize, True]]"
            ),
            "List[4, List[Equation[x, 2], Equation[y, 0]]]",
        ),
        (
            &format!(
                "NLPSolve[Plus[x, y], {rows}, Equation[Assume, NLPNonNegative], \
                 Equation[NLPIntegerVariables, List[x]], Maximize]"
            ),
            "List[11/4, List[Equation[x, 1], Equation[y, 7/4]]]",
        ),
        (
            &format!(
                "NLPSolve[Plus[Times[2, x], y], {rows}, Equation[Assume, NLPBinary], Maximize]"
            ),
            "List[3, List[Equation[x, 1], Equation[y, 1]]]",
        ),
        (
            "NLPSolve[Plus[Times[2, x], y], List[LessEqual[Plus[x, y], 3/2]], \
             Equation[x, Range[0, 5]], Equation[y, Range[0, 5]], \
             Equation[NLPBinaryVariables, List[x, y]], Maximize]",
            "List[2, List[Equation[x, 1], Equation[y, 0]]]",
        ),
        (
            "NLPSolve[Plus[Times[2, x], Times[3, y]], List[GreaterEqual[Plus[x, y], 1]], \
             Equation[Assume, NLPNonNegative]]",
            "List[2, List[Equation[x, 1], Equation[y, 0]]]",
        ),
        // Decimals made exact.
        (
            "NLPSolve[Plus[x, Times[-1, y]], List[LessEqual[Plus[x, y], 2.5], \
             GreaterEqual[x, 0], GreaterEqual[y, 0]], Equation[NLPMaximize, True]]",
            "List[5/2, List[Equation[x, 5/2], Equation[y, 0]]]",
        ),
        // LPSolve's variables come in first appearance, y then x: sorted.
        (
            "NLPSolve[y, List[GreaterEqual[Plus[x, y], 1], LessEqual[x, 1/2]], \
             Equation[Assume, NLPNonNegative]]",
            "List[1/2, List[Equation[x, 1/2], Equation[y, 1/2]]]",
        ),
        // Unbounded: LPSolve's `{Infinity, {}}`.
        (
            "NLPSolve[Plus[x, y], List[GreaterEqual[Plus[x, y], 1]], \
             Equation[Assume, NLPNonNegative], Maximize]",
            "List[Infinity, List[]]",
        ),
    ]);
}
