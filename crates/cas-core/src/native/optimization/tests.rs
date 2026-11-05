//! Tests of the optimization builtins (CAS-415 …): giac's examples, higher
//! derivatives, saddles, constraints, and what stays.

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
            "NLPSolve[Power[Plus[1, x], 2], Equation[nlp_nonnegative, True]]",
            1.0,
            &[("x", 0.0)],
        ),
        // Hock and Schittkowski's problem 19, from its initial point.
        (
            "NLPSolve[Plus[Power[Plus[-10, x1], 3], Power[Plus[-20, x2], 3]], \
             List[GreaterEqual[Plus[Power[Plus[-5, x1], 2], Power[Plus[-5, x2], 2]], 100], \
             LessEqual[Plus[Power[Plus[-6, x1], 2], Power[Plus[-5, x2], 2]], 8281/100]], \
             Equation[nlp_initialpoint, List[Equation[x1, 201/10], Equation[x2, 146/25]]]]",
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
    // Unbounded, infeasible, with a parameter: the call stays.
    for input in [
        "NLPSolve[Power[x, 3], List[]]",
        "NLPSolve[x, List[GreaterEqual[x, 1], LessEqual[x, 0]]]",
        "NLPSolve[Times[a, Power[x, 2]], List[]]",
    ] {
        check(&[(input, input)]);
    }
}
