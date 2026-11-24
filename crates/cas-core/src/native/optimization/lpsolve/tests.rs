//! Tests of the linear programs (CAS-417, CAS-418, CAS-600 to CAS-603):
//! giac's examples and options, unbounded and infeasible programs,
//! integers, binaries and decimals, branch and bound's strategies against
//! enumeration, its limits and gap.

use cas_num::{Integer, Rational};

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

/// A test program over `x1, x2, …`: the objective's coefficients and
/// constant, the rows `(a, relation, (numerator, denominator))`, each
/// variable's box, which are integers, and whether it maximizes.
struct Case {
    c: &'static [i64],
    constant: i64,
    rows: &'static [(&'static [i64], &'static str, (i64, i64))],
    boxes: &'static [(i64, i64)],
    integer: &'static [bool],
    maximum: bool,
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// `Σ a_j x_j`, plus `k` when not 0, in full form.
fn sum(a: &[i64], k: i64) -> String {
    let mut terms: Vec<String> = a
        .iter()
        .enumerate()
        .map(|(j, a)| format!("Times[{a}, x{}]", j + 1))
        .collect();
    if k != 0 {
        terms.push(k.to_string());
    }
    format!("Plus[{}]", terms.join(", "))
}

impl Case {
    /// The call with the boxes and the integers, then `options`.
    fn call(&self, boxes: &[(i64, i64)], integer: bool, options: &str) -> String {
        let rows: Vec<String> = self
            .rows
            .iter()
            .map(|(a, rel, (n, d))| format!("{rel}[{}, {n}/{d}]", sum(a, 0)))
            .collect();
        let mut o: Vec<String> = boxes
            .iter()
            .enumerate()
            .map(|(j, (l, h))| format!("Equation[x{}, Range[{l}, {h}]]", j + 1))
            .collect();
        let integers: Vec<String> = (0..self.c.len())
            .filter(|&j| integer && self.integer[j])
            .map(|j| format!("x{}", j + 1))
            .collect();
        if !integers.is_empty() {
            o.push(format!(
                "Equation[LPIntegerVariables, List[{}]]",
                integers.join(", ")
            ));
        }
        if self.maximum {
            o.push("Maximize".to_owned());
        }
        if !options.is_empty() {
            o.push(options.to_owned());
        }
        format!(
            "LPSolve[{}, List[{}], {}]",
            sum(self.c, self.constant),
            rows.join(", "),
            o.join(", ")
        )
    }

    /// Whether `x` keeps the rows, the boxes and the integers.
    fn feasible(&self, x: &[Rational]) -> bool {
        let rows = self.rows.iter().all(|(a, rel, (n, d))| {
            let ax = a
                .iter()
                .zip(x)
                .fold(Rational::ZERO, |s, (a, v)| &s + &(&q(*a, 1) * v));
            let b = q(*n, *d);
            match *rel {
                "LessEqual" => ax <= b,
                "GreaterEqual" => ax >= b,
                _ => ax == b,
            }
        });
        let boxes = x
            .iter()
            .zip(self.boxes)
            .zip(self.integer)
            .all(|((v, (l, h)), is)| *v >= q(*l, 1) && *v <= q(*h, 1) && (!is || v.is_integer()));
        rows && boxes
    }

    fn value(&self, x: &[Rational]) -> Rational {
        self.c
            .iter()
            .zip(x)
            .fold(q(self.constant, 1), |s, (c, v)| &s + &(&q(*c, 1) * v))
    }

    /// The best value by enumeration of the integer variables' boxes, the
    /// continuous ones by `LPSolve` with the integers fixed.
    fn best(&self, ctx: &Context) -> Rational {
        let ints: Vec<usize> = (0..self.c.len()).filter(|&j| self.integer[j]).collect();
        let mut point: Vec<i64> = ints.iter().map(|&j| self.boxes[j].0).collect();
        let mut best: Option<Rational> = None;
        loop {
            let mut boxes = self.boxes.to_vec();
            for (&j, &v) in ints.iter().zip(&point) {
                boxes[j] = (v, v);
            }
            let value = if ints.len() == self.c.len() {
                let x: Vec<Rational> = point.iter().map(|&v| q(v, 1)).collect();
                self.feasible(&x).then(|| self.value(&x))
            } else {
                let found = ctx.eval(&read(&self.call(&boxes, false, ""))).unwrap();
                found.args().first().and_then(Expr::to_rational)
            };
            if let Some(v) = value
                && best
                    .as_ref()
                    .is_none_or(|b| if self.maximum { v > *b } else { v < *b })
            {
                best = Some(v);
            }
            // The next point, the first variable fastest.
            let mut i = 0;
            while i < ints.len() && point[i] == self.boxes[ints[i]].1 {
                point[i] = self.boxes[ints[i]].0;
                i += 1;
            }
            if i == ints.len() {
                break;
            }
            point[i] += 1;
        }
        best.expect("a feasible integer point")
    }
}

/// The integer programs of [`strategies`], each in a finite box: giac's
/// lattice, the binary knapsack, `x + y` past fractions, three bounded
/// programs, and two mixed ones.
const CASES: &[Case] = &[
    Case {
        c: &[-5, -7],
        constant: 0,
        rows: &[
            (&[7, 1], "LessEqual", (35, 1)),
            (&[-1, 3], "LessEqual", (6, 1)),
        ],
        boxes: &[(-10, 10), (-10, 10)],
        integer: &[true, true],
        maximum: false,
    },
    Case {
        c: &[8, 11, 6, 4],
        constant: 0,
        rows: &[(&[5, 7, 4, 3], "LessEqual", (14, 1))],
        boxes: &[(0, 1), (0, 1), (0, 1), (0, 1)],
        integer: &[true, true, true, true],
        maximum: true,
    },
    Case {
        c: &[1, 1],
        constant: 0,
        rows: &[
            (&[1, 0], "GreaterEqual", (1, 3)),
            (&[0, 1], "GreaterEqual", (1, 2)),
        ],
        boxes: &[(-5, 5), (-5, 5)],
        integer: &[true, true],
        maximum: false,
    },
    Case {
        c: &[5, 4, 3],
        constant: 0,
        rows: &[
            (&[2, 3, 1], "LessEqual", (5, 1)),
            (&[4, 1, 2], "LessEqual", (11, 1)),
            (&[3, 4, 2], "LessEqual", (8, 1)),
        ],
        boxes: &[(0, 3), (0, 3), (0, 3)],
        integer: &[true, true, true],
        maximum: true,
    },
    Case {
        c: &[3, 2, 4],
        constant: 0,
        rows: &[
            (&[3, 2, 5], "LessEqual", (37, 2)),
            (&[2, 5, 1], "LessEqual", (11, 1)),
            (&[1, -1, 3], "GreaterEqual", (-2, 1)),
        ],
        boxes: &[(0, 6), (0, 6), (0, 6)],
        integer: &[true, true, true],
        maximum: true,
    },
    Case {
        c: &[-1, -1],
        constant: 0,
        rows: &[
            (&[2, 2], "LessEqual", (7, 1)),
            (&[-2, 2], "LessEqual", (3, 1)),
        ],
        boxes: &[(-4, 4), (-4, 4)],
        integer: &[true, true],
        maximum: false,
    },
    Case {
        c: &[7, -3],
        constant: 1,
        rows: &[
            (&[3, 5], "GreaterEqual", (7, 2)),
            (&[4, -1], "LessEqual", (9, 1)),
        ],
        boxes: &[(-3, 5), (-3, 5)],
        integer: &[true, true],
        maximum: false,
    },
    Case {
        c: &[8, 11, 6, 4],
        constant: 0,
        rows: &[(&[5, 7, 4, 3], "LessEqual", (14, 1))],
        boxes: &[(0, 1), (0, 1), (0, 1), (0, 1)],
        integer: &[true, true, false, false],
        maximum: true,
    },
    Case {
        c: &[-5, -7],
        constant: 0,
        rows: &[
            (&[7, 1], "LessEqual", (35, 1)),
            (&[-1, 3], "LessEqual", (6, 1)),
        ],
        boxes: &[(-10, 10), (-10, 10)],
        integer: &[true, false],
        maximum: false,
    },
];

/// Every node selection, variable selection, `lp_maxcuts` 0 and 5, and
/// `lp_heuristic` on and off (CAS-602) reach the best value found by
/// enumeration, at a point that keeps the rows, the boxes and the integers
/// (C0040: every failure, one assert).
#[test]
fn strategies() {
    let ctx = Context::new();
    let nodes = [
        "LPDepthFirst",
        "LPBreadthFirst",
        "LPBestLocalBound",
        "LPBestProjection",
        "LPHybrid",
    ];
    let vars = [
        "LPFirstFractional",
        "LPLastFractional",
        "LPMostFractional",
        "LPPseudoCost",
    ];
    let mut wrong = Vec::new();
    for case in CASES {
        let best = case.best(&ctx);
        for node in nodes {
            for var in vars {
                for cuts in [0, 5] {
                    for heuristic in ["True", "False"] {
                        let options = format!(
                            "Equation[LPNodeSelect, {node}], Equation[LPVarSelect, {var}], \
                             Equation[LPMaxCuts, {cuts}], Equation[LPHeuristic, {heuristic}]"
                        );
                        let input = case.call(case.boxes, true, &options);
                        let found = ctx.eval(&read(&input)).unwrap();
                        let point: Option<Vec<Rational>> = found.args().get(1).and_then(|p| {
                            p.args()
                                .iter()
                                .map(|e| e.args().get(1).and_then(Expr::to_rational))
                                .collect()
                        });
                        let value = found.args().first().and_then(Expr::to_rational);
                        let right = matches!((&value, &point), (Some(v), Some(x))
                            if *v == best && x.len() == case.c.len() && case.feasible(x)
                                && case.value(x) == best);
                        if !right {
                            wrong.push(format!("{input}: {found:?}, best {best:?}"));
                        }
                    }
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The binary knapsack of [`knapsack`], maximized, with the options `o`.
fn binary(o: &str) -> String {
    knapsack(&format!("Equation[Assume, LPBinary], Maximize, {o}"))
}

/// giac's `x1 + x2` over `1867x1 + 1913x2 = 3618894` in nonnegative
/// integers with the options `o`.
fn coins(o: &str) -> String {
    format!(
        "LPSolve[Plus[x1, x2], List[Equation[Plus[Times[1867, x1], Times[1913, x2]], 3618894]], \
         Equation[Assume, LPNonNegInt], {o}]"
    )
}

/// The answer of `input` and its one warning, which names `option`; `None`
/// for no warning expected.
fn warns(cases: &[(String, &str, Option<&str>)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected, option) in cases {
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        let found = found.unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        let named = match option {
            Some(name) => {
                report.warnings.len() == 1 && format!("{:?}", report.warnings[0]).contains(name)
            }
            None => report.warnings.is_empty(),
        };
        if found != expected || !named {
            wrong.push(format!("{input}: {found:?} {:?}", report.warnings));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A limit or the gap stops the search (CAS-603): the incumbent, or `{}`
/// when there is none, with a warning that names the option, as giac's
/// `[]` (probe1, probe5); the points from 1.5/T3's Python prototype.
#[test]
fn limits() {
    let point = |x: [i64; 4]| {
        format!(
            "List[Equation[x1, {}], Equation[x2, {}], Equation[x3, {}], Equation[x4, {}]]",
            x[0], x[1], x[2], x[3]
        )
    };
    let fourteen = format!("List[14, {}]", point([1, 0, 1, 0]));
    let twelve = format!("List[12, {}]", point([1, 0, 0, 1]));
    let eighteen = format!("List[18, {}]", point([1, 0, 1, 1]));
    let best = format!("List[21, {}]", point([0, 1, 1, 1]));
    let program = |o: &str| {
        format!(
            "LPSolve[Plus[Times[2, x], y, Times[-1, z], 4], List[LessEqual[x, 1], \
             GreaterEqual[y, 2], Equation[Plus[x, Times[3, y], Times[-1, z]], 2], \
             LessEqual[Plus[Times[2, x], Times[-1, y], z], 8], \
             LessEqual[Plus[Times[-1, x], y], 5]], {o}]"
        )
    };
    let (node, depth, iteration, gap) = (
        Some("lp_nodelimit"),
        Some("lp_depthlimit"),
        Some("lp_iterationlimit"),
        Some("lp_gaptolerance"),
    );
    warns(&[
        (binary("Equation[LPNodeLimit, 1]"), &fourteen, node),
        (binary("Equation[LPDepthLimit, 1]"), &fourteen, depth),
        (binary("Equation[LPIterationLimit, 17]"), &twelve, iteration),
        (binary("Equation[LPIterationLimit, 1]"), "List[]", iteration),
        (binary("Equation[LPGapTolerance, 1/2]"), &eighteen, gap),
        (binary("Equation[LPGapTolerance, 0x1p-1]"), &eighteen, gap),
        // No limit: Infinity, 0, and a gap of 0.
        (binary("Equation[LPNodeLimit, Infinity]"), &best, None),
        (binary("Equation[LPDepthLimit, 0]"), &best, None),
        (binary("Equation[LPGapTolerance, 0]"), &best, None),
        // No integer point yet.
        (coins("Equation[LPNodeLimit, 5]"), "List[]", node),
        (coins("Equation[LPDepthLimit, 2]"), "List[]", depth),
        (coins("Equation[LPIterationLimit, 10]"), "List[]", iteration),
        // A program with no integer: giac's [] too.
        (
            program("Equation[LPIterationLimit, 1]"),
            "List[]",
            iteration,
        ),
        (
            program("Equation[LPIterationLimit, 1000]"),
            "List[-4, List[Equation[x, 0], Equation[y, 5], Equation[z, 13]]]",
            None,
        ),
    ]);
}

/// `lp_timelimit` is read and changes nothing (ours has no clock); the
/// keys' values of another kind leave the call (CAS-602, CAS-603).
#[test]
fn keys() {
    let int = "List[-41, List[Equation[x, 4], Equation[y, 3]]]";
    let integer = |o: &str| lattice(&format!("Equation[Assume, integer], {o}"));
    warns(&[
        (integer("Equation[LPTimeLimit, 1]"), int, None),
        (integer("Equation[LPTimeLimit, Infinity]"), int, None),
        (integer("Equation[LPMaxCuts, Infinity]"), int, None),
        (integer("Equation[LPNodeSelect, LPHybrid]"), int, None),
        (
            integer("Equation[LPVarSelect, LPMostFractional]"),
            int,
            None,
        ),
        (integer("Equation[LPHeuristic, False]"), int, None),
    ]);
    stays(&[
        integer("Equation[LPNodeSelect, foo]"),
        integer("Equation[LPNodeSelect, LPPseudoCost]"),
        integer("Equation[LPVarSelect, LPDepthFirst]"),
        integer("Equation[LPVarSelect, 1]"),
        integer("Equation[LPMaxCuts, -1]"),
        integer("Equation[LPMaxCuts, 1/2]"),
        integer("Equation[LPNodeLimit, x]"),
        integer("Equation[LPDepthLimit, 0x1.8p+0]"),
        integer("Equation[LPIterationLimit, Times[-1, Infinity]]"),
        integer("Equation[LPTimeLimit, foo]"),
        integer("Equation[LPHeuristic, 1]"),
        integer("Equation[LPGapTolerance, -1/2]"),
        integer("Equation[LPGapTolerance, x]"),
    ]);
}
