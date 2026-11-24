//! Tests of `FrankWolfe` (CAS-609): giac's probe rows against their closed
//! minima, the `exact` and `rand` options, the `maxiter` warning, and calls
//! that stay.

use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr, Rational};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// A number as a double.
fn double(e: &Expr) -> Option<f64> {
    if let View::Ball(b) = e.view() {
        return (b.rad().to_float()?.to_f64() == 0.0).then(|| b.mid().to_f64());
    }
    let q: Rational = e.to_rational()?;
    Some(
        cas_num::ball::Float::from_rational(&q, 53, cas_num::ball::Round::Nearest)
            .0
            .to_f64(),
    )
}

/// The answer's coordinates (exact numbers) and the warnings; `None` for a
/// call that stays.
fn run(ctx: &Context, input: &str) -> (Option<Vec<f64>>, usize) {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    let found = found.unwrap_or_else(|err| panic!("{input}: {err:?}"));
    if found.head() == e.head() {
        return (None, report.warnings.len());
    }
    let x = found
        .args()
        .iter()
        .map(double)
        .collect::<Option<Vec<_>>>()
        .unwrap_or_else(|| panic!("{input}: exact coordinates, not {found:?}"));
    (Some(x), report.warnings.len())
}

const ROW1: &str = "Plus[Times[1/3, Power[Plus[x1, -2/3], 2]], Times[1/5, Power[Plus[x2, -2/3], 2]]], \
                    List[List[1, 0], List[0, 1], List[-1, 0], List[0, -1]]";
const ROW2: &str = "Plus[Power[Plus[x1, -1], 2], Power[Plus[x2, -1], 2]], \
                    List[List[0, 0], List[1, 0], List[0, 1]]";
const ROW3: &str = "Plus[Power[x1, 2], Power[x2, 2]], List[List[1, 1], List[2, 1], List[1, 2]]";

fn f1(x: &[f64]) -> f64 {
    (x[0] - 2.0 / 3.0).powi(2) / 3.0 + (x[1] - 2.0 / 3.0).powi(2) / 5.0
}
fn f2(x: &[f64]) -> f64 {
    (x[0] - 1.0).powi(2) + (x[1] - 1.0).powi(2)
}
fn f3(x: &[f64]) -> f64 {
    x[0] * x[0] + x[1] * x[1]
}
fn in1(x: &[f64]) -> bool {
    x[0].abs() + x[1].abs() <= 1.0 + 1e-12
}
fn in2(x: &[f64]) -> bool {
    x[0] >= -1e-12 && x[1] >= -1e-12 && x[0] + x[1] <= 1.0 + 1e-12
}
fn in3(x: &[f64]) -> bool {
    x[0] >= 1.0 - 1e-12 && x[1] >= 1.0 - 1e-12 && x[0] + x[1] <= 3.0 + 1e-12
}

type Row = (&'static str, fn(&[f64]) -> f64, f64, fn(&[f64]) -> bool);

const ROWS: [Row; 3] = [
    (ROW1, f1, 1.0 / 72.0, in1),
    (ROW2, f2, 0.5, in2),
    (ROW3, f3, 2.0, in3),
];

/// Each row with `options`: an answer in the hull with `f − f* ≤ tol`, no
/// warning.
fn check(options: &str, tol: f64) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (row, f, min, inside) in ROWS {
        let input = format!("FrankWolfe[{row}{options}]");
        match run(&ctx, &input) {
            (Some(x), 0) if x.len() == 2 && inside(&x) && f(&x) - min <= tol => {}
            other => wrong.push(format!("{input}: {other:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's probe rows (probe3.out 19-21), giac's defaults.
#[test]
fn closed_forms() {
    check("", 1e-3);
}

/// `tol = 2^−14`, `maxiter`, and a start `x0`.
#[test]
fn options() {
    check(", 0x1p-14", 6.103_515_625e-5);
    check(", 0x1p-10, 5000", 9.765_625e-4);
    let ctx = Context::new();
    let (x, warnings) = run(&ctx, &format!("FrankWolfe[{ROW2}, List[0x1p-2, 0x1p-2]]"));
    let x = x.unwrap_or_else(|| panic!("{warnings} warnings"));
    assert!(in2(&x) && f2(&x) - 0.5 <= 1e-3 && warnings == 0, "{x:?}");
}

/// `exact` (global search) and `exact = 10` (local).
#[test]
fn exact() {
    check(", Exact", 1e-3);
    check(", Equation[Exact, 10]", 1e-3);
}

/// `rand = 3`: the best of three random starts, by value.
#[test]
fn random_starts() {
    check(", Equation[Random, 3]", 1e-3);
}

/// `maxiter = 3` on rows 1 and 2: no certified answer, a warning, the call
/// stays (giac answers its "Unable to minimize" list).
#[test]
fn maxiter_warns() {
    let ctx = Context::new();
    for row in [ROW1, ROW2] {
        assert_eq!(
            run(&ctx, &format!("FrankWolfe[{row}, 3]")),
            (None, 1),
            "{row}"
        );
    }
}

/// Shapes that leave the call, without a warning.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "FrankWolfe[Plus[x1, x2]]".to_owned(),
        "FrankWolfe[Power[x, 2], List[List[0], List[1]]]".to_owned(),
        "FrankWolfe[Plus[x1, x2], List[List[0, 0]]]".to_owned(),
        "FrankWolfe[Plus[x1, x2], x]".to_owned(),
        "FrankWolfe[Plus[x1, x2], List[List[0, 0], List[1]]]".to_owned(),
        "FrankWolfe[Plus[x1, x2], List[List[0, 0], List[1, y]]]".to_owned(),
        "FrankWolfe[Plus[x1, x2], List[List[0, 0], List[1, Infinity]]]".to_owned(),
        "FrankWolfe[List[x1, x2], List[List[0, 0], List[1, 1]]]".to_owned(),
        "FrankWolfe[Abs[Plus[x1, x2]], List[List[0, 0], List[1, 1]]]".to_owned(),
        format!("FrankWolfe[{ROW2}, 0]"),
        format!("FrankWolfe[{ROW2}, -1/1000]"),
        format!("FrankWolfe[{ROW2}, y]"),
        format!("FrankWolfe[{ROW2}, List[0, 0, 0]]"),
        format!("FrankWolfe[{ROW2}, Equation[Exact, 0]]"),
        format!("FrankWolfe[{ROW2}, Equation[Random, 1]]"),
        format!("FrankWolfe[{ROW2}, Equation[y, 2]]"),
        format!("FrankWolfe[{ROW2}, 10, 0x1p-10, Exact, Random]"),
    ] {
        let e = read(&input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}
