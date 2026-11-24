//! `frank_wolfe` against giac (CAS-609): closed forms first (C0006), our
//! point in the hull of `A` with `f(ours) − f* ≤ tol`, giac's point within
//! `tol` of the minimum by value too, giac's defaults and with `exact`; a
//! `rand=5` row on ours only (C0044: giac's starts are random); a
//! departure, a `maxiter` too small: giac answers its "Unable to minimize"
//! list as caseval prints it (C0007), ours warns and stays `FrankWolfe[…]`.

use cas_core::{Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// A number as a double: a decimal's midpoint, or a small rational.
fn double(e: &Expr) -> Option<f64> {
    if let View::Ball(b) = e.view() {
        return Some(b.mid().to_f64());
    }
    let q = e.to_rational()?;
    let small = |k: &cas_core::Integer| k.to_i64().and_then(|k| i32::try_from(k).ok());
    Some(f64::from(small(q.numer())?) / f64::from(small(q.denom())?))
}

/// A list of two numbers as doubles.
fn point(e: &Expr) -> Option<[f64; 2]> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) {
        return None;
    }
    match e.args() {
        [a, b] => Some([double(a)?, double(b)?]),
        _ => None,
    }
}

const ROW1: &str = "(x1-2/3)^2/3+(x2-2/3)^2/5,[[1,0],[0,1],[-1,0],[0,-1]]";
const ROW2: &str = "(x1-1)^2+(x2-1)^2,[[0,0],[1,0],[0,1]]";
const ROW3: &str = "x1^2+x2^2,[[1,1],[2,1],[1,2]]";

type Row = (&'static str, fn([f64; 2]) -> f64, f64, fn([f64; 2]) -> bool);

/// The rows: `f`, its minimum on the hull (`1/72` at `[13/24, 11/24]`,
/// `1/2` at `[1/2, 1/2]`, `2` at `[1, 1]`), the hull.
const ROWS: [Row; 3] = [
    (
        ROW1,
        |x| (x[0] - 2.0 / 3.0).powi(2) / 3.0 + (x[1] - 2.0 / 3.0).powi(2) / 5.0,
        1.0 / 72.0,
        |x| x[0].abs() + x[1].abs() <= 1.0 + 1e-12,
    ),
    (
        ROW2,
        |x| (x[0] - 1.0).powi(2) + (x[1] - 1.0).powi(2),
        0.5,
        |x| x[0] >= -1e-12 && x[1] >= -1e-12 && x[0] + x[1] <= 1.0 + 1e-12,
    ),
    (
        ROW3,
        |x| x[0] * x[0] + x[1] * x[1],
        2.0,
        |x| x[0] >= 1.0 - 1e-12 && x[1] >= 1.0 - 1e-12 && x[0] + x[1] <= 3.0 + 1e-12,
    ),
];

/// Each row with `options`: ours in the hull and within `10⁻³` of the
/// minimum by value; giac's within it too when `giac_too`.
fn check(options: &str, giac_too: bool) -> Vec<String> {
    let tol = 1e-3;
    let mut wrong = Vec::new();
    for (row, f, min, inside) in ROWS {
        let input = format!("frank_wolfe({row}{options})");
        let found = ours(&input).ok();
        match found.as_ref().and_then(point) {
            Some(x) if inside(x) && f(x) - min <= tol => {}
            _ => wrong.push(format!("{input}: ours {found:?}")),
        }
        if !giac_too {
            continue;
        }
        if let Some(giac) = giac() {
            let theirs = giac.value(&input);
            let ok = theirs
                .as_ref()
                .ok()
                .and_then(point)
                .is_some_and(|x| f(x) - min <= tol);
            if !ok {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
    }
    wrong
}

/// giac's probe rows (probe3.out 19-21), giac's defaults.
#[test]
fn closed_forms() {
    let wrong = check("", true);
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The line search of `find_minimum` (`exact`).
#[test]
fn exact() {
    let wrong = check(",exact", true);
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Five random starts, ours only (giac's starts are random too).
#[test]
fn random_starts() {
    let wrong = check(",rand=5", false);
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departure: `maxiter = 3` leaves no certified answer; giac's list as
/// caseval prints it, ours stays `FrankWolfe[…]`.
#[test]
fn departures() {
    let mut wrong = Vec::new();
    for (input, answer) in [
        (
            format!("frank_wolfe({ROW2},3)"),
            "[\"Unable to minimize at given precision, last value\",[0.476908050292,0.354109911339]]",
        ),
        (
            format!("frank_wolfe({ROW1},3)"),
            "[\"Unable to minimize at given precision, last value\",[0.485764101255,0.110312297746]]",
        ),
    ] {
        if let Some(giac) = giac() {
            let text = giac.eval(&input).unwrap_or_default();
            if text != answer {
                wrong.push(format!("{input}: giac prints {text}, recorded {answer}"));
            }
        }
        let found = ours(&input);
        if found.as_ref().ok().and_then(Expr::head) != Some(&Expr::from(Symbol::FRANK_WOLFE)) {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
