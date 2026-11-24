//! `simulated_annealing` against giac (CAS-610), on the closed form first
//! (C0006): `x² + sin(20x)/5` from 2, its minimum `x* = −0.0766237535726`
//! (`f* = −0.193981966498`). giac unseeded and seeded with 42 and 7 in the
//! same input (C0014), as caseval prints them, each within `10⁻³` of `x*`;
//! ours under three seeds within `10⁻³` of `x*` too, never compared with
//! giac's (C0044: the draws differ). In its own file: giac's state is
//! shared.

use cas_core::{Expr, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// The last value of giac's printed sequence, `srand`'s seed before it.
fn last(answer: &str) -> &str {
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => start = i + 1,
            _ => {}
        }
    }
    &answer[start..]
}

/// A number as a double: a decimal's midpoint.
fn double(e: &Expr) -> Option<f64> {
    match e.view() {
        View::Ball(b) => Some(b.mid().to_f64()),
        _ => None,
    }
}

const ROW: &str = "simulated_annealing(2,x->x^2+sin(20x)/5,(x,y)->abs(x-y),(x,m)->rand(x-m,x+m),\
                   [1,0.1,1.002,1e-6],[5,20],1)";
/// The minimiser of `x² + sin(20x)/5`.
const X_STAR: f64 = -0.076_623_753_572_6;

/// giac's three answers (unseeded `-0.0766196874902`, seeded with 42
/// `-0.0766106601804`, with 7 `-0.0766229843721`) within `10⁻³` of `x*`.
#[test]
fn giac_near_the_minimum() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        ROW.to_owned(),
        format!("srand(42);{ROW}"),
        format!("srand(7);{ROW}"),
    ] {
        let text = giac.eval(&input).unwrap_or_default();
        let x = last(&text).trim().parse::<f64>().ok();
        if !x.is_some_and(|x| (x - X_STAR).abs() < 1e-3) {
            wrong.push(format!("{input}: giac prints {text}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Ours, seeded in the same input under three seeds, within `10⁻³` of
/// `x*`.
#[test]
fn ours_near_the_minimum() {
    let mut wrong = Vec::new();
    for seed in [1, 42, 7] {
        let input = format!("[srand({seed}),{ROW}]");
        let found = ours(&input).ok();
        let x = found
            .as_ref()
            .and_then(|e| e.args().get(1))
            .and_then(double);
        if !x.is_some_and(|x| (x - X_STAR).abs() < 1e-3) {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
