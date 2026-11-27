//! Tests of `FindMinimum` (CAS-608): giac's probe rows against their
//! closed forms, a minimum at an end, the tie of `sin` on `[0, 4π]`, the
//! double search alone, and calls that stay.

use super::{Stop, search};
use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr, Rational};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The answer's midpoint and radius as doubles, and the warnings.
fn run(ctx: &Context, input: &str) -> (f64, f64, usize) {
    let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
    let found = found.unwrap_or_else(|e| panic!("{input}: {e:?}"));
    let (mid, rad) = if let View::Ball(b) = found.view() {
        (
            b.mid().to_f64(),
            b.rad().to_float().map_or(f64::INFINITY, |r| r.to_f64()),
        )
    } else {
        let q: Rational = found
            .to_rational()
            .unwrap_or_else(|| panic!("{input}: a number, not {found:?}"));
        let (d, _) = cas_num::ball::Float::from_rational(&q, 53, cas_num::ball::Round::Nearest);
        (d.to_f64(), 0.0)
    };
    (mid, rad, report.warnings.len())
}

/// The root of `4x³ − 6x + 1` near `−1.3008395659`, by Newton's method.
fn quartic_minimiser() -> f64 {
    let mut x = -1.3_f64;
    for _ in 0..8 {
        x -= (4.0 * x.powi(3) - 6.0 * x + 1.0) / (12.0 * x * x - 6.0);
    }
    x
}

/// giac's probe rows (probe3.out 10-18) that compile: the closed form
/// within the answer's ball (up to the double of the closed form), its
/// radius at most `eps`, no warning.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let three_halves_pi = 3.0 * std::f64::consts::FRAC_PI_2;
    let mut wrong = Vec::new();
    for (input, closed, eps) in [
        ("FindMinimum[Cosh[x], -1, 1]", 0.0, 1e-8),
        // 1e-5 as a hex float (C0047).
        (
            "FindMinimum[Cosh[x], -1, 1, 0x1.4f8b588e368f1p-17, 20]",
            0.0,
            1e-5,
        ),
        (
            "FindMinimum[Plus[Power[x, 2], Times[-2, x]], 0, 3]",
            1.0,
            1e-8,
        ),
        (
            "FindMinimum[Power[Plus[x, -1/3], 2], 0, 1]",
            1.0 / 3.0,
            1e-8,
        ),
        (
            "FindMinimum[Function[List[x], Plus[Power[x, 4], Times[-3, Power[x, 2]], x]], -3, 3]",
            quartic_minimiser(),
            1e-8,
        ),
        ("FindMinimum[Sin[x], 0, 10]", three_halves_pi, 1e-8),
        ("FindMinimum[Abs[Plus[x, -1/2]], 0, 1]", 0.5, 1e-8),
    ] {
        let (mid, rad, warnings) = run(&ctx, input);
        if (mid - closed).abs() > rad + 1e-14 || rad > eps || warnings > 0 {
            wrong.push(format!("{input}: {mid} ± {rad}, {warnings} warnings"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `x` on `[0, 1]`: increasing, the exact `0`; `x²` on `[0, 1]`: `0` within
/// `eps`.
#[test]
fn ends() {
    let ctx = Context::new();
    assert_eq!(run(&ctx, "FindMinimum[x, 0, 1]"), (0.0, 0.0, 0));
    let (mid, rad, warnings) = run(&ctx, "FindMinimum[Times[-1, x], 0, 1]");
    assert_eq!((mid, rad, warnings), (1.0, 0.0, 0));
    let (mid, rad, warnings) = run(&ctx, "FindMinimum[Power[x, 2], 0, 1]");
    assert!(
        mid.abs() <= rad && rad <= 1e-8 && warnings == 0,
        "{mid} ± {rad}"
    );
}

/// `sin` on `[0, 4π]`: minima `−1` at `3π/2` and `7π/2`; the leftmost, with
/// a warning.
#[test]
fn tie_warns() {
    let ctx = Context::new();
    let (mid, rad, warnings) = run(&ctx, "FindMinimum[Sin[x], 0, Times[4, Pi]]");
    let closed = 3.0 * std::f64::consts::FRAC_PI_2;
    assert!((mid - closed).abs() <= rad + 1e-14, "{mid} ± {rad}");
    assert!(rad <= 1e-8, "{rad}");
    assert_eq!(warnings, 1);
}

/// The double search: giac's global queue finds the global minimum of the
/// quartic (a local one at `1.13`), and its limits.
#[test]
fn search_alone() {
    let quartic = |x: f64| x.powi(4) - 3.0 * x * x + x;
    let (arg, _) = search(&quartic, -3.0, 3.0, 1e-8, 100, true, 1000).unwrap();
    assert!((arg - quartic_minimiser()).abs() < 1e-6, "{arg}");
    let (arg, value) = search(&f64::sin, 0.0, 10.0, 1e-8, 100, true, 1000).unwrap();
    assert!(
        (arg - 3.0 * std::f64::consts::FRAC_PI_2).abs() < 1e-6,
        "{arg}"
    );
    assert!((value + 1.0).abs() < 1e-12, "{value}");
    assert_eq!(
        search(&quartic, -3.0, 3.0, 1e-8, 100, true, 1),
        Err(Stop::Limits)
    );
    assert_eq!(
        search(&f64::ln, -1.0, 1.0, 1e-8, 100, true, 1000),
        Err(Stop::NotReal)
    );
    let (arg, _) = search(&|x: f64| (x - 2.0).powi(2), 0.0, 3.0, 1e-8, 100, false, 1).unwrap();
    assert!((arg - 2.0).abs() < 1e-6, "{arg}");
}

/// Shapes that leave the call, without a warning.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "FindMinimum[Sin[x], 0]",
        "FindMinimum[Sin[x], 1, 0]",
        "FindMinimum[Sin[x], 1, 1]",
        "FindMinimum[Sin[x], 0, Infinity]",
        "FindMinimum[Sin[x], 0, y]",
        "FindMinimum[Plus[x, y], 0, 1]",
        "FindMinimum[3, 0, 1]",
        "FindMinimum[Function[List[x], 3], 0, 1]",
        "FindMinimum[Function[List[x, y], x], 0, 1]",
        "FindMinimum[Gamma[x], 1, 2]",
        "FindMinimum[Sin[x], 0, 1, 0]",
        "FindMinimum[Sin[x], 0, 1, -1]",
        "FindMinimum[Sin[x], 0, 1, 3/2]",
        "FindMinimum[Sin[x], 0, 1, y]",
        "FindMinimum[Sin[x], 0, 1, 2, 3, 4]",
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// CAS-650: an end holding a held `FindMinimum` 30 deep is refused at
/// once (`minimax::ball`), not approximated again at each level.
#[test]
fn nested_end() {
    let mut text = String::from("y");
    for _ in 0..30 {
        text = format!("FindMinimum[Sin[x], 0, {text}]");
    }
    let ctx = Context::new();
    let e = read(&text);
    let start = std::time::Instant::now();
    let (found, _) = ctx.eval_with(&e, &EvalOptions::default());
    let elapsed = start.elapsed();
    assert_eq!(found.unwrap().head(), e.head());
    assert!(elapsed.as_secs_f64() < 2.0, "{elapsed:?}");
}
