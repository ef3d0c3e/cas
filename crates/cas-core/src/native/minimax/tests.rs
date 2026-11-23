//! Tests of `MinimaxApproximation` (CAS-599): closed forms, the alternation
//! of the error on the reference, `rand` after a seed, what stays, and the
//! warning at the iteration limit.

use super::{Options, chebyshev, remez};
use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn double(e: &Expr) -> f64 {
    if let View::Ball(b) = e.view() {
        return b.mid().to_f64();
    }
    let q = e.to_rational().unwrap_or_else(|| panic!("a number: {e:?}"));
    super::to_f64(&q)
}

/// The coefficients of a polynomial in `x`, the constant first.
fn coefficients(p: &Expr) -> Vec<f64> {
    let x = Expr::symbol("x");
    let terms = match p.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => p.args().to_vec(),
        _ => vec![p.clone()],
    };
    let mut out = Vec::new();
    for t in &terms {
        let factors = match t.head().and_then(Expr::as_symbol) {
            Some(Symbol::TIMES) => t.args().to_vec(),
            _ => vec![t.clone()],
        };
        let (mut c, mut k) = (1.0, 0);
        for f in &factors {
            if *f == x {
                k = 1;
            } else if f.head().and_then(Expr::as_symbol) == Some(Symbol::POWER) && f.args()[0] == x
            {
                k = usize::try_from(f.args()[1].as_i64().unwrap()).unwrap();
            } else {
                c *= double(f);
            }
        }
        if out.len() <= k {
            out.resize(k + 1, 0.0);
        }
        out[k] += c;
    }
    out
}

/// The error ball's midpoint and radius, and `p`'s coefficients.
fn answer(ctx: &Context, input: &str) -> (f64, f64, Vec<f64>) {
    let found = ctx.eval(&read(input)).unwrap();
    let [error, p] = found.args() else {
        panic!("{input}: {found:?}");
    };
    let View::Ball(b) = error.view() else {
        panic!("{input}: {found:?}");
    };
    let radius = b.rad().to_float().unwrap().to_f64();
    (b.mid().to_f64(), radius, coefficients(p))
}

/// The closed forms: coefficients within 1e−9, the error ball within 1e−6
/// of the value, relatively.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let e = std::f64::consts::E;
    let level = (1.0 - (e - 1.0) * (1.0 - (e - 1.0).ln())) / 2.0;
    for (input, error, p) in [
        (
            "MinimaxApproximation[Power[x, 3], Equation[x, Range[-1, 1]], 2]",
            0.25,
            vec![0.0, 0.75],
        ),
        (
            "MinimaxApproximation[Abs[x], Equation[x, Range[-1, 1]], 2]",
            0.125,
            vec![0.125, 0.0, 1.0],
        ),
        (
            "MinimaxApproximation[Abs[x], Equation[x, Range[-1, 1]], 3]",
            0.125,
            vec![0.125, 0.0, 1.0],
        ),
        (
            "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 1]",
            level,
            vec![1.0 - level, e - 1.0],
        ),
    ] {
        let (mid, radius, found) = answer(&ctx, input);
        assert!(
            (mid - error).abs() <= 1e-6 * error && radius <= 1e-6 * error,
            "{input}: {mid} ± {radius}"
        );
        assert!(found.len() <= p.len() + 1, "{input}: {found:?}");
        for k in 0..found.len().max(p.len()) {
            let (u, v) = (
                found.get(k).copied().unwrap_or(0.0),
                p.get(k).copied().unwrap_or(0.0),
            );
            assert!((u - v).abs() <= 1e-9, "{input}: {found:?}");
        }
    }
    // |x| − x² − 1/8 alternates at −1, −1/2, 0, 1/2, 1.
    let (_, _, p) = answer(
        &ctx,
        "MinimaxApproximation[Abs[x], Equation[x, Range[-1, 1]], 2]",
    );
    for (i, x) in [-1.0, -0.5, 0.0, 0.5, 1.0_f64].into_iter().enumerate() {
        let value: f64 = p.iter().rev().fold(0.0, |s, c| s * x + c);
        let sign = if i % 2 == 0 { -1.0 } else { 1.0 };
        assert!((x.abs() - value - sign * 0.125).abs() <= 1e-9, "{x}: {p:?}");
    }
}

/// A function on `[−1, 1]` and the degrees.
type Case = (&'static dyn Fn(f64) -> f64, usize, usize);

/// The error alternates in sign at the `N` reference points, `|e|` equal
/// within `2^−20`.
#[test]
fn reference_alternates() {
    let cases: [Case; 4] = [
        (&|t: f64| t.abs(), 2, 0),
        (&|t: f64| t.powi(3), 2, 0),
        (&|t: f64| f64::midpoint(t, 1.0).exp(), 3, 0),
        (&|t: f64| f64::midpoint(t, 1.0).sin(), 2, 1),
    ];
    for (k, (f, n, m)) in cases.into_iter().enumerate() {
        let r = remez(f, n, m, &Options::default(), &mut || 0).expect("converges");
        assert_eq!(r.reference.len(), n + m + 2, "case {k}");
        assert!(r.iterations <= Options::default().limit, "case {k}");
        let e: Vec<f64> = r
            .reference
            .iter()
            .map(|&t| f(t) - chebyshev(&r.p, t) / chebyshev(&r.q, t))
            .collect();
        for w in e.windows(2) {
            assert!(w[0] * w[1] < 0.0, "case {k}: {e:?}");
        }
        for v in &e {
            assert!(
                (v.abs() - r.emax).abs() <= 2f64.powi(-20) * r.emax,
                "case {k}: {e:?}"
            );
        }
    }
}

/// `rand = k` and `rand` draw from the generator: the same after the same
/// seed.
#[test]
fn rand_after_a_seed() {
    let ctx = Context::new();
    for input in [
        "MinimaxApproximation[Times[x, Exp[Times[-1, x]]], Equation[x, Range[0, 4]], List[2, 1], \
         Equation[Random, 3]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 3, Random]",
    ] {
        let mut found = Vec::new();
        for _ in 0..2 {
            ctx.eval(&read("SeedRandom[7]")).unwrap();
            found.push(ctx.eval(&read(input)).unwrap());
        }
        assert_eq!(
            found[0].head().and_then(Expr::as_symbol),
            Some(Symbol::LIST),
            "{input}: {:?}",
            found[0]
        );
        assert_eq!(found[0], found[1], "{input}");
    }
}

/// Options, threshold and limit, and bad shapes: the call stays.
#[test]
fn options_and_what_stays() {
    let ctx = Context::new();
    let (_, _, p) = answer(
        &ctx,
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Equation[Threshold, 2], \
         Equation[Limit, 50]]",
    );
    assert_eq!(p.len(), 3, "{p:?}");
    for input in [
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 0]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], -2]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], List[1, 0]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], List[1, 2, 3]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[1, 0]], 2]",
        "MinimaxApproximation[Exp[x], Equation[1, Range[0, 1]], 2]",
        "MinimaxApproximation[Exp[x], Equation[x, 1], 2]",
        "MinimaxApproximation[Exp[x], x, 2]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Equation[Random, 0]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Equation[Limit, 0]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Equation[Threshold, 1]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Equation[Foo, 2]]",
        "MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 2, Foo]",
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// Not converged within `limit = 1`: a warning, and the call stays.
#[test]
fn limit_warns() {
    let ctx = Context::new();
    let e = read("MinimaxApproximation[Exp[x], Equation[x, Range[0, 1]], 3, Equation[Limit, 1]]");
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head());
    assert!(!report.warnings.is_empty());
}

/// `Min` through `Abs`: giac's `minimax(min(1/cosh(3*sin(x)),sin(9*x/10)),
/// x=-3..4,30)` evaluates, its error within `10⁻³` of `0.0319105` (a
/// 200001-point sampling gives `0.0319148`), relatively.
#[test]
fn min_compiles() {
    let ctx = Context::new();
    let input = "MinimaxApproximation[Min[Power[Cosh[Times[3, Sin[x]]], -1], \
                 Sin[Times[9, Power[10, -1], x]]], Equation[x, Range[-3, 4]], 30]";
    let (mid, _, p) = answer(&ctx, input);
    assert!((mid - 0.031_910_5).abs() <= 1e-3 * 0.031_910_5, "{mid}");
    assert!(p.len() <= 31, "{p:?}");
}
