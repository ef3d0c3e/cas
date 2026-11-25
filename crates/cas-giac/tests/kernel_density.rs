//! Xcas's kernel density estimation against giac (CAS-629):
//! `kernel_density` and `kde`, in its exact, list and piecewise forms
//! (giac's answers in .run/1.7/giac/probe3.out and .run/1.7/T5/extra.out).
//! The exact form compares by its values at sample points within giac's 12
//! digits, ours exact (no ball) on exact data with an exact bandwidth.
//! giac convolves by FFT, whose lists carry ~10⁻⁹ relative noise (C0084):
//! list entries, and the piecewise form's values at sample points (its
//! breakpoints by sampling, C0049), compare within it; ours sums directly.
//! A seeded row compares the list's size (C0014). Invalid calls stay with
//! giac's message as the warning (C0022, C0077); departures: giac's text,
//! whitespace collapsed, and ours staying without a warning.

use cas_core::{Context, EvalOptions, Expr, Symbol, View, Warning};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, decimal_balls, lower, parse};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// A number's value in `f64`: a ball's midpoint, a small integer's.
fn float(e: &Expr) -> Option<f64> {
    match e.view() {
        View::Ball(b) => Some(b.mid().to_f64()),
        View::Int(i) => i32::try_from(i).ok().map(f64::from),
        _ => None,
    }
}

/// A list's entries in `f64`.
fn floats(e: &Expr) -> Option<Vec<f64>> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) {
        return None;
    }
    e.args().iter().map(float).collect()
}

/// Whether a ball appears anywhere in `e`.
fn has_ball(e: &Expr) -> bool {
    match e.view() {
        View::Ball(_) | View::ComplexBall(_) => true,
        View::Apply { head, args } => has_ball(head) || args.iter().any(has_ball),
        _ => false,
    }
}

/// Whether `ours` is within giac's FFT noise of `theirs` (C0084): each
/// entry within `10⁻⁹ max|entry| + 1.5·10⁻¹² |giac|`.
fn within_noise(mine: &[f64], theirs: &[f64]) -> bool {
    let scale = theirs.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    mine.len() == theirs.len()
        && mine
            .iter()
            .zip(theirs)
            .all(|(a, b)| (a - b).abs() <= 1e-9 * scale + 1.5e-12 * b.abs())
}

/// `f` at `var = v`, numeric, from both sides.
fn sampled(giac: &Giac, f: &str, var: &str, points: &[&str]) -> (Vec<f64>, Vec<f64>) {
    let (mut mine, mut theirs) = (Vec::new(), Vec::new());
    for v in points {
        let input = format!("evalf(subst({f},{var}={v}))");
        theirs.push(
            giac.value(&input)
                .ok()
                .as_ref()
                .and_then(float)
                .unwrap_or(f64::NAN),
        );
        mine.push(
            ours(&input)
                .ok()
                .as_ref()
                .and_then(float)
                .unwrap_or(f64::NAN),
        );
    }
    (mine, theirs)
}

/// The exact form: values at `x = 1/2, 2, 13/4` (or the row's points)
/// within giac's 12 digits; ours exact where the data and the bandwidth are.
#[test]
fn exact_forms() {
    let Some(giac) = giac() else {
        return;
    };
    let usual: &[&str] = &["1/2", "2", "13/4"];
    let rows: [(&str, &str, &[&str], bool); 9] = [
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,exact)",
            "x",
            usual,
            true,
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,exact,x)",
            "x",
            usual,
            true,
        ),
        ("kde([1,2,3,2],bandwidth=1/2,exact)", "x", usual, true),
        ("kernel_density([1,2,3,2],exact)", "x", usual, false),
        (
            "kernel_density([1.0,2.5,3.0,2.0,4.5],exact)",
            "x",
            usual,
            false,
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=gauss,exact)",
            "x",
            usual,
            false,
        ),
        (
            "kernel_density([1.0,2.5,3.0,2.0,4.5],bandwidth=gauss,exact)",
            "x",
            usual,
            false,
        ),
        (
            "kernel_density(seq(irem(k*k,37),k,1,1001),exact)",
            "x",
            &["5", "17", "61/2"],
            false,
        ),
        ("kde([1,2,3,2],bandwidth=1/2,exact,t)", "t", usual, true),
    ];
    let mut wrong = Vec::new();
    for (input, var, points, exact) in rows {
        let (mine, theirs) = sampled(giac, input, var, points);
        let close = mine
            .iter()
            .zip(&theirs)
            .all(|(a, b)| (a - b).abs() <= 1.5e-11 * b.abs());
        if !close {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
        if exact && !ours(input).is_ok_and(|e| !has_ball(&e)) {
            wrong.push(format!("{input}: ours {:?}, not exact", ours(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The list form: each entry within giac's FFT noise ([`within_noise`]).
#[test]
fn lists() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "kernel_density([1,2,3,2],bandwidth=1/4)",
        "kernel_density([1,2,3,2])",
        "kernel_density([1,2,3,2],bandwidth=1/4,bins=4)",
        "kernel_density([1,2,3,2],bandwidth=1/4,x)",
        "kernel_density([1,2,3,2],bandwidth=1/4,exact,output=list)",
        "kernel_density([1,2,3,2],bandwidth=1/4,bins=3,output=list)",
        "kernel_density([1,2,3,2],bandwidth=1/4,range=0..4,bins=4)",
        "kernel_density([1,2,3,2],bandwidth=select)",
        "kernel_density([1,2,3,2],bandwidth=gauss)",
        "kernel_density([1.0,2.5,3.0,2.0,4.5],bandwidth=normal,bins=5)",
        "kernel_density([1,2,3,2],bandwidth=1/4,bins=2)",
        "kernel_density([1,2,3,2],bandwidth=1/4,range=[0,4],bins=4)",
        "kernel_density([1,2,3,2],bins=7)",
        "kernel_density([1.0,2.5,3.0,2.0,4.5],bins=6)",
        "kernel_density(seq(irem(k*k,37),k,1,1001),bins=12)",
        "kernel_density(seq(evalf(sin(k)),k,1,1500),bins=9)",
        "kernel_density(seq(evalf(sin(k)),k,1,1000),bins=9)",
    ] {
        let theirs = giac.value(input).ok().as_ref().and_then(floats);
        let found = ours(input).ok();
        let mine = found.as_ref().and_then(floats);
        let balls = found
            .as_ref()
            .is_some_and(|e| e.args().iter().all(|x| matches!(x.view(), View::Ball(_))));
        match (&mine, &theirs) {
            (Some(m), Some(t)) if balls && within_noise(m, t) => {}
            _ => wrong.push(format!("{input}: ours {found:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The piecewise form: values at sample points inside and outside the
/// grid, away from the breakpoints, within giac's FFT noise (C0049).
#[test]
fn piecewise_forms() {
    let Some(giac) = giac() else {
        return;
    };
    let points = ["-1", "1/2", "1", "2", "3", "7/2", "5"];
    let mut wrong = Vec::new();
    for (input, var) in [
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,bins=4,piecewise)",
            "x",
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,range=[0,4],bins=4,piecewise)",
            "x",
        ),
        (
            "kernel_density([1,2,3,2],t=0..4,bandwidth=1/4,bins=4,piecewise)",
            "t",
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,bins=5,interp=1,output=piecewise)",
            "x",
        ),
    ] {
        let (mine, theirs) = sampled(giac, input, var, &points);
        if !within_noise(&mine, &theirs) {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `seq`'s last value after `srand(42)` (C0014).
fn seeded(giac: &Giac, input: &str) -> Option<Expr> {
    let answer = giac.eval(&format!("srand(42);{input}")).ok()?;
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
    ours(&answer[start..].replace("list[", "[")).ok()
}

/// A sample's density on 50 bins has 50 entries, giac seeded inside the
/// input (C0014); ours draws its own sample.
#[test]
fn seeded_sizes() {
    let Some(giac) = giac() else {
        return;
    };
    let theirs = seeded(
        giac,
        "X:=randvar(gammad,mean=5,stddev=2):;L:=kernel_density(sample(X,500),bins=50):;size(L)",
    );
    let mine =
        ours("size(kernel_density(sample(randvar(gammad,mean=5,stddev=2),500),bins=50))").ok();
    assert_eq!(theirs, Some(Expr::from(50)), "giac");
    assert_eq!(mine, Some(Expr::from(50)), "ours");
}

/// `text` read as [`ours`] reads it, and its value with the warnings of
/// its evaluation.
fn evaluated(text: &str) -> Option<(Expr, Expr, Vec<Warning>)> {
    let parsed = parse(text);
    if !parsed.errors().is_empty() {
        return None;
    }
    let lowered = lower(&parsed, Names::english());
    let [statement] = &lowered.statements[..] else {
        return None;
    };
    let e = decimal_balls(&statement.expr, cas_core::DEFAULT_WORKING_PRECISION);
    let mut ctx = Context::new();
    ctx.set_language(Some(Hook::english()));
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    Some((e, found.ok()?, report.warnings))
}

/// The message giac's error carries: the sentence before `. Error:`, or
/// giac's error itself.
fn message(error: &str) -> &str {
    let text = error.strip_prefix("GIAC_ERROR: ").unwrap_or(error);
    match text.split_once(". Error:") {
        Some((m, _)) => m,
        None => text.strip_prefix("Error: ").unwrap_or(text),
    }
}

/// Whether ours of `input` stays, its head the input's, with a warning of
/// `expected` (or none for `None`).
fn stays_with(input: &str, expected: Option<&str>) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && match expected {
                Some(text) => warnings
                    .iter()
                    .any(|w| matches!(w, Warning::Native { message, .. } if message == text)),
                None => warnings.is_empty(),
            }
    })
}

/// Giac's text, whitespace collapsed and a lone `|` dropped (C0077).
fn collapsed(s: &str) -> String {
    s.split_whitespace()
        .filter(|w| *w != "|")
        .collect::<Vec<_>>()
        .join(" ")
}

/// Each `(input, giac's answer)`: giac's text, whitespace collapsed, and
/// ours staying, with giac's message as its warning when `warned`.
fn check_stays(rows: &[(&str, &str)], warned: bool) {
    let mut wrong = Vec::new();
    for (input, answer) in rows {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let expected = warned.then(|| message(answer));
        if !stays_with(input, expected) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {expected:?}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's errors (.run/1.7/T5/extra.out) stay with its message; `bins=1`
/// with giac's FFT's "Invalid dimension".
#[test]
fn invalid_calls_stay() {
    const BAD: &str = "Error: Bad Argument Value";
    let rows = [
        (
            "kernel_density(5)",
            "GIAC_ERROR: Expected a list or a sequence of arguments. Error: Bad Argument Type"
                .to_owned(),
        ),
        (
            "kernel_density([1],exact)",
            format!("GIAC_ERROR: Data too small. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=-1)",
            format!("GIAC_ERROR: Invalid bandwidth specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=foo)",
            format!("GIAC_ERROR: Invalid bandwidth specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],bins=0)",
            format!("GIAC_ERROR: Number of bins must be positive. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],bins=2.5)",
            format!("GIAC_ERROR: Number of bins must be positive. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],range=[4,0])",
            format!("GIAC_ERROR: Invalid range specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],range=3)",
            format!("GIAC_ERROR: Invalid range specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],4..0)",
            format!("GIAC_ERROR: Invalid range specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],t=3)",
            format!("GIAC_ERROR: Invalid range specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],output=foo)",
            format!("GIAC_ERROR: Invalid method specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],interp=0)",
            format!("GIAC_ERROR: Invalid interpolation method specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],spline=0)",
            format!("GIAC_ERROR: Invalid interpolation method specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],3)",
            format!("GIAC_ERROR: Unrecognized option. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],\"a\")",
            format!("GIAC_ERROR: Unrecognized option. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],foo=[1,2])",
            format!("GIAC_ERROR: Invalid range specification. {BAD}"),
        ),
        (
            "kernel_density([1,2,a,2])",
            format!("GIAC_ERROR: Expected a real constant. {BAD}"),
        ),
        (
            "kernel_density([1,2,3,2],bandwidth=1/4,bins=1)",
            "GIAC_ERROR: Error: Invalid dimension".to_owned(),
        ),
    ];
    let rows: Vec<(&str, &str)> = rows.iter().map(|(i, a)| (*i, a.as_str())).collect();
    check_stays(&rows, true);
}

/// Departures: giac's `eval=v` (a value of the piecewise form, or its
/// errors) is not ported, the call stays without a warning.
#[test]
fn departures() {
    check_stays(
        &[
            (
                "kernel_density([1,2,3,2],bandwidth=1/4,bins=4,piecewise,eval=2)",
                "0.797899452895",
            ),
            (
                "kernel_density([1,2,3,2],bandwidth=1/4,bins=4,output=piecewise,eval=5)",
                "0",
            ),
            (
                "kernel_density([1,2,3,2],eval=[1])",
                "GIAC_ERROR: Expected a real constant. Error: Bad Argument Value",
            ),
            (
                "kernel_density([1,2,3,2],eval=2)",
                "GIAC_ERROR: Invalid method. Error: Bad Argument Value",
            ),
        ],
        false,
    );
}
