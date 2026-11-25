//! Xcas's distribution fitting against giac (CAS-629): `fitdistr` (giac's
//! answers in .run/1.7/giac/probe3.out and .run/1.7/T5/extra.out). giac's
//! mean divides by `N − 1` (optimization.cc:7481), ours by `N`, so the laws
//! whose answer uses the mean are departures: giac's text, whitespace
//! collapsed, and ours, exact on exact data (the same as ours of its form)
//! or within 10⁻¹² of the prototype's root (.run/1.7/T5/proto.py). The
//! Weibull, Cauchy and beta fits do not depend on the mean but through
//! their starting points: each parameter within 10⁻⁸ of giac's, which stops
//! its iteration at a step of 10⁻⁵. Invalid calls stay with giac's message
//! as the warning (C0022, C0077). Seeded rows (C0014) check both answers
//! against the law sampled, within four standard errors.

use cas_core::{Context, EvalOptions, Expr, View, Warning};
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

/// A law's parameters in `f64`.
fn parameters(e: &Expr) -> Option<Vec<f64>> {
    e.args().iter().map(float).collect()
}

/// Whether `mine` is `theirs`' law, each parameter within `tolerance` of
/// it, relatively.
fn near(mine: &Expr, theirs: &Expr, tolerance: f64) -> bool {
    mine.head() == theirs.head()
        && match (parameters(mine), parameters(theirs)) {
            (Some(m), Some(t)) => {
                m.len() == t.len()
                    && m.iter()
                        .zip(&t)
                        .all(|(a, b)| (a - b).abs() <= tolerance * b.abs())
            }
            _ => false,
        }
}

/// The Weibull, Cauchy and beta fits: each parameter within 10⁻⁸ of
/// giac's; ours certified balls.
#[test]
fn fits() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "fitdistr([1.2,2.3,0.7,3.1,1.9,2.2],weibull)",
        "fitdistr([1.2,2.3,0.7,3.1,1.9,2.2],cauchy)",
        "fitdistr([0.2,0.3,0.7,0.1,0.9,0.4],betad)",
        "fitdistr([1,2,3,4,5],cauchyd)",
        "fitdistr([0.2,0.5,0.6,0.9],Beta)",
        "fitdistr([0.5,1.5,2.0,3.5,4.0,7.0,9.5],cauchy)",
    ] {
        let theirs = giac.value(input);
        let mine = ours(input);
        let balls = mine
            .as_ref()
            .is_ok_and(|e| e.args().iter().all(|x| matches!(x.view(), View::Ball(_))));
        match (&mine, &theirs) {
            (Ok(m), Ok(t)) if balls && near(m, t, 1e-8) => {}
            _ => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whitespace collapsed and a lone `|` dropped (C0077).
fn collapsed(s: &str) -> String {
    s.split_whitespace()
        .filter(|w| *w != "|")
        .collect::<Vec<_>>()
        .join(" ")
}

/// giac prints `answer` for `input`, when giac is present.
fn giac_prints(input: &str, answer: &str) -> Option<String> {
    let printed = giac()?.eval(input).unwrap_or_default();
    (collapsed(&printed) != collapsed(answer))
        .then(|| format!("{input}: giac prints {printed}, recorded {answer}"))
}

/// Departures, the mean `Σx/N` (giac's `Σx/(N − 1)`): giac's answer and
/// ours, exact rows the same as ours of the form, decimal rows within
/// 10⁻¹² of the prototype's (proto.py, the roots at 40 digits).
#[test]
fn departures() {
    let rows = [
        (
            "fitdistr([1,2,3,4,5],normal)",
            "normald(3.75,1.78972763291)",
            "normald(3,sqrt(10)/2)",
        ),
        (
            "fitdistr([1,2,3,4,5,2,3,3],normal)",
            "normald(3.28571428571,1.32149751613)",
            "normald(23/8,sqrt(1218)/28)",
        ),
        (
            "fitdistr([1.2,2.3,0.7,3.1,1.9,2.2],normal)",
            "normald(2.28,0.947248647399)",
            "normald(1.9,0.850881895447306)",
        ),
        (
            "fitdistr([1,2,3,4,5,2,3,3],poisson)",
            "poisson(3.28571428571)",
            "poisson(23/8)",
        ),
        (
            "fitdistr([1,2,3,4,5,2,3,3],exponential)",
            "exponentiald(0.304347826087)",
            "exponentiald(8/23)",
        ),
        (
            "fitdistr([1,2,3,4,5,2,3,3],geometric)",
            "geometric(0.304347826087)",
            "geometric(8/23)",
        ),
        (
            "fitdistr([1.2,2.3,0.7,3.1,1.9,2.2],gammad)",
            "gammad(1.8889697753,0.828495515483)",
            "gammad(4.91602524944337,2.58738171023335)",
        ),
        (
            "fitdistr([1,2,3,4,5],weibulld)",
            "weibulld(2.29380667075,3.39430387666)",
            "weibulld(2.29380667071232,3.39429071810949)",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, form) in rows {
        wrong.extend(giac_prints(input, answer));
        let mine = ours(input);
        let expected = ours(form);
        let exact = !form.contains('.');
        let fine = match (&mine, &expected) {
            (Ok(m), Ok(f)) if exact => m == f,
            (Ok(m), Ok(f)) => near(m, f, 1e-12),
            _ => false,
        };
        if !fine {
            wrong.push(format!("{input}: ours {mine:?}, wanted {form}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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
/// the error itself after giac's ` | ` (C0007).
fn message(error: &str) -> &str {
    let text = error.strip_prefix("GIAC_ERROR: ").unwrap_or(error);
    match text.split_once('|') {
        Some((_, rest)) => rest.trim().strip_prefix("Error: ").unwrap_or(rest.trim()),
        None => text.split_once(". Error:").map_or(text, |(m, _)| m),
    }
}

/// Whether ours of `input` stays, its head the input's, with a warning of
/// `expected`.
fn stays_with(input: &str, expected: &str) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && warnings
                .iter()
                .any(|w| matches!(w, Warning::Native { message, .. } if message == expected))
    })
}

/// giac's errors (probe3.out, extra.out) stay with its message, the four
/// laws giac does not fit among them.
#[test]
fn invalid_calls_stay() {
    const BAD: &str = "Error: Bad Argument Value";
    let unrecognized = format!("GIAC_ERROR: Unrecognized probability distribution. {BAD}");
    let rows = [
        ("fitdistr([1,2,3,4,5,2,3,3],binomial)", unrecognized.clone()),
        ("fitdistr([1,2,3,4,5,2,3,3],uniform)", unrecognized.clone()),
        (
            "fitdistr([1.2,2.3,0.7,3.1,1.9,2.2],lognormal)",
            unrecognized.clone(),
        ),
        ("fitdistr([1,2,3],negbinomial)", unrecognized.clone()),
        ("fitdistr([1,2,3],foo)", unrecognized),
        (
            "fitdistr([1,2,3])",
            "GIAC_ERROR: fitdistr([1,2,3])  |  Error: Bad Argument Type".to_owned(),
        ),
        (
            "fitdistr([1,2,3],normal,1)",
            "GIAC_ERROR: Expected a sequence of two arguments. Error: Invalid dimension".to_owned(),
        ),
        (
            "fitdistr(3,normal)",
            format!("GIAC_ERROR: Expected a list. {BAD}"),
        ),
        (
            "fitdistr([1],normal)",
            format!("GIAC_ERROR: At least two samples are required. {BAD}"),
        ),
        (
            "fitdistr([1,a],normal)",
            format!("GIAC_ERROR: Expected a real constant. {BAD}"),
        ),
        (
            "fitdistr([1,2.5],poisson)",
            format!("GIAC_ERROR: Expected a nonnegative integer. {BAD}"),
        ),
        (
            "fitdistr([1,-2],poisson)",
            format!("GIAC_ERROR: Expected a nonnegative integer. {BAD}"),
        ),
        (
            "fitdistr([1,-2],exponential)",
            format!("GIAC_ERROR: Expected a positive number. {BAD}"),
        ),
        (
            "fitdistr([0,2],geometric)",
            format!("GIAC_ERROR: Expected a positive integer. {BAD}"),
        ),
        (
            "fitdistr([1,-2],gammad)",
            format!("GIAC_ERROR: Expected a positive number. {BAD}"),
        ),
        (
            "fitdistr([0.5,1.5],betad)",
            format!("GIAC_ERROR: Expected a number in [0,1]. {BAD}"),
        ),
        (
            "fitdistr([0,0],weibull)",
            format!("GIAC_ERROR: Expected a nonzero variance. {BAD}"),
        ),
        (
            "fitdistr([1,-2],weibull)",
            format!("GIAC_ERROR: Expected a positive number. {BAD}"),
        ),
    ];
    let mut wrong = Vec::new();
    for (input, error) in &rows {
        wrong.extend(giac_prints(input, error));
        if !stays_with(input, message(error)) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {:?}",
                evaluated(input),
                message(error)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures where giac answers and ours stays with a warning: Weibull
/// of equal data (giac's variance about its mean `Σx/(N − 1)` is not 0,
/// and its iteration fails, `undef`), and beta with a 0 among the data
/// (giac answers its start, ln 0 being −∞).
#[test]
fn departures_that_stay() {
    let mut wrong = Vec::new();
    for (input, answer, warning) in [
        (
            "fitdistr([2,2,2],weibull)",
            "undef",
            "Expected a nonzero variance",
        ),
        (
            "fitdistr([0,0.5],betad)",
            "betad(-0.0,-0.0)",
            "Unable to compute parameters for Beta distribution",
        ),
    ] {
        wrong.extend(giac_prints(input, answer));
        if !stays_with(input, warning) {
            wrong.push(format!("{input}: ours {:?}", evaluated(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The last value of giac's answer to `srand(42);input` (C0014).
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
    ours(&answer[start..]).ok()
}

/// Ours of `input` seeded with 42 in the same input.
fn ours_seeded(input: &str) -> Option<Expr> {
    ours(&format!("[srand(42),{input}]"))
        .ok()
        .and_then(|e| e.args().get(1).cloned())
}

/// Whether `e` is `head`'s law, each parameter within its bound of the
/// sampled law's.
fn within(e: Option<&Expr>, head: &str, law: [f64; 2], bounds: [f64; 2]) -> bool {
    let Some(e) = e else {
        return false;
    };
    let Ok(expected) = ours(&format!("{head}(1,1)")) else {
        return false;
    };
    e.head() == expected.head()
        && parameters(e)
            .is_some_and(|p| p.len() == 2 && (0..2).all(|i| (p[i] - law[i]).abs() <= bounds[i]))
}

/// Seeded samples of 1000 (probe3's last rows): both answers within four
/// standard errors of the law sampled. Weibull `(1/2, 1)`: `sd(k) ≈ 0.78
/// k/√n`, `sd(λ) ≈ 1.05 λ/(k √n)`, bounds 0.05 and 0.27 (giac `(0.4968,
/// 0.9409)`); the logarithms of a log-normal sample, normal `(5, √2)`:
/// `sd(m) = σ/√n`, `sd(s) ≈ σ/√(2n)`, bounds 0.18 and 0.127 (giac `(5.0334,
/// 1.4515)`).
#[test]
fn seeded_rows() {
    let Some(giac) = giac() else {
        return;
    };
    let weibull = "fitdistr(randvector(1000,weibulld,1/2,1),weibull)";
    let (law, bounds) = ([0.5, 1.0], [0.05, 0.27]);
    let theirs = seeded(giac, weibull);
    let mine = ours_seeded(weibull);
    assert!(
        within(theirs.as_ref(), "weibulld", law, bounds),
        "giac's {theirs:?}"
    );
    assert!(
        within(mine.as_ref(), "weibulld", law, bounds),
        "ours {mine:?}"
    );
    let theirs = seeded(
        giac,
        "X:=randvar(normal,mean=5,variance=2):;S:=sample(exp(X),1000):;fitdistr(log(S),normal)",
    );
    let mine =
        ours_seeded("fitdistr(log(sample(exp(randvar(normal,mean=5,variance=2)),1000)),normal)");
    let (law, bounds) = ([5.0, 2f64.sqrt()], [0.18, 0.127]);
    assert!(
        within(theirs.as_ref(), "normald", law, bounds),
        "giac's {theirs:?}"
    );
    assert!(
        within(mine.as_ref(), "normald", law, bounds),
        "ours {mine:?}"
    );
}
