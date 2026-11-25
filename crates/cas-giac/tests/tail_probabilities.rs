//! Xcas's upper tails `UTPN`, `UTPT`, `UTPC` and `UTPF` (CAS-625), and
//! `cdf`, `icdf` and `mgf` of a law by its name, of a random variable's
//! object or of a list (CAS-626), and `logistic` and `logistic_regression`
//! (CAS-627, the regression's six entries one by one), against giac:
//! giac's answers read back are this CAS's value of the same input, the
//! difference expanded 0 or,
//! with symbols, the two equal at a rational point within 10⁻²⁰; giac's
//! decimals within its 12 digits, ours exact where the input is (C0006,
//! C0032, R67); invalid calls stay with giac's message as the warning
//! (libgiac's text, C0007, C0077). Departures: giac's answer as caseval
//! prints it, whitespace collapsed, and ours.

use cas_core::{Context, EvalOptions, Expr, Symbol, Warning};
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

/// Whether `e` evaluates to `True`.
fn holds(ctx: &Context, e: &Expr) -> bool {
    ctx.eval(e).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The rational point the symbols of the rows take: `t` inside every
/// law's domain of `mgf`.
fn point() -> Vec<(Expr, Expr)> {
    let q = |n: i64, d: i64| Expr::from(n) / Expr::from(d);
    [
        ("x", q(7, 10)),
        ("t", q(1, 5)),
        ("a", q(3, 2)),
        ("b", q(2, 1)),
        ("k", q(3, 1)),
        ("l", q(2, 1)),
        ("lambda", q(3, 2)),
        ("mu", q(1, 3)),
        ("n", q(4, 1)),
        ("p", q(1, 3)),
        ("sigma", q(2, 3)),
    ]
    .into_iter()
    .map(|(name, value)| (Expr::symbol(name), value))
    .collect()
}

/// Whether `a` and `b` agree within a relative `10^-digits`, at 30 digits.
fn close(ctx: &Context, a: &Expr, b: &Expr, digits: i64) -> bool {
    let numeric = |e: &Expr| Expr::call(Symbol::N, [e.clone(), Expr::from(30)]);
    let gap = Expr::call(Symbol::ABS, [numeric(a) - numeric(b)]);
    let scale = Expr::call(
        Symbol::MAX,
        [Expr::from(1), Expr::call(Symbol::ABS, [numeric(b)])],
    );
    holds(
        ctx,
        &Expr::call(Symbol::LESS, [gap, scale * Expr::from(10).pow(-digits)]),
    )
}

/// Whether two answers are one value: the same expression (`undef`, a
/// call that stays), their difference expanded 0, or the two within 10⁻²⁰
/// at [`point`] (`erf(x/sqrt(2))` against `erf(x*sqrt(2)/2)`).
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    if ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0)) {
        return true;
    }
    let at = point();
    close(ctx, &a.substitute(&at), &b.substitute(&at), 20)
}

/// Each input: giac's answer, read back, is this CAS's value.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if same_value(&ctx, &mine, &theirs) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input whose giac answer has decimals: ours within a relative
/// 10⁻⁹ of it, giac's 12 digits (C0006), at [`point`] if it has symbols.
fn check_near(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let at = point();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine))
                if close(&ctx, &mine.substitute(&at), &theirs.substitute(&at), 9) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
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

/// The message giac's error carries, after giac's ` | ` (C0007).
fn message(error: &str) -> &str {
    let text = error.strip_prefix("GIAC_ERROR: ").unwrap_or(error);
    match text.split_once('|') {
        Some((_, rest)) => rest.trim().strip_prefix("Error: ").unwrap_or(rest.trim()),
        None => text,
    }
}

/// Giac's text, whitespace collapsed and a lone `|` dropped: caseval.py
/// prints libgiac's line break as ` | ` (C0077).
fn collapsed(s: &str) -> String {
    s.split_whitespace()
        .filter(|w| *w != "|")
        .collect::<Vec<_>>()
        .join(" ")
}

/// Each `(input, giac's error)`: giac's text, and ours staying with giac's
/// message as its warning.
fn check_stays(rows: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, error) in rows {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
        let stays = evaluated(input).is_some_and(|(e, found, warnings)| {
            found.head() == e.head()
                && warnings
                    .iter()
                    .any(|w| matches!(w, Warning::Native { message: m, .. } if m == message(error)))
        });
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {:?}",
                evaluated(input),
                message(error)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The upper tails, exact and symbolic where giac's are; the calls giac
/// leaves, which ours leaves too (probe1, probe7, probe13, T3 probe).
#[test]
fn upper_tails() {
    check(&[
        "UTPN(0,1,0)",
        "UTPN(0,1,x)",
        "UTPN(x)",
        "UTPN(1/2)",
        "UTPN(0,4,2)",
        "UTPC(3,0)",
        "UTPC(-1,0)",
        "UTPC(2,-3)",
        "UTPC(0.5,2)",
        "UTPC(0.5,-1)",
        "UTPF(2,4,0)",
        "UTPF(2,2,-1)",
        "UTPF(2.5,2,1)",
    ]);
}

/// The upper tails where giac's answers are decimals, `UTPC(2,6)` among
/// them, ours `exp(-3)` (probe1, probe7, probe13, T3 probe).
#[test]
fn upper_tails_near() {
    check_near(&[
        "UTPN(1.96)",
        "UTPN(1,4,2.96)",
        "UTPN(1,4,2.96)+normal_cdf(1,2,2.96)",
        "UTPT(3,2.35)",
        "UTPT(3,-2.35)",
        "UTPT(3,0)",
        "UTPT(2,1)",
        "UTPT(1,1)",
        "UTPT(1,0.5)",
        "UTPT(3.0,2.35)",
        "UTPC(2,6.1)",
        "UTPC(4,6.1)",
        "UTPC(1,1/2)",
        "UTPC(2,6)",
        "UTPC(2,2)",
        "UTPC(5,3.0)",
        "UTPC(2.0,6)",
        "UTPF(4,10,3.5)",
        "UTPF(4,4,2.1)",
        "UTPF(2,2,1)",
        "UTPF(2,2,2)",
        "UTPF(1,3,2.0)",
        "UTPF(4.0,10,3.5)",
    ]);
}

/// `cdf`, `icdf` and `mgf` by a law's name (giac's `exp` among them) and on
/// `randvar`'s objects, exact or symbolic where giac's are; `mgf`'s
/// `undef` for the laws giac has no form for (probe1, probe7, probe13).
#[test]
fn laws() {
    check(&[
        "cdf(normald,0,1,x)",
        "cdf(normald,0,1)",
        "cdf(binomial,10,1/2,3)",
        "cdf(uniform,0,1,1/2)",
        "cdf(exponential,2,1)",
        "cdf(exp,2,1)",
        "cdf(geometric,1/4,3)",
        "cdf(weibull,2,3,1)",
        "cdf(cauchy,0,1,1)",
        "cdf(randvar(binomial,10,1/2),3)",
        "icdf(uniform,0,1,1/4)",
        "icdf(exponential,2,1/2)",
        "icdf(exp,2,1/2)",
        "icdf([1,2,3,4],0)",
        "icdf([1,2,3,4],1)",
        "icdf([1,2,3,4],2)",
        "icdf([1,2,3,4],-0.5)",
        "icdf([1,2,3,4],x)",
        "mgf(normald,1,0)",
        "mgf(normald,0,1)",
        "mgf(normald,mu,sigma)",
        "mgf(normald(0,1))",
        "mgf(poisson,5)",
        "mgf(poisson,lambda)",
        "mgf(randvar(poisson,3))",
        "mgf(binomial,n,p)",
        "mgf(negbinomial,3,1/2)",
        "mgf(gammad,a,b)",
        "mgf(chisquare,k)",
        "mgf(geometric,p)",
        "mgf(uniformd,a,b)",
        "mgf(uniform,0,1)",
        "mgf(exponentiald,l)",
        "mgf(exponential,2)",
        "mgf(exp,2)",
        "mgf(betad,2,3)",
        "mgf(student,3)",
        "mgf(fisher,2,3)",
        "mgf(weibull,2,3)",
        "mgf(cauchy,0,1)",
    ]);
}

/// `cdf` and `icdf` where giac's answers are decimals; `icdf` of a list,
/// ours its exact element (probe7, probe13, T3 probe).
#[test]
fn laws_near() {
    check_near(&[
        "cdf(normald,0,1,1.96)",
        "cdf(normald(0,1),1.96)",
        "cdf(poisson,3,2)",
        "cdf(chisquare,2,6.1)",
        "cdf(student,3,2.35)",
        "cdf(fisher,4,10,3.5)",
        "cdf(randvar(normal,mean=0,stddev=1),1.96)",
        "cdf(randvar(poisson,3),2)",
        "icdf(normald,0,1,0.975)",
        "icdf(normald(0,1),0.975)",
        "icdf(binomial,10,1/2,0.3)",
        "icdf(poisson,3,0.5)",
        "icdf(chisquare,2,0.95)",
        "icdf(randvar(normal,mean=0,stddev=1),0.975)",
        "icdf(randvar(poisson,3),0.5)",
        "icdf([1,2,3,4],0.25)",
        "icdf([1,2,3,4],1/2)",
        "icdf([1,2,3,4],0.6)",
        "icdf([3,1,2,2],0.5)",
        "mgf(randvar(normal,mean=0,stddev=1))",
    ]);
}

/// giac's errors: the call stays, giac's message the warning.
#[test]
fn invalid_calls_stay() {
    const BAD: &str = "Error: Bad Argument Value";
    const DIMENSION: &str = "Error: Invalid dimension";
    let rows: Vec<(String, String)> = [
        ("UTPN(1,0,2)", "UTPN(1,0,2)", BAD),
        ("UTPN(1,2)", "UTPN(1,2)", BAD),
        ("UTPT(3)", "UTPT(3)", BAD),
        ("UTPT(1.5,1)", "UTPT(1.5,1)", BAD),
        ("UTPT(0,1)", "UTPT(0,1)", DIMENSION),
        ("UTPT(-2,1)", "UTPT(-2,1)", DIMENSION),
        ("UTPC(-1,2)", "UTPC(-1,2)", BAD),
        ("UTPF(0,2,1)", "UTPF(0,2,1)", DIMENSION),
        ("UTPF(400,2,1)", "UTPF(400,2,1)", DIMENSION),
        ("cdf(normal,0,1,1)", "cdf('normal',0,1,1)", BAD),
        ("cdf(normald,0,1,1,2)", "cdf('normald',0,1,1,2)", BAD),
        ("cdf([1,2,3,4],x)", "cdf([1,2,3,4],x)", BAD),
        ("cdf(1,2)", "cdf(1,2)", BAD),
        ("cdf(x)", "cdf(x)", BAD),
        ("icdf([1,2,3,4])", "icdf([1,2,3,4])", BAD),
        (
            "icdf(randvar([1,2,3]),0.5)",
            "icdf(randvar([1,2,3]),0.5)",
            BAD,
        ),
        ("mgf(randvar([1,2,3]))", "mgf(randvar([1,2,3]))", BAD),
        ("mgf(normald)", "mgf('normald')", BAD),
        ("mgf(normald,0)", "mgf('normald',0)", BAD),
        ("mgf(foo,1)", "mgf(foo,1)", BAD),
    ]
    .into_iter()
    .map(|(input, call, error)| (input.to_owned(), format!("GIAC_ERROR: {call} | {error}")))
    .collect();
    let rows: Vec<(&str, &str)> = rows.iter().map(|(a, b)| (a.as_str(), b.as_str())).collect();
    check_stays(&rows);
}

/// giac's answers where ours differs, each with its reason:
/// - a list's `cdf`: giac's dichotomy answers an index, `(#{Lᵢ ≤ x} −
///   1)/n` (`−1/n` below the least value) and steps `[x₍ᵢ₎, i/n]`; ours the
///   empirical distribution function, `#{Lᵢ ≤ x}/n`, steps `(i + 1)/n`;
/// - `icdf(normald,0,1,1/2)`: giac's rounding, ours exactly 0;
/// - the gamma and beta laws' `cdf`: giac keeps `igamma` and `Beta`, ours
///   `gammad_cdf` and `betad_cdf` evaluate;
/// - the negative binomial: giac has no `cdf` or `icdf` (`undef(…)`), ours
///   `negbinomial_cdf` and `negbinomial_icdf`;
/// - `UTPT(3,x)` and `UTPC(2,x)`: giac stays, ours `1 − F` expanded.
#[test]
fn departures() {
    let cases = [
        ("cdf([1,2,3,4],2)", "1/4", "1/2"),
        ("cdf([1,2,3,4],1)", "0", "1/4"),
        ("cdf([1,2,3,4],2.5)", "1/4", "1/2"),
        ("cdf([1,2,3,4],4)", "3/4", "1"),
        ("cdf([1,2,3,4],5)", "3/4", "1"),
        ("cdf([1,2,3,4],0)", "-1/4", "0"),
        ("cdf([3,1,2,2],2)", "1/2", "3/4"),
        (
            "cdf([1,2,3,4])",
            "[[1,0],[2,1/4],[3,1/2],[4,3/4]]",
            "[[1,1/4],[2,1/2],[3,3/4],[4,1]]",
        ),
        (
            "cdf([3,1,2,2])",
            "[[1,0],[2,1/4],[2,1/2],[3,3/4]]",
            "[[1,1/4],[2,1/2],[2,3/4],[3,1]]",
        ),
        ("icdf(normald,0,1,1/2)", "-1.88884303598e-19", "0"),
        ("cdf(gammad,2,3,1)", "igamma(2,3,1)", "gammad_cdf(2,3,1)"),
        (
            "cdf(betad,2,3,1/2)",
            "12*Beta(2,3,1/2)",
            "betad_cdf(2,3,1/2)",
        ),
        (
            "cdf(negbinomial,3,1/2,2)",
            "undef(3,1/2,2)",
            "negbinomial_cdf(3,1/2,2)",
        ),
        (
            "icdf(negbinomial,3,1/2,0.5)",
            "undef(3,1/2,0.5)",
            "negbinomial_icdf(3,1/2,0.5)",
        ),
        ("UTPT(3,x)", "UTPT(3,x)", "expand(1-student_cdf(3,x))"),
        ("UTPC(2,x)", "UTPC(2,x)", "expand(1-chisquare_cdf(2,x))"),
    ];
    let mut wrong = Vec::new();
    for (input, answer, mine) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours(input);
        if found.is_err() || found.as_ref().ok() != ours(mine).ok().as_ref() {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The logistic (CAS-627): giac's exact and symbolic answers, its decimals
/// within its 12 digits, and its Invalid dimension (T3 rows).
#[test]
fn logistic() {
    check(&[
        "logistic(0)",
        "logistic(1,0,2,3)",
        "logistic(x)",
        "logistic(x,1,2,3)",
    ]);
    check_near(&["logistic(0.5)", "logistic(0.0,0.5,10,2)", "logistic(1.8)"]);
    check_stays(&[
        (
            "logistic()",
            "GIAC_ERROR: logistic(seq[]) | Error: Invalid dimension",
        ),
        (
            "logistic(1,2,3,4,5)",
            "GIAC_ERROR: logistic(1,2,3,4,5) | Error: Invalid dimension",
        ),
    ]);
}

/// Each input: giac's six entries and ours, entry by entry: a rational
/// entry the same value (C0032), any other within giac's 12 digits at
/// several values of `x` (C0006).
fn check_entries(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let q = |n: i64, d: i64| Expr::from(n) / Expr::from(d);
    let points: Vec<Vec<(Expr, Expr)>> = [q(0, 1), q(7, 10), q(5, 1), q(25, 2)]
        .into_iter()
        .map(|v| vec![(Expr::symbol("x"), v)])
        .collect();
    let list = Expr::from(Symbol::LIST);
    let mut wrong = Vec::new();
    for input in inputs {
        let agree = match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) => {
                theirs.head() == Some(&list)
                    && mine.head() == Some(&list)
                    && theirs.args().len() == 6
                    && mine.args().len() == 6
                    && theirs.args().iter().zip(mine.args()).all(|(t, m)| {
                        if t.to_rational().is_some() {
                            same_value(&ctx, m, t)
                        } else {
                            points
                                .iter()
                                .all(|at| close(&ctx, &m.substitute(at), &t.substitute(at), 9))
                        }
                    })
            }
            _ => false,
        };
        if !agree {
            wrong.push(format!(
                "{input}: ours {:?}, giac's {:?}",
                ours(input),
                giac.value(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's logistic regression (CAS-627; probe7.out, probe10.out, T3 rows):
/// exact total and peak where giac's are, the fitted curves within its
/// digits; constant quotients, giac's `[undef,undef,infinity,infinity,
/// undef,undef]`; a time that is no number stays with its message.
#[test]
fn logistic_regressions() {
    check_entries(&[
        "logistic_regression([1,2,4,7,11],0,1)",
        "logistic_regression([11,7,4,2,1],0,1)",
        "logistic_regression([0.5,1,2,3.5,4.5,4],1,1)",
        "logistic_regression([1,2,4,7,11],0,1.0)",
        "logistic_regression([1.0,2,4,7,11],0,1)",
    ]);
    check(&[
        "logistic_regression([1,2,4],0,1)",
        "logistic_regression([1,2],0,1)",
    ]);
    check_stays(&[(
        "logistic_regression([1,2,4,7,11],x,1)",
        "GIAC_ERROR: logistic_regression([1,2,4,7,11],x,1) | Error: Bad Argument Value",
    )]);
}
