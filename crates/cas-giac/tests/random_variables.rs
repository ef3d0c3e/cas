//! Xcas's random variables against giac (CAS-624): `randvar` and the
//! objects it answers, `sample` and `randvector` drawing from them, and
//! `randmultinomial`. `randvar` draws nothing, so its objects compare with
//! giac's, the parameters exact (`normal(ours − giac)`, C0032) or within
//! giac's error where giac prints decimals (C0006, C0076); invalid calls stay,
//! giac's message the warning (libgiac's text, C0007, C0018). giac's generator
//! is not ours, so a draw is a property both answers have, each side seeded in
//! the same input (C0014), never giac's digits. Departures: giac's answer as
//! caseval prints it, whitespace collapsed, and ours.

use cas_core::{Context, EvalOptions, Expr, Integer, Rational, Symbol, View, Warning};
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

/// Whether two parameters agree: `normal(a − b)` is 0, or the two within a
/// relative 10⁻⁹ of each other, giac's 12 digits (C0006).
fn same_parameter(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]);
    if ctx.eval(&difference).ok() == Some(Expr::from(0)) {
        return true;
    }
    let numeric = |e: &Expr| Expr::call(Symbol::N, [e.clone(), Expr::from(30)]);
    let gap = Expr::call(Symbol::ABS, [numeric(a) - numeric(b)]);
    let scale = Expr::call(
        Symbol::MAX,
        [Expr::from(1), Expr::call(Symbol::ABS, [numeric(a)])],
    );
    holds(
        ctx,
        &Expr::call(Symbol::LESS, [gap, scale * Expr::from(10).pow(-9)]),
    )
}

/// Each input: giac's law, read back, is ours, the same head and each
/// parameter the same ([`same_parameter`]).
fn check_laws(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine))
                if theirs.head() == mine.head()
                    && theirs.args().len() == mine.args().len()
                    && theirs
                        .args()
                        .iter()
                        .zip(mine.args())
                        .all(|(a, b)| same_parameter(&ctx, a, b)) => {}
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
/// `expected` (C0022, C0056).
fn stays_with(input: &str, expected: &str) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && warnings
                .iter()
                .any(|w| matches!(w, Warning::Native { message, .. } if message == expected))
    })
}

/// Giac's text, whitespace collapsed and a lone `|` dropped: the recorded
/// texts come from scripts/caseval.py, which prints giac's line break as
/// ` | `, while libgiac's own text keeps the break (CI run 45565458930).
fn collapsed(s: &str) -> String {
    s.split_whitespace()
        .filter(|w| *w != "|")
        .collect::<Vec<_>>()
        .join(" ")
}

/// Each `(input, giac's error)`: giac's text, whitespace collapsed, and ours
/// staying with giac's message as its warning.
fn check_stays(rows: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, error) in rows {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
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

/// Objects (CAS-624; giac's answers in crates/cas-core/proto/randvar/probe.cas,
/// crates/cas-core/proto/giac/probe2.cas): parameters given, from `mean=`, `stddev=`
/// and `variance=` (exact where giac's are: `gammad(36,3)`,
/// `uniformd(10-2*sqrt(3),10+2*sqrt(3))`, `geometric(1/4)` from giac's
/// unsimplified radical), giac's decimals for decimal input
/// (`normald(5.0,1.5)`) within its error.
#[test]
fn objects() {
    check_laws(&[
        "randvar(fisher,2,3)",
        "randvar(gammad,mean=12,variance=4)",
        "randvar(uniform,mean=10,stddev=2)",
        "randvar(normal,mean=5,stddev=1.5)",
        "randvar(normal,mean=5,variance=2,stddev=1)",
        "randvar(exp,2)",
        "randvar(poisson,3)",
        "randvar(geometric,1/4)",
        "randvar(cauchy,0,1)",
        "randvar(weibull,2,3)",
        "randvar(normal,5,2)",
        "randvar(betad,mean=2/5,variance=1/25)",
        "randvar(negbinomial,mean=6,stddev=3)",
        "randvar(geometric,variance=12)",
        "randvar(student,3)",
        "randvar(chisquare,3)",
    ]);
}

/// giac's checks (graphtheory.cc, moyal.cc): the call stays, giac's message
/// the warning; `randmultinomial` with weights not summing to 1, or values
/// not matching them, is giac's Bad Argument Value.
#[test]
fn invalid_calls_stay() {
    const BAD: &str = "Error: Bad Argument Value";
    check_stays(&[
        (
            "randvar(normal,mean=5,stddev=-2)",
            "GIAC_ERROR: Standard deviation must be positive. Error: Bad Argument Value",
        ),
        (
            "randvar(fisher,2)",
            &format!("GIAC_ERROR: randvar('fisher',2) | {BAD}"),
        ),
        (
            "randvar(normal,mean=5,variance=-1)",
            "GIAC_ERROR: Variance must be nonnegative. Error: Bad Argument Value",
        ),
        (
            "randvar(normal,mean=5,range=x)",
            "GIAC_ERROR: Invalid range specification. Error: Bad Argument Value",
        ),
        (
            "randvar(normal,mean=5,foo=2)",
            "GIAC_ERROR: Unknown distribution parameter. Error: Bad Argument Value",
        ),
        (
            "randvar(normal,mean=5,\"x\")",
            "GIAC_ERROR: Invalid distribution specification. Error: Bad Argument Value",
        ),
        (
            "randvar(fisher,mean=2,stddev=1)",
            "GIAC_ERROR: Specifying moments is not supported for this distribution. \
             Error: Bad Argument Value",
        ),
        (
            "randvar([1,2],[a,b,c])",
            "GIAC_ERROR: Items and weights sizes do not match. Error: Invalid dimension",
        ),
        (
            "randvar([1,-2,3])",
            &format!("GIAC_ERROR: randvar([1,-2,3]) | {BAD}"),
        ),
        (
            "randvar([[1,2,3],[4,5,6]])",
            "GIAC_ERROR: Expected a two-column matrix. Error: Invalid dimension",
        ),
        (
            "randvar(uniform,3..1)",
            "GIAC_ERROR: Invalid distribution parameter(s). Error: Bad Argument Value",
        ),
        (
            "randmultinomial([1/2,1/2,1/2])",
            &format!("GIAC_ERROR: randmultinomial([1/2,1/2,1/2]) | {BAD}"),
        ),
        (
            "randmultinomial([1/2,1/3])",
            &format!("GIAC_ERROR: randmultinomial([1/2,1/3]) | {BAD}"),
        ),
        (
            "randmultinomial([1/2,1/2],[a])",
            &format!("GIAC_ERROR: randmultinomial([1/2,1/2],[a]) | {BAD}"),
        ),
    ]);
}

/// Ours only: more draws than the limits allow stays (C0019).
#[test]
fn too_many_draws_stay() {
    let input = "sample(randvar(poisson,3),10^7)";
    assert!(
        stays_with(input, "too many draws"),
        "{input}: ours {:?}",
        evaluated(input)
    );
}

/// giac's objects where ours differs, each with its reason:
/// - `binomial(n,p)` from `mean=` and `stddev=`: giac's law prints as the
///   coefficient's name; our `binomial` is the coefficient, so the law is
///   `randvar(binomial,n,p)`;
/// - giac's `discreted(…)` alias table (Walker's, decimals) for weights with
///   values; ours `randvar([w…],[v…])`, exact, `[[v,w],…]` normalized to it;
/// - `string(randvar(normal))`: giac's decimals, ours exact, `normald(0,1)`;
/// - `uniform` over a range: giac's `'uniformd',a,b` sequence in decimals,
///   ours `uniformd(a,b)` exact;
/// - `randvar(multinomial,P)`: giac's sequence, which does not sample
///   (`sample` of it is Bad Argument Value); ours the weights form, which
///   draws indices ([`departures_that_draw`]).
#[test]
fn departures() {
    let cases = [
        (
            "randvar(binomial,mean=18,stddev=4)",
            "binomial(162,1/9)",
            "randvar(binomial,162,1/9)",
        ),
        (
            "randvar(binomial,mean=10,stddev=2)",
            "binomial(17,10/17)",
            "randvar(binomial,17,10/17)",
        ),
        (
            "randvar([[\"apple\",1/3],[\"orange\",1/4],[\"pear\",1/5],[\"plum\",13/60]])",
            "discreted(4,1.0,1,0.8,0.866666666667,1,0,0,0,\"apple\",\"orange\",\"pear\",\"plum\")",
            "randvar([1/3,1/4,1/5,13/60],[\"apple\",\"orange\",\"pear\",\"plum\"])",
        ),
        (
            "randvar([3,1,2,5],[alpha,beta,gamma,delta])",
            "discreted(4,0.454545454545,0.363636363636,0.727272727273,1,3,0,3,0,\
             alpha,beta,gamma,delta)",
            "randvar([3,1,2,5],[alpha,beta,gamma,delta])",
        ),
        (
            "string(randvar(normal))",
            "\"normald(0.0,1.0)\"",
            "string(normald(0,1))",
        ),
        (
            "randvar(uniform,e..pi)",
            "'uniformd',2.71828182846,3.14159265359",
            "uniformd(e,pi)",
        ),
        (
            "randvar(uniform,range=1..3)",
            "'uniformd',1.0,3.0",
            "uniformd(1,3)",
        ),
        (
            "randvar(multinomial,[1/2,1/2])",
            "'multinomial',[1/2,1/2]",
            "randvar([1/2,1/2])",
        ),
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

/// A number as a rational: itself, or a ball's midpoint.
fn value(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Ball(b) => b.mid().checked_to_rational(),
        _ => e.to_rational(),
    }
}

fn q(n: i64) -> Rational {
    Rational::from(Integer::from(n))
}

/// An integer in `[low, high]`.
fn integer_in(low: i64, high: i64) -> impl Fn(&Expr) -> bool {
    move |e| value(e).is_some_and(|v| v.is_integer() && q(low) <= v && v <= q(high))
}

/// A non-negative integer.
fn natural(e: &Expr) -> bool {
    value(e).is_some_and(|v| v.is_integer() && q(0) <= v)
}

/// A decimal.
fn decimal(e: &Expr) -> bool {
    matches!(e.view(), View::Ball(_))
}

/// A real decimal.
fn real(e: &Expr) -> bool {
    decimal(e) && value(e).is_some()
}

/// A positive decimal.
fn positive(e: &Expr) -> bool {
    decimal(e) && value(e).is_some_and(|v| q(0) < v)
}

/// A decimal strictly between 0 and 1.
fn unit_open(e: &Expr) -> bool {
    decimal(e) && value(e).is_some_and(|v| q(0) < v && v < q(1))
}

/// A list of `n` items that each pass `item`.
fn list_of(n: usize, item: impl Fn(&Expr) -> bool) -> impl Fn(&Expr) -> bool {
    move |e| {
        e.head() == Some(&Expr::from(Symbol::LIST))
            && e.args().len() == n
            && e.args().iter().all(&item)
    }
}

/// One of the items of `expected`.
fn one_of(expected: &str) -> impl Fn(&Expr) -> bool {
    let items = ours(expected).expect("readable").args().to_vec();
    move |e| items.contains(e)
}

/// A property of an answer.
type Property<'p> = &'p dyn Fn(&Expr) -> bool;

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

/// Ours of `input` seeded with 42 in the same input.
fn ours_seeded(input: &str) -> Option<Expr> {
    ours(&format!("[srand(42),{input}]"))
        .ok()
        .and_then(|e| e.args().get(1).cloned())
}

/// Each input seeded with 42 on both sides (C0014), giac's answer and ours
/// each with its property.
fn check(cases: &[(&str, Property<'_>)]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, property) in cases {
        let theirs = giac
            .eval(&format!("srand(42);{input}"))
            .map_err(|e| e.to_string())
            .and_then(|t| ours(last(&t)).map_err(|e| e.to_string()));
        let mine = ours_seeded(input);
        match (&theirs, &mine) {
            (Ok(t), Some(m)) if property(t) && property(m) => {}
            _ => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Draws from the objects (giac's in crates/cas-core/proto/giac/probe2.cas, probe14.cas,
/// crates/cas-core/proto/randvar/c2): supports and sizes; giac's weights form draws indices
/// from 0, and its parameter-only density call is itself an object.
#[test]
fn draws() {
    let fruit = "[\"apple\",\"orange\",\"pear\",\"plum\"]";
    check(&[
        ("sample(randvar(fisher,2,3),5)", &list_of(5, positive)),
        (
            "sample(randvar([[\"apple\",1/3],[\"orange\",1/4],[\"pear\",1/5],[\"plum\",13/60]]),6)",
            &list_of(6, one_of(fruit)),
        ),
        (
            "sample(randvar(k->1-(k/10)^2,range=-10..10),5)",
            &list_of(5, integer_in(-10, 10)),
        ),
        (
            "sample(randvar([3,1,2,5],[alpha,beta,gamma,delta]),8)",
            &list_of(8, one_of("[alpha,beta,gamma,delta]")),
        ),
        (
            "sample(randvar(multinomial,[1/2,1/3,1/6],[a,b,c]),6)",
            &list_of(6, one_of("[a,b,c]")),
        ),
        (
            "sample(randvar([1,2,3]),10)",
            &list_of(10, integer_in(0, 2)),
        ),
        ("sample(randvar(poisson,3),10)", &list_of(10, natural)),
        (
            "sample(randvar(binomial,10,1/2),10)",
            &list_of(10, integer_in(0, 10)),
        ),
        ("sample(randvar(betad,2,3),3)", &list_of(3, unit_open)),
        ("sample(randvar(exponential,2),3)", &list_of(3, positive)),
        ("sample(randvar(cauchy,0,1),3)", &list_of(3, real)),
        ("randvector(5,normald,0,1)", &list_of(5, real)),
        (
            "randvector(3,randvar(normal,mean=5,stddev=1.5))",
            &list_of(3, real),
        ),
        ("randvector(4,randvar(poisson,3))", &list_of(4, natural)),
        (
            "randvector(4,randvar([1,2,3]))",
            &list_of(4, integer_in(0, 2)),
        ),
        ("randvector(3,normald(5,1.5))", &list_of(3, real)),
        ("sample(normald(0,1),3)", &list_of(3, real)),
        (
            "sample(exp(randvar(normal,mean=5,stddev=2)),3)",
            &list_of(3, positive),
        ),
        (
            "randvector(3,multinomial,[1/2,1/3,1/6],[a,b,c])",
            &list_of(3, one_of("[a,b,c]")),
        ),
        ("randmultinomial([1/2,1/3,1/6])", &integer_in(0, 2)),
        ("randmultinomial([1/2,1/3,1/6],[a,b,c])", &one_of("[a,b,c]")),
    ]);
}

/// Whether `e` is a number within `half` of `center`.
fn within(center: Rational, half: Rational) -> impl Fn(&Expr) -> bool {
    move |e| {
        value(e).is_some_and(|v| {
            center.clone() - half.clone() <= v && v <= center.clone() + half.clone()
        })
    }
}

/// Whether `e` is a list of counts, one per weight of `p`, summing to `n`,
/// each count `c_k` within 4 standard deviations of its binomial mean:
/// `|c_k − n·p_k| ≤ 4·√(n·p_k·(1 − p_k))`.
fn counts(n: u32, p: &'static [f64]) -> impl Fn(&Expr) -> bool {
    move |e| {
        let found: Option<Vec<i32>> = e
            .args()
            .iter()
            .map(|c| {
                value(c)
                    .filter(Rational::is_integer)
                    .and_then(|v| v.numer().to_i64())
                    .and_then(|v| i32::try_from(v).ok())
            })
            .collect();
        e.head() == Some(&Expr::from(Symbol::LIST))
            && found.is_some_and(|found| {
                let total = f64::from(n);
                found.len() == p.len()
                    && found.iter().map(|&c| i64::from(c)).sum::<i64>() == i64::from(n)
                    && found.iter().zip(p).all(|(&c, p)| {
                        (f64::from(c) - total * p).abs() <= 4.0 * (total * p * (1.0 - p)).sqrt()
                    })
            })
    }
}

/// Means and frequencies over 10,000 draws, seeded on both sides:
/// - the mean of 10,000 draws of `normald(5,1.5)` within 4 standard errors
///   of 5, `|m − 5| ≤ 4·1.5/√10000 = 3/50` (giac 4.97611995074);
/// - giac's `randvector(n,multinomial,P)` is the counts of each index over
///   `n` draws (`[13,4,3]` for `n = 20`), each within 4 standard deviations
///   of `n·p`, `4·√(n·p(1−p))`: 200, 188.6 and 149.1 about 5000, 3333.3 and
///   1666.7.
#[test]
fn frequencies() {
    check(&[
        (
            "mean(sample(randvar(normal,mean=5,stddev=1.5),10000))",
            &within(q(5), q(3) / q(50)),
        ),
        (
            "randvector(10000,multinomial,[1/2,1/3,1/6])",
            &counts(10_000, &[1.0 / 2.0, 1.0 / 3.0, 1.0 / 6.0]),
        ),
    ]);
}

/// Ours only, giac's sample failing: `randvar(student,3)` and
/// `randvar(chisquare,3)` (giac's `sample(randvar('student',3),3)` and
/// `chisquared(seq[3])` Bad Argument Value), and `randvar(multinomial,P)`,
/// which giac leaves a sequence; ours draw real, positive decimals and
/// indices from 0.
#[test]
fn departures_that_draw() {
    let cases: [(&str, &str, Property<'_>); 3] = [
        (
            "sample(randvar(student,3),3)",
            "GIAC_ERROR: sample(randvar('student',3),3) | Error: Bad Argument Value",
            &list_of(3, real),
        ),
        (
            "sample(randvar(chisquare,3),3)",
            "GIAC_ERROR: chisquared(seq[3]) | Error: Bad Argument Value",
            &list_of(3, positive),
        ),
        (
            "sample(randvar(multinomial,[1/2,1/2]),4)",
            "GIAC_ERROR: sample(randvar('multinomial',[1/2,1/2]),4) | Error: Bad Argument Value",
            &list_of(4, integer_in(0, 1)),
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, property) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours_seeded(input);
        if !found.as_ref().is_some_and(property) {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Ours only (C0044): the same seed draws the same values again.
#[test]
fn seeded_draws_repeat() {
    let mut wrong = Vec::new();
    for x in [
        "randvar(normal,mean=5,stddev=1.5)",
        "randvar(poisson,3)",
        "randvar([1/2,1/3,1/6],[a,b,c])",
    ] {
        let input = format!("[srand(7),sample({x},5),srand(7),sample({x},5)]");
        let found = ours(&input).ok();
        let repeats = found.as_ref().is_some_and(|e| {
            let a = e.args();
            a.len() == 4 && a[1] == a[3] && list_of(5, |_| true)(&a[1])
        });
        if !repeats {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
