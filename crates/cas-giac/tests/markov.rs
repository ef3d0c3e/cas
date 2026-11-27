//! Xcas's Markov chains against giac (CAS-628): `markov`, its four lists
//! bracketed into one (C0023), and `randmarkov`. Exact rows agree with
//! giac's; where giac prints decimals (decimal input, or its iteration's
//! `1.0` and float-iterated rationals) ours is within giac's 12 digits, or
//! within `eps` where giac iterates to it (C0006), and exact where the
//! input is (R67). Invalid calls stay, giac's message the warning (libgiac's
//! text, C0007, C0077). giac's generator is not ours: draws are properties
//! both answers have, each side seeded in the same input (C0014).
//! Departures: giac's answer as caseval prints it, whitespace collapsed,
//! and ours (crates/cas-core/proto/giac/probe1.cas, probe10.cas,
//! crates/cas-core/proto/markov/).

use cas_core::{Context, EvalOptions, Expr, Rational, Symbol, View, Warning};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, decimal_balls, lower, parse};
use proptest::prelude::*;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// A number as a rational: itself, or a ball's midpoint.
fn value(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Ball(b) => b.mid().checked_to_rational(),
        _ => e.to_rational(),
    }
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(n.into(), d.into())
}

/// Whether `ours` is `theirs`: lists item by item, numbers the same or
/// within `tolerance` relative to `max(1, |theirs|)`.
fn near(ours: &Expr, theirs: &Expr, tolerance: &Rational) -> bool {
    if ours == theirs {
        return true;
    }
    let list = Expr::from(Symbol::LIST);
    if ours.head() == Some(&list) || theirs.head() == Some(&list) {
        return ours.head() == theirs.head()
            && ours.args().len() == theirs.args().len()
            && ours
                .args()
                .iter()
                .zip(theirs.args())
                .all(|(a, b)| near(a, b, tolerance));
    }
    match (value(ours), value(theirs)) {
        (Some(a), Some(b)) => {
            let scale = b.abs().max(q(1, 1));
            (&a - &b).abs() <= &scale * tolerance
        }
        _ => false,
    }
}

/// Whether `e` holds a decimal.
fn has_ball(e: &Expr) -> bool {
    matches!(e.view(), View::Ball(_)) || e.args().iter().any(has_ball)
}

/// Each input: giac's answer, read back, is ours (`==`, C0032).
#[test]
fn exact_chains() {
    let Some(giac) = giac() else {
        return;
    };
    let inputs = [
        "[markov([[1/2,1/2],[1/3,2/3]])]",
        "[markov([[1/2,1/2],[1/3,2/3]],1e-6)]",
        "[markov([[1/2,1/2],[1/3,2/3]],1/2)]",
        "[markov([[1/2,1/2],[1/3,2/3]],sqrt(2))]",
        "[markov([[1/2,1/2],[1/3,2/3]],-1)]",
        "[markov([[1/2,1/3],[1/2,2/3]])]",
        "[markov([[0,1],[1,0]])]",
        "[markov([[1/2,1/2],[1/2,1/2]])]",
        "[markov([[0,1,0],[0,0,1],[1,0,0]])]",
        "[markov([[1/2,1/2,0,0],[1/2,1/2,0,0],[0,0,1/3,2/3],[1/4,1/4,1/4,1/4]])]",
        "randmarkov([[1/2,1/2],[1/3,2/3]],0,0)",
        "randmarkov([[1,0],[0,1]],1,5)",
        "randmarkov([[1/2,1/2],[1/3,2/3]],[],4)",
    ];
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each `(input, tolerance)`: giac's decimals and ours within it (giac's
/// 12 digits, or the `eps` it iterates to), ours exact on exact input.
#[test]
fn decimal_chains() {
    let Some(giac) = giac() else {
        return;
    };
    let digits = q(1, 1_000_000_000);
    let eps = q(1, 1000);
    let cases = [
        // Exact input, giac's iteration: `1.0`, 3812798742493/7625597484987.
        (
            "[markov([[0,0,1/2,0,1/2],[0,0,1,0,0],[1/4,1/4,0,1/4,1/4],[0,0,1/2,0,1/2],[0,0,0,0,1]])]",
            &digits,
        ),
        ("[markov([[1,0,0],[1/3,1/3,1/3],[0,0,1]])]", &digits),
        ("[markov([[1]])]", &digits),
        ("[markov([[1/2,1/2,0],[1/2,1/2,0],[0,0,1]])]", &digits),
        (
            "[markov([[1,0,0,0],[0,0,1,0],[0,1,0,0],[1/3,1/3,1/3,0]])]",
            &digits,
        ),
        ("[markov([[1,0,0],[1/4,1/2,1/4],[0,0,1]],1e-3)]", &eps),
        // Decimal input: balls.
        ("[markov([[0.5,0.5],[0.2,0.8]])]", &digits),
        ("[markov([[0.5,0.5+1e-13],[0.2,0.8]])]", &digits),
        ("[markov([[0.3,0.7,0],[0.2,0.3,0.5],[0,0,1]])]", &digits),
        ("[markov([[1,0,0],[0.25,0.5,0.25],[0,0,1]])]", &digits),
        ("[markov([[0.5,0.5],[0,1]])]", &digits),
        ("[markov([[0.9,0.1],[0.5,0.5]],1e-3)]", &eps),
    ];
    let mut wrong = Vec::new();
    for (input, tolerance) in cases {
        let exact = !input.contains('.');
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine))
                if near(&mine, &theirs, tolerance) && (!exact || !has_ball(&mine)) => {}
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

/// The message giac's error carries: after giac's ` | ` and `Error: `, or
/// the sentence before ` Error:` (C0007).
fn message(error: &str) -> &str {
    let text = error.strip_prefix("GIAC_ERROR: ").unwrap_or(error);
    if let Some((_, rest)) = text.split_once('|') {
        return rest.trim().strip_prefix("Error: ").unwrap_or(rest.trim());
    }
    if let Some(rest) = text.strip_prefix("Error: ") {
        return rest;
    }
    text.split_once(" Error:").map_or(text, |(m, _)| m)
}

/// Giac's text, whitespace collapsed and a lone `|` dropped (C0077).
fn collapsed(s: &str) -> String {
    s.split_whitespace()
        .filter(|w| *w != "|")
        .collect::<Vec<_>>()
        .join(" ")
}

/// Whether ours of `input` stays, its head the input's, with a warning of
/// `expected`, or with none for `None`.
fn stays_with(input: &str, expected: Option<&str>) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && match expected {
                Some(expected) => warnings
                    .iter()
                    .any(|w| matches!(w, Warning::Native { message, .. } if message == expected)),
                None => warnings.is_empty(),
            }
    })
}

/// Invalid calls stay with giac's message as the warning (giac's text in
/// `crates/cas-core/proto/giac/probe1.cas`, `probe10.cas`,
/// `crates/cas-core/proto/markov/probe_a.cas`,
/// `probe_b.cas`).
#[test]
fn invalid_calls_stay() {
    const STOCHASTIC: &str = "GIAC_ERROR: Not a stochastic matrix! Error: Bad Argument Value";
    let rows = [
        ("markov([[1/2,1/3],[1/3,2/3]])", STOCHASTIC),
        ("markov([[1/2,1/2],[x,1-x]])", STOCHASTIC),
        ("markov([1,2,3])", STOCHASTIC),
        ("markov([[1/2,1/2],[-1/3,4/3]])", STOCHASTIC),
        ("markov([[1/2,1/2]])", STOCHASTIC),
        ("markov([[0.5,0.5+1e-10],[0.2,0.8]])", STOCHASTIC),
        ("markov([[0.5,0.5+1e-6],[0.2,0.8]])", STOCHASTIC),
        (
            "markov([[1/2,1/2],[1/3,2/3]],x)",
            "GIAC_ERROR: markov([[1/2,1/2],[1/3,2/3]],x)  |  Error: Bad Argument Value",
        ),
        (
            "markov([[1/2,1/2],[1/3,2/3]],1+i)",
            "GIAC_ERROR: markov([[1/2,1/2],[1/3,2/3]],1+i)  |  Error: Bad Argument Value",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],0)",
            "GIAC_ERROR: randmarkov([[1/2,1/2],[1/3,2/3]],0)  |  Error: Bad Argument Value",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],2,10)",
            "GIAC_ERROR: randmarkov([[1/2,1/2],[1/3,2/3]],2,10)  |  Error: Invalid dimension",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],[0,2],4)",
            "GIAC_ERROR: randmarkov([[1/2,1/2],[1/3,2/3]],[0,2],4)  |  Error: Invalid dimension",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],-1,4)",
            "GIAC_ERROR: randmarkov([[1/2,1/2],[1/3,2/3]],-1,4)  |  Error: Invalid dimension",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],0,-1)",
            "GIAC_ERROR: randmarkov([[1/2,1/2],[1/3,2/3]],0,-1)  |  Error: Bad Argument Value",
        ),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],1/2,4)",
            "GIAC_ERROR: Error: Bad Argument Value",
        ),
        ("randmarkov([[1/2,1/2],[x,1-x]],0,4)", STOCHASTIC),
        (
            "randmarkov([0,2],1)",
            "GIAC_ERROR: randmarkov([0,2],1)  |  Error: Invalid dimension",
        ),
        (
            "randmarkov([2],-1)",
            "GIAC_ERROR: randmarkov([2],-1)  |  Error: Bad Argument Value",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, error) in rows {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
        if !stays_with(input, Some(message(error))) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {:?}",
                evaluated(input),
                message(error)
            ));
        }
    }
    // Ours only: past the limits, and a plot (giac a graphic) without a
    // warning.
    for (input, expected) in [
        ("randmarkov([1000],1)", Some("too many states")),
        (
            "randmarkov([[1/2,1/2],[1/3,2/3]],0,10^6)",
            Some("too many draws"),
        ),
        ("randmarkov([[1/2,1/2],[1/3,2/3]],0,4,plot)", None),
    ] {
        if !stays_with(input, expected) {
            wrong.push(format!("{input}: ours {:?}", evaluated(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures: giac's answer as recorded, and ours, the value of the last
/// form: giac's floats where ours is exact; giac's never-ending iteration
/// for a transient state feeding a periodic class (ours only).
#[test]
fn departures() {
    let cases = [
        (
            "[markov([[1]])]",
            "[[[0]],[[1.0]],[],[[1]]]",
            "[[[0]],[[1]],[],[[1]]]",
        ),
        (
            "[markov([[1,0,0],[1/3,1/3,1/3],[0,0,1]])]",
            "[[[0],[2]],[[1.0,0,0],[0,0,1.0]],[[1]],\
             [[1,0],[3812798742493/7625597484987,3812798742493/7625597484987],[0,1]]]",
            "[[[0],[2]],[[1,0,0],[0,0,1]],[[1]],[[1,0],[1/2,1/2],[0,1]]]",
        ),
        (
            "[markov([[1,0,0],[1/4,1/2,1/4],[0,0,1]],1e-3)]",
            "[[[0],[2]],[[1.0,0,0],[0,0,1.0]],[[1]],[[1,0],[2047/4096,2047/4096],[0,1]]]",
            "[[[0],[2]],[[1,0,0],[0,0,1]],[[1]],[[1,0],[1/2,1/2],[0,1]]]",
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
    // giac's loop never ends here (crates/cas-core/proto/markov/probe_b.cas, a 15 s
    // timeout): ours only.
    let periodic = ours("[markov([[1,0,0,0],[0,0,1,0],[0,1,0,0],[1/2,1/6,1/3,0]])]");
    let expected =
        ours("[[[0],[2,1]],[[1,0,0,0],[0,1/2,1/2,0]],[[3]],[[1,0],[0,1],[0,1],[1/2,1/2]]]");
    if periodic.is_err() || periodic.as_ref().ok() != expected.as_ref().ok() {
        wrong.push(format!("periodic: ours {periodic:?}"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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

/// Ours of `input` seeded with 42 in the same input.
fn ours_seeded(input: &str) -> Option<Expr> {
    ours(&format!("[srand(42),{input}]"))
        .ok()
        .and_then(|e| e.args().get(1).cloned())
}

/// A property of an answer.
type Property<'p> = &'p dyn Fn(&Expr) -> bool;

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

/// The states of a path, `None` unless each is a small integer.
fn states(e: &Expr) -> Option<Vec<i64>> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) {
        return None;
    }
    e.args().iter().map(Expr::as_i64).collect()
}

/// A path of `n` states from `start`, each step allowed by `allowed`.
fn path(n: usize, start: i64, allowed: &[&[i64]]) -> impl Fn(&Expr) -> bool {
    move |e| {
        states(e).is_some_and(|s| {
            s.len() == n
                && s[0] == start
                && s.windows(2).all(|w| {
                    usize::try_from(w[0])
                        .ok()
                        .and_then(|i| allowed.get(i))
                        .is_some_and(|next| next.contains(&w[1]))
                })
        })
    }
}

/// A stochastic matrix of `n` states in `randmarkov(v, nt)`'s layout:
/// decimals in `[0, 1]`, rows summing to 1 within giac's printed digits.
fn layout(n: usize) -> impl Fn(&Expr) -> bool {
    move |m| {
        m.args().len() == n
            && m.args().iter().all(|row| {
                let entries: Option<Vec<Rational>> = row.args().iter().map(value).collect();
                row.args().len() == n
                    && row.args().iter().all(|x| matches!(x.view(), View::Ball(_)))
                    && entries.is_some_and(|x| {
                        let sum = x.iter().fold(q(0, 1), |s, x| &s + x);
                        x.iter().all(|x| q(0, 1) <= *x && *x <= q(1, 1))
                            && (&sum - &q(1, 1)).abs() <= q(1, 1_000_000_000)
                    })
            })
    }
}

/// `markov`'s answer, bracketed: recurrent classes of the sizes in `v`
/// (in any order) and `nt` transient states.
fn classes(v: &'static [usize], nt: usize) -> impl Fn(&Expr) -> bool {
    move |e| {
        let Some([recurrent, _, transient, _]) = e.args().get(..4) else {
            return false;
        };
        let mut sizes: Vec<usize> = recurrent.args().iter().map(|c| c.args().len()).collect();
        sizes.sort_unstable();
        let mut expected = v.to_vec();
        expected.sort_unstable();
        sizes == expected
            && transient
                .args()
                .iter()
                .map(|c| c.args().len())
                .sum::<usize>()
                == nt
    }
}

/// Seeded draws (giac's in `crates/cas-core/proto/giac/probe1.cas`, `probe10.cas`,
/// `crates/cas-core/proto/markov/probe_f.cas`): paths by length, start and steps of positive
/// probability; random chains by layout and by `markov`'s classes.
#[test]
fn draws() {
    let aide = "[[0,0,1/2,0,1/2],[0,0,1,0,0],[1/4,1/4,0,1/4,1/4],[0,0,1/2,0,1/2],[0,0,0,0,1]]";
    let aide_path = format!("randmarkov({aide},2,20)");
    let aide_steps: [&[i64]; 5] = [&[2, 4], &[2], &[0, 1, 3, 4], &[2, 4], &[4]];
    let two: [&[i64]; 2] = [&[0, 1], &[0, 1]];
    let halves = path(11, 1, &two);
    let aide_walk = path(21, 2, &aide_steps);
    let long = |e: &Expr| states(e).is_some_and(|s| s.len() == 1001);
    // State 1's share over 10,000 steps near 3/5: the chain mixes fast
    // (second eigenvalue 1/6), the standard error is about 0.006, the bound
    // 0.03 (giac 6043/10001).
    let share = |e: &Expr| {
        states(e).is_some_and(|s| {
            let ones = s.iter().filter(|&&x| x == 1).count();
            s.len() == 10001 && (5700..=6300).contains(&ones)
        })
    };
    let rows = [
        ("randmarkov([1,2,1,3],4)", layout(11)),
        ("randmarkov([1,2],1)", layout(4)),
        ("randmarkov([2],0)", layout(2)),
        ("randmarkov([1],2)", layout(3)),
    ];
    let recovered = [
        (
            "[markov(randmarkov([1,2,1,3],4))]",
            classes(&[1, 2, 1, 3], 4),
        ),
        ("[markov(randmarkov([1,2],1))]", classes(&[1, 2], 1)),
        ("[markov(randmarkov([2],0))]", classes(&[2], 0)),
        ("[markov(randmarkov([1],2))]", classes(&[1], 2)),
    ];
    let mut cases: Vec<(&str, Property<'_>)> = vec![
        (aide_path.as_str(), &aide_walk),
        ("randmarkov([[1/2,1/2],[1/3,2/3]],1,10)", &halves),
        ("randmarkov([[0.5,0.5],[0.25,0.75]],1,10)", &halves),
        ("randmarkov([[1/2,1/2],[1/3,2/3]],0,1000)", &long),
        ("randmarkov([[1/2,1/2],[1/3,2/3]],0,10000)", &share),
    ];
    for (input, property) in &rows {
        cases.push((input, property));
    }
    for (input, property) in &recovered {
        cases.push((input, property));
    }
    check(&cases);
}

/// Ours only: a seed repeats a path and a random chain.
#[test]
fn seeded_draws_repeat() {
    for input in [
        "[srand(7),randmarkov([[1/2,1/2],[1/3,2/3]],0,30),srand(7),randmarkov([[1/2,1/2],[1/3,2/3]],0,30)]",
        "[srand(7),randmarkov([1,2],1),srand(7),randmarkov([1,2],1)]",
    ] {
        let found = ours(input).unwrap_or_else(|e| panic!("{input}: {e}"));
        assert_eq!(found.args()[1], found.args()[3], "{input}");
    }
}

/// giac's answer to `input`, read back, is ours; else the failure names
/// the input, giac's answer and ours.
fn agree(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    giac.agrees(input)
        .map_err(|e| TestCaseError::fail(e.to_string()))
}

/// Row weights of a chain of 2 to 5 states: each weight 0 half the time
/// (so chains with transient states and several classes come up), else 1
/// to 4; each row's sum positive.
fn weights() -> impl Strategy<Value = Vec<Vec<u32>>> {
    (2usize..=5).prop_flat_map(|n| {
        prop::collection::vec(
            prop::collection::vec(prop_oneof![Just(0u32), 1u32..=4], n)
                .prop_filter("a row sums to more than 0", |r| r.iter().sum::<u32>() > 0),
            n,
        )
    })
}

/// Whether giac answers the chain of row weights `w` exactly, read off its
/// graph: no absorbing state, and one recurrent class or no transient
/// state (see `random_markov_chains`).
fn exact_in_giac(w: &[Vec<u32>]) -> bool {
    let n = w.len();
    let mut reach: Vec<Vec<bool>> = (0..n)
        .map(|i| (0..n).map(|j| i == j || w[i][j] > 0).collect())
        .collect();
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if reach[i][k] && reach[k][j] {
                    reach[i][j] = true;
                }
            }
        }
    }
    // A recurrent state reaches only its own class.
    let recurrent: Vec<usize> = (0..n)
        .filter(|&i| (0..n).all(|j| !reach[i][j] || reach[j][i]))
        .collect();
    let classes = recurrent
        .iter()
        .filter(|&&i| (0..i).all(|j| !(reach[i][j] && reach[j][i])))
        .count();
    recurrent
        .iter()
        .all(|&i| (0..n).filter(|&j| reach[i][j]).count() > 1)
        && (classes == 1 || recurrent.len() == n)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `markov` of a random chain (`weights`), entries `w/sum` exact,
    /// bracketed (C0023). Skipped by the drawn chain's graph (C0116) where
    /// giac departs (`departures`): an absorbing state (giac's invariant
    /// vector `1.0`, `markov([[1,0],[1/2,1/2]])`); two or more recurrent
    /// classes with a transient state (giac iterates the absorption: `0.0`,
    /// `2199023255551/2199023255552`, `markov([[1/2,1/2,0],[0,1,0],[0,0,1]])`;
    /// or never ends where a transient state feeds a periodic class,
    /// `markov([[0,0,0,0,1],[0,1,0,0,0],[0,0,0,1,0],[0,0,1,0,0],
    /// [2/11,0,4/11,1/11,4/11]])`, 15 s). In 7000 draws run in one libgiac
    /// (T4.c1) exactly the skipped ones had floats in giac's answer; about
    /// a third are skipped.
    #[test]
    fn random_markov_chains(w in weights()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        if !exact_in_giac(&w) {
            return Ok(());
        }
        let rows: Vec<String> = w
            .iter()
            .map(|r| {
                let sum: u32 = r.iter().sum();
                let entries: Vec<String> = r.iter().map(|x| format!("{x}/{sum}")).collect();
                format!("[{}]", entries.join(","))
            })
            .collect();
        agree(giac, &format!("[markov([{}])]", rows.join(",")))?;
    }
}
