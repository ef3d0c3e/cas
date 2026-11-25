//! Tests of Markov chains (CAS-628): `markov`'s four lists on the rows of
//! the prototype (crates/cas-core/proto/markov/markov.py, giac's answers made exact), its
//! balls on decimals, each warning, and seeded `randmarkov` paths and
//! matrices by their properties.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr, Integer, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression,
/// with no warning.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        let found = found.unwrap();
        let expected = eval(&ctx, expected);
        if found != expected || !report.warnings.is_empty() {
            wrong.push(format!(
                "{input}: {found:?}, expected {expected:?}, {:?}",
                report.warnings
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head(), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

/// A ball's midpoint in `f64`, or a small integer's value.
fn float(e: &Expr) -> f64 {
    if let View::Ball(b) = e.view() {
        return b.mid().to_f64();
    }
    let i = e
        .as_i64()
        .unwrap_or_else(|| panic!("{e:?}: no small number"));
    i32::try_from(i).map(f64::from).expect("a small integer")
}

fn is_ball(e: &Expr) -> bool {
    matches!(e.view(), View::Ball(_))
}

const HALVES: &str = "List[List[1/2, 1/2], List[1/3, 2/3]]";

/// The prototype's rows, exact: giac's `1.0` entries are 1, its iterated
/// absorption `3812798742493/7625597484987` is 1/2.
#[test]
fn exact_chains() {
    let row = |m: &str, answer: &str| (format!("List[MarkovChain[{m}]]"), answer.to_owned());
    let cases = [
        row(
            HALVES,
            "List[List[List[1, 0]], List[List[2/5, 3/5]], List[], List[List[1], List[1]]]",
        ),
        (
            format!("List[MarkovChain[{HALVES}, 1/1000000]]"),
            "List[List[List[1, 0]], List[List[2/5, 3/5]], List[], List[List[1], List[1]]]".into(),
        ),
        // Columns summing to 1: transposed.
        row(
            "List[List[1/2, 1/3], List[1/2, 2/3]]",
            "List[List[List[1, 0]], List[List[2/5, 3/5]], List[], List[List[1], List[1]]]",
        ),
        row(
            "List[List[0, 1], List[1, 0]]",
            "List[List[List[1, 0]], List[List[1/2, 1/2]], List[], List[List[1], List[1]]]",
        ),
        row(
            "List[List[1/2, 1/2], List[1/2, 1/2]]",
            "List[List[List[1, 0]], List[List[1/2, 1/2]], List[], List[List[1], List[1]]]",
        ),
        row(
            "List[List[0, 1, 0], List[0, 0, 1], List[1, 0, 0]]",
            "List[List[List[2, 1, 0]], List[List[1/3, 1/3, 1/3]], List[], \
             List[List[1], List[1], List[1]]]",
        ),
        row(
            "List[List[1/2, 1/2, 0, 0], List[1/2, 1/2, 0, 0], List[0, 0, 1/3, 2/3], \
             List[1/4, 1/4, 1/4, 1/4]]",
            "List[List[List[1, 0]], List[List[1/2, 1/2, 0, 0]], List[List[3, 2]], \
             List[List[1], List[1], List[1], List[1]]]",
        ),
        row(
            "List[List[0, 0, 1/2, 0, 1/2], List[0, 0, 1, 0, 0], List[1/4, 1/4, 0, 1/4, 1/4], \
             List[0, 0, 1/2, 0, 1/2], List[0, 0, 0, 0, 1]]",
            "List[List[List[4]], List[List[0, 0, 0, 0, 1]], List[List[3, 1, 2, 0]], \
             List[List[1], List[1], List[1], List[1], List[1]]]",
        ),
        row(
            "List[List[1, 0, 0], List[1/3, 1/3, 1/3], List[0, 0, 1]]",
            "List[List[List[0], List[2]], List[List[1, 0, 0], List[0, 0, 1]], List[List[1]], \
             List[List[1, 0], List[1/2, 1/2], List[0, 1]]]",
        ),
        row(
            "List[List[1]]",
            "List[List[List[0]], List[List[1]], List[], List[List[1]]]",
        ),
        row(
            "List[List[1/2, 1/2, 0], List[1/2, 1/2, 0], List[0, 0, 1]]",
            "List[List[List[1, 0], List[2]], List[List[1/2, 1/2, 0], List[0, 0, 1]], List[], \
             List[List[1, 0], List[1, 0], List[0, 1]]]",
        ),
        // A transient state feeding a periodic class: giac's iteration
        // never ends.
        row(
            "List[List[1, 0, 0, 0], List[0, 0, 1, 0], List[0, 1, 0, 0], List[1/2, 1/6, 1/3, 0]]",
            "List[List[List[0], List[2, 1]], List[List[1, 0, 0, 0], List[0, 1/2, 1/2, 0]], \
             List[List[3]], List[List[1, 0], List[0, 1], List[0, 1], List[1/2, 1/2]]]",
        ),
    ];
    let cases: Vec<(&str, &str)> = cases
        .iter()
        .map(|(a, b)| (a.as_str(), b.as_str()))
        .collect();
    check(&cases);
}

/// Decimal chains answer balls (giac's floats, within their error): the
/// classes exact, the probabilities balls, the padding exact 0.
#[test]
fn decimal_chains() {
    let ctx = Context::new();
    let answer = |m: &str| {
        eval(&ctx, &format!("List[MarkovChain[{m}]]"))
            .args()
            .to_vec()
    };
    // [[0.5,0.5],[0.2,0.8]]: π = [2/7, 5/7].
    let found =
        answer("List[List[0x1p-1, 0x1p-1], List[0x1.999999999999ap-3, 0x1.999999999999ap-1]]");
    assert_eq!(found[0], eval(&ctx, "List[List[1, 0]]"));
    let pi = found[1].args()[0].args();
    assert!(pi.iter().all(is_ball), "{pi:?}");
    assert!((float(&pi[0]) - 2.0 / 7.0).abs() < 1e-12, "{pi:?}");
    assert!((float(&pi[1]) - 5.0 / 7.0).abs() < 1e-12, "{pi:?}");
    assert_eq!(found[2], eval(&ctx, "List[]"));
    assert_eq!(found[3], eval(&ctx, "List[List[1], List[1]]"));
    // [[0.3,0.7,0],[0.2,0.3,0.5],[0,0,1]]: giac's [0,0,1.0].
    let found = answer(
        "List[List[0x1.3333333333333p-2, 0x1.6666666666666p-1, 0], \
         List[0x1.999999999999ap-3, 0x1.3333333333333p-2, 0x1p-1], List[0, 0, 1]]",
    );
    assert_eq!(found[0], eval(&ctx, "List[List[2]]"));
    let pi = found[1].args()[0].args();
    assert_eq!(&pi[..2], &[Expr::from(0), Expr::from(0)], "{pi:?}");
    assert!(
        is_ball(&pi[2]) && (float(&pi[2]) - 1.0).abs() < 1e-12,
        "{pi:?}"
    );
    assert_eq!(found[2], eval(&ctx, "List[List[1, 0]]"));
    // [[1,0,0],[0.25,0.5,0.25],[0,0,1]]: giac's [0.5,0.5].
    let found = answer("List[List[1, 0, 0], List[0x1p-2, 0x1p-1, 0x1p-2], List[0, 0, 1]]");
    assert_eq!(found[0], eval(&ctx, "List[List[0], List[2]]"));
    let b = found[3].args();
    assert_eq!(b[0], eval(&ctx, "List[1, 0]"));
    assert!(
        b[1].args()
            .iter()
            .all(|x| is_ball(x) && (float(x) - 0.5).abs() < 1e-12),
        "{b:?}"
    );
    assert_eq!(b[2], eval(&ctx, "List[0, 1]"));
}

/// giac's errors stay with its message; ours only, the limits.
#[test]
fn warnings() {
    let ctx = Context::new();
    let cases = [
        (
            "MarkovChain[List[List[1/2, 1/3], List[1/3, 2/3]]]",
            "Not a stochastic matrix!",
        ),
        (
            "MarkovChain[List[List[1/2, 1/2], List[x, Plus[1, Times[-1, x]]]]]",
            "Not a stochastic matrix!",
        ),
        ("MarkovChain[List[1, 2, 3]]", "Not a stochastic matrix!"),
        (
            "MarkovChain[List[List[1/2, 1/2], List[-1/3, 4/3]]]",
            "Not a stochastic matrix!",
        ),
        (
            "MarkovChain[List[List[1/2, 1/2]]]",
            "Not a stochastic matrix!",
        ),
        (
            "MarkovChain[List[List[0x1p-1, 0x1.00000000dbe70p-1], List[0x1.999999999999ap-3, 0x1.999999999999ap-1]]]",
            "Not a stochastic matrix!",
        ),
        (
            "MarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], x]",
            "Bad Argument Value",
        ),
        (
            "MarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], Plus[1, I]]",
            "Bad Argument Value",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 0]",
            "Bad Argument Value",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 2, 10]",
            "Invalid dimension",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], List[0, 2], 4]",
            "Invalid dimension",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], -1, 4]",
            "Invalid dimension",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 0, -1]",
            "Bad Argument Value",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 1/2, 4]",
            "Bad Argument Value",
        ),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[x, Plus[1, Times[-1, x]]]], 0, 4]",
            "Not a stochastic matrix!",
        ),
        ("RandomMarkovChain[List[0, 2], 1]", "Invalid dimension"),
        ("RandomMarkovChain[List[2], -1]", "Bad Argument Value"),
        ("RandomMarkovChain[List[1000], 1]", "too many states"),
        (
            "RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 0, 1000000]",
            "too many draws",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = warning(&ctx, input);
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // A chain of more entries than `terms`.
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits { terms: 3, ..limits });
    assert_eq!(
        warning(&small, &format!("MarkovChain[{HALVES}]")),
        "too many states"
    );
    // A plot argument stays without a warning.
    let e = read("RandomMarkovChain[List[List[1/2, 1/2], List[1/3, 2/3]], 0, 4, plot]");
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap(), e);
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);
}

/// Seeded paths: their length, states of positive probability, absorbing
/// states, and state 1's share near 3/5.
#[test]
fn paths() {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[42]");
    let states = |input: &str| -> Vec<i64> {
        let found = eval(&ctx, input);
        found
            .args()
            .iter()
            .map(|s| s.as_i64().unwrap_or_else(|| panic!("{input}: {found:?}")))
            .collect()
    };
    assert_eq!(states(&format!("RandomMarkovChain[{HALVES}, 0, 0]")), [0]);
    assert_eq!(
        states("RandomMarkovChain[List[List[1, 0], List[0, 1]], 1, 5]"),
        [1, 1, 1, 1, 1, 1]
    );
    let walk = states(
        "RandomMarkovChain[List[List[0, 0, 1/2, 0, 1/2], List[0, 0, 1, 0, 0], \
         List[1/4, 1/4, 0, 1/4, 1/4], List[0, 0, 1/2, 0, 1/2], List[0, 0, 0, 0, 1]], 2, 20]",
    );
    assert_eq!(walk.len(), 21);
    assert_eq!(walk[0], 2);
    let allowed: [&[i64]; 5] = [&[2, 4], &[2], &[0, 1, 3, 4], &[2, 4], &[4]];
    for pair in walk.windows(2) {
        let from = usize::try_from(pair[0]).expect("a state");
        assert!(allowed[from].contains(&pair[1]), "{walk:?}");
    }
    let decimal =
        states("RandomMarkovChain[List[List[0x1p-1, 0x1p-1], List[0x1p-2, 0x1.8p-1]], 0, 50]");
    assert_eq!(decimal.len(), 51);
    assert!(decimal.iter().all(|s| (0..2).contains(s)), "{decimal:?}");
    let long = states(&format!("RandomMarkovChain[{HALVES}, 0, 10000]"));
    assert_eq!(long.len(), 10001);
    let ones = long.iter().filter(|&&s| s == 1).count();
    let share = f64::from(u32::try_from(ones).expect("a count")) / 10001.0;
    // The chain mixes fast (second eigenvalue 1/6): 0.03 is many standard
    // errors (about 0.006).
    assert!((share - 0.6).abs() < 0.03, "{share}");
    let rows = eval(&ctx, &format!("RandomMarkovChain[{HALVES}, List[0, 1], 4]"));
    assert_eq!(rows.args().len(), 2);
    for (row, start) in rows.args().iter().zip([0, 1]) {
        assert_eq!(row.args().len(), 5);
        assert_eq!(row.args()[0], Expr::from(start));
    }
    // A seed repeats a path.
    let twice = eval(
        &ctx,
        &format!(
            "List[SeedRandom[7], RandomMarkovChain[{HALVES}, 0, 30], \
             SeedRandom[7], RandomMarkovChain[{HALVES}, 0, 30]]"
        ),
    );
    assert_eq!(twice.args()[1], twice.args()[3]);
}

/// The fuzz harness's way (CAS-628): `Context::seed_random` before each of
/// two evaluations repeats the draws of `randmarkov`.
#[test]
fn reseeded_draws_repeat() {
    let ctx = Context::new();
    let input = read("Table[RandomMarkovChain[x, 1], x, -2, 3]");
    ctx.seed_random(&Integer::from(0));
    let first = ctx.eval(&input).unwrap();
    ctx.seed_random(&Integer::from(0));
    assert_eq!(ctx.eval(&input).unwrap(), first);
}

/// `randmarkov(v, nt)`: rows of decimals in `[0, 1]` summing to 1, and
/// `markov` of the answer finds classes of the sizes in `v` and `nt`
/// transient states.
#[test]
fn random_chains() {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[42]");
    for (v, nt, sizes) in [
        ("List[1, 2, 1, 3]", 4, &[1usize, 2, 1, 3][..]),
        ("List[1, 2]", 1, &[1, 2]),
        ("List[2]", 0, &[2]),
        ("List[1]", 2, &[1]),
        ("3", 1, &[3]),
    ] {
        let m = eval(&ctx, &format!("RandomMarkovChain[{v}, {nt}]"));
        let n: usize = sizes.iter().sum::<usize>() + nt;
        assert_eq!(m.args().len(), n, "{m:?}");
        for row in m.args() {
            assert_eq!(row.args().len(), n);
            assert!(
                row.args()
                    .iter()
                    .all(|x| is_ball(x) && (0.0..=1.0).contains(&float(x)))
            );
            let sum: f64 = row.args().iter().map(float).sum();
            assert!((sum - 1.0).abs() < 1e-12, "{row:?}");
        }
        let chain = Expr::call(
            Symbol::LIST,
            [Expr::call(Symbol::MARKOV_CHAIN, [m.clone()])],
        );
        let found = ctx.eval(&chain).unwrap();
        let mut recurrent: Vec<usize> = found.args()[0]
            .args()
            .iter()
            .map(|c| c.args().len())
            .collect();
        recurrent.sort_unstable();
        let mut expected = sizes.to_vec();
        expected.sort_unstable();
        assert_eq!(recurrent, expected, "{found:?}");
        let transient: usize = found.args()[2].args().iter().map(|c| c.args().len()).sum();
        assert_eq!(transient, nt, "{found:?}");
    }
}
