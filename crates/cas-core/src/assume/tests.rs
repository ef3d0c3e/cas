//! Tests of assumptions and the number domain: tables in both domains, and a
//! property test that what they decide never contradicts exact evaluation at
//! points that satisfy them.

use std::time::{Duration, Instant};

use proptest::prelude::*;

use super::{Domain, Sign, Truth};
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form, evaluated in `ctx`.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn complex() -> Context {
    let mut ctx = Context::new();
    ctx.set_domain(Domain::Complex);
    ctx
}

fn assuming(ctx: &Context, facts: &[&str]) -> Context {
    facts.iter().fold(ctx.clone(), |c, f| c.assuming(read(f)))
}

#[test]
fn the_real_domain_is_the_default() {
    let ctx = Context::new();
    assert_eq!(ctx.domain(), Domain::Real);
    check(
        &ctx,
        &[
            // Symbols are real.
            ("Positive[Plus[Power[x, 2], 1]]", "True"),
            ("NonNegative[Power[x, 2]]", "True"),
            ("Negative[Power[x, 2]]", "False"),
            ("Positive[x]", "Positive[x]"),
            ("Positive[Times[x, y]]", "Positive[Times[x, y]]"),
            ("RealQ[x]", "True"),
            ("RealQ[Plus[x, Times[2, y]]]", "True"),
            ("IntegerQ[x]", "False"),
            ("ComplexQ[x]", "True"),
            ("RealQ[f[x]]", "False"),
            ("Element[x, Reals]", "True"),
            ("Element[x, Integers]", "Element[x, Integers]"),
            // Powers.
            ("Power[Power[x, 2], 1/2]", "Abs[x]"),
            ("Power[Power[x, 3], 1/3]", "x"),
            ("Power[Power[x, 2], 1/3]", "Power[x, 2/3]"),
            ("Power[Power[x, 4], 1/2]", "Power[x, 2]"),
            ("Power[Power[x, 6], 1/4]", "Power[Abs[x], 3/2]"),
            // Imaginary for a negative `x`: it stays.
            ("Power[Power[x, 3], 1/2]", "Power[Power[x, 3], 1/2]"),
            ("Power[Abs[x], 2]", "Power[x, 2]"),
            ("Power[Abs[x], 3]", "Power[Abs[x], 3]"),
            // Odd roots of negatives are real (O7).
            ("Power[-8, 1/3]", "-2"),
            ("Power[-8, 2/3]", "4"),
            ("Power[-8, 4/3]", "16"),
            ("Power[-2, 1/3]", "Times[-1, Power[2, 1/3]]"),
            ("Power[-1/27, 1/3]", "-1/3"),
            ("Power[-4, 1/2]", "Complex[0, 2]"),
            ("Power[-4, 1/4]", "Power[-4, 1/4]"),
            // The constants lie between rationals.
            ("Greater[Pi, 3]", "True"),
            ("Less[E, 3]", "True"),
            ("Greater[Plus[Pi, -3], 0]", "True"),
            ("Positive[Plus[E, -2]]", "True"),
        ],
    );
}

#[test]
fn the_complex_domain_assumes_nothing() {
    let ctx = complex();
    check(
        &ctx,
        &[
            ("Power[Power[x, 2], 1/2]", "Power[Power[x, 2], 1/2]"),
            ("Power[Power[x, 3], 1/3]", "Power[Power[x, 3], 1/3]"),
            ("Power[-8, 1/3]", "Power[-8, 1/3]"),
            ("Power[-4, 1/2]", "Complex[0, 2]"),
            (
                "Positive[Plus[Power[x, 2], 1]]",
                "Positive[Plus[1, Power[x, 2]]]",
            ),
            ("RealQ[x]", "False"),
            ("ComplexQ[x]", "True"),
            ("Element[x, Reals]", "Element[x, Reals]"),
            ("Abs[Power[x, 2]]", "Abs[Power[x, 2]]"),
            // Numbers are what they are.
            ("Greater[Pi, 3]", "True"),
            ("Positive[-1/2]", "False"),
            ("Positive[Complex[0, 1]]", "False"),
        ],
    );
    // Facts still make a value real.
    let ctx = assuming(&ctx, &["Greater[x, 0]"]);
    check(
        &ctx,
        &[
            ("Power[Power[x, 2], 1/2]", "x"),
            ("RealQ[x]", "True"),
            ("Positive[Plus[Power[x, 2], 1]]", "True"),
        ],
    );
}

#[test]
fn abs() {
    check(
        &Context::new(),
        &[
            ("Abs[-3]", "3"),
            ("Abs[-3/4]", "3/4"),
            ("Abs[Complex[3, 4]]", "5"),
            ("Abs[Complex[1, 1]]", "Power[2, 1/2]"),
            ("Abs[Complex[0, -2]]", "2"),
            ("Abs[Times[-2, x]]", "Times[2, Abs[x]]"),
            ("Abs[Abs[x]]", "Abs[x]"),
            ("Abs[Power[x, 2]]", "Power[x, 2]"),
            ("Abs[Plus[Power[x, 2], 1]]", "Plus[1, Power[x, 2]]"),
            ("Abs[Pi]", "Pi"),
            ("Abs[Times[-1, Pi]]", "Pi"),
            ("Abs[Infinity]", "Infinity"),
            ("Abs[DirectedInfinity[-1]]", "Infinity"),
            ("Abs[x]", "Abs[x]"),
            ("Abs[List[-1, 2]]", "List[1, 2]"),
        ],
    );
}

#[test]
fn bounds_decide_comparisons() {
    let ctx = assuming(&Context::new(), &["Greater[x, 0]"]);
    check(
        &ctx,
        &[
            ("Positive[x]", "True"),
            ("NonNegative[x]", "True"),
            ("Negative[x]", "False"),
            ("Greater[x, -1]", "True"),
            ("Less[x, 0]", "False"),
            ("LessEqual[x, 0]", "False"),
            ("Equal[x, 0]", "False"),
            ("Unequal[x, 0]", "True"),
            ("Greater[x, 1]", "Greater[x, 1]"),
            ("Abs[x]", "x"),
            ("Abs[Times[-3, x]]", "Times[3, x]"),
            ("Power[Power[x, 2], 1/2]", "x"),
            ("Positive[Plus[x, 1]]", "True"),
            ("Positive[Times[x, Pi]]", "True"),
            ("Positive[Power[x, -1]]", "True"),
            ("Positive[Power[x, 1/2]]", "True"),
            ("Element[x, Reals]", "True"),
            // A coefficient of known sign turns an infinity.
            ("Times[x, Infinity]", "Infinity"),
            ("Times[-2, x, Infinity]", "DirectedInfinity[-1]"),
            ("Times[x, y, Infinity]", "Times[Infinity, x, y]"),
        ],
    );
    let ctx = assuming(&Context::new(), &["Greater[x, 2]"]);
    check(
        &ctx,
        &[
            ("Greater[Plus[x, -1], 0]", "True"),
            ("Greater[x, 1]", "True"),
            ("Greater[x, 2]", "True"),
            ("GreaterEqual[x, 2]", "True"),
            ("Greater[x, 3]", "Greater[x, 3]"),
            ("Less[Times[-1, x], -2]", "True"),
            ("Greater[Times[3, x], 6]", "True"),
            ("Greater[Power[x, 2], 0]", "True"),
        ],
    );
    // Bounds on both sides, as one assumption.
    let ctx = assuming(
        &Context::new(),
        &["And[GreaterEqual[x, 0], LessEqual[x, 1]]"],
    );
    check(
        &ctx,
        &[
            ("LessEqual[x, 1]", "True"),
            ("Greater[x, 2]", "False"),
            ("NonNegative[x]", "True"),
            ("Positive[x]", "Positive[x]"),
            ("Less[Plus[x, -1], 1]", "True"),
            ("Less[0, x, 2]", "Less[0, x, 2]"),
            ("LessEqual[-1, x, 2]", "True"),
        ],
    );
}

/// A ball's range is its interval, with exact ends; one whose midpoint or
/// radius would be a fraction of more than `BALL_RANGE_BITS` bits has the
/// range of the sign it proves, where its exact ends took fractions of 10^9
/// bits at each decision (CAS-171).
#[test]
fn balls_range_by_their_ends_or_far_by_their_signs() {
    let start = Instant::now();
    let ctx = Context::new();
    let tiny = "Ball[0x1p-1000000000, 0x0p+0]";
    let huge = "Ball[0x3p1000000000, 0x1p1000000000]";
    assert_eq!(ctx.sign(&read(tiny)), Sign::Positive);
    assert_eq!(
        ctx.sign(&read(&format!("Times[-1, {huge}]"))),
        Sign::Negative
    );
    assert_eq!(
        ctx.sign(&read("Ball[0x1p-1000000000, 0x1p-999999999]")),
        Sign::Unknown
    );
    let positive = assuming(&ctx, &["Greater[x, 0]"]);
    for ball in [tiny, huge] {
        let query = format!("Greater[Plus[x, {ball}], 0]");
        assert_eq!(positive.decide(&read(&query)), Truth::True, "{query}");
    }
    let elapsed = start.elapsed();
    assert!(elapsed < Duration::from_secs(5), "{elapsed:?}");
    // 0 < x < 1/2 gives x + 1/4 < 1 by the ball's exact upper end.
    let bounded = assuming(&ctx, &["Greater[x, 0]", "Less[x, 1/2]"]);
    let query = "Less[Plus[x, Ball[0x1p-2, 0x0p+0]], 1]";
    assert_eq!(bounded.decide(&read(query)), Truth::True);
}

#[test]
fn negative_and_negated_facts() {
    let ctx = assuming(
        &Context::new(),
        &[
            "Less[x, 0]",
            "Not[Greater[y, 1]]",
            "Less[z, 0]",
            "Not[Less[w, 0]]",
        ],
    );
    check(
        &ctx,
        &[
            ("Abs[x]", "Times[-1, x]"),
            ("Negative[x]", "True"),
            ("Positive[Times[x, x]]", "True"),
            ("Negative[Power[x, 3]]", "True"),
            ("Times[x, Infinity]", "DirectedInfinity[-1]"),
            ("LessEqual[y, 1]", "True"),
            ("Greater[y, 2]", "False"),
            ("Greater[y, 0]", "Greater[y, 0]"),
            // Two negatives make a positive.
            ("Positive[Times[x, z]]", "True"),
            ("Negative[Times[x, z, Pi]]", "False"),
            ("NonNegative[w]", "True"),
            ("Negative[w]", "False"),
        ],
    );
}

#[test]
fn facts_about_differences() {
    // Written in any form, it is a bound on `x - y`.
    for fact in [
        "Greater[x, Plus[y, 1]]",
        "Greater[Plus[Times[2, x], Times[-2, y]], 2]",
        "Less[Plus[y, Times[-1, x]], -1]",
    ] {
        let ctx = assuming(&Context::new(), &[fact]);
        check(
            &ctx,
            &[
                ("Greater[x, y]", "True"),
                ("Less[y, x]", "True"),
                ("Greater[Times[2, x], Times[2, y]]", "True"),
                ("Greater[Plus[x, Times[-1, y]], 1]", "True"),
                ("Greater[Plus[x, 1], y]", "True"),
                ("Equal[x, y]", "False"),
                ("Greater[x, Plus[y, 2]]", "Greater[x, Plus[2, y]]"),
            ],
        );
    }
    // A difference that cancels decides by itself.
    check(
        &Context::new(),
        &[
            ("Greater[Plus[x, 1], x]", "True"),
            ("Less[Plus[x, 1], x]", "False"),
            ("Equal[Plus[x, 1], x]", "False"),
        ],
    );
}

#[test]
fn exclusions_and_number_sets() {
    let ctx = assuming(&Context::new(), &["Unequal[x, 2]", "Element[n, Integers]"]);
    check(
        &ctx,
        &[
            ("Unequal[x, 2]", "True"),
            ("Equal[x, 2]", "False"),
            ("Equal[Plus[x, -2], 0]", "False"),
            ("Equal[x, 3]", "Equal[x, 3]"),
            ("IntegerQ[n]", "True"),
            ("IntegerQ[Plus[n, 1]]", "True"),
            ("IntegerQ[Times[2, n]]", "True"),
            ("IntegerQ[Times[1/2, n]]", "False"),
            ("RationalQ[Times[1/2, n]]", "True"),
            ("RealQ[n]", "True"),
            ("Element[Power[n, 2], Integers]", "True"),
            ("Element[Power[n, -1], Rationals]", "True"),
            ("IntegerQ[Power[n, -1]]", "False"),
            ("Element[n, Integers]", "True"),
            ("Element[1/2, Integers]", "False"),
            ("Element[Complex[0, 1], Reals]", "False"),
        ],
    );
    // The language's `n :: integer` follows them.
    let mut ctx = ctx;
    ctx.define_delayed(
        &read("g[PatternTest[Pattern[k, Blank[]], IntegerQ]]"),
        &read("int[k]"),
    )
    .unwrap();
    check(&ctx, &[("g[n]", "int[n]"), ("g[x]", "g[x]")]);
}

#[test]
fn other_assumptions_hold_as_written() {
    let ctx = assuming(&Context::new(), &["f[x]", "Not[h[y]]"]);
    assert_eq!(ctx.decide(&read("f[x]")), Truth::True);
    assert_eq!(ctx.decide(&read("Not[f[x]]")), Truth::False);
    assert_eq!(ctx.decide(&read("h[y]")), Truth::False);
    assert_eq!(ctx.decide(&read("f[y]")), Truth::Unknown);
    assert_eq!(ctx.decide(&read("And[f[x], Greater[Pi, 3]]")), Truth::True);
    assert_eq!(ctx.decide(&read("Or[f[y], Less[Pi, 3]]")), Truth::Unknown);
    assert_eq!(ctx.decide(&read("Or[f[y], Greater[Pi, 3]]")), Truth::True);
}

#[test]
fn scopes_and_the_api() {
    let parent = Context::new();
    let child = parent.assuming(read("Greater[x, 0]"));
    let grandchild = child.assuming(read("Less[x, 1]"));
    assert_eq!(parent.sign(&read("x")), Sign::Unknown);
    assert_eq!(child.sign(&read("x")), Sign::Positive);
    assert_eq!(grandchild.decide(&read("Less[x, 2]")), Truth::True);
    assert_eq!(child.decide(&read("Less[x, 2]")), Truth::Unknown);
    assert_eq!(parent.assumptions(), &[] as &[Expr]);
    assert_eq!(grandchild.assumptions().len(), 2);
    // Assumptions are kept as written, unevaluated.
    let held = parent.assuming(read("Greater[Plus[x, x], 0]"));
    assert_eq!(held.assumptions(), &[read("Greater[Plus[x, x], 0]")]);
    assert_eq!(held.sign(&read("x")), Sign::Positive);
    // The sign of what evaluates.
    assert_eq!(parent.sign(&read("Plus[x, Times[-1, x]]")), Sign::Zero);
    assert_eq!(parent.sign(&read("Plus[Power[x, 2], 1]")), Sign::Positive);
    assert_eq!(complex().sign(&read("Plus[Power[x, 2], 1]")), Sign::Unknown);
    assert_eq!(parent.sign(&read("Complex[0, 1]")), Sign::Nonzero);
    // A new domain is a new state.
    let mut ctx = Context::new();
    let id = ctx.id();
    ctx.set_domain(Domain::Complex);
    assert_ne!(ctx.id(), id);
    assert_eq!(ctx.domain(), Domain::Complex);
}

#[test]
fn rule_conditions_see_the_facts() {
    let mut ctx = Context::new();
    ctx.define_delayed(
        &read("g[Pattern[x, Blank[]]]"),
        &read("Condition[positive, Positive[x]]"),
    )
    .unwrap();
    check(&ctx, &[("g[a]", "g[a]"), ("g[2]", "positive")]);
    let inside = ctx.assuming(read("Greater[a, 0]"));
    check(&inside, &[("g[a]", "positive"), ("g[b]", "g[b]")]);
}

/// A rational with a small numerator and denominator.
fn small() -> impl Strategy<Value = (i64, i64)> {
    (-6i64..=6, 1i64..=3)
}

fn rational((n, d): (i64, i64)) -> Expr {
    Expr::from(cas_num::Rational::new(n.into(), d.into()))
}

/// Polynomials in `x` and `y`, with `Abs` and powers that may be negative.
fn expression() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        3 => prop::sample::select(vec!["x", "y"]).prop_map(|s| Expr::from(Symbol::new(s))),
        2 => small().prop_map(rational),
        1 => Just(Expr::from(Symbol::PI)),
    ];
    leaf.prop_recursive(3, 16, 3, |inner| {
        prop_oneof![
            prop::collection::vec(inner.clone(), 2..4)
                .prop_map(|args| Expr::apply(&Symbol::PLUS.into(), &args)),
            prop::collection::vec(inner.clone(), 2..4)
                .prop_map(|args| Expr::apply(&Symbol::TIMES.into(), &args)),
            (inner.clone(), prop::sample::select(vec![-1i64, 2, 3]))
                .prop_map(|(b, n)| Expr::apply(&Symbol::POWER.into(), &[b, n.into()])),
            inner.prop_map(|x| Expr::apply(&Symbol::ABS.into(), &[x])),
        ]
    })
}

/// An assumption about `v` that `value` satisfies: a bound on either side,
/// strict or not, both, or a value it is not.
fn fact_about() -> impl Strategy<Value = (i64, i64, usize)> {
    (0i64..=3, 1i64..=2, 0usize..6)
}

fn fact(v: &str, value: (i64, i64), (gap, d, kind): (i64, i64, usize)) -> Option<Expr> {
    let x = Expr::from(Symbol::new(v));
    let point = cas_num::Rational::new(value.0.into(), value.1.into());
    let delta = cas_num::Rational::new(gap.into(), d.into());
    let (below, above) = (&point - &delta, &point + &delta);
    let compare = |head: Symbol, bound: cas_num::Rational| {
        Expr::apply(&head.into(), &[x.clone(), Expr::from(bound)])
    };
    Some(match kind {
        0 => compare(Symbol::GREATER_EQUAL, below),
        1 if !delta.is_zero() => compare(Symbol::GREATER, below),
        2 => compare(Symbol::LESS_EQUAL, above),
        3 if !delta.is_zero() => compare(Symbol::LESS, above),
        4 if !delta.is_zero() => compare(Symbol::UNEQUAL, above),
        5 => Expr::apply(
            &Symbol::AND.into(),
            &[
                compare(Symbol::GREATER_EQUAL, below),
                compare(Symbol::LESS_EQUAL, above),
            ],
        ),
        _ => return None,
    })
}

/// Queries about an expression: comparisons with a number or with another
/// expression, and the sign predicates.
fn query(e: Expr, f: Expr, c: Expr, kind: usize) -> Expr {
    let heads = [
        Symbol::GREATER,
        Symbol::GREATER_EQUAL,
        Symbol::LESS,
        Symbol::LESS_EQUAL,
        Symbol::EQUAL,
        Symbol::UNEQUAL,
    ];
    match kind {
        0..6 => Expr::apply(&heads[kind].into(), &[e, c]),
        6..12 => Expr::apply(&heads[kind - 6].into(), &[e, f]),
        12 => Expr::apply(&Symbol::POSITIVE.into(), &[e]),
        13 => Expr::apply(&Symbol::NON_NEGATIVE.into(), &[e]),
        14 => Expr::apply(&Symbol::NEGATIVE.into(), &[e]),
        _ => Expr::apply(&Symbol::NON_POSITIVE.into(), &[e]),
    }
}

/// `e` with `x` and `y` replaced by their values, evaluated without facts:
/// `Some(truth)` when it comes out `True` or `False`.
fn at(e: &Expr, x: (i64, i64), y: (i64, i64)) -> Option<bool> {
    let values = [
        (Expr::from(Symbol::new("x")), rational(x)),
        (Expr::from(Symbol::new("y")), rational(y)),
    ];
    let mut plain = Context::new();
    plain.set_domain(Domain::Complex);
    match plain.eval(&e.substitute(&values)).ok()?.as_symbol()? {
        Symbol::TRUE => Some(true),
        Symbol::FALSE => Some(false),
        _ => None,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// What `decide` says holds at every point that satisfies the
    /// assumptions, in both domains: tried at the point the assumptions were
    /// made around, and at their bounds.
    #[test]
    fn decide_never_contradicts_evaluation(
        point in (small(), small()),
        facts in (fact_about(), fact_about()),
        e in expression(),
        f in expression(),
        c in small(),
        kind in 0usize..16,
        complex_domain in any::<bool>(),
    ) {
        let (px, py) = point;
        let mut ctx = Context::new();
        if complex_domain {
            ctx.set_domain(Domain::Complex);
        }
        let assumptions: Vec<Expr> = [fact("x", px, facts.0), fact("y", py, facts.1)]
            .into_iter()
            .flatten()
            .collect();
        for a in &assumptions {
            prop_assert_eq!(at(a, px, py), Some(true), "{:?} at the point", a);
            ctx = ctx.assuming(a.clone());
        }
        let q = query(e, f, rational(c), kind);
        let answer = ctx.decide(&q);
        // Points that satisfy the assumptions: the one they were made
        // around, and others nearby.
        let candidates = [(px, py), ((px.0 * 2 + 1, px.1 * 2), py), (px, (py.0 * 3 - 1, py.1 * 3))];
        for (x, y) in candidates {
            if !assumptions.iter().all(|a| at(a, x, y) == Some(true)) {
                continue;
            }
            match (answer, at(&q, x, y)) {
                (Truth::True, Some(false)) | (Truth::False, Some(true)) => {
                    prop_assert!(false, "{:?} is {:?} under {:?}, but not at x = {:?}, y = {:?}", q, answer, assumptions, x, y);
                }
                _ => {}
            }
        }
        // `sign` agrees with the value at the point.
        if let Some(v) = ctx.eval(&q.args()[0]).ok().map(|e| e.substitute(&[])) {
            let sign = ctx.sign(&v);
            let positive = at(&Expr::apply(&Symbol::GREATER.into(), &[v.clone(), 0.into()]), px, py);
            let zero = at(&Expr::apply(&Symbol::EQUAL.into(), &[v.clone(), 0.into()]), px, py);
            let ok = match sign {
                Sign::Positive => positive != Some(false) && zero != Some(true),
                Sign::Negative => positive != Some(true) && zero != Some(true),
                Sign::Zero => zero != Some(false),
                Sign::NonNegative => !(positive == Some(false) && zero == Some(false)),
                Sign::NonPositive => positive != Some(true),
                Sign::Nonzero => zero != Some(true),
                Sign::Unknown => true,
            };
            prop_assert!(ok, "sign {:?} of {:?} under {:?}", sign, v, assumptions);
        }
    }
}
