//! Tests of the evaluation report: warnings, the steps and what took them,
//! the levels, completeness and verification.

use std::fmt::Write as _;
use std::sync::Arc;

use proptest::prelude::*;

use super::{Origin, Step, StepKind, Trace, Warning};
use crate::{
    Attributes, Context, Domain, EvalError, EvalOptions, Expr, Limits, Symbol, TraceLevel,
};

fn sym(name: &str) -> Symbol {
    Symbol::new(name)
}

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn delayed(ctx: &mut Context, lhs: &str, rhs: &str) {
    ctx.define_delayed(&read(lhs), &read(rhs)).unwrap();
}

fn up(ctx: &mut Context, tag: &str, lhs: &str, rhs: &str) {
    ctx.define_up_delayed(sym(tag), &read(lhs), &read(rhs))
        .unwrap();
}

fn traced(
    ctx: &Context,
    input: &str,
    level: TraceLevel,
) -> (Result<Expr, EvalError>, super::Report) {
    ctx.eval_with(&read(input), &EvalOptions { trace: level })
}

/// What took a step, in a word or a few.
fn kind(k: &StepKind) -> String {
    let name = |s: &Symbol| format!("{:?}", Expr::from(*s));
    match k {
        StepKind::Rule(id) => {
            let origin = match id.origin {
                Origin::User => String::new(),
                Origin::System { name, .. } => format!(" {name}"),
            };
            format!(
                "rule {:?} {} {}{origin}",
                id.kind,
                name(&id.symbol),
                id.index
            )
        }
        StepKind::Native(s) => format!("native {}", name(s)),
        StepKind::Algorithm(note) => format!("algorithm {note}"),
        StepKind::Memo { first } => format!("memo {first:?}"),
        StepKind::Normalise => "normalise".into(),
        StepKind::Evaluation => "evaluation".into(),
        StepKind::Attempts(n) => format!("attempts {n}"),
    }
}

/// The trace as a tree, a line per step and its children indented under
/// it: `position kind: before -> after`.
fn show(trace: &Trace) -> String {
    fn line(trace: &Trace, s: &Step, depth: usize, out: &mut String) {
        let indent = "  ".repeat(depth);
        let (position, kind) = (&s.position, kind(&s.kind));
        writeln!(
            out,
            "{indent}{position:?} {kind}: {:?} -> {:?}",
            s.before, s.after
        )
        .unwrap();
        for c in trace.children(s) {
            line(trace, c, depth + 1, out);
        }
    }
    let mut out = String::new();
    for s in trace.top() {
        line(trace, s, 0, &mut out);
    }
    out
}

/// Each `(input, trace)` at `level` in `ctx`, the trace as [`show`] writes
/// it.
fn check(ctx: &Context, level: TraceLevel, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (result, report) = traced(ctx, input, level);
        let found = show(&report.trace);
        if let Err(e) = consistent(&report.trace) {
            wrong.push(format!("{input}: {e}\n{found}"));
        } else if result.is_err() || found.trim_end() != expected.trim() {
            wrong.push(format!(
                "{input} ({result:?}):\n{found}expected:\n{expected}\n"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The context of the table tests: a value, exact, patterned, up- and
/// sub-rules, and a user's native function.
fn rules_context() -> Context {
    let mut ctx = Context::new();
    ctx.set(sym("x"), read("2")).unwrap();
    ctx.set(sym("pv"), read("True")).unwrap();
    ctx.define(&read("g[0]"), &read("zero")).unwrap();
    delayed(&mut ctx, "g[Pattern[n, Blank[]]]", "h[n]");
    delayed(&mut ctx, "m[Pattern[n, Blank[]]]", "Plus[n, 1]");
    up(&mut ctx, "t", "Power[t, 2]", "tt");
    delayed(
        &mut ctx,
        "s[Pattern[a, Blank[]]][Pattern[b, Blank[]]]",
        "Plus[a, b]",
    );
    delayed(
        &mut ctx,
        "c[Pattern[y, Blank[]]]",
        "Condition[y, Greater[y, 0]]",
    );
    delayed(&mut ctx, "c[Pattern[y, Blank[]]]", "Times[-1, y]");
    ctx.define_native(
        sym("nr"),
        Arc::new(|ev, _| {
            ev.eval(&Expr::apply(&Symbol::PLUS.into(), &[1.into(), 1.into()]))?;
            Ok(None)
        }),
    )
    .unwrap();
    delayed(&mut ctx, "nr[Pattern[y, Blank[]]]", "done");
    ctx.define_native(
        sym("count"),
        Arc::new(|_, e| {
            let [arg] = e.args() else { return Ok(None) };
            Ok(Some(Expr::from(i64::try_from(arg.args().len()).unwrap())))
        }),
    )
    .unwrap();
    ctx
}

#[test]
fn steps_name_what_took_them() {
    let ctx = rules_context();
    check(
        &ctx,
        TraceLevel::Rules,
        &[
            // A value, then the native of `Plus`.
            (
                "Plus[x, 1]",
                "[1] rule Value x 0: x -> 2\n[] native Plus: Plus[1, 2] -> 3",
            ),
            // An exact rule, then a patterned one, at their positions.
            (
                "List[g[0], g[1]]",
                "[1] rule Exact g 0: g[0] -> zero\n[2] rule Down g 0: g[1] -> h[1]",
            ),
            ("Power[t, 2]", "[] rule Up t 0: Power[t, 2] -> tt"),
            (
                "s[1][2]",
                "[] rule Sub s 0: s[1][2] -> Plus[1, 2]\n[] native Plus: Plus[1, 2] -> 3",
            ),
            (
                "Function[y, Plus[y, 1]][2]",
                "[] native Function: Function[y, Plus[y, 1]][2] -> Plus[2, 1]\n\
                 [] native Plus: Plus[1, 2] -> 3",
            ),
            // The value of a builtin symbol.
            (
                "Times[2, I]",
                "[2] native I: I -> Complex[0, 1]\n\
                 [] native Times: Times[2, Complex[0, 1]] -> Complex[0, 2]",
            ),
            ("count[w[a, b]]", "[] native count: count[w[a, b]] -> 2"),
            // A native's evaluations are its children, located from what
            // they evaluate.
            (
                "And[pv, qq]",
                "[] native And: And[pv, qq] -> qq\n  [] evaluation: pv -> True\n    \
                 [] rule Value pv 0: pv -> True",
            ),
            // The memo points to the first step of the earlier evaluation.
            (
                "List[m[1], m[1]]",
                "[1] rule Down m 0: m[1] -> Plus[1, 1]\n[1] native Plus: Plus[1, 1] -> 2\n\
                 [2] memo Some(0): m[1] -> 2",
            ),
            // A condition is not recorded at this level, nor the rule
            // whose condition failed, nor normalisation: `Times` gets its
            // arguments sorted.
            ("c[3]", "[] rule Down c 0: c[3] -> 3"),
            // Nor what a native that did not apply evaluated, before a rule
            // that did, or alone.
            ("nr[1]", "[] rule Down nr 0: nr[1] -> done"),
            // Steps after those dropped are numbered again, their children
            // too.
            (
                "List[nr[1], And[pv, qq]]",
                "[1] rule Down nr 0: nr[1] -> done\n[2] native And: And[pv, qq] -> qq\n  \
                 [] evaluation: pv -> True\n    [] rule Value pv 0: pv -> True",
            ),
            ("Greater[Times[2, a], a]", ""),
            (
                "c[-2]",
                "[] rule Down c 1: c[-2] -> Times[-1, -2]\n[] native Times: Times[-2, -1] -> 2",
            ),
        ],
    );
}

#[test]
fn verbose_adds_normalisation_conditions_and_attempts() {
    let ctx = rules_context();
    check(
        &ctx,
        TraceLevel::Verbose,
        &[
            ("Plus[b, a]", "[] normalise: Plus[b, a] -> Plus[a, b]"),
            (
                "nr[1]",
                "[] attempts 0: nr[1] -> nr[1]\n  [] evaluation: Plus[1, 1] -> 2\n    \
                 [] native Plus: Plus[1, 1] -> 2\n[] rule Down nr 0: nr[1] -> done",
            ),
            (
                "Greater[Times[2, a], a]",
                "[] attempts 0: Greater[Times[2, a], a] -> Greater[Times[2, a], a]\n  \
                 [] evaluation: Plus[Times[2, a], Times[-1, a]] -> a\n    \
                 [] normalise: Plus[Times[2, a], Times[-1, a]] -> Plus[Times[-1, a], Times[2, a]]\n    \
                 [] native Plus: Plus[Times[-1, a], Times[2, a]] -> a",
            ),
            (
                "c[3]",
                "[] rule Down c 0: c[3] -> 3\n  [] evaluation: Greater[3, 0] -> True\n    \
                 [] native Greater: Greater[3, 0] -> True",
            ),
            (
                "c[-2]",
                "[] attempts 1: c[-2] -> c[-2]\n  [] evaluation: Greater[-2, 0] -> False\n    \
                 [] native Greater: Greater[-2, 0] -> False\n\
                 [] rule Down c 1: c[-2] -> Times[-1, -2]\n\
                 [] normalise: Times[-1, -2] -> Times[-2, -1]\n\
                 [] native Times: Times[-2, -1] -> 2",
            ),
        ],
    );
    // Off records nothing, and Rules no normalisation.
    for (level, steps) in [
        (TraceLevel::Off, 0),
        (TraceLevel::Rules, 0),
        (TraceLevel::Verbose, 1),
    ] {
        let (_, report) = traced(&ctx, "Plus[b, a]", level);
        assert_eq!(report.trace.steps.len(), steps, "{level:?}");
    }
}

#[test]
fn natives_report_steps_and_warnings_of_their_own() {
    let mut ctx = Context::new();
    // `dream[a, b]` is `a^2 + b^2`, wrongly reported as `(a + b)^2`.
    ctx.define_native(
        sym("dream"),
        Arc::new(|ev, e| {
            let [a, b] = e.args() else { return Ok(None) };
            let square = |x: &Expr| Expr::apply(&Symbol::POWER.into(), &[x.clone(), 2.into()]);
            let sum = Expr::apply(&Symbol::PLUS.into(), &[square(a), square(b)]);
            if ev.tracing() {
                let whole = square(&Expr::apply(&Symbol::PLUS.into(), &[a.clone(), b.clone()]));
                ev.step("expand", &whole, &sum);
            }
            ev.warn(Warning::Native {
                message: "dreaming".into(),
                expr: e.clone(),
            });
            Ok(Some(sum))
        }),
    )
    .unwrap();
    let (result, report) = traced(&ctx, "dream[a, b]", TraceLevel::Rules);
    assert_eq!(result, Ok(read("Plus[Power[a, 2], Power[b, 2]]")));
    assert_eq!(
        show(&report.trace),
        "[] native dream: dream[a, b] -> Plus[Power[a, 2], Power[b, 2]]\n  \
         [] algorithm expand: Power[Plus[a, b], 2] -> Plus[Power[a, 2], Power[b, 2]]\n"
    );
    assert_eq!(report.warnings[0].to_string(), "dreaming: dream[a, b]");
    // Verification finds the step the native reported wrong.
    let verification = report.trace.verify(&ctx);
    let wrong = verification.wrong.expect("a wrong step");
    assert_eq!(wrong.step.kind, StepKind::Algorithm("expand".into()));
    assert_eq!(wrong.point.len(), 2);
    assert_ne!(wrong.before, wrong.after);
    // Without a trace, the native skips its step and still warns.
    let (_, report) = traced(&ctx, "dream[a, b]", TraceLevel::Off);
    assert!(report.trace.steps.is_empty());
    assert_eq!(report.warnings.len(), 1);
}

/// A part that warned in its own rewrite is attempted again where it is
/// met; one whose parts warned is a stamp hit in later calls, which takes no
/// step but gives the warning again (CAS-149).
#[test]
fn a_warned_part_met_again_warns_again() {
    let mut ctx = Context::new();
    delayed(
        &mut ctx,
        "f[PatternTest[Pattern[n, Blank[]], IntegerQ]]",
        "Power[n, 2]",
    );
    check(
        &ctx,
        TraceLevel::Verbose,
        &[(
            "List[f[a], g[f[a]]]",
            "
[1] attempts 1: f[a] -> f[a]
  [] evaluation: IntegerQ[a] -> False
    [] native IntegerQ: IntegerQ[a] -> False
[2, 1] attempts 1: f[a] -> f[a]
  [] evaluation: IntegerQ[a] -> False
    [] memo Some(0): IntegerQ[a] -> False",
        )],
    );
    for input in ["List[f[a], g[f[a]]]", "g[f[a]]"] {
        let (_, report) = traced(&ctx, input, TraceLevel::Verbose);
        assert_eq!(
            report.warnings,
            [Warning::NoDefinitionMatches(read("f[a]"))],
            "{input}"
        );
        assert!(report.trace.steps.is_empty(), "{input}");
    }
}

#[test]
fn warnings() {
    let mut ctx = Context::new();
    delayed(
        &mut ctx,
        "f[PatternTest[Pattern[n, Blank[]], IntegerQ]]",
        "Power[n, 2]",
    );
    up(&mut ctx, "t", "Power[t, 2]", "tt");
    let cases: &[(&str, &[Warning])] = &[
        (
            "Times[1, Power[0, -1]]",
            &[Warning::Infinite(read("Power[0, -1]"))],
        ),
        (
            "Power[0, 0]",
            &[Warning::Indeterminate(read("Power[0, 0]"))],
        ),
        (
            "Times[0, Infinity]",
            &[Warning::Indeterminate(read("Times[0, Infinity]"))],
        ),
        (
            "Plus[Infinity, Times[-1, Infinity]]",
            &[Warning::Indeterminate(read(
                "Plus[Infinity, DirectedInfinity[-1]]",
            ))],
        ),
        (
            "Power[Infinity, 0]",
            &[Warning::Indeterminate(read("Power[Infinity, 0]"))],
        ),
        // From an argument already indeterminate or infinite: nothing new.
        ("Power[Indeterminate, 2]", &[]),
        ("Plus[Indeterminate, 1]", &[]),
        ("Times[ComplexInfinity, 2]", &[]),
        (
            "Plus[f[a], f[3]]",
            &[Warning::NoDefinitionMatches(read("f[a]"))],
        ),
        ("f[1, 2]", &[Warning::NoDefinitionMatches(read("f[1, 2]"))]),
        // Each once, in order.
        (
            "List[f[a], f[b], f[a]]",
            &[
                Warning::NoDefinitionMatches(read("f[a]")),
                Warning::NoDefinitionMatches(read("f[b]")),
            ],
        ),
        // No down-rules to match: an up-rule's tag, a builtin.
        ("t[1]", &[]),
        ("Length[1]", &[]),
        // A pole of the factorial.
        ("Factorial[-1]", &[Warning::Infinite(read("Factorial[-1]"))]),
    ];
    for (input, expected) in cases {
        for level in [TraceLevel::Off, TraceLevel::Verbose] {
            let (_, report) = traced(&ctx, input, level);
            assert_eq!(report.warnings, *expected, "{input} at {level:?}");
        }
    }
    // Plain evaluation drops them, but stamps what warned with them: they
    // come again.
    assert_eq!(ctx.eval(&read("g[f[a]]")), Ok(read("g[f[a]]")));
    for _ in 0..2 {
        let (_, report) = traced(&ctx, "g[f[a]]", TraceLevel::Off);
        assert_eq!(
            report.warnings,
            [Warning::NoDefinitionMatches(read("f[a]"))]
        );
    }
    let messages = [
        (
            Warning::Infinite(read("Power[0, -1]")),
            "infinite expression 0^(-1) encountered",
        ),
        (
            Warning::Indeterminate(read("Power[0, 0]")),
            "indeterminate expression 0^0 encountered",
        ),
        (
            Warning::NoDefinitionMatches(read("f[1, 2]")),
            "no definition of f matches f[1, 2]",
        ),
    ];
    for (w, message) in messages {
        assert_eq!(w.to_string(), message);
    }
    // A large expression is not written out: this one has 2^40 parts.
    let mut big = read("f[x]");
    for _ in 0..40 {
        big = Expr::apply(&read("g"), &[big.clone(), big]);
    }
    let w = Warning::NoDefinitionMatches(big);
    assert_eq!(w.to_string(), "no definition of g matches g[…]");
}

/// A number applied as a function stays, with a warning (CAS-DOC-7 §3.1):
/// `two(x+1)` after `two := 2`.
#[test]
fn numbers_applied_as_functions_warn() {
    let mut ctx = Context::new();
    ctx.set(Symbol::new("two"), 2.into()).unwrap();
    let cases: &[(&str, &str, &[Warning])] = &[
        (
            "two[Plus[x, 1]]",
            "2[Plus[1, x]]",
            &[Warning::NumberApplied(read("2[Plus[1, x]]"))],
        ),
        (
            "Times[1/2, 4][x]",
            "2[x]",
            &[Warning::NumberApplied(read("2[x]"))],
        ),
        (
            "Complex[0, 1][]",
            "Complex[0, 1][]",
            &[Warning::NumberApplied(read("Complex[0, 1][]"))],
        ),
        ("g[a][x]", "g[a][x]", &[]),
    ];
    for (input, result, expected) in cases {
        for level in [TraceLevel::Off, TraceLevel::Verbose] {
            let (found, report) = traced(&ctx, input, level);
            assert_eq!(found, Ok(read(result)), "{input} at {level:?}");
            assert_eq!(report.warnings, *expected, "{input} at {level:?}");
        }
    }
    assert_eq!(
        Warning::NumberApplied(read("2[Plus[1, x]]")).to_string(),
        "the number 2 applied as a function in 2[1 + x]"
    );
}

#[test]
fn a_failed_evaluation_reports_the_steps_in_progress() {
    let mut ctx = Context::new();
    delayed(&mut ctx, "loop[Pattern[n, Blank[]]]", "loop[Plus[n, 1]]");
    let mut limits = ctx.limits();
    limits.iterations = 5;
    ctx.set_limits(limits);
    let (result, report) = traced(&ctx, "loop[0]", TraceLevel::Rules);
    assert_eq!(result, Err(EvalError::IterationLimit));
    let rules = report
        .trace
        .top()
        .filter(|s| matches!(s.kind, StepKind::Rule(_)))
        .count();
    assert!(rules >= 4, "{}", show(&report.trace));
    // Deep inside the evaluations of a native, each still open closes.
    let mut ctx = Context::new();
    delayed(
        &mut ctx,
        "deep[Pattern[n, Blank[]]]",
        "And[deep[Plus[n, 1]]]",
    );
    limits = ctx.limits();
    limits.recursion = 40;
    ctx.set_limits(limits);
    let (result, report) = traced(&ctx, "deep[0]", TraceLevel::Rules);
    assert_eq!(result, Err(EvalError::RecursionLimit));
    let last = report.trace.top().last().expect("a step");
    assert_eq!(last.kind, StepKind::Attempts(0));
    assert_eq!(last.before, read("And[deep[Plus[0, 1]]]"));
    let mut depth = 0;
    let mut s = last;
    while let Some(&c) = s.children.last() {
        s = &report.trace.steps[c];
        depth += 1;
    }
    assert!(depth > 4, "{}", show(&report.trace));
}

#[test]
fn verify_locates_a_wrong_rule() {
    let mut ctx = Context::new();
    // Wrong: `x^2` is not `2 x`.
    up(&mut ctx, "x", "Power[x, 2]", "Times[2, x]");
    // Right: `t^2` is `(t + 1)^2 - 2 t - 1`.
    up(
        &mut ctx,
        "t",
        "Power[t, 2]",
        "Plus[Power[Plus[t, 1], 2], Times[-2, t], -1]",
    );
    let (result, report) = traced(&ctx, "Plus[Power[t, 2], Power[x, 2]]", TraceLevel::Rules);
    assert_eq!(
        result,
        Ok(read(
            "Plus[-1, Times[-2, t], Times[2, x], Power[Plus[1, t], 2]]"
        ))
    );
    let verification = report.trace.verify(&ctx);
    let wrong = verification.wrong.expect("a wrong step");
    assert_eq!(wrong.step.position, [2]);
    assert_eq!(wrong.step.before, read("Power[x, 2]"));
    assert!(matches!(wrong.step.kind, StepKind::Rule(id) if id.symbol == sym("x")));
    let [(s, v)] = &wrong.point[..] else {
        panic!("{:?}", wrong.point)
    };
    assert_eq!(*s, sym("x"));
    assert_eq!(
        wrong.before,
        ctx.eval(&read("Power[x, 2]").substitute(&[(sym("x").into(), v.clone())]))
            .unwrap()
    );
    // The steps before it agreed: those of `t`'s rule.
    assert!(verification.passed >= 1);
    assert!(verification.unchecked.is_empty());
    // The same again: the points are drawn the same way.
    assert_eq!(report.trace.verify(&ctx).wrong.unwrap().point, wrong.point);
    // Without the wrong step, everything agrees.
    let (_, report) = traced(&ctx, "Power[t, 2]", TraceLevel::Rules);
    let verification = report.trace.verify(&ctx);
    assert_eq!(verification.wrong, None);
    assert_eq!(verification.passed, report.trace.steps.len());
}

#[test]
fn verify_draws_points_that_satisfy_the_assumptions() {
    let mut ctx = Context::new();
    // Right only where `u >= 0`, `w` is real, `q` is an integer.
    up(&mut ctx, "u", "Abs[u]", "u");
    up(&mut ctx, "w", "Element[w, Reals]", "True");
    up(&mut ctx, "q", "Element[q, Integers]", "True");
    let mut complex = ctx.clone();
    complex.set_domain(Domain::Complex);
    // What verification says of the one step `input` takes: points that
    // satisfy the assumptions must be found for it to pass.
    let outcome = |ctx: &Context, input: &str| {
        let (_, report) = traced(ctx, input, TraceLevel::Rules);
        let verification = report.trace.verify(ctx);
        match (
            &verification.wrong,
            verification.passed,
            verification.unchecked.len(),
        ) {
            (Some(_), 0, 0) => "wrong",
            (None, 1, 0) => "passed",
            (None, 0, 1) => "unchecked",
            other => panic!("{input}: {other:?}"),
        }
    };
    let cases: &[(&Context, &str, &str)] = &[
        (&ctx, "Abs[u]", "wrong"),
        (&ctx.assuming(read("Greater[u, 0]")), "Abs[u]", "passed"),
        // A fact that ties `u` to `v`: both get values.
        (&ctx.assuming(read("Greater[u, v]")), "Abs[u]", "wrong"),
        (
            &ctx.assuming(read("And[Greater[u, v], Greater[v, 0]]")),
            "Abs[u]",
            "passed",
        ),
        // Real points in the Real domain, complex ones in the Complex
        // domain.
        (&ctx, "Element[w, Reals]", "passed"),
        (&complex, "Element[w, Reals]", "wrong"),
        (
            &complex.assuming(read("Element[w, Reals]")),
            "Element[w, Reals]",
            "passed",
        ),
        // Integer points where the facts say integer.
        (&ctx, "Element[q, Integers]", "wrong"),
        (
            &ctx.assuming(read("Element[q, Integers]")),
            "Element[q, Integers]",
            "passed",
        ),
    ];
    for (ctx, input, expected) in cases {
        assert_eq!(outcome(ctx, input), *expected, "{input}");
    }
}

/// Points evaluate with small numbers only (`Limits::number_bits`): two
/// sides equal once computed, but not in form, tell nothing where their
/// powers would have over 2^12 bits.
#[test]
fn verify_leaves_large_powers_alone() {
    let mut ctx = Context::new();
    // `3^(200000 x)` reported as `9^(100000 x)`: equal, and at every point
    // (`|x| >= 1/12`) over 2^12 bits, though at most 2^24.
    ctx.define_native(
        sym("big"),
        Arc::new(|ev, e| {
            let [x] = e.args() else { return Ok(None) };
            let power = |b: i64, k: i64| {
                let exponent = Expr::apply(&Symbol::TIMES.into(), &[k.into(), x.clone()]);
                Expr::apply(&Symbol::POWER.into(), &[b.into(), exponent])
            };
            ev.step("square", &power(3, 200_000), &power(9, 100_000));
            Ok(Some(power(9, 100_000)))
        }),
    )
    .unwrap();
    let (_, report) = traced(&ctx, "big[x]", TraceLevel::Rules);
    let verification = report.trace.verify(&ctx);
    assert_eq!(verification.wrong, None);
    let algorithm: Vec<&StepKind> = verification.unchecked.iter().map(|s| &s.kind).collect();
    assert_eq!(algorithm, [&StepKind::Algorithm("square".into())]);
    // What the fuzz target found (CAS-60): its points made powers of
    // millions of bits, minutes of verification.
    let (_, report) = traced(
        &ctx,
        "Power[Power[Power[Power[Power[-3/4, -3], 7], Power[Plus[c, c, 4], 5]], -3], -3]",
        TraceLevel::Rules,
    );
    assert_eq!(report.trace.verify(&ctx).wrong, None);
}

/// Balls at a point compare by overlap (CAS-85): the sides of a right step
/// agree however differently they round, and a ball apart from the other
/// side's, or without its exact value, makes the step wrong. Such steps went
/// unchecked.
#[test]
fn verify_compares_balls_by_overlap() {
    let mut ctx = Context::new();
    // Right: `x^2` is `(x + 0.5)^2 - x - 0.25`, a ball at each point that
    // holds the exact square.
    up(
        &mut ctx,
        "x",
        "Power[x, 2]",
        "Plus[Power[Plus[x, Decimal[5, -1]], 2], Times[-1, x], Decimal[-25, -2]]",
    );
    // Wrong: the `- 0.25` left out.
    up(
        &mut ctx,
        "y",
        "Power[y, 2]",
        "Plus[Power[Plus[y, Decimal[5, -1]], 2], Times[-1, y]]",
    );
    let (_, report) = traced(&ctx, "Power[x, 2]", TraceLevel::Rules);
    let verification = report.trace.verify(&ctx);
    assert_eq!(verification.wrong, None);
    assert!(
        !verification
            .unchecked
            .iter()
            .any(|s| matches!(s.kind, StepKind::Rule(_))),
        "{}",
        show(&report.trace)
    );
    let (_, report) = traced(&ctx, "Power[y, 2]", TraceLevel::Rules);
    let wrong = report.trace.verify(&ctx).wrong.expect("a wrong step");
    assert_eq!(wrong.step.before, read("Power[y, 2]"));
    // Steps between balls, and lists of them.
    let tenth = "Times[Decimal[1, -1], z]";
    let cases = [
        (
            "double",
            format!("Plus[{tenth}, {tenth}]"),
            "Times[Decimal[2, -1], z]".to_owned(),
            true,
        ),
        (
            "triple",
            format!("Plus[{tenth}, {tenth}]"),
            "Times[Decimal[3, -1], z]".to_owned(),
            false,
        ),
        (
            "pair",
            format!("List[{tenth}, 1]"),
            format!("List[{tenth}, 2]"),
            false,
        ),
    ];
    for (name, before, after, right) in cases {
        let (before, after) = (read(&before), read(&after));
        let returned = after.clone();
        ctx.define_native(
            sym(name),
            Arc::new(move |ev, _| {
                ev.step(name, &before, &after);
                Ok(Some(returned.clone()))
            }),
        )
        .unwrap();
        let (_, report) = traced(&ctx, &format!("{name}[z]"), TraceLevel::Rules);
        let verification = report.trace.verify(&ctx);
        let algorithm = StepKind::Algorithm(name.into());
        if right {
            assert_eq!(verification.wrong, None, "{name}");
            assert!(
                !verification.unchecked.iter().any(|s| s.kind == algorithm),
                "{name}"
            );
        } else {
            let wrong = verification.wrong.expect("a wrong step");
            assert_eq!(wrong.step.kind, algorithm, "{name}");
        }
    }
}

#[test]
fn verify_flags_what_it_cannot_check() {
    let mut ctx = Context::new();
    delayed(&mut ctx, "k[Pattern[n, Blank[]]]", "Sin[n]");
    // A binder: no values can be put in.
    let (_, report) = traced(&ctx, "Function[y, Plus[y, 1]][z]", TraceLevel::Rules);
    let verification = report.trace.verify(&ctx);
    assert_eq!(verification.unchecked.len(), 1);
    assert_eq!(
        verification.unchecked[0].kind,
        StepKind::Native(Symbol::FUNCTION)
    );
    // `k[a]` and `Sin[a]` are equal at each point though not numbers.
    let (_, report) = traced(&ctx, "k[a]", TraceLevel::Rules);
    assert_eq!(report.trace.verify(&ctx).passed, 1);
    // No point satisfies a contradiction.
    let contradiction = ctx.assuming(read("And[Greater[a, 0], Less[a, 0]]"));
    let (_, report) = traced(&contradiction, "k[a]", TraceLevel::Rules);
    assert_eq!(report.trace.verify(&contradiction).unchecked.len(), 1);
}

/// The context of the completeness test: user heads with each attribute,
/// values, and rules of each kind, with conditions, that terminate.
fn property_context() -> Context {
    let mut ctx = Context::new();
    for (name, attributes) in [
        ("p_flat", Attributes::FLAT),
        ("p_orderless", Attributes::ORDERLESS),
        ("p_listable", Attributes::LISTABLE),
        ("p_hold", Attributes::HOLD_ALL),
    ] {
        ctx.set_attributes(sym(name), attributes).unwrap();
    }
    ctx.set(sym("p_one"), read("1")).unwrap();
    ctx.set(sym("p_list"), read("List[1, p_one]")).unwrap();
    ctx.define(&read("f[0]"), &read("1")).unwrap();
    delayed(
        &mut ctx,
        "f[Pattern[a, Blank[]], Pattern[b, Blank[]]]",
        "Plus[a, b]",
    );
    delayed(
        &mut ctx,
        "f[Pattern[a, Blank[]]]",
        "Condition[Times[2, a], Greater[a, 0]]",
    );
    delayed(&mut ctx, "f[Pattern[a, Blank[List]]]", "Length[a]");
    up(&mut ctx, "u", "g[u, Pattern[a, Blank[]]]", "List[a]");
    delayed(
        &mut ctx,
        "s[Pattern[a, Blank[]]][Pattern[b, Blank[]]]",
        "f[a, b]",
    );
    delayed(
        &mut ctx,
        "k[Pattern[n, Blank[]]]",
        "Condition[k[Plus[n, -1]], Greater[n, 0]]",
    );
    // `q` evaluates `a + a` for each argument `a`, then does not apply; a
    // rule takes two arguments.
    ctx.define_native(
        sym("q"),
        Arc::new(|ev, e| {
            for a in e.args() {
                ev.eval(&Expr::apply(&Symbol::PLUS.into(), &[a.clone(), a.clone()]))?;
            }
            Ok(None)
        }),
    )
    .unwrap();
    delayed(&mut ctx, "q[Pattern[a, Blank[]], Pattern[b, Blank[]]]", "b");
    ctx
}

/// Random expressions over the symbols and rules of [`property_context`]
/// and builtins whose natives evaluate, with templates that reach each kind
/// of rule, the memo and nested evaluations often.
fn expr() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "b", "u", "p_one", "p_list", "True", "False"])
            .prop_map(|n| Expr::from(sym(n))),
        (-2i64..4).prop_map(Expr::from),
    ];
    let apply = |head: &str, args: &[Expr]| Expr::apply(&sym(head).into(), args);
    leaf.prop_recursive(4, 32, 4, move |inner| {
        prop_oneof![
            6 => (
                prop::sample::select(vec![
                    "p_flat",
                    "p_orderless",
                    "p_listable",
                    "p_hold",
                    "f",
                    "q",
                    "Plus",
                    "Times",
                    "Power",
                    "List",
                    "Sequence",
                    "Hold",
                    "Not",
                    "Greater",
                    "If",
                ]),
                prop::collection::vec(inner.clone(), 0..4),
            )
                .prop_map(move |(head, args)| apply(head, &args)),
            // The exact, patterned and conditional rules of `f` and `k`.
            3 => (
                prop::sample::select(vec!["f", "k"]),
                prop_oneof![3 => inner.clone(), 1 => Just(Expr::from(0))],
            )
                .prop_map(move |(head, a)| apply(head, &[a])),
            1 => inner.clone().prop_map(move |a| apply("g", &[sym("u").into(), a])),
            1 => (inner.clone(), inner.clone())
                .prop_map(move |(a, b)| Expr::apply(&apply("s", &[a]), &[b])),
            1 => inner.clone().prop_map(|a| Expr::apply(&read("Function[y, Plus[y, 1]]"), &[a])),
            1 => inner.clone().prop_map(move |a| apply("ReplaceAll", &[a, read("Rule[a, b]")])),
            // A part twice, for the memo.
            1 => inner.clone().prop_map(move |a| apply("List", &[a.clone(), a])),
            // Natives that evaluate held parts.
            1 => (prop::sample::select(vec!["And", "Or"]), inner.clone(), inner)
                .prop_map(move |(head, a, b)| apply(head, &[a, b])),
        ]
    })
}

/// Replays `steps`, located in `start`, and checks each evaluation among
/// their children: the expression they make, or what went wrong.
fn replay<'t>(
    trace: &'t Trace,
    start: &Expr,
    steps: impl Iterator<Item = &'t Step>,
) -> Result<Expr, String> {
    let mut e = start.clone();
    for s in steps {
        match s.kind {
            StepKind::Algorithm(_) => continue,
            StepKind::Evaluation => return Err(format!("an evaluation located in {e:?}: {s:?}")),
            _ => {}
        }
        let part = e
            .part(&s.position)
            .ok_or_else(|| format!("no part {s:?} in {e:?}"))?;
        if *part != s.before {
            return Err(format!(
                "{part:?} at {:?} in {e:?} is not {:?}",
                s.position, s.before
            ));
        }
        e = e
            .replace_part(&s.position, s.after.clone())
            .expect("the part exists");
        for c in trace.children(s) {
            if c.kind == StepKind::Evaluation {
                let r = replay(trace, &c.before, trace.children(c))?;
                if r != c.after {
                    return Err(format!(
                        "{:?} replays to {r:?}, not {:?}",
                        c.before, c.after
                    ));
                }
            }
        }
    }
    Ok(e)
}

/// The tree a trace describes: parents and children agree, ids are places.
fn consistent(trace: &Trace) -> Result<(), String> {
    for (i, s) in trace.steps.iter().enumerate() {
        if s.id != i {
            return Err(format!("step {i} has id {}", s.id));
        }
        for &c in &s.children {
            if c >= i || trace.steps[c].parent != Some(i) {
                return Err(format!("child {c} of {i}"));
            }
        }
        if let StepKind::Memo { first: Some(f) } = s.kind {
            if f >= i {
                return Err(format!("memo {i} points at {f}"));
            }
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4000))]

    /// Tracing changes no result; at the Verbose level the top steps
    /// replayed from the input give the result, and each evaluation
    /// recorded replays to its value.
    #[test]
    fn a_verbose_trace_replays_to_the_result(e in expr()) {
        let mut ctx = property_context();
        ctx.set_limits(Limits {
            recursion: 200,
            iterations: 200,
            ..Limits::default()
        });
        let plain = ctx.eval(&e);
        for level in [TraceLevel::Rules, TraceLevel::Verbose] {
            let (result, report) = ctx.eval_with(&e, &EvalOptions { trace: level });
            prop_assert_eq!(&result, &plain, "{:?} at {:?}", e, level);
            prop_assert!(consistent(&report.trace).is_ok(), "{:?}", consistent(&report.trace));
            if level == TraceLevel::Rules {
                prop_assert!(report.trace.steps.iter().all(|s| s.kind != StepKind::Normalise));
            }
            if let (TraceLevel::Verbose, Ok(result)) = (level, &result) {
                let replayed = replay(&report.trace, &e, report.trace.top());
                prop_assert_eq!(replayed.as_ref(), Ok(result), "{:?}\n{}", e, show(&report.trace));
            }
        }
    }
}
