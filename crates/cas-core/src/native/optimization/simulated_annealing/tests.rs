//! Tests of `SimulatedAnnealing` (CAS-610): a seeded run repeats, a convex
//! cost's minimum, the user's functions applied, the schedule's checks and
//! calls that stay.

use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The answer and the warnings; `None` for a call that stays.
fn run(ctx: &Context, input: &str) -> (Option<Expr>, usize) {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    let found = found.unwrap_or_else(|err| panic!("{input}: {err:?}"));
    if found.head() == e.head() {
        return (None, report.warnings.len());
    }
    (Some(found), report.warnings.len())
}

/// A number as a double: a decimal's midpoint, or a rational.
fn double(e: &Expr) -> f64 {
    if let View::Ball(b) = e.view() {
        return b.mid().to_f64();
    }
    let q = e
        .to_rational()
        .unwrap_or_else(|| panic!("a number, not {e:?}"));
    cas_num::ball::Float::from_rational(&q, 53, cas_num::ball::Round::Nearest)
        .0
        .to_f64()
}

/// `(x − 1)²` from 5, steps drawn in `[x − m, x + m]`, `m = 1/2`; 74
/// temperatures of 20 steps (1480 at most).
fn convex(schedule: &str, iterations: &str) -> String {
    format!(
        "SimulatedAnnealing[5, Function[List[x], Power[Plus[x, -1], 2]], \
         Function[List[x, y], Abs[Plus[x, Times[-1, y]]]], \
         Function[List[x, m], Random[Plus[x, Times[-1, m]], Plus[x, m]]], \
         {schedule}, {iterations}, 1/2]"
    )
}

const SCHEDULE: &str = "List[1, 1, 11/10, 1/1000]";
const ITERATIONS: &str = "List[1, 20]";

/// Seeded twice, the same answer.
#[test]
fn seeded_run_repeats() {
    let ctx = Context::new();
    let input = convex(SCHEDULE, ITERATIONS);
    let mut answers = Vec::new();
    for _ in 0..2 {
        ctx.eval(&read("SeedRandom[5]")).unwrap();
        let (found, warnings) = run(&ctx, &input);
        assert_eq!(warnings, 0, "{input}");
        answers.push(found.unwrap_or_else(|| panic!("{input}: stays")));
    }
    assert_eq!(answers[0], answers[1]);
}

/// The minimum 1 within `10⁻²` under three seeds (a prototype of the loop:
/// within `2.5·10⁻³` under 200 seeds).
#[test]
fn convex_minimum() {
    let ctx = Context::new();
    let input = convex(SCHEDULE, ITERATIONS);
    for seed in [1, 2, 3] {
        ctx.eval(&read(&format!("SeedRandom[{seed}]"))).unwrap();
        let (found, warnings) = run(&ctx, &input);
        let x = double(&found.unwrap_or_else(|| panic!("{input}: stays")));
        assert_eq!(warnings, 0, "{input}");
        assert!((x - 1.0).abs() < 1e-2, "seed {seed}: {x}");
    }
}

/// The user's step and cost applied: steps of −1 (or −`maxstep`) always
/// lower the cost `x`, 2 temperatures (1, 1/2) of 3 steps each.
#[test]
fn functions_applied() {
    let ctx = Context::new();
    for (input, answer) in [
        (
            "SimulatedAnnealing[0, Function[List[x], x], Function[List[x, y], 0], \
             Function[List[x], Plus[x, -1]], List[1, 1, 2, 1/2], List[1, 3]]",
            "-6",
        ),
        (
            "SimulatedAnnealing[0, Function[x, x], Function[List[x, y], 0], \
             Function[List[x, m], Plus[x, Times[-1, m]]], List[1, 1, 2, 1/2], List[1, 3], 2]",
            "-12",
        ),
    ] {
        assert_eq!(run(&ctx, input), (Some(read(answer)), 0), "{input}");
    }
}

/// The schedule's checks: a warning and the call stays.
#[test]
fn bad_schedules_warn() {
    let ctx = Context::new();
    for (schedule, iterations) in [
        ("List[0, 1, 11/10, 1/1000]", ITERATIONS),
        ("List[-1, 1, 11/10, 1/1000]", ITERATIONS),
        ("List[1, 1, 1, 1/1000]", ITERATIONS),
        ("List[1, 1, 1/2, 1/1000]", ITERATIONS),
        ("List[1, 1, 11/10, 1]", ITERATIONS),
        ("List[1, 1, 11/10, 2]", ITERATIONS),
        ("List[1, 1, 11/10, 0]", ITERATIONS),
        ("List[1, 1, 11/10, -1/1000]", ITERATIONS),
        // 74 temperatures of 10⁶ steps: past `terms`.
        (SCHEDULE, "List[1, 1000000]"),
    ] {
        let input = convex(schedule, iterations);
        assert_eq!(run(&ctx, &input), (None, 1), "{input}");
    }
    // A cost that is no real number.
    let input = "SimulatedAnnealing[0, Function[List[x], y], Function[List[x, y], 0], \
                 Function[List[x], Plus[x, -1]], List[1, 1, 2, 1/2], List[1, 3]]";
    assert_eq!(run(&ctx, input), (None, 1), "{input}");
}

/// Malformed calls stay without a warning.
#[test]
fn bad_shapes_stay() {
    let ctx = Context::new();
    let cost = "Function[List[x], Power[x, 2]]";
    let dist = "Function[List[x, y], 0]";
    let step1 = "Function[List[x], Plus[x, 1]]";
    let step2 = "Function[List[x, m], Plus[x, m]]";
    for args in [
        format!("1, {cost}, {dist}, {step1}"),
        format!("1, {dist}, {dist}, {step1}, {SCHEDULE}"),
        format!("1, {cost}, {cost}, {step1}, {SCHEDULE}"),
        format!("1, {cost}, {dist}, {step2}, {SCHEDULE}"),
        format!("1, {cost}, {dist}, {step1}, {SCHEDULE}, {ITERATIONS}, 1"),
        format!("1, {cost}, {dist}, {step2}, {SCHEDULE}, {ITERATIONS}, 0"),
        format!("1, {cost}, {dist}, {step1}, {SCHEDULE}, List[1, 0]"),
        format!("1, {cost}, {dist}, {step1}, {SCHEDULE}, List[0, 1]"),
        format!("1, {cost}, {dist}, {step1}, List[1, 1, 11/10]"),
        format!("1, {cost}, {dist}, {step1}, List[1, a, 11/10, 1/1000]"),
        format!("1, Function[Power[Slot[1], 2]], {dist}, {step1}, {SCHEDULE}"),
        format!("1, {cost}, {dist}, {step2}, {SCHEDULE}, {ITERATIONS}, 1, 1"),
    ] {
        let input = format!("SimulatedAnnealing[{args}]");
        assert_eq!(run(&ctx, &input), (None, 0), "{input}");
    }
}
