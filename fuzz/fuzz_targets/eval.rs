//! The cas-core evaluator on random expressions, under random rules,
//! attributes, defaults and assumptions, with small limits (CAS-60,
//! CAS-DOC-4 §10). Evaluation never panics, overflows the stack (each
//! input runs on a thread with just the stack its recursion limit allows)
//! or runs past its time budget; a result evaluates to itself; a trace at the
//! Verbose level changes no result, describes a consistent tree, and its
//! steps replayed from the input give the result; verification and the
//! warnings' messages never panic.
//!
//! The time budget follows from the limits only if rules cannot branch:
//! `f[x_] := f[g[x]] + f[h[x]]` makes `2^depth` distinct calls, and the
//! evaluator has no budget per call (CAS-DOC-4 R8). So a right-hand side
//! applies a user head at most once and passes none as an argument, and
//! has none of the natives that apply something to each part (`Map`,
//! `Apply`, `ReplaceAll`, `Function`, and `Unapply`, which makes one); the
//! expression evaluated may. It
//! uses each pattern variable at most once, and a rule on a `Flat` head
//! has no tests or conditions (see [`pattern`]).
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{
    Attributes, Context, Domain, EvalOptions, Expr, Limits, STACK_PER_LEVEL, Step, StepKind,
    Symbol, Trace, TraceLevel,
};
use libfuzzer_sys::fuzz_target;

/// The heads rules are defined on.
const HEADS: [&str; 4] = ["f", "g", "h", "k"];

/// Atoms besides numbers.
const ATOMS: [&str; 12] = [
    "a", "b", "c", "x", "y", "u", "I", "Pi", "E", "Infinity", "True", "False",
];

/// Builtins anywhere.
const BUILTINS: [&str; 133] = [
    "Plus",
    "Times",
    "Power",
    "List",
    "Sequence",
    "Hold",
    "And",
    "Or",
    "Not",
    "If",
    "Greater",
    "Less",
    "Equal",
    "SameQ",
    "Length",
    "Part",
    "Abs",
    "Element",
    "Positive",
    "ElementwiseTimes",
    "Sqrt",
    "DirectedInfinity",
    "IntegerQ",
    "FreeQ",
    "D",
    // What the language lowers to (CAS-123).
    "At",
    "ReplaceAt",
    "Range",
    "FiniteSet",
    "Factorial",
    "Equation",
    "ElementwiseDivide",
    "ElementwisePower",
    // Balls (CAS-83), decimal literals and `exact` (CAS-148).
    "Ball",
    "Decimal",
    "Exact",
    // `N` and the elementary functions of balls (CAS-84); `N`'s precision
    // stays within the limit on numbers.
    "N",
    "Exp",
    "Log",
    "Sin",
    "Cos",
    "ArcSin",
    "ArcTan",
    // The polynomial builtins (CAS-140, CAS-156), within a small limit on
    // terms.
    "Expand",
    "Together",
    "GCD",
    "LCM",
    "PolynomialQuotient",
    "PolynomialRemainder",
    "PolynomialQuotientRemainder",
    "Exponent",
    "Coefficient",
    "LeadingCoefficient",
    "Numerator",
    "Denominator",
    "Content",
    "PrimitivePart",
    "Resultant",
    "Discriminant",
    // Factoring (CAS-185, CAS-187), its searches within the limit on terms.
    "Factor",
    "FactorRational",
    "FactorList",
    "FactorSquareFree",
    // Integrals (CAS-194), series and limits (CAS-203, CAS-205), and
    // equations (CAS-210, CAS-214).
    "Integrate",
    "Series",
    "Limit",
    "Solve",
    "CSolve",
    "FSolve",
    // Differential equations (CAS-259).
    "DSolve",
    // Number theory (CAS-267 to CAS-271).
    "PrimeQ",
    "NextPrime",
    "PreviousPrime",
    "IntegerFactor",
    "IntegerFactorList",
    "Divisors",
    "EulerPhi",
    "PowerMod",
    "ExtendedGCD",
    "ChineseRemainder",
    "JacobiSymbol",
    "LegendreSymbol",
    "BernoulliB",
    // Linear algebra (CAS-218 to CAS-222).
    "Transpose",
    "IdentityMatrix",
    "Tr",
    "Dot",
    "Cross",
    "Det",
    "Inverse",
    "RowReduce",
    "MatrixRank",
    "NullSpace",
    "ColumnSpace",
    "LinSolve",
    "CharacteristicPolynomial",
    "Eigenvalues",
    "Eigenvectors",
    // Sums and products (CAS-224, CAS-225).
    "Sum",
    "Product",
    "Binomial",
    // Simplification (CAS-229, CAS-233).
    "Simplify",
    "TrigExpand",
    "TrigReduce",
    "TrigToExp",
    "ExpToTrig",
    // Hyperbolic functions (CAS-238).
    "Sinh",
    "Cosh",
    "Tanh",
    "ArcSinh",
    "ArcCosh",
    "ArcTanh",
    // Reciprocal trigonometric functions (CAS-239).
    "Sec",
    "Csc",
    "Cot",
    "ArcSec",
    "ArcCsc",
    "ArcCot",
    // Complex parts (CAS-240).
    "Re",
    "Im",
    "Conjugate",
    "Arg",
    // Rounding, sign and order (CAS-241).
    "Floor",
    "Ceiling",
    "Round",
    "Sign",
    "Max",
    "Min",
    // Special functions (CAS-243).
    "Gamma",
    "Beta",
    "Zeta",
    "Erf",
    "Erfc",
];

/// Builtins that apply something to parts, only in the expression
/// evaluated, and `Unapply`, which makes a `Function`. Not
/// `ReplaceRepeated`: a rule that matches what it makes, such as an integer
/// becoming `-3 + y`, grows the expression at each of its rounds.
const APPLYING: [&str; 5] = ["Map", "Apply", "ReplaceAll", "MatchQ", "Unapply"];

/// The names pattern variables take.
const VARIABLES: [&str; 3] = ["n", "m", "s"];

/// What a thread of its own keeps besides its levels, as cas-core's
/// `STACK_RESERVE`.
const RESERVE: usize = 256 << 10;

fn sym(name: &str) -> Expr {
    Symbol::new(name).into()
}

fn apply(head: &str, args: &[Expr]) -> Expr {
    Expr::apply(&sym(head), args)
}

/// A number, an atom, or a pattern variable when `room` has one.
fn leaf(u: &mut Unstructured<'_>, room: &mut Room) -> Result<Expr> {
    Ok(match u.int_in_range(0u8..=9)? {
        0..=3 => Expr::from(i64::from(u.int_in_range(-3i8..=9)?)),
        4 => Expr::parse_full_form(*u.choose(&[
            "1/2",
            "-3/4",
            "Complex[1, 2]",
            "2/3",
            // Balls (CAS-83): 0.1 at 53 bits, an exact 2.0, an inexact
            // ball around 0, a complex ball.
            "Ball[0x1.999999999999ap-4, 0x1p-57]",
            "Ball[0x1p+1, 0x0p+0]",
            "Ball[0x0p+0, 0x1p-50]",
            "Complex[Ball[0x1p-1, 0x1p-60], Ball[0x1p+0, 0x0p+0]]",
        ])?)
        .expect("a number"),
        5 if room.variables => {
            let i = u.choose_index(VARIABLES.len())?;
            match &mut room.used {
                Some(used) if used[i] => sym(u.choose(&ATOMS)?),
                Some(used) => {
                    used[i] = true;
                    sym(VARIABLES[i])
                }
                None => sym(VARIABLES[i]),
            }
        }
        _ => sym(u.choose(&ATOMS)?),
    })
}

/// What an expression may contain: the calls to user heads left (none
/// once spent; `None` for any number), whether natives that apply to
/// parts are allowed, whether pattern variables appear as atoms, and, in
/// a right-hand side, those used: each at most once, since `k[n___] :=
/// k[n, n, n]` triples its arguments at each rewrite.
#[derive(Clone, Copy)]
struct Room {
    calls: Option<u8>,
    applying: bool,
    variables: bool,
    used: Option<[bool; VARIABLES.len()]>,
}

/// No calls, no natives that apply to parts, variables as often as wanted:
/// the tests of conditions.
const TESTS: Room = Room {
    calls: Some(0),
    applying: false,
    variables: true,
    used: None,
};

fn expr(u: &mut Unstructured<'_>, depth: u8, room: &mut Room) -> Result<Expr> {
    if depth == 0 || u.ratio(1, 3)? {
        return leaf(u, room);
    }
    let count = u.int_in_range(0u8..=3)?;
    let args = |u: &mut Unstructured<'_>, room: &mut Room| -> Result<Vec<Expr>> {
        (0..count).map(|_| expr(u, depth - 1, room)).collect()
    };
    Ok(match u.int_in_range(0u8..=11)? {
        0..=5 => {
            let head = *u.choose(&BUILTINS)?;
            apply(head, &args(u, room)?)
        }
        6 | 7 if room.calls != Some(0) => {
            room.calls = room.calls.map(|c| c - 1);
            let head = *u.choose(&HEADS)?;
            apply(head, &args(u, room)?)
        }
        // A sub-rule's subject, `s[…][…]`.
        8 if room.calls != Some(0) => {
            room.calls = room.calls.map(|c| c - 1);
            let inner = Expr::apply(&sym("s"), &args(u, room)?);
            Expr::apply(&inner, &[expr(u, depth - 1, room)?])
        }
        9 if room.applying => {
            let head = *u.choose(&APPLYING)?;
            let subject = expr(u, depth - 1, room)?;
            let second = match head {
                "Map" | "Apply" => sym(u.choose(&HEADS)?),
                "MatchQ" => pattern(u, depth - 1, true)?,
                "Unapply" => sym(u.choose(&ATOMS)?),
                _ => apply(
                    "Rule",
                    &[pattern(u, depth - 1, true)?, expr(u, depth - 1, room)?],
                ),
            };
            match head {
                "Map" | "Apply" => apply(head, &[second, subject]),
                _ => apply(head, &[subject, second]),
            }
        }
        10 if room.applying => {
            let body = expr(u, depth - 1, room)?;
            let function = if u.ratio(1, 3)? {
                apply("Unapply", &[body, sym(u.choose(&["x", "y", "u"])?)])
            } else {
                apply("Function", &[sym("n"), body])
            };
            Expr::apply(&function, &args(u, room)?)
        }
        // A number applied as a function, which warns.
        11 if u.ratio(1, 4)? => {
            let number = Expr::from(i64::from(u.int_in_range(-3i8..=9)?));
            Expr::apply(&number, &args(u, room)?)
        }
        _ => apply("List", &args(u, room)?),
    })
}

/// A pattern: blanks, named or with a head, sequences, defaults,
/// alternatives, a literal, and with `tests` tests and conditions without
/// calls. Not under a `Flat` head: there a blank may bind a new group
/// `k[a, b]`, whose test evaluates it, which matches the rule again with a
/// budget of its own (R8), and so on.
fn pattern(u: &mut Unstructured<'_>, depth: u8, tests: bool) -> Result<Expr> {
    let blank = |u: &mut Unstructured<'_>| -> Result<Expr> {
        let kind = *u.choose(&["Blank", "BlankSequence", "BlankNullSequence"])?;
        Ok(if u.ratio(1, 3)? {
            // Heads of applications, and of kinds of atoms (CAS-150).
            let heads = [
                "List", "f", "Plus", "Integer", "Rational", "Real", "Complex", "Symbol",
            ];
            apply(kind, &[sym(u.choose(&heads)?)])
        } else {
            apply(kind, &[])
        })
    };
    let mut no_calls = TESTS;
    Ok(match u.int_in_range(0u8..=9)? {
        0..=3 => apply("Pattern", &[sym(u.choose(&VARIABLES)?), blank(u)?]),
        4 => blank(u)?,
        5 => apply(
            "Optional",
            &[apply(
                "Pattern",
                &[sym(u.choose(&VARIABLES)?), apply("Blank", &[])],
            )],
        ),
        6 if depth > 0 => apply(
            "Alternatives",
            &[pattern(u, depth - 1, tests)?, pattern(u, depth - 1, tests)?],
        ),
        7 if depth > 0 && tests => apply(
            "Condition",
            &[
                pattern(u, depth - 1, tests)?,
                expr(u, depth - 1, &mut no_calls)?,
            ],
        ),
        8 if tests => apply(
            "PatternTest",
            &[
                apply(
                    "Pattern",
                    &[sym(u.choose(&VARIABLES)?), apply("Blank", &[])],
                ),
                sym(u.choose(&["IntegerQ", "ListQ", "NumberQ"])?),
            ],
        ),
        _ => expr(u, depth.min(1), &mut no_calls)?,
    })
}

/// A context with random attributes, defaults, rules, domain and
/// assumptions, and small limits, `recursion` levels deep, evaluating on
/// the caller's thread.
fn context(u: &mut Unstructured<'_>, recursion: u32) -> Result<Context> {
    let mut ctx = Context::new();
    // First: definitions and assumptions evaluate too.
    let mut limits: Limits = ctx.limits();
    limits.recursion = recursion;
    limits.iterations = u.int_in_range(16..=128)?;
    limits.match_steps = u.int_in_range(1_000..=10_000)?;
    limits.number_bits = 1 << u.int_in_range(6u32..=14)?;
    limits.terms = 1 << u.int_in_range(4u32..=12)?;
    limits.stack = 0;
    ctx.set_limits(limits);
    let mut attributes = [Attributes::NONE; HEADS.len()];
    for a in &mut attributes {
        for one in [
            Attributes::FLAT,
            Attributes::ORDERLESS,
            Attributes::ONE_IDENTITY,
            Attributes::LISTABLE,
            Attributes::HOLD_ALL,
            Attributes::HOLD_FIRST,
            Attributes::HOLD_REST,
        ] {
            if u.ratio(1, 6)? {
                *a = a.union(one);
            }
        }
    }
    let flat = |head: &str| {
        matches!(head, "Plus" | "Times")
            || HEADS
                .iter()
                .zip(&attributes)
                .any(|(h, a)| *h == head && a.contains(Attributes::FLAT))
    };
    let set_attributes = |ctx: &mut Context| {
        for (head, a) in HEADS.iter().zip(attributes) {
            let _ = ctx.set_attributes(Symbol::new(head), a);
        }
    };
    // Attributes before the rules or after them, which recompiles them.
    let before = u.arbitrary::<bool>()?;
    if before {
        set_attributes(&mut ctx);
    }
    if u.ratio(1, 4)? {
        let _ = ctx.set_default(Symbol::new(u.choose(&HEADS)?), leaf(u, &mut { TESTS })?);
    }
    for _ in 0..u.int_in_range(0u8..=5)? {
        // A down-rule, a sub-rule `s[…][…]` or an up-rule.
        let kind = u.int_in_range(0u8..=4)?;
        let head = match kind {
            1 => "s",
            2 => *u.choose(&["g", "Plus", "Times", "Power", "Abs"])?,
            _ => *u.choose(&HEADS)?,
        };
        let tests = !flat(head);
        let count = u.int_in_range(0u8..=3)?;
        let mut args: Vec<Expr> = (0..count)
            .map(|_| pattern(u, 2, tests))
            .collect::<Result<_>>()?;
        let mut room = Room {
            calls: Some(1),
            applying: false,
            variables: true,
            used: Some([false; VARIABLES.len()]),
        };
        let rhs = expr(u, 3, &mut room)?;
        let rhs = if tests && u.ratio(1, 4)? {
            apply("Condition", &[rhs, expr(u, 2, &mut { TESTS })?])
        } else {
            rhs
        };
        let _ = match kind {
            0 => ctx.define(&apply(head, &args), &rhs),
            1 => {
                let inner = Expr::apply(&sym("s"), &args);
                let arg = pattern(u, 1, tests)?;
                ctx.define_delayed(&Expr::apply(&inner, &[arg]), &rhs)
            }
            2 => {
                let tag = *u.choose(&["a", "b", "u"])?;
                args.push(sym(tag));
                ctx.define_up_delayed(Symbol::new(tag), &apply(head, &args), &rhs)
            }
            _ => ctx.define_delayed(&apply(head, &args), &rhs),
        };
    }
    if !before {
        set_attributes(&mut ctx);
    }
    if u.ratio(1, 5)? {
        ctx.set_domain(Domain::Complex);
    }
    for _ in 0..u.int_in_range(0u8..=2)? {
        let assumption = match u.int_in_range(0u8..=4)? {
            0 => "Greater[x, 0]",
            1 => "Element[a, Integers]",
            2 => "Less[b, x]",
            3 => "Unequal[y, 0]",
            _ => "And[Positive[u], Greater[u, b]]",
        };
        ctx = ctx.assuming(Expr::parse_full_form(assumption).expect("an assumption"));
    }
    Ok(ctx)
}

/// Replays `steps`, located in `start`, checking each evaluation among
/// their children, as the trace tests do: the expression they make.
fn replay<'t>(trace: &'t Trace, start: &Expr, steps: impl Iterator<Item = &'t Step>) -> Expr {
    let mut e = start.clone();
    for s in steps {
        match s.kind {
            StepKind::Algorithm(_) => continue,
            StepKind::Evaluation => panic!("an evaluation located in {e:?}: {s:?}"),
            _ => {}
        }
        let part = e
            .part(&s.position)
            .unwrap_or_else(|| panic!("no part {s:?} in {e:?}"));
        assert_eq!(*part, s.before, "at {:?} in {e:?}", s.position);
        e = e
            .replace_part(&s.position, s.after.clone())
            .expect("the part exists");
        for c in trace.children(s) {
            if c.kind == StepKind::Evaluation {
                assert_eq!(replay(trace, &c.before, trace.children(c)), c.after);
            }
        }
    }
    e
}

/// The tree a trace describes: ids are places, parents and children agree.
fn check_tree(trace: &Trace) {
    for (i, s) in trace.steps.iter().enumerate() {
        assert_eq!(s.id, i);
        for &c in &s.children {
            assert!(
                c < i && trace.steps[c].parent == Some(i),
                "child {c} of {i}"
            );
        }
        if let StepKind::Memo { first: Some(f) } = s.kind {
            assert!(f < i, "memo {i} points at {f}");
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(recursion) = u.int_in_range(16u32..=128) else {
        return;
    };
    // A thread with just the stack the levels may take, so that a path
    // through a level larger than `STACK_PER_LEVEL` overflows it.
    let stack = RESERVE + recursion as usize * STACK_PER_LEVEL;
    std::thread::scope(|scope| {
        let thread = std::thread::Builder::new()
            .stack_size(stack)
            .spawn_scoped(scope, move || check(u, recursion))
            .expect("the thread starts");
        if let Err(panic) = thread.join() {
            std::panic::resume_unwind(panic);
        }
    });
});

/// The checks on one input.
fn check(mut u: Unstructured<'_>, recursion: u32) {
    let u = &mut u;
    let Ok(ctx) = context(u, recursion) else {
        return;
    };
    let mut room = Room {
        calls: None,
        applying: true,
        variables: false,
        used: None,
    };
    let Ok(e) = expr(u, 5, &mut room) else {
        return;
    };
    let plain = ctx.eval(&e);
    if let Ok(r) = &plain {
        assert_eq!(ctx.eval(r).as_ref(), Ok(r), "{e:?} gave {r:?}");
    }
    let (traced, report) = ctx.eval_with(
        &e,
        &EvalOptions {
            trace: TraceLevel::Verbose,
        },
    );
    assert_eq!(
        traced.is_ok(),
        plain.is_ok(),
        "{e:?}: {traced:?} and {plain:?}"
    );
    check_tree(&report.trace);
    if let (Ok(r), Ok(t)) = (&plain, &traced) {
        assert_eq!(t, r, "{e:?}");
        assert_eq!(&replay(&report.trace, &e, report.trace.top()), r, "{e:?}");
    }
    for w in &report.warnings {
        let _ = w.to_string();
    }
    if report.trace.steps.len() <= 64 {
        let _ = report.trace.verify(&ctx);
    }
}
