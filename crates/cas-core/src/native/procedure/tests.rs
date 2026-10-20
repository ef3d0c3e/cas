//! Tests of procedures in full form; programs in the language are in
//! cas-lang's session scripts.

use proptest::prelude::*;

use crate::{Context, EvalError, Expr, Interrupt, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn statements() {
    check(&[
        // The value of the last statement; of an assignment, the value.
        ("Block[List[x], Set[x, 2], Times[x, x]]", "4"),
        ("Block[List[x], Set[x, Plus[1, 1]]]", "2"),
        ("Block[List[]]", "Null"),
        (
            "Block[List[x, y], Set[List[x, y], List[1, 2]], List[y, x]]",
            "List[2, 1]",
        ),
        // A name that is not a local stays a global.
        ("Block[List[x], Set[x, a], Plus[x, b]]", "Plus[a, b]"),
        // Branches in the same frame.
        (
            "Block[List[x], Set[x, 1], If[Greater[x, 0], Statements[Set[x, 5]], Statements[Set[x, 6]]], x]",
            "5",
        ),
        (
            "Block[List[x], Set[x, 1], If[Less[x, 0], Statements[Set[x, 5]]], x]",
            "1",
        ),
        // Loops: a range, a list, a set, an empty range.
        (
            "Block[List[s, k], Set[s, 0], For[k, Range[1, 4], Statements[Set[s, Plus[s, k]]]], s]",
            "10",
        ),
        (
            "Block[List[s, k], Set[s, 1], For[k, List[2, 3], Statements[Set[s, Times[s, k]]]], s]",
            "6",
        ),
        (
            "Block[List[s, k], Set[s, 0], For[k, FiniteSet[5, 5, 1], Statements[Set[s, Plus[s, k]]]], s]",
            "6",
        ),
        (
            "Block[List[s, k], Set[s, 0], For[k, Range[3, 1], Statements[Set[s, 1]]], s]",
            "0",
        ),
        (
            "Block[List[m, n], Set[m, 10], Set[n, 0], While[Greater[m, 1], Statements[Set[m, Times[1/2, m]], Set[n, Plus[n, 1]]]], List[m, n]]",
            "List[5/8, 4]",
        ),
        // Jumps.
        (
            "Block[List[k], For[k, Range[1, 10], Statements[If[Greater[k, 3], Statements[Return[k]]]]], 0]",
            "4",
        ),
        (
            "Block[List[s, k], Set[s, 0], For[k, Range[1, 10], Statements[If[Greater[k, 3], Statements[Break[]]], Set[s, Plus[s, k]]]], s]",
            "6",
        ),
        (
            "Block[List[s, k], Set[s, 0], For[k, Range[1, 5], Statements[If[EvenQ[k], Statements[Continue[]]], Set[s, Plus[s, k]]]], s]",
            "9",
        ),
        // `break` ends the loop: later rounds would add more.
        (
            "Block[List[s, k], Set[s, 0], For[k, Range[1, 10], Statements[Set[s, Plus[s, k]], If[Greater[s, 5], Statements[Break[]]]]], s]",
            "6",
        ),
        ("Block[List[], Return[]]", "Null"),
        // A loop's value is `Null`.
        ("Block[List[k], For[k, Range[1, 2], Statements[k]]]", "Null"),
    ]);
}

/// The values of locals go into the statements before they are evaluated,
/// except where a name is bound again.
#[test]
fn substitution_respects_binders() {
    check(&[
        // A lambda captures the value when it is made.
        (
            "Block[List[y, g], Set[y, 3], Set[g, Function[List[t], Plus[t, y]]], Set[y, 100], g[1]]",
            "4",
        ),
        // Its own parameter is not replaced.
        (
            "Block[List[t], Set[t, 5], Function[List[t], Times[2, t]][1]]",
            "2",
        ),
        // A nested block's locals are its own.
        (
            "Block[List[t], Set[t, 5], Block[List[t], Set[t, 1], Plus[t, 1]]]",
            "2",
        ),
        (
            "Block[List[a], Set[a, 2], Block[List[c], Set[c, Plus[a, 1]], Times[c, 2]]]",
            "6",
        ),
        // A pattern variable in a rule is not the local of the same name.
        (
            "Block[List[x], Set[x, 7], ReplaceAll[f[1, 2], RuleDelayed[f[Pattern[x, Blank[]], y], g[x]]]]",
            "f[1, 2]",
        ),
        (
            "Block[List[x], Set[x, 7], ReplaceAll[f[1, y], RuleDelayed[f[Pattern[x, Blank[]], y], g[x]]]]",
            "g[1]",
        ),
        // A literal name in a rule takes the local's value.
        (
            "Block[List[x], Set[x, 7], ReplaceAll[f[7], Rule[x, 0]]]",
            "f[0]",
        ),
    ]);
}

#[test]
fn errors() {
    let ctx = Context::new();
    for (input, error) in [
        (
            "Block[List[], If[Greater[z, 0], Statements[1]]]",
            EvalError::UndecidedCondition,
        ),
        (
            "Block[List[], While[z, Statements[]]]",
            EvalError::UndecidedCondition,
        ),
        (
            "Block[List[k], For[k, z, Statements[]]]",
            EvalError::NotIterable,
        ),
        (
            "Block[List[k], For[k, Range[1, n], Statements[]]]",
            EvalError::NotIterable,
        ),
        (
            "Block[List[k], For[k, Range[1/2, 3], Statements[]]]",
            EvalError::NotIterable,
        ),
        (
            "Block[List[a, b], Set[List[a, b], List[1, 2, 3]]]",
            EvalError::Shape,
        ),
        // Only locals are written.
        ("Block[List[], Set[g, 1]]", EvalError::Shape),
    ] {
        assert_eq!(ctx.eval(&read(input)), Err(error), "{input}");
    }
    // What is not a block stays.
    for input in ["Block[x, 1]", "Block[List[1], 1]", "Block[]"] {
        assert_eq!(ctx.eval(&read(input)), Ok(read(input)), "{input}");
    }
}

/// A block takes a level of `Limits::recursion`, and a group of its
/// statements another one, as the stack test counts them: here three per
/// call, with the evaluation of the call.
#[test]
fn blocks_take_levels() {
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.recursion = 600;
    ctx.set_limits(limits);
    ctx.define_delayed(&read("p[0]"), &read("done")).unwrap();
    ctx.define_delayed(
        &read("p[Pattern[n, Blank[]]]"),
        &read("Block[List[r], While[True, Statements[Set[r, p[Plus[n, -1]]], Break[]]], r]"),
    )
    .unwrap();
    assert_eq!(ctx.eval(&read("p[150]")), Ok(read("done")));
    assert_eq!(ctx.eval(&read("p[250]")), Err(EvalError::RecursionLimit));
}

/// A loop without a rewrite still stops at an interrupt.
#[test]
fn an_endless_loop_stops_at_an_interrupt() {
    let mut ctx = Context::new();
    let interrupt = Interrupt::new();
    ctx.set_interrupt(interrupt.clone());
    interrupt.interrupt();
    for input in [
        "Block[List[], While[True, Statements[]]]",
        "Block[List[k], For[k, Range[1, 1000000000000], Statements[]]]",
    ] {
        assert_eq!(
            ctx.eval(&read(input)),
            Err(EvalError::Interrupted),
            "{input}"
        );
    }
}

/// Expressions over `a`, `b`, small integers and the names `locals`.
fn expr(locals: &'static [&'static str]) -> impl Strategy<Value = Expr> {
    let names: Vec<&str> = ["a", "b"].iter().chain(locals).copied().collect();
    let leaf = prop_oneof![
        prop::sample::select(names).prop_map(|n| Expr::from(Symbol::new(n))),
        (-3i64..4).prop_map(Expr::from),
    ];
    leaf.prop_recursive(3, 16, 3, |inner| {
        (
            prop::sample::select(vec!["Plus", "Times", "List", "f", "Power", "Hold"]),
            prop::collection::vec(inner, 0..3),
        )
            .prop_map(|(head, args)| Expr::apply(&Symbol::new(head).into(), &args))
    })
}

/// A block that assigns `x` and then `y`, each from what comes before it,
/// with a loop that adds to `y` sometimes, and ends with an expression of
/// both.
fn block() -> impl Strategy<Value = Expr> {
    (
        expr(&[]),
        expr(&["x"]),
        expr(&["x", "y", "k"]),
        expr(&["x", "y"]),
        any::<bool>(),
    )
        .prop_map(|(e1, e2, step, last, looping)| {
            let s = |name: &str| Expr::from(Symbol::new(name));
            let apply = |head: Symbol, args: &[Expr]| Expr::apply(&head.into(), args);
            let mut statements = vec![
                apply(Symbol::LIST, &[s("x"), s("y"), s("k")]),
                apply(Symbol::SET, &[s("x"), e1]),
                apply(Symbol::SET, &[s("y"), e2]),
            ];
            if looping {
                let body = apply(
                    Symbol::STATEMENTS,
                    &[apply(
                        Symbol::SET,
                        &[s("y"), apply(Symbol::PLUS, &[s("y"), step])],
                    )],
                );
                let range = apply(Symbol::RANGE, &[1.into(), 3.into()]);
                statements.push(apply(Symbol::FOR, &[s("k"), range, body]));
            }
            statements.push(last);
            apply(Symbol::BLOCK, &statements)
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    /// A procedure's result evaluates to itself (CAS-53).
    #[test]
    fn procedures_are_idempotent(e in block()) {
        let ctx = Context::new();
        if let Ok(r) = ctx.eval(&e) {
            prop_assert!(ctx.eval(&r) == Ok(r.clone()), "{:?} gave {:?}", e, r);
        }
    }

    /// Stamps and the memo change no result of a procedure (CAS-53).
    #[test]
    fn caches_change_nothing_in_procedures(e in block()) {
        let ctx = Context::new();
        let uncached = ctx.eval_uncached(&e);
        prop_assert!(ctx.eval(&e) == uncached, "{:?}", e);
        prop_assert!(ctx.eval(&e) == uncached, "{:?}", e);
    }
}
