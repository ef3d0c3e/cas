//! Tests of the evaluation loop. Native functions for tests live in a global
//! table, so each test gives them symbols of its own: tests run in parallel.

use std::collections::HashMap;
use std::hint::black_box;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{LazyLock, Mutex, PoisonError};

use proptest::prelude::*;

use super::{Eval, NativeFn, STACK_PER_LEVEL};
use crate::{Attributes, Context, EvalError, Expr, Limits, Symbol};

static NATIVES: LazyLock<Mutex<HashMap<Symbol, NativeFn>>> = LazyLock::new(Mutex::default);

/// The native function a test gave `s`, which the evaluator consults in
/// tests.
pub(super) fn native(s: Symbol) -> Option<NativeFn> {
    NATIVES
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .get(&s)
        .copied()
}

/// Gives the symbol `name` the native function `f`.
pub(crate) fn define_native(name: &str, f: NativeFn) -> Symbol {
    let s = Symbol::new(name);
    NATIVES
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .insert(s, f);
    s
}

fn sym(name: &str) -> Symbol {
    Symbol::new(name)
}

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

fn eval_err(ctx: &Context, input: &str) -> EvalError {
    ctx.eval(&read(input)).expect_err(input)
}

/// `f[f[…f[leaf]…]]`, `depth` deep.
fn nested(f: Symbol, leaf: Expr, depth: usize) -> Expr {
    (0..depth).fold(leaf, |e, _| Expr::apply(&f.into(), &[e]))
}

#[test]
fn atoms_and_values() {
    let mut ctx = Context::new();
    check(
        &ctx,
        &[
            ("42", "42"),
            ("-3/2", "-3/2"),
            ("Complex[0, 1]", "Complex[0, 1]"),
            (r#""s""#, r#""s""#),
            ("atoms_x", "atoms_x"),
            ("Plus", "Plus"),
        ],
    );
    ctx.set(sym("atoms_x"), read("atoms_y")).unwrap();
    ctx.set(sym("atoms_y"), read("atoms_f[atoms_z, 2]"))
        .unwrap();
    ctx.set(sym("atoms_z"), read("3")).unwrap();
    check(
        &ctx,
        &[
            ("atoms_x", "atoms_f[3, 2]"),
            ("List[atoms_z, atoms_w]", "List[3, atoms_w]"),
            // A head is evaluated too.
            ("atoms_z[1]", "3[1]"),
        ],
    );
}

#[test]
fn held_arguments_are_left_alone() {
    let mut ctx = Context::new();
    let x = sym("holds_x");
    ctx.set(x, read("1")).unwrap();
    for (name, attributes) in [
        ("holds_first", Attributes::HOLD_FIRST),
        ("holds_rest", Attributes::HOLD_REST),
        ("holds_all", Attributes::HOLD_ALL),
    ] {
        ctx.set_attributes(sym(name), attributes).unwrap();
    }
    check(
        &ctx,
        &[
            ("holds_first[holds_x, holds_x]", "holds_first[holds_x, 1]"),
            (
                "holds_rest[holds_x, holds_x, holds_x]",
                "holds_rest[1, holds_x, holds_x]",
            ),
            ("holds_all[holds_x, holds_x]", "holds_all[holds_x, holds_x]"),
            ("Hold[holds_x]", "Hold[holds_x]"),
            ("If[holds_x, holds_x, holds_x]", "If[1, holds_x, holds_x]"),
            (
                "Function[List[holds_x], holds_x]",
                "Function[List[holds_x], holds_x]",
            ),
        ],
    );
    // A held part is not evaluated, so not stamped.
    let part = read("holds_g[holds_x]");
    ctx.eval(&Expr::apply(
        &Symbol::HOLD.into(),
        std::slice::from_ref(&part),
    ))
    .unwrap();
    assert!(!ctx.is_fixpoint(&part));
}

#[test]
fn statements_are_not_run() {
    let mut ctx = Context::new();
    let x = sym("pure_x");
    ctx.set(x, read("1")).unwrap();
    let id = ctx.id();
    check(
        &ctx,
        &[
            ("Set[pure_y, pure_x]", "Set[pure_y, 1]"),
            (
                "SetDelayed[pure_f[pure_x], pure_x]",
                "SetDelayed[pure_f[pure_x], pure_x]",
            ),
            ("Clear[pure_x]", "Clear[pure_x]"),
        ],
    );
    assert_eq!(ctx.id(), id);
    assert_eq!(ctx.value(sym("pure_y")), None);
    assert_eq!(ctx.value(x), Some(&read("1")));
}

#[test]
fn sequences_are_spliced() {
    let mut ctx = Context::new();
    ctx.set(sym("seq_x"), read("2")).unwrap();
    check(
        &ctx,
        &[
            ("seq_f[a, Sequence[b, c], d]", "seq_f[a, b, c, d]"),
            ("seq_f[Sequence[]]", "seq_f[]"),
            ("seq_f[Sequence[a, Sequence[b]], c]", "seq_f[a, b, c]"),
            ("seq_f[Sequence[seq_x]]", "seq_f[2]"),
            // Spliced in held positions too, unevaluated.
            ("Hold[Sequence[seq_x, a], b]", "Hold[seq_x, a, b]"),
            ("Hold[Sequence[Sequence[seq_x]]]", "Hold[seq_x]"),
            ("Sequence[Sequence[a], b]", "Sequence[a, b]"),
        ],
    );
}

#[test]
fn flat_heads_flatten() {
    let mut ctx = Context::new();
    ctx.set_attributes(sym("flat_f"), Attributes::FLAT).unwrap();
    check(
        &ctx,
        &[
            ("flat_f[a, flat_f[b, flat_f[c]], d]", "flat_f[a, b, c, d]"),
            ("flat_f[flat_f[]]", "flat_f[]"),
            // Held arguments flatten at any depth in one step.
            ("And[p, And[q, And[r, s]]]", "And[p, q, r, s]"),
            ("Plus[a, Plus[b, c]]", "Plus[a, b, c]"),
            ("flat_g[flat_g[a]]", "flat_g[flat_g[a]]"),
        ],
    );
}

#[test]
fn listable_heads_thread_over_lists() {
    let mut ctx = Context::new();
    ctx.set_attributes(sym("list_g"), Attributes::LISTABLE)
        .unwrap();
    ctx.set(sym("list_v"), read("List[1, 2]")).unwrap();
    check(
        &ctx,
        &[
            ("list_g[List[1, 2]]", "List[list_g[1], list_g[2]]"),
            (
                "list_g[List[1, 2], List[a, b]]",
                "List[list_g[1, a], list_g[2, b]]",
            ),
            ("list_g[List[]]", "List[]"),
            (
                "list_g[List[List[1, 2], 3]]",
                "List[List[list_g[1], list_g[2]], list_g[3]]",
            ),
            ("list_g[list_v]", "List[list_g[1], list_g[2]]"),
            ("list_g[a]", "list_g[a]"),
            (
                "Plus[List[1, 2], List[a, b]]",
                "List[Plus[1, a], Plus[2, b]]",
            ),
            // `*` is linear algebra: `Times` and `Power` do not thread (R13).
            (
                "Times[List[1, 2], List[3, 4]]",
                "Times[List[1, 2], List[3, 4]]",
            ),
            ("Power[List[1, 2], 2]", "Power[List[1, 2], 2]"),
        ],
    );
    for input in [
        "list_g[List[1], List[2, 3]]",
        "list_g[List[1, 2], 3]",
        // A scalar plus a list is an error (O8).
        "Plus[List[1, 2], 1]",
        "Plus[List[1, 2], List[1]]",
    ] {
        assert_eq!(eval_err(&ctx, input), EvalError::Shape, "{input}");
    }
}

#[test]
fn orderless_heads_sort() {
    let mut ctx = Context::new();
    ctx.set_attributes(sym("sort_f"), Attributes::ORDERLESS)
        .unwrap();
    check(
        &ctx,
        &[
            ("sort_f[c, b, a]", "sort_f[a, b, c]"),
            ("sort_f[x, 2, 1/2]", "sort_f[1/2, 2, x]"),
            ("Plus[c, Plus[b, a], 1]", "Plus[1, a, b, c]"),
            // In `Times`, lists keep their order, after the other factors.
            (
                "Times[List[2], b, List[1], a]",
                "Times[a, b, List[2], List[1]]",
            ),
            ("Times[List[2], List[1]]", "Times[List[2], List[1]]"),
        ],
    );
}

#[test]
fn natives_rewrite_and_their_results_are_evaluated() {
    // `peel[h[x]]` is `x`.
    define_native("native_peel", |_, e| {
        Ok(match e.args() {
            [a] if !a.is_atom() => a.args().first().cloned(),
            _ => None,
        })
    });
    // `same[x]` returns its input: no change, no loop.
    define_native("native_same", |_, e| Ok(Some(e.clone())));
    // `twice[x]` evaluates `peel[x]` itself, a level deeper.
    define_native("native_twice", |ev, e| {
        let inner = Expr::apply(&sym("native_peel").into(), e.args());
        ev.eval(&inner).map(Some)
    });
    let mut ctx = Context::new();
    ctx.set(sym("native_x"), read("native_h[native_k[7]]"))
        .unwrap();
    check(
        &ctx,
        &[
            ("native_peel[native_h[a]]", "a"),
            ("native_peel[native_x]", "native_k[7]"),
            ("native_peel[native_peel[native_x]]", "7"),
            ("native_peel[a]", "native_peel[a]"),
            ("native_same[a]", "native_same[a]"),
            ("native_twice[native_x]", "native_k[7]"),
        ],
    );
}

#[test]
fn the_memo_evaluates_a_repeated_part_once_per_call() {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    define_native("memo_counted", |_, e| {
        CALLS.fetch_add(1, Ordering::Relaxed);
        Ok(Some(Expr::apply(&sym("memo_done").into(), e.args())))
    });
    let ctx = Context::new();
    let e = read("memo_f[memo_counted[1], memo_counted[1], memo_g[memo_counted[1]]]");
    let r = ctx.eval(&e).unwrap();
    assert_eq!(
        format!("{r:?}"),
        "memo_f[memo_done[1], memo_done[1], memo_g[memo_done[1]]]"
    );
    assert_eq!(CALLS.load(Ordering::Relaxed), 1);
    // The memo lasts one call; the input was not a fixpoint.
    ctx.eval(&e).unwrap();
    assert_eq!(CALLS.load(Ordering::Relaxed), 2);
    // Without caches, every occurrence is evaluated.
    ctx.eval_uncached(&e).unwrap();
    assert_eq!(CALLS.load(Ordering::Relaxed), 5);
}

#[test]
fn stamps_make_evaluating_a_result_again_free() {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    // Applies to nothing, but counts its calls.
    define_native("stamp_inert", |_, _| {
        CALLS.fetch_add(1, Ordering::Relaxed);
        Ok(None)
    });
    let mut ctx = Context::new();
    let e = read("stamp_f[stamp_inert[1], stamp_inert[2]]");
    let r = ctx.eval(&e).unwrap();
    assert_eq!(r, e);
    assert!(ctx.is_fixpoint(&r));
    assert_eq!(CALLS.load(Ordering::Relaxed), 2);
    assert_eq!(ctx.eval(&r).unwrap(), r);
    assert_eq!(
        CALLS.load(Ordering::Relaxed),
        2,
        "a stamp hit calls nothing"
    );
    // Another state: the stamps miss.
    ctx.set(sym("stamp_x"), read("1")).unwrap();
    assert_eq!(ctx.eval(&r).unwrap(), r);
    assert_eq!(CALLS.load(Ordering::Relaxed), 4);
}

#[test]
fn iteration_limit() {
    let mut ctx = Context::new();
    ctx.set(sym("iter_x"), read("iter_y")).unwrap();
    ctx.set(sym("iter_y"), read("iter_x")).unwrap();
    assert_eq!(eval_err(&ctx, "iter_x"), EvalError::IterationLimit);
    define_native("iter_ping", |_, e| {
        Ok(Some(Expr::apply(&sym("iter_pong").into(), e.args())))
    });
    define_native("iter_pong", |_, e| {
        Ok(Some(Expr::apply(&sym("iter_ping").into(), e.args())))
    });
    assert_eq!(eval_err(&ctx, "iter_ping[1]"), EvalError::IterationLimit);
    // A chain of values within the limit.
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.iterations = 9;
    ctx.set_limits(limits);
    for i in 0..10 {
        ctx.set(
            sym(&format!("iter_v{i}")),
            read(&format!("iter_v{}", i + 1)),
        )
        .unwrap();
    }
    assert_eq!(eval_err(&ctx, "iter_v0"), EvalError::IterationLimit);
    check(&ctx, &[("iter_v1", "iter_v10")]);
}

#[test]
fn recursion_limit_on_the_callers_thread() {
    let mut ctx = Context::new();
    let f = sym("rec_f");
    let limit = ctx.limits().recursion as usize;
    // Each level of the expression is one level of evaluation.
    let x = Expr::from(sym("rec_x"));
    assert!(ctx.eval(&nested(f, x.clone(), limit - 1)).is_ok());
    let deep = nested(f, x, limit);
    assert_eq!(ctx.eval(&deep), Err(EvalError::RecursionLimit));
    // Far deeper than the stack would allow without the limit.
    let deeper = nested(f, Expr::from(1), 100_000);
    assert_eq!(ctx.eval(&deeper), Err(EvalError::RecursionLimit));
    // A value that contains itself.
    ctx.set(sym("rec_x"), read("rec_f[rec_x]")).unwrap();
    assert_eq!(eval_err(&ctx, "rec_x"), EvalError::RecursionLimit);
}

#[test]
fn deep_evaluation_runs_on_a_thread_of_its_own() {
    let mut ctx = Context::new();
    ctx.set_limits(Limits::deep());
    let levels = Limits::deep().recursion as usize;
    assert!(levels >= 80_000, "{levels}");
    let f = sym("deep_f");
    let e = nested(f, Expr::from(1), levels - 1);
    // Compared without `assert_eq!`, which would print 10⁵ levels.
    let same = ctx.eval(&e).unwrap() == e;
    assert!(same, "a deep fixpoint evaluates to itself");
    assert_eq!(
        ctx.eval(&nested(f, Expr::from(2), levels)),
        Err(EvalError::RecursionLimit)
    );
    // A recursion limit beyond what the stack holds is capped by the stack.
    let mut limits = Limits::deep();
    limits.recursion = u32::MAX;
    ctx.set_limits(limits);
    assert_eq!(
        ctx.eval(&nested(f, Expr::from(3), levels + 10)),
        Err(EvalError::RecursionLimit)
    );
}

/// A way to nest one more level around an expression.
type Shape = fn(Symbol, Expr) -> Expr;

/// `f[e]`, as a layer of a pattern or of its subject.
fn matcher_layer(f: Symbol, e: Expr) -> Expr {
    Expr::apply(&f.into(), &[e])
}

/// The stack each level of evaluation takes, at most, through arguments,
/// through heads, through a native that evaluates its arguments (`And` and
/// `Or`, which alternate so they do not flatten), through the matcher, and
/// through rules whose conditions recurse: measured from the address of a
/// local in a native at the bottom.
#[test]
fn a_level_takes_less_than_stack_per_level() {
    static BOTTOM: AtomicUsize = AtomicUsize::new(usize::MAX);
    let probe = define_native("stack_probe", |_, _| {
        let here = black_box(0u8);
        BOTTOM.fetch_min(std::ptr::addr_of!(here) as usize, Ordering::Relaxed);
        Ok(None)
    });
    let f = sym("stack_f");
    // Each with the levels one layer takes.
    let shapes: [(&str, Shape, usize); 5] = [
        ("arguments", |f, e| Expr::apply(&f.into(), &[e]), 1),
        ("heads", |_, e| Expr::apply(&e, &[]), 1),
        (
            "natives",
            |_, e| {
                let or = Expr::apply(&Symbol::OR.into(), &[e]);
                Expr::apply(&Symbol::AND.into(), &[or])
            },
            2,
        ),
        // `MatchQ[f[…f[a]…], f[…f[Condition[_, probe[]]]…]]`: the matcher
        // nests one node per level, and its condition evaluates the probe.
        ("matcher", |f, e| matcher_layer(f, e), 1),
        // `chain[n]`, whose rule's condition evaluates `chain[n - 1]`: the
        // rule's evaluation, its four match nodes and the condition's
        // evaluation per step.
        ("rules", |_, e| e, 6),
    ];
    let chain = sym("stack_chain");
    let depth = 2000;
    let measured = std::thread::Builder::new()
        .stack_size(64 << 20)
        .spawn(move || {
            let mut ctx = Context::new();
            let mut limits = ctx.limits();
            limits.recursion = 10 * depth;
            ctx.set_limits(limits);
            let rules = [
                ("stack_chain[0]", "stack_probe[]"),
                (
                    "stack_chain[Pattern[n, Blank[]]]",
                    "Condition[done, SameQ[stack_chain[Plus[n, -1]], done]]",
                ),
            ];
            for (lhs, rhs) in rules {
                ctx.define_delayed(&read(lhs), &read(rhs)).unwrap();
            }
            shapes.map(|(name, shape, levels)| {
                let e = if name == "rules" {
                    Expr::apply(&chain.into(), &[i64::from(depth).into()])
                } else if name == "matcher" {
                    let condition = Expr::apply(
                        &Symbol::CONDITION.into(),
                        &[
                            Expr::apply(&Symbol::BLANK.into(), &[]),
                            Expr::apply(&probe.into(), &[]),
                        ],
                    );
                    let pattern = (0..depth).fold(condition, |e, _| shape(f, e));
                    let subject = (0..depth).fold(Expr::from(sym("stack_a")), |e, _| shape(f, e));
                    let pattern = Expr::apply(&Symbol::HOLD_PATTERN.into(), &[pattern]);
                    Expr::apply(&Symbol::MATCH_Q.into(), &[subject, pattern])
                } else {
                    (0..depth).fold(Expr::apply(&probe.into(), &[]), |e, _| shape(f, e))
                };
                BOTTOM.store(usize::MAX, Ordering::Relaxed);
                let top = black_box(0u8);
                let top = std::ptr::addr_of!(top) as usize;
                ctx.eval_uncached(&e).unwrap();
                (
                    name,
                    (top - BOTTOM.load(Ordering::Relaxed)) / (levels * depth as usize),
                )
            })
        })
        .unwrap()
        .join()
        .unwrap();
    eprintln!("stack per level: {measured:?}");
    for (name, bytes) in measured {
        assert!(bytes <= STACK_PER_LEVEL, "{name}: {bytes} bytes per level");
    }
}

#[test]
fn eval_errors_print() {
    assert_eq!(EvalError::Shape.to_string(), "lists of different shapes");
    let ctx = Context::new();
    let ev = Eval::new(&ctx, 1);
    assert!(format!("{ev:?}").starts_with("Eval { context: StateId("));
}

/// The context of the property tests: symbols with each attribute, and
/// values that do not refer to themselves.
fn property_context() -> Context {
    let mut ctx = Context::new();
    for (name, attributes) in [
        ("p_flat", Attributes::FLAT),
        ("p_orderless", Attributes::ORDERLESS),
        ("p_listable", Attributes::LISTABLE),
        ("p_ac", Attributes::FLAT.union(Attributes::ORDERLESS)),
        ("p_first", Attributes::HOLD_FIRST),
        ("p_rest", Attributes::HOLD_REST),
        ("p_all", Attributes::HOLD_ALL),
    ] {
        ctx.set_attributes(sym(name), attributes).unwrap();
    }
    for (name, value) in [
        ("p_one", "1"),
        ("p_list", "List[1, p_one]"),
        ("p_seq", "Sequence[a, p_one]"),
        ("p_sum", "p_ac[p_one, b, p_ac[a]]"),
    ] {
        ctx.set(sym(name), read(value)).unwrap();
    }
    ctx
}

/// Random expressions over the symbols of [`property_context`] and some
/// builtins.
fn expr() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec![
            "a",
            "b",
            "p_one",
            "p_list",
            "p_seq",
            "p_sum",
            "p_flat",
            "p_listable",
        ])
        .prop_map(|n| Expr::from(sym(n))),
        (-3i64..3).prop_map(Expr::from),
    ];
    leaf.prop_recursive(4, 32, 4, |inner| {
        (
            prop::sample::select(vec![
                "p_flat",
                "p_orderless",
                "p_listable",
                "p_ac",
                "p_first",
                "p_rest",
                "p_all",
                "f",
                "Plus",
                "Times",
                "List",
                "Sequence",
                "Hold",
            ]),
            prop::collection::vec(inner, 0..4),
        )
            .prop_map(|(head, args)| Expr::apply(&sym(head).into(), &args))
    })
}

proptest! {
    /// Evaluating a result again gives it back.
    #[test]
    fn evaluation_is_idempotent(e in expr()) {
        let ctx = property_context();
        if let Ok(r) = ctx.eval(&e) {
            prop_assert!(ctx.eval(&r) == Ok(r.clone()), "{:?} gave {:?}", e, r);
            prop_assert!(ctx.eval_uncached(&r) == Ok(r.clone()), "{:?} gave {:?}", e, r);
        }
    }

    /// Stamps and the memo change no result.
    #[test]
    fn caches_change_nothing(e in expr()) {
        let ctx = property_context();
        let uncached = ctx.eval_uncached(&e);
        prop_assert!(ctx.eval(&e) == uncached, "{:?}", e);
        // Again, now that parts are stamped.
        prop_assert!(ctx.eval(&e) == uncached, "{:?}", e);
    }

    /// Permuting the arguments of an `Orderless` head changes nothing.
    #[test]
    fn argument_order_of_orderless_heads_changes_nothing(
        args in prop::collection::vec(expr(), 0..5),
        perm in any::<prop::sample::Index>(),
        // Not `Times`, where lists do not commute.
        head in prop::sample::select(vec!["p_orderless", "p_ac", "Plus"]),
    ) {
        let ctx = property_context();
        let head = Expr::from(sym(head));
        let mut shuffled = args.clone();
        let mid = perm.index(shuffled.len().max(1)) % shuffled.len().max(1);
        shuffled.rotate_left(mid);
        shuffled.reverse();
        let original = ctx.eval(&Expr::apply(&head, &args));
        prop_assert!(ctx.eval(&Expr::apply(&head, &shuffled)) == original);
    }
}
