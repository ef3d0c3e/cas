//! Tests of contexts: protection, scopes, state ids and fixpoint stamps.

use std::collections::HashSet;
use std::sync::Barrier;

use proptest::prelude::*;

use super::*;

fn sym(name: &str) -> Symbol {
    Symbol::new(name)
}

/// A node (atoms have no stamp) private to the calling test.
fn node(name: &str) -> Expr {
    Expr::apply(&sym(name).into(), &[1.into()])
}

#[test]
fn context_is_send_and_sync() {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<Context>();
}

#[test]
fn builtins_cannot_change() {
    let mut ctx = Context::new();
    let id = ctx.id();
    let err = Err(ContextError::Protected(Symbol::PLUS));
    assert_eq!(ctx.set_attributes(Symbol::PLUS, Attributes::NONE), err);
    assert_eq!(ctx.unprotect(Symbol::PLUS), err);
    assert_eq!(ctx.set(Symbol::PLUS, 1.into()), err);
    assert_eq!(ctx.unset(Symbol::PLUS), err);
    assert_eq!(ctx.id(), id, "a rejected change is no change");
    assert_eq!(
        Some(ctx.attributes(Symbol::PLUS)),
        Symbol::PLUS.builtin_attributes()
    );
    assert_eq!(
        ContextError::Protected(Symbol::PLUS).to_string(),
        "symbol Plus is protected"
    );
}

#[test]
fn user_attributes_and_protection() {
    use Attributes as A;
    let (f, x) = (sym("ctx_attr_f"), sym("ctx_attr_x"));
    let mut ctx = Context::new();
    assert_eq!(ctx.attributes(f), A::NONE);
    ctx.set_attributes(f, A::ORDERLESS).unwrap();
    assert_eq!(ctx.attributes(f), A::ORDERLESS);
    ctx.set_attributes(f, A::NONE).unwrap();
    assert_eq!(ctx.attributes(f), A::NONE);

    ctx.set_attributes(x, A::FLAT | A::PROTECTED).unwrap();
    assert_eq!(
        ctx.set_attributes(x, A::NONE),
        Err(ContextError::Protected(x))
    );
    assert_eq!(ctx.set(x, 1.into()), Err(ContextError::Protected(x)));
    ctx.unprotect(x).unwrap();
    assert_eq!(ctx.attributes(x), A::FLAT);
    ctx.set(x, 1.into()).unwrap();
    ctx.unprotect(x).unwrap(); // unprotected already: fine
    ctx.set_attributes(x, A::PROTECTED).unwrap();
    ctx.unprotect(x).unwrap();
    assert_eq!(ctx.attributes(x), A::NONE);
}

#[test]
fn values_and_ids() {
    let x = sym("ctx_values_x");
    let mut ctx = Context::new();
    let ids = &mut vec![ctx.id()];
    assert_eq!(ctx.value(x), None);
    ctx.set(x, 2.into()).unwrap();
    ids.push(ctx.id());
    assert_eq!(ctx.value(x), Some(&2.into()));
    let copy = ctx.clone();
    assert_eq!(copy.id(), ctx.id(), "a clone is the same state");
    ctx.unset(x).unwrap();
    ids.push(ctx.id());
    assert_eq!(ctx.value(x), None);
    assert_eq!(copy.value(x), Some(&2.into()), "copy on write");
    let id = ctx.id();
    ctx.unset(x).unwrap();
    assert_eq!(ctx.id(), id, "removing nothing changes nothing");
    ctx.set_limits(Limits {
        recursion: 10,
        ..Limits::default()
    });
    ids.push(ctx.id());
    assert_eq!(ctx.limits().recursion, 10);
    assert_eq!(copy.limits(), Limits::default());
    assert_eq!(ids.iter().collect::<HashSet<_>>().len(), ids.len());
    let e = node("ctx_values_e");
    assert_eq!(ctx.eval(&e), Ok(e.clone()));
}

#[test]
fn child_scopes_do_not_leak() {
    let (x, y, f) = (sym("ctx_scope_x"), sym("ctx_scope_y"), sym("ctx_scope_f"));
    let (a1, a2) = (node("ctx_scope_a1"), node("ctx_scope_a2"));
    let mut parent = Context::new();
    parent.set(x, 1.into()).unwrap();
    let before = parent.id();

    let mut child = parent.assuming(a1.clone());
    child.set(x, 2.into()).unwrap();
    child.set(y, 3.into()).unwrap();
    child.set_attributes(f, Attributes::FLAT).unwrap();
    let grandchild = child.assuming(a2.clone());

    assert_eq!(grandchild.assumptions(), &[a1.clone(), a2]);
    assert_eq!(child.assumptions(), &[a1]);
    assert_eq!(grandchild.value(x), Some(&2.into()));
    assert_eq!(grandchild.attributes(f), Attributes::FLAT);

    assert_eq!(parent.id(), before);
    assert!(parent.assumptions().is_empty());
    assert_eq!(parent.value(x), Some(&1.into()));
    assert_eq!(parent.value(y), None);
    assert_eq!(parent.attributes(f), Attributes::NONE);
    drop((child, grandchild));
    assert_eq!(parent.value(x), Some(&1.into()));
}

/// `forgetting` drops the assumptions that mention a symbol and keeps the
/// others, in order, as Xcas's `purge` does.
#[test]
fn forgetting_drops_the_assumptions_on_a_symbol() {
    let read = |s: &str| Expr::parse_full_form(s).unwrap();
    let mut ctx = Context::new();
    ctx.set(sym("v"), 3.into()).unwrap();
    let ctx = ctx
        .assuming(read("Greater[x, 0]"))
        .assuming(read("Greater[y, 0]"))
        .assuming(read("Less[x, y]"));
    let forgot = ctx.forgetting(sym("x"));
    assert_eq!(forgot.assumptions(), &[read("Greater[y, 0]")]);
    assert_ne!(forgot.id(), ctx.id());
    assert_eq!(forgot.value(sym("v")), Some(&Expr::from(3)));
    let positive = |c: &Context, e: &str| c.eval(&read(&format!("Positive[{e}]"))).unwrap();
    assert_eq!(positive(&ctx, "x"), read("True"));
    assert_eq!(positive(&forgot, "x"), read("Positive[x]"));
    assert_eq!(positive(&forgot, "y"), read("True"));
    // Nothing to forget: the same scope.
    let same = forgot.forgetting(sym("z"));
    assert_eq!(same.id(), forgot.id());
    assert_eq!(ctx.assumptions().len(), 3, "self is unchanged");
}

/// An `Interrupt` stops an evaluation from another thread: here a branching
/// recursion with 2^25 distinct calls (`2x` and `2x + 1`, which the memo
/// cannot share), minutes of work that no limit stops.
#[test]
fn an_interrupt_stops_an_evaluation_from_another_thread() {
    let read = |s: &str| Expr::parse_full_form(s).unwrap();
    let mut ctx = Context::new();
    // Rules are tried in definition order: the base case first.
    ctx.define_delayed(&read("g[0, Pattern[x, Blank[]]]"), &read("x"))
        .unwrap();
    ctx.define_delayed(
        &read("g[Pattern[n, Blank[]], Pattern[x, Blank[]]]"),
        &read("Plus[g[Plus[n, -1], Times[2, x]], g[Plus[n, -1], Plus[Times[2, x], 1]]]"),
    )
    .unwrap();
    let interrupt = Interrupt::new();
    ctx.set_interrupt(interrupt.clone());
    let started = std::time::Instant::now();
    std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(100));
            interrupt.interrupt();
        });
        let (result, _) = ctx.eval_with(&read("g[25, 0]"), &crate::EvalOptions::default());
        assert_eq!(result, Err(EvalError::Interrupted));
    });
    assert!(started.elapsed() < std::time::Duration::from_secs(20));
    // Until it is reset, every evaluation stops; then they run again.
    assert_eq!(ctx.eval(&read("g[2, 0]")), Err(EvalError::Interrupted));
    interrupt.reset();
    assert_eq!(ctx.eval(&read("g[2, 0]")), Ok(read("6")));
    assert_eq!(EvalError::Interrupted.to_string(), "interrupted");
}

/// The matcher checks the flag too, every 1024 steps, so a long match stops
/// before its budget runs out.
#[test]
fn an_interrupt_stops_a_long_match() {
    let read = |s: &str| Expr::parse_full_form(s).unwrap();
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.match_steps = 2000;
    ctx.set_limits(limits);
    let args: Vec<String> = (1..=30).map(|i| i.to_string()).collect();
    let e = read(&format!(
        "MatchQ[f[{}], f[Pattern[a, BlankNullSequence[]], Pattern[b, BlankNullSequence[]], \
         Pattern[c, BlankNullSequence[]], 99]]",
        args.join(", ")
    ));
    assert_eq!(ctx.eval(&e), Err(EvalError::MatchLimit));
    let interrupt = Interrupt::new();
    ctx.set_interrupt(interrupt.clone());
    interrupt.interrupt();
    assert_eq!(ctx.eval(&e), Err(EvalError::Interrupted));
}

#[test]
fn ids_are_unique_across_threads() {
    const THREADS: usize = 8;
    let barrier = Barrier::new(THREADS);
    let x = sym("ctx_ids_x");
    let ids: Vec<StateId> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..THREADS)
            .map(|_| {
                s.spawn(|| {
                    barrier.wait();
                    let mut ids = Vec::new();
                    for i in 0..1000 {
                        let mut ctx = Context::new();
                        ids.push(ctx.id());
                        ctx.set(x, i.into()).unwrap();
                        ids.push(ctx.id());
                        ids.push(ctx.assuming(Expr::from(i)).id());
                    }
                    ids
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });
    assert_eq!(ids.iter().collect::<HashSet<_>>().len(), ids.len());
}

#[test]
fn stamps() {
    let e = node("ctx_stamps_e");
    let mut a = Context::new();
    let b = Context::new();
    assert!(!a.is_fixpoint(&e), "unmarked");
    a.mark_fixpoint(&e);
    assert!(a.is_fixpoint(&e) && a.clone().is_fixpoint(&e));
    assert!(!b.is_fixpoint(&e));
    // The stamp is on the shared node: an equal expression sees it.
    assert!(a.is_fixpoint(&node("ctx_stamps_e")));
    let old = a.clone();
    a.set(sym("ctx_stamps_x"), 1.into()).unwrap();
    assert!(!a.is_fixpoint(&e), "the state changed");
    b.mark_fixpoint(&e);
    assert!(
        b.is_fixpoint(&e) && !old.is_fixpoint(&e),
        "one slot: last mark wins"
    );
    // Atoms have no stamp.
    let x = Expr::from(sym("ctx_stamps_x"));
    a.mark_fixpoint(&x);
    a.mark_fixpoint(&Expr::from(3));
    assert!(!a.is_fixpoint(&x) && !a.is_fixpoint(&Expr::from(3)));
    // Visible from another thread.
    std::thread::scope(|s| {
        s.spawn(|| a.mark_fixpoint(&e));
    });
    assert!(a.is_fixpoint(&e) && !b.is_fixpoint(&e));
}

/// An operation on context `i` (modulo the number of contexts) or node `n`.
#[derive(Clone, Debug)]
enum Op {
    Clone(usize),
    Set(usize),
    Assume(usize),
    SetAttributes(usize),
    Mark(usize, usize),
    Check(usize, usize),
}

fn op() -> impl Strategy<Value = Op> {
    let i = 0..8usize;
    let n = 0..4usize;
    prop_oneof![
        i.clone().prop_map(Op::Clone),
        i.clone().prop_map(Op::Set),
        i.clone().prop_map(Op::Assume),
        i.clone().prop_map(Op::SetAttributes),
        (i.clone(), n.clone()).prop_map(|(i, n)| Op::Mark(i, n)),
        (i, n).prop_map(|(i, n)| Op::Check(i, n)),
    ]
}

proptest! {
    /// Against a model where each state is a token, a clone shares its
    /// parent's token and every change takes a new one: ids are equal iff
    /// tokens are, and a node is a fixpoint in a context iff it was last
    /// marked in that context's state.
    #[test]
    fn stamps_match_only_their_state(ops in prop::collection::vec(op(), 1..64)) {
        // Fresh nodes per case, so earlier cases leave no stamps.
        static CASE: AtomicU64 = AtomicU64::new(0);
        let case = CASE.fetch_add(1, Ordering::Relaxed);
        let nodes: Vec<Expr> = (0..4).map(|n| node(&format!("ctx_prop_{case}_{n}"))).collect();
        let (x, f) = (sym("ctx_prop_x"), sym("ctx_prop_f"));
        let mut tokens = 0u64;
        let mut ctxs = vec![(Context::new(), 0u64)];
        let mut marked: Vec<Option<u64>> = vec![None; nodes.len()];
        for op in ops {
            let len = ctxs.len();
            match op {
                Op::Clone(i) => {
                    let c = ctxs[i % len].clone();
                    ctxs.push(c);
                }
                Op::Set(i) | Op::Assume(i) | Op::SetAttributes(i) => {
                    let (ctx, token) = &mut ctxs[i % len];
                    match op {
                        Op::Set(_) => ctx.set(x, Expr::from(i64::try_from(i).unwrap())).unwrap(),
                        Op::Assume(_) => *ctx = ctx.assuming(Expr::from(x)),
                        _ => ctx.set_attributes(f, Attributes::FLAT).unwrap(),
                    }
                    tokens += 1;
                    *token = tokens;
                }
                Op::Mark(i, n) => {
                    let (ctx, token) = &ctxs[i % len];
                    ctx.mark_fixpoint(&nodes[n]);
                    marked[n] = Some(*token);
                }
                Op::Check(i, n) => {
                    let (ctx, token) = &ctxs[i % len];
                    prop_assert_eq!(ctx.is_fixpoint(&nodes[n]), marked[n] == Some(*token));
                }
            }
        }
        for (a, ta) in &ctxs {
            for (b, tb) in &ctxs {
                prop_assert_eq!(a.id() == b.id(), ta == tb);
            }
            for (n, e) in nodes.iter().enumerate() {
                prop_assert_eq!(a.is_fixpoint(e), marked[n] == Some(*ta));
            }
        }
    }
}
