//! Tests of the rule store and dispatch. Rules live in contexts, so tests
//! share symbol names freely.

use std::sync::Arc;

use proptest::prelude::*;

use super::*;
use crate::Limits;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn sym(name: &str) -> Symbol {
    Symbol::new(name)
}

/// `(input, result)` in full form, evaluated in `ctx` with and without the
/// caches.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        for found in [ctx.eval(&read(input)), ctx.eval_uncached(&read(input))] {
            match found {
                Ok(found) if format!("{found:?}") == *expected => {}
                found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn define(ctx: &mut Context, lhs: &str, rhs: &str) {
    ctx.define(&read(lhs), &read(rhs))
        .unwrap_or_else(|e| panic!("{lhs} = {rhs}: {e}"));
}

fn delayed(ctx: &mut Context, lhs: &str, rhs: &str) {
    ctx.define_delayed(&read(lhs), &read(rhs))
        .unwrap_or_else(|e| panic!("{lhs} := {rhs}: {e}"));
}

/// The definitions of `s` in the order they are tried, `lhs -> rhs` in full
/// form: the exact ones in definition order, then the down-rules, the
/// sub-rules and the up-rules.
fn definitions(ctx: &Context, s: Symbol) -> Vec<String> {
    let Some(own) = ctx.rules().symbols.get(s) else {
        return Vec::new();
    };
    let mut exact: Vec<(&Expr, &Exact)> = own.exact.iter().collect();
    exact.sort_by_key(|(_, x)| x.rank);
    exact
        .into_iter()
        .map(|(lhs, x)| format!("{lhs:?} -> {:?}", x.rhs))
        .chain(
            [&own.down, &own.sub, &own.up]
                .into_iter()
                .flat_map(|list| &list.rules)
                .map(|r| format!("{:?} -> {:?}", r.lhs, r.rule.rhs())),
        )
        .collect()
}

#[test]
fn exact_definitions_are_keyed_by_the_evaluated_left_hand_side() {
    let mut ctx = Context::new();
    define(&mut ctx, "f[Plus[1, 1]]", "two");
    define(&mut ctx, "f[a, Times[2, 3]]", "six");
    check(
        &ctx,
        &[
            ("f[2]", "two"),
            ("f[Plus[3, -1]]", "two"),
            ("f[a, 6]", "six"),
            ("f[3]", "f[3]"),
            ("f[6, a]", "f[6, a]"),
        ],
    );
    assert_eq!(
        definitions(&ctx, sym("f")),
        ["f[2] -> two", "f[a, 6] -> six"]
    );
}

#[test]
fn orderless_left_hand_sides_are_sorted() {
    let mut ctx = Context::new();
    ctx.set_attributes(sym("f"), Attributes::ORDERLESS).unwrap();
    define(&mut ctx, "f[b, a]", "1");
    delayed(&mut ctx, "f[g[Pattern[x, Blank[]]], c]", "x");
    check(
        &ctx,
        &[
            ("f[a, b]", "1"),
            ("f[b, a]", "1"),
            ("f[c, g[1]]", "1"),
            ("f[g[2], c]", "2"),
            ("f[g[2], d]", "f[d, g[2]]"),
        ],
    );
    assert_eq!(
        definitions(&ctx, sym("f")),
        ["f[a, b] -> 1", "f[c, g[Pattern[x, Blank[]]]] -> x"]
    );
}

#[test]
fn patterned_rules_apply_in_definition_order_after_exact_ones() {
    let int = "f[PatternTest[Pattern[x, Blank[]], IntegerQ]]";
    let any = "f[Pattern[x, Blank[]]]";
    let mut ctx = Context::new();
    delayed(&mut ctx, int, "int");
    delayed(&mut ctx, any, "any");
    define(&mut ctx, "f[0]", "zero");
    check(&ctx, &[("f[1]", "int"), ("f[a]", "any"), ("f[0]", "zero")]);
    // Not by specificity: the general rule first hides the other.
    let mut ctx = Context::new();
    delayed(&mut ctx, any, "any");
    delayed(&mut ctx, int, "int");
    check(&ctx, &[("f[1]", "any"), ("f[a]", "any")]);
}

#[test]
fn recursion_through_rules() {
    let mut ctx = Context::new();
    define(&mut ctx, "fact[0]", "1");
    delayed(
        &mut ctx,
        "fact[Pattern[n, Blank[]]]",
        "Times[n, fact[Plus[n, -1]]]",
    );
    let product = |n: i64| {
        let factors: Vec<Expr> = (1..=n).map(Expr::from).collect();
        Context::new()
            .eval(&Expr::apply(&Symbol::TIMES.into(), &factors))
            .unwrap()
    };
    let fact = |n: i64| Expr::apply(&sym("fact").into(), &[n.into()]);
    assert_eq!(ctx.eval(&fact(20)), Ok(product(20)));
    // Each step takes a few levels: too many for the caller's thread...
    assert_eq!(ctx.eval(&fact(1000)), Err(EvalError::RecursionLimit));
    // ...but not for a thread of its own, on which definitions run too.
    ctx.set_limits(Limits::deep());
    define(&mut ctx, "fact[1]", "1");
    assert_eq!(ctx.eval(&fact(1000)), Ok(product(1000)));
}

#[test]
fn a_redefinition_replaces_in_place() {
    let f = sym("f");
    let mut ctx = Context::new();
    delayed(&mut ctx, "f[Pattern[x, Blank[]]]", "first[x]");
    delayed(
        &mut ctx,
        "f[PatternTest[Pattern[x, Blank[]], IntegerQ]]",
        "int",
    );
    // The first rule, up to the name of its variable: still first.
    delayed(&mut ctx, "f[Pattern[y, Blank[]]]", "second[y]");
    check(&ctx, &[("f[1]", "second[1]"), ("f[a]", "second[a]")]);
    assert_eq!(
        definitions(&ctx, f),
        [
            "f[Pattern[y, Blank[]]] -> second[y]",
            "f[PatternTest[Pattern[x, Blank[]], IntegerQ]] -> int"
        ]
    );
    // Exact definitions too, by their evaluated left-hand sides.
    define(&mut ctx, "f[0]", "zero");
    define(&mut ctx, "f[Plus[1, -1]]", "nought");
    check(&ctx, &[("f[0]", "nought")]);
    assert_eq!(definitions(&ctx, f).len(), 3);
}

#[test]
fn rules_with_other_conditions_are_other_rules() {
    let mut ctx = Context::new();
    delayed(
        &mut ctx,
        "f[Pattern[x, Blank[]]]",
        "Condition[pos, Greater[x, 0]]",
    );
    delayed(
        &mut ctx,
        "f[Pattern[y, Blank[]]]",
        "Condition[neg, Less[y, 0]]",
    );
    check(
        &ctx,
        &[
            ("f[1]", "pos"),
            ("f[-1]", "neg"),
            ("f[0]", "f[0]"),
            ("f[a]", "f[a]"),
        ],
    );
    // The first rule again: the same pattern and condition, on the left.
    delayed(
        &mut ctx,
        "Condition[f[Pattern[z, Blank[]]], Greater[z, 0]]",
        "positive",
    );
    check(&ctx, &[("f[1]", "positive"), ("f[-1]", "neg")]);
    assert_eq!(
        definitions(&ctx, sym("f")),
        [
            "Condition[f[Pattern[z, Blank[]]], Greater[z, 0]] -> positive",
            "f[Pattern[y, Blank[]]] -> Condition[neg, Less[y, 0]]"
        ]
    );
}

#[test]
fn a_literal_that_shares_a_variable_name_is_not_renamed() {
    let mut ctx = Context::new();
    delayed(&mut ctx, "f[Pattern[x, Blank[]], x]", "one");
    delayed(&mut ctx, "f[Pattern[y, Blank[]], y]", "two");
    check(
        &ctx,
        &[
            ("f[a, x]", "one"),
            ("f[a, y]", "two"),
            ("f[a, a]", "f[a, a]"),
        ],
    );
    assert_eq!(definitions(&ctx, sym("f")).len(), 2);
}

#[test]
fn set_evaluates_now_and_set_delayed_at_each_application() {
    let mut ctx = Context::new();
    let y = sym("y");
    ctx.set(y, 1.into()).unwrap();
    define(&mut ctx, "now[]", "y");
    delayed(&mut ctx, "later[]", "y");
    define(&mut ctx, "a", "y");
    delayed(&mut ctx, "b", "y");
    ctx.set(y, 2.into()).unwrap();
    check(
        &ctx,
        &[("now[]", "1"), ("later[]", "2"), ("a", "1"), ("b", "2")],
    );
}

#[test]
fn sub_rules_attach_to_the_innermost_head() {
    let mut ctx = Context::new();
    delayed(
        &mut ctx,
        "f[Pattern[a, Blank[]]][Pattern[b, Blank[]]]",
        "List[a, b]",
    );
    // The head is normalised too: the rule is for `g[2][…]`.
    delayed(&mut ctx, "g[Plus[1, 1]][Pattern[x, Blank[]]]", "x");
    check(
        &ctx,
        &[
            ("f[1][2]", "List[1, 2]"),
            ("f[1]", "f[1]"),
            ("f[1][2][3]", "List[1, 2][3]"),
            ("g[2][5]", "5"),
            ("g[Plus[1, 1]][5]", "5"),
            ("g[3][5]", "g[3][5]"),
        ],
    );
    assert_eq!(
        definitions(&ctx, sym("f")),
        ["f[Pattern[a, Blank[]]][Pattern[b, Blank[]]] -> List[a, b]"]
    );
}

#[test]
fn up_rules_extend_builtins() {
    let mut ctx = Context::new();
    let (phi, vec, f) = (sym("phi"), sym("vec"), sym("f"));
    // On an atomic argument.
    ctx.define_up(phi, &read("Power[phi, 2]"), &read("Plus[phi, 1]"))
        .unwrap();
    // On an argument's head, before the native function of `Length`.
    ctx.define_up_delayed(
        vec,
        &read("Length[vec[BlankNullSequence[]]]"),
        &read("dimension"),
    )
    .unwrap();
    // On the head a blank requires.
    ctx.define_up_delayed(f, &read("g[Pattern[x, Blank[f]]]"), &read("gf[x]"))
        .unwrap();
    check(
        &ctx,
        &[
            ("Times[phi, phi]", "Plus[1, phi]"),
            ("Power[phi, 3]", "Power[phi, 3]"),
            ("Length[vec[1, 2]]", "dimension"),
            ("Length[List[1, 2]]", "2"),
            ("g[f[1]]", "gf[f[1]]"),
            ("g[h[1]]", "g[h[1]]"),
            ("Plus[f[1], phi, vec[]]", "Plus[phi, f[1], vec[]]"),
        ],
    );
    assert_eq!(ctx.rules().up_count, 3);
    for s in [phi, vec, f] {
        ctx.clear(s).unwrap();
    }
    assert_eq!((ctx.rules().up_count, ctx.rules().up_bits), (0, 0));
    check(&ctx, &[("Length[vec[1, 2]]", "2")]);
}

#[test]
fn up_rules_and_values_exclude_each_other() {
    let mut ctx = Context::new();
    let (t, v) = (sym("t"), sym("v"));
    ctx.define_up(t, &read("Sin[t]"), &read("half")).unwrap();
    assert_eq!(ctx.set(t, 1.into()), Err(ContextError::HasUpRules(t)));
    assert_eq!(
        ctx.define(&read("t"), &read("1")),
        Err(DefineError::HasUpRules(t))
    );
    ctx.set(v, 1.into()).unwrap();
    assert_eq!(
        ctx.define_up(v, &read("Sin[v]"), &read("half")),
        Err(DefineError::HasValue(v))
    );
    // The tag must be at level 1.
    assert_eq!(
        ctx.define_up(t, &read("Sin[Cos[t]]"), &read("1")),
        Err(DefineError::NotATag(t))
    );
    assert_eq!(
        ctx.define_up(t, &read("t"), &read("1")),
        Err(DefineError::NotATag(t))
    );
    check(&ctx, &[("Sin[t]", "half")]);
    // Once cleared, a value is fine.
    ctx.clear(t).unwrap();
    ctx.set(t, 1.into()).unwrap();
    check(&ctx, &[("Sin[t]", "Sin[1]")]);
}

#[test]
fn protected_symbols_take_no_definitions() {
    let mut ctx = Context::new();
    let id = ctx.id();
    let plus = Symbol::PLUS;
    assert_eq!(
        ctx.define(&read("Plus[Pattern[a, Blank[]], b]"), &read("0")),
        Err(DefineError::Protected(plus))
    );
    assert_eq!(
        ctx.define_delayed(&read("Pi"), &read("3")),
        Err(DefineError::Protected(Symbol::PI))
    );
    assert_eq!(
        ctx.define_delayed(&read("List[a][Pattern[x, Blank[]]]"), &read("x")),
        Err(DefineError::Protected(Symbol::LIST))
    );
    assert_eq!(
        ctx.define_up(plus, &read("f[Plus]"), &read("0")),
        Err(DefineError::Protected(plus))
    );
    assert_eq!(
        ctx.define_native(plus, Arc::new(|_, _| Ok(None))),
        Err(ContextError::Protected(plus))
    );
    assert_eq!(
        ctx.set_default(plus, 0.into()),
        Err(ContextError::Protected(plus))
    );
    assert_eq!(ctx.clear(plus), Err(ContextError::Protected(plus)));
    assert_eq!(ctx.id(), id, "a rejected definition is no change");
    // A user symbol, protected, then not.
    let p = sym("p");
    ctx.set_attributes(p, Attributes::PROTECTED).unwrap();
    assert_eq!(
        ctx.define(&read("p[1]"), &read("1")),
        Err(DefineError::Protected(p))
    );
    ctx.unprotect(p).unwrap();
    define(&mut ctx, "p[1]", "one");
    check(&ctx, &[("p[1]", "one")]);
}

#[test]
fn clear_keeps_attributes_and_default() {
    let mut ctx = Context::new();
    let f = sym("f");
    ctx.set_attributes(f, Attributes::ORDERLESS).unwrap();
    ctx.set_default(f, 7.into()).unwrap();
    define(&mut ctx, "f[1]", "one");
    delayed(&mut ctx, "f[Pattern[x, Blank[]], b]", "rule");
    ctx.define_native(f, Arc::new(|_, _| Ok(Some(Expr::from(0)))))
        .unwrap();
    check(&ctx, &[("f[1]", "0")]);
    ctx.clear(f).unwrap();
    check(&ctx, &[("f[1]", "f[1]"), ("f[b, a]", "f[a, b]")]);
    assert!(definitions(&ctx, f).is_empty());
    assert_eq!(ctx.attributes(f), Attributes::ORDERLESS);
    delayed(
        &mut ctx,
        "f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]]]]",
        "List[x, y]",
    );
    check(&ctx, &[("f[1]", "List[1, 7]")]);
    // Clearing a symbol without anything is no change.
    let id = ctx.id();
    ctx.clear(sym("nothing")).unwrap();
    assert_eq!(ctx.id(), id);
}

#[test]
fn user_natives_run_before_rules_and_may_evaluate() {
    let mut ctx = Context::new();
    let (f, twice) = (sym("f"), sym("twice"));
    ctx.define_native(
        f,
        Arc::new(|_, e| {
            Ok(match e.args() {
                [n] if n.as_i64().is_some() => {
                    Some(Expr::apply(&Symbol::PLUS.into(), &[n.clone(), 1.into()]))
                }
                _ => None,
            })
        }),
    )
    .unwrap();
    delayed(&mut ctx, "f[Pattern[x, Blank[]]]", "rule[x]");
    define(&mut ctx, "f[2]", "exact");
    ctx.define_native(
        twice,
        Arc::new(|ev, e| {
            let [x] = e.args() else { return Ok(None) };
            ev.eval(&Expr::apply(&Symbol::PLUS.into(), &[x.clone(), x.clone()]))
                .map(Some)
        }),
    )
    .unwrap();
    check(
        &ctx,
        &[
            ("f[1]", "2"),
            ("f[2]", "3"),
            ("f[a]", "rule[a]"),
            ("twice[a]", "Times[2, a]"),
        ],
    );
}

#[test]
fn defaults_fill_missing_arguments_and_recompile() {
    let mut ctx = Context::new();
    let (f, h) = (sym("f"), sym("h"));
    ctx.set_default(f, 7.into()).unwrap();
    delayed(
        &mut ctx,
        "f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]]]]",
        "List[x, y]",
    );
    check(&ctx, &[("f[1]", "List[1, 7]"), ("f[1, 2]", "List[1, 2]")]);
    ctx.set_default(f, 8.into()).unwrap();
    check(&ctx, &[("f[1]", "List[1, 8]")]);
    // A rule of `g` that mentions `h`, which has no rules.
    ctx.set_default(h, 8.into()).unwrap();
    delayed(&mut ctx, "g[h[Optional[Pattern[y, Blank[]]]]]", "y");
    check(&ctx, &[("g[h[]]", "8"), ("g[h[3]]", "3")]);
    ctx.set_default(h, 9.into()).unwrap();
    check(&ctx, &[("g[h[]]", "9")]);
    // Without a default, `x_.` must be there.
    delayed(&mut ctx, "k[Optional[Pattern[y, Blank[]]]]", "y");
    check(&ctx, &[("k[]", "k[]"), ("k[4]", "4")]);
}

#[test]
fn an_attribute_change_recompiles_the_rules_that_mention_the_symbol() {
    let mut ctx = Context::new();
    let f = sym("f");
    // Rules of `g` and `h` that mention `f`, and one of `f`.
    define(&mut ctx, "g[f[b, a]]", "exact");
    delayed(&mut ctx, "h[f[Pattern[x, Blank[]], b]]", "x");
    define(&mut ctx, "f[d, c]", "own");
    check(
        &ctx,
        &[
            ("g[f[b, a]]", "exact"),
            ("h[f[c, b]]", "c"),
            ("h[f[b, c]]", "h[f[b, c]]"),
            ("f[c, d]", "f[c, d]"),
        ],
    );
    ctx.set_attributes(f, Attributes::ORDERLESS).unwrap();
    // Subjects are sorted now, and so are the rules.
    check(
        &ctx,
        &[
            ("g[f[b, a]]", "exact"),
            ("g[f[a, b]]", "exact"),
            ("h[f[c, b]]", "c"),
            ("h[f[b, c]]", "c"),
            ("f[c, d]", "own"),
            ("f[d, c]", "own"),
        ],
    );
    assert_eq!(definitions(&ctx, sym("g")), ["g[f[a, b]] -> exact"]);
    assert_eq!(
        definitions(&ctx, sym("h")),
        ["h[f[b, Pattern[x, Blank[]]]] -> x"]
    );
    // Flat as well: `f[b, f[a, c]]` is `f[a, b, c]`.
    delayed(&mut ctx, "k[f[a, f[b, Pattern[x, Blank[]]]]]", "x");
    ctx.set_attributes(f, Attributes::ORDERLESS | Attributes::FLAT)
        .unwrap();
    check(&ctx, &[("k[f[c, a, b]]", "c"), ("k[f[b, f[c, a]]]", "c")]);
}

#[test]
fn held_left_hand_sides() {
    let mut ctx = Context::new();
    ctx.set_attributes(sym("f"), Attributes::HOLD_ALL).unwrap();
    define(&mut ctx, "f[Plus[1, 1]]", "held");
    check(&ctx, &[("f[Plus[1, 1]]", "held"), ("f[2]", "f[2]")]);
    // `HoldPattern` keeps the left-hand side as written, which no
    // evaluated subject is.
    delayed(&mut ctx, "HoldPattern[k[Plus[1, 1]]]", "never");
    check(&ctx, &[("k[Plus[1, 1]]", "k[2]")]);
    assert_eq!(
        definitions(&ctx, sym("k")),
        ["HoldPattern[k[Plus[1, 1]]] -> never"]
    );
}

#[test]
fn define_many_is_all_or_nothing() {
    let pairs = |defs: &[(&str, &str)]| -> Vec<(Expr, Expr)> {
        defs.iter().map(|(l, r)| (read(l), read(r))).collect()
    };
    let mut ctx = Context::new();
    let library = pairs(&[
        ("fact[0]", "1"),
        ("fact[Pattern[n, Blank[]]]", "Times[n, fact[Plus[n, -1]]]"),
        // Each definition sees those before it: this is `m[2]`.
        ("k[1]", "2"),
        ("m[k[1]]", "three"),
    ]);
    ctx.define_many(library.iter().map(|(l, r)| (l, r)))
        .unwrap();
    check(&ctx, &[("fact[4]", "24"), ("m[2]", "three")]);
    let id = ctx.id();
    let bad = pairs(&[("g[1]", "one"), ("Plus[a]", "1")]);
    assert_eq!(
        ctx.define_many(bad.iter().map(|(l, r)| (l, r))),
        Err(DefineError::Protected(Symbol::PLUS))
    );
    assert_eq!(ctx.id(), id);
    check(&ctx, &[("g[1]", "g[1]")]);
}

#[test]
fn a_definition_changes_only_its_context() {
    let mut ctx = Context::new();
    define(&mut ctx, "f[1]", "one");
    let before = ctx.clone();
    define(&mut ctx, "f[2]", "two");
    delayed(&mut ctx, "f[Pattern[x, Blank[]]]", "any");
    assert_ne!(ctx.id(), before.id());
    check(
        &before,
        &[("f[1]", "one"), ("f[2]", "f[2]"), ("f[3]", "f[3]")],
    );
    check(&ctx, &[("f[1]", "one"), ("f[2]", "two"), ("f[3]", "any")]);
}

#[test]
fn malformed_definitions() {
    let mut ctx = Context::new();
    for lhs in ["1", "\"s\"", "HoldPattern[x]"] {
        assert_eq!(
            ctx.define(&read(lhs), &read("2")),
            Err(DefineError::NoHead(read(lhs)))
        );
    }
    assert!(matches!(
        ctx.define_delayed(&read("f[Pattern[1, Blank[]]]"), &read("x")),
        Err(DefineError::Pattern(_))
    ));
    let mut limits = ctx.limits();
    limits.recursion = 4;
    ctx.set_limits(limits);
    assert_eq!(
        ctx.define(&read("f[x]"), &read("g[g[g[g[g[1]]]]]")),
        Err(DefineError::Eval(EvalError::RecursionLimit))
    );
    let messages = [
        (
            DefineError::Protected(Symbol::PLUS),
            "symbol Plus is protected",
        ),
        (
            DefineError::HasValue(sym("v")),
            "symbol v has a value, which would hide its up-rules",
        ),
        (
            DefineError::HasUpRules(sym("t")),
            "symbol t has up-rules, which a value would hide",
        ),
        (
            DefineError::NotATag(sym("t")),
            "symbol t is not an argument of the left-hand side or the head of one",
        ),
        (DefineError::NoHead(read("1")), "1 has no symbol to define"),
        (
            DefineError::Eval(EvalError::RecursionLimit),
            "recursion limit exceeded",
        ),
    ];
    for (e, message) in messages {
        assert_eq!(e.to_string(), message);
    }
}

#[test]
fn system_rules_apply_to_builtin_heads() {
    let system = system(&[
        (read("Factorial[0]"), read("1")),
        (
            read("Factorial[Pattern[n, Blank[]]]"),
            read("Times[n, Factorial[Plus[n, -1]]]"),
        ),
        (
            read("Range[Pattern[n, Blank[]]][Pattern[k, Blank[]]]"),
            read("Plus[n, k]"),
        ),
    ]);
    let mut ctx = Context::new();
    define(&mut ctx, "fact[0]", "1");
    let mut ev = Eval::new(&ctx, 64);
    let factorial = Symbol::FACTORIAL;
    let own = system.symbols.get(factorial).unwrap();
    assert_eq!(
        own_rules(own, &mut ev, &read("Factorial[0]")),
        Ok(Some(read("1")))
    );
    assert_eq!(
        own_rules(own, &mut ev, &read("Factorial[3]")),
        Ok(Some(read("Times[3, Factorial[Plus[3, -1]]]")))
    );
    assert_eq!(
        sub_rules(&system, &mut ev, &read("Range[2][3]"), Symbol::RANGE),
        Ok(Some(read("Plus[2, 3]")))
    );
    // A user's symbol has the rules of the context.
    assert_eq!(
        down(&mut ev, &read("fact[0]"), sym("fact")),
        Ok(Some(read("1")))
    );
    // The shared set has none yet.
    assert_eq!(down(&mut ev, &read("Factorial[0]"), factorial), Ok(None));
    assert_eq!(sub(&mut ev, &read("Range[2][3]")), Ok(None));
}

#[test]
fn rules_print_their_size() {
    let mut ctx = Context::new();
    define(&mut ctx, "f[1]", "1");
    delayed(&mut ctx, "f[Pattern[x, Blank[]]]", "x");
    assert_eq!(
        format!("{:?}", ctx.rules()),
        "Rules { symbols: 1, up: 0, .. }"
    );
    assert_eq!(
        format!("{:?}", ctx.rules().symbols.get(sym("f")).unwrap()),
        "SymbolRules { exact: 1, down: 1, up: 0, sub: 0, native: false, default: None, .. }"
    );
}

/// Atoms, and with `patterns` pattern variables and blanks, under `f`.
fn argument(patterns: bool) -> impl Strategy<Value = Expr> {
    let atom = prop::sample::select(vec!["a", "b", "c"]).prop_map(|s| Expr::from(sym(s)));
    let leaf = if patterns {
        prop_oneof![
            2 => atom,
            1 => prop::sample::select(vec![
                "Pattern[x, Blank[]]",
                "Pattern[y, Blank[]]",
                "Blank[]",
                "Pattern[x, BlankSequence[]]",
                "BlankNullSequence[]",
                "Optional[Pattern[y, Blank[]], b]",
                "Pattern[x, Blank[f]]",
                "PatternTest[Pattern[y, Blank[]], SymbolQ]",
            ])
            .prop_map(read),
        ]
        .boxed()
    } else {
        atom.boxed()
    };
    leaf.prop_recursive(2, 8, 3, |inner| {
        prop::collection::vec(inner, 0..3).prop_map(|args| Expr::apply(&sym("f").into(), &args))
    })
}

fn g(args: &[Expr]) -> Expr {
    Expr::apply(&sym("g").into(), args)
}

/// Definitions of `g`, the right-hand side of the `i`-th `r[i, x, y]`.
fn define_all(ctx: &mut Context, lhs: &[Vec<Expr>]) {
    for (i, args) in (0i64..).zip(lhs) {
        let rhs = Expr::apply(
            &sym("r").into(),
            &[i.into(), sym("x").into(), sym("y").into()],
        );
        ctx.define_delayed(&g(args), &rhs).unwrap();
    }
}

proptest! {
    /// Setting attributes after the definitions gives the rules and the
    /// results that setting them before gives: recompilation is defining
    /// again under the new attributes.
    #[test]
    fn attributes_set_later_recompile(
        lhs in prop::collection::vec(prop::collection::vec(argument(true), 0..3), 1..8),
        subjects in prop::collection::vec(prop::collection::vec(argument(false), 0..3), 0..8),
        f_attributes in prop::sample::select(vec![
            Attributes::FLAT,
            Attributes::ORDERLESS,
            Attributes::FLAT | Attributes::ORDERLESS,
        ]),
        g_attributes in prop::sample::select(vec![
            Attributes::NONE,
            Attributes::ORDERLESS,
        ]),
    ) {
        let (f, g_) = (sym("f"), sym("g"));
        let mut later = Context::new();
        define_all(&mut later, &lhs);
        later.set_attributes(f, f_attributes).unwrap();
        later.set_attributes(g_, g_attributes).unwrap();
        let mut first = Context::new();
        first.set_attributes(f, f_attributes).unwrap();
        first.set_attributes(g_, g_attributes).unwrap();
        define_all(&mut first, &lhs);
        prop_assert_eq!(definitions(&later, g_), definitions(&first, g_));
        for args in &subjects {
            prop_assert_eq!(later.eval(&g(args)), first.eval(&g(args)));
        }
    }

    /// The filters skip only rules that would not match, whatever the
    /// attributes of the head.
    #[test]
    fn filters_change_nothing(
        lhs in prop::collection::vec(prop::collection::vec(argument(true), 0..4), 1..8),
        subjects in prop::collection::vec(prop::collection::vec(argument(false), 0..4), 1..8),
        attributes in prop::sample::select(vec![
            Attributes::NONE,
            Attributes::ORDERLESS,
            Attributes::FLAT,
            Attributes::FLAT | Attributes::ORDERLESS,
            Attributes::ORDERLESS | Attributes::ONE_IDENTITY,
            Attributes::FLAT | Attributes::ORDERLESS | Attributes::ONE_IDENTITY,
        ]),
    ) {
        let mut ctx = Context::new();
        ctx.set_attributes(sym("g"), attributes).unwrap();
        define_all(&mut ctx, &lhs);
        let mut unfiltered = ctx.clone();
        let own = unfiltered.rules_mut().symbols.entry(sym("g"), SymbolRules::default);
        for stored in &mut own.down.rules {
            stored.filter = Filter::default();
        }
        for args in &subjects {
            prop_assert_eq!(ctx.eval(&g(args)), unfiltered.eval(&g(args)));
        }
    }
}
