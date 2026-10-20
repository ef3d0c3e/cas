//! Tests of patterns, rules, `Function` and the `Replace` family.

use std::collections::HashMap;

use proptest::prelude::*;

use super::{Pattern, PatternErrorKind as K, Rule};
use crate::{Context, EvalError, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form, each evaluated in `ctx`.
fn check_in(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn check(cases: &[(&str, &str)]) {
    check_in(&Context::new(), cases);
}

#[test]
fn blanks_and_bindings() {
    check(&[
        ("MatchQ[x, Blank[]]", "True"),
        ("MatchQ[f[x], Blank[f]]", "True"),
        ("MatchQ[g[x], Blank[f]]", "False"),
        ("MatchQ[1, Blank[f]]", "False"),
        ("MatchQ[f[1, 2], f[Blank[], Blank[]]]", "True"),
        ("MatchQ[f[1], f[Blank[], Blank[]]]", "False"),
        ("MatchQ[f[1, 2], g[Blank[], Blank[]]]", "False"),
        ("MatchQ[f[1, 2], Blank[][Blank[], Blank[]]]", "True"),
        // A repeated name binds equal expressions.
        (
            "MatchQ[f[a, a], f[Pattern[x, Blank[]], Pattern[x, Blank[]]]]",
            "True",
        ),
        (
            "MatchQ[f[a, b], f[Pattern[x, Blank[]], Pattern[x, Blank[]]]]",
            "False",
        ),
        (
            "ReplaceAll[h[f[1, 2]], RuleDelayed[f[Pattern[x, Blank[]], Pattern[y, Blank[]]], g[y, x]]]",
            "h[g[2, 1]]",
        ),
        // Literal parts match only themselves.
        ("MatchQ[f[1, x], f[1, Blank[]]]", "True"),
        ("MatchQ[f[2, x], f[1, Blank[]]]", "False"),
    ]);
}

#[test]
fn sequences_take_the_fewest_first() {
    check(&[
        (
            "Replace[f[1, 2, 3], RuleDelayed[f[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]]], g[List[x], List[y]]]]",
            "g[List[1], List[2, 3]]",
        ),
        (
            "Replace[f[1, 2, 3], RuleDelayed[f[Pattern[x, BlankSequence[]], Pattern[y, Blank[]]], g[List[x], y]]]",
            "g[List[1, 2], 3]",
        ),
        (
            "Replace[f[], RuleDelayed[f[Pattern[x, BlankNullSequence[]]], g[List[x]]]]",
            "g[List[]]",
        ),
        ("MatchQ[f[], f[BlankSequence[]]]", "False"),
        ("MatchQ[f[a, b], f[BlankNullSequence[], b]]", "True"),
        // A sequence with a head: every element has it.
        ("MatchQ[f[g[1], g[2]], f[BlankSequence[g]]]", "True"),
        ("MatchQ[f[g[1], h[2]], f[BlankSequence[g]]]", "False"),
        // Bindings agree across nested parts: the first split of the inner
        // sequence fails at the outer one, the second succeeds.
        (
            "Replace[f[g[1, 2, 3], 1, 2], RuleDelayed[f[g[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]]], Pattern[x, BlankSequence[]]], List[List[x], List[y]]]]",
            "List[List[1, 2], List[3]]",
        ),
        // A sequence variable in the right-hand side is spliced.
        (
            "Replace[f[1, 2], RuleDelayed[f[Pattern[x, BlankSequence[]]], g[0, x]]]",
            "g[0, 1, 2]",
        ),
    ]);
}

#[test]
fn alternatives_optional_conditions_tests() {
    check(&[
        ("MatchQ[b, Alternatives[a, b]]", "True"),
        ("MatchQ[c, Alternatives[a, b]]", "False"),
        (
            "Replace[f[b], RuleDelayed[f[Alternatives[Pattern[x, a], Pattern[y, b]]], List[x, y]]]",
            "List[b]",
        ),
        // `Optional`: the default when the argument is missing.
        (
            "Replace[f[1], RuleDelayed[f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], 5]], List[x, y]]]",
            "List[1, 5]",
        ),
        (
            "Replace[f[1, 2], RuleDelayed[f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], 5]], List[x, y]]]",
            "List[1, 2]",
        ),
        // Defaults of the builtins: 1 for an exponent.
        (
            "Replace[Hold[Power[x]], RuleDelayed[Hold[Power[Pattern[b, Blank[]], Optional[Pattern[n, Blank[]]]]], List[b, n]]]",
            "List[x, 1]",
        ),
        ("MatchQ[f[1], f[Blank[], Optional[Blank[]]]]", "False"),
        // `Condition` and `PatternTest`.
        (
            "MatchQ[5, Condition[Pattern[x, Blank[]], Greater[x, 0]]]",
            "True",
        ),
        (
            "MatchQ[-5, Condition[Pattern[x, Blank[]], Greater[x, 0]]]",
            "False",
        ),
        (
            "MatchQ[x, Condition[Pattern[x, Blank[]], Greater[x, 0]]]",
            "False",
        ),
        ("MatchQ[2, PatternTest[Blank[], IntegerQ]]", "True"),
        ("MatchQ[1/2, PatternTest[Blank[], IntegerQ]]", "False"),
        // `HoldPattern`: matched without evaluating the pattern.
        (
            "ReplaceAll[Hold[Plus[1, 2]], RuleDelayed[HoldPattern[Plus[1, 2]], three]]",
            "Hold[three]",
        ),
    ]);
}

#[test]
fn replacing() {
    check(&[
        ("Replace[f[a], RuleDelayed[a, b]]", "f[a]"),
        ("Replace[a, RuleDelayed[a, b]]", "b"),
        ("ReplaceAll[f[a, g[a]], RuleDelayed[a, b]]", "f[b, g[b]]"),
        // The first rule that applies, at the highest part.
        (
            "ReplaceAll[f[a], List[RuleDelayed[a, b], RuleDelayed[a, c]]]",
            "f[b]",
        ),
        (
            "ReplaceAll[f[a], List[RuleDelayed[f[a], 1], RuleDelayed[a, c]]]",
            "1",
        ),
        // A rewritten part is not looked into again.
        ("ReplaceAll[x, RuleDelayed[x, f[x]]]", "f[x]"),
        (
            "ReplaceRepeated[x, List[RuleDelayed[x, y], RuleDelayed[y, z]]]",
            "z",
        ),
        // `Rule` evaluates its right-hand side, `RuleDelayed` after replacing.
        ("ReplaceAll[x, Rule[x, Plus[1, 1]]]", "2"),
        (
            "ReplaceAll[f[2], RuleDelayed[f[Pattern[n, Blank[]]], Times[n, n]]]",
            "4",
        ),
        // A condition on the right-hand side: the match goes on to the next
        // split while it fails.
        (
            "Replace[f[1, 2, 3], RuleDelayed[f[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]]], Condition[g[List[x], List[y]], Equal[Length[List[x]], 2]]]]",
            "g[List[1, 2], List[3]]",
        ),
        (
            "Replace[f[1], RuleDelayed[f[Pattern[x, Blank[]]], Condition[x, False]]]",
            "f[1]",
        ),
        // Not rules: left unevaluated.
        ("ReplaceAll[x, y]", "ReplaceAll[x, y]"),
        (
            "MatchQ[x, Pattern[1, Blank[]]]",
            "MatchQ[x, Pattern[1, Blank[]]]",
        ),
    ]);
    // A rule that never stops: each round counts against the iterations.
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.iterations = 64;
    ctx.set_limits(limits);
    assert_eq!(
        ctx.eval(&read("ReplaceRepeated[x, RuleDelayed[x, f[x]]]")),
        Err(EvalError::IterationLimit)
    );
}

#[test]
fn functions() {
    check(&[
        ("Function[List[x], Power[x, 2]][3]", "9"),
        ("Function[List[x, y], Plus[x, y]][1, 2]", "3"),
        ("Function[x, Times[2, x]][5]", "10"),
        ("Function[List[], 7][]", "7"),
        ("Function[Plus[Slot[1], Slot[2]]][1, 2]", "3"),
        ("Function[Times[Slot[], 2]][4]", "8"),
        // Too few arguments: left alone.
        ("Function[List[x, y], x][1]", "Function[List[x, y], x][1]"),
        ("Function[Slot[2]][1]", "Function[Slot[2]][1]"),
        // A nested function binds its own parameter.
        (
            "Function[List[x], Function[List[x], x]][1]",
            "Function[List[x], x]",
        ),
        (
            "Function[List[x], Function[List[y], List[x, y]]][1][2]",
            "List[1, 2]",
        ),
        // Slots of a nested pure function are its own; a function with named
        // parameters has none.
        ("Function[Function[Slot[1]]][1]", "Function[Slot[1]]"),
        (
            "Function[Function[List[y], List[Slot[1], y]]][1][2]",
            "List[1, 2]",
        ),
        (
            "Map[Function[List[x], Power[x, 2]], List[1, 2, 3]]",
            "List[1, 4, 9]",
        ),
        ("Apply[Function[List[x, y], Times[x, y]], List[3, 4]]", "12"),
    ]);
}

/// `f[a_] :> x -> a x`: instantiating the rule substitutes `a` and leaves the
/// lambda's own `x` alone, also when a pattern variable has the same name.
#[test]
fn closures_keep_their_bound_parameter() {
    check(&[
        (
            "Replace[f[2], RuleDelayed[f[Pattern[a, Blank[]]], Function[List[x], Times[a, x]]]]",
            "Function[List[x], Times[2, x]]",
        ),
        (
            "Replace[f[2], RuleDelayed[f[Pattern[a, Blank[]]], Function[List[x], Times[a, x]]]][3]",
            "6",
        ),
        (
            "Replace[f[2], RuleDelayed[f[Pattern[x, Blank[]]], Function[List[x], x]]]",
            "Function[List[x], x]",
        ),
        (
            "Replace[f[2], RuleDelayed[f[Pattern[x, Blank[]]], List[x, Function[x, x]]]]",
            "List[2, Function[x, x]]",
        ),
    ]);
}

/// A rule compiled once gives what compiling it again gives, in contexts
/// with other values: its condition reads the context it runs in.
#[test]
fn a_rule_compiled_once_works_in_any_context() {
    let rule = Rule::new(
        &read("f[Pattern[x, Blank[]]]"),
        &read("Condition[g[x], Greater[x, limit]]"),
    )
    .unwrap();
    let e = read("List[f[1], f[5], f[10]]");
    let mut contexts = vec![Context::new()];
    for limit in [0, 4, 9, 20] {
        let mut ctx = Context::new();
        ctx.set(Symbol::new("limit"), Expr::from(limit)).unwrap();
        contexts.push(ctx);
    }
    for ctx in &contexts {
        let again = Rule::new(
            &read("f[Pattern[x, Blank[]]]"),
            &read("Condition[g[x], Greater[x, limit]]"),
        )
        .unwrap();
        assert_eq!(
            ctx.replace_all(&e, std::slice::from_ref(&rule)),
            ctx.replace_all(&e, &[again])
        );
    }
    let mut ctx = Context::new();
    ctx.set(Symbol::new("limit"), Expr::from(4)).unwrap();
    assert_eq!(
        format!("{:?}", ctx.replace_all(&e, &[rule]).unwrap()),
        "List[f[1], g[5], g[10]]"
    );
}

#[test]
fn pattern_errors() {
    let cases = [
        ("Pattern[1, Blank[]]", K::NotAName),
        ("f[Pattern[x]]", K::Arity),
        ("Blank[a, b]", K::Arity),
        ("Condition[Blank[]]", K::Arity),
        ("Optional[Blank[], 1, 2]", K::Arity),
        ("HoldPattern[a, b]", K::Arity),
    ];
    for (text, kind) in cases {
        let error = Pattern::new(&read(text)).expect_err(text);
        assert_eq!(error.kind, kind, "{text}");
        let _ = error.to_string();
    }
    let p = Pattern::new(&read(
        "f[Pattern[x, Blank[]], g[Pattern[y, Blank[]]], Pattern[x, Blank[]]]",
    ))
    .unwrap();
    assert_eq!(p.variables(), [Symbol::new("x"), Symbol::new("y")]);
    assert_eq!(p.heads(), [Symbol::new("f"), Symbol::new("g")]);
}

// A reference matcher: every way of matching, by brute force, in the order
// the matcher tries them.

type Env = HashMap<Symbol, Expr>;

fn head_symbol(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn bind(env: &Env, x: Symbol, value: Expr) -> Option<Env> {
    match env.get(&x) {
        Some(bound) if *bound != value => None,
        _ => {
            let mut env = env.clone();
            env.insert(x, value);
            Some(env)
        }
    }
}

fn reference(p: &Expr, s: &Expr, env: &Env) -> Vec<Env> {
    let View::Apply { head, args } = p.view() else {
        return if p == s { vec![env.clone()] } else { vec![] };
    };
    let blank = |args: &[Expr]| match args {
        [] => vec![env.clone()],
        [h] if s.head() == Some(h) => vec![env.clone()],
        _ => vec![],
    };
    match head.as_symbol() {
        Some(Symbol::BLANK | Symbol::BLANK_SEQUENCE | Symbol::BLANK_NULL_SEQUENCE) => blank(args),
        Some(Symbol::PATTERN) => {
            let x = args[0].as_symbol().unwrap();
            bind(env, x, s.clone())
                .map(|env| reference(&args[1], s, &env))
                .unwrap_or_default()
        }
        Some(Symbol::ALTERNATIVES) => args.iter().flat_map(|a| reference(a, s, env)).collect(),
        _ => {
            let View::Apply { head: sh, args: sa } = s.view() else {
                return vec![];
            };
            reference(head, sh, env)
                .iter()
                .flat_map(|env| reference_args(args, sa, env))
                .collect()
        }
    }
}

/// `(variable, least length, head)` of a sequence pattern.
fn reference_sequence(p: &Expr) -> Option<(Option<Symbol>, usize, Option<Expr>)> {
    let (x, blank) = if head_symbol(p) == Some(Symbol::PATTERN) {
        (p.args()[0].as_symbol(), &p.args()[1])
    } else {
        (None, p)
    };
    let min = match head_symbol(blank) {
        Some(Symbol::BLANK_SEQUENCE) => 1,
        Some(Symbol::BLANK_NULL_SEQUENCE) => 0,
        _ => return None,
    };
    Some((x, min, blank.args().first().cloned()))
}

fn reference_args(ps: &[Expr], ss: &[Expr], env: &Env) -> Vec<Env> {
    let Some((first, rest)) = ps.split_first() else {
        return if ss.is_empty() {
            vec![env.clone()]
        } else {
            vec![]
        };
    };
    if let Some((x, min, head)) = reference_sequence(first) {
        let mut out = Vec::new();
        for len in min..=ss.len() {
            if head
                .as_ref()
                .is_some_and(|h| ss[..len].iter().any(|s| s.head() != Some(h)))
            {
                continue;
            }
            let run = Expr::apply(&Symbol::SEQUENCE.into(), &ss[..len]);
            let env = match x {
                Some(x) => bind(env, x, run),
                None => Some(env.clone()),
            };
            if let Some(env) = env {
                out.extend(reference_args(rest, &ss[len..], &env));
            }
        }
        return out;
    }
    let Some((s, more)) = ss.split_first() else {
        return vec![];
    };
    reference(first, s, env)
        .iter()
        .flat_map(|env| reference_args(rest, more, env))
        .collect()
}

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn app(head: &str, args: impl AsRef<[Expr]>) -> Expr {
    Expr::apply(&sym(head), args.as_ref())
}

/// Random patterns over `Blank`, sequences, variables `x` and `y` (so
/// names repeat), `Alternatives` and applications.
fn pattern() -> impl Strategy<Value = Expr> {
    let blank = prop::option::of(prop::sample::select(vec!["f", "g"]))
        .prop_map(|h| app("Blank", h.map(sym).into_iter().collect::<Vec<_>>()));
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "b"]).prop_map(sym),
        (1i64..3).prop_map(Expr::from),
        blank.clone(),
        (prop::sample::select(vec!["x", "y"]), blank)
            .prop_map(|(x, b)| app("Pattern", vec![sym(x), b])),
    ];
    leaf.prop_recursive(3, 16, 4, |inner| {
        let sequence = (
            prop::option::of(prop::sample::select(vec!["x", "y"])),
            prop::sample::select(vec!["BlankSequence", "BlankNullSequence"]),
        )
            .prop_map(|(x, blank)| {
                let blank = app(blank, vec![]);
                match x {
                    Some(x) => app("Pattern", vec![sym(x), blank]),
                    None => blank,
                }
            });
        prop_oneof![
            (
                prop::sample::select(vec!["f", "g"]),
                prop::collection::vec(prop_oneof![inner.clone(), sequence], 0..4),
            )
                .prop_map(|(h, args)| app(h, args)),
            prop::collection::vec(inner, 2..3).prop_map(|alts| app("Alternatives", alts)),
        ]
    })
}

/// Random subjects over the same symbols.
fn subject() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "b"]).prop_map(sym),
        (1i64..3).prop_map(Expr::from),
    ];
    leaf.prop_recursive(3, 16, 4, |inner| {
        (
            prop::sample::select(vec!["f", "g"]),
            prop::collection::vec(inner, 0..4),
        )
            .prop_map(|(h, args)| app(h, args))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// The matcher finds a match exactly when the reference does, and the
    /// same first one.
    #[test]
    fn matches_like_the_reference(p in pattern(), s in subject()) {
        let ctx = Context::new();
        let compiled = Pattern::new(&p).unwrap();
        let expected = reference(&p, &s, &Env::new());
        prop_assert_eq!(ctx.matches(&s, &compiled).unwrap(), !expected.is_empty(), "{:?} on {:?}", p, s);
        if let Some(first) = expected.first() {
            let vars: Vec<Expr> = ["x", "y"].map(sym).to_vec();
            let rule = Rule::new(&p, &app("List", vars.clone())).unwrap();
            let found = ctx.eval(&Expr::apply(
                &Symbol::REPLACE.into(),
                &[s.clone(), Expr::apply(&Symbol::RULE_DELAYED.into(), &[p.clone(), app("List", vars.clone())])],
            )).unwrap();
            // The pattern's own variables take their bindings, or the empty
            // sequence when the alternative taken leaves them out; other
            // symbols stay.
            let expected: Vec<Expr> = vars
                .iter()
                .map(|v| {
                    let x = v.as_symbol().unwrap();
                    if !compiled.variables().contains(&x) {
                        return v.clone();
                    }
                    first
                        .get(&x)
                        .cloned()
                        .unwrap_or_else(|| Expr::apply(&Symbol::SEQUENCE.into(), &[]))
                })
                .collect();
            let expected = ctx.eval(&app("List", expected)).unwrap();
            prop_assert_eq!(&found, &expected, "{:?} on {:?}", p, s);
            prop_assert_eq!(ctx.replace_all(&s, &[rule]).unwrap() == expected, true);
        }
    }

    /// The filter of literal symbols only rejects subjects that cannot match.
    #[test]
    fn the_symbol_filter_changes_nothing(p in pattern(), s in subject()) {
        let ctx = Context::new();
        let compiled = Pattern::new(&p).unwrap();
        let mut unfiltered = compiled.clone();
        unfiltered.required = 0;
        prop_assert_eq!(ctx.matches(&s, &compiled).unwrap(), ctx.matches(&s, &unfiltered).unwrap());
    }
}
