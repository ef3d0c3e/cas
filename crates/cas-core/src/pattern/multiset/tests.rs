//! Tests of matching under `Flat`, `Orderless` and `OneIdentity` (CAS-57),
//! against tables and against a brute-force oracle.

use proptest::prelude::*;

use crate::{Attributes, Context, EvalError, Expr, Limits, Pattern, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn sym(name: &str) -> Symbol {
    Symbol::new(name)
}

/// A context where `o` is `Orderless`, `fl` `Flat`, `ac` both, `aci` both
/// and `OneIdentity`, `oi` `Orderless` and `OneIdentity`, and `pl` plain.
fn context() -> Context {
    let mut ctx = Context::new();
    let a = |list: &[Attributes]| list.iter().fold(Attributes::NONE, |x, &y| x | y);
    let (f, o, i) = (
        Attributes::FLAT,
        Attributes::ORDERLESS,
        Attributes::ONE_IDENTITY,
    );
    for (name, attributes) in [
        ("o", a(&[o])),
        ("fl", a(&[f])),
        ("ac", a(&[f, o])),
        ("aci", a(&[f, o, i])),
        ("oi", a(&[o, i])),
    ] {
        ctx.set_attributes(sym(name), attributes).unwrap();
    }
    ctx
}

/// The first match of `pattern` against `subject`, both in full form, the
/// subject normalised by evaluation: the bindings by name, in the order of
/// the variables.
fn first(
    ctx: &Context,
    pattern: &str,
    subject: &str,
) -> Result<Option<Vec<(String, String)>>, EvalError> {
    let p = Pattern::with_context(&read(pattern), ctx).unwrap();
    let subject = ctx.eval(&read(subject)).unwrap();
    first_match(ctx, &p, &subject)
}

fn first_match(
    ctx: &Context,
    p: &Pattern,
    subject: &Expr,
) -> Result<Option<Vec<(String, String)>>, EvalError> {
    let mut found = None;
    let mut ev = crate::eval::Eval::new(ctx, ctx.limits().recursion);
    p.match_with(&mut ev, subject, &mut |_, bound| {
        found = Some(
            p.variables()
                .iter()
                .zip(bound)
                .map(|(x, v)| (x.name().to_owned(), format!("{v:?}")))
                .collect(),
        );
        Ok(true)
    })?;
    Ok(found)
}

/// `(pattern, subject, bindings)` in full form, the bindings as
/// `x=value` separated by spaces, or `None` for no match.
fn check(ctx: &Context, cases: &[(&str, &str, Option<&str>)]) {
    let mut wrong = Vec::new();
    for &(pattern, subject, expected) in cases {
        let found = first(ctx, pattern, subject).map(|m| {
            m.map(|b| {
                b.iter()
                    .map(|(x, v)| format!("{x}={v}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
        });
        if found.as_ref().ok().map(Option::as_deref) != Some(expected) {
            wrong.push(format!(
                "{pattern} vs {subject}: {found:?}, expected {expected:?}"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn orderless_arguments_match_in_any_order() {
    let ctx = context();
    check(
        &ctx,
        &[
            ("o[Pattern[x, Blank[]], a]", "o[b, a]", Some("x=b")),
            ("o[a, Pattern[x, Blank[]]]", "o[a, b]", Some("x=b")),
            ("o[g[Pattern[x, Blank[]]], b]", "o[b, g[1]]", Some("x=1")),
            (
                "o[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "o[b, a]",
                Some("x=a y=b"),
            ),
            (
                "o[Pattern[x, Blank[]], Pattern[x, Blank[]]]",
                "o[a, a]",
                Some("x=a"),
            ),
            (
                "o[Pattern[x, Blank[]], Pattern[x, Blank[]]]",
                "o[a, b]",
                None,
            ),
            // A fixed head first: it binds `x`, which the blank then must be.
            (
                "o[Pattern[x, Blank[]], g[Pattern[x, Blank[]]]]",
                "o[g[1], 1]",
                Some("x=1"),
            ),
            (
                "o[Pattern[x, Blank[]], g[Pattern[x, Blank[]]]]",
                "o[g[1], 2]",
                None,
            ),
            (
                "o[Pattern[x, BlankSequence[]], a]",
                "o[c, a, b]",
                Some("x=Sequence[b, c]"),
            ),
            (
                "o[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]]]",
                "o[a, b, c]",
                Some("x=Sequence[a] y=Sequence[b, c]"),
            ),
            (
                "o[Pattern[x, BlankNullSequence[]], a, b]",
                "o[b, a]",
                Some("x=Sequence[]"),
            ),
            (
                "o[Pattern[x, Blank[g]], Pattern[y, Blank[]]]",
                "o[g[2], a, g[1]]",
                None,
            ),
            (
                "o[Pattern[x, BlankSequence[g]], Pattern[y, Blank[]]]",
                "o[g[2], a, g[1]]",
                Some("x=Sequence[g[1], g[2]] y=a"),
            ),
            (
                "o[a, Optional[Pattern[x, Blank[]], z]]",
                "o[a]",
                Some("x=z"),
            ),
            (
                "o[a, Optional[Pattern[x, Blank[]], z]]",
                "o[b, a]",
                Some("x=b"),
            ),
            // Without `Flat`, one node takes one argument.
            (
                "o[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "o[a, b, c]",
                None,
            ),
        ],
    );
}

#[test]
fn flat_arguments_group() {
    let ctx = context();
    check(
        &ctx,
        &[
            (
                "fl[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "fl[a, b]",
                Some("x=a y=b"),
            ),
            (
                "fl[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "fl[a, b, c]",
                Some("x=a y=fl[b, c]"),
            ),
            (
                "fl[Pattern[x, Blank[]], c]",
                "fl[a, b, c]",
                Some("x=fl[a, b]"),
            ),
            ("fl[b, Pattern[x, Blank[]]]", "fl[a, b, c]", None),
            (
                "fl[Pattern[x, Blank[]], Pattern[x, Blank[]]]",
                "fl[a, b, a, b]",
                Some("x=fl[a, b]"),
            ),
            // Sequences bind sequences, not groups.
            (
                "fl[Pattern[x, BlankSequence[]], c]",
                "fl[a, b, c]",
                Some("x=Sequence[a, b]"),
            ),
            // A blank with another head takes one argument only.
            (
                "fl[Pattern[x, Blank[g]], Pattern[y, Blank[]]]",
                "fl[g[1], a, b]",
                Some("x=g[1] y=fl[a, b]"),
            ),
            ("fl[Pattern[x, Blank[g]]]", "fl[g[1], g[2]]", None),
            // Order still counts without `Orderless`.
            ("fl[c, Pattern[x, Blank[]]]", "fl[a, b, c]", None),
        ],
    );
}

#[test]
fn flat_orderless_arguments_group_as_multisets() {
    let ctx = context();
    check(
        &ctx,
        &[
            (
                "ac[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "ac[a, b, c]",
                Some("x=a y=ac[b, c]"),
            ),
            (
                "ac[Pattern[x, Blank[]], b]",
                "ac[c, b, a]",
                Some("x=ac[a, c]"),
            ),
            (
                "ac[Pattern[x, Blank[]], Pattern[x, Blank[]]]",
                "ac[a, b, a, b]",
                Some("x=ac[a, b]"),
            ),
            (
                "ac[Pattern[x, Blank[]], Pattern[x, Blank[]]]",
                "ac[a, a, b]",
                None,
            ),
            (
                "ac[g[Pattern[x, Blank[]]], Pattern[y, Blank[]]]",
                "ac[a, g[1], b]",
                Some("x=1 y=ac[a, b]"),
            ),
            (
                "ac[Pattern[x, Blank[]], Pattern[y, BlankNullSequence[]]]",
                "ac[a, b]",
                Some("x=a y=Sequence[b]"),
            ),
            // The last node takes all that is left.
            (
                "ac[a, Pattern[x, Blank[]]]",
                "ac[a, b, c, d]",
                Some("x=ac[b, c, d]"),
            ),
        ],
    );
}

#[test]
fn one_identity_lets_a_lone_expression_stand_for_the_application() {
    // The builtins: `Plus` and `Times` are `Flat`, `Orderless` and
    // `OneIdentity`, `Power` is `OneIdentity`.
    let ctx = Context::new();
    check(
        &ctx,
        &[
            (
                "Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]]",
                "y",
                Some("x=y n=1"),
            ),
            (
                "Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]]",
                "Power[y, 3]",
                Some("x=y n=3"),
            ),
            (
                "Times[Optional[Pattern[c, Blank[]]], Pattern[x, Blank[]]]",
                "y",
                Some("c=1 x=y"),
            ),
            (
                "Times[Optional[Pattern[c, Blank[]]], Pattern[x, Blank[]]]",
                "Times[3, y]",
                Some("c=3 x=y"),
            ),
            (
                "Times[Optional[Pattern[c, Blank[]]], Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]]]",
                "Times[3, Power[y, 2]]",
                Some("c=3 x=y n=2"),
            ),
            (
                "Times[Optional[Pattern[c, Blank[]]], Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]]]",
                "Power[y, 2]",
                Some("c=1 x=y n=2"),
            ),
            (
                "Times[Optional[Pattern[c, Blank[]]], Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]]]",
                "y",
                Some("c=1 x=y n=1"),
            ),
            (
                "Plus[Pattern[x, Blank[]], Optional[Pattern[a, Blank[]]]]",
                "y",
                Some("x=y a=0"),
            ),
            // Two arguments that must be there need two.
            ("Plus[Pattern[x, Blank[]], Pattern[y, Blank[]]]", "y", None),
            // An application of the head matches only as one.
            (
                "Power[Pattern[x, Blank[g]], Optional[Pattern[n, Blank[]]]]",
                "Power[y, 2]",
                None,
            ),
        ],
    );
    // A user head with `OneIdentity` and a default.
    let mut ctx = context();
    ctx.set_default(sym("aci"), read("e")).unwrap();
    check(
        &ctx,
        &[
            (
                "aci[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]]]]",
                "q",
                Some("x=q y=e"),
            ),
            (
                "oi[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], w]]",
                "q",
                Some("x=q y=w"),
            ),
            // Without `OneIdentity`, no.
            (
                "ac[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], w]]",
                "q",
                None,
            ),
        ],
    );
}

#[test]
fn lists_in_times_keep_their_order() {
    let ctx = Context::new();
    check(
        &ctx,
        &[
            (
                "Times[Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "Times[List[1], List[2]]",
                Some("x=List[1] y=List[2]"),
            ),
            (
                "Times[Pattern[y, Blank[List]], Pattern[x, Blank[List]]]",
                "Times[List[1], List[2]]",
                Some("y=List[1] x=List[2]"),
            ),
            (
                "Times[List[2], Pattern[x, Blank[]]]",
                "Times[List[1], List[2]]",
                None,
            ),
            (
                "Times[Pattern[x, Blank[]], List[2]]",
                "Times[List[1], List[2]]",
                Some("x=List[1]"),
            ),
            (
                "Times[Pattern[c, Blank[]], Pattern[x, Blank[]], Pattern[y, Blank[]]]",
                "Times[2, List[1], List[2]]",
                Some("c=2 x=List[1] y=List[2]"),
            ),
        ],
    );
    // A pattern that matches only lists keeps its place in a product, as a
    // list does, so a rule for `y x` binds `y` to the left factor.
    let mut ctx = Context::new();
    ctx.define_delayed(
        &read("f[Times[Pattern[y, Blank[List]], Pattern[x, Blank[List]]]]"),
        &read("left[y]"),
    )
    .unwrap();
    assert_eq!(
        ctx.eval(&read("f[Times[List[1], List[2]]]")).unwrap(),
        read("left[List[1]]")
    );
}

#[test]
fn rules_on_attributed_heads_dispatch_through_the_filter() {
    let mut ctx = context();
    for (lhs, rhs) in [
        ("fl[Pattern[x, Blank[]], Pattern[y, Blank[]]]", "pair[x, y]"),
        ("o[g[Pattern[x, Blank[]]], b]", "first[x]"),
        // Stored as `o[b, x_]`: a match need not have `b` first.
        ("o[b, Pattern[x, Blank[]]]", "second[x]"),
        (
            "aci[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], w]]",
            "lone[x, y]",
        ),
    ] {
        ctx.define_delayed(&read(lhs), &read(rhs)).unwrap();
    }
    for (input, expected) in [
        // The group `fl[b, c]` is a pair too, once evaluated.
        ("fl[a, b, c]", "pair[a, pair[b, c]]"),
        ("o[g[1], b]", "first[1]"),
        ("o[b, g[1]]", "first[1]"),
        ("o[b, a]", "second[a]"),
        // The stored left-hand side is sorted: `y_.` comes first.
        ("aci[q, r]", "lone[r, q]"),
    ] {
        assert_eq!(ctx.eval(&read(input)).unwrap(), read(expected), "{input}");
    }
}

#[test]
fn the_budget_stops_a_match_deterministically() {
    // Three sequences over twelve arguments and a condition that never
    // holds: every way of splitting is tried.
    let mut ctx = context();
    let pattern = read(
        "Condition[o[Pattern[x, BlankSequence[]], Pattern[y, BlankSequence[]], Pattern[z, BlankSequence[]]], False]",
    );
    let args: Vec<Expr> = (0..12).map(Expr::from).collect();
    let subject = Expr::apply(&sym("o").into(), &args);
    let p = Pattern::with_context(&pattern, &ctx).unwrap();
    assert_eq!(first_match(&ctx, &p, &subject), Err(EvalError::MatchLimit));
    assert_eq!(first_match(&ctx, &p, &subject), Err(EvalError::MatchLimit));
    // With room for every split, no match.
    let mut limits = ctx.limits();
    limits.match_steps = 10_000_000;
    ctx.set_limits(limits);
    let p = Pattern::with_context(&pattern, &ctx).unwrap();
    assert_eq!(first_match(&ctx, &p, &subject), Ok(None));
    // A tight budget stops a small match too, and the error says so.
    ctx.set_limits(Limits {
        match_steps: 3,
        ..Limits::default()
    });
    assert_eq!(
        first(
            &ctx,
            "o[Pattern[x, Blank[]], Pattern[y, Blank[]], c]",
            "o[a, b, c]"
        ),
        Err(EvalError::MatchLimit)
    );
    assert_eq!(
        EvalError::MatchLimit.to_string(),
        "match step limit exceeded"
    );
}

#[test]
fn pruning_keeps_easy_matches_within_a_small_budget() {
    let mut ctx = context();
    let mut limits = ctx.limits();
    limits.match_steps = 50;
    ctx.set_limits(limits);
    let many = |head: &str, extra: &[&str]| -> Expr {
        let mut args: Vec<Expr> = (0..1000).map(Expr::from).collect();
        args.extend(extra.iter().map(|e| read(e)));
        ctx.eval(&Expr::apply(&sym(head).into(), &args)).unwrap()
    };
    // Three nodes cannot take a thousand arguments: no search at all.
    let p = Pattern::with_context(
        &read("o[Pattern[x, Blank[]], Pattern[y, Blank[]], Pattern[z, Blank[]]]"),
        &ctx,
    )
    .unwrap();
    assert_eq!(first_match(&ctx, &p, &many("o", &[])), Ok(None));
    // `x` is bound by `g[x_]` before the blank: the blank takes the equal
    // subject at once, though it is the last of a thousand.
    let p = Pattern::with_context(
        &read("o[g[Pattern[x, Blank[]]], Pattern[x, Blank[]], Pattern[r, BlankNullSequence[]]]"),
        &ctx,
    )
    .unwrap();
    let found = first_match(&ctx, &p, &many("o", &["g[999]"])).unwrap();
    assert_eq!(found.map(|b| b[0].1.clone()), Some("999".to_owned()));
}

/// Bindings by name, in the order they were made.
type Bindings = Vec<(Symbol, Expr)>;

/// A brute-force matcher over pattern expressions: every match, in the
/// order the matcher defines, with no pruning, no shortcut for known values
/// and every copy of equal subjects tried.
struct Oracle<'c> {
    ctx: &'c Context,
}

fn bind(b: &Bindings, x: Symbol, value: &Expr) -> Option<Bindings> {
    if let Some((_, v)) = b.iter().find(|(y, _)| *y == x) {
        return (v == value).then(|| b.clone());
    }
    let mut b = b.clone();
    b.push((x, value.clone()));
    Some(b)
}

/// The bindings of an optional argument left out: its variable takes the
/// default.
fn bind_missing(inner: &Expr, d: &Expr, b: &Bindings) -> Option<Bindings> {
    match (head_symbol(inner), inner.args()) {
        (Some(Symbol::PATTERN), [x, _]) => bind(b, x.as_symbol().unwrap(), d),
        _ => Some(b.clone()),
    }
}

fn head_symbol(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// `(variable, least, head)` for a sequence pattern.
fn sequence(p: &Expr) -> Option<(Option<Symbol>, usize, Option<Expr>)> {
    match (head_symbol(p), p.args()) {
        (Some(Symbol::BLANK_SEQUENCE), args) => Some((None, 1, args.first().cloned())),
        (Some(Symbol::BLANK_NULL_SEQUENCE), args) => Some((None, 0, args.first().cloned())),
        (Some(Symbol::PATTERN), [x, q]) => sequence(q)
            .filter(|(v, ..)| v.is_none())
            .map(|(_, least, h)| (x.as_symbol(), least, h)),
        _ => None,
    }
}

/// The inner pattern and default of `Optional[q, d]`.
fn optional(p: &Expr) -> Option<(&Expr, &Expr)> {
    match (head_symbol(p), p.args()) {
        (Some(Symbol::OPTIONAL), [q, d]) => Some((q, d)),
        _ => None,
    }
}

impl Oracle<'_> {
    fn lone(&self, h: Symbol, qs: &[Expr]) -> bool {
        let least: usize = qs
            .iter()
            .map(|q| match sequence(q) {
                Some((_, least, _)) => least,
                None if optional(q).is_some() => 0,
                None => 1,
            })
            .sum();
        self.ctx.attributes(h).contains(Attributes::ONE_IDENTITY) && least <= 1 && !qs.is_empty()
    }

    fn fixed_head(&self, q: &Expr) -> Option<Expr> {
        if !q.has_pattern() {
            return q.head().cloned();
        }
        match (head_symbol(q), q.args()) {
            (Some(Symbol::BLANK), [h]) => Some(h.clone()),
            (Some(Symbol::PATTERN), [_, r]) => self.fixed_head(r),
            (Some(s), _) if !crate::builtin::is_pattern_head(s) && !self.lone(s, q.args()) => {
                q.head().cloned()
            }
            _ => None,
        }
    }

    fn category(&self, q: &Expr) -> u8 {
        if sequence(q).is_some() {
            return 4;
        }
        if !q.has_pattern() {
            return 0;
        }
        if let Some((r, _)) = optional(q) {
            return self.category(r);
        }
        if self.stands_alone(q) {
            3
        } else if self.fixed_head(q).is_some() {
            1
        } else {
            2
        }
    }

    fn stands_alone(&self, q: &Expr) -> bool {
        match (head_symbol(q), q.args()) {
            (Some(Symbol::PATTERN), [_, r]) => self.stands_alone(r),
            (Some(s), args) if !crate::builtin::is_pattern_head(s) => self.lone(s, args),
            _ => false,
        }
    }

    /// Tries the matches of `p` against `s` that extend `b`, in order,
    /// until `found` takes one; whether it did.
    fn matches(&self, p: &Expr, s: &Expr, b: &Bindings, found: Found<'_>) -> bool {
        if !p.has_pattern() {
            return p == s && found(b);
        }
        let args = p.args();
        match head_symbol(p) {
            Some(Symbol::BLANK | Symbol::BLANK_SEQUENCE | Symbol::BLANK_NULL_SEQUENCE) => {
                args.first().is_none_or(|h| s.head() == Some(h)) && found(b)
            }
            Some(Symbol::PATTERN) => match bind(b, args[0].as_symbol().unwrap(), s) {
                Some(b) => self.matches(&args[1], s, &b, found),
                None => false,
            },
            Some(Symbol::OPTIONAL) => self.matches(&args[0], s, b, found),
            _ => {
                let h = head_symbol(p);
                if let (Some(sh), false) = (s.head(), s.is_atom()) {
                    if self.matches(p.head().unwrap(), sh, b, &mut |b| {
                        self.arguments(h, args, s.args(), b, found)
                    }) {
                        return true;
                    }
                    if h.is_some_and(|h| sh.as_symbol() == Some(h)) {
                        return false;
                    }
                }
                match h.filter(|&h| self.lone(h, args)) {
                    Some(h) => self.arguments(Some(h), args, std::slice::from_ref(s), b, found),
                    None => false,
                }
            }
        }
    }

    fn arguments(
        &self,
        h: Option<Symbol>,
        qs: &[Expr],
        ss: &[Expr],
        b: &Bindings,
        found: Found<'_>,
    ) -> bool {
        let attributes = h.map_or(Attributes::NONE, |h| self.ctx.attributes(h));
        let flat = attributes.contains(Attributes::FLAT);
        if attributes.contains(Attributes::ORDERLESS) {
            let mut order: Vec<usize> = (0..qs.len()).collect();
            order.sort_by_key(|&i| self.category(&qs[i]));
            let search = Search {
                h: h.unwrap(),
                flat,
                qs,
                ss,
            };
            self.multiset(&search, &order, &mut vec![None; ss.len()], b, found)
        } else {
            self.sequential(h, flat, qs, ss, b, found)
        }
    }

    #[expect(
        clippy::many_single_char_names,
        reason = "the matcher's names: pattern q, subject s, bindings b, head h"
    )]
    fn sequential(
        &self,
        h: Option<Symbol>,
        flat: bool,
        qs: &[Expr],
        ss: &[Expr],
        b: &Bindings,
        found: Found<'_>,
    ) -> bool {
        let Some((q, more)) = qs.split_first() else {
            return ss.is_empty() && found(b);
        };
        if let Some((x, least, head)) = sequence(q) {
            for len in least..=ss.len() {
                let run = &ss[..len];
                if head
                    .as_ref()
                    .is_some_and(|h| run.iter().any(|s| s.head() != Some(h)))
                {
                    continue;
                }
                let value = Expr::apply(&Symbol::SEQUENCE.into(), run);
                if let Some(b) = x.map_or(Some(b.clone()), |x| bind(b, x, &value)) {
                    if self.sequential(h, flat, more, &ss[len..], &b, found) {
                        return true;
                    }
                }
            }
            return false;
        }
        let (inner, default) = match optional(q) {
            Some((r, d)) => (r, Some(d)),
            None => (q, None),
        };
        let most = if flat { ss.len() } else { ss.len().min(1) };
        for len in 1..=most {
            let item = match &ss[..len] {
                [one] => one.clone(),
                run => Expr::apply(&h.unwrap().into(), run),
            };
            if self.matches(inner, &item, b, &mut |b| {
                self.sequential(h, flat, more, &ss[len..], b, found)
            }) {
                return true;
            }
        }
        match default.and_then(|d| bind_missing(inner, d, b)) {
            Some(b) => self.sequential(h, flat, more, ss, &b, found),
            None => false,
        }
    }

    /// The nodes of `order` take subjects in turn: every subset of the free
    /// ones, by size and in increasing order, each copy of equal subjects
    /// tried.
    #[expect(
        clippy::many_single_char_names,
        reason = "the matcher's names: pattern q, subject s, bindings b, node i"
    )]
    fn multiset(
        &self,
        search: &Search<'_>,
        order: &[usize],
        taker: &mut Vec<Option<usize>>,
        b: &Bindings,
        found: Found<'_>,
    ) -> bool {
        let Some((&i, later)) = order.split_first() else {
            return taker.iter().all(Option::is_some)
                && (search.h != Symbol::TIMES || lists_in_order(search.ss, taker))
                && found(b);
        };
        let free: Vec<usize> = (0..search.ss.len())
            .filter(|&j| taker[j].is_none())
            .collect();
        let q = &search.qs[i];
        // Gives `chosen` to node `i`, then goes on.
        let mut next = |chosen: &[usize], b: &Bindings, taker: &mut Vec<Option<usize>>| {
            for &j in chosen {
                taker[j] = Some(i);
            }
            let done = self.multiset(search, later, taker, b, found);
            for &j in chosen {
                taker[j] = None;
            }
            done
        };
        if let Some((x, least, head)) = sequence(q) {
            for size in least..=free.len() {
                for chosen in combinations(&free, size) {
                    let run: Vec<Expr> = chosen.iter().map(|&j| search.ss[j].clone()).collect();
                    if head
                        .as_ref()
                        .is_some_and(|h| run.iter().any(|s| s.head() != Some(h)))
                    {
                        continue;
                    }
                    let value = Expr::apply(&Symbol::SEQUENCE.into(), &run);
                    if let Some(b) = x.map_or(Some(b.clone()), |x| bind(b, x, &value)) {
                        if next(&chosen, &b, taker) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
        let (inner, default) = match optional(q) {
            Some((r, d)) => (r, Some(d)),
            None => (q, None),
        };
        let most = if search.flat {
            free.len()
        } else {
            free.len().min(1)
        };
        for size in 1..=most {
            for chosen in combinations(&free, size) {
                let item = if let [j] = chosen.as_slice() {
                    search.ss[*j].clone()
                } else {
                    let group: Vec<Expr> = chosen.iter().map(|&j| search.ss[j].clone()).collect();
                    Expr::apply(&search.h.into(), &group)
                };
                let owners = &mut *taker;
                if self.matches(inner, &item, b, &mut |b| next(&chosen, b, owners)) {
                    return true;
                }
            }
        }
        match default.and_then(|d| bind_missing(inner, d, b)) {
            Some(b) => next(&[], &b, taker),
            None => false,
        }
    }
}

/// What receives each complete match, and says whether to stop there.
type Found<'a> = &'a mut dyn FnMut(&Bindings) -> bool;

/// What a multiset search matches: the arguments `qs` of a head `h` against
/// the subjects `ss`.
struct Search<'a> {
    h: Symbol,
    flat: bool,
    qs: &'a [Expr],
    ss: &'a [Expr],
}

fn lists_in_order(ss: &[Expr], taker: &[Option<usize>]) -> bool {
    let takers: Vec<usize> = ss
        .iter()
        .zip(taker)
        .filter(|(s, _)| crate::eval::is_list_like(s))
        .map(|(_, t)| t.unwrap())
        .collect();
    takers.windows(2).all(|w| w[0] <= w[1])
}

/// The `size`-subsets of `items`, in increasing lexicographic order.
fn combinations(items: &[usize], size: usize) -> Vec<Vec<usize>> {
    fn go(items: &[usize], size: usize, prefix: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if size == 0 {
            out.push(prefix.clone());
            return;
        }
        for (k, &x) in items.iter().enumerate() {
            if items.len() - k < size {
                break;
            }
            prefix.push(x);
            go(&items[k + 1..], size - 1, prefix, out);
            prefix.pop();
        }
    }
    let mut out = Vec::new();
    go(items, size, &mut Vec::new(), &mut out);
    out
}

const HEADS: [&str; 6] = ["o", "fl", "ac", "aci", "oi", "pl"];

/// Subjects: atoms with repeats, under the test heads, two levels deep.
fn subject() -> impl Strategy<Value = Expr> {
    let atom = prop::sample::select(vec!["a", "b", "c"]).prop_map(|s| Expr::from(sym(s)));
    atom.prop_recursive(2, 12, 5, |inner| {
        (
            prop::sample::select(HEADS.to_vec()),
            prop::collection::vec(inner, 0..5),
        )
            .prop_map(|(h, args)| Expr::apply(&sym(h).into(), &args))
    })
}

/// The arguments of a pattern: literals, blanks, variables (repeated),
/// fixed heads, sequences, optional arguments and nested applications, among
/// them one that a lone expression may stand for.
fn pattern_arguments() -> impl Strategy<Value = Vec<Expr>> {
    let leaf = prop::sample::select(vec![
        "a",
        "b",
        "Blank[]",
        "aci[Pattern[x, Blank[]], Optional[Pattern[z, Blank[]], b]]",
        "Pattern[x, Blank[]]",
        "Pattern[y, Blank[]]",
        "Pattern[x, Blank[o]]",
        "Pattern[x, BlankSequence[]]",
        "Pattern[y, BlankNullSequence[]]",
        "BlankSequence[]",
        "Optional[Pattern[z, Blank[]], a]",
    ])
    .prop_map(read);
    let argument = leaf.prop_recursive(1, 6, 3, |inner| {
        (
            prop::sample::select(HEADS.to_vec()),
            prop::collection::vec(inner, 0..3),
        )
            .prop_map(|(h, args)| Expr::apply(&sym(h).into(), &args))
    });
    prop::collection::vec(argument, 0..4)
}

/// A pattern and a subject, mostly with the same head so that most cases
/// get past it.
fn case() -> impl Strategy<Value = (Expr, Expr)> {
    let head = prop::sample::select(HEADS.to_vec());
    (
        head.clone(),
        head,
        prop::bool::weighted(0.85),
        pattern_arguments(),
        prop::collection::vec(subject(), 0..6),
    )
        .prop_map(|(h, other, same, pattern, subject)| {
            let s = if same { h } else { other };
            (
                Expr::apply(&sym(h).into(), &pattern),
                Expr::apply(&sym(s).into(), &subject),
            )
        })
}

/// A subject and a pattern made from it: each argument kept, replaced by a
/// blank, a variable (two names, so they may clash), a blank with its head
/// or a sequence, or dropped, so that another node must take it; the
/// pattern's arguments in order or reversed. Most of these match.
fn generalized() -> impl Strategy<Value = (Expr, Expr)> {
    (
        prop::sample::select(HEADS.to_vec()),
        prop::collection::vec(subject(), 0..6),
    )
        .prop_flat_map(|(h, args)| {
            let n = args.len();
            (
                Just(h),
                Just(args),
                prop::collection::vec(0u8..7, n),
                any::<bool>(),
            )
        })
        .prop_map(|(h, args, choices, reverse)| {
            let sequence = read("Pattern[w, BlankSequence[]]");
            let mut pattern: Vec<Expr> = Vec::new();
            for (a, c) in args.iter().zip(choices) {
                match c {
                    0 => pattern.push(a.clone()),
                    1 => pattern.push(read("Blank[]")),
                    2 => pattern.push(read("Pattern[x, Blank[]]")),
                    3 => pattern.push(read("Pattern[y, Blank[]]")),
                    4 if pattern.last() != Some(&sequence) => pattern.push(sequence.clone()),
                    4 | 5 => {}
                    _ => pattern.push(match a.head() {
                        Some(g) => Expr::apply(
                            &Symbol::PATTERN.into(),
                            &[
                                sym("x").into(),
                                Expr::apply(&Symbol::BLANK.into(), std::slice::from_ref(g)),
                            ],
                        ),
                        None => read("Blank[]"),
                    }),
                }
            }
            if reverse {
                pattern.reverse();
            }
            (
                Expr::apply(&sym(h).into(), &pattern),
                Expr::apply(&sym(h).into(), &args),
            )
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// The matcher finds a match exactly when the oracle does, and the same
    /// first one.
    #[test]
    fn the_matcher_agrees_with_the_oracle((p, s) in prop_oneof![case(), generalized()]) {
        let ctx = context();
        let s = ctx.eval(&s).unwrap();
        let compiled = Pattern::with_context(&p, &ctx).unwrap();
        let found = first_match(&ctx, &compiled, &s).unwrap();
        let oracle = Oracle { ctx: &ctx };
        let mut first = None;
        oracle.matches(&p, &s, &Vec::new(), &mut |b| {
            first = Some(b.clone());
            true
        });
        let expected = first.map(|b| {
            compiled
                .variables()
                .iter()
                .map(|&x| {
                    let v = b.iter().find(|(y, _)| *y == x).map_or_else(
                        || Expr::apply(&Symbol::SEQUENCE.into(), &[]),
                        |(_, v)| v.clone(),
                    );
                    (x.name().to_owned(), format!("{v:?}"))
                })
                .collect::<Vec<_>>()
        });
        prop_assert_eq!(found, expected, "{:?} against {:?}", p, s);
    }
}
