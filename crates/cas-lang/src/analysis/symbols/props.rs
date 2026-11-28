//! The index against lowering, and what stays stable under an edit
//! elsewhere (CAS-DOC-60 §4).

use std::collections::BTreeSet;
use std::ops::Range;

use cas_core::{Expr, Symbol};
use proptest::prelude::*;

use super::{ScopeId, ScopeKind, SymbolKind, Symbols, symbols};
use crate::lower::{Spans, lower, pattern_variables, range};
use crate::names::{Names, unescape_name};
use crate::parser::{Parse, parse};

/// A scope that binds names: its kind, start, end and the names it binds.
type Bound = (String, usize, usize, BTreeSet<String>);

/// The user names of `symbols`, unescaped.
fn unescaped(symbols: impl IntoIterator<Item = Symbol>) -> BTreeSet<String> {
    symbols
        .into_iter()
        .map(|s| unescape_name(s.name()).to_owned())
        .collect()
}

/// The `Block`s, rules and `Function`s in `e`, whose spans are `spans` at
/// `path`: a block's locals, a rule's pattern variables, a function's
/// parameters.
fn lowered_scopes(e: &Expr, spans: &Spans, path: &mut Vec<usize>, out: &mut Vec<Bound>) {
    let Some(head) = e.head() else {
        return;
    };
    let args = e.args();
    let listed = || {
        args.first()
            .map(|l| unescaped(l.args().iter().filter_map(Expr::as_symbol)))
            .unwrap_or_default()
    };
    let bound = match head.as_symbol() {
        Some(Symbol::BLOCK) => Some((ScopeKind::Block, listed())),
        Some(Symbol::FUNCTION) => Some((ScopeKind::Lambda, listed())),
        Some(Symbol::SET_DELAYED | Symbol::UP_SET_DELAYED) => {
            let mut variables = Vec::new();
            if let Some(lhs) = args.first() {
                pattern_variables(lhs, &mut variables);
            }
            Some((ScopeKind::Rule, unescaped(variables)))
        }
        _ => None,
    };
    if let Some((kind, names)) = bound {
        let r = spans.find(path);
        out.push((format!("{kind:?}"), r.start, r.end, names));
    }
    if head.as_symbol().is_none() {
        path.push(0);
        lowered_scopes(head, spans, path, out);
        path.pop();
    }
    for (i, a) in args.iter().enumerate() {
        path.push(i + 1);
        lowered_scopes(a, spans, path, out);
        path.pop();
    }
}

/// The index's scopes but the document: a block's `Local`s, a rule's or a
/// lambda's `Parameter`s.
fn index_scopes(s: &Symbols) -> Vec<Bound> {
    s.scopes()
        .iter()
        .enumerate()
        .filter(|(_, c)| c.kind != ScopeKind::Document)
        .map(|(id, c)| {
            let kind = if c.kind == ScopeKind::Block {
                SymbolKind::Local
            } else {
                SymbolKind::Parameter
            };
            let names = s
                .symbols()
                .iter()
                .filter(|i| i.scope == ScopeId(id) && i.kind == kind)
                .map(|i| i.name.clone())
                .collect();
            (format!("{:?}", c.kind), c.range.start, c.range.end, names)
        })
        .collect()
}

/// The scopes of `text` as the index and lowering have them, sorted, when
/// it parses and lowers without a diagnostic.
fn agreement(text: &str) -> Option<(Vec<Bound>, Vec<Bound>)> {
    let parse = parse(text);
    let names = Names::english();
    let lowered = lower(&parse, names);
    if !parse.errors().is_empty() || !lowered.diagnostics.is_empty() {
        return None;
    }
    let mut ours = index_scopes(&symbols(&parse, names));
    let mut theirs = Vec::new();
    for s in &lowered.statements {
        lowered_scopes(&s.expr, &s.spans, &mut Vec::new(), &mut theirs);
    }
    ours.sort();
    theirs.sort();
    Some((ours, theirs))
}

/// Every binding form, each program without a diagnostic (C0098: capital
/// names lower to `$A`).
#[test]
fn agrees_with_lowering() {
    let programs = [
        "x := 1\n[a, b] := [1, 2]\nv[k] := 0\nFill(7, A)",
        "f(x) := x^2 when x > 0",
        "f(n) := n*f(n - 1) when n > 0",
        "f(x, n = 1) := x^n",
        "f(x, xs...) := g(xs...) + h(xs)",
        "f(n :: integer) := n",
        "f(a)(b) := a b",
        "f(g(x)) := x",
        "twice(g, x) := g(g(x))",
        "phi^2 := phi + 1\nsin(myangle) := 1/2\nf(a) + f(b) := f(a + b)",
        "g := (x, y) -> x y\nh := t -> t^2\nadd := x -> y -> x + y\n(x -> x^2)(3)",
        "begin x := 1; x end",
        "begin end",
        "f(n) := begin m := n; while m > 1 do m := m - 1 end; m end",
        "begin s := 0; for k in 1..3 do s := s + k end; s end",
        "begin if a then x := 1 elif b then x := 2 else x := 3 end; x end",
        "begin for k in 1..3 do t := k end; 0 end",
        "begin local u, v := 2; u := v; [a, b] := [u, v]; a[0] := 1; a end",
        "begin a := 1; b := begin c := a; c end; t -> t + a end",
        "f(x) := if x > 0 then return 1 else return -1 end",
        "f(n) := if n > 0 then r := 1; r else 0 end",
        "begin g := (y -> y^2); y := 3; g(y) end",
        "f(x) := begin g := (x -> x^2); g(x) end",
        "A := begin B := 2; Total := B + 1; Total end\nF(X) := begin Acc := X; Acc end",
        "outer := begin p := 1; q := begin pp := p; begin ppp := pp; ppp end end; q end",
    ];
    let mut wrong = Vec::new();
    for text in programs {
        match agreement(text) {
            None => {
                let parse = parse(text);
                let lowered = lower(&parse, Names::english());
                wrong.push(format!(
                    "{text:?}: {:?} {:?}",
                    parse.errors(),
                    lowered.diagnostics
                ));
            }
            Some((ours, theirs)) if ours != theirs => {
                wrong.push(format!("{text:?}: index {ours:?}, lowering {theirs:?}"));
            }
            Some(_) => {}
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Top-level statements: finished ones of every form, and half-typed ones.
const STATEMENTS: &[&str] = &[
    "x := 1",
    "A := x + 1",
    "[a, b] := [1, 2]",
    "v[k] := 0",
    "f(x) := x^2 when x > 0",
    "f(x, n = 1) := x^n",
    "f(x, xs...) := g(xs...) + h(xs)",
    "f(n :: integer) := n",
    "f(a)(b) := a b",
    "f(g(x)) := x",
    "phi^2 := phi + 1",
    "sin(myangle) := 1/2",
    "f(a) + f(b) := f(a + b)",
    "g := (x, y) -> x y",
    "h := t -> t + A",
    "r := begin s := 0; for k in 1..3 do s := s + k end; s end",
    "begin local u, v := 2; u := v; [a, b] := [u, v]; a[0] := 1; a end",
    "f(n) := begin m := n; while m > 1 do m := m - 1 end; m end",
    "f(x) := if x > 0 then return 1 else return -1 end",
    "begin a := 1; b := begin c := a; c end; t -> t + a end",
    "f(x) := begin g := (x -> x^2); g(x) end",
    "sin(x) + Digits*pi",
    "diff(x^2, x)",
    "Fill(7, A)",
    "f(x, ",
    "y := ",
    "begin z := 1",
    "g(a) := a +",
    "k -> ",
];

/// Short strings an edit inserts, multi-byte ones included.
const PIECES: &[&str] = &[
    "(", ")", ",", ":=", "begin ", " end", "x", "A", "k", "π", "é", "²", " ", "->", "...",
    " :: real", "1", "; ",
];

/// What an occurrence stands for, as §4 compares it.
#[derive(Debug, PartialEq, Eq)]
enum Meaning {
    /// A reserved target or a type: its kind and spelling.
    Reserved(SymbolKind, String),
    /// The global of a key.
    Global(String),
    /// A parameter or a local: its kind, key, and its definitions' name and
    /// construct bytes from the statement's start, `None` for one outside
    /// the statement.
    Bound(SymbolKind, String, Vec<Option<[usize; 4]>>),
}

/// The occurrences inside `statement`: range from its start, binding flag,
/// meaning.
fn occurrences_in(s: &Symbols, statement: &Range<usize>) -> Vec<(Range<usize>, bool, Meaning)> {
    let start = statement.start;
    let inside = |r: &Range<usize>| statement.start <= r.start && r.end <= statement.end;
    s.occurrences()
        .iter()
        .filter(|o| inside(&o.range))
        .map(|o| {
            let info = &s.symbols()[o.symbol.0];
            let meaning = if info.kind.is_user() {
                if info.scope == ScopeId(0) {
                    Meaning::Global(info.name.clone())
                } else {
                    let definitions = info
                        .definitions
                        .iter()
                        .map(|d| {
                            if inside(&d.name) && inside(&d.node) {
                                Some(
                                    [d.name.start, d.name.end, d.node.start, d.node.end]
                                        .map(|b| b - start),
                                )
                            } else {
                                None
                            }
                        })
                        .collect();
                    Meaning::Bound(info.kind, info.name.clone(), definitions)
                }
            } else {
                Meaning::Reserved(info.kind, info.name.clone())
            };
            (
                o.range.start - start..o.range.end - start,
                o.definition,
                meaning,
            )
        })
        .collect()
}

/// The top-level statements' bytes.
fn statements(parse: &Parse) -> Vec<Range<usize>> {
    parse.syntax().children().map(|n| range(&n)).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// An edit inside one statement S that keeps the statement partition
    /// leaves every other statement's occurrences as they were, shifted
    /// (CAS-DOC-60 §4); both texts, when they lower without a diagnostic,
    /// agree with lowering.
    #[test]
    fn stable_under_an_edit_elsewhere(
        lines in prop::collection::vec(prop::sample::select(STATEMENTS), 2..6),
        which in any::<prop::sample::Index>(),
        at in any::<prop::sample::Index>(),
        op in 0..3u8,
        cut in 1..4usize,
        piece in prop::sample::select(PIECES),
    ) {
        let old = lines.join("\n");
        let before = parse(&old);
        let parts = statements(&before);
        prop_assume!(!parts.is_empty());
        let s = which.index(parts.len());
        let r = parts[s].clone();
        // Char boundaries of S: insert (0), delete (1) or replace (2).
        let mut bounds: Vec<usize> = old[r.clone()].char_indices().map(|(i, _)| r.start + i).collect();
        bounds.push(r.end);
        let k = at.index(bounds.len());
        let pos = bounds[k];
        let end = if op == 0 { pos } else { bounds[(k + cut).min(bounds.len() - 1)] };
        let insert = if op == 1 { "" } else { piece };
        let new = format!("{}{insert}{}", &old[..pos], &old[end..]);
        prop_assume!(new != old);
        let after = parse(&new);
        let edited = statements(&after);
        // The partition: as many statements, the others' texts unchanged
        // and shifted when after S.
        prop_assume!(edited.len() == parts.len());
        for (t, (a, b)) in parts.iter().zip(&edited).enumerate() {
            if t != s {
                prop_assume!(old[a.clone()] == new[b.clone()]);
                if t < s {
                    prop_assume!(a == b);
                } else {
                    prop_assume!(a.start + new.len() == b.start + old.len());
                }
            }
        }
        let names = Names::english();
        let (x, y) = (symbols(&before, names), symbols(&after, names));
        for (t, (a, b)) in parts.iter().zip(&edited).enumerate() {
            if t == s {
                continue;
            }
            let was = occurrences_in(&x, a);
            let own = was
                .iter()
                .all(|(_, _, m)| !matches!(m, Meaning::Bound(_, _, d) if d.contains(&None)));
            prop_assert!(own, "statement {} of {:?}: {:?}", t, old, was);
            prop_assert_eq!(was, occurrences_in(&y, b), "statement {} of {:?} edited to {:?}", t, old, new);
        }
        for text in [&old, &new] {
            if let Some((ours, theirs)) = agreement(text) {
                prop_assert_eq!(ours, theirs, "{:?}", text);
            }
        }
    }
}
