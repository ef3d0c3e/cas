use cas_core::Symbol;

use super::{Occurrence, ScopeKind, SymbolInfo, SymbolKind, Symbols, symbols};
use crate::names::{Names, SyntaxForm, TypeName};
use crate::parser::parse;

fn index(text: &str) -> Symbols {
    symbols(&parse(text), Names::english())
}

/// Each occurrence as `text kind scope range`, `text=` for a binding site,
/// the scope being the symbol's.
fn describe(text: &str) -> Vec<String> {
    let s = index(text);
    s.occurrences()
        .iter()
        .map(|o| {
            let info = &s.symbols()[o.symbol.0];
            let scope = &s.scopes()[info.scope.0];
            let kind = match info.kind {
                SymbolKind::Builtin(_) => "Builtin".to_owned(),
                kind => format!("{kind:?}"),
            };
            let mark = if o.definition { "=" } else { "" };
            let name = &text[o.range.clone()];
            format!("{name}{mark} {kind} {:?} {:?}", scope.kind, scope.range)
        })
        .collect()
}

/// The symbol of the first occurrence spelled `word`.
fn symbol<'s>(s: &'s Symbols, text: &str, word: &str) -> &'s SymbolInfo {
    let o = s
        .occurrences()
        .iter()
        .find(|o| &text[o.range.clone()] == word);
    &s.symbols()[o.expect(word).symbol.0]
}

/// The construct ranges of the definitions of the first occurrence `word`,
/// as `(start, end)`.
fn sites(text: &str, word: &str) -> Vec<(usize, usize)> {
    let s = index(text);
    let info = symbol(&s, text, word);
    info.definitions
        .iter()
        .map(|d| (d.node.start, d.node.end))
        .collect()
}

/// Whether the occurrences at `a` and `b` stand for one symbol.
fn same(text: &str, a: usize, b: usize) -> bool {
    let s = index(text);
    let o = s.occurrences();
    o[a].symbol == o[b].symbol
}

#[test]
fn assignment() {
    assert_eq!(describe("x := 1"), ["x= Variable Document 0..6"]);
    assert_eq!(sites("x := 1", "x"), [(0, 6)]);
    assert_eq!(describe("x := "), ["x= Variable Document 0..5"]);
    // The next line becomes the value (CAS-665's tree).
    assert_eq!(
        describe("x := \ny := 1\nz := 2"),
        [
            "x= Variable Document 0..19",
            "y Variable Document 0..19",
            "z= Variable Document 0..19"
        ]
    );
}

#[test]
fn rule_parameters() {
    assert_eq!(
        describe("f(x) := x"),
        [
            "f= Function Document 0..9",
            "x= Parameter Rule 0..9",
            "x Parameter Rule 0..9"
        ]
    );
    assert_eq!(sites("f(x) := x", "f"), [(0, 9)]);
    assert_eq!(sites("f(x) := x", "x"), [(2, 3)]);
    assert_eq!(
        describe("f(x) := "),
        ["f= Function Document 0..8", "x= Parameter Rule 0..7"]
    );
    // Half-typed, a call: its names are uses.
    assert_eq!(
        describe("f(x, "),
        ["f Variable Document 0..5", "x Variable Document 0..5"]
    );
}

#[test]
fn parameter_forms() {
    let text = "g(a, n = k, xs..., m :: real) := n";
    assert_eq!(
        describe(text),
        [
            "g= Function Document 0..34",
            "a= Parameter Rule 0..34",
            "n= Parameter Rule 0..34",
            "k Variable Document 0..34",
            "xs= Parameter Rule 0..34",
            "m= Parameter Rule 0..34",
            "real Type(Real) Document 0..34",
            "n Parameter Rule 0..34"
        ]
    );
    // The default, the spread and the typed name, each whole.
    assert_eq!(sites(text, "n"), [(5, 10)]);
    assert_eq!(sites(text, "xs"), [(12, 17)]);
    assert_eq!(sites(text, "m"), [(19, 28)]);
}

#[test]
fn sub_rule() {
    assert_eq!(
        describe("f(a)(b) := a + b"),
        [
            "f= Function Document 0..16",
            "a= Parameter Rule 0..16",
            "b= Parameter Rule 0..16",
            "a Parameter Rule 0..16",
            "b Parameter Rule 0..16"
        ]
    );
}

#[test]
fn up_rule() {
    // Tagged names are parameters, plain user names globals the up-rule
    // attaches to, and the builtin head is no definition.
    let text = "sin(x :: real, y) := x + y";
    assert_eq!(
        describe(text),
        [
            "sin Builtin Document 0..26",
            "x= Parameter Rule 0..26",
            "real Type(Real) Document 0..26",
            "y Variable Document 0..26",
            "x Parameter Rule 0..26",
            "y Variable Document 0..26"
        ]
    );
    assert_eq!(sites(text, "x"), [(4, 13)]);
    assert_eq!(sites(text, "sin"), []);
}

#[test]
fn list_and_element_assignment() {
    assert_eq!(
        describe("[a, b] := [1, 2]"),
        ["a= Variable Document 0..16", "b= Variable Document 0..16"]
    );
    assert_eq!(sites("[a, b] := [1, 2]", "b"), [(0, 16)]);
    // `v[k] := e` modifies `v`.
    assert_eq!(
        describe("v[k] := 0"),
        ["v= Variable Document 0..9", "k Variable Document 0..9"]
    );
    assert_eq!(describe("[1, "), Vec::<String>::new());
}

#[test]
fn lambda_parameters() {
    assert_eq!(
        describe("g := (x, y) -> x*y + z"),
        [
            "g= Variable Document 0..22",
            "x= Parameter Lambda 5..22",
            "y= Parameter Lambda 5..22",
            "x Parameter Lambda 5..22",
            "y Parameter Lambda 5..22",
            "z Variable Document 0..22"
        ]
    );
    assert_eq!(sites("h := t -> t", "t"), [(5, 6)]);
}

#[test]
fn block_locals() {
    // By assignment, `local` and the loop variable.
    let text = "r := begin x := 1; local y; for k in v do y := k end end";
    assert_eq!(
        describe(text),
        [
            "r= Variable Document 0..56",
            "x= Local Block 5..56",
            "y= Local Block 5..56",
            "k= Local Block 5..56",
            "v Variable Document 0..56",
            "y= Local Block 5..56",
            "k Local Block 5..56"
        ]
    );
    assert_eq!(sites(text, "y"), [(25, 26), (42, 48)]);
    assert_eq!(sites(text, "k"), [(28, 52)]);
    assert_eq!(describe("begin x := 1"), ["x= Local Block 0..12"]);
}

#[test]
fn locals_leave_out_parameters() {
    // A user rule's parameter stays the parameter, as lowering has it.
    assert_eq!(
        describe("f(x) := begin x := 1; y := x end"),
        [
            "f= Function Document 0..32",
            "x= Parameter Rule 0..32",
            "x= Parameter Rule 0..32",
            "y= Local Block 8..32",
            "x Parameter Rule 0..32"
        ]
    );
    assert_eq!(
        describe("g := x -> begin x := 1 end"),
        [
            "g= Variable Document 0..26",
            "x= Parameter Lambda 5..26",
            "x= Parameter Lambda 5..26"
        ]
    );
    // An up-rule leaves lowering's parameters alone: `x` is a local.
    assert_eq!(
        describe("sin(x :: real) := begin x := 1; x end"),
        [
            "sin Builtin Document 0..37",
            "x= Parameter Rule 0..37",
            "real Type(Real) Document 0..37",
            "x= Local Block 18..37",
            "x Local Block 18..37"
        ]
    );
}

#[test]
fn nested_scopes() {
    let text = "begin x := 1; begin x := 2 end; y := x -> x end";
    assert_eq!(
        describe(text),
        [
            "x= Local Block 0..47",
            "x= Local Block 14..30",
            "y= Local Block 0..47",
            "x= Parameter Lambda 37..43",
            "x Parameter Lambda 37..43"
        ]
    );
    let s = index(text);
    let scopes: Vec<_> = s
        .scopes()
        .iter()
        .map(|c| (c.kind, c.range.clone(), c.parent.map(|p| p.0)))
        .collect();
    assert_eq!(
        scopes,
        [
            (ScopeKind::Document, 0..47, None),
            (ScopeKind::Block, 0..47, Some(0)),
            (ScopeKind::Block, 14..30, Some(1)),
            (ScopeKind::Lambda, 37..43, Some(1))
        ]
    );
}

#[test]
fn if_body_is_a_block() {
    assert_eq!(
        describe("f(n) := if n > 0 then r := 1; r else 0 end"),
        [
            "f= Function Document 0..42",
            "n= Parameter Rule 0..42",
            "n Parameter Rule 0..42",
            "r= Local Block 8..42",
            "r Local Block 8..42"
        ]
    );
}

#[test]
fn fill() {
    let text = "Fill(v, A)";
    assert_eq!(
        describe(text),
        [
            "Fill Syntax(Fill) Document 0..10",
            "v Variable Document 0..10",
            "A= Variable Document 0..10"
        ]
    );
    assert_eq!(sites(text, "A"), [(0, 10)]);
    // The source spelling, not lowering's escaped `$A`.
    assert_eq!(symbol(&index(text), text, "A").name, "A");
}

#[test]
fn reserved_in_binding_position_is_a_use() {
    assert_eq!(describe("sin := 3"), ["sin Builtin Document 0..8"]);
    assert_eq!(sites("sin := 3", "sin"), []);
}

#[test]
fn reserved_targets() {
    let text = "π + pi";
    let s = index(text);
    assert!(same(text, 0, 1));
    assert_eq!(symbol(&s, text, "π").kind, SymbolKind::Builtin(Symbol::PI));
    assert_eq!(symbol(&s, text, "π").name, "pi");

    let text = "ln(x) + log(x)";
    let s = index(text);
    assert!(same(text, 0, 2));
    assert_eq!(
        symbol(&s, text, "log").kind,
        SymbolKind::Builtin(Symbol::LOG)
    );
    assert_eq!(symbol(&s, text, "log").name, "ln");

    let text = "Digits := 20";
    let s = index(text);
    assert_eq!(
        symbol(&s, text, "Digits").kind,
        SymbolKind::Builtin(Symbol::DIGITS)
    );
    assert_eq!(describe(text), ["Digits Builtin Document 0..12"]);

    let text = "i^2";
    assert_eq!(symbol(&index(text), text, "i").kind, SymbolKind::Imaginary);
    // `Planned`: en.txt's [planned] section is empty, so there is no name
    // to test.

    let text = "call(f, x); assume(x > 0); purge(x); Fill(v, A)";
    let s = index(text);
    for (word, form) in [
        ("call", SyntaxForm::Call),
        ("assume", SyntaxForm::Assume),
        ("purge", SyntaxForm::Purge),
        ("Fill", SyntaxForm::Fill),
    ] {
        assert_eq!(
            symbol(&s, text, word).kind,
            SymbolKind::Syntax(form),
            "{word}"
        );
    }
}

#[test]
fn types() {
    let text = "f(x :: real, y :: foo) := x";
    assert_eq!(
        describe(text),
        [
            "f= Function Document 0..27",
            "x= Parameter Rule 0..27",
            "real Type(Real) Document 0..27",
            "y= Parameter Rule 0..27",
            "x Parameter Rule 0..27"
        ]
    );
    assert_eq!(
        symbol(&index(text), text, "real").kind,
        SymbolKind::Type(TypeName::Real)
    );
}

#[test]
fn globals() {
    // No binding site: a variable without definitions.
    let text = "diff(x^2, x)";
    let s = index(text);
    let x = symbol(&s, text, "x");
    assert_eq!(
        (x.kind, x.definitions.len(), x.uses.len()),
        (SymbolKind::Variable, 0, 2)
    );
    // One global per key, whatever the order; a rule definition makes it a
    // function.
    let text = "y := f(2); f := 1; f(x) := x";
    let s = index(text);
    let f = symbol(&s, text, "f");
    assert_eq!(
        (f.kind, f.definitions.len(), f.uses.len()),
        (SymbolKind::Function, 2, 1)
    );
    assert!(same(text, 1, 2) && same(text, 1, 3));
}

#[test]
fn local_read_before_assignment() {
    let text = "begin y := x; x := 1 end";
    assert_eq!(
        describe(text),
        [
            "y= Local Block 0..24",
            "x Local Block 0..24",
            "x= Local Block 0..24"
        ]
    );
    assert!(same(text, 1, 2));
}

#[test]
fn syntax_errors() {
    // The statement keeps what its tree shows; the names in the error node
    // (`:= 1 z := 2`) are none.
    assert_eq!(
        describe("f(g(x, \ny := 1\nz := 2"),
        [
            "f Variable Document 0..21",
            "g Variable Document 0..21",
            "x Variable Document 0..21",
            "y Variable Document 0..21"
        ]
    );
}

#[test]
fn at_boundaries() {
    let s = index("sin(x)");
    let at = |o: usize| s.at(o).map(|o: &Occurrence| o.range.clone());
    assert_eq!(at(3), Some(0..3)); // sin|(
    assert_eq!(at(4), Some(4..5)); // (|x
    assert_eq!(at(1), Some(0..3));
    assert_eq!(at(6), None); // after `)`, the end
    assert_eq!(at(usize::MAX), None);
    assert_eq!(index("x").at(1).map(|o| o.range.clone()), Some(0..1));
    assert_eq!(index("").at(0), None);
    assert_eq!(index("a + b").at(2), None);
}
