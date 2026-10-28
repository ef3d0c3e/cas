use proptest::prelude::*;

use super::{LowerErrorKind as E, Lowered, Severity, lower};
use crate::names::Names;
use crate::parser::{MAX_DEPTH, ParseErrorKind, SyntaxErrorKind, parse};

fn lowered(text: &str) -> Lowered {
    lower(&parse(text), Names::english())
}

/// The full form of each statement of `text`, which must lower without
/// errors.
fn forms(text: &str) -> Vec<String> {
    let lowered = lowered(text);
    assert!(!lowered.has_errors(), "{text:?}: {:?}", lowered.diagnostics);
    lowered
        .statements
        .iter()
        .map(|s| format!("{:?}", s.expr))
        .collect()
}

/// `(text, full form)`, one statement each.
fn check(cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (text, expected) in cases {
        let found = forms(text);
        if found != [*expected] {
            wrong.push(format!("{text:?}: found {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The errors of `text`, as `(covered text, kind)`.
fn errors(text: &str) -> Vec<(&str, E)> {
    lowered(text)
        .diagnostics
        .into_iter()
        .filter(|d| d.severity == Severity::Error)
        .map(|d| (&text[d.range], d.kind))
        .collect()
}

#[test]
fn literals_and_names() {
    check(&[
        ("42", "42"),
        // Decimals, balls when their statement runs (CAS-148).
        ("0.1", "Decimal[1, -1]"),
        ("1.5e-3", "Decimal[15, -4]"),
        ("-2.50", "Decimal[-250, -2]"),
        ("2E10", "Decimal[2, 10]"),
        (r#""s\n""#, r#""s\n""#),
        ("x", "x"),
        ("θ", "θ"),
        // A capitalised name of the user's is escaped (§5).
        ("E", "$E"),
        ("Delta", "$Delta"),
        ("pi", "Pi"),
        ("π", "Pi"),
        ("e", "E"),
        ("i", "Complex[0, 1]"),
        ("true", "True"),
        ("undef", "Indeterminate"),
        ("infinity", "ComplexInfinity"),
        ("inf", "Infinity"),
        ("∞", "Infinity"),
        ("+infinity", "Infinity"),
        ("-infinity", "DirectedInfinity[-1]"),
        ("sin", "Sin"),
    ]);
}

#[test]
fn arithmetic() {
    check(&[
        ("x + 2y", "Plus[x, Times[2, y]]"),
        ("a + b - c + d", "Plus[a, b, Times[-1, c], d]"),
        ("x - 2", "Plus[x, -2]"),
        ("-x", "Times[-1, x]"),
        ("-2", "-2"),
        ("-(x + y)", "Times[-1, Plus[x, y]]"),
        ("a * b / c", "Times[a, b, Power[c, -1]]"),
        ("1/2x", "Times[1, Power[Times[2, x], -1]]"),
        ("2x y", "Times[2, x, y]"),
        ("x^2", "Power[x, 2]"),
        ("x²", "Power[x, 2]"),
        ("x⁻¹", "Power[x, -1]"),
        ("2^3^2", "Power[2, Power[3, 2]]"),
        ("√x", "Power[x, 1/2]"),
        ("∛8", "Power[8, 1/3]"),
        ("a .* b .* c", "ElementwiseTimes[a, b, c]"),
        ("a ./ b", "ElementwiseDivide[a, b]"),
        // Runs of one binding power, folded from the left (CAS-164).
        ("a - b - c", "Plus[a, Times[-1, b], Times[-1, c]]"),
        ("a / b * c", "Times[a, Power[b, -1], c]"),
        ("a * b .* c", "ElementwiseTimes[Times[a, b], c]"),
        ("a .* b * c", "Times[ElementwiseTimes[a, b], c]"),
        ("2x * y .* z", "ElementwiseTimes[Times[2, x, y], z]"),
        (
            "a ./ b ./ c",
            "ElementwiseDivide[ElementwiseDivide[a, b], c]",
        ),
        (
            "a * b .* c / d",
            "Times[ElementwiseTimes[Times[a, b], c], Power[d, -1]]",
        ),
        (
            "a .* b ./ c .* d",
            "ElementwiseTimes[ElementwiseDivide[ElementwiseTimes[a, b], c], d]",
        ),
        ("a and b and c or d", "Or[And[a, b, c], d]"),
        ("a .^ 2", "ElementwisePower[a, 2]"),
        ("n!", "Factorial[n]"),
        ("3x^2(x - 1)", "Times[3, Power[x, 2], Plus[x, -1]]"),
    ]);
}

#[test]
fn comparisons_and_logic() {
    check(&[
        ("a = b", "Equation[a, b]"),
        ("a == b", "Equal[a, b]"),
        ("a != b", "Unequal[a, b]"),
        ("0 < x <= 1", "And[Less[0, x], LessEqual[x, 1]]"),
        ("x >= 2", "GreaterEqual[x, 2]"),
        ("x in S", "Element[x, $S]"),
        ("a and b and c", "And[a, b, c]"),
        ("a or b", "Or[a, b]"),
        ("not p", "Not[p]"),
        ("1..n", "Range[1, n]"),
    ]);
}

/// The set operators (CAS-DOC-32 §3): `union` and `minus` from the left,
/// `intersect` tighter, all between the comparisons and the ranges.
#[test]
fn set_operators() {
    check(&[
        ("a union b", "Union[a, b]"),
        ("a union b union c", "Union[a, b, c]"),
        ("a ∪ b", "Union[a, b]"),
        ("a intersect b ∩ c", "Intersection[a, b, c]"),
        ("a minus b", "Complement[a, b]"),
        ("a minus b minus c", "Complement[a, b, c]"),
        ("a ∖ b", "Complement[a, b]"),
        ("a union b minus c", "Complement[Union[a, b], c]"),
        ("a minus b union c", "Union[Complement[a, b], c]"),
        ("a union b intersect c", "Union[a, Intersection[b, c]]"),
        ("x in a union b", "Element[x, Union[a, b]]"),
        ("1..3 union 5..7", "Union[Range[1, 3], Range[5, 7]]"),
        (
            "{1, 2} union {2, 3}",
            "Union[FiniteSet[1, 2], FiniteSet[2, 3]]",
        ),
    ]);
}

#[test]
fn structure() {
    check(&[
        ("[1, 2]", "List[1, 2]"),
        ("{a, b}", "FiniteSet[a, b]"),
        ("v[0]", "At[v, 0]"),
        ("m[0, k]", "At[m, 0, k]"),
        ("v[1..2]", "At[v, Range[1, 2]]"),
        ("x -> x^2", "Function[List[x], Power[x, 2]]"),
        ("(x, y) -> x y", "Function[List[x, y], Times[x, y]]"),
        ("if c then a else b end", "If[c, a, b]"),
        ("if c then a end", "If[c, a]"),
        ("if c then end", "If[c, Null]"),
        (
            "if c then a elif d then b else f end",
            "If[c, a, If[d, b, f]]",
        ),
    ]);
}

/// Calls and products: giac's rule (CAS-DOC-7 §3.1).
#[test]
fn calls_and_products() {
    check(&[
        ("f(x)", "f[x]"),
        ("f(a, b)", "f[a, b]"),
        ("f()", "f[]"),
        ("sin(x)", "Sin[x]"),
        ("y(x + 1)", "y[Plus[x, 1]]"),
        ("x(x + 1)", "Times[x, Plus[x, 1]]"),
        ("n(n + 1)/2", "Times[n, Plus[n, 1], Power[2, -1]]"),
        ("x(y(x))", "Times[x, y[x]]"),
        ("pi(x + 1)", "Times[Pi, Plus[x, 1]]"),
        ("i(x + 1)", "Times[Complex[0, 1], Plus[x, 1]]"),
        ("2(x + 1)", "Times[2, Plus[x, 1]]"),
        ("(a + b)(a - b)", "Times[Plus[a, b], Plus[a, Times[-1, b]]]"),
        ("2x(x + 1)", "Times[2, x, Plus[x, 1]]"),
        ("f(a)(b)", "f[a][b]"),
        ("x(x)(y)", "Times[Times[x, x], y]"),
        ("g(g(u))", "g[g[u]]"),
        ("call(x, x + 1)", "x[Plus[x, 1]]"),
        ("(x -> x^2)(3)", "Function[List[x], Power[x, 2]][3]"),
        ("v[0](x)", "At[v, 0][x]"),
    ]);
    // The product reading of a name in its argument comes with a note.
    let lowered = lowered("x(x + 1)");
    assert_eq!(lowered.diagnostics.len(), 1);
    let note = &lowered.diagnostics[0];
    assert_eq!(note.severity, Severity::Note);
    assert_eq!(note.to_string(), "read as x*(x + 1)");
    assert!(self::lowered("pi(x + 1)").diagnostics.is_empty());
}

#[test]
fn assignments_and_definitions() {
    check(&[
        ("x := 1", "Set[x, 1]"),
        // The working precision, the one builtin an assignment sets (CAS-151).
        ("Digits := 30", "Set[Digits, 30]"),
        ("Digits", "Digits"),
        ("[a, b] := [1, 2]", "Set[List[a, b], List[1, 2]]"),
        ("v[0] := 1", "Set[v, ReplaceAt[v, 0, 1]]"),
        (
            "f(x) := x^2",
            "SetDelayed[f[Pattern[x, Blank[]]], Power[x, 2]]",
        ),
        ("f(0) := 1", "SetDelayed[f[0], 1]"),
        (
            "f(x, x) := x",
            "SetDelayed[f[Pattern[x, Blank[]], Pattern[x, Blank[]]], x]",
        ),
        ("f(pi) := 0", "SetDelayed[f[Pi], 0]"),
        (
            "f(n) := n*f(n - 1) when n > 0",
            "SetDelayed[f[Pattern[n, Blank[]]], Condition[Times[n, f[Plus[n, -1]]], Greater[n, 0]]]",
        ),
        (
            "f(n :: integer) := n",
            "SetDelayed[f[PatternTest[Pattern[n, Blank[]], IntegerQ]], n]",
        ),
        ("f(x :: any) := x", "SetDelayed[f[Pattern[x, Blank[]]], x]"),
        (
            "f(x, xs...) := g(xs...) + h(xs)",
            "SetDelayed[f[Pattern[x, Blank[]], Pattern[xs, BlankNullSequence[]]], \
             Plus[g[xs], h[List[xs]]]]",
        ),
        (
            "f(x, n = 1) := x^n",
            "SetDelayed[f[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]], 1]], Power[x, n]]",
        ),
        ("f(g(x)) := x", "SetDelayed[f[g[Pattern[x, Blank[]]]], x]"),
        (
            "f(a)(b) := a b",
            "SetDelayed[f[Pattern[a, Blank[]]][Pattern[b, Blank[]]], Times[a, b]]",
        ),
        (
            "twice(g, x) := g(g(x))",
            "SetDelayed[twice[Pattern[g, Blank[]], Pattern[x, Blank[]]], g[g[x]]]",
        ),
        // A builtin head makes an up-rule; names are literal there.
        (
            "phi^2 := phi + 1",
            "UpSetDelayed[Power[phi, 2], Plus[phi, 1]]",
        ),
        (
            "sin(myangle) := 1/2",
            "UpSetDelayed[Sin[myangle], Times[1, Power[2, -1]]]",
        ),
        (
            "f(a) + f(b) := f(a + b)",
            "UpSetDelayed[Plus[f[a], f[b]], f[Plus[a, b]]]",
        ),
        ("assume(x > 0)", "Assume[Greater[x, 0]]"),
        ("purge(f, g)", "Clear[f, g]"),
    ]);
}

/// Standalone rules are equations; `::` marks their variables (§6).
#[test]
fn rules() {
    check(&[
        (
            "replace(p, sin(a :: any)^2 = 1 - cos(a)^2)",
            "ReplaceAll[p, RuleDelayed[Power[Sin[Pattern[a, Blank[]]], 2], \
             Plus[1, Times[-1, Power[Cos[a], 2]]]]]",
        ),
        ("subst(p, x = 2)", "ReplaceAll[p, RuleDelayed[x, 2]]"),
        (
            "subst(p, [x = 1, y = 2])",
            "ReplaceAll[p, List[RuleDelayed[x, 1], RuleDelayed[y, 2]]]",
        ),
        ("subst(p, r)", "ReplaceAll[p, r]"),
    ]);
}

#[test]
fn statements() {
    let lowered = lowered("a; b\nc;\n// done");
    let quiet: Vec<_> = lowered.statements.iter().map(|s| s.quiet).collect();
    assert_eq!(quiet, [true, false, true]);
    // A statement with a parse error is left out; the parser reported it.
    let lowered = self::lowered("x := 1\nf(x\ny := 2");
    assert_eq!(lowered.statements.len(), 1);
    assert!(lowered.diagnostics.is_empty());
}

#[test]
fn errors_with_fixes() {
    let text = "2 sin x";
    let result = lowered(text);
    let error = &result.diagnostics[0];
    assert_eq!(error.kind, E::FunctionWithoutParentheses("sin".into()));
    let fix = error.fix.as_ref().unwrap();
    assert_eq!(
        (&text[fix.range.clone()], fix.replacement.as_str()),
        ("sin x", "sin(x)")
    );
    for (text, covered, replacement) in [
        ("sin^2(x)", "sin^2(x)", "sin(x)^2"),
        ("sin²x", "sin²x", "sin(x)^2"),
        ("2 cos^3 (y + 1)", "cos^3 (y + 1)", "cos(y + 1)^3"),
    ] {
        let lowered = self::lowered(text);
        let error = &lowered.diagnostics[0];
        assert!(matches!(error.kind, E::PowerOfFunction(_)), "{text}");
        let fix = error.fix.as_ref().unwrap();
        assert_eq!(&text[fix.range.clone()], covered, "{text}");
        assert_eq!(fix.replacement, replacement, "{text}");
    }
    assert!(self::lowered("sin^2").diagnostics[0].fix.is_none());
    assert!(self::lowered("x sin").diagnostics[0].fix.is_none());
}

#[test]
fn errors_leave_statements_out() {
    let cases: &[(&str, &[(&str, E)])] = &[
        ("pi := 3", &[("pi", E::ReservedName("pi".into()))]),
        ("sin := 3", &[("sin", E::ReservedName("sin".into()))]),
        ("i -> i^2", &[("i", E::ReservedName("i".into()))]),
        ("[x, e] := l", &[("e", E::ReservedName("e".into()))]),
        ("purge(pi)", &[("pi", E::ReservedName("pi".into()))]),
        ("x :: integer", &[("x :: integer", E::TypeOutsidePattern)]),
        (
            "f(n :: natural) := 1",
            &[("natural", E::UnknownType("natural".into()))],
        ),
        ("g(xs...)", &[("xs...", E::MisplacedSpread)]),
        ("x := 1 when x > 0", &[("x > 0", E::GuardOnAssignment)]),
        ("2 := 3", &[("2", E::InvalidTarget)]),
        ("pi(x) := 1", &[("pi(x)", E::InvalidTarget)]),
        (
            "sin(x :: any) := 0",
            &[("sin(x :: any)", E::NoUserName("sin".into()))],
        ),
        ("call", &[("call", E::MisplacedForm("call".into()))]),
        (
            "f(assume(x))",
            &[("assume(x)", E::MisplacedForm("assume".into()))],
        ),
        // Procedures (§8): statements outside blocks, jumps outside loops.
        ("local x", &[("local x", E::OutsideBlock("local"))]),
        (
            "if a then b; c end",
            &[("b; c", E::OutsideBlock("an if with statements in a branch"))],
        ),
        (
            "for k in r do k end",
            &[("for k in r do k end", E::OutsideBlock("a loop"))],
        ),
        ("return 1", &[("return 1", E::OutsideBlock("return"))]),
        ("break", &[("break", E::OutsideLoop("break"))]),
        (
            "begin continue end",
            &[("continue", E::OutsideLoop("continue"))],
        ),
        // A parameter is not assigned, nor a global; no definitions inside.
        (
            "f(n) := begin n := 1 end",
            &[("n", E::AssignedParameter("n".into()))],
        ),
        (
            "g() := begin h(y) := y end",
            &[("h(y) := y", E::DefinitionInBlock)],
        ),
        // Nor the working precision, which changes the context (CAS-151).
        (
            "g() := begin Digits := 20 end",
            &[("Digits := 20", E::DefinitionInBlock)],
        ),
        (
            "p() := begin t := t + 1 end",
            &[("t", E::ReadBeforeAssignment("t".into()))],
        ),
        (
            "p(c) := begin if c then t := 1 end; t end",
            &[("t", E::ReadBeforeAssignment("t".into()))],
        ),
        (
            "p() := begin while c do t := 1 end; t end",
            &[("t", E::ReadBeforeAssignment("t".into()))],
        ),
        (
            "p() := begin local t; t end",
            &[("t", E::ReadBeforeAssignment("t".into()))],
        ),
        (
            "begin y := if c then a := 1; a else 2 end end",
            &[("a := 1; a", E::StatementsInValue)],
        ),
        (
            "1e2000000",
            &[("1e2000000", E::Number(crate::NumberError::ExponentTooLarge))],
        ),
    ];
    let mut wrong = Vec::new();
    for (text, expected) in cases {
        let found = errors(text);
        if found != *expected {
            wrong.push(format!("{text:?}: found {found:?}"));
        }
        assert!(lowered(text).statements.is_empty(), "{text:?} is left out");
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's option `assume = v` of a call is `Equation[Assume, v]`, as
/// `lpsolve`'s `assume=integer` (CAS-418); elsewhere `assume` stays a
/// misplaced form.
#[test]
fn assume_options() {
    check(&[
        (
            "lpsolve(x, [x >= 1], assume = integer)",
            "LPSolve[x, List[GreaterEqual[x, 1]], Equation[Assume, integer]]",
        ),
        ("f(assume = lp_binary)", "f[Equation[Assume, lp_binary]]"),
    ]);
    assert_eq!(
        errors("y := assume"),
        [("assume", E::MisplacedForm("assume".into()))]
    );
}

/// A planned name, reserved before its builtin exists, is not available
/// to call or to read: `names/en.txt` may list none, so the test adds one.
#[test]
fn planned_names_are_not_available() {
    let names = Names::with_planned("Pending: pending");
    for (text, covered) in [
        ("pending(x^2)", "pending(x^2)"),
        ("y := pending", "pending"),
    ] {
        let lowered = lower(&parse(text), &names);
        let errors: Vec<(&str, E)> = lowered
            .diagnostics
            .into_iter()
            .filter(|d| d.severity == Severity::Error)
            .map(|d| (&text[d.range], d.kind))
            .collect();
        assert_eq!(
            errors,
            [(covered, E::NotAvailable("pending".into()))],
            "{text:?}"
        );
        assert!(lowered.statements.is_empty(), "{text:?} is left out");
    }
}

/// Blocks, their statements and their locals (CAS-DOC-7 §8).
#[test]
fn procedures() {
    check(&[
        ("begin x := 1; x end", "Block[List[x], Set[x, 1], x]"),
        ("begin end", "Block[List[]]"),
        (
            "f(n) := begin m := n; while m > 1 do m := m - 1 end; m end",
            "SetDelayed[f[Pattern[n, Blank[]]], Block[List[m], Set[m, n], \
             While[Greater[m, 1], Statements[Set[m, Plus[m, -1]]]], m]]",
        ),
        (
            "begin s := 0; for k in 1..3 do s := s + k end; s end",
            "Block[List[s, k], Set[s, 0], For[k, Range[1, 3], Statements[Set[s, Plus[s, k]]]], s]",
        ),
        // A statement `if` has groups for branches, `elif` nested.
        (
            "begin if a then x := 1 elif b then x := 2 else x := 3 end; x end",
            "Block[List[x], If[a, Statements[Set[x, 1]], If[b, Statements[Set[x, 2]], \
             Statements[Set[x, 3]]]], x]",
        ),
        // A name first assigned in a loop's body is a local too.
        (
            "begin for k in 1..3 do t := k end; 0 end",
            "Block[List[k, t], For[k, Range[1, 3], Statements[Set[t, k]]], 0]",
        ),
        // `local`, with and without a value; destructuring; indexing.
        (
            "begin local u, v := 2; u := v; [a, b] := [u, v]; a[0] := 1; a end",
            "Block[List[u, v, a, b], Set[v, 2], Set[u, v], Set[List[a, b], List[u, v]], \
             Set[a, ReplaceAt[a, 0, 1]], a]",
        ),
        // Jumps.
        (
            "begin for k in r do if k > 2 then break end; continue end; return end",
            "Block[List[k], For[k, r, Statements[If[Greater[k, 2], Statements[Break[]]], \
             Continue[]]], Return[]]",
        ),
        // A nested block and a lambda are scopes of their own.
        (
            "begin a := 1; b := begin c := a; c end; t -> t + a end",
            "Block[List[a, b], Set[a, 1], Set[b, Block[List[c], Set[c, a], c]], \
             Function[List[t], Plus[t, a]]]",
        ),
        // An `if` with statements is a block of its own as a body.
        (
            "f(x) := if x > 0 then return 1 else return -1 end",
            "SetDelayed[f[Pattern[x, Blank[]]], Block[List[], If[Greater[x, 0], \
             Statements[Return[1]], Statements[Return[-1]]]]]",
        ),
        // A lambda's parameter is not the local of the same name, assigned
        // later.
        (
            "begin g := (y -> y^2); y := 3; g(y) end",
            "Block[List[g, y], Set[g, Function[List[y], Power[y, 2]]], Set[y, 3], g[y]]",
        ),
        // A lambda's parameter is bound, not assigned.
        (
            "f(x) := begin g := (x -> x^2); g(x) end",
            "SetDelayed[f[Pattern[x, Blank[]]], Block[List[g], \
             Set[g, Function[List[x], Power[x, 2]]], g[x]]]",
        ),
    ]);
}

#[test]
fn spans_lead_back_to_the_text() {
    let text = "f(x) + 2 y - z";
    let lowered = lowered(text);
    let spans = &lowered.statements[0].spans;
    let part = |path: &[usize]| &text[spans.find(path)];
    assert_eq!(part(&[]), text);
    assert_eq!(part(&[1]), "f(x)");
    assert_eq!(part(&[1, 0]), "f");
    assert_eq!(part(&[1, 1]), "x");
    assert_eq!(part(&[2]), "2 y");
    assert_eq!(part(&[2, 2]), "y");
    // `-z` is `Times[-1, z]`, made for `z`.
    assert_eq!(part(&[3]), "z");
    assert_eq!(part(&[3, 2]), "z");
    // A path past the parts stops at the nearest one.
    assert_eq!(part(&[2, 2, 5]), "y");
}

/// Pieces of the language, for random programs.
fn fragment() -> impl Strategy<Value = &'static str> {
    prop::sample::select(vec![
        "x", "f", "2", "sin", "pi", "i", "E", "(", ")", "[", "]", ",", "+", "-", "*", "/", "^",
        "²", "!", "...", "::", ":=", "->", "=", "==", "<", "..", " ", "\n", ";", "and", "not",
        "if", "then", "else", "end", "when", "integer", "call", "assume", "subst", "√",
    ])
}

proptest! {
    /// Lowering never panics, and every range it reports lies in the text.
    #[test]
    fn lowering_is_total(fragments in prop::collection::vec(fragment(), 0..30)) {
        let text = fragments.concat();
        let lowered = lowered(&text);
        for d in &lowered.diagnostics {
            prop_assert!(d.range.end <= text.len(), "{:?} in {:?}", d, text);
            if let Some(fix) = &d.fix {
                prop_assert!(fix.range.end <= text.len());
            }
        }
        for s in &lowered.statements {
            prop_assert!(s.spans.range().end <= text.len());
        }
    }
}

/// A fuzz finding (CAS-163): thousands of `!` after a name nested as deep as
/// the text was long, past the parser's depth limit, and lowering's recursion
/// overflowed the stack. The text now nests too deeply to parse, and lowers
/// what parsed.
#[test]
fn long_postfix_chains_stay_within_the_depth_limit() {
    for text in [
        "h0h".to_owned() + &"!".repeat(4000),
        "f".to_owned() + &"(1)".repeat(2000),
        "a".to_owned() + &"[1]".repeat(2000),
    ] {
        let parse = parse(&text);
        assert!(
            parse
                .errors()
                .iter()
                .any(|e| e.kind == SyntaxErrorKind::Parse(ParseErrorKind::TooDeep)),
            "{}",
            &text[..20]
        );
        let _ = lower(&parse, Names::english());
    }
}

/// A run of binary operators of one binding power is one node, which
/// lowering takes in a loop (CAS-164): long sums and products lower flat,
/// where they nested once per operator, and lowering recursed as deep as
/// they were long, and dropping the tree too. `./` after `./`, and `*`
/// alternating with `.*`, nest once lowered, so they are levels of nesting:
/// within the depth limit they lower in a loop, past it they do not parse.
#[test]
fn binary_runs_lower_in_a_loop() {
    let lowered_one = |text: &str| {
        let lowered = lowered(text);
        assert!(!lowered.has_errors(), "{:?}", lowered.diagnostics);
        lowered.statements[0].expr.clone()
    };
    for op in [" + ", " * "] {
        let flat = lowered_one(&vec!["x"; 100_000].join(op));
        assert_eq!(flat.args().len(), 100_000, "{op}");
    }
    // Each `./` divides all before it, and each run of `*` or `.*` takes
    // all before it as its first operand.
    let n = MAX_DEPTH / 2;
    for text in [
        vec!["a"; n + 1].join(" ./ "),
        "a".to_owned() + &" * a .* a".repeat(n / 2),
    ] {
        let mut depth = 0;
        let mut inner = lowered_one(&text);
        while let [first, _] = inner.args() {
            inner = first.clone();
            depth += 1;
        }
        assert_eq!(depth, n, "{}", &text[..20]);
    }
    for text in [
        vec!["a"; 4001].join(" ./ "),
        "a".to_owned() + &" * a .* a".repeat(2000),
    ] {
        let parse = parse(&text);
        assert!(
            parse
                .errors()
                .iter()
                .any(|e| e.kind == SyntaxErrorKind::Parse(ParseErrorKind::TooDeep)),
            "{}",
            &text[..20]
        );
        let _ = lower(&parse, Names::english());
    }
}

/// Superscripts, `=` and `..` wrap what came before them too: past the depth
/// limit, the text does not parse, and lowering what did stays within the
/// stack (CAS-168).
#[test]
fn long_chains_of_other_forms_stay_within_the_depth_limit() {
    for text in [
        "x".to_owned() + &"² ".repeat(2000),
        "a".to_owned() + &" = a".repeat(2000),
        "a".to_owned() + &"..a".repeat(2000),
        "a".to_owned() + &" < a = a".repeat(1000),
    ] {
        let parse = parse(&text);
        assert!(
            parse
                .errors()
                .iter()
                .any(|e| e.kind == SyntaxErrorKind::Parse(ParseErrorKind::TooDeep)),
            "{}",
            text.chars().take(20).collect::<String>()
        );
        let _ = lower(&parse, Names::english());
    }
}
