use std::cmp::Ordering;
use std::fmt::Write as _;

use proptest::prelude::*;
use rowan::NodeOrToken;

use super::{MAX_DEPTH, ParseErrorKind as E, SyntaxErrorKind, parse};
use crate::syntax::{SyntaxKind as K, SyntaxNode};

/// The tree without trivia: `(Kind child …)` for a node, the text of a
/// token, and only the text for a name or a literal.
fn sexp(node: &SyntaxNode) -> String {
    let mut out = String::new();
    write_sexp(node, &mut out);
    out
}

fn write_sexp(node: &SyntaxNode, out: &mut String) {
    if matches!(node.kind(), K::Name | K::Literal) {
        write!(out, "{}", node.text()).unwrap();
        return;
    }
    write!(out, "({:?}", node.kind()).unwrap();
    for child in node.children_with_tokens() {
        match child {
            NodeOrToken::Node(n) => {
                out.push(' ');
                write_sexp(&n, out);
            }
            NodeOrToken::Token(t) if !t.kind().is_trivia() && t.kind() != K::Newline => {
                write!(out, " {}", t.text()).unwrap();
            }
            NodeOrToken::Token(_) => {}
        }
    }
    out.push(')');
}

/// The parse errors of `text`, as `(covered text, kind)`.
fn parse_errors(text: &str) -> Vec<(&str, E)> {
    parse(text)
        .errors()
        .iter()
        .filter_map(|e| match e.kind {
            SyntaxErrorKind::Parse(kind) => Some((&text[e.range.clone()], kind)),
            SyntaxErrorKind::Lex(_) => None,
        })
        .collect()
}

/// The tree of `text`, which must have no errors and keep every byte.
fn tree(text: &str) -> String {
    let parse = parse(text);
    assert!(parse.errors().is_empty(), "{text:?}: {:?}", parse.errors());
    assert_eq!(parse.syntax().text().to_string(), text);
    sexp(&parse.syntax())
}

/// The tree of the one expression statement of `text`.
fn expr(text: &str) -> String {
    let tree = tree(text);
    tree.strip_prefix("(Root (ExprStmt ")
        .and_then(|t| t.strip_suffix("))"))
        .unwrap_or_else(|| panic!("{text:?} is not one expression: {tree}"))
        .to_owned()
}

fn check_exprs(cases: &[(&str, &str)]) {
    for (text, expected) in cases {
        assert_eq!(expr(text), *expected, "{text:?}");
    }
}

fn check_trees(cases: &[(&str, &str)]) {
    for (text, expected) in cases {
        assert_eq!(tree(text), *expected, "{text:?}");
    }
}

#[test]
fn primaries() {
    check_exprs(&[
        ("x", "x"),
        ("42", "42"),
        ("0.5", "0.5"),
        (r#""s""#, r#""s""#),
        ("(x)", "(Paren ( x ))"),
        ("[]", "(ListExpr [ ])"),
        ("[1, 2, 3]", "(ListExpr [ 1 , 2 , 3 ])"),
        ("[1, 2,]", "(ListExpr [ 1 , 2 , ])"),
        ("{a, b}", "(SetExpr { a , b })"),
        ("π", "π"),
        ("∞", "∞"),
    ]);
}

#[test]
fn postfix_forms() {
    check_exprs(&[
        ("f(x)", "(Application f (ArgList ( x )))"),
        ("f()", "(Application f (ArgList ( )))"),
        ("f(a, b)", "(Application f (ArgList ( a , b )))"),
        ("f (x)", "(Application f (ArgList ( x )))"),
        (
            "f(a)(b)",
            "(Application (Application f (ArgList ( a ))) (ArgList ( b )))",
        ),
        ("v[0]", "(IndexExpr v (ArgList [ 0 ]))"),
        ("m[i, j]", "(IndexExpr m (ArgList [ i , j ]))"),
        (
            "v[0](x)",
            "(Application (IndexExpr v (ArgList [ 0 ])) (ArgList ( x )))",
        ),
        ("n!", "(Factorial n !)"),
        ("f(x)!", "(Factorial (Application f (ArgList ( x ))) !)"),
        ("xs...", "(Spread xs ...)"),
        ("n :: integer", "(Typed n :: integer)"),
        (
            "(x -> x^2)(3)",
            "(Application (Paren ( (Lambda x -> (PowerExpr x ^ 2)) )) (ArgList ( 3 )))",
        ),
    ]);
}

/// A number or a parenthesised expression before `(…)` makes a product, not
/// an application, so it binds below `^` (CAS-DOC-7 §3.1, L4).
#[test]
fn products_with_parentheses() {
    check_exprs(&[
        (
            "2(x + 1)",
            "(ImplicitProduct 2 (Paren ( (BinaryExpr x + 1) )))",
        ),
        (
            "(a + b)(a - b)",
            "(ImplicitProduct (Paren ( (BinaryExpr a + b) )) (Paren ( (BinaryExpr a - b) )))",
        ),
        (
            "2(x + 1)(x - 1)",
            "(ImplicitProduct 2 (Paren ( (BinaryExpr x + 1) )) (Paren ( (BinaryExpr x - 1) )))",
        ),
        (
            "x^2(x + 1)",
            "(ImplicitProduct (PowerExpr x ^ 2) (Paren ( (BinaryExpr x + 1) )))",
        ),
        (
            "3x^2(x - 1)",
            "(ImplicitProduct 3 (PowerExpr x ^ 2) (Paren ( (BinaryExpr x - 1) )))",
        ),
        (
            "(a + b)x",
            "(ImplicitProduct (Paren ( (BinaryExpr a + b) )) x)",
        ),
        // After a name it is an application; lowering decides (§3.1).
        (
            "x(x + 1)",
            "(Application x (ArgList ( (BinaryExpr x + 1) )))",
        ),
        (
            "2x(x + 1)",
            "(ImplicitProduct 2 (Application x (ArgList ( (BinaryExpr x + 1) ))))",
        ),
        ("e^f(x)", "(PowerExpr e ^ (Application f (ArgList ( x ))))"),
    ]);
}

#[test]
fn implicit_products_and_powers() {
    check_exprs(&[
        ("2x", "(ImplicitProduct 2 x)"),
        ("x y z", "(ImplicitProduct x y z)"),
        ("1/2x", "(BinaryExpr 1 / (ImplicitProduct 2 x))"),
        ("a/2b", "(BinaryExpr a / (ImplicitProduct 2 b))"),
        ("x^2y", "(ImplicitProduct (PowerExpr x ^ 2) y)"),
        ("2^3x", "(ImplicitProduct (PowerExpr 2 ^ 3) x)"),
        ("2x^2", "(ImplicitProduct 2 (PowerExpr x ^ 2))"),
        ("2^3^2", "(PowerExpr 2 ^ (PowerExpr 3 ^ 2))"),
        ("x^-2", "(PowerExpr x ^ (PrefixExpr - 2))"),
        ("e^-x", "(PowerExpr e ^ (PrefixExpr - x))"),
        (
            "x^-2y",
            "(ImplicitProduct (PowerExpr x ^ (PrefixExpr - 2)) y)",
        ),
        ("a .^ 2", "(PowerExpr a .^ 2)"),
        ("x²", "(PowerExpr x ²)"),
        ("x²y", "(ImplicitProduct (PowerExpr x ²) y)"),
        ("(a+b)¹⁰", "(PowerExpr (Paren ( (BinaryExpr a + b) )) ¹⁰)"),
        ("x^2!", "(PowerExpr x ^ (Factorial 2 !))"),
        ("x!^2", "(PowerExpr (Factorial x !) ^ 2)"),
        (
            "sin(x)cos(x)",
            "(ImplicitProduct (Application sin (ArgList ( x ))) (Application cos (ArgList ( x ))))",
        ),
        ("2πr", "(ImplicitProduct 2 π r)"),
    ]);
}

#[test]
fn prefix_operators() {
    check_exprs(&[
        ("-x", "(PrefixExpr - x)"),
        ("-x^2", "(PrefixExpr - (PowerExpr x ^ 2))"),
        ("-2x", "(ImplicitProduct (PrefixExpr - 2) x)"),
        ("--x", "(PrefixExpr - (PrefixExpr - x))"),
        ("−x", "(PrefixExpr − x)"),
        ("-x!", "(PrefixExpr - (Factorial x !))"),
        ("a -b", "(BinaryExpr a - b)"),
        ("a - -b", "(BinaryExpr a - (PrefixExpr - b))"),
        ("√x", "(PrefixExpr √ x)"),
        ("√x^2", "(PrefixExpr √ (PowerExpr x ^ 2))"),
        ("√2x", "(ImplicitProduct (PrefixExpr √ 2) x)"),
        ("2√x", "(ImplicitProduct 2 (PrefixExpr √ x))"),
        ("∛8", "(PrefixExpr ∛ 8)"),
        ("not p", "(PrefixExpr not p)"),
        ("not p and q", "(BinaryExpr (PrefixExpr not p) and q)"),
        ("not a < b", "(PrefixExpr not (Comparison a < b))"),
        ("not not p", "(PrefixExpr not (PrefixExpr not p))"),
    ]);
}

#[test]
fn comparisons_equations_and_ranges() {
    check_exprs(&[
        ("a = b", "(Equation a = b)"),
        (
            "x^2 + 1 = 0",
            "(Equation (BinaryExpr (PowerExpr x ^ 2) + 1) = 0)",
        ),
        ("a == b", "(Comparison a == b)"),
        ("0 < x <= 1", "(Comparison 0 < x <= 1)"),
        ("a != b", "(Comparison a != b)"),
        ("x in S", "(Comparison x in S)"),
        ("x ∈ S", "(Comparison x ∈ S)"),
        ("1 ≤ x", "(Comparison 1 ≤ x)"),
        ("0..n-1", "(RangeExpr 0 .. (BinaryExpr n - 1))"),
        ("x = 0..1", "(Equation x = (RangeExpr 0 .. 1))"),
        ("a and b or c", "(BinaryExpr (BinaryExpr a and b) or c)"),
        ("a or b and c", "(BinaryExpr a or (BinaryExpr b and c))"),
        (
            "x > 0 and x < 1",
            "(BinaryExpr (Comparison x > 0) and (Comparison x < 1))",
        ),
    ]);
}

#[test]
fn lambdas() {
    check_exprs(&[
        ("x -> x^2", "(Lambda x -> (PowerExpr x ^ 2))"),
        (
            "(x, y) -> x + y",
            "(Lambda (ParamList ( x , y )) -> (BinaryExpr x + y))",
        ),
        ("() -> 1", "(Lambda (ParamList ( )) -> 1)"),
        (
            "x -> y -> x y",
            "(Lambda x -> (Lambda y -> (ImplicitProduct x y)))",
        ),
        ("x → 2x", "(Lambda x → (ImplicitProduct 2 x))"),
        (
            "map(x -> x^2, l)",
            "(Application map (ArgList ( (Lambda x -> (PowerExpr x ^ 2)) , l )))",
        ),
        // Not a lambda: `(a b)` is a product.
        ("(a b)", "(Paren ( (ImplicitProduct a b) ))"),
    ]);
}

/// Binary forms from the loosest to the tightest (CAS-DOC-7 §3): the
/// operator, the node, the binding power and how it chains.
const LEVELS: &[(&str, &str, u8, Chain)] = &[
    ("or", "BinaryExpr", 1, Chain::Left),
    ("and", "BinaryExpr", 3, Chain::Left),
    ("<", "Comparison", 7, Chain::Flat),
    ("..", "RangeExpr", 9, Chain::Error),
    ("+", "BinaryExpr", 11, Chain::Left),
    ("*", "BinaryExpr", 13, Chain::Left),
    (" ", "ImplicitProduct", 15, Chain::Flat),
    ("^", "PowerExpr", 20, Chain::Right),
];

#[derive(Clone, Copy, PartialEq)]
enum Chain {
    Left,
    Right,
    /// One node for the whole chain.
    Flat,
    /// Chaining is an error.
    Error,
}

/// The tree of `lhs op rhs` for an operator of `LEVELS`.
fn form(node: &str, op: &str, lhs: &str, rhs: &str) -> String {
    if op == " " {
        format!("({node} {lhs} {rhs})")
    } else {
        format!("({node} {lhs} {op} {rhs})")
    }
}

/// Every pair of levels, in both orders: `a op1 b op2 c` groups by binding
/// power, and within a level by how the level chains.
#[test]
fn every_pair_of_levels() {
    for &(op1, node1, bp1, chain) in LEVELS {
        for &(op2, node2, bp2, _) in LEVELS {
            let join = |op: &str| {
                if op == " " {
                    " ".to_owned()
                } else {
                    format!(" {op} ")
                }
            };
            let text = format!("a{}b{}c", join(op1), join(op2));
            let left = form(node2, op2, &form(node1, op1, "a", "b"), "c");
            let right = form(node1, op1, "a", &form(node2, op2, "b", "c"));
            match (bp1.cmp(&bp2), chain) {
                (Ordering::Greater, _) | (Ordering::Equal, Chain::Left) => {
                    assert_eq!(expr(&text), left, "{text}");
                }
                (Ordering::Less, _) | (Ordering::Equal, Chain::Right) => {
                    assert_eq!(expr(&text), right, "{text}");
                }
                (Ordering::Equal, Chain::Flat) => {
                    let expected = if op1 == " " {
                        format!("({node1} a b c)")
                    } else {
                        format!("({node1} a {op1} b {op2} c)")
                    };
                    assert_eq!(expr(&text), expected, "{text}");
                }
                (Ordering::Equal, Chain::Error) => {
                    assert_eq!(parse_errors(&text), [("..", E::ChainedRange)], "{text}");
                }
            }
        }
    }
}

#[test]
fn statements_and_definitions() {
    check_trees(&[
        ("", "(Root)"),
        ("x := 1", "(Root (Definition x := 1))"),
        (
            "f(x) := x^2",
            "(Root (Definition (Application f (ArgList ( x ))) := (PowerExpr x ^ 2)))",
        ),
        (
            "f(n) := n*f(n-1) when n > 0",
            "(Root (Definition (Application f (ArgList ( n ))) := (BinaryExpr n * \
             (Application f (ArgList ( (BinaryExpr n - 1) )))) (WhenClause when (Comparison n > 0))))",
        ),
        (
            "[a, b] := e",
            "(Root (Definition (ListExpr [ a , b ]) := e))",
        ),
        (
            "v[0] := 1",
            "(Root (Definition (IndexExpr v (ArgList [ 0 ])) := 1))",
        ),
        ("a; b", "(Root (ExprStmt a) ; (ExprStmt b))"),
        ("a;", "(Root (ExprStmt a) ;)"),
        (";;a", "(Root ; ; (ExprStmt a))"),
        ("a\nb", "(Root (ExprStmt a) (ExprStmt b))"),
        ("a\n\n\nb\n", "(Root (ExprStmt a) (ExprStmt b))"),
        ("// note\nx // more", "(Root (ExprStmt x))"),
    ]);
}

/// A line break ends a complete statement, unless the next line starts with
/// a token that cannot start a statement (CAS-DOC-7 §4).
#[test]
fn line_breaks() {
    check_trees(&[
        // Incomplete: continues.
        ("a +\nb", "(Root (ExprStmt (BinaryExpr a + b)))"),
        (
            "f(a,\nb)",
            "(Root (ExprStmt (Application f (ArgList ( a , b )))))",
        ),
        ("[1,\n2]", "(Root (ExprStmt (ListExpr [ 1 , 2 ])))"),
        // Inside brackets, line breaks do not separate.
        (
            "(a\nb)",
            "(Root (ExprStmt (Paren ( (ImplicitProduct a b) ))))",
        ),
        // Tokens that cannot start a statement continue it.
        ("a\n* b", "(Root (ExprStmt (BinaryExpr a * b)))"),
        ("a\nand b", "(Root (ExprStmt (BinaryExpr a and b)))"),
        ("x\n:= 1", "(Root (Definition x := 1))"),
        (
            "f(n) := 1\n  when n == 0",
            "(Root (Definition (Application f (ArgList ( n ))) := 1 (WhenClause when (Comparison n == 0))))",
        ),
        // Tokens that can start one end it: no product, call or sum across
        // lines.
        ("a\nb", "(Root (ExprStmt a) (ExprStmt b))"),
        ("f\n(x)", "(Root (ExprStmt f) (ExprStmt (Paren ( x ))))"),
        ("a\n- b", "(Root (ExprStmt a) (ExprStmt (PrefixExpr - b)))"),
        ("a\n[0]", "(Root (ExprStmt a) (ExprStmt (ListExpr [ 0 ])))"),
        ("return\nx", "(Root (ReturnStmt return) (ExprStmt x))"),
    ]);
}

#[test]
fn blocks_and_procedures() {
    check_trees(&[
        (
            "if a then b end",
            "(Root (ExprStmt (IfExpr if a then (Body (ExprStmt b)) end)))",
        ),
        (
            "if a then b elif c then d else e end",
            "(Root (ExprStmt (IfExpr if a then (Body (ExprStmt b)) elif c then (Body (ExprStmt d)) else (Body (ExprStmt e)) end)))",
        ),
        (
            "y := if x > 0 then x else -x end",
            "(Root (Definition y := (IfExpr if (Comparison x > 0) then (Body (ExprStmt x)) else (Body (ExprStmt (PrefixExpr - x))) end)))",
        ),
        (
            "if a then\n  b\n  c\nend",
            "(Root (ExprStmt (IfExpr if a then (Body (ExprStmt b) (ExprStmt c)) end)))",
        ),
        (
            "for k in 0..n do s := s + k end",
            "(Root (ForStmt for k in (RangeExpr 0 .. n) do (Body (Definition s := (BinaryExpr s + k))) end))",
        ),
        (
            "while m != 1 do m := m/2 end",
            "(Root (WhileStmt while (Comparison m != 1) do (Body (Definition m := (BinaryExpr m / 2))) end))",
        ),
        (
            "begin local y := 1, z; return y end",
            "(Root (ExprStmt (BlockExpr begin (Body (LocalStmt local (LocalItem y := 1) , (LocalItem z)) ; (ReturnStmt return y)) end)))",
        ),
        (
            "while c do break; continue end",
            "(Root (WhileStmt while c do (Body (BreakStmt break) ; (ContinueStmt continue)) end))",
        ),
        // A block inside brackets still separates its statements by line.
        (
            "f(begin\n a\n b\nend)",
            "(Root (ExprStmt (Application f (ArgList ( (BlockExpr begin (Body (ExprStmt a) (ExprStmt b)) end) )))))",
        ),
        (
            "if a then end",
            "(Root (ExprStmt (IfExpr if a then (Body) end)))",
        ),
    ]);
    let collatz = "collatz(n) := begin
  m := n
  steps := 0
  while m != 1 do
    m := if even(m) then m/2 else 3m + 1 end
    steps := steps + 1
  end
  steps
end";
    assert_eq!(
        tree(collatz),
        "(Root (Definition (Application collatz (ArgList ( n ))) := (BlockExpr begin (Body \
         (Definition m := n) (Definition steps := 0) (WhileStmt while (Comparison m != 1) do (Body \
         (Definition m := (IfExpr if (Application even (ArgList ( m ))) then (Body (ExprStmt \
         (BinaryExpr m / 2))) else (Body (ExprStmt (BinaryExpr (ImplicitProduct 3 m) + 1))) end)) \
         (Definition steps := (BinaryExpr steps + 1))) end) (ExprStmt steps)) end)))"
    );
}

#[test]
fn errors_and_recovery() {
    let cases: &[(&str, &[(&str, E)])] = &[
        ("a = b = c", &[("=", E::ChainedEquation)]),
        ("a = b < c", &[("<", E::ChainedEquation)]),
        ("a < b = c", &[("=", E::ChainedEquation)]),
        ("1..2..3", &[("..", E::ChainedRange)]),
        ("(a, b)", &[(",", E::Tuple)]),
        ("x 2", &[("2", E::NumberInProduct)]),
        ("2 3", &[("3", E::NumberInProduct)]),
        ("a +", &[("", E::ExpectedExpression)]),
        (
            "a + )",
            &[(")", E::ExpectedExpression), (")", E::ExpectedSeparator)],
        ),
        ("f(a b c", &[("", E::Expected(K::RParen))]),
        (
            "f(a ; b)",
            &[(";", E::Expected(K::RParen)), (")", E::ExpectedSeparator)],
        ),
        ("f(* , a)", &[("*", E::ExpectedExpression)]),
        (
            "[1, 2)",
            &[(")", E::Expected(K::RBracket)), (")", E::ExpectedSeparator)],
        ),
        ("()", &[(")", E::ExpectedExpression)]),
        ("a b c d e )", &[(")", E::ExpectedSeparator)]),
        ("x := ", &[("", E::ExpectedExpression)]),
        ("if a then b", &[("", E::Expected(K::EndKw))]),
        // `a b` is a condition, a product.
        ("if a b end", &[("end", E::Expected(K::ThenKw))]),
        ("if a\n  b\nend", &[("b", E::Expected(K::ThenKw))]),
        (
            "for 1 in s do end",
            &[("1", E::ExpectedName), ("1", E::Expected(K::InKw))],
        ),
        ("end", &[("end", E::Unexpected(K::EndKw))]),
        (")", &[(")", E::Unexpected(K::RParen))]),
        (
            "n :: 2",
            &[("2", E::ExpectedName), ("2", E::NumberInProduct)],
        ),
        ("(if a then b)", &[(")", E::Expected(K::EndKw))]),
        ("a b; c d", &[]),
    ];
    let mut wrong = Vec::new();
    for (text, expected) in cases {
        let found = parse_errors(text);
        if found != *expected {
            wrong.push(format!("{text:?}: found {found:?}"));
        }
        assert_eq!(parse(text).syntax().text().to_string(), *text);
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // Lexing errors come along, in order, and do not cascade.
    let parse = parse("f'(x) + |y|");
    let kinds: Vec<_> = parse.errors().iter().map(|e| e.kind).collect();
    assert!(
        kinds.iter().all(|k| matches!(k, SyntaxErrorKind::Lex(_))),
        "{kinds:?}"
    );
    assert_eq!(kinds.len(), 3);
}

/// A bad statement does not disturb the lines around it.
#[test]
fn recovery_is_local() {
    let text = "a := 1\nb := ) 2\nc := 3";
    assert_eq!(
        sexp(&parse(text).syntax()),
        "(Root (Definition a := 1) (Definition b :=) (ErrorNode ) 2) (Definition c := 3))"
    );
    assert_eq!(
        parse_errors(text),
        [(")", E::ExpectedExpression), (")", E::ExpectedSeparator)]
    );
}

#[test]
fn nesting_up_to_the_limit() {
    // The test threads have small stacks, and debug frames are large.
    let deep = |n: usize| format!("{}x{}", "(".repeat(n), ")".repeat(n));
    let ok = deep(MAX_DEPTH - 1);
    assert!(parse(&ok).errors().is_empty());
    for text in [
        deep(MAX_DEPTH + 10),
        "-".repeat(MAX_DEPTH + 10) + "x",
        "not ".repeat(MAX_DEPTH + 10) + "p",
        "x^".repeat(MAX_DEPTH + 10) + "x",
        "x -> ".repeat(MAX_DEPTH + 10) + "x",
        "begin ".repeat(MAX_DEPTH + 10) + &"end ".repeat(MAX_DEPTH + 10),
        "[".repeat(MAX_DEPTH + 10),
        // Postfix forms nest what comes before them (CAS-163).
        "x".to_owned() + &"!".repeat(MAX_DEPTH + 10),
        "f".to_owned() + &"(1)".repeat(MAX_DEPTH + 10),
        "a".to_owned() + &"[1]".repeat(MAX_DEPTH + 10),
        "x".to_owned() + &"...".repeat(MAX_DEPTH + 10),
        "x".to_owned() + &"::t".repeat(MAX_DEPTH + 10),
    ] {
        let parse = parse(&text);
        assert_eq!(parse.syntax().text().to_string(), text);
        let errors: Vec<_> = parse.errors().iter().map(|e| e.kind).collect();
        assert_eq!(
            errors,
            [SyntaxErrorKind::Parse(E::TooDeep)],
            "{}",
            &text[..20]
        );
    }
}

/// Pieces of the language, so that random concatenations exercise the
/// parser's decisions.
fn fragment() -> impl Strategy<Value = &'static str> {
    prop::sample::select(vec![
        "x", "f", "2", "0.5", "\"s\"", " ", "\n", ";", "(", ")", "[", "]", "{", "}", ",", "+", "-",
        "*", "/", "^", "²", "!", "...", "::", ":=", "->", "=", "==", "<", "..", "and", "or", "not",
        "if", "then", "elif", "else", "end", "for", "in", "do", "while", "begin", "local",
        "return", "break", "when", "√", "'", "|", "@",
    ])
}

/// Whether every byte of `text` is in the tree and every error lies in it.
fn check_lossless(text: &str) {
    let parse = parse(text);
    let root = parse.syntax();
    assert_eq!(root.text().to_string(), text);
    assert_eq!(root.kind(), K::Root);
    let mut last = 0;
    for e in parse.errors() {
        assert!(
            last <= e.range.start && e.range.start <= e.range.end,
            "{text:?}"
        );
        assert!(e.range.end <= text.len());
        last = e.range.start;
    }
    // In a text without errors, nodes span from their first to their last
    // significant token: no node but the root starts or ends with trivia.
    // (Where an operand is missing, trivia already added can end a node.)
    if !parse.errors().is_empty() {
        return;
    }
    for node in root.descendants() {
        let first = node.first_token();
        let last = node.last_token();
        if let (Some(first), Some(last)) = (first, last) {
            if node.kind() != K::Root {
                assert!(
                    !first.kind().is_trivia() && first.kind() != K::Newline,
                    "{text:?}"
                );
                assert!(
                    !last.kind().is_trivia() && last.kind() != K::Newline,
                    "{text:?}"
                );
            }
        }
    }
}

/// Lines from a pool of statements, joined by line breaks.
fn program() -> impl Strategy<Value = Vec<&'static str>> {
    prop::collection::vec(
        prop::sample::select(vec![
            "x := 1",
            "f(x) := x^2 when x > 0",
            "y := [1, 2, 3]",
            "a + b*c",
            "z := if a then b else c end",
            "g := (x, y) -> x y",
            "n! + 2x^2(x + 1)",
            "s := sum(k^2, k, 1, n)",
        ]),
        1..6,
    )
}

proptest! {
    #[test]
    fn lossless(fragments in prop::collection::vec(fragment(), 0..40)) {
        check_lossless(&fragments.concat());
    }

    #[test]
    fn lossless_any_text(text in any::<String>()) {
        check_lossless(&text);
    }

    /// A stray token in one line gives errors in that line only, and the
    /// other lines parse as they do alone.
    #[test]
    fn errors_stay_in_their_line(
        lines in program(),
        line in any::<prop::sample::Index>(),
        at in any::<prop::sample::Index>(),
        junk in prop::sample::select(vec![")", "]", "}", "*", ":=", "end", "then", ",", "@", "'", "else"]),
    ) {
        let bad = line.index(lines.len());
        let mut text = String::new();
        let mut range = 0..0;
        for (i, l) in lines.iter().enumerate() {
            if i > 0 {
                text.push('\n');
            }
            if i == bad {
                // Between two tokens of the line: not inside one, and not at
                // its ends, where a token that cannot start a statement
                // joins the line to its neighbour by design (CAS-DOC-7 §4).
                let lexed = crate::lex(l);
                let cuts: Vec<usize> = (1..lexed.len()).map(|i| lexed.range(i).start).collect();
                let cut = cuts[at.index(cuts.len())];
                range.start = text.len();
                write!(text, "{} {junk} {}", &l[..cut], &l[cut..]).unwrap();
                range.end = text.len();
            } else {
                text.push_str(l);
            }
        }
        let parse = parse(&text);
        // A stray token may make sense where it lands (`f * (x)`); if it does
        // not, the errors are in its line.
        for e in parse.errors() {
            prop_assert!(
                range.start <= e.range.start && e.range.end <= range.end,
                "{:?} outside line {} of {:?}", e, bad, text
            );
        }
        // Every other line keeps its tree.
        let good: Vec<String> = lines
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != bad)
            .map(|(_, l)| sexp(&parse_one(l)))
            .collect();
        let statements: Vec<String> = parse.syntax().children().map(|n| sexp(&n)).collect();
        for g in &good {
            prop_assert!(statements.contains(g), "{} missing from {:?}", g, statements);
        }
    }
}

/// The statement node of a line that parses alone.
fn parse_one(line: &str) -> SyntaxNode {
    parse(line).syntax().children().next().unwrap()
}
