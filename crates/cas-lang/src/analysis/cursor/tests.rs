use std::fmt::Write as _;

use super::{Call, Cursor, CursorKind, cursor, token_at};
use crate::parser::parse;

/// The cursor as `Kind range`, then `call callee arguments argument` when
/// there is a call, the callee as its text.
fn describe(text: &str, cursor: &Cursor) -> String {
    let mut found = match &cursor.kind {
        CursorKind::Name { range } => format!("Name {range:?}"),
        CursorKind::Type { range } => format!("Type {range:?}"),
        kind => format!("{kind:?}"),
    };
    if let Some(call) = &cursor.call {
        let callee = &text[call.callee.clone()];
        write!(
            found,
            " call {callee} {:?} {}",
            call.arguments, call.argument
        )
        .unwrap();
    }
    found
}

/// The cursor at the `|` of `marked`, which is not part of the text.
fn at(marked: &str) -> String {
    let offset = marked.find('|').unwrap();
    let text = marked.replacen('|', "", 1);
    let found = cursor(&parse(&text), offset);
    assert_eq!(found.offset, offset, "{marked:?}");
    describe(&text, &found)
}

/// The token's text at the `|` of `marked`.
fn token(marked: &str) -> String {
    let offset = marked.find('|').unwrap();
    let text = marked.replacen('|', "", 1);
    token_at(&parse(&text), offset).unwrap().text().to_owned()
}

/// The cursor at every offset of T4's half-typed input (CAS-665), alone and
/// with lines below it: where the tokens start and end and between them.
/// After an open `(` the call runs to the first significant token past its
/// arguments; in `f(g(x, ` with lines below, `f`'s arguments run to the end.
#[test]
fn half_typed_input() {
    let mut found = String::new();
    for start in [
        "f(x, ",
        "sin(",
        "x := ",
        "begin x := 1",
        "[1, ",
        "f(x) := ",
        "f(g(x, ",
    ] {
        for text in [start.to_owned(), format!("{start}\ny := 1\nz := 2")] {
            let parse = parse(&text);
            writeln!(found, "{text:?}").unwrap();
            for offset in 0..=text.len() {
                let c = cursor(&parse, offset);
                writeln!(found, "  {offset} {}", describe(&text, &c)).unwrap();
            }
        }
    }
    assert_eq!(found, HALF_TYPED);
}

/// The cursor at each offset of each text of [`half_typed_input`].
const HALF_TYPED: &str = r#""f(x, "
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..4 0
  3 Name 2..3 call f 1..4 0
  4 Name 4..4 call f 1..4 1
  5 Name 5..5 call f 1..4 1
"f(x, \ny := 1\nz := 2"
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..4 0
  3 Name 2..3 call f 1..4 0
  4 Name 4..4 call f 1..4 1
  5 Name 5..5 call f 1..4 1
  6 Name 6..7 call f 1..4 1
  7 Name 6..7
  8 Name 8..8
  9 Other
  10 Name 10..10
  11 Number
  12 Number
  13 Name 13..14
  14 Name 13..14
  15 Name 15..15
  16 Other
  17 Name 17..17
  18 Number
  19 Number
"sin("
  0 Name 0..3
  1 Name 0..3
  2 Name 0..3
  3 Name 0..3
  4 Name 4..4 call sin 3..4 0
"sin(\ny := 1\nz := 2"
  0 Name 0..3
  1 Name 0..3
  2 Name 0..3
  3 Name 0..3
  4 Name 4..4 call sin 3..4 0
  5 Name 5..6 call sin 3..4 0
  6 Name 5..6
  7 Name 7..7
  8 Other
  9 Name 9..9
  10 Number
  11 Number
  12 Name 12..13
  13 Name 12..13
  14 Name 14..14
  15 Other
  16 Name 16..16
  17 Number
  18 Number
"x := "
  0 Name 0..1
  1 Name 0..1
  2 Name 2..2
  3 Other
  4 Name 4..4
  5 Name 5..5
"x := \ny := 1\nz := 2"
  0 Name 0..1
  1 Name 0..1
  2 Name 2..2
  3 Other
  4 Name 4..4
  5 Name 5..5
  6 Name 6..7
  7 Name 6..7
  8 Name 8..8
  9 Other
  10 Name 10..10
  11 Number
  12 Number
  13 Name 13..14
  14 Name 13..14
  15 Name 15..15
  16 Other
  17 Name 17..17
  18 Number
  19 Number
"begin x := 1"
  0 Name 0..5
  1 Name 0..5
  2 Name 0..5
  3 Name 0..5
  4 Name 0..5
  5 Name 0..5
  6 Name 6..7
  7 Name 6..7
  8 Name 8..8
  9 Other
  10 Name 10..10
  11 Number
  12 Number
"begin x := 1\ny := 1\nz := 2"
  0 Name 0..5
  1 Name 0..5
  2 Name 0..5
  3 Name 0..5
  4 Name 0..5
  5 Name 0..5
  6 Name 6..7
  7 Name 6..7
  8 Name 8..8
  9 Other
  10 Name 10..10
  11 Number
  12 Number
  13 Name 13..14
  14 Name 13..14
  15 Name 15..15
  16 Other
  17 Name 17..17
  18 Number
  19 Number
  20 Name 20..21
  21 Name 20..21
  22 Name 22..22
  23 Other
  24 Name 24..24
  25 Number
  26 Number
"[1, "
  0 Name 0..0
  1 Number
  2 Number
  3 Name 3..3
  4 Name 4..4
"[1, \ny := 1\nz := 2"
  0 Name 0..0
  1 Number
  2 Number
  3 Name 3..3
  4 Name 4..4
  5 Name 5..6
  6 Name 5..6
  7 Name 7..7
  8 Other
  9 Name 9..9
  10 Number
  11 Number
  12 Name 12..13
  13 Name 12..13
  14 Name 14..14
  15 Other
  16 Name 16..16
  17 Number
  18 Number
"f(x) := "
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..4 0
  3 Name 2..3 call f 1..4 0
  4 Name 4..4
  5 Name 5..5
  6 Other
  7 Name 7..7
  8 Name 8..8
"f(x) := \ny := 1\nz := 2"
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..4 0
  3 Name 2..3 call f 1..4 0
  4 Name 4..4
  5 Name 5..5
  6 Other
  7 Name 7..7
  8 Name 8..8
  9 Name 9..10
  10 Name 9..10
  11 Name 11..11
  12 Other
  13 Name 13..13
  14 Number
  15 Number
  16 Name 16..17
  17 Name 16..17
  18 Name 18..18
  19 Other
  20 Name 20..20
  21 Number
  22 Number
"f(g(x, "
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..6 0
  3 Name 2..3 call f 1..6 0
  4 Name 4..5 call g 3..6 0
  5 Name 4..5 call g 3..6 0
  6 Name 6..6 call g 3..6 1
  7 Name 7..7 call g 3..6 1
"f(g(x, \ny := 1\nz := 2"
  0 Name 0..1
  1 Name 0..1
  2 Name 2..3 call f 1..21 0
  3 Name 2..3 call f 1..21 0
  4 Name 4..5 call g 3..6 0
  5 Name 4..5 call g 3..6 0
  6 Name 6..6 call g 3..6 1
  7 Name 7..7 call g 3..6 1
  8 Name 8..9 call g 3..6 1
  9 Name 8..9 call f 1..21 0
  10 Name 10..10 call f 1..21 0
  11 Other call f 1..21 0
  12 Name 12..12 call f 1..21 0
  13 Number call f 1..21 0
  14 Number call f 1..21 0
  15 Name 15..16 call f 1..21 0
  16 Name 15..16 call f 1..21 0
  17 Name 17..17 call f 1..21 0
  18 Other call f 1..21 0
  19 Name 19..19 call f 1..21 0
  20 Number call f 1..21 0
  21 Number call f 1..21 0
"#;

#[test]
fn calls() {
    // CAS-DOC-60 §3.
    assert_eq!(at("f(x, |"), "Name 5..5 call f 1..4 1");
    assert_eq!(at("sin(|"), "Name 4..4 call sin 3..4 0");
    // The innermost call; commas of a list are not the call's.
    assert_eq!(at("f(g(x, |"), "Name 7..7 call g 3..6 1");
    assert_eq!(at("f([1, 2], |"), "Name 10..10 call f 1..9 1");
    assert_eq!(at("f(a)(b, |"), "Name 8..8 call f(a) 4..7 1");
    assert_eq!(at("call(f, |"), "Name 8..8 call call 4..7 1");
    // An index is not a call.
    assert_eq!(at("v[1, |"), "Name 5..5");
    // Closed: inside up to the `)`, outside after it and before the `(`.
    assert_eq!(at("f(x|)"), "Name 2..3 call f 1..4 0");
    assert_eq!(at("f(x)|"), "Name 4..4");
    assert_eq!(at("f|(x)"), "Name 0..1");
    // Inside a string.
    assert_eq!(at("f(\"a|"), "String call f 1..4 0");
}

#[test]
fn kinds() {
    assert_eq!(at("x := |"), "Name 5..5");
    assert_eq!(at("begin x := 1|"), "Number");
    assert_eq!(at("x|²"), "Name 0..1");
    assert_eq!(at("x²|"), "Number");
    // A type after `::`, empty or begun.
    assert_eq!(at("x :: |"), "Type 5..5");
    assert_eq!(at("x :: re|"), "Type 5..7");
    // A keyword is a name.
    assert_eq!(at("in|"), "Name 0..2");
    // Inside an operator.
    assert_eq!(at(":|="), "Other");
    // Unterminated strings and comments hold the end; closed ones do not.
    assert_eq!(at("\"ab|"), "String");
    assert_eq!(at("\"a|b\""), "String");
    assert_eq!(at("\"ab\"|"), "Name 4..4");
    assert_eq!(at("|\"ab\""), "Name 0..0");
    assert_eq!(at("/* a /* b */ |"), "Comment");
    assert_eq!(at("/* a */|"), "Name 7..7");
    assert_eq!(at("// c|"), "Comment");
    assert_eq!(at("// c|\nx"), "Comment");
}

#[test]
fn tokens() {
    assert_eq!(token("sin|("), "sin");
    assert_eq!(token("(|x"), "x");
    assert_eq!(token("a|+b"), "a");
    assert_eq!(token("π|r"), "π");
    assert_eq!(token("|x"), "x");
    assert_eq!(token("x + y|"), "y");
    assert_eq!(token("x |+"), "+");
    assert_eq!(token("x|"), "x");
    assert!(token_at(&parse(""), 0).is_none());
}

#[test]
fn offsets() {
    let name = |offset, range| Cursor {
        offset,
        kind: CursorKind::Name { range },
        call: None,
    };
    let x = parse("x");
    for offset in [2, 99, usize::MAX] {
        assert_eq!(cursor(&x, offset), name(1, 0..1));
        assert_eq!(token_at(&x, offset).unwrap().text(), "x");
    }
    // Inside `é`: its start.
    let e = parse("é");
    assert_eq!(cursor(&e, 1), name(0, 0..2));
    assert_eq!(token_at(&e, 1).unwrap().text(), "é");
    assert_eq!(cursor(&e, 3), name(2, 0..2));
    let empty = parse("");
    assert_eq!(cursor(&empty, usize::MAX), name(0, 0..0));
    let f = parse("f(é");
    assert_eq!(
        cursor(&f, 3),
        Cursor {
            offset: 2,
            kind: CursorKind::Name { range: 2..4 },
            call: Some(Call {
                callee: 0..1,
                arguments: 1..4,
                argument: 0
            }),
        }
    );
}
