//! Tests of Xcas's string commands (CAS-611, CAS-613): every row of the
//! phase 1.6 prototype, giac's answers where giac answers, and the calls
//! that stay.

use crate::trace::{EvalOptions, Warning};
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input in full form keeps its head, with no warning (C0056).
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

#[test]
fn codes() {
    check(&[
        ("Char[65]", "\"A\""),
        ("Char[List[65, 66, 67]]", "\"ABC\""),
        ("Char[List[]]", "\"\""),
        ("Char[233]", "\"\u{e9}\""),
        ("Asc[\"bonjour\"]", "List[98, 111, 110, 106, 111, 117, 114]"),
        ("Asc[\"A\"]", "List[65]"),
        ("Asc[\"\"]", "List[]"),
        ("Asc[\"abc\", 1]", "98"),
        ("Asc[\"\u{e9}\"]", "List[233]"),
        ("Ord[\"A\"]", "65"),
        ("Ord[\"ABC\"]", "65"),
        ("Ord[List[\"a\", \"b\", \"c\"]]", "List[97, 98, 99]"),
        ("Ord[List[]]", "List[]"),
        ("Ord[List[\"ab\", List[\"c\"]]]", "List[97, List[99]]"),
        ("Ord[\"\u{e9}\"]", "233"),
    ]);
    stays(&[
        "Char[0x1.06p+6]",
        "Char[-1]",
        "Char[55296]",
        "Char[1114112]",
        "Char[List[65, x]]",
        "Char[65, 66]",
        "Asc[5]",
        "Asc[\"abc\", 3]",
        "Asc[\"abc\", -1]",
        "Ord[\"\"]",
        "Ord[List[\"a\", \"\"]]",
        "Ord[x]",
    ]);
}

#[test]
fn positions() {
    check(&[
        ("InString[\"abcd\", \"b\"]", "1"),
        ("InString[\"abcd\", \"e\"]", "-1"),
        ("InString[\"abc\", \"bc\"]", "1"),
        ("InString[\"\", \"a\"]", "-1"),
        ("InString[\"abcb\", \"b\"]", "1"),
        ("InString[\"abc\", \"\", 5]", "-1"),
        ("InString[\"abc\", \"\"]", "0"),
        ("InString[\"abcb\", \"b\", 2]", "3"),
        ("InString[\"abcb\", \"b\", -1]", "-1"),
        ("InString[\"abc\", \"c\", 3]", "-1"),
        ("InString[\"abc\", \"\", 3]", "3"),
        ("InString[\"h\u{e9}llo\", \"l\"]", "2"),
    ]);
    stays(&["InString[\"abc\", x]", "InString[\"abc\", \"b\", 0x1p+0]"]);
}

#[test]
fn containing() {
    check(&[
        ("Contains[List[0, 1, 2, 3], 4]", "0"),
        ("Contains[List[1, 2, 1], 1]", "1"),
        ("Contains[\"abc\", \"b\"]", "2"),
        ("Contains[List[List[1, 2]], List[1, 2]]", "1"),
        ("Contains[List[Range[1, 2]], 0x1.8p+0]", "1"),
        ("Contains[List[Range[1, 2]], 3]", "0"),
        ("Contains[List[Range[1, 2]], 1]", "1"),
        ("Contains[List[Range[1, 2]], \"x\"]", "0"),
        ("Contains[List[Range[1, 2]], x]", "0"),
        ("Contains[List[1, Range[2, 3]], 0x1.4p+1]", "0"),
        ("Contains[\"abc\", \"\"]", "1"),
        ("Contains[List[1, 2], List[1, 2]]", "0"),
        ("Contains[List[Range[1, 2]], 0x1.b333333333333p+0]", "1"),
        ("Contains[List[Range[0, 2]], 0x1.8p+0]", "1"),
        ("Contains[List[Range[1, 3]], 0x1.4p+1]", "1"),
        ("Contains[List[5, Range[1, 2]], 0x1.8p+0]", "0"),
        ("Contains[List[Range[1, 2]], 2]", "1"),
        ("Contains[List[Range[1, 2]], 0x1.4p+1]", "0"),
        ("Contains[List[Range[1, 2]], 3/2]", "1"),
        ("Contains[List[\"x\"], \"x\"]", "1"),
        ("Contains[List[x], x]", "1"),
        ("Contains[List[Range[a, b]], Range[a, b]]", "1"),
    ]);
    stays(&[
        "Contains[5, 1]",
        "Contains[List[Range[1, 2]], List[1]]",
        "Contains[List[Range[1, 2]], Range[1, 2]]",
    ]);
}

#[test]
fn pieces() {
    check(&[
        ("Mid[\"azertyuiop\", 2, 4]", "\"erty\""),
        ("Mid[\"azertyuiop\", 2]", "\"ertyuiop\""),
        ("Mid[List[0, 1, 2, 3, 4, 5, 6], 2, 3]", "List[2, 3, 4]"),
        ("Mid[List[0, 1, 2, 3, 4, 5, 6], 2]", "List[2, 3, 4, 5, 6]"),
        (
            "Mid[List[List[1, 2], List[3, 4], List[5, 6]], 1]",
            "List[List[3, 4], List[5, 6]]",
        ),
        ("Mid[\"abc\", 5]", "\"\""),
        ("Mid[\"abc\", -1]", "\"\""),
        ("Mid[\"abc\", 1, 10]", "\"bc\""),
        ("Mid[\"abc\", 0, 0]", "\"\""),
        ("Mid[List[1, 2, 3], 1, -1]", "List[]"),
        ("Mid[\"abc\", 1, -1]", "\"\""),
        ("Mid[\"abc\", 3]", "\"\""),
        ("Mid[List[1, 2, 3], 3]", "List[]"),
        ("Mid[\"abc\", 1, 1]", "\"b\""),
        ("Left[\"abcdefg\", 3]", "\"abc\""),
        ("Right[\"abcdefg\", 3]", "\"efg\""),
        ("Left[\"abc\", 0]", "\"\""),
        ("Left[\"abc\", 10]", "\"abc\""),
        ("Left[\"abc\", -1]", "\"abc\""),
        ("Right[\"abc\", 10]", "\"abc\""),
        ("Left[List[1, 2, 3, 4], 2]", "List[1, 2]"),
        ("Right[List[1, 2, 3, 4], 2]", "List[3, 4]"),
        ("Left[Equation[\"a\", \"b\"]]", "\"a\""),
        ("Right[Equation[Plus[1, Power[x, 2]], 5]]", "5"),
        ("Left[Equal[a, b]]", "a"),
        ("Left[Range[1, 5]]", "1"),
        ("Right[Range[1, 5]]", "5"),
        ("Left[List[1, 2], -1]", "List[]"),
        ("Right[List[1, 2, 3], -1]", "List[]"),
        ("Right[\"abc\", -1]", "\"\""),
        ("Right[\"abc\", 0]", "\"\""),
    ]);
    stays(&[
        "Right[Range[1, 5], 2]",
        "Left[List[1, 2, 3]]",
        "Left[\"abc\", 0x1p+1]",
        "Mid[\"abc\", 0x1p+0]",
        "Mid[x, 1]",
        "Right[x]",
    ]);
}

#[test]
fn distances() {
    check(&[
        ("Levenshtein[\"kitten\", \"sitting\"]", "3"),
        ("Levenshtein[List[1, 2, 3], List[1, 3]]", "1"),
        ("Levenshtein[\"abc\", \"abc\"]", "0"),
        ("Levenshtein[\"\", \"abc\"]", "3"),
        ("Levenshtein[List[], List[]]", "0"),
        ("Levenshtein[List[1, 2], List[2, 1]]", "2"),
        ("Hamming[\"cats\", \"dogs\"]", "3"),
        ("Hamming[List[1, 0, 1, 1], List[1, 1, 1, 0]]", "2"),
    ]);
    stays(&[
        "Hamming[\"ab\", \"abc\"]",
        "Hamming[\"ab\", List[1, 2]]",
        "Levenshtein[x, \"abc\"]",
    ]);
    // A string and a list: giac's message as a warning.
    let ctx = Context::new();
    let e = read("Levenshtein[\"abc\", List[1, 2]]");
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap(), e);
    assert!(
        matches!(
            &report.warnings[..],
            [Warning::Native { message, .. }]
                if message == "Both arguments must be either lists or strings"
        ),
        "{:?}",
        report.warnings
    );
}
