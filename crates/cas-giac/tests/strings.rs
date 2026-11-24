//! Xcas's string commands against giac (CAS-611, CAS-613): giac's answers
//! read back are this CAS's value of the same input; departures, giac's
//! answer as caseval prints it, whitespace collapsed (C0007, C0047), and
//! this CAS's own.

use cas_core::{Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Each input: giac's answer, read back, is this CAS's value.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `char`, `asc` and `ord` on ASCII (probe1, probe5, probe8).
#[test]
fn codes() {
    check(&[
        "char(65)",
        "chr(65)",
        "char([65,66,67])",
        "char([])",
        "asc(\"bonjour\")",
        "asc(\"A\")",
        "asc(\"\")",
        "asc(\"abc\",1)",
        "asc(5)",
        "ord(\"A\")",
        "ord(\"ABC\")",
        "ord([\"a\",\"b\",\"c\"])",
        "ord([])",
        "ord([\"ab\",[\"c\"]])",
    ]);
}

/// `inString` and `contains`, with giac's interval `[lo..hi]` (probe1,
/// probe5, probe8).
#[test]
fn positions() {
    check(&[
        "inString(\"abcd\",\"b\")",
        "inString(\"abcd\",\"e\")",
        "inString(\"abc\",\"bc\")",
        "inString(\"\",\"a\")",
        "inString(\"abcb\",\"b\")",
        "inString(\"abc\",\"\",5)",
        "inString(\"abc\",\"\")",
        "inString(\"abcb\",\"b\",2)",
        "inString(\"abcb\",\"b\",-1)",
        "inString(\"abc\",\"c\",3)",
        "inString(\"abc\",\"\",3)",
        "contains([0,1,2,3],4)",
        "contains([1,2,1],1)",
        "contains(\"abc\",\"b\")",
        "contains([[1,2]],[1,2])",
        "contains([1..2],1.5)",
        "contains([1..2],3)",
        "contains([1..2],1)",
        "contains([1..2],x)",
        "contains([1,2..3],2.5)",
        "contains([1,2..3],5/2)",
        "contains(\"abc\",\"\")",
        "contains([1,2],[1,2])",
        "contains([1..2],1.7)",
        "contains([0..2],1.5)",
        "contains([1..3],2.5)",
        "contains([5,1..2],1.5)",
        "contains([1..2],2)",
        "contains([1..2],2.5)",
        "contains([1..2],3/2)",
        "contains([x],x)",
    ]);
}

/// `mid`, `left`, `right` and their aliases `lhs`, `rhs` (probe1, probe5,
/// probe8, probe10).
#[test]
fn pieces() {
    check(&[
        "mid(\"azertyuiop\",2,4)",
        "mid(\"azertyuiop\",2)",
        "mid([0,1,2,3,4,5,6],2,3)",
        "mid([0,1,2,3,4,5,6],2)",
        "mid([[1,2],[3,4],[5,6]],1)",
        "mid(\"abc\",5)",
        "mid(\"abc\",-1)",
        "mid(\"abc\",1,10)",
        "mid(\"abc\",0,0)",
        "mid([1,2,3],1,-1)",
        "mid(\"abc\",1,-1)",
        "mid(\"abc\",3)",
        "mid([1,2,3],3)",
        "mid(\"abc\",1,1)",
        "left(\"abcdefg\",3)",
        "right(\"abcdefg\",3)",
        "left(\"abc\",0)",
        "left(\"abc\",10)",
        "left(\"abc\",-1)",
        "right(\"abc\",10)",
        "left([1,2],-1)",
        "right([1,2,3],-1)",
        "right(\"abc\",-1)",
        "right(\"abc\",0)",
        "left(a=b)",
        "right(a=b)",
        "lhs(a=b)",
        "rhs(a=b)",
        "left(x^2+1=5)",
        "right(x^2+1=5)",
        "left(1..5)",
        "right(1..5)",
    ]);
}

/// `levenshtein` (probe5, probe8), and `hamming` where giac and ours stay.
#[test]
fn distances() {
    check(&[
        "levenshtein(\"kitten\",\"sitting\")",
        "levenshtein(\"\",\"abc\")",
        "levenshtein([],[])",
        "levenshtein([1,2],[2,1])",
        "hamming(\"ab\",\"abc\")",
    ]);
}

/// Departures, giac's answer as caseval prints it (probe5, probe9,
/// probe10) and ours: giac counts bytes of UTF-8 where ours counts Unicode
/// scalar values; giac's errors stay here; giac has no `hamming`.
#[test]
fn departures() {
    let cases = [
        ("asc(\"\u{e9}\")", "[195,169]", Some("[233]")),
        ("ord(\"\u{e9}\")", "-61", Some("233")),
        ("char(233)", "\"\u{fffd}\"", Some("\"\u{e9}\"")),
        (
            "char(65.5)",
            "GIAC_ERROR: char(65.5) \n Error: Bad Argument Value",
            None,
        ),
        (
            "ord(\"\")",
            "GIAC_ERROR: ord(\"\") \n Error: Bad Argument Value",
            None,
        ),
        (
            "levenshtein(\"abc\",[1,2])",
            "GIAC_ERROR: Both arguments must be either lists or strings. Error: Bad Argument Type",
            None,
        ),
        (
            "hamming(\"cats\",\"dogs\")",
            "hamming(\"cats\",\"dogs\")",
            Some("3"),
        ),
        (
            "hamming([1,0,1,1],[1,1,1,0])",
            "hamming([1,0,1,1],[1,1,1,0])",
            Some("2"),
        ),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, mine) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours(input);
        let ok = match mine {
            Some(mine) => found.as_ref().ok() == ours(mine).ok().as_ref(),
            None => found
                .as_ref()
                .ok()
                .and_then(|e| e.head().and_then(Expr::as_symbol))
                .is_some_and(|head| {
                    [Symbol::CHAR, Symbol::ORD, Symbol::LEVENSHTEIN].contains(&head)
                }),
        };
        if !ok {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
