//! Xcas's string commands against giac (CAS-611, CAS-612, CAS-613), and
//! strings in `size`, `concat`, `[]`, the comparisons and `sort`: giac's answers
//! read back are this CAS's value of the same input; departures, giac's
//! answer as caseval prints it, whitespace collapsed (C0007, C0047), and
//! this CAS's own; and random ASCII strings a run to each command but
//! `strip`, `size`, `concat`, `[]`, the comparisons and `sort`.

use cas_core::{Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

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

/// `split`, `join`, `trim` and `strip` (probe1, probe5, probe8, probe9),
/// with giac's `split(e, [x, y])` of an expression.
#[test]
fn words() {
    check(&[
        "split(\"This is a phrase\",\" \")",
        "split(\"a,,b\",\",\")",
        "split(\"\",\",\")",
        "split(\"a,\",\",\")",
        "split(\",a\",\",\")",
        "split(\"aXYbXY\",\"XY\")",
        "split(x^3*y^2-y^2+x^3-1,[x,y])",
        "split(x^3*y^2-y^2+x^3+1,[x,y])",
        "split(2*x*y,[x,y])",
        "split(x+y,[x,y])",
        "split(3,[x,y])",
        "split(x^2*sin(y),[x,y])",
        "join(\" \",[\"This\",\"is\",\"a\",\"phrase\"])",
        "join(\",\",[])",
        "join(\"\",[\"a\",\"b\"])",
        "trim(\"   an apple     \")",
        "trim(\" a b \")",
        "trim(\"\")",
        "trim(\"ab\",left)",
        "trim(\" ab \",left)",
        "trim(\" ab \",right)",
        "trim([0,0,1,2,3,0,1e-16,0],right)",
        "trim([0.01,0.023,0.005,0.12,0.27,0.82,0.53,0.21,0.05,0.001],0.1)",
        "trim([0,0],right)",
        "trim([1e-13,1,0])",
        "trim([1e-11,1,0])",
        "trim([0,1,0],left)",
        "trim([0,1,0])",
        "trim([0,x,0])",
        "trim([0,0.0,1])",
        "trim([1,0.5,2],1)",
        "strip(\"\ta\t\")",
        "strip(\"xxaxx\",\"x\")",
        "strip(\"xxx\",\"x\")",
        "strip(\"  \")",
    ]);
}

/// Strings in `size`/`length`, `concat`, `[]`, the comparisons and `sort`
/// (probe5, probe8, probe9).
#[test]
fn elsewhere() {
    check(&[
        "size(\"abc\")",
        "size(\"\")",
        "length(\"abc\")",
        "concat(\"ab\",\"cd\")",
        "concat(\"a\",\"b\",\"c\")",
        "concat([\"a\"],\"b\")",
        // Non-strings printed (CAS-614).
        "join(\",\",[1,2])",
        "concat(\"a\",1)",
        "\"abc\"[1]",
        "\"abc\"[0]",
        "\"abc\"[-1]",
        "\"abc\"[1..2]",
        "\"abc\"==\"abc\"",
        "\"a\"==\"b\"",
        "\"b\">\"a\"",
        "\"a\"<\"b\"",
        "\"B\"<\"a\"",
        "\"b\">=\"a\"",
        "sort([\"b\",\"a\",\"c\"])",
        "sort([\"b\",\"A\",\"a\",\"aa\",\"\"])",
    ]);
}

/// An index out of a string's range: giac's error (probe8, probe9), and
/// ours; a slice past the end, which giac cuts short, is ours too, as for
/// lists.
#[test]
fn index_errors() {
    let error = "GIAC_ERROR: Gen [int] Error: Bad Argument Type";
    let cases = [
        ("\"abc\"[5]", error),
        ("\"abc\"[-4]", error),
        ("\"abc\"[0..5]", "\"abc\""),
    ];
    for (input, answer) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            let collapsed = printed.split_whitespace().collect::<Vec<_>>().join(" ");
            assert_eq!(collapsed, answer, "{input}");
        }
        assert!(ours(input).is_err(), "{input}");
    }
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
        // CAS-612: ours counts characters, giac bytes.
        ("length(\"h\u{e9}llo\")", "6", Some("5")),
        (
            "split(\"abc\",\"\")",
            "GIAC_ERROR: split(\"abc\",\"\") \n Error: Invalid dimension",
            None,
        ),
        (
            "trim([1,0.5,2],-1)",
            "GIAC_ERROR: Invalid threshold specification. Error: Bad Argument Value",
            None,
        ),
        (
            "trim([1,0.5,2],x)",
            "GIAC_ERROR: Invalid threshold specification. Error: Bad Argument Value",
            None,
        ),
        // `+` does not join strings.
        ("\"ab\"+\"cd\"", "\"abcd\"", None),
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
                    [
                        Symbol::CHAR,
                        Symbol::ORD,
                        Symbol::LEVENSHTEIN,
                        Symbol::SPLIT,
                        Symbol::TRIM,
                        Symbol::PLUS,
                    ]
                    .contains(&head)
                }),
        };
        if !ok {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Printable ASCII, often `a`, `b`, a space, `"` or `\`, at most `max`
/// characters.
fn ascii(max: usize) -> impl Strategy<Value = String> {
    let c = prop_oneof![0x20u8..=0x7e, prop::sample::select(&b"ab \"\\"[..])];
    prop::collection::vec(c, 0..=max).prop_map(|v| v.into_iter().map(char::from).collect())
}

/// `s` as a string literal: `"` and `\` escaped.
fn literal(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

/// Whether giac reads `s` written as [`literal`]: its lexer takes a string
/// ending in `\` for an unclosed one (`size("a\\")` is `undef`, caseval).
fn writable(s: &str) -> bool {
    !s.ends_with('\\')
}

/// Each input agrees ([`Giac::agrees`]), or the failure names it.
fn agree(giac: &Giac, inputs: &[String]) -> Result<(), TestCaseError> {
    for input in inputs {
        giac.agrees(input)
            .map_err(|e| TestCaseError::fail(format!("{input}: {e}")))?;
    }
    Ok(())
}

// giac prints a string's `"` doubled and its `\` bare (caseval:
// `"a\"b"` is `"a""b"`), which do not read back: commands answering a
// string are compared through `asc`.
proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `asc`, `ord` and `char` (through `asc`) on ASCII; `ord("")` and
    /// non-ASCII are departures.
    #[test]
    fn random_codes(
        s in ascii(12),
        t in ascii(4),
        codes in prop::collection::vec(0x20i64..=0x7e, 0..=4),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        if !writable(&s) || !writable(&t) {
            return Ok(());
        }
        let (s, t) = (literal(&s), literal(&t));
        let codes: Vec<String> = codes.iter().map(i64::to_string).collect();
        let mut inputs = vec![
            format!("asc({s})"),
            format!("asc(char([{}]))", codes.join(",")),
        ];
        if let Some(c) = codes.first() {
            inputs.push(format!("asc(char({c}))"));
        }
        if s != "\"\"" {
            inputs.push(format!("ord({s})"));
            if t != "\"\"" {
                inputs.push(format!("ord([{s},{t}])"));
            }
        }
        agree(giac, &inputs)?;
    }

    /// `inString`, with and without a start, and `contains`, of a random
    /// needle or a piece of the string.
    #[test]
    fn random_positions(
        s in ascii(12),
        other in ascii(2),
        piece in any::<bool>(),
        at in 0usize..=12,
        n in 0usize..=3,
        start in -1i64..=13,
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let t = if piece {
            s.chars().skip(at).take(n).collect()
        } else {
            other
        };
        if !writable(&s) || !writable(&t) {
            return Ok(());
        }
        let (s, t) = (literal(&s), literal(&t));
        agree(giac, &[
            format!("inString({s},{t})"),
            format!("inString({s},{t},{start})"),
            format!("contains({s},{t})"),
        ])?;
    }

    /// `mid`, `left` and `right` (through `asc`).
    #[test]
    fn random_pieces(s in ascii(12), a in -1i64..=13, n in -1i64..=13) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        if !writable(&s) {
            return Ok(());
        }
        let s = literal(&s);
        agree(giac, &[
            format!("asc(mid({s},{a}))"),
            format!("asc(mid({s},{a},{n}))"),
            format!("asc(left({s},{n}))"),
            format!("asc(right({s},{n}))"),
        ])?;
    }

    /// `levenshtein` of two strings.
    #[test]
    fn random_levenshtein(s in ascii(12), t in ascii(12)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        if !writable(&s) || !writable(&t) {
            return Ok(());
        }
        agree(giac, &[format!("levenshtein({},{})", literal(&s), literal(&t))])?;
    }

    /// `join` and `trim` (through `asc`), and `split` by a nonempty
    /// separator: `split(s, "")` is giac's error (a departure), and `split`
    /// answers a list of strings, so its strings hold no `"` or `\`.
    #[test]
    fn random_words(
        s in ascii(12),
        sep in ascii(2),
        words in prop::collection::vec(ascii(4), 0..=3),
        plain in "[ ab,x]{0,12}",
        by in "[ ab,x]{1,2}",
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        if !writable(&s) || !writable(&sep) || !words.iter().all(|w| writable(w)) {
            return Ok(());
        }
        let (s, sep) = (literal(&s), literal(&sep));
        let words: Vec<String> = words.iter().map(|w| literal(w)).collect();
        agree(giac, &[
            format!("asc(join({sep},[{}]))", words.join(",")),
            format!("asc(trim({s}))"),
            format!("asc(trim({s},left))"),
            format!("asc(trim({s},right))"),
            format!("split({},{})", literal(&plain), literal(&by)),
        ])?;
    }
}
