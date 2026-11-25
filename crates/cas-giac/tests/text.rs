//! Expressions as text and back against giac (CAS-614): `string`/`str`,
//! `expr`/`execute`, `cat`, and non-strings in `join` and `concat`. A text
//! is giac's input form there and this CAS's here (`x+1` against `1 + x`),
//! so text rows compare what the two texts read back to (C0032); `expr`
//! rows compare values. Departures: giac's answer as libgiac prints it,
//! whitespace collapsed (C0007, C0047), and this CAS's own.

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

/// The text of a string answer: giac's quotes removed and its doubled
/// inner quotes undone; ours, the string itself.
fn unquoted(answer: &str) -> Option<String> {
    let inner = answer.trim().strip_prefix('"')?.strip_suffix('"')?;
    Some(inner.replace("\"\"", "\""))
}

/// Whether two texts say the same: equal without their spaces, or read
/// back by this CAS to the same value (`==`, or `normal(a - b)` is 0).
fn same_text(giac: &str, mine: &str) -> bool {
    let bare = |s: &str| s.split_whitespace().collect::<String>();
    if bare(giac) == bare(mine) {
        return true;
    }
    let (Ok(a), Ok(b)) = (ours(giac), ours(mine)) else {
        return false;
    };
    a == b
        || ours(&format!("({giac}) == ({mine})")).is_ok_and(|e| e.as_symbol() == Some(Symbol::TRUE))
        || ours(&format!("normal(({giac}) - ({mine}))")).is_ok_and(|e| e == Expr::from(0))
}

/// Each input answers a string: giac's text and ours say the same.
fn check_text(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in inputs {
        let answer = giac.eval(input).unwrap_or_default();
        let mine = ours(input);
        let ok = match (unquoted(&answer), mine.as_ref().ok().and_then(Expr::as_str)) {
            (Some(theirs), Some(mine)) => same_text(&theirs, mine),
            _ => false,
        };
        if !ok {
            wrong.push(format!("{input}: giac {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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

/// `string` and `str` (probe1, probe5, probe11).
#[test]
fn strings_of_expressions() {
    check_text(&[
        "string(1.23)",
        "string(x^2+1)",
        "string([1,2])",
        "string(\"a\")",
        "string(1/3)",
        "string(sqrt(2))",
        "string(-x)",
        "string(2*x^2-3*x+1/2)",
        // `a` is assigned by giac's `expr("a:=3")` below.
        "string(u=v)",
        "string(sin(x)/2)",
        "string(-1/2*x)",
        "string(3.5e-20)",
        "string(pi)",
        "string(i)",
        "string(1+2*i)",
        "str(x*y)",
    ]);
}

/// `expr` and `execute` (probe1, probe5, probe11).
#[test]
fn expressions_of_strings() {
    check(&[
        "expr(\"123\")",
        "expr(\"1+1\")",
        "expr(sin,x)",
        "expr(\"\")",
        "expr(\"   \")",
        "expr(\"x^2+2*x+1\")",
        "expr(\"[1,2,3]\")",
        "expr(\"\\\"ab\\\"\")",
        "expr(f,x,y)",
        "expr(string(x^3-1))",
        "execute(\"2+3\")",
        "expr(\"1.5\")",
    ]);
    // giac's product reads back multiplied out; ours is `ifactor`'s held
    // product, compared by its value (as tests/number_theory.rs does).
    let Some(giac) = giac() else {
        return;
    };
    let input = "expr(\"ifactor(54)\")";
    let held = ours(input).ok();
    let value = held.as_ref().and_then(|h| match h.args() {
        [inner] if h.head() == Some(&Expr::from(Symbol::HOLD)) => {
            cas_core::Context::new().eval(inner).ok()
        }
        _ => None,
    });
    let theirs = giac.value(input).ok();
    assert!(
        value.is_some() && value == theirs,
        "{input}: ours {held:?}, giac's {theirs:?}"
    );
}

/// `cat`, and non-strings in `join` and `concat` (probe1, probe5,
/// probe11).
#[test]
fn concatenations() {
    check_text(&[
        "cat(\"aaa\",c,12*3)",
        "cat()",
        "cat(\"a\")",
        "cat(\"a\",\"b\")",
        "cat(1/2,x^2,[1,2])",
        "cat(1.5,\"x\")",
        "cat(\"s\",-x)",
        "join(\"\",[1/2,[1,2]])",
    ]);
    check(&[
        "join(\",\",[1,2])",
        "concat(\"a\",1)",
        "concat(1,\"a\")",
        "concat(\"a\",[1,2])",
    ]);
}

/// Departures: giac's answer, and this CAS's value (`None`: the call
/// stays).
#[test]
fn departures() {
    let cases: &[(&str, &str, Option<&str>)] = &[
        // Our grammar has no octal integers.
        ("expr(\"0123\")", "83", Some("123")),
        // Statements are not read: they stay, with the reader's warning.
        ("expr(\"a:=3\")", "3", None),
        ("expr(\"1;2\")", "1,2", None),
        (
            "expr(\"x+\")",
            "GIAC_ERROR: Syntax compatibility mode xcas. Parse error line 2 column 2 at end of \
             input Error: Bad Argument Value",
            None,
        ),
        // Our input form puts the number first and spaces sums: "1 + x",
        // giac's "x+1".
        ("size(string(x+1))", "3", Some("5")),
        ("cat(\"a\",x+1)", "\"ax+1\"", Some("\"a1 + x\"")),
        ("concat(\"a\",x+1)", "\"ax+1\"", Some("\"a1 + x\"")),
        (
            "join(\",\",[1,x+1,\"z\"])",
            "\"1,x+1,z\"",
            Some("\"1,1 + x,z\""),
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
            None => {
                found
                    .as_ref()
                    .ok()
                    .and_then(|e| e.head().and_then(Expr::as_symbol))
                    == Some(Symbol::TO_EXPRESSION)
            }
        };
        if !ok {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    if let Some(giac) = giac() {
        giac.eval("purge(a)").ok();
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
