//! `increasing_power` against giac (CAS-618): with the flag off, a sum of
//! powers of one symbol prints by decreasing degree in `string`, `latex`,
//! `export_mathml` and `mathml`, giac's (its evaluation sorts so, gen.cc
//! :10094) and ours (`cas_lang::power_order` before every printer). giac's
//! input sets the flag and back to giac's default, 0, within one call, as
//! its instance is shared by the tests of this file (C0011); ours runs in a
//! session with `increasing_power := 0`.
//!
//! Rows: `string`, giac's text is ours with its spaces removed (the term
//! order is the point), and both read back to one value; the markup forms,
//! giac's text, its quotes removed, is ours, whitespace collapsed (C0047),
//! in the default `export_mathml` form the annotation (each one's input
//! form) compared without spaces. Departures: a rational coefficient
//! (giac `1/2*x^2-x+1/2`, ours `x^2/2 - x + 1/2`, same order), and sums of
//! two symbols, which ours leaves in its canonical order (giac
//! `x^2+y^2+2*x*y`): same value.

use cas_giac::{Giac, ours};
use cas_lang::Session;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The last of the values giac prints for a sequence of statements; commas
/// inside strings are not separators.
fn last(answer: &str) -> &str {
    let (mut depth, mut quoted, mut start) = (0i32, false, 0);
    for (i, c) in answer.char_indices() {
        match c {
            '"' => quoted = !quoted,
            '(' | '[' if !quoted => depth += 1,
            ')' | ']' if !quoted => depth -= 1,
            ',' if depth == 0 && !quoted => start = i + 1,
            _ => {}
        }
    }
    &answer[start..]
}

/// The text of a string answer: giac's quotes removed and its doubled
/// inner quotes undone.
fn unquoted(answer: &str) -> Option<String> {
    let inner = answer.trim().strip_prefix('"')?.strip_suffix('"')?;
    Some(inner.replace("\"\"", "\""))
}

fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// giac's text of `call` with the flag as `flag` says, set back to giac's
/// default, 0, after it.
fn giac_text(giac: &Giac, call: &str, flag: u8) -> String {
    let input = format!("increasing_power:={flag};a:={call};increasing_power:=0;a");
    let answer = giac
        .eval(&input)
        .unwrap_or_else(|e| panic!("{input}: {e:?}"));
    unquoted(last(&answer)).unwrap_or_else(|| panic!("{input}: {answer}"))
}

/// Our text of `call`, a string, in a session with `increasing_power` as
/// `flag` says.
fn our_text(call: &str, flag: u8) -> String {
    let mut session = Session::new();
    let outputs = session.run(&format!("increasing_power := {flag}; {call}"));
    let output = outputs.last().expect("a statement");
    let result = output
        .result
        .as_ref()
        .unwrap_or_else(|| panic!("{call}: {:?}", output.messages));
    let value = ours(result.text()).unwrap_or_else(|e| panic!("{call}: {e:?}"));
    value
        .as_str()
        .unwrap_or_else(|| panic!("{call}: {value:?}"))
        .to_owned()
}

/// Whether two input forms read back to one value.
fn same_value(a: &str, b: &str) -> bool {
    matches!((ours(a), ours(b)), (Ok(a), Ok(b)) if a == b)
}

#[test]
fn string_by_decreasing_degree() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for e in [
        "expand((x+1)^3)",
        "expand(1-x^2)",
        "expand(2*x^2-x+5)",
        "expand(x^5-x^3+x)",
        "expand((t-2)^3)",
        "expand(3*(x+1))",
        "expand((x+i)^2)",
        "sin(expand((x+1)^2))",
        "[expand((x+1)^2),expand((t-2)^3)]",
        "expand((x+1)^3)=0",
    ] {
        let call = format!("string({e})");
        let theirs = giac_text(giac, &call, 0);
        let mine = our_text(&call, 0);
        if mine.replace(' ', "") != theirs || !same_value(&theirs, &mine) {
            wrong.push(format!("{call}: giac {theirs}, ours {mine}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn markup_by_decreasing_degree() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for call in [
        "latex(expand((x+1)^3))",
        "latex(expand(1-x^2))",
        "latex(expand(2*x^2-x+5))",
        "latex(expand(-x^5+x^3-1))",
        "latex(sin(expand((x+1)^2)))",
        "latex(2*x^2-3*x+1)",
        "latex(expand((x-1)^2/2))",
        "export_mathml(expand((x+1)^3),display)",
        "export_mathml(expand((x+1)^3),content)",
        "export_mathml(expand(1-x^2),display)",
        "mathml(expand((x+1)^3))",
        "mathml(expand(1-x^2))",
    ] {
        let theirs = giac_text(giac, call, 0);
        let mine = our_text(call, 0);
        if collapsed(&theirs) != collapsed(&mine) {
            wrong.push(format!("{call}: giac {theirs}, ours {mine}"));
        }
    }
    // The default form: the formula as giac's, the annotation giac's
    // string without spaces.
    let call = "export_mathml(expand((x+1)^3))";
    let (theirs, mine) = (giac_text(giac, call, 0), our_text(call, 0));
    let split = |t: &str| {
        let (head, rest) = t.split_once("<annotation encoding='giac'>")?;
        Some((head.to_owned(), rest.replace(' ', "")))
    };
    if split(&theirs).is_none() || split(&theirs) != split(&mine) {
        wrong.push(format!("{call}: giac {theirs}, ours {mine}"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures, same values: a rational coefficient spelled otherwise, and
/// sums of two symbols, which ours leaves as with the flag on.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let call = "string(expand((x-1)^2/2))";
    let theirs = giac_text(giac, call, 0);
    assert_eq!(theirs, "1/2*x^2-x+1/2");
    let mine = our_text(call, 0);
    assert_eq!(mine, "x^2/2 - x + 1/2");
    assert!(same_value(&theirs, &mine), "{call}");
    for (call, giacs) in [
        ("string(expand((x+y)^2))", "x^2+y^2+2*x*y"),
        ("string(x^2+y+1)", "x^2+y+1"),
    ] {
        let theirs = giac_text(giac, call, 0);
        assert_eq!(theirs, giacs);
        let mine = our_text(call, 0);
        assert_eq!(mine, our_text(call, 1), "{call}");
        assert!(
            same_value(&theirs, &mine),
            "{call}: giac {theirs}, ours {mine}"
        );
    }
}

/// With the flag back on, ours by increasing degree as giac's (probe7).
#[test]
fn increasing_degree_with_the_flag_on() {
    let Some(giac) = giac() else {
        return;
    };
    for call in ["latex(expand((x+1)^3))", "string(expand((x+1)^3))"] {
        let theirs = giac_text(giac, call, 1);
        let mine = our_text(call, 1);
        assert_eq!(mine.replace(' ', ""), theirs.replace(' ', ""), "{call}");
    }
}
