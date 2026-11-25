//! Tests of `ToString`, `ToExpression`, `Cat` (CAS-614), `LaTeX`
//! (CAS-615) and `ExportMathML` (CAS-616), and of the
//! non-strings in `StringJoin` and `Join`: what each does with no language
//! hook, and the logic (caps, warnings) through a stub language that prints
//! and reads full form. cas-lang's session tests run them with the
//! language itself.

use std::sync::Arc;

use crate::context::Limits;
use crate::language::{Form, Language};
use crate::trace::{EvalOptions, Warning};
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Full form both ways, refusing text past `max_len` bytes.
struct Stub;

impl Language for Stub {
    fn print(&self, e: &Expr, _: Form, _: &Context, max_len: usize) -> Result<String, String> {
        let text = e.full_form().to_string();
        if text.len() > max_len {
            return Err("too long".to_owned());
        }
        Ok(text)
    }

    fn read(&self, text: &str, _: &Context) -> Result<Expr, String> {
        Expr::parse_full_form(text).map_err(|e| e.to_string())
    }
}

/// Full form after the name of the form asked for: `Mathml:f[x]`.
struct Forms;

impl Language for Forms {
    fn print(&self, e: &Expr, form: Form, _: &Context, _: usize) -> Result<String, String> {
        Ok(format!("{form:?}:{}", e.full_form()))
    }

    fn read(&self, text: &str, _: &Context) -> Result<Expr, String> {
        Expr::parse_full_form(text).map_err(|e| e.to_string())
    }
}

fn stub() -> Context {
    let mut ctx = Context::new();
    ctx.set_language(Some(Arc::new(Stub)));
    ctx
}

/// `(input, expected)` in full form, both evaluated in `ctx`: the same
/// expression.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
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

/// Each input keeps its head in `ctx`, with the warnings `warned` says.
fn stays(ctx: &Context, inputs: &[&str], warned: bool) {
    for input in inputs {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert_eq!(
            !report.warnings.is_empty(),
            warned,
            "{input}: {:?}",
            report.warnings
        );
    }
}

/// The one warning of `input` in `ctx`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head(), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

#[test]
fn without_a_hook() {
    let ctx = Context::new();
    stays(
        &ctx,
        &[
            "ToString[x]",
            "ToString[\"a\"]",
            "LaTeX[x]",
            "ExportMathML[x]",
            "ExportMathML[x, display]",
            "ToExpression[\"1\"]",
            "Cat[\"a\", x]",
            "Cat[1]",
            "StringJoin[\",\", List[1, 2]]",
            "Join[\"a\", 1]",
            "Join[1, \"a\"]",
        ],
        false,
    );
    check(
        &ctx,
        &[
            ("Cat[]", "\"\""),
            ("Cat[\"a\"]", "\"a\""),
            ("Cat[\"a\", \"b\", \"\"]", "\"ab\""),
            ("StringJoin[\",\", List[\"a\", \"b\"]]", "\"a,b\""),
            ("Join[\"a\", \"b\"]", "\"ab\""),
            ("Join[\"a\", List[1, 2]]", "List[\"a\", 1, 2]"),
            ("ToExpression[\"\"]", "Indeterminate"),
            ("ToExpression[\" \\n\"]", "Indeterminate"),
            ("ToExpression[Sin, 0]", "0"),
            ("ToExpression[Function[List[y], Plus[y, 1]], 2]", "3"),
        ],
    );
}

#[test]
fn shapes_that_stay() {
    let ctx = stub();
    stays(
        &ctx,
        &[
            "ToString[]",
            "ToString[x, y]",
            "LaTeX[]",
            "LaTeX[x, \"file\"]",
            "ExportMathML[]",
            "ExportMathML[x, y]",
            "ExportMathML[x, \"display\"]",
            "ExportMathML[x, display, y]",
            "ToExpression[]",
            "ToExpression[x]",
            "ToExpression[\"1\", 2]",
            "ToExpression[f, x, y]",
            "ToExpression[Sin]",
            "StringJoin[x, List[1]]",
            "StringJoin[\",\", x]",
            "Join[x, 1]",
        ],
        false,
    );
}

#[test]
fn through_the_hook() {
    let ctx = stub();
    check(
        &ctx,
        &[
            ("ToString[f[x]]", "\"f[x]\""),
            ("LaTeX[f[x]]", "\"f[x]\""),
            ("ToString[\"a\"]", "\"\\\"a\\\"\""),
            ("ToExpression[\"Plus[1, 1]\"]", "2"),
            ("ToExpression[ToString[f[x, 2]]]", "f[x, 2]"),
            ("Cat[\"aaa\", c, Times[12, 3]]", "\"aaac36\""),
            ("Cat[f[x], \"b\"]", "\"f[x]b\""),
            ("StringJoin[\",\", List[1, x, \"z\"]]", "\"1,x,z\""),
            ("Join[\"a\", 1]", "\"a1\""),
            ("Join[1, \"a\"]", "\"1a\""),
        ],
    );
    // The reader's message as the warning.
    assert!(!warning(&ctx, "ToExpression[\"Plus[\"]").is_empty());
}

/// `ExportMathML`'s second argument picks the form: none the default,
/// `display` presentation, `content` (or `Content`, its builtin) content.
#[test]
fn export_mathml_forms() {
    let mut ctx = Context::new();
    ctx.set_language(Some(Arc::new(Forms)));
    check(
        &ctx,
        &[
            ("ExportMathML[f[x]]", "\"Mathml:f[x]\""),
            ("ExportMathML[f[x], display]", "\"MathmlDisplay:f[x]\""),
            ("ExportMathML[f[x], content]", "\"MathmlContent:f[x]\""),
            ("ExportMathML[f[x], Content]", "\"MathmlContent:f[x]\""),
            ("LaTeX[f[x]]", "\"Latex:f[x]\""),
        ],
    );
}

/// A symbol with a value or rules of its own is applied; one without stays.
#[test]
fn applies_user_functions() {
    let mut ctx = stub();
    ctx.define(&read("g[2]"), &read("5")).unwrap();
    ctx.set(crate::Symbol::new("h"), read("Sin")).unwrap();
    check(
        &ctx,
        &[("ToExpression[g, 2]", "5"), ("ToExpression[h, 0]", "0")],
    );
    stays(&ctx, &["ToExpression[k, 2]"], false);
}

/// Text past the `terms` limit stays: a print, a concatenation, a string to
/// read (with a warning).
#[test]
fn the_terms_limit() {
    let mut ctx = stub();
    let limits = ctx.limits();
    ctx.set_limits(Limits { terms: 4, ..limits });
    check(
        &ctx,
        &[
            ("ToString[f[x]]", "\"f[x]\""),
            ("Cat[\"ab\", \"cd\"]", "\"abcd\""),
            ("StringJoin[\"\", List[\"ab\", g]]", "\"abg\""),
        ],
    );
    assert_eq!(warning(&ctx, "ToString[f[xy]]"), "too long");
    stays(
        &ctx,
        &[
            "Cat[\"abc\", \"de\"]",
            "Cat[\"ab\", f[x]]",
            "StringJoin[\",\", List[\"ab\", \"cd\"]]",
        ],
        false,
    );
    assert_eq!(
        warning(&ctx, "ToExpression[\"12345\"]"),
        "the text is longer than the terms limit"
    );
}
