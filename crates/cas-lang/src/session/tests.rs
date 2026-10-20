//! Tests of the session; scripted sessions are in `tests/sessions`.

use std::time::{Duration, Instant};

use super::*;

fn results(outputs: &[Output]) -> Vec<Option<&str>> {
    outputs
        .iter()
        .map(|o| o.result.as_ref().map(Shown::text))
        .collect()
}

/// An evaluation stopped from another thread ends its statement with an
/// error, and the session goes on with the next one: here a branching
/// recursion with 2^30 distinct calls.
#[test]
fn an_interrupted_statement_ends_and_the_next_one_runs() {
    let mut session = Session::new();
    let interrupt = session.interrupt();
    let text = "g(0, x) := x\ng(n, x) := g(n - 1, 2x) + g(n - 1, 2x + 1)\ng(30, 0)\ng(3, 0)";
    let started = Instant::now();
    let outputs = std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(Duration::from_millis(200));
            interrupt.interrupt();
        });
        session.run(text)
    });
    assert!(started.elapsed() < Duration::from_secs(30));
    assert_eq!(results(&outputs), [None, None, None, Some("28")]);
    let stopped = &outputs[2];
    let start = text.find("g(30, 0)").unwrap();
    assert_eq!(
        stopped.messages,
        [Message {
            level: Level::Error,
            range: start..start + 8,
            text: "interrupted".to_owned(),
            fix: None,
        }]
    );
    assert_eq!(stopped.range, start..start + 8);
}

/// Each statement gives an output, with the diagnostics inside its text;
/// statements that do not lower give their diagnostics only.
#[test]
fn outputs_follow_the_statements() {
    let mut session = Session::new();
    let text = "1 + 1; 2 + 2\nsin x\nx(x + 1)";
    let outputs = session.run(text);
    let texts: Vec<&str> = outputs.iter().map(|o| &text[o.range.clone()]).collect();
    assert_eq!(texts, ["1 + 1", "2 + 2", "sin x", "x(x + 1)"]);
    assert_eq!(
        results(&outputs),
        [None, Some("4"), None, Some("x*(1 + x)")]
    );
    assert_eq!(outputs[2].messages.len(), 1);
    assert_eq!(outputs[2].messages[0].level, Level::Error);
    assert!(outputs[2].messages[0].fix.is_some());
    assert_eq!(outputs[3].messages[0].level, Level::Note);
    assert!(session.run("").is_empty());
}

/// What a statement defines or assumes stays for the next `run`.
#[test]
fn the_context_keeps_what_statements_did() {
    let mut session = Session::new();
    session.run("f(x) := 2x\nassume(a > 0)");
    assert_eq!(session.context().assumptions().len(), 1);
    assert_eq!(
        results(&session.run("f(a)\nabs(a)")),
        [Some("2*a"), Some("a")]
    );
    session.run("purge(a, f)");
    assert!(session.context().assumptions().is_empty());
    assert_eq!(results(&session.run("f(a)")), [Some("f(a)")]);
}

/// A locale's names read and print.
#[test]
fn a_locale_names_input_and_output() {
    let names = Names::with_locale("[builtins]\nLength: taille\nSin: sinus\n").unwrap();
    let mut session = Session::with_names(names);
    let outputs = session.run("taille([1, 2, 3])\nsinus(x)\nsin(y)");
    assert_eq!(
        results(&outputs),
        [Some("3"), Some("sinus(x)"), Some("sinus(y)")]
    );
    let spelling = session
        .names()
        .spelling(&crate::Target::Builtin(Symbol::LENGTH));
    assert_eq!(spelling, Some("taille"));
}

/// A result the language cannot write shows in the evaluator's notation,
/// with a note.
#[test]
fn results_the_language_cannot_write() {
    let mut session = Session::new();
    let outputs = session.run("subst(f(2), x = 3)");
    assert_eq!(results(&outputs), [Some("f(2)")]);
    let mut session = Session::new();
    session
        .context
        .set(
            Symbol::new("r"),
            Expr::parse_full_form("Rule[a, b]").unwrap(),
        )
        .unwrap();
    let outputs = session.run("r");
    assert_eq!(
        outputs[0].result,
        Some(Shown::Evaluator("Rule[a, b]".to_owned()))
    );
    assert_eq!(outputs[0].messages[0].level, Level::Note);
    assert_eq!(
        outputs[0].messages[0].text,
        "shown in the evaluator's notation: the language has no syntax for Rule with 2 arguments"
    );
}
