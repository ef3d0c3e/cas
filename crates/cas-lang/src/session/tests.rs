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

/// A loop that never ends, and rewrites nothing, stops at an interrupt too.
#[test]
fn an_endless_loop_is_interrupted() {
    let mut session = Session::new();
    let interrupt = session.interrupt();
    let outputs = std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(Duration::from_millis(100));
            interrupt.interrupt();
        });
        session.run("begin while true do end end\n1 + 1")
    });
    assert_eq!(results(&outputs), [None, Some("2")]);
    assert_eq!(outputs[0].messages[0].text, "interrupted");
}

/// Decimals are balls made when their statement runs, at the session's
/// working precision then (CAS-148).
#[test]
fn decimals_take_the_precision_of_their_statement() {
    let mut session = Session::new();
    assert_eq!(session.working_precision(), 53);
    let shown = |session: &mut Session, text: &str| -> Vec<String> {
        session
            .run(text)
            .iter()
            .filter_map(|o| o.result.as_ref().map(|r| r.text().to_owned()))
            .collect()
    };
    assert_eq!(shown(&mut session, "1/3.0"), ["0.333333333333333"]);
    assert_eq!(shown(&mut session, "g(t) := t + 0.1"), Vec::<String>::new());
    session.set_working_precision(100);
    assert_eq!(
        shown(&mut session, "1/3.0"),
        ["0.333333333333333333333333333333"]
    );
    // An inexact literal is read at 100 bits too.
    assert_eq!(
        shown(&mut session, "1/3.1"),
        ["0.322580645161290322580645161290"]
    );
    // A literal made before keeps its precision: `g` read 0.1 at 53 bits.
    assert_eq!(shown(&mut session, "g(1)"), ["1.1"]);
    // The 53-bit 0.1 of `g` against a 100-bit 1.1: no digit in common.
    let [difference] = shown(&mut session, "g(1) - 1.1").try_into().unwrap();
    assert!(difference.starts_with("ball(0.0, "), "{difference}");
}

/// Xcas's settings `all_trig_solutions` and `with_sqrt` (CAS-545): set by
/// assignment or giac's call, read back by name, 0 or 1 only.
#[test]
fn settings_are_set_both_ways() {
    let mut session = Session::new();
    assert_eq!(
        results(&session.run("all_trig_solutions\nwith_sqrt")),
        [Some("0"), Some("1")]
    );
    assert_eq!(
        results(&session.run("all_trig_solutions := 1\nsolve(sin(x) = 0, x)")),
        [Some("1"), Some("[n_0*pi]")]
    );
    assert!(session.context().all_trig_solutions());
    assert_eq!(
        results(&session.run("all_trig_solutions(0);\nsolve(sin(x) = 0, x)")),
        [None, Some("[0, pi]")]
    );
    assert_eq!(
        results(&session.run("with_sqrt(false)\nfactor(x^2 - 2)\nwith_sqrt")),
        [Some("false"), Some("-2 + x^2"), Some("0")]
    );
    assert!(!session.context().with_sqrt());
    let outputs = session.run("with_sqrt := 2");
    assert_eq!(results(&outputs), [None]);
    assert_eq!(
        outputs[0].messages[0].text,
        "with_sqrt takes 0 or 1, false or true"
    );
    assert!(!session.context().with_sqrt());
}

/// Xcas's settings of phase 1.6 (CAS-617, CAS-618): read by name or by
/// `name()`, set by assignment or giac's call, with giac's answers;
/// `autosimplify` and `approx_mode` act on the results; `cas_setup` sets
/// them in giac's order.
#[test]
fn phase_1_6_settings() {
    let mut session = Session::new();
    assert_eq!(
        results(&session.run(
            "approx_mode\ncomplex_mode\ncomplex_variables\nangle_radian\nkeep_algext\n\
             increasing_power\neval_level\nautosimplify"
        )),
        [
            Some("0"),
            Some("0"),
            Some("0"),
            Some("1"),
            Some("0"),
            Some("1"),
            Some("25"),
            Some("0")
        ]
    );
    assert_eq!(
        results(&session.run("approx_mode(1)\napprox_mode()\n1/3\napprox_mode := 0\n1/3")),
        [
            Some("1"),
            Some("1"),
            Some("0.333333333333333"),
            Some("0"),
            Some("1/3")
        ]
    );
    let expanded = results(&session.run("expand((x + 1)^2)"))[0].map(str::to_owned);
    let outputs = session.run("autosimplify(expand)\n(x + 1)^2\nautosimplify()\nautosimplify(0)");
    assert_eq!(
        results(&outputs),
        [Some("1"), expanded.as_deref(), Some("expand"), Some("1")]
    );
    assert!(session.context().autosimplify().is_none());
    assert_eq!(
        results(
            &session.run("complex_mode := 1\nkeep_algext(true)\neval_level(3)\nangle_radian(0)")
        ),
        [Some("1"), Some("true"), Some("3"), Some("0")]
    );
    let ctx = session.context();
    assert!(ctx.complex_mode() && ctx.keep_algext());
    assert_eq!((ctx.eval_level(), ctx.angle()), (3, Angle::Degree));
    session.run("epsilon := 1e-6");
    assert_eq!(session.context().epsilon().to_bits(), 1e-6_f64.to_bits());
    for (input, message) in [
        (
            "angle_radian := 3",
            "angle_radian takes 0, 1 or 2: degrees, radians or grads",
        ),
        ("epsilon(-1)", "epsilon takes a positive number"),
        ("eval_level := 0", "eval_level takes a positive integer"),
        (
            "approx_mode := 2",
            "approx_mode takes 0 or 1, false or true",
        ),
    ] {
        let outputs = session.run(input);
        assert_eq!(results(&outputs), [None], "{input}");
        assert_eq!(outputs[0].messages[0].text, message, "{input}");
    }
    assert_eq!(session.context().angle(), Angle::Degree);
    assert_eq!(
        results(&session.run(
            "rur_certify(1)\nrur_certify()\nrur_certify(-3)\ngbasis_max_pairs(100)\n\
             gbasis_max_pairs()\ngbasis_simult_primes(0)\ngbasis_simult_primes()"
        )),
        [
            Some("0"),
            Some("1"),
            Some("3"),
            Some("32768"),
            Some("1073741824"),
            Some("20"),
            Some("1")
        ]
    );
    session.run("cas_setup(1, 1, 0, 2, 0, [1e-3, 1e-15], 20, [1, 50, 0, 10], 0, 0, 1);");
    let ctx = session.context();
    assert!(ctx.approx_mode() && ctx.complex_variables() && !ctx.complex_mode());
    assert_eq!((ctx.angle(), ctx.eval_level()), (Angle::Grad, 10));
    assert_eq!(ctx.epsilon().to_bits(), 1e-3_f64.to_bits());
    assert_eq!(ctx.working_precision(), cas_core::digits_bits(20));
    assert!(!ctx.increasing_power() && !ctx.with_sqrt() && ctx.all_trig_solutions());
    assert_eq!(
        results(&session.run("approx_mode(0)\ncas_setup([0, 0, 0])")),
        [Some("0"), Some("cas_setup([0, 0, 0])")]
    );
    let outputs = session.run("cas_setup()");
    let setup = results(&outputs)[0].unwrap_or_default();
    assert!(setup.starts_with("[0, 1, 0, 2, 0, ["), "{setup}");
}

/// `string`, `expr` and `cat` through the session's language hook
/// (CAS-614): the input form, read back to the same value.
#[test]
fn expressions_as_text_and_back() {
    let mut session = Session::new();
    assert_eq!(
        results(&session.run("size(string(x + 1))\ncat(\"aaa\", c, 12*3)")),
        [Some("5"), Some("\"aaac36\"")]
    );
    for text in [
        "x^2 + 2*x + 1",
        "[1, x, [2]]",
        "3/4",
        "sin(x)/2",
        "\"ab\"",
        "f(x, -y)",
    ] {
        let outputs = session.run(&format!("expr(string({text}))\n{text}"));
        let found = results(&outputs);
        assert!(
            found[0].is_some() && found[0] == found[1],
            "{text}: {found:?}"
        );
    }
    // A statement is not read: it stays, with the reader's warning, and
    // assigns nothing.
    let outputs = session.run("expr(\"a:=3\")\na");
    assert_eq!(outputs[0].messages.len(), 1, "{:?}", outputs[0].messages);
    assert_eq!(outputs[0].messages[0].level, Level::Warning);
    assert!(results(&outputs)[0].is_some_and(|t| t.starts_with("expr(")));
    assert_eq!(results(&outputs)[1], Some("a"));
}
