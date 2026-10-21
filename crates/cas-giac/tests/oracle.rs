//! The oracle against giac itself: answers recorded in CAS-DOC-8, read back
//! in the language; agreement with this CAS where it computes already; the
//! fast entry point against the slow one; the skip path.

use cas_giac::{Giac, OracleError, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac 1.9's answers recorded in CAS-DOC-8 §3, §4, §7 and §8, as it
/// prints them; each reads back in the language and evaluates to the same
/// value as the recorded text.
#[test]
fn recorded_answers() {
    let Some(giac) = giac() else {
        return;
    };
    let recorded = [
        // §3: roots of rationals.
        ("gcd(x^2-2, x-sqrt(2))", "x-sqrt(2)"),
        ("normal((x^2-2)/(x-sqrt(2)))", "x+sqrt(2)"),
        ("normal(1/(sqrt(2)-1))", "sqrt(2)+1"),
        ("normal(1/(sqrt(2)+sqrt(3)))", "sqrt(3)-sqrt(2)"),
        ("normal(sqrt(6)-sqrt(2)*sqrt(3))", "0"),
        ("normal(1/(2^(1/3)-1))", "(2^(1/3)+1)*2^(1/3)+1"),
        // §4: factor.
        ("factor(x^2-2)", "(x-sqrt(2))*(x+sqrt(2))"),
        ("factor(x^2-3*x+1)", "(x+(-sqrt(5)-3)/2)*(x+(sqrt(5)-3)/2)"),
        ("factor(x^4+1)", "(x^2+(-sqrt(2))*x+1)*(x^2+sqrt(2)*x+1)"),
        ("factor(x^2+1)", "x^2+1"),
        ("factor(2*x^2-2)", "2*(x-1)*(x+1)"),
        ("factor(-x^2+1)", "-(x-1)*(x+1)"),
        ("factor(6*x+12)", "6*(x+2)"),
        ("factor(x^4-y^4)", "(x-y)*(x+y)*(x^2+y^2)"),
        // §7: expand.
        ("expand(sin((x+1)^2))", "sin(x^2+2*x+1)"),
        ("expand(exp(a+b))", "exp(a+b)"),
        // §8: normal.
        ("normal(1/x+1/(x+1))", "(2*x+1)/(x^2+x)"),
        ("normal((x^2-1)/(x^2+2*x+1))", "(x-1)/(x+1)"),
        ("normal((x^2-a^2)/(x-a))", "a+x"),
        ("normal(1/(1-x))", "-1/(x-1)"),
    ];
    for (input, answer) in recorded {
        assert_eq!(giac.eval(input).unwrap(), answer, "{input}");
        giac.answers(input, answer)
            .unwrap_or_else(|e| panic!("{e}"));
    }
}

/// Where this CAS computes already, its value and giac's agree.
#[test]
fn agreement() {
    let Some(giac) = giac() else {
        return;
    };
    for input in [
        "1/3+1/6",
        "(2/3)^3",
        "[10,20,30][1]",
        "diff(sin(x^2), x)",
        "diff(x^3+2*x, x)",
        "sqrt(8)",
        "exact(0.1)",
        "evalf(pi, 30)",
        // giac's default is 12 digits (Xcas's Digits), ours 15.
        "evalf(1/3, 15)",
        "2^100",
        "irem(17, 5)",
    ] {
        giac.agrees(input).unwrap_or_else(|e| panic!("{e}"));
    }
    // A difference shows both sides.
    let err = giac.answers("1/3+1/6", "1/3").unwrap_err();
    assert!(
        matches!(&err, OracleError::Differs { giac, .. } if giac == "1/2"),
        "{err}"
    );
    assert_eq!(
        err.to_string(),
        "1/3+1/6: giac gives 1/2 (1/2 here), this CAS 1/3"
    );
}

/// `nws_caseval`, formatting undone, against `caseval`.
#[test]
fn fast_answers_are_the_slow_ones() {
    let Some(giac) = giac() else {
        return;
    };
    for input in [
        "factor(x^2-2)",
        "normal(1/(sqrt(2)+sqrt(3)))",
        "[1/2,sqrt(2)]",
        "[[1,2],[3,4]]",
        "[[1/2,2],[3,4]]",
        "1/3",
        "-sqrt(2)",
        "2^100",
        "1/2^100",
        "x^2+2*x+1",
        "sin(1)",
    ] {
        assert_eq!(
            giac.eval_fast(input).unwrap(),
            giac.eval(input).unwrap(),
            "{input}"
        );
    }
}

/// giac's state persists until `restart`. A name no other test uses keeps
/// the parallel tests apart.
#[test]
fn restart_forgets() {
    let Some(giac) = giac() else {
        return;
    };
    assert_eq!(giac.eval("oracle_restart:=3").unwrap(), "3");
    assert_eq!(giac.eval("oracle_restart+1").unwrap(), "4");
    giac.restart();
    assert_eq!(giac.eval("oracle_restart+1").unwrap(), "oracle_restart+1");
}

/// Without the library, loading says why, and checks skip.
#[test]
fn a_missing_library_skips() {
    let err = Giac::load("/nonexistent/libgiac.so.0").unwrap_err();
    assert!(err.contains("/nonexistent/libgiac.so.0"), "{err}");
    assert!(Giac::load("lib\0giac").is_err());
}

#[test]
fn ours_reads_the_language() {
    assert_eq!(ours("1/3+1/6"), ours("1/2"));
    // Decimals are balls in held parts too, as the session reads them.
    let held = ours("hold(0.5)").unwrap();
    assert!(
        matches!(held.args()[0].view(), cas_core::View::Ball(_)),
        "{held:?}"
    );
    assert!(matches!(ours("sin x"), Err(OracleError::Unreadable(..))));
    assert!(matches!(ours("1; 2"), Err(OracleError::Unreadable(..))));
}
