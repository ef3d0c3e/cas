//! Xcas's session commands against giac (CAS-619, CAS-620): `VARS()`,
//! `rm_all_vars()`, `rm_a_z()`, `restart` and `nodisp`. giac's instance is
//! shared, and these commands change its whole state, so one test runs
//! giac's inputs in turn, each after giac's `restart` (C0011); ours runs
//! each in a new session, statement by statement, and its values are
//! compared with giac's sequence in order. The departures record giac's
//! text (C0007, C0015) and check ours alone (C0044): its history, which
//! `caseval` does not keep, `time`'s list, `monotonic`'s clock, its
//! version, and its session commands nested in an expression.

use cas_core::{Context, Expr, Symbol, View};
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

/// The values giac prints for a sequence of statements, in order.
fn values(answer: &str) -> Vec<&str> {
    let mut depth = 0i32;
    let mut start = 0;
    let mut values = Vec::new();
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => {
                values.push(&answer[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    values.push(&answer[start..]);
    values
}

/// This CAS's values of the statements of `input` that give one, giac's `;`
/// read as separate statements, run in one session.
fn session_values(input: &str) -> Vec<Expr> {
    let mut session = Session::new();
    session
        .run(&input.replace(';', "\n"))
        .iter()
        .filter_map(|output| {
            let Some(result) = &output.result else {
                assert!(output.messages.is_empty(), "{input}: {:?}", output.messages);
                return None;
            };
            Some(ours(result.text()).unwrap_or_else(|e| panic!("{input}: {e:?}")))
        })
        .collect()
}

/// Whether two values are one: equal, or decimals within giac's printed
/// digits (C0006).
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let tolerance = ours("1e-10").unwrap();
    let difference = Expr::call(Symbol::ABS, [a.clone() - b.clone()]);
    let bound = tolerance * Expr::call(Symbol::ABS, [a.clone()]);
    ctx.eval(&Expr::call(Symbol::LESS, [difference, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// A nonnegative decimal.
fn nonnegative_decimal(e: &Expr) -> bool {
    matches!(e.view(), View::Ball(b) if b.mid().to_f64() >= 0.0)
}

/// Inputs whose values are giac's, one by one; `true` where our
/// definitions show nothing and giac's show a function: the last values
/// only are compared.
#[test]
fn session_commands_as_giac() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, last_only) in [
        ("a:=1;b:=2;VARS()", false),
        ("b:=1;a:=2;VARS()", false),
        ("zz:=1;f(x):=x^2;a:=3;VARS()", true),
        ("a:=5;rm_all_vars();a", false),
        ("a:=5;ab:=6;rm_a_z();[a,ab]", false),
        ("restart", false),
        ("a:=5;restart()", false),
        ("angle_radian:=0;restart;angle_radian", false),
    ] {
        giac.restart();
        let answer = giac
            .eval(input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let mut texts = values(&answer);
        let mut mine = session_values(input);
        if last_only {
            texts = texts.split_off(texts.len().saturating_sub(1));
            mine = mine.split_off(mine.len().saturating_sub(1));
        }
        let theirs: Vec<Expr> = texts
            .into_iter()
            .map(|v| ours(v).unwrap_or_else(|e| panic!("{input}: {e:?}")))
            .collect();
        if theirs.len() != mine.len() || theirs.iter().zip(&mine).any(|(a, b)| !same(&ctx, a, b)) {
            wrong.push(format!("{input}: giac's {answer}, ours {mine:?}"));
        }
    }
    // A string, compared raw: giac's without its quotes.
    giac.restart();
    let answer = giac.eval("nodisp(1+1)").unwrap();
    let theirs = answer.trim_matches('"');
    let mine = session_values("nodisp(1+1)");
    if mine.len() != 1 || mine[0].as_str() != Some(theirs) {
        wrong.push(format!("nodisp(1+1): giac's {answer}, ours {mine:?}"));
    }
    giac.restart();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures: giac's text recorded, ours checked alone (C0044).
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    // giac's caseval keeps no history: `undef`. Ours answers from the
    // session's history.
    for (input, answer, decided) in [
        ("1+1;ans()", "2,undef", "2"),
        ("1+1;quest()", "2,undef", "\"1+1\""),
        // giac answers 1: its restart() inside an expression is a
        // sequence of nothing; ours stays, with a warning.
        ("1+restart()", "1", "1+restart()"),
    ] {
        giac.restart();
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let mut session = Session::new();
        let outputs = session.run(&input.replace(';', "\n"));
        let found = outputs
            .last()
            .and_then(|o| o.result.as_ref())
            .map(|r| ours(r.text()).unwrap_or_else(|e| panic!("{input}: {e:?}")));
        let expected = ours(decided).unwrap();
        if found.as_ref() != Some(&expected) {
            wrong.push(format!("{input}: ours is {found:?}, not {decided}"));
        }
    }
    giac.restart();
    // giac's time(1+1) is [real, cpu] (probe: [2.6e-08,2.5644606e-08]),
    // its monotonic() a clock in seconds (1154478.57914), its version()
    // "giac 2.0.0, (c) B. Parisse and R. De Graeve, …": ours one
    // nonnegative decimal, the seconds since the session began, and "cas "
    // with this crate's version.
    for input in ["time(1+1)", "monotonic()"] {
        let found = session_values(input).pop();
        if !found.as_ref().is_some_and(nonnegative_decimal) {
            wrong.push(format!("{input}: ours is {found:?}"));
        }
    }
    let version = session_values("version()").pop();
    if !version
        .as_ref()
        .and_then(Expr::as_str)
        .is_some_and(|v| v.starts_with("cas "))
    {
        wrong.push(format!("version(): ours is {version:?}"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
