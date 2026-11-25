//! Xcas's settings of phase 1.6 against giac (CAS-617, CAS-618): their
//! values read by name or by `name()`, set by `name(v)` or `name := v`, giac's
//! answers of `rur_certify`, `gbasis_max_pairs` and `gbasis_simult_primes`,
//! and `cas_setup()` element by element where ours is giac's. giac's instance
//! is shared by the tests of this file, so each input sets and resets a
//! setting within one call (C0011); ours runs in a session, statement by
//! statement, and its values are compared with giac's sequence in order.
//! The departures record giac's text (C0007, C0015): its bare names of
//! settings without a value, its evaluation levels and `regroup` by default,
//! its `approx_mode` that leaves exact results, and its mode switches that
//! are not builtins here (CAS-DOC-7's one grammar).

use cas_core::{Context, Expr, Symbol};
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

/// Inputs whose values are giac's, one by one.
#[test]
fn settings_as_giac_reads_and_sets_them() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "approx_mode",
        "complex_variables",
        "angle_radian",
        "epsilon",
        "epsilon:=1e-6;epsilon;epsilon:=1e-12",
        "approx_mode(1);approx_mode();approx_mode(0)",
        "complex_mode(1);complex_mode();complex_mode(0)",
        "keep_algext(1);keep_algext();keep_algext(0)",
        "increasing_power(1);increasing_power();increasing_power(0)",
        "angle_radian(0);angle_radian();angle_radian(1)",
        // giac's answers (solve.cc): the stored value, the old one.
        "rur_certify(1)",
        "gbasis_max_pairs(100);gbasis_max_pairs(32768)",
        "gbasis_simult_primes(4);gbasis_simult_primes(20)",
    ] {
        let answer = giac
            .eval(input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let theirs: Vec<Expr> = values(&answer)
            .into_iter()
            .map(|v| ours(v).unwrap_or_else(|e| panic!("{input}: {e:?}")))
            .collect();
        let mine = session_values(input);
        if theirs.len() != mine.len() || theirs.iter().zip(&mine).any(|(a, b)| !same(&ctx, a, b)) {
            wrong.push(format!("{input}: giac's {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `cas_setup()`: giac's list element by element, but its `Digits` (12,
/// ours the working precision's 15), its threads and recursion limit
/// (`[12,100,0,25]`, ours `[1, recursion, 0, 25]`) and its
/// `increasing_power` (0, ours 1): departures. `cas_setup` of fewer than 7
/// elements is giac's error; ours stays.
#[test]
fn cas_setup_lists_the_settings() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    // giac's threads (element 7's first) are the host's cores: 12 here, 4
    // on CI; so its text is not asserted.
    let answer = giac.eval("cas_setup()").unwrap();
    let theirs = ours(&answer).unwrap();
    let [mine] = &session_values("cas_setup()")[..] else {
        panic!("one value");
    };
    assert_eq!(theirs.args().len(), 12, "{answer}");
    assert_eq!(theirs.args().len(), mine.args().len(), "{mine:?}");
    for i in [0, 1, 2, 3, 4, 9, 10, 11] {
        assert_eq!(theirs.args()[i], mine.args()[i], "element {i}: {mine:?}");
    }
    let (a, b) = (&theirs.args()[5], &mine.args()[5]);
    assert!(
        a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(&ctx, a, b)),
        "element 5: {mine:?}"
    );
    for i in [6, 7, 8] {
        assert_ne!(theirs.args()[i], mine.args()[i], "element {i} agrees now");
    }
    assert_eq!(
        giac.eval("cas_setup([0,0,0])").unwrap(),
        "GIAC_ERROR: cas_setup([0,0,0]) \n Error: Invalid dimension"
    );
    assert_eq!(
        session_values("cas_setup([0,0,0])"),
        [ours("cas_setup([0,0,0])").unwrap()]
    );
}

/// Departures: `(giac's input, its answer, our input, our last value
/// decided)`, ours run in a session.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac's bare names of settings it reads by their call only; ours
        // read the value.
        ("complex_mode", "'complex_mode'", "complex_mode", "0"),
        ("keep_algext", "'keep_algext'", "keep_algext", "0"),
        (
            "increasing_power",
            "'increasing_power'",
            "increasing_power",
            "1",
        ),
        ("eval_level", "'eval_level'", "eval_level", "25"),
        ("autosimplify", "'autosimplify'", "autosimplify", "0"),
        // giac's library leaves exact results with approx_mode on; ours
        // shows decimals.
        (
            "approx_mode:=1;1/3+sqrt(2);approx_mode:=0",
            "1,1/3+sqrt(2),0",
            "approx_mode:=1;1/3+sqrt(2)",
            "evalf(1/3+sqrt(2))",
        ),
        // giac's autosimplify(factor) leaves x^2-1; ours factors.
        (
            "autosimplify(factor);x^2-1;autosimplify(regroup)",
            "1,x^2-1,1",
            "autosimplify(factor);x^2-1",
            "factor(x^2-1)",
        ),
        // giac evaluates a once at level 1; ours to a fixpoint.
        (
            "eval_level(1);a:=b;b:=2;a;eval_level(25);purge(a,b)",
            "1,b,2,b,25,b,2",
            "eval_level(1);a:=b;b:=2;a",
            "2",
        ),
        // giac's default command is regroup; ours none.
        ("autosimplify()", "'regroup'", "autosimplify()", "0"),
        // giac's mode switches; not builtins here (CAS-DOC-7's one
        // grammar): the calls stay.
        (
            "xcas_mode(0)",
            "\"Warning: some commands like subs might change arguments order\"",
            "xcas_mode(0)",
            "xcas_mode(0)",
        ),
        (
            "maple_mode(0)",
            "\"Warning: some commands like subs might change arguments order\"",
            "maple_mode(0)",
            "maple_mode(0)",
        ),
        (
            "python_compat(0)",
            "0",
            "python_compat(0)",
            "python_compat(0)",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, mine, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let found = session_values(mine).pop();
        let expected = session_values(decided).pop();
        if found != expected {
            wrong.push(format!("{mine}: ours is {found:?}, not {decided}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
