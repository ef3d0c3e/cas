//! Xcas's `complex_variables` against giac (CAS-617): with it on, symbols
//! without facts are complex, this CAS's `Domain::Complex`. giac's instance
//! is shared by the tests of this file, so each input sets and resets the
//! setting within one call (C0011); ours runs in a session, statement by
//! statement. As reading a value back in a new, real context would make
//! `conj(x)` `x`, the values are read unevaluated and compared as they are
//! or evaluated in a complex context. The departures record giac's text
//! (C0007, C0015): the complex domain changes more than giac's flag does.

use cas_core::{Context, Expr};
use cas_giac::Giac;
use cas_lang::{Names, Session, lower, parse};

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

/// `text`, one statement, parsed and lowered with the English names, not
/// evaluated.
fn read(text: &str) -> Expr {
    let parsed = parse(text);
    assert!(parsed.errors().is_empty(), "{text}: {:?}", parsed.errors());
    let lowered = lower(&parsed, Names::english());
    assert!(!lowered.has_errors(), "{text}: {:?}", lowered.diagnostics);
    let [statement] = &lowered.statements[..] else {
        panic!("{text}: not one statement");
    };
    statement.expr.clone()
}

/// This CAS's values of the statements of `input` that give one, giac's `;`
/// read as separate statements, run in one session, read unevaluated.
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
            Some(read(result.text()))
        })
        .collect()
}

/// Whether two values are one: equal as read, or evaluated in `complex`.
fn same(complex: &Context, a: &Expr, b: &Expr) -> bool {
    a == b
        || complex
            .eval(a)
            .is_ok_and(|a| complex.eval(b).ok() == Some(a))
}

/// A context whose symbols are complex.
fn complex() -> Context {
    let mut complex = Context::new();
    complex.set_complex_variables(true);
    complex
}

#[test]
fn complex_symbols_as_giac() {
    let Some(giac) = giac() else {
        return;
    };
    let complex = complex();
    let mut wrong = Vec::new();
    for input in [
        "complex_variables:=1;conj(x);re(x);complex_variables:=0;conj(x)",
        "complex_variables:=1;abs(x);im(x);complex_variables:=0;abs(x);im(x);re(x)",
        "complex_variables:=1;factor(x^2+1);complex_variables:=0",
    ] {
        let answer = giac
            .eval(input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let theirs: Vec<Expr> = values(&answer).into_iter().map(read).collect();
        let mine = session_values(input);
        if theirs.len() != mine.len()
            || theirs.iter().zip(&mine).any(|(a, b)| !same(&complex, a, b))
        {
            wrong.push(format!("{input}: giac's {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures: `(giac's input, its answer, our input, our last value
/// decided)`, ours run in a session, decided in a complex context.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let complex = complex();
    let cases = [
        // giac's abs(x)^2 is x^2 whatever the flag; ours, of a complex x,
        // stays.
        (
            "complex_variables:=1;abs(x)^2;complex_variables:=0",
            "1,x^2,0",
            "complex_variables:=1;abs(x)^2",
            "abs(x)^2",
        ),
        // giac's solve stays real under the flag; ours solves in the
        // complex domain.
        (
            "complex_variables:=1;solve(x^2+1=0,x);complex_variables:=0",
            "1,[],0",
            "complex_variables:=1;solve(x^2+1=0,x)",
            "[i,-i]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, mine, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let found = session_values(mine).pop().expect("a value");
        if !same(&complex, &found, &read(decided)) {
            wrong.push(format!("{mine}: ours is {found:?}, not {decided}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
