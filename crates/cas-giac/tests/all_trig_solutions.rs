//! `all_trig_solutions` against giac (CAS-545): with the flag on, giac's
//! general solutions of trigonometric equations, read back, are this CAS's
//! in giac's order, element by element equal once expanded (giac writes
//! `nπ + π/2` as `-1/2*(-2*n*pi-pi)`). giac's parameters `n_k`, counted
//! across its session, read as `n_0`, `n_1`, … in their order. giac's input
//! sets the flag and resets it within one call, as its instance is shared by
//! the tests of this file; ours runs in a session with the flag set.

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

/// The parameters `n_5`, `n_7`, … of `answer` as `n_0`, `n_1`, … in the
/// order of their numbers.
fn renamed(answer: &str) -> String {
    let bytes = answer.as_bytes();
    let mut found: Vec<(usize, usize, u64)> = Vec::new();
    let mut i = 0;
    while i + 2 < bytes.len() {
        if bytes[i] == b'n' && bytes[i + 1] == b'_' {
            let digits = bytes[i + 2..]
                .iter()
                .take_while(|b| b.is_ascii_digit())
                .count();
            if digits > 0 {
                let end = i + 2 + digits;
                found.push((i, end, answer[i + 2..end].parse().expect("digits")));
                i = end;
                continue;
            }
        }
        i += 1;
    }
    let mut numbers: Vec<u64> = found.iter().map(|&(_, _, k)| k).collect();
    numbers.sort_unstable();
    numbers.dedup();
    let mut out = String::with_capacity(answer.len());
    let mut last = 0;
    for (start, end, k) in found {
        let index = numbers.binary_search(&k).expect("collected");
        out.push_str(&answer[last..start]);
        out.push_str("n_");
        out.push_str(&index.to_string());
        last = end;
    }
    out.push_str(&answer[last..]);
    out
}

/// The last of the values giac prints for a sequence of statements.
fn last(answer: &str) -> &str {
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => start = i + 1,
            _ => {}
        }
    }
    &answer[start..]
}

/// giac's general solutions of `solve(equation,x)`, in the language.
fn theirs(giac: &Giac, equation: &str) -> Expr {
    let input = format!("all_trig_solutions(1); a:=solve({equation},x); all_trig_solutions(0); a");
    let answer = giac
        .eval(&input)
        .unwrap_or_else(|e| panic!("{input}: {e:?}"));
    let text = renamed(&last(&answer).replace("list[", "["));
    ours(&text).unwrap_or_else(|e| panic!("{input}: {text}: {e:?}"))
}

/// This CAS's `solve(equation,x)` with `all_trig_solutions` on.
fn with_flag(equation: &str) -> Expr {
    let mut session = Session::new();
    let input = format!("all_trig_solutions(1); solve({equation},x)");
    let outputs = session.run(&input);
    let output = outputs.last().expect("a statement");
    let result = output
        .result
        .as_ref()
        .unwrap_or_else(|| panic!("{input}: {:?}", output.messages));
    ours(result.text()).unwrap_or_else(|e| panic!("{input}: {e:?}"))
}

/// Whether two lists are one, element by element, once their differences
/// are expanded.
fn same(a: &Expr, b: &Expr) -> bool {
    let ctx = Context::new();
    a.args().len() == b.args().len()
        && a.args().iter().zip(b.args()).all(|(a, b)| {
            let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
            ctx.eval(&difference).ok() == Some(Expr::from(0))
        })
}

#[test]
fn general_solutions() {
    let Some(giac) = giac() else {
        return;
    };
    for equation in [
        "sin(x)=0",
        "sin(x)=1/2",
        "sin(x)=-1/2",
        "sin(x)=1",
        "cos(x)=1/2",
        "cos(x)=1",
        "cos(x)=0",
        "cos(x)=-1",
        "tan(x)=1",
        "tan(x)=0",
        "sin(2*x)=0",
        "cos(2*x)=1/2",
        "sin(x)^2=1/4",
        "exp(x)=2",
    ] {
        let (theirs, ours) = (theirs(giac, equation), with_flag(equation));
        assert!(
            same(&theirs, &ours),
            "solve({equation},x): giac's {theirs:?}, ours {ours:?}"
        );
    }
}

/// The flag off, as by default, the principal solutions as before.
#[test]
fn principal_by_default() {
    let mut session = Session::new();
    let outputs = session.run("all_trig_solutions(1); all_trig_solutions(0); solve(sin(x)=0,x)");
    let result = outputs.last().and_then(|o| o.result.as_ref());
    assert_eq!(result.map(cas_lang::Shown::text), Some("[0, pi]"));
}
