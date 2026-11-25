//! Xcas's `epsilon` against giac (CAS-617): the default bound of
//! `epsilon2zero` and of `trim`. giac's instance is shared by the tests of
//! this file, so each input sets and resets the setting within one call
//! (C0011); ours runs in a session, statement by statement, and its values
//! are compared with giac's sequence in order.

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

/// Whether two values are one: equal, their difference expanded and
/// together 0 (C0032), decimals within giac's printed digits (C0006), or
/// lists of the same values in any order (C0023).
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let is_list = |e: &Expr| e.head() == Some(&Expr::from(Symbol::LIST));
    if is_list(a) && is_list(b) {
        let mut rest: Vec<&Expr> = b.args().iter().collect();
        return a.args().len() == rest.len()
            && a.args().iter().all(|x| {
                let found = rest.iter().position(|y| same(ctx, x, y));
                found.map(|i| rest.remove(i)).is_some()
            });
    }
    let difference = a.clone() - b.clone();
    let normal = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [difference.clone()])],
    );
    if ctx.eval(&normal).ok() == Some(Expr::from(0)) {
        return true;
    }
    let tolerance = ours("1e-10").unwrap();
    let bound = tolerance * Expr::call(Symbol::ABS, [a.clone()]);
    let size = Expr::call(Symbol::ABS, [difference]);
    ctx.eval(&Expr::call(Symbol::LESS, [size, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// Inputs whose values are giac's, one by one, each run whole by giac and
/// in one session here.
fn check(giac: &Giac, inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let answer = giac
            .eval(input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let theirs: Vec<Expr> = values(&answer)
            .into_iter()
            // giac's `list[…]` is a list, as `Giac::value` reads it.
            .map(|v| ours(&v.replace("list[", "[")).unwrap_or_else(|e| panic!("{input}: {e:?}")))
            .collect();
        let mine = session_values(input);
        if theirs.len() != mine.len() || theirs.iter().zip(&mine).any(|(a, b)| !same(&ctx, a, b)) {
            wrong.push(format!("{input}: giac's {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn the_default_bound_as_giac() {
    let Some(giac) = giac() else {
        return;
    };
    check(
        giac,
        &[
            "epsilon:=1e-3;epsilon2zero(0.0001+x);epsilon:=1e-12",
            "epsilon:=1e-3;epsilon2zero(0.01+x);epsilon:=1e-12",
            "epsilon:=1e-3;trim([0.0001,1,0.0001,0.01]);epsilon:=1e-12",
        ],
    );
}
