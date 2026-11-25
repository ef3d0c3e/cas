//! Xcas's `angle_radian` against giac (CAS-617): in degrees and grads, the
//! trigonometric functions of a number and the inverse ones in that unit,
//! a call still held and other arguments kept. giac's instance is shared by
//! the tests of this file, so each input sets and resets the setting within
//! one call (C0011); ours runs in a session, statement by statement, and
//! its values, read back as giac's are in radians, are compared with
//! giac's sequence in order.

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
fn degrees_and_grads_as_giac() {
    let Some(giac) = giac() else {
        return;
    };
    check(
        giac,
        &[
            "angle_radian:=0;sin(30);cos(60);asin(1/2);atan(1);sin(x);angle_radian:=1;sin(30)",
            "angle_radian:=0;sin(30.0);angle_radian:=1",
            "angle_radian:=0;sin(pi);angle_radian:=1",
            "angle_radian:=0;asin(x);angle_radian:=1",
            "angle_radian:=0;sin(x+30);angle_radian:=1",
            "angle_radian:=0;tan(45);sec(60);angle_radian:=1",
            "angle_radian:=0;acos(0);atan(x);asin(0.5);angle_radian:=1",
            "angle_radian:=0;cos(1);angle_radian:=1",
            "angle_radian:=0;acot(0);sin(390);cos(1.0);angle_radian:=1",
            // The hyperbolic functions and exp are the same in every mode.
            "angle_radian:=0;sinh(1);exp(i*30);angle_radian:=1",
            "angle_radian:=2;sin(100);asin(1);angle_radian:=1",
        ],
    );
}

/// Departure: giac writes `sin(15)` in degrees as `(sqrt(6)-sqrt(2))/4`;
/// ours has no closed form of `sin(pi/12)`, so the call stays.
#[test]
fn sin_15_degrees_stays() {
    let Some(giac) = giac() else {
        return;
    };
    let input = "angle_radian:=0;sin(15);angle_radian:=1";
    assert_eq!(giac.eval(input).unwrap(), "0,(sqrt(6)-sqrt(2))/4,1");
    let mine = session_values(input);
    assert_eq!(mine.get(1), Some(&ours("sin(15)").unwrap()), "{mine:?}");
}
