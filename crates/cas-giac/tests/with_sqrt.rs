//! `with_sqrt` against giac (CAS-545): with the flag off, giac's `factor`,
//! read back, is this CAS's, over ℚ without square roots. giac refuses
//! `with_sqrt:=0` (a reserved word), so its input calls `with_sqrt(0)` and
//! `with_sqrt(1)` within one call, as its instance is shared by the tests
//! of this file; ours runs in a session with the flag off.

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

/// This CAS's value of the last statement of `input`, run in a session.
fn session_value(input: &str) -> Expr {
    let mut session = Session::new();
    let outputs = session.run(input);
    let output = outputs.last().expect("a statement");
    let result = output
        .result
        .as_ref()
        .unwrap_or_else(|| panic!("{input}: {:?}", output.messages));
    ours(result.text()).unwrap_or_else(|e| panic!("{input}: {e:?}"))
}

#[test]
fn factor_without_square_roots() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for p in [
        "x^2-2",
        "3*x^2-6",
        "(x^2-2)*(x+1)",
        "x^4-4",
        "x^2-3*x+1",
        "x^2-4",
        "x^2+1",
    ] {
        let input = format!("with_sqrt(0); a:=factor({p}); with_sqrt(1); a");
        let answer = giac
            .eval(&input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let theirs = ours(last(&answer)).unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let ours = session_value(&format!("with_sqrt(0); factor({p})"));
        let equal = Expr::call(Symbol::EQUAL, [theirs.clone(), ours.clone()]);
        assert!(
            theirs == ours
                || ctx.eval(&equal).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE),
            "factor({p}): giac's {answer}, ours {ours:?}"
        );
    }
}

/// The flag set back, `factor` splits over square roots again, as by
/// default.
#[test]
fn square_roots_by_default() {
    let split = session_value("with_sqrt(0); with_sqrt(1); factor(x^2-2)");
    assert_eq!(split, session_value("factor(x^2-2)"));
    assert_ne!(split, session_value("with_sqrt(0); factor(x^2-2)"));
}
