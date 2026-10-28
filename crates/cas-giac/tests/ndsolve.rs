//! `NDSolve` against giac's `odesolve` (CAS-DOC-35 §5, CAS-PLAN-29):
//! giac integrates in doubles by a Runge–Kutta method, whose answers to
//! equations with closed forms err by up to about `2·10⁻⁹`; ours enclose
//! the solution, and the two agree within `10⁻⁸`, relative. giac 2.0.0
//! refuses the systems of its manual, `odesolve([z,-y],[t,[y,z]],…)`.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether two answers, lists of one number, agree within `10⁻⁸ (1 +
/// |b|)`.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    let (Some(a), Some(b)) = (
        (a.head() == Some(&list))
            .then(|| a.args().first())
            .flatten(),
        (b.head() == Some(&list))
            .then(|| b.args().first())
            .flatten(),
    ) else {
        return false;
    };
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(100_000_000_i64);
    let gap = Expr::call(Symbol::ABS, [a.clone() - b.clone()]);
    ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

#[test]
fn values() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "odesolve(sin(t*y),[t,y],[0,1],2)",
        "odesolve(y,t=0..1,y,1)",
        "odesolve(-2*t*y,t=0..3,y,1)",
        "odesolve(y*(1-y),t=0..5,y,1/2)",
        "odesolve(t^2-y,[t,y],[0,1],1)",
        "odesolve(y,[t,y],[1,1],0)",
        "odesolve(cos(t)*y+1,[t,y],[0,0],2)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if close(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
