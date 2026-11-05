//! Piecewise functions against giac (CAS-DOC-53, CAS-PLAN-47): giac's
//! answers, read back, are this CAS's.

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

/// Whether two answers are one: the same, or equal by `==`.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    a == b
        || ctx
            .eval(&Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]))
            .ok()
            .and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
}

/// giac's answers and this CAS's, the same by [`same`].
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Evaluation (CAS-513): by substitution, decided conditions and defaults.
/// No case applying is giac's string and `undef` here, so not compared.
#[test]
fn evaluation() {
    check(&[
        "subst(piecewise(x<0,-1,x<1,0),x=1/2)",
        "subst(piecewise(x<0,-1,x<1,0),x=-2)",
        "subst(piecewise(x<0,1,x>=0 and x<2,x,4),x=1)",
        "subst(piecewise(x<0,1,x>=0 and x<2,x,4),x=7)",
        "subst(piecewise(x<0,-x,x),x=-3)",
        "subst(piecewise(x^2<1,0,1),x=1/2)",
        "subst(piecewise(x<0,-1,1),x=0)",
        "piecewise(true,1,2)",
        "piecewise(false,1,2)",
        "piecewise(1)",
        "piecewise(-2<0,-1,1)",
    ]);
}
