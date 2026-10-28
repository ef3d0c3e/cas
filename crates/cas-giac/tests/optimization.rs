//! Optimization against giac (CAS-DOC-38, CAS-PLAN-32): giac's answers to
//! its examples of `extrema`, `minimize`, `maximize` and `lpsolve`, read
//! back, are this CAS's by value, lists entry by entry.

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists and equations part by part, the
/// same, or numbers within 10⁻¹⁰.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    if a.head().is_some()
        && a.head() == b.head()
        && a.head()
            .and_then(Expr::as_symbol)
            .is_some_and(|h| matches!(h, Symbol::LIST | Symbol::EQUATION))
    {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [a.clone() - b.clone()])) else {
        return false;
    };
    let parts = match difference.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![difference],
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts.into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
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

/// Local extrema (CAS-415).
#[test]
fn extrema() {
    check(&[
        "extrema(x^3-3*x,x)",
        "extrema(x^2+y^2-2*x,[x,y])",
        "extrema(x^3+y^3-3*x*y,[x,y])",
        "extrema(x^4,x)",
        "extrema(x*y*z,x+y+z=1,[x,y,z])",
        "extrema(x^2+y^2,x*y=1,[x=0..inf,y=0..inf])",
    ]);
}

/// Global optima (CAS-416).
#[test]
fn optima() {
    check(&[
        "minimize(x^2-2*x+3,x)",
        "maximize(sin(x),x=0..4)",
        "minimize(x^3-3*x,x=-3..3)",
        "minimize(x^2+y^2,x+y=1,[x,y])",
        "minimize(2x^2-y^2+6y,x^2+y^2<=16,[x,y])",
        "maximize(2x^2-y^2+6y,x^2+y^2<=16,[x,y])",
        "minimize(x^2-3x+y^2+3y+3,[x=2..4,y=-4..-2],point)",
    ]);
}

/// Linear programs (CAS-417, CAS-418).
#[test]
fn linear_programs() {
    check(&[
        "lpsolve(2x+y-z+4,[x<=1,y>=2,x+3y-z=2,2x-y+z<=8,-x+y<=5])",
        "lpsolve(-7x+2y,[4x-12y<=20,-x+3y<=3],x=-5..5,y=0..inf,maximize)",
        "lpsolve(-x-y,[y<=3x+1/2,y<=-5x+2],assume=lp_nonnegative)",
        "lpsolve(-6x+4y+z,[5x-10y<=20,2z-3y=6,-x+3y<=3],x=1..20,y=0..inf)",
        "lpsolve(-5x-7y,[7x+y<=35,-x+3y<=6],assume=integer)",
        "lpsolve(8x1+11x2+6x3+4x4,[5x1+7x2+4x3+3x4<=14],assume=lp_binary,maximize)",
        "lpsolve(x+y,[x>=1/3,y>=1/2],assume=integer)",
    ]);
}
