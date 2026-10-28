//! Interpolation and orthogonal polynomials against giac (CAS-DOC-27 §6,
//! CAS-PLAN-21): polynomials and rational functions compared after
//! `Expand` and `Together` of the difference, lists element by element,
//! since giac's Newton forms nest otherwise.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// Whether `a − b` is 0 once expanded and put together, lists element by
/// element.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| same(ctx, x, y));
    }
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

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

/// `lagrange` and `interp` in their forms (CAS-328).
#[test]
fn lagrange() {
    check(&[
        "lagrange([1,2,3],[4,8,14])",
        "lagrange([0,1,2,3],[1,0,2,1])",
        "lagrange([1,2],[a,b])",
        "lagrange([1/2,1,3],[1,-1,2],t)",
        "lagrange([[1,2,3],[1,4,9]])",
        "interp([1,2,3],[4,8,14])",
    ]);
}

/// Natural splines of degrees 1 to 5 (CAS-329).
#[test]
fn splines() {
    check(&[
        "spline([0,1,2],[1,3,0],x,3)",
        "spline([0,1,2,3],[0,1,0,1],x,1)",
        "spline([0,1,2],[1,3,0],x,2)",
        "spline([0,1,2,3],[1,0,2,1],x,2)",
        "spline([0,1,2,3,4],[1,0,2,1,3],x,4)",
        "spline([0,1,2,3,4],[1,0,2,1,3],x,5)",
        "spline([0,1,3],[1,2,0],t,3)",
        "spline([-1,0,1/2,2],[2,-1,0,3],x,3)",
    ]);
}

/// Thiele's continued fractions (CAS-330).
#[test]
fn thiele() {
    check(&[
        "thiele([1,2,3],[1,4,9],x)",
        "thiele([0,1,2,3],[1,2,5,10],x)",
        "thiele([1,2,3,4],[1,1/2,1/3,1/4],x)",
        "thiele([0,1,2],[1,3,2],t)",
    ]);
}

/// The orthogonal polynomials, and values (CAS-331). `laguerre` with its
/// parameter: giac's two-argument form keeps a symbol `a`.
#[test]
fn orthogonal_polynomials() {
    check(&[
        "legendre(5,x)",
        "hermite(6,x)",
        "laguerre(4,x,0)",
        "laguerre(3,x,b)",
        "tchebyshev1(7,x)",
        "tchebyshev2(6,x)",
        "legendre(10,1/3)",
        "hermite(5,3/2)",
        "tchebyshev1(9,y)",
    ]);
}
