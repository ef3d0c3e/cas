//! Polynomial tools against giac (CAS-DOC-29 §4, CAS-PLAN-23):
//! `partfrac` and `poly2symb` compared after `Expand` and `Together` of the
//! difference, the others exactly. Not `symb2poly`, whose `poly1[…]`
//! reads otherwise here (`coeffs` gives the same list); not `realroot`,
//! whose intervals the unit tests check; not `sturmab` of a polynomial
//! with a multiple root, which giac counts otherwise.

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

#[test]
fn polynomial_tools() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "partfrac(x/(x^2-3*x+2))",
        "partfrac((x^2+1)/(x^3-x))",
        "partfrac((x^3+1)/(x*(x-1)^2))",
        "partfrac(1/(x^3+1))",
        "partfrac((2*x+3)/((x+1)^2*(x^2+1)))",
        "partfrac((x^4+1)/(x^2-1))",
        "partfrac(1/(x^2-2))",
        "partfrac(1/((x-a)*(x-b)),x)",
        "partfrac(1/(x^2+1)/(x-1)^2,x)",
        "coeffs(x^3-2*x+1)",
        "coeffs(x^2*y+y,x)",
        "poly2symb([1,0,-2,1],x)",
        "poly2symb([2,0,1],t)",
        "tcoeff(3*x^4+2*x^2)",
        "valuation(3*x^4+2*x^2)",
        "sturmab(x^3-2*x,x,-2,2)",
        "sturmab(x^2+1,x,-5,5)",
        "sturmab(x^4-5*x^2+4,x,0,3)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
