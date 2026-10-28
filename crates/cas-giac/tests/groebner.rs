//! Gröbner bases against giac (CAS-DOC-21 §5, CAS-PLAN-15): reduced bases
//! compared as sets of polynomials, each up to a constant factor, as giac
//! normalizes neither the signs of its elements nor their order; `greduce`'s
//! remainders exactly.

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

/// Whether `a` is a non-zero constant multiple of `b`: their quotient, put
/// together, a number.
fn proportional(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let quotient = Expr::call(Symbol::TOGETHER, [a.clone() / b.clone()]);
    ctx.eval(&quotient)
        .is_ok_and(|q| q.is_number() && q != Expr::from(0))
}

/// Whether two answers are one basis: lists of as many elements, each of
/// one a constant multiple of one of the other's.
fn same_basis(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    a.head() == Some(&list)
        && b.head() == Some(&list)
        && a.args().len() == b.args().len()
        && a.args()
            .iter()
            .all(|p| b.args().iter().any(|q| proportional(ctx, p, q)))
}

/// Whether two answers are one expression, their difference expanded 0.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0))
}

/// giac's answers and this CAS's to `inputs`, alike by `alike`.
fn check(inputs: &[&str], alike: fn(&Context, &Expr, &Expr) -> bool) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if alike(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Reduced bases under each order (CAS-288).
#[test]
fn bases() {
    check(
        &[
            "gbasis([x^2+y^2-1,x-y],[x,y])",
            "gbasis([x^2+y^2-1,x-y],[x,y],plex)",
            "gbasis([x^2+y^2-1,x-y],[x,y],tdeg)",
            "gbasis([x^2+y^2=1,x=y],[x,y])",
            "gbasis([x^3-2*x*y,x^2*y-2*y^2+x],[x,y])",
            "gbasis([x^3-2*x*y,x^2*y-2*y^2+x],[x,y],plex)",
            "gbasis([2*x+4*y,6*x*y],[x,y])",
            "gbasis([x^2-y,x^3-z],[x,y,z],plex)",
            "gbasis([x^2+y,x*y+1],[x,y],tdeg)",
            "gbasis([x-1/2,y],[x,y])",
            "gbasis([x^2*y-1,x*y^2-x],[x,y])",
            "gbasis([x+y+z,x*y+y*z+z*x,x*y*z-1],[x,y,z],plex)",
            "gbasis([x*y-z^2,y^2-x*z,x^2*y-z^3],[x,y,z])",
            "gbasis([x^2+y^2+z^2-1,x*y-z,x-y+z],[x,y,z],plex)",
            "gbasis([x^3-y,y^2-x*z],[x,y,z],tdeg)",
            "gbasis([x^2-1,x-2],[x])",
            "gbasis([x^2-2],[x,y])",
            "gbasis([0],[x])",
        ],
        same_basis,
    );
}

/// Remainders of division by lists (CAS-288).
#[test]
fn remainders() {
    check(
        &[
            "greduce(x^3+y^2,[2*y^2-1,x-y],[x,y])",
            "greduce(x*y,[2*y^2-1,x-y],[x,y])",
            "greduce(x^2,[x-y],[x,y],plex)",
            "greduce(x^2*y,[x*y-1],[x,y])",
        ],
        same,
    );
}
