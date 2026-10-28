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

/// Whether two lists of polynomials generate one ideal: the reduced basis
/// of giac's, under `revlex` in the symbols of both in canonical order, as
/// `Eliminate` orders them, is ours. giac mostly eliminates by resultants,
/// which it says "may reduce further", and normalizes nothing (CAS-292).
fn same_ideal(ctx: &Context, ours: &Expr, theirs: &Expr) -> bool {
    let mut symbols: Vec<Symbol> = [ours, theirs]
        .iter()
        .flat_map(|e| e.symbols())
        .filter(|s| !s.is_builtin())
        .collect();
    symbols.sort_by_cached_key(|&s| Expr::from(s));
    symbols.dedup();
    if symbols.is_empty() {
        return ours == theirs;
    }
    let vars = Expr::call(Symbol::LIST, symbols.into_iter().map(Expr::from));
    let basis = Expr::call(Symbol::GROEBNER_BASIS, [theirs.clone(), vars]);
    ctx.eval(&basis).is_ok_and(|b| same_basis(ctx, ours, &b))
}

/// Elimination ideals (CAS-292).
#[test]
fn eliminations() {
    check(
        &[
            "eliminate([x=t^2,y=t^3],t)",
            "eliminate([x=t^2,y=t^3],[t])",
            "eliminate(x=t^2,t)",
            "eliminate([x^2+y^2+z^2=1,x+y+z=0,x*y=z],[z])",
            "eliminate([x^2+y^2+z^2=1,x+y+z=0,x*y=z],[y,z])",
            "eliminate([x+y=1,x-y=3],[y])",
            "eliminate([x+y=1,x+y=2],[y])",
            "eliminate([x*y-1,x^2-2],[x])",
            "eliminate([x=a*t,y=b*t^2],t)",
            "eliminate([x+y+z,x*y+y*z+z*x,x*y*z-1],[x,y])",
            "eliminate([x^2-y,y^2-z],[y])",
            "eliminate([x-y,y-z],[x,y,z])",
            "eliminate([x^2-1],[y])",
            "eliminate([x^2+1,y-x],[x])",
            "eliminate([x^2+y^2-1,x-y,x+z],[x])",
            "eliminate([x*y-1,y*z-1,x*z-1],[x])",
            "eliminate([x^3-y,x^2-z],[x])",
            "eliminate([u=x+y,v=x*y,x-y=1],[x,y])",
            "eliminate([x-t,y-t^2,z-t^3],t)",
            "eliminate([x=u*v,y=u^2,z=v^2],[u,v])",
            "eliminate([x^2+y^2=r^2,y=m*x],[y])",
            "eliminate([a*x+b,c*x+d],x)",
        ],
        same_ideal,
    );
}

/// Whether two numbers agree to 10⁻²⁰, real and imaginary parts apart.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// Whether two answers are one set of tuples: as many, each of one close
/// to one of the other's, value by value.
fn same_tuples(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    a.head() == Some(&list)
        && b.head() == Some(&list)
        && a.args().len() == b.args().len()
        && a.args().iter().all(|s| {
            b.args().iter().any(|t| {
                s.args().len() == t.args().len()
                    && s.args().iter().zip(t.args()).all(|(u, v)| close(ctx, u, v))
            })
        })
}

/// Polynomial systems by lexicographic bases (CAS-289), as sets of tuples.
#[test]
fn systems() {
    check(
        &[
            "solve([x^2+y^2+z^2=1,x-y=0,y-z=0],[x,y,z])",
            "solve([x*y=1,y*z=1,x*z=1],[x,y,z])",
            "solve([x+y+z=1,x^2+y^2+z^2=1,x^3+y^3+z^3=1],[x,y,z])",
            "solve([x^2-1,y^2-1,x-y],[x,y])",
            "solve([x+y=1,x+y=2,x*y=3],[x,y])",
            "solve([x+y+z=0,x*y+y*z+z*x=0,x*y*z=1],[x,y,z])",
            "csolve([x+y+z=0,x*y+y*z+z*x=0,x*y*z=1],[x,y,z])",
            "solve([x^2+y^2=2,x*y=1,x+z=0],[x,y,z])",
            "solve([x^2+y+z=1,x+y^2+z=1,x+y+z^2=1],[x,y,z])",
        ],
        same_tuples,
    );
}
