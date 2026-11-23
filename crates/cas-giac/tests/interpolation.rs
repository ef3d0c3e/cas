//! Interpolation and orthogonal polynomials against giac (CAS-DOC-27 §6,
//! CAS-PLAN-21): polynomials and rational functions compared after
//! `Expand` and `Together` of the difference, lists element by element,
//! since giac's Newton forms nest otherwise; `ratinterp` and `triginterp`
//! too (CAS-594), decimal `triginterp` by value at three points.

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

/// Floater and Hormann's rational interpolants (CAS-594): the default `d`,
/// `d = 1` and `d = 0`, the matrix form, values at points.
#[test]
fn rational_interpolants() {
    check(&[
        "ratinterp([1,3,5,8],[2,-1,3,4],x)",
        "ratinterp([1,3,5,8],[2,-1,3,4],x,1)",
        "ratinterp([[1,2],[3,-1],[5,3],[8,4]],x)",
        "ratinterp([0,1,2,3,4],[1,2,1,2,1],x,1)",
        "ratinterp([0,1,2,3,4],[1,2,1,2,1],x,0)",
        "ratinterp([0,1,2,3],[1,0,1,0],x,1)",
        "ratinterp([1,2,3],[1,4,9])",
        "ratinterp([1,2,3],[1,4,9],[0,1/2,4])",
        "ratinterp([[1,2],[3,-1],[5,3],[8,4]],[2,4],1)",
    ]);
}

/// Trigonometric interpolation (CAS-594), exact in its two forms.
#[test]
fn trigonometric_interpolants() {
    check(&[
        "triginterp([1,2,3],x=0..2)",
        "triginterp([11,10,17,24,32,26,23,19],x=0..21)",
        "triginterp([11,10,17,24,32,26,23,19],0,21,x)",
        "triginterp([1,0,1,0],x=0..3)",
        "triginterp([0,1,0,-1],x=0..3)",
    ]);
}

/// Trigonometric interpolation of decimals, giac's doubles: compared by
/// value at three points of `x`, within `10⁻¹⁰` (C0006).
#[test]
fn trigonometric_interpolants_of_decimals() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let mut wrong = Vec::new();
    for input in [
        "triginterp([1.0,2,3],x=0..2)",
        "triginterp([1,2,3],x=0..2.0)",
    ] {
        let (Ok(theirs), Ok(ours)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        for point in [Expr::from(3) / 10, Expr::from(17) / 10, Expr::from(29) / 10] {
            let at = |e: &Expr| e.substitute(&[(x.clone(), point.clone())]);
            let gap = Expr::call(
                Symbol::ABS,
                [Expr::call(Symbol::N, [at(&ours) - at(&theirs)])],
            );
            let claim = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-10)]);
            if ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) != Some(Symbol::TRUE) {
                wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
                break;
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where giac errs or answers `undef` and this CAS leaves the call: giac's
/// answer as caseval prints it, whitespace collapsed (C0007, C0047).
#[test]
fn interpolant_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        (
            "ratinterp([1,2,3],[1,4,9],x,5)",
            "GIAC_ERROR: Expected an integer 0<=d<=|X|. Error: Bad Argument Value",
            Symbol::RATIONAL_INTERPOLATION,
        ),
        (
            "ratinterp([1],[1],x)",
            "GIAC_ERROR: At least two sample points are required. Error: Bad Argument Value",
            Symbol::RATIONAL_INTERPOLATION,
        ),
        (
            "ratinterp([1,2,3],[1,4,9],x,3)",
            "undef",
            Symbol::RATIONAL_INTERPOLATION,
        ),
        (
            "ratinterp([1,1,3],[1,4,9],x)",
            "undef",
            Symbol::RATIONAL_INTERPOLATION,
        ),
        (
            "triginterp([1],x=0..2)",
            "GIAC_ERROR: At least two data points are required. Error: Bad Argument Value",
            Symbol::TRIGONOMETRIC_INTERPOLATION,
        ),
        (
            "triginterp([1,2,3],x=0..0)",
            "undef",
            Symbol::TRIGONOMETRIC_INTERPOLATION,
        ),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, head) in cases {
        let printed = giac.eval(input).unwrap();
        if collapsed(&printed) != collapsed(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let stays = ours(input)
            .ok()
            .and_then(|e| e.head().and_then(Expr::as_symbol));
        if stays != Some(head) {
            wrong.push(format!("{input}: ours {:?} does not stay", ours(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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
