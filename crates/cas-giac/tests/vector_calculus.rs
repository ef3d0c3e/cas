//! Vector calculus and curves against giac (CAS-DOC-37, CAS-PLAN-31):
//! giac's answers, read back, are this CAS's by value at a point, lists
//! entry by entry. Vector potentials differ by gradients, so theirs is
//! checked through its curl.

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

/// The symbols' values where answers are compared.
const POINT: [(&str, i64, i64); 8] = [
    ("p", 3, 7),
    ("r", 5, 11),
    ("s", 7, 5),
    ("v", 2, 3),
    ("t", 5, 7),
    ("x", 2, 9),
    ("y", 7, 3),
    ("z", -3, 4),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The parts of a number, a complex ball's real and imaginary ones apart.
fn parts(e: &Expr) -> Vec<Expr> {
    match e.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![e.clone()],
    }
}

/// Whether two values are one: lists entry by entry, the same, or each
/// part of their difference within 10⁻¹⁰ once approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    if a == b {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = POINT
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let difference = (a.clone() - b.clone()).substitute(&rules);
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [difference])) else {
        return false;
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts(&difference).into_iter().all(|p| {
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

/// The differential operators and the scalar potential (CAS-408,
/// CAS-409).
#[test]
fn operators() {
    check(&[
        "grad(2*x^2*y-x*z^3,[x,y,z])",
        "divergence([x^2+y,x+z+y,z^3+x^2],[x,y,z])",
        "divergence([x*y,y*z],[x,y])",
        "curl([2*x*y,x*z,y*z],[x,y,z])",
        "curl([x*y,y*z,z*x],[x,y,z])",
        "laplacian(exp(z)*cos(x*y),[x,y,z])",
        "laplacian(3)",
        "hessian(2*x^2*y-x*z,[x,y,z])",
        "hessian(x*y*z,[x,y,z])",
        "potential([2*x*y+3,x^2-4*z,-4*y],[x,y,z])",
        "potential([2*x*y,x^2+3*y^2],[x,y])",
        "potential([y,x+z,y],[x,y,z])",
    ]);
}

/// Vector potentials: giac's and this CAS's differ by a gradient, so each
/// one's curl is the field.
#[test]
fn vector_potentials() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for (field, vars) in [
        ("[2*x*y+3,x^2-4*z,-2*y*z]", "[x,y,z]"),
        ("[x,y,-2*z]", "[x,y,z]"),
    ] {
        let input = format!("curl(vpotential({field},{vars}),{vars})");
        let theirs = giac.value(&input);
        let ours = ours(&input);
        match (theirs, ours) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => panic!("{input}: ours {ours:?}, giac's {theirs:?}"),
        }
    }
}

/// Curves (CAS-410, CAS-411). giac's `frenet` answers a sequence, read as
/// the list it makes in brackets.
#[test]
fn curves() {
    check(&[
        "curvature([t,t^2],t)",
        "curvature([t,-t^2],t)",
        "curvature([t,t^2],t,1)",
        "curvature(x^2,x)",
        "curvature([5*cos(t),5*sin(t)],t)",
        "curvature([3*exp(t/2)*cos(t),3*exp(t/2)*sin(t)],t)",
        "curvature([cos(t),sin(t),t],t,0)",
        "arcLen(t^2,t,1,2)",
        "arcLen([t,t^2],t,1,2)",
        "arcLen([cos(t),sin(t)],t,1,2)",
        "arcLen([cos(t),sin(t),t],t,0,1)",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for input in [
        "frenet([t,t^2],t,1)",
        "frenet([t,-t^2],t,1)",
        "frenet([cos(t),sin(t),t],t,0)",
    ] {
        match (giac.value(&format!("[{input}]")), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => panic!("{input}: ours {ours:?}, giac's {theirs:?}"),
        }
    }
}

/// Implicit differentiation (CAS-412).
#[test]
fn implicit_derivatives() {
    check(&[
        "implicitdiff(x^2+y^2=1,y,x)",
        "implicitdiff(x^2+y^2=1,x,y)",
        "implicitdiff(x^2*y+y^2=1,y,x)",
        "implicitdiff(r=p*v/s,p,s)",
        "implicitdiff(x^2+y^2=1,y,x,x)",
        "implicitdiff(x^2*y+y^3=2,y,x,x)",
    ]);
}
