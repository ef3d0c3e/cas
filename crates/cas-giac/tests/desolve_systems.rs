//! Systems of differential equations against giac (CAS-DOC-20 §5,
//! CAS-PLAN-14). giac solves a linear system with constant coefficients
//! in matrix form, `desolve(y'=A*y+b,y)`, and answers `[]` to the same as
//! a list of scalar equations, which this CAS's language writes: so each
//! of this CAS's systems is compared with giac's matrix form of it,
//! component by component, at points where the constants have values. A
//! higher-order system is compared with giac's first-order form, the
//! unknowns' components picked from it; both take the initial values at 0
//! as the constants.

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// The points where answers are compared, the constants' values among
/// them: before and after the steps at 1.
const POINTS: [[(&str, i64, i64); 7]; 2] = [
    [
        ("x", 2, 9),
        ("a", 3, 4),
        ("c_0", 7, 5),
        ("c_1", -2, 3),
        ("c_2", 5, 11),
        ("c_3", 3, 13),
        ("c_4", -9, 7),
    ],
    [
        ("x", 7, 3),
        ("a", 3, 4),
        ("c_0", 7, 5),
        ("c_1", -2, 3),
        ("c_2", 5, 11),
        ("c_3", 3, 13),
        ("c_4", -9, 7),
    ],
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

/// Whether two numbers agree: each part of `a − b` within `10⁻¹⁰ (1 + |b|)`.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let Ok(difference) = ctx.eval(&(a.clone() - b.clone())) else {
        return false;
    };
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(10_000_000_000_i64);
    parts(&difference).into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, bound.clone()]))
    })
}

/// `e` approximated at `point`.
fn at(ctx: &Context, e: &Expr, point: &[(&str, i64, i64)]) -> Option<Expr> {
    let rules: Vec<(Expr, Expr)> = point
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
        .ok()
}

/// The components of an answer `{{v₁, …, vₙ}}`.
fn components(answer: &Expr) -> Option<&[Expr]> {
    let list = Expr::from(Symbol::LIST);
    match answer.args() {
        [solution] if answer.head() == Some(&list) && solution.head() == Some(&list) => {
            Some(solution.args())
        }
        _ => None,
    }
}

/// Whether this CAS's answer is giac's components `picked`, at every
/// point.
fn same(ctx: &Context, ours: &Expr, theirs: &Expr, picked: &[usize]) -> bool {
    let (Some(ours), Some(theirs)) = (components(ours), components(theirs)) else {
        return false;
    };
    ours.len() == picked.len()
        && ours.iter().zip(picked).all(|(a, &k)| {
            theirs.get(k).is_some_and(|b| {
                POINTS.iter().all(|point| {
                    matches!(
                        (at(ctx, a, point), at(ctx, b, point)),
                        (Some(a), Some(b)) if close(ctx, &a, &b)
                    )
                })
            })
        })
}

/// A system: this CAS's input, giac's matrix form, and the components of
/// giac's answer that are this CAS's, in order.
struct Case(&'static str, &'static str, &'static [usize]);

/// Each case's answers, this CAS's and giac's, the same by [`same`].
fn check(cases: &[Case]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for &Case(input, matrix, picked) in cases {
        match (ours(input), giac.value(matrix)) {
            (Ok(ours), Ok(theirs)) if same(&ctx, &ours, &theirs, picked) => {}
            (ours, theirs) => wrong.push(format!(
                "{input}: ours {ours:?}, giac's {matrix}: {theirs:?}"
            )),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// First-order systems, homogeneous and forced, and a higher-order one
/// (CAS-283). The oracle's libgiac answers the system forced by
/// `Heaviside(x-1)` with `erf(1)` and `sqrt(pi)`, wrongly, where giac 1.9
/// agrees with this CAS: the unit tests'.
#[test]
fn constant_coefficients() {
    check(&[
        Case(
            "desolve([diff(y(x),x)=z(x),diff(z(x),x)=-y(x)],[y(x),z(x)])",
            "desolve(y'=[[0,1],[-1,0]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=y(x)+z(x),diff(z(x),x)=y(x)-z(x)],[y(x),z(x)])",
            "desolve(y'=[[1,1],[1,-1]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=2*y(x),diff(z(x),x)=3*z(x)],[y,z])",
            "desolve(y'=[[2,0],[0,3]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=3*y(x)-2*z(x),diff(z(x),x)=2*y(x)-2*z(x)],x,[y,z])",
            "desolve(y'=[[3,-2],[2,-2]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=y(x)+z(x),diff(z(x),x)=z(x)],[y(x),z(x)])",
            "desolve(y'=[[1,1],[0,1]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=y(x)-z(x),diff(z(x),x)=y(x)+z(x)],[y(x),z(x)])",
            "desolve(y'=[[1,-1],[1,1]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=y(x)+2*z(x),diff(z(x),x)=3*y(x)+4*z(x)],[y(x),z(x)])",
            "desolve(y'=[[1,2],[3,4]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=a*y(x),diff(z(x),x)=y(x)],[y(x),z(x)])",
            "desolve(y'=[[a,0],[1,0]]*y,y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=z(x),diff(z(x),x)=w(x),diff(w(x),x)=y(x)],[y(x),z(x),w(x)])",
            "desolve(y'=[[0,1,0],[0,0,1],[1,0,0]]*y,y)",
            &[0, 1, 2],
        ),
        // Forced.
        Case(
            "desolve([diff(y(x),x)=z(x),diff(z(x),x)=-y(x)+1],[y(x),z(x)])",
            "desolve(y'=[[0,1],[-1,0]]*y+[0,1],y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=z(x)+x,diff(z(x),x)=-y(x)+exp(x)],[y(x),z(x)])",
            "desolve(y'=[[0,1],[-1,0]]*y+[x,exp(x)],y)",
            &[0, 1],
        ),
        // y'' = −z, z' = y against Y = (y, y', z).
        Case(
            "desolve([diff(y(x),x,x)=-z(x),diff(z(x),x)=y(x)],[y(x),z(x)])",
            "desolve(y'=[[0,1,0],[0,0,-1],[1,0,0]]*y,y)",
            &[0, 2],
        ),
    ]);
}

/// Conditions fix the constants (CAS-284).
#[test]
fn conditions() {
    check(&[
        Case(
            "desolve([diff(y(x),x)=z(x),diff(z(x),x)=-y(x),y(0)=1,z(0)=0],[y(x),z(x)])",
            "desolve([y'=[[0,1],[-1,0]]*y,y(0)=[1,0]],y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=y(x)+z(x),diff(z(x),x)=y(x)-z(x),y(0)=2,z(0)=-1],[y(x),z(x)])",
            "desolve([y'=[[1,1],[1,-1]]*y,y(0)=[2,-1]],y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=2*y(x),diff(z(x),x)=3*z(x),y(1)=1,z(1)=1],[y(x),z(x)])",
            "desolve([y'=[[2,0],[0,3]]*y,y(1)=[1,1]],y)",
            &[0, 1],
        ),
        Case(
            "desolve([diff(y(x),x)=z(x)+x,diff(z(x),x)=-y(x)+exp(x),y(0)=0,z(0)=1],[y(x),z(x)])",
            "desolve([y'=[[0,1],[-1,0]]*y+[x,exp(x)],y(0)=[0,1]],y)",
            &[0, 1],
        ),
        // y'' = -z, z' = y against Y = (y, y', z), y'(0) by subst.
        Case(
            "desolve([diff(y(x),x,x)=-z(x),diff(z(x),x)=y(x),y(0)=1,\
             subst(diff(y(x),x),x=0)=0,z(0)=2],[y(x),z(x)])",
            "desolve([y'=[[0,1,0],[0,0,-1],[1,0,0]]*y,y(0)=[1,0,2]],y)",
            &[0, 2],
        ),
    ]);
}
