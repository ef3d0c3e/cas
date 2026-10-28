//! Laplace transforms against giac (CAS-DOC-18, CAS-PLAN-12): giac's
//! answers, read back with `Heaviside`, `Dirac` and `Gamma` named, are this
//! CAS's by value at a point.

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
const POINT: [(&str, i64, i64); 6] = [
    ("x", 7, 3),
    ("s", 11, 4),
    ("t", 13, 5),
    ("a", 3, 4),
    ("b", 2, 7),
    ("n", 5, 2),
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

/// Whether two answers are one: the same, equal by `==`, or close once
/// approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let equal = |a: &Expr, b: &Expr| holds(ctx, &Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]));
    if a == b || equal(a, b) {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = POINT
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let at = |e: &Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
            .ok()
    };
    match (at(a), at(b)) {
        (Some(a), Some(b)) => a == b || equal(&a, &b) || close(ctx, &a, &b),
        _ => false,
    }
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

/// The step and the delta (CAS-273).
#[test]
fn heaviside_and_dirac() {
    check(&[
        "Heaviside(0)",
        "Heaviside(2)",
        "Heaviside(-1)",
        "Dirac(1)",
        "diff(Heaviside(x),x)",
        "integrate(Dirac(x-1)*sin(x),x,0,3)",
        "integrate(Dirac(x),x,-1,1)",
    ]);
}

/// The transform (CAS-274). giac's `laplace(t^(1/2),t,s)` is 0 and its
/// transforms of derivatives stay as integrals: the unit tests'.
#[test]
fn transforms() {
    check(&[
        "laplace(1,t,s)",
        "laplace(t^3,t,s)",
        "laplace(exp(-2*t),t,s)",
        "laplace(sin(3*t),t,s)",
        "laplace(cos(3*t),t,s)",
        "laplace(t*exp(t),t,s)",
        "laplace(exp(t)*sin(t),t,s)",
        "laplace(sinh(2*t),t,s)",
        "laplace(cosh(t),t,s)",
        "laplace(t*sin(t),t,s)",
        "laplace(a*t+b,t,s)",
        "laplace(sin(t)^2,t,s)",
        "laplace(t^n,t,s)",
        "laplace(exp(a*t),t,s)",
        "laplace(Heaviside(t-1),t,s)",
        "laplace(Dirac(t-2),t,s)",
        "laplace(t*Heaviside(t-1),t,s)",
        "laplace(sin(t-1)*Heaviside(t-1),t,s)",
        "laplace(exp(-t)*t^2,t,s)",
        "laplace(t*cos(2*t),t,s)",
        "laplace(1)",
        "laplace(exp(-x))",
    ]);
}
