//! Exponential integrals, polygamma and Lambert W against giac (CAS-DOC-31
//! §8, CAS-PLAN-25): giac's answers, read back, are this CAS's by value,
//! the symbols at a point, decimals within giac's printed digits. Not
//! `series(LambertW(x),x=0,4)`, which giac writes without its remainder,
//! nor `Ei` and `Ci` at 0, where giac's series fail.

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
const POINT: [(&str, i64, i64); 3] = [("x", 2, 9), ("y", 7, 3), ("n", 5, 1)];

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

/// Whether two numbers agree to the digits giac prints: each part of
/// `a − b` within `10⁻¹⁰ (1 + |b|)`.
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

/// A series' terms and its remainder `x^m OrderSize[x]` apart: the
/// remainder compared as it is, the terms by value, as `N` leaves
/// `OrderSize` and the remainders of two forms do not cancel.
fn remainder(e: &Expr) -> (Expr, Vec<Expr>) {
    let order_size = Expr::from(Symbol::new("OrderSize"));
    if !e.contains(&order_size) {
        return (e.clone(), Vec::new());
    }
    let terms: Vec<Expr> = if e.head() == Some(&Expr::from(Symbol::PLUS)) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    };
    let (rest, mut kept): (Vec<Expr>, Vec<Expr>) =
        terms.into_iter().partition(|t| !t.contains(&order_size));
    kept.sort();
    (Expr::call(Symbol::PLUS, rest), kept)
}

/// Whether two values are one: lists part by part, the same, equal by `==`,
/// or equal by `==` or close once approximated at [`POINT`]; series by
/// [`remainder`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let ((a_terms, a_rest), (b_terms, b_rest)) = (remainder(a), remainder(b));
    if !a_rest.is_empty() || !b_rest.is_empty() {
        return a_rest == b_rest && same(ctx, &a_terms, &b_terms);
    }
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

/// Values and derivatives (CAS-354 to CAS-356).
#[test]
fn values() {
    check(&[
        "Ei(1.0)",
        "Ei(-2.0)",
        "Ei(0)",
        "Si(1.0)",
        "Si(0)",
        "Si(-x)",
        "Ci(1.0)",
        "Ci(-1.0)",
        "Li(2.0)",
        "Li(0.5)",
        "Li(x)",
        "Psi(1)",
        "Psi(3)",
        "Psi(1/2)",
        "Psi(5/2)",
        "Psi(-1/2)",
        "Psi(1/3)",
        "Psi(2/3)",
        "Psi(1/4)",
        "Psi(1/6)",
        "Psi(7/3)",
        "Psi(2,1)",
        "Psi(1,3)",
        "Psi(1/2,1)",
        "Psi(1,2)",
        "Psi(2.5)",
        "Psi(2.5,2)",
        "LambertW(0)",
        "LambertW(e)",
        "LambertW(-1/e)",
        "LambertW(1.0)",
        "LambertW(10.0)",
        "LambertW(-0.2,-1)",
        "evalf(euler_gamma)",
        "diff(Ei(x),x)",
        "diff(Si(x),x)",
        "diff(Ci(x),x)",
        "diff(Li(x),x)",
        "diff(Gamma(x),x)",
        "diff(Psi(x),x)",
        "diff(Psi(x,2),x)",
        "diff(LambertW(x),x)",
        "diff(Beta(x,y),x)",
    ]);
}

/// Antiderivatives and definite integrals (CAS-357, CAS-358).
#[test]
fn integrals() {
    check(&[
        "int(exp(x)/x,x)",
        "int(sin(x)/x,x)",
        "int(cos(x)/x,x)",
        "int(1/ln(x),x)",
        "int(exp(2*x)/x,x)",
        "int(exp(x)/x^2,x)",
        "int(sin(2*x)/x,x)",
        "int(sinh(x)/x,x)",
        "int(exp(-x)/x,x)",
        "int(sin(x)/x^2,x)",
        "int(cos(2*x+1)/x,x)",
        "int(exp(x)*(1+x)/x^2,x)",
        "int(x/ln(x),x)",
        "int(1/ln(2*x),x)",
        "int(exp(x^2)/x,x)",
        "int(Ei(x),x)",
        "int(Si(x),x)",
        "int(Ci(x),x)",
        "int(erf(x),x)",
        "int(LambertW(x),x)",
        "int(exp(-x^2),x)",
        "int(x^2*exp(-x^2),x)",
        "int(exp(-x^2-2*x),x)",
        "int(exp(-2*x^2+x),x)",
        "int(exp(-x^2)*x^4,x)",
        "int(sin(x)/x,x,0,inf)",
        "int(exp(-x^2/2),x,-inf,inf)",
        "int(exp(-x)/x,x,1,inf)",
        "int(exp(x)/x,x,1,2)",
        "int(sin(x)/x,x,0,1)",
    ]);
}

/// Series and limits (CAS-358).
#[test]
fn series_and_limits() {
    check(&[
        "series(Si(x),x=0,6)",
        "series(erf(x),x=0,5)",
        "series(Psi(1+x),x=0,3)",
        "series(Psi(x),x=0,2)",
        "series(Gamma(x),x=0,2)",
        "series(Si(x),x=1,2)",
        "series(Psi(x),x=-1,1)",
        "series(Gamma(x),x=-1,1)",
        "series(LambertW(x),x=1,2)",
        "limit(Si(x),x,inf)",
        "limit(Ci(x),x,inf)",
        "limit(Ei(x),x,-inf)",
        "limit(erf(x),x,inf)",
        "limit(Psi(x)-ln(x),x,inf)",
        "limit(Si(x)/x,x,0)",
        "limit(Gamma(x)*x,x,0)",
        "limit(LambertW(x)/x,x,0)",
    ]);
}

/// Sums of rational functions (CAS-359) and equations (CAS-360).
#[test]
fn sums_and_equations() {
    check(&[
        "sum(1/k,k,1,n)",
        "sum(1/(k+1),k,0,n)",
        "sum(1/(2*k+1),k,0,n)",
        "sum(1/k^2,k,1,n)",
        "sum(1/(k+1/3),k,0,n)",
        "sum(1/(k*(2*k+1)),k,1,inf)",
        "sum(1/(k*(k+1/2)),k,1,inf)",
        "sum(1/(2*k-1)^2,k,1,inf)",
        "sum(1/(3*k+1)^2,k,0,inf)",
        "solve(x*exp(x)=1,x)",
        "solve(x*exp(x)=y,x)",
        "solve((x+1)*exp(x)=2,x)",
        "solve(x=exp(-x),x)",
        "solve(x*exp(x)=-1/4,x)",
        "solve(x*exp(x)=-1,x)",
    ]);
}
