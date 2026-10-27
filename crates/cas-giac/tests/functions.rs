//! Elementary and special functions against giac (CAS-DOC-15, CAS-PLAN-9):
//! giac's answers, read back, are this CAS's by value, the symbols at a
//! point, decimals within giac's printed digits.

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
const POINT: [(&str, i64, i64); 2] = [("x", 2, 9), ("y", 7, 3)];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The parts of a number, a complex ball's real and imaginary ones apart:
/// the modulus of a complex ball about 0 is unbounded.
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

/// Whether two values are one: lists part by part, the same, equal by `==`,
/// or equal by `==` or close once approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
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

/// The hyperbolic functions (CAS-238). giac writes the inverses of symbols
/// as logarithms; `sinh(ln(2))` stays there, and is 3/4 here.
#[test]
fn hyperbolic() {
    check(&[
        "sinh(0)",
        "cosh(0)",
        "tanh(0)",
        "sinh(-x)",
        "cosh(-x)",
        "tanh(-x)",
        "sinh(1.0)",
        "cosh(2.0)",
        "tanh(0.5)",
        "asinh(1.0)",
        "acosh(2.0)",
        "atanh(0.5)",
        "asinh(1)",
        "acosh(2)",
        "atanh(1/2)",
        "sinh(i*x)",
        "cosh(i*x)",
        "sinh(asinh(x))",
        "diff(sinh(x),x)",
        "diff(cosh(x),x)",
        "diff(tanh(x),x)",
        "diff(asinh(x),x)",
        "diff(acosh(x),x)",
        "diff(atanh(x),x)",
        "limit(sinh(x)/x,x,0)",
        "limit(tanh(x),x,inf)",
        "integrate(sinh(x),x)",
        "integrate(cosh(2*x),x)",
        "sinh(ln(2))",
    ]);
    assert_eq!(ours("sinh(ln(2))").unwrap(), ours("3/4").unwrap());
}

/// The reciprocal trigonometric functions (CAS-239).
#[test]
fn reciprocal_trigonometric() {
    check(&[
        "sec(x)",
        "csc(x)",
        "cot(x)",
        "sec(pi/3)",
        "cot(pi/4)",
        "acot(1)",
        "acot(-1)",
        "acot(0)",
        "asec(2)",
        "asec(-2)",
        "acsc(2)",
        "acsc(-1)",
    ]);
}
