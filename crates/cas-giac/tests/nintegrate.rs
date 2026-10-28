//! Certified numerical integration against giac (CAS-DOC-26 §6,
//! CAS-PLAN-20): giac's `romberg`, `gaussquad` and `nInt` answer doubles,
//! which this CAS's balls must agree with to `10⁻¹⁰`, relatively.

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

/// Whether `a` is within `10⁻¹⁰ max(1, |b|)` of `b`.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(Symbol::N, [a.clone() - b.clone()])],
    );
    let scale = Expr::call(
        Symbol::MAX,
        [
            Expr::from(1),
            Expr::call(Symbol::ABS, [Expr::call(Symbol::N, [b.clone()])]),
        ],
    );
    let claim = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-10) * scale]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Smooth integrands, a switch, products. Not `romberg(sqrt(x),x,0,1)`,
/// where giac gives up and answers its last two approximations.
#[test]
fn quadratures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "romberg(sin(sin(x)),x,0,1)",
        "gaussquad(exp(-x^2),x,0,1)",
        "nInt(1/(1+x^4),x,0,1)",
        "romberg(ln(1+x)/(1+x^2),x,0,1)",
        "gaussquad(cos(x^2),x,0,2)",
        "nInt(exp(x)*sin(5*x),x,-1,1)",
        "romberg(abs(sin(x)),x,0,3)",
        "nInt(atan(x)/(1+x^2),x,-2,1)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
