//! The Bessel and Airy functions against giac (CAS-DOC-34 §7,
//! CAS-PLAN-28): giac evaluates `BesselJ`, `BesselY`, `Airy_Ai` and
//! `Airy_Bi` in double precision, exact arguments too, and knows nothing
//! else of them: ours, approximated, are within its printed digits.

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

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two numbers agree to the digits giac prints: `a − b` within
/// `10⁻¹⁰ (1 + |b|)`, both approximated.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let at = |e: &Expr| ctx.eval(&Expr::call(Symbol::N, [e.clone()])).ok();
    let (Some(a), Some(b)) = (at(a), at(b)) else {
        return false;
    };
    if !matches!(a.view(), View::Ball(_)) && a.to_rational().is_none() {
        return false;
    }
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(10_000_000_000_i64);
    let gap = Expr::call(Symbol::ABS, [a - b]);
    holds(ctx, &Expr::call(Symbol::LESS, [gap, bound]))
}

/// giac's values and ours, close.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if close(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Values of `J` and `Y`, exact arguments among them, a negative one for
/// `J`, orders below and above the argument.
#[test]
fn bessel_values() {
    check(&[
        "BesselJ(0,1.0)",
        "BesselJ(1,2.5)",
        "BesselJ(2,1.0)",
        "BesselJ(5,10.0)",
        "BesselJ(0,1)",
        "BesselJ(3,0.5)",
        "BesselJ(0,-1.0)",
        "BesselJ(0,0)",
        "BesselJ(1,40.5)",
        "BesselY(0,1.0)",
        "BesselY(2,3.5)",
        "BesselY(1,0.25)",
        "BesselY(4,30.0)",
    ]);
}

/// Values of `Ai` and `Bi` on both sides of 0, and at 0.
#[test]
fn airy_values() {
    check(&[
        "Airy_Ai(0)",
        "Airy_Bi(0)",
        "Airy_Ai(1.0)",
        "Airy_Bi(1.0)",
        "Airy_Ai(2.5)",
        "Airy_Bi(2.5)",
        "Airy_Ai(-5.0)",
        "Airy_Bi(-5.0)",
        "Airy_Ai(-0.5)",
    ]);
}
