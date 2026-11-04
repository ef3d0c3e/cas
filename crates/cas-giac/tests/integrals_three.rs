//! Integration III against giac (CAS-DOC-49 §4, CAS-PLAN-43): giac's
//! definite integrals into the incomplete gamma and beta functions, and
//! its antiderivatives by their differences between two points (the
//! constants of integration apart), against this CAS's by value within
//! 10⁻¹⁵. giac's `// ∫ ~=` notes are not part of its answers.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether `a` and `b` are one number within 10⁻¹⁵, by balls.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let numeric = |e: &Expr| ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]));
    let (Ok(a), Ok(b)) = (numeric(a), numeric(b)) else {
        return false;
    };
    let gap = Expr::call(Symbol::ABS, [a - b]);
    let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-15)]);
    ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Definite integrals (CAS-488).
#[test]
fn definite() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "int(x^(3/2)*exp(-x),x,0,inf)",
        "int(x^(1/3)*exp(-2*x),x,0,inf)",
        "int(x^(1/3)*(1-x)^(2/3),x,0,1)",
        "int(x^(1/2)*(1-x)^(1/2),x,0,1)",
        "int(x^(2/3)*exp(-x),x,1,2)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(t), Ok(o)) if near(&ctx, &o, &t) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Antiderivatives (CAS-487), by `F(2) − F(1)`.
#[test]
fn antiderivatives() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for integrand in ["x^(5/2)*exp(-2*x)", "x^(1/3)*exp(-2*x)", "x^(-1/2)*exp(-x)"] {
        let input = format!("int({integrand},x)");
        let difference = |f: &Expr| {
            let at = |p: i64| f.substitute(&[(Expr::symbol("x"), Expr::from(p))]);
            at(2) - at(1)
        };
        match (giac.value(&input), ours(&input)) {
            (Ok(t), Ok(o))
                if !o.contains(&Expr::from(Symbol::INTEGRATE))
                    && near(&ctx, &difference(&o), &difference(&t)) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
