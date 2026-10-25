//! Integrals against giac (CAS-DOC-9, CAS-PLAN-3): an antiderivative here
//! and giac's differ at most by a constant, and definite integrals agree.
//! Each step of the plan adds its cases.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, OracleError, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether `a − b` is constant in `x`: its values at two points, as balls
/// where they are not exact, are equal.
fn differ_by_a_constant(ours: &Expr, theirs: &Expr) -> bool {
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let difference = ours.clone() - theirs.clone();
    let at = |v: Expr| {
        let value = difference.substitute(&[(x.clone(), v)]);
        ctx.eval(&Expr::call(Symbol::N, [value]))
    };
    let (Ok(first), Ok(second)) = (
        at(Expr::from(3) / Expr::from(7)),
        at(Expr::from(5) / Expr::from(2)),
    ) else {
        return false;
    };
    ctx.eval(&Expr::call(Symbol::EQUAL, [first, second]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// giac's answer to `input` and this CAS's, the same or, for an
/// antiderivative in `x`, a constant apart.
fn agrees(giac: &Giac, input: &str) -> Result<(), String> {
    match giac.agrees(input) {
        Ok(()) => return Ok(()),
        Err(OracleError::Differs { .. }) => {}
        Err(e) => return Err(e.to_string()),
    }
    let theirs = giac.value(input).map_err(|e| e.to_string())?;
    let ours = ours(input).map_err(|e| e.to_string())?;
    if differ_by_a_constant(&ours, &theirs) {
        Ok(())
    } else {
        Err(format!("{input}: ours {ours:?}, giac's {theirs:?}"))
    }
}

/// The integrals of CAS-194: what is free of `x`, sums, constant factors,
/// powers of `x`, polynomials, lists, and definite integrals of polynomials.
#[test]
fn polynomials_and_powers() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "integrate(x^2,x)",
        "integrate(3*x^2+2*x+1,x)",
        "integrate(5,x)",
        "integrate(a,x)",
        "integrate(a*x,x)",
        "integrate(1/x,x)",
        "integrate(1/x^2,x)",
        "integrate(x^(1/3),x)",
        "integrate(sqrt(x),x)",
        "integrate(x^n,x)",
        "integrate((x+1)^2,x)",
        "integrate(x*(x+y),x)",
        "integrate((x^2-1)/(y+2),x)",
        "integrate(i*x^3,x)",
        "integrate([x,x^2],x)",
        "integrate(x+f(x),x)",
        "integrate(2*f(x),x)",
        "integrate(x^2,x,0,1)",
        "integrate((x+1)^2,x,-1,2)",
        "integrate(x^3-x,x,a,b)",
    ] {
        if let Err(e) = agrees(giac, input) {
            wrong.push(e);
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
