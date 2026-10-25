//! Series against giac (CAS-DOC-10, CAS-PLAN-4): giac's `series` and
//! `taylor`, read back, are this CAS's, or differ from it by 0 once over
//! one denominator (giac writes `series(1/x,x=1,3)` as `1-x+1+…`).

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

/// giac's answer to `input` and this CAS's: the same, or their difference
/// 0 by `Together`.
fn agrees(giac: &Giac, input: &str) -> Result<(), String> {
    match giac.agrees(input) {
        Ok(()) => return Ok(()),
        Err(OracleError::Differs { .. }) => {}
        Err(e) => return Err(e.to_string()),
    }
    let theirs = giac.value(input).map_err(|e| e.to_string())?;
    let ours = ours(input).map_err(|e| e.to_string())?;
    let difference = Expr::call(Symbol::TOGETHER, [ours.clone() - theirs.clone()]);
    if Context::new().eval(&difference) == Ok(Expr::from(0)) {
        Ok(())
    } else {
        Err(format!("{input}: ours {ours:?}, giac's {theirs:?}"))
    }
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| agrees(giac, input).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The expansions of CAS-204: at 0, at other points, at infinity, Laurent
/// and Puiseux terms, a series that ends, parameters.
#[test]
fn expansions() {
    check(&[
        "series(exp(x),x=0,5)",
        "series(sin(x)/x,x=0,6)",
        "series(1/(1-x),x=0,4)",
        "series(ln(1+x),x=0,4)",
        "series(exp(x),x=1,3)",
        "series(1/x,x=1,3)",
        "series(1/sin(x),x=0,3)",
        "series(sqrt(1+x),x=0,3)",
        "series(tan(x),x=0,7)",
        "series(x/(exp(x)-1),x=0,4)",
        "taylor(cos(x),x=0,6)",
        "series(1/x+x,x=0,3)",
        "series(exp(1/x),x=inf,3)",
        "series(exp(a*x),x=0,3)",
        "series(atan(x),x=0,5)",
        "series(asin(x),x=0,5)",
        "series(cos(x)^2,x=0,4)",
        "series(1/(x^2*(1-x)),x=0,2)",
        "series(sin(sin(x)),x=0,5)",
        "series(ln(x),x=1,3)",
        "series(sqrt(x+x^2),x=0,2)",
    ]);
}
