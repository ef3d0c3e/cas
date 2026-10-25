//! Limits against giac (CAS-DOC-10, CAS-PLAN-4): giac's `limit`, read
//! back, is this CAS's (`infinity` is `ComplexInfinity`, `+infinity`
//! `Infinity`).

use cas_giac::Giac;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The limits of CAS-205, by series.
#[test]
fn limits_by_series() {
    check(&[
        "limit(sin(x)/x,x=0)",
        "limit((x^2-1)/(x-1),x=1)",
        "limit(1/x,x=0)",
        "limit(1/x,x=0,1)",
        "limit(1/x,x=0,-1)",
        "limit(x*ln(x),x=0,1)",
        "limit((1-cos(x))/x^2,x=0)",
        "limit((exp(x)-1-x)/x^2,x=0)",
        "limit(sqrt(x^2+x)-x,x=inf)",
        "limit(atan(x),x=-inf)",
        "limit((2*x^2+1)/(3*x^2-x),x=inf)",
        "limit(x^x,x=0,1)",
        "limit(tan(x),x=pi/2)",
        "limit(tan(x),x=pi/2,-1)",
        "limit(x/abs(x),x=0,1)",
        "limit(sqrt(x),x=0,1)",
        "limit((1+1/x)^x,x=inf)",
        "limit(ln(x)/x,x=inf)",
        "limit(x^(1/x),x=inf)",
        "limit(ln(ln(x))/ln(x),x=inf)",
    ]);
}

/// The limits of CAS-206, by Gruntz's algorithm.
#[test]
fn limits_by_gruntz() {
    check(&[
        "limit(exp(x)/x^10,x=inf)",
        "limit(x^10*exp(-x),x=inf)",
        "limit(exp(x)*x,x=-inf)",
        "limit(exp(x)/exp(2*x),x=inf)",
        "limit(exp(1/x),x=0,1)",
        "limit(exp(1/x),x=0,-1)",
        "limit(x/exp(sqrt(x)),x=inf)",
        "limit(exp(x+exp(-x))-exp(x),x=inf)",
    ]);
}
