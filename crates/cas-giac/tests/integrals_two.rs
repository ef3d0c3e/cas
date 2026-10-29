//! Integration II against giac (CAS-DOC-40, CAS-PLAN-34): for giac's
//! integrals of rational functions of exponentials, of sines and cosines,
//! and of roots, this CAS's antiderivative is whole and differs from
//! giac's by a constant, their difference's derivative 0 at points where
//! both are real; definite integrals by value.

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

/// Whether `e`, a number, is within 10⁻¹⁰ of 0, its parts apart.
fn near_zero(ctx: &Context, e: &Expr) -> bool {
    let Ok(value) = ctx.eval(&Expr::call(Symbol::N, [e.clone()])) else {
        return false;
    };
    let parts = match value.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![value],
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts.into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// For each `int(f,x)`: this CAS's antiderivative is whole, and its
/// difference with giac's has the derivative 0 at the points `p/q`.
fn check_antiderivatives(inputs: &[&str], points: &[(i64, i64)]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let mut wrong = Vec::new();
    for input in inputs {
        let (theirs, ours) = match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) => (theirs, ours),
            (theirs, ours) => {
                wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
                continue;
            }
        };
        if ours.contains(&Expr::from(Symbol::INTEGRATE)) {
            wrong.push(format!("{input}: ours stays, {ours:?}"));
            continue;
        }
        let difference = Expr::call(Symbol::D, [ours.clone() - theirs.clone(), x.clone()]);
        let Ok(difference) = ctx.eval(&difference) else {
            wrong.push(format!("{input}: no derivative"));
            continue;
        };
        for &(p, q) in points {
            let at = difference.substitute(&[(x.clone(), Expr::from(p) / Expr::from(q))]);
            if !near_zero(&ctx, &at) {
                wrong.push(format!(
                    "{input}: ours {ours:?} and giac's {theirs:?} differ in slope at {p}/{q}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Rational functions of exponentials (CAS-433).
#[test]
fn exponentials() {
    check_antiderivatives(
        &[
            "int(1/(1+exp(x)),x)",
            "int(exp(x)/(1+exp(2*x)),x)",
            "int(1/(exp(x)+exp(-x)),x)",
            "int(1/cosh(x),x)",
            "int(tanh(x)^2,x)",
            "int(exp(2*x)/(exp(x)+1),x)",
            "int(sinh(x)^2,x)",
        ],
        &[(1, 3), (7, 5), (-5, 2)],
    );
}

/// Rational functions of sines and cosines (CAS-434).
#[test]
fn sines_and_cosines() {
    check_antiderivatives(
        &[
            "int(1/(2+cos(x)),x)",
            "int(1/sin(x),x)",
            "int(1/cos(x),x)",
            "int(tan(x)^2,x)",
            "int(1/(1+sin(x)),x)",
            "int(sin(x)/(1+cos(x)),x)",
            "int(sin(x)^3/cos(x)^2,x)",
            "int(1/(3+5*sin(x)),x)",
            "int(tan(x)^3,x)",
            "int(1/(1+sin(x)^2),x)",
            "int(cos(x)^2/sin(x),x)",
        ],
        &[(1, 3), (7, 5), (5, 2)],
    );
}

/// Roots of linear functions and Euler's substitutions (CAS-435).
#[test]
fn roots() {
    check_antiderivatives(
        &[
            "int(x*sqrt(x+1),x)",
            "int(1/(1+sqrt(x)),x)",
            "int(sqrt(x)/(1+x),x)",
            "int((x+1)^(1/3)*x,x)",
            "int(1/(x*sqrt(1+x^2)),x)",
            "int(x/sqrt(x^2+2*x+5),x)",
            "int(sqrt(x^2+1)/x,x)",
            "int(1/(x^2*sqrt(x^2+4)),x)",
        ],
        &[(1, 3), (7, 5), (5, 2)],
    );
    check_antiderivatives(&["int(1/(x*sqrt(x^2-1)),x)"], &[(3, 2), (7, 4), (5, 2)]);
    check_antiderivatives(&["int(x^2*sqrt(1-x^2),x)"], &[(1, 3), (-1, 2), (3, 4)]);
}

/// Definite integrals (CAS-436): giac's values.
#[test]
fn definite() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "int(1/(2+cos(x)),x,0,2*pi)",
        "int(1/(2+cos(x)),x,0,pi)",
        "int(1/(1+sin(x)^2),x,0,pi)",
        "int(1/(1+exp(x)),x,0,1)",
        "int(1/cosh(x),x,-1,1)",
        "int(x*sqrt(x+1),x,0,3)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours))
                if !ours.contains(&Expr::from(Symbol::INTEGRATE))
                    && near_zero(&ctx, &(ours.clone() - theirs.clone())) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
