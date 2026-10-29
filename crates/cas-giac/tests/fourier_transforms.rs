//! Fourier transforms against giac (CAS-DOC-41, CAS-PLAN-35): giac's
//! `fourier` and `ifourier`, read back, are this CAS's by value at points,
//! where a delta's argument is not 0 and so the delta is.

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

/// Whether `ours` and `theirs` agree within 10⁻¹⁰ at `variable =
/// numerator/denominator`, parts apart.
fn agree(
    ctx: &Context,
    ours: &Expr,
    theirs: &Expr,
    variable: &str,
    (numerator, denominator): (i64, i64),
) -> bool {
    let point = Expr::from(numerator) / Expr::from(denominator);
    let at = (ours.clone() - theirs.clone()).substitute(&[(Expr::symbol(variable), point)]);
    let Ok(value) = ctx.eval(&Expr::call(Symbol::N, [at])) else {
        return false;
    };
    let parts = match value.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![value],
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts.into_iter().all(|part| {
        let gap = Expr::call(Symbol::ABS, [part]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// giac's answers and this CAS's agree in `v` at the points.
fn check(inputs: &[&str], v: &str) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours))
                if !ours.contains(&Expr::from(Symbol::FOURIER_TRANSFORM))
                    && !ours.contains(&Expr::from(Symbol::INVERSE_FOURIER_TRANSFORM))
                    && [(1, 3), (7, 5), (-5, 2)]
                        .into_iter()
                        .all(|point| agree(&ctx, &ours, &theirs, v, point)) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The table and the rules (CAS-439, CAS-440).
#[test]
fn table() {
    check(
        &[
            "fourier(exp(-x^2),x,s)",
            "fourier(exp(-2*x^2),x,s)",
            "fourier(exp(-(x-1)^2),x,s)",
            "fourier(x^2*exp(-x^2),x,s)",
            "fourier(x*exp(-x^2),x,s)",
            "fourier(exp(-x^2+2*x),x,s)",
            "fourier(exp(-x^2)*cos(x),x,s)",
            "fourier(exp(-x^2)*sin(x),x,s)",
            "fourier(exp(i*x)*exp(-x^2),x,s)",
            "fourier(exp(-abs(x)),x,s)",
            "fourier(exp(-abs(x-1)),x,s)",
            "fourier(exp(-abs(2*x)),x,s)",
            "fourier(exp(-abs(x))*cos(x),x,s)",
            "fourier(Heaviside(x)*exp(-x),x,s)",
            "fourier(Heaviside(x)*exp(-2*x),x,s)",
            "fourier(Heaviside(-x)*exp(x),x,s)",
            "fourier(Heaviside(x)*exp(-x)*cos(x),x,s)",
            "fourier(Heaviside(x),x,s)",
            "fourier(Heaviside(x-1),x,s)",
            "fourier(sign(x),x,s)",
            "fourier(sign(x-1),x,s)",
            "fourier(rect(x),x,s)",
            "fourier(rect(2*x+1),x,s)",
            "fourier(1/x,x,s)",
            "fourier(sin(x)/x,x,s)",
            "fourier(Dirac(x-2),x,s)",
        ],
        "s",
    );
}

/// Rational functions by residues (CAS-441).
#[test]
fn rational_functions() {
    check(
        &[
            "fourier(1/(1+x^2),x,s)",
            "fourier(1/(x^2+4),x,s)",
            "fourier(x/(x^2+1),x,s)",
            "fourier(1/(x^2+2*x+2),x,s)",
            "fourier(1/(1+x^4),x,s)",
        ],
        "s",
    );
}

/// Inverse transforms (CAS-442).
#[test]
fn inverse_transforms() {
    check(
        &[
            "ifourier(2/(1+s^2),s,x)",
            "ifourier(pi*exp(-abs(s)),s,x)",
            "ifourier(sqrt(pi)*exp(-s^2/4),s,x)",
        ],
        "x",
    );
}
