//! Integrals against giac (CAS-DOC-9, CAS-PLAN-3): an antiderivative here
//! and giac's differ at most by a constant, and definite integrals have the
//! same value. Each step of the plan adds its cases.

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

/// The parameters of the inputs, besides `x`, and their values where the
/// two answers are compared.
const PARAMETERS: [(&str, i64, i64); 4] = [("a", 2, 3), ("b", 5, 4), ("y", 7, 5), ("n", 3, 2)];

/// `e` with the parameters at their values and `x` at `x`, as balls where
/// it is not exact.
fn at(ctx: &Context, e: &Expr, x: Option<Expr>) -> Option<Expr> {
    let mut rules: Vec<(Expr, Expr)> = PARAMETERS
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    if let Some(x) = x {
        rules.push((Expr::symbol("x"), x));
    }
    ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
        .ok()
}

/// Whether `left == right` is `True`: exact numbers equal, balls
/// overlapping.
fn equal(ctx: &Context, left: Expr, right: Expr) -> bool {
    ctx.eval(&Expr::call(Symbol::EQUAL, [left, right]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// Whether `difference` is constant in `x`, at two points.
fn constant(ctx: &Context, difference: &Expr) -> bool {
    let first = at(ctx, difference, Some(Expr::from(3) / Expr::from(7)));
    let second = at(ctx, difference, Some(Expr::from(5) / Expr::from(2)));
    match (first, second) {
        (Some(first), Some(second)) => equal(ctx, first, second),
        _ => false,
    }
}

/// giac's answer to `input` and this CAS's: the same, or with the
/// parameters at their values, for an antiderivative a constant apart, for
/// a definite integral the same value.
fn agrees(giac: &Giac, input: &str, definite: bool) -> Result<(), String> {
    match giac.agrees(input) {
        Ok(()) => return Ok(()),
        Err(OracleError::Differs { .. }) => {}
        Err(e) => return Err(e.to_string()),
    }
    let theirs = giac.value(input).map_err(|e| e.to_string())?;
    let ours = ours(input).map_err(|e| e.to_string())?;
    let ctx = Context::new();
    let difference = ours.clone() - theirs.clone();
    let agree = if definite {
        at(&ctx, &difference, None).is_some_and(|d| equal(&ctx, d, Expr::from(0)))
    } else {
        constant(&ctx, &difference)
    };
    if agree {
        Ok(())
    } else {
        Err(format!("{input}: ours {ours:?}, giac's {theirs:?}"))
    }
}

/// Each of `inputs` against giac.
fn check(inputs: &[&str], definite: bool) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| agrees(giac, input, definite).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The antiderivatives of CAS-194: what is free of `x`, sums, constant
/// factors, powers of `x`, polynomials, lists, and parts that stay.
#[test]
fn polynomials_and_powers() {
    check(
        &[
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
        ],
        false,
    );
}

/// The definite integrals of CAS-194, of polynomials.
#[test]
fn definite_integrals_of_polynomials() {
    check(
        &[
            "integrate(x^2,x,0,1)",
            "integrate((x+1)^2,x,-1,2)",
            "integrate(x^3-x,x,a,b)",
        ],
        true,
    );
}

/// The table at linear arguments of CAS-195 (CAS-DOC-9 §4).
#[test]
fn the_table_at_linear_arguments() {
    check(
        &[
            "integrate(exp(x),x)",
            "integrate(exp(2*x+1),x)",
            "integrate(exp(-x),x)",
            "integrate(exp(a*x),x)",
            "integrate(2^x,x)",
            "integrate(3^(2*x),x)",
            "integrate(sin(x),x)",
            "integrate(cos(3*x),x)",
            "integrate(sin(a*x+b),x)",
            "integrate(cos(x/2),x)",
            "integrate(tan(x),x)",
            "integrate(tan(2*x),x)",
            "integrate(ln(x),x)",
            "integrate(ln(2*x+1),x)",
            "integrate(1/(2*x+3),x)",
            "integrate(1/(x+a),x)",
            "integrate((2*x+1)^(-3),x)",
            "integrate(sqrt(2*x+1),x)",
            "integrate(1/sqrt(1-x^2),x)",
            "integrate(1/sqrt(1-4*x^2),x)",
            "integrate(1/sqrt(4-x^2),x)",
            "integrate(1/(1+x^2),x)",
            "integrate(1/(1+(2*x+1)^2),x)",
            "integrate(1/(1+4*x^2),x)",
            "integrate(3*cos(x)+x^2,x)",
        ],
        false,
    );
}
