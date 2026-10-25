//! Integrals against giac (CAS-DOC-9, CAS-PLAN-3): an antiderivative here
//! and giac's differ at most by a constant, and definite integrals have the
//! same value. Each step of the plan adds its cases.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, OracleError, ours};
use proptest::prelude::*;

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

/// Whether `ours` and `theirs` differ by a constant in `x`: they change
/// as much from `x = 3/7` to `x = 5/2`. Not their difference at the two
/// points, a ball about 0 when they are the same function written apart,
/// which `Equal` cannot tell from 0.
fn constant(ctx: &Context, ours: &Expr, theirs: &Expr) -> bool {
    let x = Expr::symbol("x");
    let change = |e: &Expr| {
        let point = |v: Expr| e.substitute(&[(x.clone(), v)]);
        let change = point(Expr::from(5) / Expr::from(2)) - point(Expr::from(3) / Expr::from(7));
        at(ctx, &change, None)
    };
    match (change(ours), change(theirs)) {
        (Some(a), Some(b)) => equal(ctx, a, b),
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
    let agree = if definite {
        // The values, not their difference: a ball about 0 that `Equal`
        // cannot tell from 0.
        match (at(&ctx, &ours, None), at(&ctx, &theirs, None)) {
            (Some(a), Some(b)) => equal(&ctx, a, b),
            _ => false,
        }
    } else {
        constant(&ctx, &ours, &theirs)
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

/// Definite integrals of CAS-198: continuous integrands between rational
/// or real bounds, the same value as giac's.
#[test]
fn definite_integrals() {
    check(
        &[
            "integrate(1/x,x,1,2)",
            "integrate(1/x,x,2,1)",
            "integrate(1/(1+x^2),x,0,1)",
            "integrate(1/(x^2+x+1),x,0,1)",
            "integrate(1/(x^2-2),x,0,1)",
            "integrate(1/(x^2-1),x,2,3)",
            "integrate(1/(x^3+1),x,0,1)",
            "integrate(x/(x^2+1),x,0,2)",
            "integrate(sin(x),x,0,pi)",
            "integrate(cos(x),x,0,pi/2)",
            "integrate(exp(x),x,0,1)",
            "integrate(sqrt(x),x,0,1)",
            "integrate(ln(x),x,1,2)",
            "integrate(tan(x),x,0,1)",
            "integrate(1/sqrt(1-x^2),x,0,1/2)",
            "integrate(x^2,x,0,sqrt(2))",
            "integrate(x/(x^2+1),x,a,b)",
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

/// Rational functions of CAS-196: Hermite's reduction, logarithms with
/// rational residues, and what is left meeting the table.
#[test]
fn rational_functions() {
    check(
        &[
            "integrate(1/(x^2-1),x)",
            "integrate(x/(x^2+1),x)",
            "integrate(1/(x*(x+1)^2),x)",
            "integrate(1/(x^2*(x+1)),x)",
            "integrate(1/(x^2+1)^2,x)",
            "integrate(1/(x^4-1),x)",
            "integrate((x^3+1)/(x-1),x)",
            "integrate(x^3/(x^2+1),x)",
            "integrate((3*x^2+1)/(x^3+x),x)",
            "integrate((2*x+1)/(x^2+x),x)",
            "integrate(1/((2*x-1)*(x+1)),x)",
            "integrate((x^4+x+1)/(x^3-x),x)",
            "integrate(1/(x^2-4*x+3)^2,x)",
        ],
        false,
    );
}

/// The heuristics of CAS-199: derivative-divides, parts, an exponential
/// times a sine or cosine, products of sines and cosines; giac's answers up
/// to a constant (`1/(x ln x)` aside: giac's `ln(ln(x))` is complex below
/// 1, where this CAS's `ln|ln x|` is real).
#[test]
fn heuristics() {
    check(
        &[
            "integrate(2*x*exp(x^2),x)",
            "integrate(x*exp(x^2),x)",
            "integrate(cos(x)*exp(sin(x)),x)",
            "integrate(sqrt(x^2+1)*x,x)",
            "integrate(ln(x)/x,x)",
            "integrate(sin(x)*cos(x)^2,x)",
            "integrate(exp(sqrt(x))/sqrt(x),x)",
            "integrate(cos(ln(x))/x,x)",
            "integrate(x*exp(x),x)",
            "integrate(x^2*sin(x),x)",
            "integrate(x*cos(2*x+1),x)",
            "integrate(x^3*exp(-x),x)",
            "integrate(x*ln(x),x)",
            "integrate(x^2*ln(x),x)",
            "integrate(ln(x)^2,x)",
            "integrate(atan(x),x)",
            "integrate(x*atan(x),x)",
            "integrate(asin(x),x)",
            "integrate(exp(x)*sin(x),x)",
            "integrate(exp(2*x)*cos(3*x),x)",
            "integrate(sin(x)*cos(x),x)",
            "integrate(cos(x)^2,x)",
            "integrate(sin(x)^2,x)",
            "integrate(sin(x)^3,x)",
            "integrate(sin(2*x)*cos(3*x),x)",
            "integrate(x*sin(x)^2,x)",
            "integrate(exp(x)*sin(x)*cos(x),x)",
        ],
        false,
    );
}

/// Pairs of quadratic residues of CAS-197: real ones joined in one
/// logarithm, complex ones as Rioboo's arctangents of polynomials, which
/// giac's answers here match up to a constant.
#[test]
fn quadratic_residues() {
    check(
        &[
            "integrate(1/(x^2+x+1),x)",
            "integrate(x/(x^2+x+1),x)",
            "integrate(1/(x^2+x+1)^2,x)",
            "integrate(1/(3*x^2+2*x+1),x)",
            "integrate((2*x+3)/(x^2+2*x+5),x)",
            "integrate(1/(x^2-x-1),x)",
            "integrate(2*x/(x^2-x-1),x)",
            "integrate(1/(x^2-2),x)",
            "integrate((x+1)/(x^2-3),x)",
            "integrate(1/(x^3+1),x)",
            "integrate(x/(x^4+1),x)",
            "integrate((8*x^2-24)/(x^4-2*x^2+9),x)",
            "integrate((x^4-3*x^2+6)/(x^6-5*x^4+5*x^2+4),x)",
        ],
        false,
    );
}

/// `p / ((x − r₁)^e₁ ⋯ (x² + b x + c))` as Xcas text: a numerator of small
/// coefficients, linear factors with multiplicities, at most one quadratic,
/// whose roots are integers when they are rational, so that neither point
/// of [`constant`] is a pole.
fn rational_function() -> impl Strategy<Value = String> {
    (
        prop::collection::vec(-5i64..6, 1..5),
        prop::collection::vec((-3i64..4, 1u32..3), 0..3),
        prop::option::of((-3i64..4, -3i64..5)),
    )
        .prop_map(|(numerator, roots, quadratic)| {
            let numerator: Vec<String> = numerator
                .iter()
                .enumerate()
                .map(|(k, c)| format!("({c})*x^{k}"))
                .collect();
            let mut factors: Vec<String> = roots
                .iter()
                .map(|(r, e)| format!("(x-({r}))^{e}"))
                .collect();
            if let Some((b, c)) = quadratic {
                factors.push(format!("(x^2+({b})*x+({c}))"));
            }
            factors.push("1".to_owned());
            format!(
                "integrate(({})/({}),x)",
                numerator.join("+"),
                factors.join("*")
            )
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Random products of CAS-197: linear factors and a quadratic, whose
    /// residues are rational or a pair, with any numerator; giac's answer
    /// up to a constant.
    #[test]
    fn random_rational_functions(input in rational_function()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agrees(giac, &input, false).map_err(TestCaseError::fail)?;
    }
}
