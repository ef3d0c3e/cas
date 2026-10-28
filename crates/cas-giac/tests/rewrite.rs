//! The rewriting forms against giac (CAS-DOC-36, CAS-PLAN-30): giac's
//! answers, read back, are this CAS's by value, the symbols at a point
//! where the formal rules hold (every symbol positive, below 1 where an
//! arcsine needs it). giac's forms are not this CAS's (`evalc` writes
//! `exp(x)` through `tan(x/2)` there), only their values are.

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

/// The symbols' values where answers are compared.
const POINT: [(&str, i64, i64); 4] = [("a", 3, 7), ("b", 5, 11), ("x", 2, 9), ("y", 7, 3)];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The parts of a number, a complex ball's real and imaginary ones apart:
/// the modulus of a complex ball about 0 is unbounded.
fn parts(e: &Expr) -> Vec<Expr> {
    match e.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![e.clone()],
    }
}

/// Whether two values are one: the same, or each part of their difference
/// within 10⁻¹⁰ once approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = POINT
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let difference = (a.clone() - b.clone()).substitute(&rules);
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [difference])) else {
        return false;
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts(&difference).into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// giac's answers and this CAS's, the same by [`same`].
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Exponentials, powers and logarithms (CAS-396).
#[test]
fn exponentials_and_logarithms() {
    check(&[
        "hyp2exp(sinh(x))",
        "hyp2exp(cosh(x)+tanh(x))",
        "hyp2exp(asinh(x))",
        "hyp2exp(acosh(y))",
        "hyp2exp(atanh(x))",
        "exp2pow(exp(3*ln(x)))",
        "exp2pow(exp(a*ln(x)+1))",
        "exp2pow(exp(ln(x)/2))",
        "exp2pow(exp(x*ln(x)))",
        "exp2pow(exp(2*ln(x)+3*ln(y)))",
        "pow2exp(x^y)",
        "pow2exp(2^x)",
        "pow2exp(x^2)",
        "pow2exp(sqrt(x))",
        "expexpand(exp(x+y))",
        "expexpand(exp(3*x))",
        "expexpand(exp(x+2*y))",
        "expexpand(exp(x-y))",
        "expexpand(exp(3*x+1))",
        "lnexpand(ln(x*y))",
        "lnexpand(ln(x^2/y))",
        "lnexpand(ln(sqrt(x)))",
        "lnexpand(ln(2*x*y^3))",
        "lnexpand(ln(2/3))",
        "lncollect(ln(x)+ln(y))",
        "lncollect(2*ln(x)-ln(y))",
        "lncollect(ln(x)/2+ln(y)/3)",
        "lncollect(a*ln(x)+ln(y))",
        "lncollect(ln(2)+ln(3))",
        "lncollect(3*ln(2))",
        "lncollect(-2*ln(x))",
    ]);
}

/// Trigonometric substitutions and Pythagorean forms (CAS-397, CAS-398).
#[test]
fn trigonometric_forms() {
    check(&[
        "halftan(sin(x))",
        "halftan(cos(x))",
        "halftan(tan(x))",
        "halftan(sin(2*x)/(1+cos(2*x)))",
        "halftan(cos(x)+sin(y))",
        "halftan_hyp2exp(sin(x)+sinh(x))",
        "tan2sincos(tan(x))",
        "sin2costan(sin(x))",
        "cos2sintan(cos(x))",
        "tan2sincos2(tan(x))",
        "tan2cossin2(tan(x))",
        "sincos(exp(i*x))",
        "trigcos(sin(x)^4+sin(x)^2)",
        "trigcos(sin(x)^3)",
        "trigcos(tan(x))",
        "trigcos(sin(x)^2+tan(x)^2)",
        "trigsin(cos(x)^4+cos(x)^2)",
        "trigsin(cos(x)^3*sin(x))",
        "trigtan(cos(x)^4+sin(x)^2)",
        "trigtan(sin(x)*cos(x))",
        "trigtan(sin(x))",
        "trigtan(sin(x)^2/cos(x)^2)",
    ]);
}

/// `tcollect` (CAS-399).
#[test]
fn collected() {
    check(&[
        "tcollect(sin(x)+cos(x))",
        "tcollect(sin(x)*cos(x))",
        "tcollect(cos(x)^2-sin(x)^2)",
        "tcollect(3*sin(x)+4*cos(x))",
        "tcollect(sin(2*x)+cos(2*x)+sin(x))",
        "tcollect(sin(x)-cos(x))",
        "tcollect(-sin(x)-cos(x))",
        "tcollect(sin(x)^2)",
        "tcollect(sin(x)+sqrt(3)*cos(x))",
    ]);
}

/// Inverse trigonometric functions (CAS-400).
#[test]
fn inverse_functions() {
    check(&[
        "acos2asin(acos(x))",
        "acos2asin(acos(x)+asin(x))",
        "acos2atan(acos(x))",
        "asin2acos(asin(x))",
        "asin2atan(asin(x))",
        "atan2acos(atan(x))",
        "atan2asin(atan(x))",
        "atrig2ln(atan(x))",
        "atrig2ln(asin(x))",
        "atrig2ln(acos(x))",
    ]);
}

/// Complex parts (CAS-401).
#[test]
fn complex_parts() {
    check(&[
        "evalc(1/(x+y*i))",
        "evalc((x+i*y)^2)",
        "evalc((x+i*y)^3)",
        "evalc(exp(x+i*y))",
        "evalc(sin(x+i*y))",
        "evalc(cos(x+i*y))",
        "evalc(tan(x+i*y))",
        "evalc(abs(x+i*y))",
        "evalc(sinh(x+i*y))",
        "evalc(conj(x+i*y))",
        "evalc(-3+4*i+exp(i))",
        "evalc(ln(x+i*y))",
        "mult_conjugate((2-sqrt(2))/(2+sqrt(2)))",
        "mult_conjugate(1/(sqrt(x)+1))",
        "mult_conjugate(sqrt(3)-sqrt(2))",
        "mult_conjugate(1/(sqrt(3)-sqrt(2)))",
        "mult_conjugate(x/(sqrt(x)-2))",
        "mult_c_conjugate(1/(3+2*i))",
        "mult_c_conjugate(1/(x+i*y))",
        "mult_c_conjugate((1+i)/(2-i))",
    ]);
}

/// `combine`, `tsimplify`, `trigsimplify` and the Heaviside forms
/// (CAS-402).
#[test]
fn combinations() {
    check(&[
        "combine(ln(x)+ln(y),ln)",
        "combine(2*ln(x)-ln(y),ln)",
        "combine(exp(x)*exp(y),exp)",
        "combine(exp(x)^2,exp)",
        "combine(sin(x)*cos(x),trig)",
        "combine(sin(x)^2,trig)",
        "tsimplify(exp(2*x)+2*exp(x)+1)",
        "tsimplify(sin(x)^2+cos(x)^2)",
        "tsimplify((exp(x)-exp(-x))/(exp(x)+exp(-x)))",
        "tsimplify(cosh(x)^2-sinh(x)^2)",
        "tsimplify(sin(x))",
        "trigsimplify(sin(x)^2+cos(x)^2)",
        "trigsimplify(3*sin(x)-4*sin(x)^3)",
        "trigsimplify(1-cos(x)^2)",
        "sign2Heaviside(sign(x))",
        "Heaviside2sign(Heaviside(x))",
    ]);
}
