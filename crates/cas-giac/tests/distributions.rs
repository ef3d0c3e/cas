//! Probability distributions against giac (CAS-DOC-23 §7, CAS-PLAN-17):
//! this CAS's exact answers against giac's, exactly where giac's are
//! exact. Where giac's are doubles (Poisson's), which it prints with 12
//! digits and `evalf(…, 30)` does not extend, within a relative 10⁻¹¹.
//! giac's answers at a symbolic point and its unevaluated incomplete `Beta`
//! are the unit tests' alone.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// Whether two numbers agree to 10⁻²⁰, real and imaginary parts apart.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// Whether two answers are one value: their difference expanded is 0,
/// `==` holds, or they agree to 10⁻²⁰ (equal values in different forms,
/// `erf(sqrt(2)/2)` against `erf(2^(-1/2))`).
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    let equal = Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0))
        || ctx.eval(&equal).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
        || close(ctx, a, b)
}

/// This CAS's answer to each `ours` input against giac's to its `theirs`.
fn check_pairs(pairs: &[(&str, &str)]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, theirs) in pairs {
        match (giac.value(theirs), ours(input)) {
            (Ok(t), Ok(o)) if same_value(&ctx, &o, &t) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {theirs} {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The same input to both.
fn check(inputs: &[&str]) {
    let pairs: Vec<(&str, &str)> = inputs.iter().map(|i| (*i, *i)).collect();
    check_pairs(&pairs);
}

/// Whether `ours` is within a relative 10⁻¹¹ of giac's double `theirs`.
fn near_double(ctx: &Context, ours: &Expr, theirs: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [ours.clone() - theirs.clone(), Expr::from(30)],
        )],
    );
    let scale = Expr::call(
        Symbol::MAX,
        [Expr::from(1), Expr::call(Symbol::ABS, [theirs.clone()])],
    );
    let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-11) * scale]);
    ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// giac's doubles against this CAS's exact answers to `inputs`.
fn check_doubles(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(t), Ok(o)) if near_double(&ctx, &o, &t) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Binomial, geometric and negative binomial: giac's exact answers
/// (CAS-301).
#[test]
fn exact_discrete() {
    check(&[
        "binomial(10,3/10,2)",
        "binomial(10,3/10,11)",
        "binomial_cdf(10,3/10,2)",
        "binomial_cdf(10,3/10,2,4)",
        "binomial_cdf(10,1/2,10)",
        "binomial_icdf(10,3/10,0.5)",
        "binomial_icdf(10,1/2,0.999)",
        "binomial_icdf(10,1/2,1)",
        "binomial_icdf(10,1/2,0)",
        "geometric(1/3,2)",
        "geometric_cdf(1/3,2)",
        "geometric_cdf(1/3,2,4)",
        "geometric_icdf(1/3,1/2)",
        "negbinomial(3,1/2,2)",
        "poisson_icdf(2,0.5)",
        "poisson_icdf(1/2,0.9)",
    ]);
}

/// Poisson's masses and cumulative probabilities: exact here, decimals
/// in giac (CAS-301).
#[test]
fn poisson_decimals() {
    check_doubles(&[
        "poisson(2,3)",
        "poisson(1/2,3)",
        "poisson_cdf(2,3)",
        "poisson_cdf(1,2)",
        "poisson_cdf(1,1,2)",
    ]);
}

/// The normal, exponential and uniform distributions, giac's exact
/// answers (CAS-303).
#[test]
fn elementary_continuous() {
    check(&[
        "normald(0,1,x)",
        "normald(1/2)",
        "normald_cdf(0,1,1)",
        "normald_cdf(1)",
        "normald_cdf(0,1,-1,1)",
        "normald_cdf(1,2)",
        "exponentiald(2,1)",
        "exponentiald_cdf(2,1)",
        "exponentiald_cdf(2,1,3)",
        "exponentiald_icdf(2,1/2)",
        "uniformd(0,2,1)",
        "uniformd_cdf(0,2,1)",
        "uniformd_icdf(0,2,1/4)",
    ]);
    check_doubles(&["normald_icdf(0,1,0.975)", "normald_icdf(0.975)"]);
}
