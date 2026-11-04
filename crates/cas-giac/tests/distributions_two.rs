//! Incomplete gamma and beta functions, and distributions II, against giac
//! (CAS-DOC-46 §5, CAS-PLAN-40): giac's exact answers by value (their
//! difference expanded, or within 10⁻²⁰), its doubles (`igamma`,
//! `ugamma` and `Beta` at decimals, the quantiles), which it prints with
//! 12 digits, within a relative 10⁻¹¹. Densities and cumulative functions
//! are taken inside their supports, where giac's formulas hold.

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

/// Whether `a` and `b` agree within `tolerance` relatively (absolutely
/// below 1), at 30 digits.
fn near(ctx: &Context, a: &Expr, b: &Expr, tolerance: i64) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [a.clone() - b.clone(), Expr::from(30)],
        )],
    );
    let scale = Expr::call(
        Symbol::MAX,
        [Expr::from(1), Expr::call(Symbol::ABS, [b.clone()])],
    );
    let small = Expr::call(Symbol::LESS, [gap, scale * Expr::from(10).pow(-tolerance)]);
    ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two answers are one value: their difference expands to 0, or
/// they agree within `tolerance`.
fn same(ctx: &Context, a: &Expr, b: &Expr, tolerance: i64) -> bool {
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0)) || near(ctx, a, b, tolerance)
}

/// This CAS's answer to each input against giac's.
fn check(inputs: &[&str], tolerance: i64) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(t), Ok(o)) if same(&ctx, &o, &t, tolerance) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The incomplete functions at decimals, giac's doubles (CAS-469).
#[test]
fn incomplete_functions() {
    check(
        &[
            "igamma(2.5,1.5)",
            "igamma(2.5,1.5,1)",
            "ugamma(2.5,1.5)",
            "ugamma(2.5,1.5,1)",
            "Gamma(2.5,1.5)",
            "igamma(0.3,0.001,1)",
            "igamma(7.25,3.5,1)",
            "ugamma(0.5,20.0,1)",
            "Beta(2.5,1.5,0.3)",
            "Beta(2.5,1.5,0.3,1)",
            "Beta(0.5,0.25,0.99,1)",
            "Beta(12.5,7.25,0.6,1)",
        ],
        11,
    );
    // Closed forms at naturals.
    check(&["igamma(2,x)", "igamma(1,3)", "igamma(3,2,1)"], 20);
}

/// The four distributions inside their supports (CAS-470).
#[test]
fn distributions() {
    check(
        &[
            "gammad(2,3,x)",
            "gammad_cdf(2,3,x)",
            "gammad_cdf(2,3,1/2)",
            "betad(2,3,x)",
            "betad_cdf(2,3,2/5)",
            "weibulld(2,3,x)",
            "weibulld(2,3,1,x)",
            "weibulld_cdf(2,3,x)",
            "weibulld_cdf(2,3,1,x)",
            "cauchyd(1,2,x)",
            "cauchyd(x)",
            "cauchyd_cdf(1,2,x)",
            "cauchyd_cdf(x)",
            "cauchyd_icdf(1,2,3/4)",
            "uniform(0,1,1/2)",
            "exponential(2,x)",
            "normal_cdf(1/2)",
        ],
        20,
    );
    check(
        &[
            "gammad(2,3,0.7)",
            "betad_cdf(2,3,0.1,0.5)",
            "gammad_cdf(2.5,1,2)",
            "gammad_icdf(2,3,0.5)",
            "betad_icdf(2,3,0.5)",
            "weibulld_icdf(2,3,0.5)",
            "weibulld_icdf(2,3,1,0.5)",
            "normal_icdf(0.9)",
        ],
        11,
    );
}
