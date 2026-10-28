//! Tests of the closed forms of the regularized incomplete gamma and beta
//! functions (CAS-304): known values, and the reflection `I_x(a, b) = 1 −
//! I_{1−x}(b, a)` at half-integers.

use super::{beta, lower_gamma};
use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The value of `found` is that of `expected`, both evaluated.
fn same(ctx: &Context, found: &Expr, expected: &str) {
    let found = ctx.eval(found).unwrap();
    let expected = ctx.eval(&read(expected)).unwrap();
    assert_eq!(found, expected);
}

/// Whether `e` is within `10^-25` of `value`.
fn near(ctx: &Context, e: &Expr, value: &str) -> bool {
    let claim = Expr::call(
        crate::symbol::Symbol::LESS,
        [
            Expr::call(
                crate::symbol::Symbol::ABS,
                [Expr::call(
                    crate::symbol::Symbol::N,
                    [e.clone() - read(value), Expr::from(30)],
                )],
            ),
            Expr::from(10).pow(-25),
        ],
    );
    ctx.eval(&claim).unwrap() == read("True")
}

#[test]
fn gamma() {
    let ctx = Context::new();
    let x = read("x");
    same(
        &ctx,
        &lower_gamma(2, &x).unwrap(),
        "Plus[1, Times[-1, Exp[Times[-1, x]]]]",
    );
    same(&ctx, &lower_gamma(1, &x).unwrap(), "Erf[Power[x, 1/2]]");
    // P(2, 1) = 1 − 2/e: χ²(4) at 2.
    same(
        &ctx,
        &lower_gamma(4, &Expr::from(1)).unwrap(),
        "Plus[1, Times[-2, Power[E, -1]]]",
    );
    // P(3/2, 1/2) = erf(1/√2) − √(2/π) e^(−1/2): χ²(3) at 1.
    let half = read("1/2");
    assert!(near(
        &ctx,
        &lower_gamma(3, &half).unwrap(),
        "Plus[Erf[Power[2, -1/2]], Times[-1, Power[Times[2, Power[Pi, -1]], 1/2], Power[E, -1/2]]]",
    ));
    assert_eq!(lower_gamma(0, &x), None);
}

#[test]
fn incomplete_beta() {
    let ctx = Context::new();
    // I_{1/3}(1, 2) = 5/9: giac's Beta(1,2,1/3,1).
    same(&ctx, &beta(2, 4, &read("1/3")).unwrap(), "5/9");
    // I_{1/2}(1/2, 1/2) = 1/2, and I_{1/4}(1/2, 3/2) = 1/3 + √3/(2π).
    same(&ctx, &beta(1, 1, &read("1/2")).unwrap(), "1/2");
    assert!(near(
        &ctx,
        &beta(1, 3, &read("1/4")).unwrap(),
        "Plus[1/3, Times[1/2, Power[3, 1/2], Power[Pi, -1]]]",
    ));
    // The reflection at half-integers and mixed parities.
    for (a, b) in [(3, 5), (1, 4), (4, 1), (5, 6), (7, 7)] {
        let sum = beta(a, b, &read("2/5")).unwrap() + beta(b, a, &read("3/5")).unwrap();
        assert!(near(&ctx, &sum, "1"), "{a}/2, {b}/2");
    }
    assert_eq!(beta(0, 2, &read("1/2")), None);
}
