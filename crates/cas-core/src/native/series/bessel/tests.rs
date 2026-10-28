//! Tests of the Bessel and Airy series (CAS-383): the power series at 0,
//! the series at regular points against their differential equations, and
//! the limits at 0 and at infinity.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// `(input, expected)`, both evaluated: the same once put together.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = eval(&ctx, &read(input));
        let expected = eval(&ctx, &read(expected));
        let difference = eval(
            &ctx,
            &Expr::call(
                Symbol::TOGETHER,
                [Expr::call(
                    Symbol::EXPAND,
                    [found.clone() - expected.clone()],
                )],
            ),
        );
        if found != expected && difference != Expr::from(0) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The power series of `J`, `I` and Airy's functions at 0.
#[test]
fn at_zero() {
    check(&[
        (
            "Series[BesselJ[0, x], List[x, 0, 6]]",
            "Plus[1, Times[-1/4, Power[x, 2]], Times[1/64, Power[x, 4]], Times[-1/2304, Power[x, 6]], \
             Times[Power[x, 7], OrderSize[x]]]",
        ),
        (
            "Series[BesselJ[1, x], List[x, 0, 5]]",
            "Plus[Times[1/2, x], Times[-1/16, Power[x, 3]], Times[1/384, Power[x, 5]], \
             Times[Power[x, 6], OrderSize[x]]]",
        ),
        (
            "Series[BesselI[0, x], List[x, 0, 4]]",
            "Plus[1, Times[1/4, Power[x, 2]], Times[1/64, Power[x, 4]], Times[Power[x, 5], OrderSize[x]]]",
        ),
        (
            "Series[BesselI[1, Times[2, x]], List[x, 0, 3]]",
            "Plus[x, Times[1/2, Power[x, 3]], Times[Power[x, 4], OrderSize[x]]]",
        ),
        (
            "Series[AiryAi[x], List[x, 0, 3]]",
            "Plus[AiryAi[0], Times[AiryAiPrime[0], x], Times[1/6, AiryAi[0], Power[x, 3]], \
             Times[Power[x, 4], OrderSize[x]]]",
        ),
        (
            "Series[AiryBiPrime[x], List[x, 0, 2]]",
            "Plus[AiryBiPrime[0], Times[1/2, AiryBi[0], Power[x, 2]], Times[Power[x, 3], OrderSize[x]]]",
        ),
    ]);
}

/// The terms of a series, its remainder dropped.
fn terms(s: &Expr) -> Expr {
    let order_size = Expr::from(Symbol::new("OrderSize"));
    let kept: Vec<Expr> = if s.head() == Some(&Expr::from(Symbol::PLUS)) {
        s.args()
            .iter()
            .filter(|t| !t.contains(&order_size))
            .cloned()
            .collect()
    } else if s.contains(&order_size) {
        Vec::new()
    } else {
        vec![s.clone()]
    };
    Expr::call(Symbol::PLUS, kept)
}

/// At regular points, the series to order 8 satisfies its equation to
/// order 5: the coefficients of `(x − x₀)^k`, `k ≤ 5`, of the equation's
/// residual vanish.
#[test]
fn regular_points_satisfy_the_equations() {
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let h = Expr::symbol("h");
    for (f, x0, equation) in [
        // x² y″ + x y′ + (x² − n²) y, and the same with −x² for I and K.
        (
            "BesselJ[0, x]",
            "1",
            "Plus[Times[Power[x, 2], D[y, x, x]], Times[x, D[y, x]], Times[Power[x, 2], y]]",
        ),
        (
            "BesselY[1, x]",
            "2",
            "Plus[Times[Power[x, 2], D[y, x, x]], Times[x, D[y, x]], Times[Plus[Power[x, 2], -1], y]]",
        ),
        (
            "BesselI[2, x]",
            "1/2",
            "Plus[Times[Power[x, 2], D[y, x, x]], Times[x, D[y, x]], Times[-1, Plus[Power[x, 2], 4], y]]",
        ),
        (
            "BesselK[0, x]",
            "3",
            "Plus[Times[Power[x, 2], D[y, x, x]], Times[x, D[y, x]], Times[-1, Power[x, 2], y]]",
        ),
        // y″ − x y.
        ("AiryAi[x]", "-1", "Plus[D[y, x, x], Times[-1, x, y]]"),
        ("AiryBi[x]", "0", "Plus[D[y, x, x], Times[-1, x, y]]"),
    ] {
        let series = eval(
            &ctx,
            &Expr::call(
                Symbol::SERIES,
                [
                    read(f),
                    Expr::call(Symbol::LIST, [x.clone(), read(x0), Expr::from(8)]),
                ],
            ),
        );
        assert_ne!(
            series.head(),
            Some(&Expr::from(Symbol::SERIES)),
            "{f} at {x0}"
        );
        let y = terms(&series);
        let residual = read(equation).substitute(&[(Expr::symbol("y"), y)]);
        let shifted = eval(&ctx, &residual).substitute(&[(x.clone(), read(x0) + h.clone())]);
        let expanded = eval(&ctx, &Expr::call(Symbol::EXPAND, [shifted]));
        for k in 0..=5i64 {
            let c = eval(
                &ctx,
                &Expr::call(
                    Symbol::TOGETHER,
                    [Expr::call(
                        Symbol::EXPAND,
                        [Expr::call(
                            Symbol::COEFFICIENT,
                            [expanded.clone(), h.clone(), Expr::from(k)],
                        )],
                    )],
                ),
            );
            assert_eq!(
                c,
                Expr::from(0),
                "{f} at {x0}: the coefficient of h^{k} is {c:?}"
            );
        }
    }
}

/// The fuzz target `bessel`'s `Bi′` at −1/2 (CAS-386): the coefficients by
/// `Bi″ = x Bi`, and the remainder at `h = 1/32` against Airy's power
/// series in 60-digit decimals, `3.3637492563587428766…·10⁻¹²`, where the
/// `h⁵` and `h⁶` terms nearly cancel.
#[test]
fn airy_bi_prime_at_minus_half() {
    let (b, bp) = ("AiryBi[-1/2]", "AiryBiPrime[-1/2]");
    let h = "Plus[1/2, x]";
    check(&[(
        "Series[AiryBiPrime[x], List[x, -1/2, 4]]",
        &format!(
            "Plus[{bp}, Times[-1/2, {b}, {h}], \
             Times[Plus[Times[1/2, {b}], Times[-1/4, {bp}]], Power[{h}, 2]], \
             Times[Plus[Times[1/24, {b}], Times[1/3, {bp}]], Power[{h}, 3]], \
             Times[Plus[Times[-1/12, {b}], Times[1/96, {bp}]], Power[{h}, 4]], \
             Times[Power[{h}, 5], OrderSize[{h}]]]"
        ),
    )]);
    let ctx = Context::new();
    let s = eval(&ctx, &read("Series[AiryBiPrime[x], List[x, -1/2, 4]]"));
    let at = read("-15/32");
    let gap = terms(&s).substitute(&[(Expr::symbol("x"), at.clone())])
        - Expr::call(Symbol::AIRY_BI_PRIME, [at]);
    let reference = read("33637492563587428766/10000000000000000000000000000000");
    let within = Expr::call(
        Symbol::LESS,
        [
            Expr::call(
                Symbol::ABS,
                [Expr::call(Symbol::N, [gap, Expr::from(40)]) - reference],
            ),
            read("1/1000000000000000000000000"),
        ],
    );
    assert_eq!(eval(&ctx, &within).as_symbol(), Some(Symbol::TRUE), "{s:?}");
}

/// Limits at 0, through the logarithmic series, and at infinity.
#[test]
fn limits() {
    check(&[
        ("Limit[BesselJ[0, x], x, Infinity]", "0"),
        ("Limit[BesselY[1, x], x, Infinity]", "0"),
        ("Limit[BesselK[0, x], x, Infinity]", "0"),
        ("Limit[AiryAi[x], x, Infinity]", "0"),
        ("Limit[AiryAi[x], x, Times[-1, Infinity]]", "0"),
        ("Limit[AiryBi[x], x, Times[-1, Infinity]]", "0"),
        ("Limit[Times[BesselJ[1, x], Power[x, -1]], x, 0]", "1/2"),
        ("Limit[Times[x, BesselK[1, x]], x, 0, 1]", "1"),
        (
            "Limit[Times[x, BesselY[1, x]], x, 0, 1]",
            "Times[-2, Power[Pi, -1]]",
        ),
        ("Limit[BesselY[0, x], x, 0, 1]", "Times[-1, Infinity]"),
    ]);
}
