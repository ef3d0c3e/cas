//! Tests of the special functions' series and limits (CAS-358): giac's
//! series at regular points and poles, ours where giac fails (`Ei` and `Ci`
//! at 0), the limits at infinity, and definite integrals through them.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same once put together.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        let difference = ctx
            .eval(&Expr::call(
                Symbol::TOGETHER,
                [Expr::call(
                    Symbol::EXPAND,
                    [found.clone() - expected.clone()],
                )],
            ))
            .unwrap();
        if found != expected && difference != Expr::from(0) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `series(Si(x),x=0,6)`, `series(erf(x),x=0,5)`,
/// `series(LambertW(x),x=0,4)`, `series(Psi(1+x),x=0,3)`,
/// `series(Psi(x),x=0,2)`, `series(Gamma(x),x=0,2)`,
/// `series(Si(x),x=1,2)`; `Ei` and `Ci` at 0, where giac fails.
#[test]
fn series_at_points() {
    check(&[
        (
            "Series[SinIntegral[x], List[x, 0, 6]]",
            "Plus[x, Times[-1/18, Power[x, 3]], Times[1/600, Power[x, 5]], Times[Power[x, 7], OrderSize[x]]]",
        ),
        (
            "Series[Erf[x], List[x, 0, 5]]",
            "Plus[Times[2, Power[Pi, -1/2], x], Times[-2/3, Power[Pi, -1/2], Power[x, 3]], Times[1/5, Power[Pi, -1/2], Power[x, 5]], Times[Power[x, 6], OrderSize[x]]]",
        ),
        (
            "Series[ExpIntegralEi[x], List[x, 0, 3]]",
            "Plus[EulerGamma, Log[x], x, Times[1/4, Power[x, 2]], Times[1/18, Power[x, 3]], Times[Power[x, 4], OrderSize[x]]]",
        ),
        (
            "Series[CosIntegral[x], List[x, 0, 4]]",
            "Plus[EulerGamma, Log[x], Times[-1/4, Power[x, 2]], Times[1/96, Power[x, 4]], Times[Power[x, 5], OrderSize[x]]]",
        ),
        (
            "Series[LambertW[x], List[x, 0, 4]]",
            "Plus[x, Times[-1, Power[x, 2]], Times[3/2, Power[x, 3]], Times[-8/3, Power[x, 4]], Times[Power[x, 5], OrderSize[x]]]",
        ),
        (
            "Series[PolyGamma[Plus[1, x]], List[x, 0, 3]]",
            "Plus[Times[-1, EulerGamma], Times[1/6, Power[Pi, 2], x], Times[-1, Zeta[3], Power[x, 2]], Times[1/90, Power[Pi, 4], Power[x, 3]], Times[Power[x, 4], OrderSize[x]]]",
        ),
        (
            "Series[PolyGamma[x], List[x, 0, 2]]",
            "Plus[Times[-1, Power[x, -1]], Times[-1, EulerGamma], Times[1/6, Power[Pi, 2], x], Times[-1, Zeta[3], Power[x, 2]], Times[Power[x, 3], OrderSize[x]]]",
        ),
        (
            "Series[Gamma[x], List[x, 0, 2]]",
            "Plus[Power[x, -1], Times[-1, EulerGamma], Times[1/12, Plus[Times[6, Power[EulerGamma, 2]], Power[Pi, 2]], x], Times[1/12, Plus[Times[-2, Power[EulerGamma, 3]], Times[-1, EulerGamma, Power[Pi, 2]], Times[-4, Zeta[3]]], Power[x, 2]], Times[Power[x, 3], OrderSize[x]]]",
        ),
        (
            "Series[SinIntegral[x], List[x, 1, 2]]",
            "Plus[SinIntegral[1], Times[Sin[1], Plus[-1, x]], Times[1/2, Plus[Cos[1], Times[-1, Sin[1]]], Power[Plus[-1, x], 2]], Times[Power[Plus[-1, x], 3], OrderSize[Plus[-1, x]]]]",
        ),
        (
            "Series[PolyGamma[x], List[x, -1, 1]]",
            "Plus[1, Times[-1, EulerGamma], Times[-1, Power[Plus[1, x], -1]], Times[Plus[1, Times[1/6, Power[Pi, 2]]], Plus[1, x]], Times[Power[Plus[1, x], 2], OrderSize[Plus[1, x]]]]",
        ),
        (
            "Series[Gamma[x], List[x, -1, 1]]",
            "Plus[-1, EulerGamma, Times[-1, Power[Plus[1, x], -1]], Times[Plus[-1, EulerGamma, Times[-1/12, Plus[Times[6, Power[EulerGamma, 2]], Power[Pi, 2]]]], Plus[1, x]], Times[Power[Plus[1, x], 2], OrderSize[Plus[1, x]]]]",
        ),
    ]);
}

/// giac's `series(LambertW(x),x=1,2)`, by Newton's steps.
#[test]
fn lambert_series() {
    check(&[(
        "Series[LambertW[x], List[x, 1, 2]]",
        "Plus[LambertW[1], Times[LambertW[1], Power[Plus[1, LambertW[1]], -1], Plus[-1, x]], Times[Plus[Times[-1, Power[LambertW[1], 3]], Times[-2, Power[LambertW[1], 2]]], Power[Plus[Times[2, Power[LambertW[1], 3]], Times[6, Power[LambertW[1], 2]], Times[6, LambertW[1]], 2], -1], Power[Plus[-1, x], 2]], Times[Power[Plus[-1, x], 3], OrderSize[Plus[-1, x]]]]",
    )]);
}

/// giac's limits at infinity and at 0.
#[test]
fn limits() {
    check(&[
        ("Limit[SinIntegral[x], x, Infinity]", "Times[1/2, Pi]"),
        (
            "Limit[SinIntegral[x], x, Times[-1, Infinity]]",
            "Times[-1/2, Pi]",
        ),
        ("Limit[CosIntegral[x], x, Infinity]", "0"),
        ("Limit[ExpIntegralEi[x], x, Times[-1, Infinity]]", "0"),
        ("Limit[Erf[x], x, Infinity]", "1"),
        ("Limit[Erfc[x], x, Infinity]", "0"),
        (
            "Limit[Plus[PolyGamma[x], Times[-1, Log[x]]], x, Infinity]",
            "0",
        ),
        ("Limit[Times[SinIntegral[x], Power[x, -1]], x, 0]", "1"),
        ("Limit[Times[Gamma[x], x], x, 0]", "1"),
        ("Limit[Times[LambertW[x], Power[x, -1]], x, 0]", "1"),
    ]);
}

/// giac's `int(sin(x)/x,x,0,inf)`, `int(exp(-x^2/2),x,-inf,inf)` and
/// `int(exp(-x)/x,x,1,inf)`.
#[test]
fn definite_integrals() {
    check(&[
        (
            "Integrate[Times[Power[x, -1], Sin[x]], List[x, 0, Infinity]]",
            "Times[1/2, Pi]",
        ),
        (
            "Integrate[Power[E, Times[-1/2, Power[x, 2]]], List[x, Times[-1, Infinity], Infinity]]",
            "Power[Times[2, Pi], 1/2]",
        ),
        (
            "Integrate[Times[Power[E, Times[-1, x]], Power[x, -1]], List[x, 1, Infinity]]",
            "Times[-1, ExpIntegralEi[-1]]",
        ),
    ]);
}
