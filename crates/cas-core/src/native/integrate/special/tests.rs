//! Tests of the integrals that give `Ei`, `Si`, `Ci` and `erf` (CAS-357):
//! giac's answers, up to the arrangement `Together` sees through, and
//! what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(integrand, antiderivative)`: `Integrate` of the first is the second,
/// once their difference is put together.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (f, expected) in cases {
        let found = ctx
            .eval(&Expr::call(Symbol::INTEGRATE, [read(f), read("x")]))
            .unwrap();
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
        if difference != Expr::from(0) {
            wrong.push(format!("∫ {f}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `int(exp(x)/x,x)`, `int(exp(2*x)/x,x)`, `int(exp(x)/x^2,x)`,
/// `int(exp(x)*(1+x)/x^2,x)`, `int(sinh(x)/x,x)`, and the shifts.
#[test]
fn exponentials() {
    check(&[
        ("Times[Power[E, x], Power[x, -1]]", "ExpIntegralEi[x]"),
        (
            "Times[Power[E, Times[2, x]], Power[x, -1]]",
            "ExpIntegralEi[Times[2, x]]",
        ),
        (
            "Times[Power[E, Plus[1, Times[2, x]]], Power[x, -1]]",
            "Times[E, ExpIntegralEi[Times[2, x]]]",
        ),
        (
            "Times[Power[E, x], Power[x, -2]]",
            "Plus[ExpIntegralEi[x], Times[-1, Power[E, x], Power[x, -1]]]",
        ),
        (
            "Times[Plus[1, x], Power[E, x], Power[x, -2]]",
            "Plus[Times[2, ExpIntegralEi[x]], Times[-1, Power[E, x], Power[x, -1]]]",
        ),
        (
            "Times[Power[E, Times[-1, x]], Power[x, -1]]",
            "ExpIntegralEi[Times[-1, x]]",
        ),
        (
            "Times[Sinh[x], Power[x, -1]]",
            "Times[1/2, Plus[ExpIntegralEi[x], Times[-1, ExpIntegralEi[Times[-1, x]]]]]",
        ),
        (
            "Times[Cosh[x], Power[x, -1]]",
            "Times[1/2, Plus[ExpIntegralEi[x], ExpIntegralEi[Times[-1, x]]]]",
        ),
        (
            "Times[Power[E, x], Power[Plus[-1, x], -1]]",
            "Times[E, ExpIntegralEi[Plus[-1, x]]]",
        ),
        (
            "Times[x, Power[E, x], Power[Plus[1, x], -1]]",
            "Plus[Power[E, x], Times[-1, Power[E, -1], ExpIntegralEi[Plus[1, x]]]]",
        ),
        // Through derivative-divides, u = x²: giac's Ei(x^2)/2.
        (
            "Times[Power[E, Power[x, 2]], Power[x, -1]]",
            "Times[1/2, ExpIntegralEi[Power[x, 2]]]",
        ),
    ]);
}

/// giac's `int(sin(x)/x,x)`, `int(cos(x)/x,x)`, `int(sin(2*x)/x,x)`,
/// `int(sin(x)/x^2,x)`; `cos(2x + 1)/x` by the addition formulas.
#[test]
fn sines_and_cosines() {
    check(&[
        ("Times[Power[x, -1], Sin[x]]", "SinIntegral[x]"),
        ("Times[Power[x, -1], Cos[x]]", "CosIntegral[x]"),
        (
            "Times[Power[x, -1], Sin[Times[2, x]]]",
            "SinIntegral[Times[2, x]]",
        ),
        (
            "Times[Power[x, -2], Sin[x]]",
            "Plus[CosIntegral[x], Times[-1, Power[x, -1], Sin[x]]]",
        ),
        (
            "Times[Power[x, -2], Cos[x]]",
            "Plus[Times[-1, SinIntegral[x]], Times[-1, Power[x, -1], Cos[x]]]",
        ),
        (
            "Times[Power[x, -1], Cos[Plus[1, Times[2, x]]]]",
            "Plus[Times[Cos[1], CosIntegral[Times[2, x]]], Times[-1, Sin[1], SinIntegral[Times[2, x]]]]",
        ),
    ]);
}

/// giac's `int(1/ln(x),x)`, `int(1/ln(2*x),x)` and `int(x/ln(x),x)`.
#[test]
fn logarithms() {
    check(&[
        ("Power[Log[x], -1]", "ExpIntegralEi[Log[x]]"),
        (
            "Power[Log[Times[2, x]], -1]",
            "Times[1/2, ExpIntegralEi[Log[Times[2, x]]]]",
        ),
        (
            "Times[x, Power[Log[x], -1]]",
            "ExpIntegralEi[Times[2, Log[x]]]",
        ),
    ]);
}

/// giac's `int(exp(-x^2),x)`, `int(x^2*exp(-x^2),x)`,
/// `int(exp(-x^2-2*x),x)`, `int(exp(-2*x^2+x),x)` and `int(x^4*exp(-x^2),x)`.
#[test]
fn gaussians() {
    check(&[
        (
            "Power[E, Times[-1, Power[x, 2]]]",
            "Times[1/2, Power[Pi, 1/2], Erf[x]]",
        ),
        (
            "Times[Power[x, 2], Power[E, Times[-1, Power[x, 2]]]]",
            "Plus[Times[-1/2, x, Power[E, Times[-1, Power[x, 2]]]], Times[1/4, Power[Pi, 1/2], Erf[x]]]",
        ),
        (
            "Power[E, Plus[Times[-1, Power[x, 2]], Times[-2, x]]]",
            "Times[1/2, E, Power[Pi, 1/2], Erf[Plus[1, x]]]",
        ),
        (
            "Power[E, Plus[x, Times[-2, Power[x, 2]]]]",
            "Times[1/2, Power[2, -1/2], Power[E, 1/8], Power[Pi, 1/2], Erf[Times[Power[2, 1/2], Plus[-1/4, x]]]]",
        ),
        (
            "Times[Power[x, 4], Power[E, Times[-1, Power[x, 2]]]]",
            "Plus[Times[Plus[Times[-1/2, Power[x, 3]], Times[-3/4, x]], Power[E, Times[-1, Power[x, 2]]]], Times[3/8, Power[Pi, 1/2], Erf[x]]]",
        ),
    ]);
}

/// giac's `int(Ei(x),x)`, `int(Si(x),x)`, `int(Ci(x),x)`,
/// `int(erf(x),x)` and `int(LambertW(x),x)`; `∫ ψ′ = ψ`.
#[test]
fn table() {
    check(&[
        (
            "ExpIntegralEi[x]",
            "Plus[Times[x, ExpIntegralEi[x]], Times[-1, Power[E, x]]]",
        ),
        ("SinIntegral[x]", "Plus[Times[x, SinIntegral[x]], Cos[x]]"),
        (
            "CosIntegral[Times[2, x]]",
            "Times[1/2, Plus[Times[2, x, CosIntegral[Times[2, x]]], Times[-1, Sin[Times[2, x]]]]]",
        ),
        (
            "Erf[x]",
            "Plus[Times[x, Erf[x]], Times[Power[E, Times[-1, Power[x, 2]]], Power[Pi, -1/2]]]",
        ),
        (
            "Erfc[x]",
            "Plus[Times[x, Erfc[x]], Times[-1, Power[E, Times[-1, Power[x, 2]]], Power[Pi, -1/2]]]",
        ),
        ("PolyGamma[x, 1]", "PolyGamma[x]"),
        ("PolyGamma[x, 3]", "PolyGamma[x, 2]"),
        (
            "LambertW[x]",
            "Times[x, Plus[-1, LambertW[x], Power[LambertW[x], -1]]]",
        ),
    ]);
}

/// What needs `Erfi`, a sign the assumptions do not give, or `LogGamma`
/// stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for f in [
        "Power[E, Power[x, 2]]",
        "Power[E, Times[-1, a, Power[x, 2]]]",
        "PolyGamma[x]",
        "Times[Power[E, x], Power[Plus[1, Power[x, 2]], -1]]",
    ] {
        let found = ctx
            .eval(&Expr::call(Symbol::INTEGRATE, [read(f), read("x")]))
            .unwrap();
        assert_eq!(
            found.head(),
            Some(&Expr::from(Symbol::INTEGRATE)),
            "∫ {f}: {found:?}"
        );
    }
}

/// The Bessel and Airy functions (CAS-384): `C_1` of a linear argument,
/// `x^(n+1) C_n(a x)` and `x^(1−n) C_n(a x)`, `x Ai(a x)`, and `Ai′`.
#[test]
fn bessel_functions() {
    check(&[
        ("BesselJ[1, x]", "Times[-1, BesselJ[0, x]]"),
        (
            "BesselI[1, Plus[Times[2, x], 1]]",
            "Times[1/2, BesselI[0, Plus[Times[2, x], 1]]]",
        ),
        ("BesselK[1, x]", "Times[-1, BesselK[0, x]]"),
        ("Times[x, BesselJ[0, x]]", "Times[x, BesselJ[1, x]]"),
        (
            "Times[Power[x, 2], BesselJ[1, Times[2, x]]]",
            "Times[1/2, Power[x, 2], BesselJ[2, Times[2, x]]]",
        ),
        (
            "Times[Power[x, -1], BesselJ[2, x]]",
            "Times[-1, Power[x, -1], BesselJ[1, x]]",
        ),
        ("Times[x, BesselY[0, x]]", "Times[x, BesselY[1, x]]"),
        ("Times[3, x, BesselK[0, x]]", "Times[-3, x, BesselK[1, x]]"),
        (
            "Times[x, BesselI[0, Times[3, x]]]",
            "Times[1/3, x, BesselI[1, Times[3, x]]]",
        ),
        (
            "Times[Power[x, -2], BesselI[3, x]]",
            "Times[Power[x, -2], BesselI[2, x]]",
        ),
        ("Times[x, AiryAi[x]]", "AiryAiPrime[x]"),
        (
            "Times[x, AiryBi[Times[2, x]]]",
            "Times[1/4, AiryBiPrime[Times[2, x]]]",
        ),
        (
            "AiryBiPrime[Times[2, x]]",
            "Times[1/2, AiryBi[Times[2, x]]]",
        ),
    ]);
}

/// The derivative of `Integrate[f, x]` at `x = 7/5` against `f`'s, within
/// `10⁻²⁰`, numbers on both sides: closed forms through `erf`.
fn check_derivatives(integrands: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for f in integrands {
        let g = ctx
            .eval(&Expr::call(Symbol::INTEGRATE, [read(f), read("x")]))
            .unwrap();
        let point = |e: Expr| {
            let at = e.substitute(&[(read("x"), read("7/5"))]);
            ctx.eval(&Expr::call(Symbol::N, [at, Expr::from(30)]))
                .unwrap()
        };
        let derivative = ctx
            .eval(&Expr::call(Symbol::D, [g.clone(), read("x")]))
            .unwrap();
        let gap = Expr::call(Symbol::ABS, [point(derivative) - point(read(f))]);
        let small = Expr::call(Symbol::LESS, [gap, read("Power[10, -20]")]);
        if g.contains(&Expr::from(Symbol::INTEGRATE))
            || ctx.eval(&small).unwrap().as_symbol() != Some(Symbol::TRUE)
        {
            wrong.push(format!("∫ {f}: {g:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Powers times exponentials into `LowerGamma`, and beta integrands into
/// `Beta` (CAS-487): giac's `igamma(a+1,b*x)/b^(a+1)`.
#[test]
fn incomplete_functions() {
    check(&[
        (
            "Times[Power[x, 1/3], Power[E, Times[-2, x]]]",
            "Times[Power[2, -4/3], LowerGamma[4/3, Times[2, x]]]",
        ),
        (
            "Times[Power[x, a], Power[E, Times[-1, x]]]",
            "LowerGamma[Plus[1, a], x]",
        ),
        (
            "Times[Power[x, 1/3], Power[Plus[1, Times[-1, x]], 2/3]]",
            "Beta[4/3, 5/3, x]",
        ),
    ]);
    // Half-integers close into erf; a shift e^(d) comes out.
    check_derivatives(&[
        "Times[Power[x, 3/2], Power[E, Times[-1, x]]]",
        "Times[Power[x, -1/2], Power[E, Plus[1, Times[-3, x]]]]",
        "Times[Power[x, 1/3], Power[E, Times[-2, x]]]",
        "Times[Power[x, 1/3], Power[Plus[1, Times[-1, x]], 2/3]]",
    ]);
}

/// `input` and `expected` are one number, within `10⁻²⁰` by balls.
fn check_numbers(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let numeric = |e: &str| {
            ctx.eval(&Expr::call(Symbol::N, [read(e), Expr::from(30)]))
                .unwrap()
        };
        let (found, want) = (numeric(input), numeric(expected));
        let gap = Expr::call(Symbol::ABS, [found.clone() - want.clone()]);
        let small = Expr::call(Symbol::LESS, [gap, read("Power[10, -20]")]);
        if ctx.eval(&small).unwrap().as_symbol() != Some(Symbol::TRUE) {
            wrong.push(format!("{input}: {found:?}, expected {want:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Limits of the incomplete functions and the definite integrals through
/// them (CAS-488): giac's `Gamma(1/3)/(6*2^(1/3))` and `Gamma(a+1)*
/// Gamma(b+1)/Gamma(a+b+2)`.
#[test]
fn incomplete_limits() {
    let ctx = Context::new();
    for (input, expected) in [
        ("Limit[LowerGamma[a, x], Rule[x, Infinity]]", "Gamma[a]"),
        ("Limit[LowerGamma[1/3, x, 1], Rule[x, Infinity]]", "1"),
        ("Limit[Gamma[1/3, x], Rule[x, Infinity]]", "0"),
        ("Limit[LowerGamma[1/3, x], Rule[x, 0]]", "0"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
    check_numbers(&[
        (
            "Integrate[Times[Power[x, 1/3], Power[E, Times[-2, x]]], List[x, 0, Infinity]]",
            "Times[Gamma[1/3], Power[Times[6, Power[2, 1/3]], -1]]",
        ),
        (
            "Integrate[Times[Power[x, 1/3], Power[Plus[1, Times[-1, x]], 2/3]], List[x, 0, 1]]",
            "Times[Gamma[4/3], Gamma[5/3], Power[Gamma[3], -1]]",
        ),
    ]);
}
