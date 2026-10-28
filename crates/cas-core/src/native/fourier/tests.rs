//! Tests of the Fourier coefficients (CAS-DOC-25 §7, CAS-319): giac's
//! answers, compared as values at sample points of `n`, `T` and `a`; the
//! closed forms for a symbolic `n`; what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a` and `b` agree to `10⁻³⁰` as numbers, with `n`, `T` and `a`
/// at `n`, `5/2` and `1/3`.
fn agree(ctx: &Context, a: &Expr, b: &Expr, n: i64) -> bool {
    let rules = [
        (Expr::symbol("n"), Expr::from(n)),
        (Expr::symbol("T"), read("5/2")),
        (Expr::symbol("a"), read("1/3")),
    ];
    let difference = Expr::call(
        Symbol::N,
        [a.substitute(&rules) - b.substitute(&rules), Expr::from(40)],
    );
    let small = Expr::call(
        Symbol::LESS,
        [
            Expr::call(Symbol::ABS, [difference]),
            Expr::from(10).pow(-30),
        ],
    );
    ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// `(input, expected)`: the same values at `n = 2 … 5`.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = read(expected);
        if found.contains(&Expr::symbol("FourierCosCoefficient"))
            || found.contains(&Expr::symbol("FourierSinCoefficient"))
            || found.contains(&Expr::symbol("FourierCoefficient"))
            || !(2..=5).all(|n| agree(&ctx, &found, &expected, n))
        {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac 1.9's `fourier_an`, `fourier_bn` and `fourier_cn` of polynomials:
/// a symbolic `n`, `n = 0`, a given `n`, the default start 0.
#[test]
fn polynomials() {
    check(&[
        // 4*(-1)^n/n^2, -2*(-1)^n/n, i*(-1)^n/n.
        (
            "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[4, Power[-1, n], Power[n, -2]]",
        ),
        (
            "FourierSinCoefficient[x, x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[-2, Power[-1, n], Power[n, -1]]",
        ),
        (
            "FourierCoefficient[x, x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[I, Power[-1, n], Power[n, -1]]",
        ),
        // pi^2/3, -4/9, 4/n^2 over [0, 2 pi], -2/n.
        (
            "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], 0, Times[-1, Pi]]",
            "Times[1/3, Power[Pi, 2]]",
        ),
        (
            "FourierCoefficient[Power[x, 2], x, Times[2, Pi], 0, Times[-1, Pi]]",
            "Times[1/3, Power[Pi, 2]]",
        ),
        (
            "FourierSinCoefficient[Power[x, 2], x, Times[2, Pi], 0, Times[-1, Pi]]",
            "0",
        ),
        (
            "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], 3, Times[-1, Pi]]",
            "-4/9",
        ),
        (
            "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], -2, Times[-1, Pi]]",
            "1",
        ),
        (
            "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], n]",
            "Times[4, Power[n, -2]]",
        ),
        (
            "FourierSinCoefficient[x, x, Times[2, Pi], n]",
            "Times[-2, Power[n, -1]]",
        ),
        // (-2*n^2*pi^2*(-1)^n+12*(-1)^n)/n^3.
        (
            "FourierSinCoefficient[Power[x, 3], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[Power[-1, n], Plus[12, Times[-2, Power[n, 2], Power[Pi, 2]]], Power[n, -3]]",
        ),
        // 8/pi for n = 1/2.
        (
            "FourierSinCoefficient[x, x, Times[2, Pi], 1/2, Times[-1, Pi]]",
            "Times[8, Power[Pi, -1]]",
        ),
    ]);
}

/// giac 1.9's coefficients of `abs`, `sign` and `|sin x|` (CAS-316,
/// CAS-317), of exponentials, at resonance (CAS-318), over other periods
/// and starts.
#[test]
fn other_functions() {
    check(&[
        // (2*(-1)^n-2)/(n^2*pi), (-2*(-1)^n+2)/(n*pi).
        (
            "FourierCosCoefficient[Abs[x], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[2, Plus[-1, Power[-1, n]], Power[n, -2], Power[Pi, -1]]",
        ),
        (
            "FourierSinCoefficient[Sign[x], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[2, Plus[1, Times[-1, Power[-1, n]]], Power[n, -1], Power[Pi, -1]]",
        ),
        (
            "FourierCoefficient[Abs[x], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[Plus[-1, Power[-1, n]], Power[n, -2], Power[Pi, -1]]",
        ),
        // ((-1)^n*exp(pi)^2-(-1)^n)/(n^2*pi*exp(pi)+pi*exp(pi)).
        (
            "FourierCosCoefficient[Exp[x], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[Power[-1, n], Plus[Exp[Pi], Times[-1, Exp[Times[-1, Pi]]]], Power[Times[Pi, Plus[1, Power[n, 2]]], -1]]",
        ),
        // (exp(2)-1)/(n^2*pi^2*exp(2)+exp(2)).
        (
            "FourierCosCoefficient[Exp[Times[-1, x]], x, 2, n]",
            "Times[Plus[1, Times[-1, Exp[-2]]], Power[Plus[1, Times[Power[n, 2], Power[Pi, 2]]], -1]]",
        ),
        // T^2/(n^2*pi^2), -1/(n*pi), i/(2*n*pi).
        (
            "FourierCosCoefficient[Power[x, 2], x, T, n]",
            "Times[Power[T, 2], Power[n, -2], Power[Pi, -2]]",
        ),
        (
            "FourierSinCoefficient[x, x, 1, n]",
            "Times[-1, Power[Times[n, Pi], -1]]",
        ),
        (
            "FourierCoefficient[x, x, 1, n]",
            "Times[I, Power[Times[2, n, Pi], -1]]",
        ),
        // giac leaves sin((2*T*n*pi+2*a*n*pi)/T): -T cos(2 pi n a/T)/(pi n).
        (
            "FourierSinCoefficient[x, x, T, n, a]",
            "Times[-1, T, Cos[Times[2, Pi, n, a, Power[T, -1]]], Power[Times[Pi, n], -1]]",
        ),
        // Given n: 1/2, and resonance (CAS-318).
        (
            "FourierCosCoefficient[Power[Cos[x], 2], x, Times[2, Pi], 2, Times[-1, Pi]]",
            "1/2",
        ),
        (
            "FourierCoefficient[Cos[x], x, Times[2, Pi], 1, Times[-1, Pi]]",
            "1/2",
        ),
        (
            "FourierCoefficient[Cos[x], x, Times[2, Pi], -1, Times[-1, Pi]]",
            "1/2",
        ),
        (
            "FourierCoefficient[Exp[Times[I, x]], x, Times[2, Pi], 1, Times[-1, Pi]]",
            "1",
        ),
        // The rectified sine (CAS-317): 2 (1 + (-1)^n)/(pi (1 - n^2)).
        (
            "FourierCosCoefficient[Abs[Sin[x]], x, Times[2, Pi], n, Times[-1, Pi]]",
            "Times[2, Plus[1, Power[-1, n]], Power[Times[Pi, Plus[1, Times[-1, Power[n, 2]]]], -1]]",
        ),
    ]);
}

/// A symbolic `n` gives closed forms in `(−1)ⁿ`, without `sin(nπ)` or
/// `cos(nπ)` (CAS-313).
#[test]
fn closed_forms() {
    let ctx = Context::new();
    for input in [
        "FourierCosCoefficient[Power[x, 2], x, Times[2, Pi], n, Times[-1, Pi]]",
        "FourierSinCoefficient[x, x, T, n]",
        "FourierCoefficient[Abs[x], x, Times[2, Pi], n, Times[-1, Pi]]",
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        let trig = [Symbol::SIN, Symbol::COS].map(|s| found.contains(&Expr::from(s)));
        assert_eq!(trig, [false, false], "{input}: {found:?}");
    }
}

/// What stays: an unknown function, an integral that stays, a period 0, a
/// variable that is no symbol, a period in the variable.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "FourierCosCoefficient[f[x], x, Times[2, Pi], n]",
        "FourierSinCoefficient[Exp[Power[x, 2]], x, Times[2, Pi], n]",
        "FourierCoefficient[x, x, 0, n]",
        "FourierCosCoefficient[x, Plus[1, x], Times[2, Pi], n]",
        "FourierCosCoefficient[x, x, x, n]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
