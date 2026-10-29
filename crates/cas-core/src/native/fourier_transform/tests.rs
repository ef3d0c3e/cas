//! Tests of the Fourier transforms (CAS-439, CAS-440): giac's transforms,
//! by value at three points `s` where they are functions, as expressions
//! where they have deltas; unknown functions; what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a − b` at `s = p/q` is within 10⁻²⁰ of 0, its parts apart.
fn agree_at(ctx: &Context, a: &Expr, b: &Expr, p: i64, q: i64) -> bool {
    let s = Expr::symbol("s");
    let point = [(s, Expr::from(p) / Expr::from(q))];
    let gap = (a.clone() - b.clone()).substitute(&point);
    let Ok(n) = ctx.eval(&Expr::call(Symbol::N, [gap, Expr::from(30)])) else {
        return false;
    };
    let tolerance = read("1/100000000000000000000");
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let size = Expr::call(Symbol::ABS, [Expr::call(part, [n.clone()])]);
        ctx.eval(&Expr::call(Symbol::LESS, [size, tolerance.clone()]))
            .is_ok_and(|t| t.as_symbol() == Some(Symbol::TRUE))
    })
}

/// `(f, F)`: `FourierTransform[f, x, s]` is `F` at `s` = 3/7, 5/2, −4/3.
fn check_values(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (f, expected) in cases {
        let input = format!("FourierTransform[{f}, x, s]");
        let found = ctx.eval(&read(&input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found.contains(&Expr::from(Symbol::FOURIER_TRANSFORM))
            || ![(3, 7), (5, 2), (-4, 3)]
                .into_iter()
                .all(|(p, q)| agree_at(&ctx, &found, &expected, p, q))
        {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `(input, expected)`, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Gaussians, two-sided and one-sided exponentials, shifted, scaled and
/// modulated: giac's answers.
#[test]
fn functions() {
    check_values(&[
        // sqrt(pi)*exp(-s^2/4); sqrt(pi)*exp(-1/8*s^2)/sqrt(2).
        (
            "Exp[Times[-1, Power[x, 2]]]",
            "Times[Power[Pi, 1/2], Exp[Times[-1/4, Power[s, 2]]]]",
        ),
        (
            "Exp[Times[-2, Power[x, 2]]]",
            "Times[Power[Pi, 1/2], Power[2, -1/2], Exp[Times[-1/8, Power[s, 2]]]]",
        ),
        // sqrt(pi)*exp(-(s^2+4*i*s)/4).
        (
            "Exp[Times[-1, Power[Plus[x, -1], 2]]]",
            "Times[Power[Pi, 1/2], Exp[Times[-1/4, Plus[Power[s, 2], Times[4, I, s]]]]]",
        ),
        // sqrt(pi)*(-s^2+2)*exp(-s^2/4)/4; (-i)*sqrt(pi)/2*s*exp(-s^2/4).
        (
            "Times[Power[x, 2], Exp[Times[-1, Power[x, 2]]]]",
            "Times[1/4, Power[Pi, 1/2], Plus[2, Times[-1, Power[s, 2]]], Exp[Times[-1/4, Power[s, 2]]]]",
        ),
        (
            "Times[x, Exp[Times[-1, Power[x, 2]]]]",
            "Times[-1/2, I, Power[Pi, 1/2], s, Exp[Times[-1/4, Power[s, 2]]]]",
        ),
        // 2/(s^2+1); 2/(s^2+1)*exp(-i*s); 4/(s^2+4).
        (
            "Exp[Times[-1, Abs[x]]]",
            "Times[2, Power[Plus[1, Power[s, 2]], -1]]",
        ),
        (
            "Exp[Times[-1, Abs[Plus[x, -1]]]]",
            "Times[2, Power[Plus[1, Power[s, 2]], -1], Exp[Times[-1, I, s]]]",
        ),
        (
            "Exp[Times[-1, Abs[Times[2, x]]]]",
            "Times[4, Power[Plus[4, Power[s, 2]], -1]]",
        ),
        // (-i)/(s-i); (-i)/(s-2*i); i/(s+i).
        (
            "Times[HeavisideTheta[x], Exp[Times[-1, x]]]",
            "Times[-1, I, Power[Plus[s, Times[-1, I]], -1]]",
        ),
        (
            "Times[HeavisideTheta[x], Exp[Times[-2, x]]]",
            "Times[-1, I, Power[Plus[s, Times[-2, I]], -1]]",
        ),
        (
            "Times[HeavisideTheta[Times[-1, x]], Exp[x]]",
            "Times[I, Power[Plus[s, I], -1]]",
        ),
        // (-2*i)*exp(-i*s)/s; 2*sin(s/2)/s; i*(exp(i/4*s)-exp(3*i/4*s))/s.
        (
            "Sign[Plus[x, -1]]",
            "Times[-2, I, Exp[Times[-1, I, s]], Power[s, -1]]",
        ),
        ("UnitBox[x]", "Times[2, Sin[Times[1/2, s]], Power[s, -1]]"),
        (
            "UnitBox[Plus[1, Times[2, x]]]",
            "Times[I, Plus[Exp[Times[1/4, I, s]], Times[-1, Exp[Times[3/4, I, s]]]], Power[s, -1]]",
        ),
        // Modulations: giac's sqrt(pi)*(exp(-(s+1)^2/4)+exp(-(s-1)^2/4))/2,
        // i*sqrt(pi)*(exp(-(s+1)^2/4)-exp(-(s-1)^2/4))/2,
        // sqrt(pi)*exp(-(s-1)^2/4), (2*s^2+4)/(s^4+4),
        // (-i)*(s-i)/(s+1-i)/(s-1-i).
        (
            "Times[Exp[Times[-1, Power[x, 2]]], Cos[x]]",
            "Times[1/2, Power[Pi, 1/2], Plus[Exp[Times[-1/4, Power[Plus[s, 1], 2]]], Exp[Times[-1/4, Power[Plus[s, -1], 2]]]]]",
        ),
        (
            "Times[Exp[Times[-1, Power[x, 2]]], Sin[x]]",
            "Times[1/2, I, Power[Pi, 1/2], Plus[Exp[Times[-1/4, Power[Plus[s, 1], 2]]], Times[-1, Exp[Times[-1/4, Power[Plus[s, -1], 2]]]]]]",
        ),
        (
            "Times[Exp[Times[I, x]], Exp[Times[-1, Power[x, 2]]]]",
            "Times[Power[Pi, 1/2], Exp[Times[-1/4, Power[Plus[s, -1], 2]]]]",
        ),
        (
            "Times[Exp[Times[-1, Abs[x]]], Cos[x]]",
            "Times[Plus[4, Times[2, Power[s, 2]]], Power[Plus[4, Power[s, 4]], -1]]",
        ),
        (
            "Times[HeavisideTheta[x], Exp[Times[-1, x]], Cos[x]]",
            "Times[-1, I, Plus[s, Times[-1, I]], Power[Plus[s, 1, Times[-1, I]], -1], Power[Plus[s, -1, Times[-1, I]], -1]]",
        ),
        // The principal value of 1/x, and sin(x)/x by it: giac's
        // (-i)*pi*sign(s) and (sign(s+1)-sign(s-1))*pi/2.
        ("Power[x, -1]", "Times[-1, I, Pi, Sign[s]]"),
        (
            "Times[Sin[x], Power[x, -1]]",
            "Times[1/2, Pi, Plus[Sign[Plus[s, 1]], Times[-1, Sign[Plus[s, -1]]]]]",
        ),
        // sqrt(pi)*exp(-(s^2+4*i*s-4)/4).
        (
            "Exp[Plus[Times[-1, Power[x, 2]], Times[2, x]]]",
            "Times[Power[Pi, 1/2], Exp[Times[-1/4, Plus[Power[s, 2], Times[4, I, s], -4]]]]",
        ),
    ]);
}

/// Rational functions by residues (CAS-441): giac's answers.
#[test]
fn rational_functions() {
    check_values(&[
        // pi*exp(-abs(s)); 1/2*pi*exp(-2*abs(s)).
        (
            "Power[Plus[1, Power[x, 2]], -1]",
            "Times[Pi, Exp[Times[-1, Abs[s]]]]",
        ),
        (
            "Power[Plus[4, Power[x, 2]], -1]",
            "Times[1/2, Pi, Exp[Times[-2, Abs[s]]]]",
        ),
        // Odd: -i π sign(s) e^(−|s|).
        (
            "Times[x, Power[Plus[1, Power[x, 2]], -1]]",
            "Times[-1, I, Pi, Sign[s], Exp[Times[-1, Abs[s]]]]",
        ),
        // Neither: giac's (-exp((1+i)*s)*Heaviside(s)+exp((1+i)*s)
        // +exp((-1+i)*s)*Heaviside(s))*pi.
        (
            "Power[Plus[2, Times[2, x], Power[x, 2]], -1]",
            "Times[Pi, Plus[Times[-1, Exp[Times[Plus[1, I], s]], HeavisideTheta[s]], \
             Exp[Times[Plus[1, I], s]], Times[Exp[Times[Plus[-1, I], s]], HeavisideTheta[s]]]]",
        ),
        // giac's sqrt(2)*exp(-sqrt(2)*abs(s)/2)*(cos(sqrt(2)*abs(s)/2)
        // +sin(sqrt(2)*abs(s)/2))*pi/2.
        (
            "Power[Plus[1, Power[x, 4]], -1]",
            "Times[1/2, Power[2, 1/2], Pi, Exp[Times[-1/2, Power[2, 1/2], Abs[s]]], \
             Plus[Cos[Times[1/2, Power[2, 1/2], Abs[s]]], Sin[Times[1/2, Power[2, 1/2], Abs[s]]]]]",
        ),
        // Double poles, by Residue: (π/2)(1 + |s|) e^(−|s|).
        (
            "Power[Plus[1, Power[x, 2]], -2]",
            "Times[1/2, Pi, Plus[1, Abs[s]], Exp[Times[-1, Abs[s]]]]",
        ),
        // A numerator of degree 1 not odd: (x + 2)/(x² + 1).
        (
            "Times[Plus[x, 2], Power[Plus[1, Power[x, 2]], -1]]",
            "Times[Pi, Plus[2, Times[-1, I, Sign[s]]], Exp[Times[-1, Abs[s]]]]",
        ),
    ]);
    // Real poles stay.
    check(&[(
        "FourierTransform[Power[Plus[-1, Power[x, 2]], -1], x, s]",
        "FourierTransform[Power[Plus[-1, Power[x, 2]], -1], x, s]",
    )]);
}

/// Distributions: constants, deltas, oscillations, the step and powers.
#[test]
fn distributions() {
    check(&[
        ("FourierTransform[1, x, s]", "Times[2, Pi, DiracDelta[s]]"),
        (
            "FourierTransform[DiracDelta[Plus[x, -2]], x, s]",
            "Exp[Times[-2, I, s]]",
        ),
        (
            "FourierTransform[Exp[Times[I, x]], x, s]",
            "Times[2, Pi, DiracDelta[Plus[s, -1]]]",
        ),
        (
            "FourierTransform[Cos[Times[3, x]], x, s]",
            "Plus[Times[Pi, DiracDelta[Plus[s, -3]]], Times[Pi, DiracDelta[Plus[s, 3]]]]",
        ),
        (
            "FourierTransform[HeavisideTheta[x], x, s]",
            "Plus[Times[Pi, DiracDelta[s]], Times[-1, I, Power[s, -1]]]",
        ),
        (
            "FourierTransform[x, x, s]",
            "Times[2, Pi, I, Derivative[1][DiracDelta][s]]",
        ),
        // A delta sifts: g(2) e^(−2is).
        (
            "FourierTransform[Times[Power[x, 2], DiracDelta[Plus[x, -2]]], x, s]",
            "Times[4, Exp[Times[-2, I, s]]]",
        ),
    ]);
}

/// Unknown functions stay, their derivatives by `(i s)ⁿ`; what no rule
/// transforms stays.
#[test]
fn unknown_and_what_stays() {
    check(&[
        (
            "FourierTransform[f[x], x, s]",
            "FourierTransform[f[x], x, s]",
        ),
        (
            "FourierTransform[D[f[x], x], x, s]",
            "Times[I, s, FourierTransform[f[x], x, s]]",
        ),
        // Growing, and not integrable.
        (
            "FourierTransform[Exp[x], x, s]",
            "FourierTransform[Exp[x], x, s]",
        ),
        (
            "FourierTransform[Exp[Power[x, 2]], x, s]",
            "FourierTransform[Exp[Power[x, 2]], x, s]",
        ),
        (
            "FourierTransform[Log[Plus[1, Power[x, 2]]], x, s]",
            "FourierTransform[Log[Plus[1, Power[x, 2]]], x, s]",
        ),
    ]);
}

/// The box's values and derivative (CAS-DOC-41 §1).
#[test]
fn unit_box() {
    check(&[
        ("UnitBox[0]", "1"),
        ("UnitBox[1/2]", "1/2"),
        ("UnitBox[-3/4]", "0"),
        ("UnitBox[2]", "0"),
        ("UnitBox[x]", "UnitBox[x]"),
        (
            "D[UnitBox[x], x]",
            "Plus[DiracDelta[Plus[x, 1/2]], Times[-1, DiracDelta[Plus[x, -1/2]]]]",
        ),
    ]);
}
