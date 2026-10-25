//! Tests of `Integrate` (CAS-194 to CAS-198): giac's antiderivatives,
//! compared as values once evaluated; definite integrals, by value, and the
//! improper ones that stay with their warnings; what stays; the trace and
//! its verification; and a proptest that the derivative of a polynomial's
//! antiderivative is the polynomial.

use proptest::prelude::*;

use crate::trace::{EvalOptions, Step, StepKind, Trace, TraceLevel};
use crate::{Context, Expr, Symbol, Warning};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a` and `b`, evaluated, are the same expression, or rational
/// functions of their generators whose difference is 0.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let difference = Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]);
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// `(input, expected)`, both evaluated: the same value.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx
            .eval(&read(input))
            .unwrap_or_else(|e| panic!("{input}: {e}"));
        let expected = ctx.eval(&read(expected)).unwrap();
        if !same(&ctx, &found, &expected) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answers (CAS-DOC-9 §2), as values.
#[test]
fn antiderivatives() {
    check(&[
        // integrate(x^2,x) = x^3/3.
        ("Integrate[Power[x, 2], x]", "Times[1/3, Power[x, 3]]"),
        // 3*x^3/3+2*x^2/2+x.
        (
            "Integrate[Plus[Times[3, Power[x, 2]], Times[2, x], 1], x]",
            "Plus[x, Power[x, 2], Power[x, 3]]",
        ),
        // ln(abs(x)) and -1/x.
        ("Integrate[Power[x, -1], x]", "Log[Abs[x]]"),
        ("Integrate[Power[x, -2], x]", "Times[-1, Power[x, -1]]"),
        // x^(n+1)/(n+1), even for a symbolic n.
        (
            "Integrate[Power[x, n], x]",
            "Times[Power[x, Plus[n, 1]], Power[Plus[n, 1], -1]]",
        ),
        // 3*x^(1/3)*x/4.
        ("Integrate[Power[x, 1/3], x]", "Times[3/4, Power[x, 4/3]]"),
        // a*x^2/2, and constants.
        ("Integrate[Times[a, x], x]", "Times[1/2, a, Power[x, 2]]"),
        ("Integrate[c, x]", "Times[c, x]"),
        ("Integrate[3, x]", "Times[3, x]"),
        ("Integrate[x, x]", "Times[1/2, Power[x, 2]]"),
        (
            "Integrate[Times[x, Power[Plus[y, 1], -1]], x]",
            "Times[1/2, Power[x, 2], Power[Plus[1, y], -1]]",
        ),
        // x^2/2+integrate(f(x),x), and a product that stays whole.
        (
            "Integrate[Plus[x, f[x]], x]",
            "Plus[Times[1/2, Power[x, 2]], Integrate[f[x], x]]",
        ),
        (
            "Integrate[Plus[x, Times[2, f[x]]], x]",
            "Plus[Times[1/2, Power[x, 2]], Integrate[Times[2, f[x]], x]]",
        ),
        // [1/2*x^2,1/3*x^3].
        (
            "Integrate[List[x, Power[x, 2]], x]",
            "List[Times[1/2, Power[x, 2]], Times[1/3, Power[x, 3]]]",
        ),
        // A power of a linear form, as giac: (x+1)^3/3; polynomials that
        // are not expanded.
        (
            "Integrate[Power[Plus[x, 1], 2], x]",
            "Times[1/3, Power[Plus[1, x], 3]]",
        ),
        (
            "Integrate[Times[x, Plus[x, y]], x]",
            "Plus[Times[1/3, Power[x, 3]], Times[1/2, Power[x, 2], y]]",
        ),
        (
            "Integrate[Times[Plus[x, 1], Plus[x, -1], Power[Plus[y, 2], -1]], x]",
            "Times[Plus[Times[1/3, Power[x, 3]], Times[-1, x]], Power[Plus[y, 2], -1]]",
        ),
        // Gaussian coefficients.
        (
            "Integrate[Times[Complex[0, 1], Power[x, 3]], x]",
            "Times[Complex[0, 1/4], Power[x, 4]]",
        ),
        // Definite integrals of polynomials: integrate(x^2,x,0,1) = 1/3.
        ("Integrate[Power[x, 2], x, 0, 1]", "1/3"),
        ("Integrate[Power[x, 2], List[x, 0, 1]]", "1/3"),
        (
            "Integrate[Plus[x, 1], x, a, b]",
            "Plus[Times[1/2, Power[b, 2]], b, Times[-1/2, Power[a, 2]], Times[-1, a]]",
        ),
        ("Integrate[Power[Plus[x, 1], 2], x, -1, 2]", "9"),
    ]);
}

/// giac's answers at linear arguments (CAS-DOC-9 §4, CAS-195), as values.
#[test]
fn the_table_at_linear_arguments() {
    check(&[
        // exp(x), exp(2*x+1)/2, -exp(-x), exp(a*x)/a, and E^x.
        ("Integrate[Exp[x], x]", "Exp[x]"),
        (
            "Integrate[Exp[Plus[Times[2, x], 1]], x]",
            "Times[1/2, Exp[Plus[1, Times[2, x]]]]",
        ),
        (
            "Integrate[Exp[Times[-1, x]], x]",
            "Times[-1, Exp[Times[-1, x]]]",
        ),
        (
            "Integrate[Exp[Times[a, x]], x]",
            "Times[Power[a, -1], Exp[Times[a, x]]]",
        ),
        ("Integrate[Power[E, x], x]", "Power[E, x]"),
        // 1/ln(2)*2^x, and 3^(2x)/(2 ln 3).
        (
            "Integrate[Power[2, x], x]",
            "Times[Power[2, x], Power[Log[2], -1]]",
        ),
        (
            "Integrate[Power[3, Times[2, x]], x]",
            "Times[1/2, Power[3, Times[2, x]], Power[Log[3], -1]]",
        ),
        // -cos(x), sin(3*x)/3, -cos(a*x+b)/a, 2*sin(x/2).
        ("Integrate[Sin[x], x]", "Times[-1, Cos[x]]"),
        (
            "Integrate[Cos[Times[3, x]], x]",
            "Times[1/3, Sin[Times[3, x]]]",
        ),
        (
            "Integrate[Sin[Plus[Times[a, x], b]], x]",
            "Times[-1, Power[a, -1], Cos[Plus[b, Times[a, x]]]]",
        ),
        (
            "Integrate[Cos[Times[1/2, x]], x]",
            "Times[2, Sin[Times[1/2, x]]]",
        ),
        // -ln(abs(cos(x))) and -ln(abs(cos(2*x)))/2.
        ("Integrate[Tan[x], x]", "Times[-1, Log[Abs[Cos[x]]]]"),
        (
            "Integrate[Tan[Times[2, x]], x]",
            "Times[-1/2, Log[Abs[Cos[Times[2, x]]]]]",
        ),
        // x*ln(x)-x and (-2*x-1+(2*x+1)*ln(2*x+1))/2.
        (
            "Integrate[Log[x], x]",
            "Plus[Times[-1, x], Times[x, Log[x]]]",
        ),
        (
            "Integrate[Log[Plus[Times[2, x], 1]], x]",
            "Times[1/2, Plus[-1, Times[-2, x], Times[Plus[1, Times[2, x]], Log[Plus[1, Times[2, x]]]]]]",
        ),
        // Powers of linear forms: 1/2*ln(abs(2*x+3)), ln(abs(x+a)),
        // -(2*x+1)^-2/4, sqrt(2*x+1)*(2*x+1)/3.
        (
            "Integrate[Power[Plus[Times[2, x], 3], -1], x]",
            "Times[1/2, Log[Abs[Plus[3, Times[2, x]]]]]",
        ),
        (
            "Integrate[Power[Plus[x, a], -1], x]",
            "Log[Abs[Plus[a, x]]]",
        ),
        (
            "Integrate[Power[Plus[Times[2, x], 1], -3], x]",
            "Times[-1/4, Power[Plus[1, Times[2, x]], -2]]",
        ),
        (
            "Integrate[Power[Plus[Times[2, x], 1], 1/2], x]",
            "Times[1/3, Power[Plus[1, Times[2, x]], 3/2]]",
        ),
        // asin(x), 1/2*asin(2*x), asin(1/2*x).
        (
            "Integrate[Power[Plus[1, Times[-1, Power[x, 2]]], -1/2], x]",
            "ArcSin[x]",
        ),
        (
            "Integrate[Power[Plus[1, Times[-4, Power[x, 2]]], -1/2], x]",
            "Times[1/2, ArcSin[Times[2, x]]]",
        ),
        (
            "Integrate[Power[Plus[4, Times[-1, Power[x, 2]]], -1/2], x]",
            "ArcSin[Times[1/2, x]]",
        ),
        // atan(x), 1/2*atan(2*x+1), and atan(2x)/2 for 1/(1+4x^2).
        ("Integrate[Power[Plus[1, Power[x, 2]], -1], x]", "ArcTan[x]"),
        (
            "Integrate[Power[Plus[1, Power[Plus[Times[2, x], 1], 2]], -1], x]",
            "Times[1/2, ArcTan[Plus[1, Times[2, x]]]]",
        ),
        (
            "Integrate[Power[Plus[1, Times[4, Power[x, 2]]], -1], x]",
            "Times[1/2, ArcTan[Times[2, x]]]",
        ),
        // With a constant factor.
        ("Integrate[Times[3, Cos[x]], x]", "Times[3, Sin[x]]"),
    ]);
}

/// giac's answers for rational functions with rational residues (CAS-196),
/// as values; an arctangent where what is left meets the table.
#[test]
fn rational_functions() {
    check(&[
        // 1/2*ln(abs(x-1))-1/2*ln(abs(x+1)).
        (
            "Integrate[Power[Plus[-1, Power[x, 2]], -1], x]",
            "Plus[Times[1/2, Log[Abs[Plus[-1, x]]]], Times[-1/2, Log[Abs[Plus[1, x]]]]]",
        ),
        // 1/2*ln(x^2+1): no real root, no absolute value.
        (
            "Integrate[Times[x, Power[Plus[1, Power[x, 2]], -1]], x]",
            "Times[1/2, Log[Plus[1, Power[x, 2]]]]",
        ),
        // ln(abs(x))-ln(abs(x+1))+1/(x+1).
        (
            "Integrate[Times[Power[x, -1], Power[Plus[1, x], -2]], x]",
            "Plus[Log[Abs[x]], Times[-1, Log[Abs[Plus[1, x]]]], Power[Plus[1, x], -1]]",
        ),
        // -1/x-ln(abs(x))+ln(abs(x+1)).
        (
            "Integrate[Times[Power[x, -2], Power[Plus[1, x], -1]], x]",
            "Plus[Times[-1, Power[x, -1]], Times[-1, Log[Abs[x]]], Log[Abs[Plus[1, x]]]]",
        ),
        // x/(2*(x^2+1))+1/2*atan(x): what Hermite's reduction leaves meets
        // the table.
        (
            "Integrate[Power[Plus[1, Power[x, 2]], -2], x]",
            "Plus[Times[1/2, x, Power[Plus[1, Power[x, 2]], -1]], Times[1/2, ArcTan[x]]]",
        ),
        // 1/4*ln(abs(x-1))-1/4*ln(abs(x+1))-1/2*atan(x).
        (
            "Integrate[Power[Plus[-1, Power[x, 4]], -1], x]",
            "Plus[Times[1/4, Log[Abs[Plus[-1, x]]]], Times[-1/4, Log[Abs[Plus[1, x]]]], Times[-1/2, ArcTan[x]]]",
        ),
        // (x^3+1)/(x-1): x^3/3+x^2/2+x+2*ln(abs(x-1)).
        (
            "Integrate[Times[Plus[1, Power[x, 3]], Power[Plus[-1, x], -1]], x]",
            "Plus[x, Times[1/2, Power[x, 2]], Times[1/3, Power[x, 3]], Times[2, Log[Abs[Plus[-1, x]]]]]",
        ),
        // x^3/(x^2+1): 1/2*(x^2-ln(x^2+1)).
        (
            "Integrate[Times[Power[x, 3], Power[Plus[1, Power[x, 2]], -1]], x]",
            "Plus[Times[1/2, Power[x, 2]], Times[-1/2, Log[Plus[1, Power[x, 2]]]]]",
        ),
        // (3*x^2+1)/(x^3+x): ln(abs(x^3+x)), one residue.
        (
            "Integrate[Times[Plus[1, Times[3, Power[x, 2]]], Power[Plus[x, Power[x, 3]], -1]], x]",
            "Log[Abs[Plus[x, Power[x, 3]]]]",
        ),
        // A residue's argument with fractions, as its primitive integer
        // multiple: 1/(2x - 1)/(x + 1).
        (
            "Integrate[Times[Power[Plus[-1, Times[2, x]], -1], Power[Plus[1, x], -1]], x]",
            "Plus[Times[1/3, Log[Abs[Plus[-1, Times[2, x]]]]], Times[-1/3, Log[Abs[Plus[1, x]]]]]",
        ),
    ]);
}

/// giac's answers for pairs of conjugate quadratic residues (CAS-197), as
/// values: real ones joined in one logarithm, complex ones as arctangents
/// of polynomials, each surd as `Power` writes it.
#[test]
fn quadratic_residues() {
    check(&[
        // 2/sqrt(3)*atan((2*x+1)/sqrt(3)).
        (
            "Integrate[Power[Plus[1, x, Power[x, 2]], -1], x]",
            "Times[2/3, Power[3, 1/2], ArcTan[Times[1/3, Power[3, 1/2], Plus[1, Times[2, x]]]]]",
        ),
        // 1/sqrt(5)*ln(abs((2*x-1-sqrt(5))/(2*x-1+sqrt(5)))).
        (
            "Integrate[Power[Plus[-1, Times[-1, x], Power[x, 2]], -1], x]",
            "Times[1/5, Power[5, 1/2], Log[Abs[Times[Plus[-1, Times[2, x], Times[-1, Power[5, 1/2]]], Power[Plus[-1, Times[2, x], Power[5, 1/2]], -1]]]]]",
        ),
        // sqrt(2)/4*ln(abs(x-sqrt(2)))-sqrt(2)/4*ln(abs(x+sqrt(2))), joined.
        (
            "Integrate[Power[Plus[-2, Power[x, 2]], -1], x]",
            "Times[1/4, Power[2, 1/2], Log[Abs[Times[Plus[x, Times[-1, Power[2, 1/2]]], Power[Plus[x, Power[2, 1/2]], -1]]]]]",
        ),
        // 1/2*atan(x^2).
        (
            "Integrate[Times[x, Power[Plus[1, Power[x, 4]], -1]], x]",
            "Times[1/2, ArcTan[Power[x, 2]]]",
        ),
        // ln(x^2+2*x+5)+1/2*atan((x+1)/2).
        (
            "Integrate[Times[Plus[3, Times[2, x]], Power[Plus[5, Times[2, x], Power[x, 2]], -1]], x]",
            "Plus[Log[Plus[5, Times[2, x], Power[x, 2]]], Times[1/2, ArcTan[Times[1/2, Plus[1, x]]]]]",
        ),
        // Bronstein's example [Bro05 §2.8]: R = 45796 (4t^2 + 1)^3,
        // atan(x)+atan(x^3)+atan((x^5-3*x^3+x)/2).
        (
            "Integrate[Times[Plus[6, Times[-3, Power[x, 2]], Power[x, 4]], Power[Plus[4, Times[5, Power[x, 2]], Times[-5, Power[x, 4]], Power[x, 6]], -1]], x]",
            "Plus[ArcTan[x], ArcTan[Power[x, 3]], ArcTan[Times[1/2, Plus[x, Times[-3, Power[x, 3]], Power[x, 5]]]]]",
        ),
        // Real residues ±√2 at the complex roots of x^4 - 2x^2 + 9: no
        // absolute value, as giac's sqrt(2)*ln(x^2-2*sqrt(2)*x+3)-....
        (
            "Integrate[Times[Plus[-24, Times[8, Power[x, 2]]], Power[Plus[9, Times[-2, Power[x, 2]], Power[x, 4]], -1]], x]",
            "Times[Power[2, 1/2], Log[Times[Plus[3, Times[-2, Power[2, 1/2], x], Power[x, 2]], Power[Plus[3, Times[2, Power[2, 1/2], x], Power[x, 2]], -1]]]]",
        ),
    ]);
}

/// A rational function whose residues are not all rational or quadratic:
/// those integrate, the rest stays under `Integrate`, and the derivative is
/// the integrand (CAS-196, CAS-197).
#[test]
fn partial_rational_functions() {
    let ctx = Context::new();
    for input in [
        // 1/((x - 1)(x^3 - 2)): -ln|x - 1| and the rest over x^3 - 2.
        "Integrate[Power[Times[Plus[-1, x], Plus[-2, Power[x, 3]]], -1], x]",
        // A pair over x^2 + 1, residues (1 ± 2i)/10, and the rest.
        "Integrate[Power[Times[Plus[1, Power[x, 2]], Plus[-2, Power[x, 3]]], -1], x]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert!(
            found.contains(&Expr::from(Symbol::LOG))
                || found.contains(&Expr::from(Symbol::ARC_TAN)),
            "{input}: {found:?}"
        );
        assert!(
            found.contains(&Expr::from(Symbol::INTEGRATE)),
            "{input}: {found:?}"
        );
        let derivative = ctx
            .eval(&Expr::call(Symbol::D, [found, Expr::symbol("x")]))
            .unwrap();
        assert!(
            same(&ctx, &derivative, &e.args()[0]),
            "{input}: {derivative:?}"
        );
    }
}

/// Whether `a` and `b` are the same, or numbers whose difference `N` takes
/// for 0, as balls: `Cos[0] - Cos[Pi]` stays as it is, and is 2.
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if same(ctx, a, b) {
        return true;
    }
    let difference = Expr::call(Symbol::N, [a.clone() - b.clone()]);
    let zero = ctx
        .eval(&difference)
        .map(|d| Expr::call(Symbol::EQUAL, [d, Expr::from(0)]));
    zero.and_then(|z| ctx.eval(&z))
        .is_ok_and(|t| t.as_symbol() == Some(Symbol::TRUE))
}

/// `(input, value)` in `ctx`: a definite integral that integrates, to the
/// same value, without a warning.
fn check_values(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        let found = found.unwrap_or_else(|e| panic!("{input}: {e}"));
        let expected = ctx.eval(&read(expected)).unwrap();
        if found.contains(&Expr::from(Symbol::INTEGRATE))
            || !same_value(ctx, &found, &expected)
            || !report.warnings.is_empty()
        {
            wrong.push(format!(
                "{input}: {found:?} {:?}, expected {expected:?}",
                report.warnings
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's definite integrals (CAS-198): `F(b) − F(a)` where the integrand
/// is continuous between the bounds, rational or real numbers, symbols, or
/// a complex bound for an integrand continuous on the plane.
#[test]
fn definite_integrals() {
    let ctx = Context::new();
    check_values(
        &ctx,
        &[
            // ln(2), pi/4, pi/(3*sqrt(3)).
            ("Integrate[Power[x, -1], x, 1, 2]", "Log[2]"),
            (
                "Integrate[Power[Plus[1, Power[x, 2]], -1], x, 0, 1]",
                "Times[1/4, Pi]",
            ),
            (
                "Integrate[Power[Plus[1, Power[x, 2]], -1], List[x, 0, 1]]",
                "Times[1/4, Pi]",
            ),
            (
                "Integrate[Power[Plus[1, x, Power[x, 2]], -1], x, 0, 1]",
                "Times[1/9, Power[3, 1/2], Pi]",
            ),
            // 1/4*(-sqrt(2)*ln(sqrt(2)+1)+sqrt(2)*ln(sqrt(2)-1)): the pole
            // at sqrt(2) is past the bound.
            (
                "Integrate[Power[Plus[-2, Power[x, 2]], -1], x, 0, 1]",
                "Times[1/4, Power[2, 1/2], Plus[Log[Plus[-1, Power[2, 1/2]]], Times[-1, Log[Plus[1, Power[2, 1/2]]]]]]",
            ),
            // 1/18*sqrt(3)*pi+1/18*(sqrt(3)*pi+6*ln(2)).
            (
                "Integrate[Power[Plus[1, Power[x, 3]], -1], x, 0, 1]",
                "Plus[Times[1/9, Power[3, 1/2], Pi], Times[1/3, Log[2]]]",
            ),
            // The bounds the other way: -ln(2).
            ("Integrate[Power[x, -1], x, 2, 1]", "Times[-1, Log[2]]"),
            // Real numbers as bounds, and the table.
            ("Integrate[Sin[x], x, 0, Pi]", "2"),
            ("Integrate[Cos[x], x, 0, Times[1/2, Pi]]", "1"),
            ("Integrate[Exp[x], x, 0, 1]", "Plus[-1, E]"),
            (
                "Integrate[Power[x, 2], x, 0, Power[2, 1/2]]",
                "Times[2/3, Power[2, 1/2]]",
            ),
            ("Integrate[Power[x, -1], x, 1, E]", "1"),
            // sqrt(x) is continuous at 0, and so is its antiderivative.
            ("Integrate[Power[x, 1/2], x, 0, 1]", "2/3"),
            ("Integrate[Log[x], x, 1, 2]", "Plus[-1, Times[2, Log[2]]]"),
            ("Integrate[Tan[x], x, 0, 1]", "Times[-1, Log[Cos[1]]]"),
            (
                "Integrate[Power[Plus[1, Times[-1, Power[x, 2]]], -1/2], x, 0, 1/2]",
                "Times[1/6, Pi]",
            ),
            // A complex bound, for a polynomial.
            ("Integrate[x, x, 0, I]", "-1/2"),
            // Symbolic bounds, where the integrand is continuous on the line.
            (
                "Integrate[x, x, a, b]",
                "Plus[Times[1/2, Power[b, 2]], Times[-1/2, Power[a, 2]]]",
            ),
            (
                "Integrate[Times[x, Power[Plus[1, Power[x, 2]], -1]], x, a, b]",
                "Plus[Times[1/2, Log[Plus[1, Power[b, 2]]]], Times[-1/2, Log[Plus[1, Power[a, 2]]]]]",
            ),
            (
                "Integrate[Exp[Times[2, x]], x, a, b]",
                "Plus[Times[1/2, Exp[Times[2, b]]], Times[-1/2, Exp[Times[2, a]]]]",
            ),
        ],
    );
    // An assumption bounds a symbolic bound: 1 - 1/a.
    let assumed = ctx.assuming(read("Greater[a, 1]"));
    check_values(
        &assumed,
        &[(
            "Integrate[Power[x, -2], x, 1, a]",
            "Plus[1, Times[-1, Power[a, -1]]]",
        )],
    );
}

/// Improper integrals and those not shown proper stay, each with its
/// warning, where giac answers `+infinity`, a principal value or `undef`
/// (CAS-198).
#[test]
fn improper_integrals_stay() {
    let ctx = Context::new();
    let unbounded = "an improper integral: the integrand is unbounded between the bounds";
    let unknown = "the integrand is not shown continuous between the bounds";
    for (input, message) in [
        ("Integrate[Power[x, -2], x, -1, 1]", unbounded),
        ("Integrate[Power[x, -1], x, -1, 1]", unbounded),
        ("Integrate[Power[Plus[-1, x], -1], x, 0, 2]", unbounded),
        (
            "Integrate[Power[Plus[-2, Power[x, 2]], -1], x, 0, 2]",
            unbounded,
        ),
        // Convergent, but at a singular point: for limits.
        ("Integrate[Power[x, -1/2], x, 0, 1]", unbounded),
        ("Integrate[Log[x], x, 0, 1]", unbounded),
        // tan has a pole at pi/2.
        ("Integrate[Tan[x], x, 0, 2]", unbounded),
        // Where the bounds may lie on either side of the pole.
        ("Integrate[Power[x, -1], x, a, b]", unknown),
        ("Integrate[Power[x, -2], x, 1, a]", unknown),
        // A complex bound and a pole.
        ("Integrate[Power[x, -1], x, 1, I]", unknown),
        // A pole at the bound sqrt(2), which balls cannot tell apart.
        (
            "Integrate[Power[Plus[-2, Power[x, 2]], -1], x, 0, Power[2, 1/2]]",
            unknown,
        ),
        (
            "Integrate[Exp[x], x, 0, Infinity]",
            "an infinite bound, which waits for limits",
        ),
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        let found = found.unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
        assert_eq!(found.args().len(), e.args().len(), "{input}: {found:?}");
        let messages: Vec<String> = report.warnings.iter().map(Warning::to_string).collect();
        assert!(
            messages.iter().any(|m| m.starts_with(message)),
            "{input}: {messages:?}"
        );
    }
}

/// What the methods of CAS-194 to CAS-197 leave as it is.
#[test]
fn integrals_that_stay() {
    let ctx = Context::new();
    for input in [
        // Arguments that are not linear, products the heuristics of CAS-199
        // take.
        "Integrate[Sin[Power[x, 2]], x]",
        "Integrate[Exp[Power[x, 2]], x]",
        "Integrate[Times[x, Exp[x]], x]",
        // Residues of degree 3 and 4.
        "Integrate[Power[Plus[-2, Power[x, 3]], -1], x]",
        "Integrate[Power[Plus[1, Power[x, 4]], -1], x]",
        // Parameters in a denominator with x (CAS-DOC-9 §3.4).
        "Integrate[Power[Plus[a, Power[x, 2]], -1], x]",
        "Integrate[f[x], x]",
        // giac's integrate(2*f(x),x) stays whole.
        "Integrate[Times[2, f[x]], x]",
        "Integrate[Power[x, x], x]",
        // Not a symbol, or the wrong arguments.
        "Integrate[x, 2]",
        "Integrate[x, x, 1]",
        "Integrate[x]",
        "Integrate[x, List[x, 1]]",
    ] {
        // The arguments evaluate, the call stays.
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
        assert_eq!(found.args().len(), e.args().len(), "{input}: {found:?}");
    }
}

/// One `"Integrate"` step a call, which `Trace::verify` checks through the
/// derivative of the antiderivative, as it checks the threading of lists.
#[test]
fn traces_verify() {
    let ctx = Context::new();
    for input in [
        "Integrate[Plus[Times[3, Power[x, 2]], Times[2, x], 1], x]",
        "Integrate[Power[x, -1], x]",
        "Integrate[Power[x, 1/2], x]",
        "Integrate[Tan[Times[2, x]], x]",
        "Integrate[Log[Plus[x, 1]], x]",
        "Integrate[Power[Plus[4, Times[-1, Power[x, 2]]], -1/2], x]",
        "Integrate[Power[3, x], x]",
        "Integrate[Power[Plus[-1, Power[x, 2]], -1], x]",
        "Integrate[Power[Plus[1, Power[x, 2]], -2], x]",
        "Integrate[Power[Plus[1, x, Power[x, 2]], -1], x]",
        "Integrate[Power[Plus[-1, Times[-1, x], Power[x, 2]], -1], x]",
        "Integrate[Power[Plus[1, Power[x, 3]], -1], x]",
        "Integrate[Times[a, Power[Plus[x, y], 2]], x]",
        "Integrate[Plus[x, f[x]], x]",
        "Integrate[List[x, Power[x, 2]], x]",
    ] {
        let options = EvalOptions {
            trace: TraceLevel::Rules,
        };
        let (result, report) = ctx.eval_with(&read(input), &options);
        result.unwrap();
        let verification = report.trace.verify(&ctx);
        assert!(
            verification.wrong.is_none(),
            "{input}: {:?}",
            verification.wrong
        );
        assert!(verification.passed > 0, "{input}");
        let steps: Vec<&Step> = report
            .trace
            .steps
            .iter()
            .filter(|s| matches!(s.kind, StepKind::Algorithm(_)))
            .collect();
        assert!(!steps.is_empty(), "{input}");
        assert_eq!(
            steps[0].before.head(),
            Some(&Expr::from(Symbol::INTEGRATE)),
            "{input}"
        );
    }
}

/// A wrong antiderivative fails the verification: the step's claim is its
/// derivative.
#[test]
fn a_wrong_antiderivative_is_caught() {
    let ctx = Context::new();
    let step = |after: &str| Step {
        id: 0,
        kind: StepKind::Algorithm("Integrate".into()),
        position: Vec::new(),
        before: read("Integrate[Power[x, 2], x]"),
        after: read(after),
        parent: None,
        children: Vec::new(),
    };
    let wrong = Trace {
        steps: vec![step("Times[1/2, Power[x, 3]]")],
    };
    assert!(wrong.verify(&ctx).wrong.is_some());
    // A constant more is right.
    let right = Trace {
        steps: vec![step("Plus[7, Times[1/3, Power[x, 3]]]")],
    };
    let verification = right.verify(&ctx);
    assert!(verification.wrong.is_none() && verification.passed == 1);
}

/// Terms in `x` and `y`, small exponents and coefficients.
fn terms() -> impl Strategy<Value = Vec<((u32, u32), i64)>> {
    prop::collection::vec(((0u32..5, 0u32..3), -9i64..9), 1..6)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// `D[Integrate[p, x], x]` is `p` for polynomials `p` in `x` and `y`.
    #[test]
    fn derivatives_of_antiderivatives(terms in terms()) {
        let ctx = Context::new();
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let p: Expr = terms
            .iter()
            .map(|&((i, j), c)| Expr::from(c) * x.pow(i64::from(i)) * y.pow(i64::from(j)))
            .sum();
        let p = ctx.eval(&p).unwrap();
        let integral = ctx.eval(&Expr::call(Symbol::INTEGRATE, [p.clone(), x.clone()])).unwrap();
        prop_assert!(!integral.contains(&Expr::from(Symbol::INTEGRATE)), "{integral:?}");
        let derivative = ctx.eval(&Expr::call(Symbol::D, [integral, x])).unwrap();
        prop_assert!(same(&ctx, &derivative, &p), "{derivative:?} against {p:?}");
    }
}
