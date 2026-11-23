//! Tests of interpolation (CAS-DOC-27 §2–§4): giac's answers compared after
//! `Together` of the difference, the forms of the calls, and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a − b` is 0 once put together, each element of lists apart.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| same(ctx, x, y));
    }
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// `(input, expected)`, both evaluated: the same.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if !same(&ctx, &found, &expected) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `(input, expected)`, both evaluated: the same expression.
fn check_exact(cases: &[(&str, &str)]) {
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

/// giac's `lagrange([1,2,3],[4,8,14])`, `(x-2+4)*(x-1)+4`, in its forms;
/// symbolic values; Newton's nested form.
#[test]
fn interpolating_polynomials() {
    let quadratic = "Plus[2, x, Power[x, 2]]";
    check(&[
        (
            "InterpolatingPolynomial[List[1, 2, 3], List[4, 8, 14]]",
            quadratic,
        ),
        (
            "InterpolatingPolynomial[List[List[1, 2, 3], List[4, 8, 14]]]",
            quadratic,
        ),
        (
            "InterpolatingPolynomial[List[List[1, 4], List[2, 8], List[3, 14]], x]",
            quadratic,
        ),
        (
            "InterpolatingPolynomial[List[1, 2, 3], List[4, 8, 14], t]",
            "Plus[2, t, Power[t, 2]]",
        ),
        // (b-a)*(x-1)+a.
        (
            "InterpolatingPolynomial[List[1, 2], List[a, b]]",
            "Plus[a, Times[Plus[b, Times[-1, a]], Plus[-1, x]]]",
        ),
        // ((-x+2+3/2)*(x-1)-1)*x+1.
        (
            "InterpolatingPolynomial[List[0, 1, 2, 3], List[1, 0, 2, 1]]",
            "Plus[1, Times[x, Plus[-1, Times[Plus[-1, x], Plus[7/2, Times[-1, x]]]]]]",
        ),
        ("InterpolatingPolynomial[List[5], List[7]]", "7"),
    ]);
    let ctx = Context::new();
    let nested = ctx
        .eval(&read(
            "InterpolatingPolynomial[List[1, 2, 3], List[4, 8, 14]]",
        ))
        .unwrap();
    assert_eq!(nested.head(), Some(&Expr::from(Symbol::PLUS)), "{nested:?}");
    assert_ne!(nested, ctx.eval(&read("Plus[2, x, Power[x, 2]]")).unwrap());
}

/// giac's natural splines of degrees 1 to 4.
#[test]
fn splines() {
    check(&[
        (
            "Spline[List[0, 1, 2], List[1, 3, 0], x, 3]",
            "List[Plus[Times[-5/4, Power[x, 3]], Times[13/4, x], 1], Plus[Times[5/4, Power[Plus[x, -1], 3]], Times[-15/4, Power[Plus[x, -1], 2]], Times[-1/2, Plus[x, -1]], 3]]",
        ),
        (
            "Spline[List[0, 1, 2, 3], List[0, 1, 0, 1], x, 1]",
            "List[x, Plus[2, Times[-1, x]], Plus[x, -2]]",
        ),
        (
            "Spline[List[0, 1, 2], List[1, 3, 0], x, 2]",
            "List[Plus[Times[2, Power[x, 2]], 1], Plus[Times[-7, Power[Plus[x, -1], 2]], Times[4, Plus[x, -1]], 3]]",
        ),
        (
            "Spline[List[0, 1, 2, 3], List[1, 0, 2, 1], x, 2]",
            "List[Plus[Times[-1, Power[x, 2]], 1], Plus[Times[4, Power[Plus[x, -1], 2]], Times[-2, Plus[x, -1]]], Plus[Times[-7, Power[Plus[x, -2], 2]], Times[6, Plus[x, -2]], 2]]",
        ),
        (
            "Spline[List[0, 1, 3], List[1, 2, 0], t, 3]",
            "List[Plus[Times[-1/3, Power[t, 3]], Times[4/3, t], 1], Plus[Times[1/6, Power[Plus[t, -1], 3]], Times[-1, Power[Plus[t, -1], 2]], Times[1/3, Plus[t, -1]], 2]]",
        ),
        (
            "Spline[List[0, 1, 2, 3, 4], List[1, 0, 2, 1, 3], x, 4]",
            "List[Plus[Times[67/200, Power[x, 4]], Times[-267/200, x], 1], Plus[Times[-271/200, Power[Plus[x, -1], 4]], Times[67/50, Power[Plus[x, -1], 3]], Times[201/100, Power[Plus[x, -1], 2]], Times[1/200, Plus[x, -1]]], Plus[Times[511/200, Power[Plus[x, -2], 4]], Times[-102/25, Power[Plus[x, -2], 3]], Times[-21/10, Power[Plus[x, -2], 2]], Times[21/8, Plus[x, -2]], 2], Plus[Times[-307/200, Power[Plus[x, -3], 4]], Times[307/50, Power[Plus[x, -3], 3]], Times[99/100, Power[Plus[x, -3], 2]], Times[-719/200, Plus[x, -3]], 1]]",
        ),
    ]);
}

/// giac's `thiele`: `(-11*x+6)/(x-6)` and `x^2+1`.
#[test]
fn thiele() {
    check(&[
        (
            "ThieleInterpolation[List[1, 2, 3], List[1, 4, 9], x]",
            "Times[Plus[6, Times[-11, x]], Power[Plus[-6, x], -1]]",
        ),
        (
            "ThieleInterpolation[List[0, 1, 2, 3], List[1, 2, 5, 10], x]",
            "Plus[1, Power[x, 2]]",
        ),
        (
            "ThieleInterpolation[List[1, 2], List[3, 5], t]",
            "Plus[1, Times[2, t]]",
        ),
        // The fraction ends early: giac's 1/x, and a constant.
        (
            "ThieleInterpolation[List[1, 2, 3, 4], List[1, 1/2, 1/3, 1/4], x]",
            "Power[x, -1]",
        ),
        ("ThieleInterpolation[List[1, 2, 3], List[5, 5, 5], x]", "5"),
    ]);
}

/// giac's `linear_interpolate` (CAS-593), exact where giac answers
/// doubles: its `6.33333333333` is `19/3`; unsorted points, a negative
/// step, the last point sampled.
#[test]
fn linear_interpolation() {
    let points = "List[List[1, 2, 6, 9], List[3, 4, 6, 7]]";
    check_exact(&[
        (
            &format!("LinearInterpolate[{points}, 1, 9, 1]"),
            "List[List[1, 2, 3, 4, 5, 6, 7, 8, 9], \
             List[3, 4, 9/2, 5, 11/2, 6, 19/3, 20/3, 7]]",
        ),
        (
            &format!("LinearInterpolate[{points}, 1, 9, 3]"),
            "List[List[1, 4, 7], List[3, 5, 19/3]]",
        ),
        (
            "LinearInterpolate[List[List[0, 2, 3], List[1, 5, 4]], 0, 3, 1/2]",
            "List[List[0, 1/2, 1, 3/2, 2, 5/2, 3], List[1, 2, 3, 4, 5, 9/2, 4]]",
        ),
        (
            "LinearInterpolate[List[List[3, 0, 2], List[4, 1, 5]], 0, 3, -1]",
            "List[List[0, 1, 2, 3], List[1, 3, 5, 4]]",
        ),
    ]);
}

/// giac's `interval2center` and `center2interval` (CAS-593): giac's
/// `[1.0 .. 3.5,3.5 .. 7.0,7.0 .. 11.0]` exact.
#[test]
fn centres_and_intervals() {
    check_exact(&[
        ("IntervalToCenter[Range[2, 5]]", "7/2"),
        ("IntervalToCenter[3]", "3"),
        ("IntervalToCenter[x]", "x"),
        (
            "IntervalToCenter[List[Range[1, 3], List[Range[5, 7]], 4]]",
            "List[2, List[6], 4]",
        ),
        (
            "List[CenterToInterval[List[2, 5, 9], 1]]",
            "List[Range[1, 7/2], Range[7/2, 7], Range[7, 11]]",
        ),
        (
            "List[CenterToInterval[List[2, 5, 8]]]",
            "List[Range[1/2, 7/2], Range[7/2, 13/2], Range[13/2, 19/2]]",
        ),
        (
            "List[CenterToInterval[List[1, 2, 4]]]",
            "List[Range[1/2, 3/2], Range[3/2, 3], Range[3, 5]]",
        ),
    ]);
}

/// giac's `ratinterp` (CAS-594): its default `d`, `d = 1` and `d = 0`, the
/// matrix form, values at points, a symbolic abscissa.
#[test]
fn rational_interpolation() {
    let (xs, ys) = ("List[1, 3, 5, 8]", "List[2, -1, 3, 4]");
    let cubic = "Plus[Times[-29/168, Power[x, 3]], Times[17/7, Power[x, 2]], \
                 Times[-1507/168, x], 61/7]";
    check(&[
        (&format!("RationalInterpolation[{xs}, {ys}, x]"), cubic),
        (
            "RationalInterpolation[List[List[1, 2], List[3, -1], List[5, 3], List[8, 4]], x]",
            cubic,
        ),
        (
            "RationalInterpolation[List[List[1, 2], List[3, -1], List[5, 3], List[8, 4]]]",
            cubic,
        ),
        (
            &format!("RationalInterpolation[{xs}, {ys}, x, 1]"),
            "Times[Plus[Times[-19, Power[x, 3]], Times[288, Power[x, 2]], Times[-1133, x], 1200], \
             Power[Plus[Times[6, Power[x, 2]], Times[-48, x], 210], -1]]",
        ),
        (
            "RationalInterpolation[List[0, 1, 2, 3, 4], List[1, 2, 1, 2, 1], x, 1]",
            "Times[Plus[Times[-1, Power[x, 4]], Times[8, Power[x, 3]], Times[-19, Power[x, 2]], \
             Times[12, x], 6], Power[Plus[Power[x, 2], Times[-4, x], 6], -1]]",
        ),
        (
            "RationalInterpolation[List[0, 1, 2, 3, 4], List[1, 2, 1, 2, 1], x, 0]",
            "Times[Plus[Times[-1, Power[x, 4]], Times[8, Power[x, 3]], Times[-15, Power[x, 2]], \
             Times[-4, x], 24], Power[Plus[Power[x, 4], Times[-8, Power[x, 3]], \
             Times[25, Power[x, 2]], Times[-36, x], 24], -1]]",
        ),
        (
            "RationalInterpolation[List[0, 1, 2, 3], List[1, 0, 1, 0], t, 1]",
            "Times[Plus[Times[-3, Power[t, 3]], Times[14, Power[t, 2]], Times[-17, t], 6], \
             Power[Plus[Power[t, 2], Times[-3, t], 6], -1]]",
        ),
        (
            "RationalInterpolation[List[1, 2, 3], List[1, 4, 9]]",
            "Power[x, 2]",
        ),
        (
            "RationalInterpolation[List[1, 2, 3], List[1, 4, 9], List[0, 1/2, 4]]",
            "List[0, 1/4, 16]",
        ),
        (
            "RationalInterpolation[List[1, 2, 3], List[1, 4, 9], 5/2, 1]",
            "25/4",
        ),
        (
            "RationalInterpolation[List[List[1, 2], List[3, -1], List[5, 3], List[8, 4]], \
             List[2, 4], 1]",
            "List[-11/23, 10/19]",
        ),
        (
            "RationalInterpolation[List[0, a], List[1, 2], x]",
            "Plus[1, Times[x, Power[a, -1]]]",
        ),
    ]);
}

/// giac's `triginterp` (CAS-594) in its two forms, odd and even `n`.
#[test]
fn trigonometric_interpolation() {
    let eight = "List[11, 10, 17, 24, 32, 26, 23, 19]";
    let answer = "Plus[81/4, \
         Times[21/8, Plus[-2, Times[-1, Sqrt[2]]], Cos[Times[1/12, Pi, x]]], \
         Times[1/8, Plus[-12, Times[-11, Sqrt[2]]], Sin[Times[1/12, Pi, x]]], \
         Times[3/4, Cos[Times[1/6, Pi, x]]], Times[-7/4, Sin[Times[1/6, Pi, x]]], \
         Times[21/8, Plus[-2, Sqrt[2]], Cos[Times[1/4, Pi, x]]], \
         Times[1/8, Plus[12, Times[-11, Sqrt[2]]], Sin[Times[1/4, Pi, x]]], \
         Times[1/2, Cos[Times[1/3, Pi, x]]]]";
    check(&[
        (
            "TrigonometricInterpolation[List[1, 2, 3], Equation[x, Range[0, 2]]]",
            "Plus[2, Times[-1, Cos[Times[2/3, Pi, x]]], \
             Times[-1/3, Sqrt[3], Sin[Times[2/3, Pi, x]]]]",
        ),
        (
            &format!("TrigonometricInterpolation[{eight}, Equation[x, Range[0, 21]]]"),
            answer,
        ),
        (
            &format!("TrigonometricInterpolation[{eight}, 0, 21, x]"),
            answer,
        ),
        (
            "TrigonometricInterpolation[List[1, 0, 1, 0], Equation[x, Range[0, 3]]]",
            "Plus[1/2, Times[1/2, Cos[Times[Pi, x]]]]",
        ),
        (
            "TrigonometricInterpolation[List[0, 1, 0, -1], Equation[x, Range[0, 3]]]",
            "Sin[Times[1/2, Pi, x]]",
        ),
        (
            "TrigonometricInterpolation[List[1, 2], Equation[t, Range[0, 1]]]",
            "Plus[3/2, Times[-1/2, Cos[Times[Pi, t]]]]",
        ),
    ]);
    // A decimal: the result numerically, as giac's evalf.
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "TrigonometricInterpolation[List[0x1p+0, 2, 3], Equation[x, Range[0, 2]]]",
        ))
        .unwrap();
    assert!(
        found.may_contain(Symbol::COS) && !found.may_contain(Symbol::PI),
        "{found:?}"
    );
}

/// giac's `numdiff` and `nDeriv` (CAS-595): orders, several at once, an
/// order past the degree, symbolic values, the difference quotient.
#[test]
fn numerical_derivatives() {
    check(&[
        (
            "NumericalDerivative[List[1, 3, 4], List[3, 0, 1], 2]",
            "-3/2",
        ),
        (
            "NumericalDerivative[List[0, 1, 2, 3], List[0, 1, 8, 27], 3/2, 2]",
            "9",
        ),
        ("NumericalDerivative[List[0, 1, 2], List[0, 1, 4], 1]", "2"),
        (
            "NumericalDerivative[List[0, 1, 2, 3], List[0, 1, 8, 27], 1, 3]",
            "6",
        ),
        ("NumericalDerivative[List[1, 2], List[3, 4], 1, 0]", "3"),
        (
            "NumericalDerivative[List[0, 1/2, 1], List[1, 2, 5], 1/4]",
            "2",
        ),
        (
            "NumericalDerivative[List[0, 1, 2], List[0, 1, 4], 1, 1, 2, 5]",
            "List[2, 2, 0]",
        ),
        (
            "NumericalDerivative[List[0, 1, 2], List[a, b, c], t]",
            "Plus[Times[1/2, Plus[Times[2, t], -3], a], Times[Plus[2, Times[-2, t]], b], \
             Times[1/2, Plus[Times[2, t], -1], c]]",
        ),
        ("NDeriv[Power[x, 2], x]", "Times[2, x]"),
        (
            "NDeriv[Power[x, 3], x, 1/2]",
            "Plus[Times[3, Power[x, 2]], 1/4]",
        ),
        (
            "NDeriv[f[x], x, h]",
            "Times[1/2, Power[h, -1], Plus[f[Plus[x, h]], Times[-1, f[Plus[x, Times[-1, h]]]]]]",
        ),
    ]);
    // A decimal: the result numerically, as giac's evalf.
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "NumericalDerivative[List[0, 1, 2], List[0, 1, 4], 0x1p+0]",
        ))
        .unwrap();
    let gap = Expr::call(Symbol::ABS, [found.clone() - Expr::from(2)]);
    let close = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-10)]);
    assert!(found != Expr::from(2), "{found:?}");
    assert_eq!(ctx.eval(&close), Ok(Expr::from(Symbol::TRUE)), "{found:?}");
}

/// What stays: a repeated abscissa, lists of unequal lengths, abscissas
/// not increasing, Thiele's differences vanishing at some points only, a
/// degree 0.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "InterpolatingPolynomial[List[1, 1], List[2, 3]]",
        "InterpolatingPolynomial[List[1, 2], List[2]]",
        "Spline[List[0, 2, 1], List[1, 2, 3], x, 3]",
        "Spline[List[0, 1, 1, 2], List[1, 2, 3, 4], x, 3]",
        "Spline[List[0, 1], List[1, 2], x, 0]",
        "ThieleInterpolation[List[1, 2, 3], List[1, 2, 1], x]",
        // An unattainable point (crash-9952273f of the fuzz target).
        "ThieleInterpolation[List[-1, -1/2, 2], List[0, -9, -9], x]",
        // Sampling (CAS-593): past the points, a step 0, xmax < xmin, a
        // repeated abscissa, a value not real, too many samples.
        "LinearInterpolate[List[List[0, 2, 3], List[1, 5, 4]], -1, 4, 1]",
        "LinearInterpolate[List[List[0, 2, 3], List[1, 5, 4]], 0, 3, 0]",
        "LinearInterpolate[List[List[0, 2, 3], List[1, 5, 4]], 2, 1, 1]",
        "LinearInterpolate[List[List[0, 0, 3], List[1, 5, 4]], 0, 3, 1]",
        "LinearInterpolate[List[List[0, 2, 3], List[1, y, 4]], 0, 3, 1]",
        "LinearInterpolate[List[List[0, 2, 3], List[1, 5, 4]], 0, 3, 1/100000]",
        "LinearInterpolate[List[List[0, 1], List[0, 1]], 0, 1, 0x1p-40]",
        "CenterToInterval[List[3]]",
        "CenterToInterval[List[1, x]]",
        "CenterToInterval[List[1, 2, 4], 5]",
        "CenterToInterval[List[2, 1]]",
        // Interpolants (CAS-594): d out of range, d = n + 1 (giac's undef),
        // a repeated abscissa, one point, a matrix not of two columns; one
        // value, a = b, a variable not a symbol.
        "RationalInterpolation[List[1, 2, 3], List[1, 4, 9], x, 5]",
        "RationalInterpolation[List[1, 2, 3], List[1, 4, 9], x, 3]",
        "RationalInterpolation[List[1, 2, 3], List[1, 4, 9], x, -1]",
        "RationalInterpolation[List[1, 1, 3], List[1, 4, 9], x]",
        "RationalInterpolation[List[1], List[1], x]",
        "RationalInterpolation[List[List[1, 2, 3], List[4, 5, 6]], x]",
        "TrigonometricInterpolation[List[1], Equation[x, Range[0, 2]]]",
        "TrigonometricInterpolation[List[1, 2, 3], Equation[x, Range[1, 1]]]",
        "TrigonometricInterpolation[List[1, 2, 3], 0, 2, 3]",
        // Numerical derivatives (CAS-595): too few arguments, unequal
        // lengths, an order not a nonnegative integer, a repeated abscissa,
        // exact and decimal; a variable not a symbol, h = 0, no variable.
        "NumericalDerivative[List[1, 2], List[3, 4]]",
        "NumericalDerivative[List[1, 2, 3], List[3, 4], 1]",
        "NumericalDerivative[List[1, 2], List[3, 4], 1, 1/2]",
        "NumericalDerivative[List[1, 2], List[3, 4], 1, -1]",
        "NumericalDerivative[List[1, 1], List[3, 4], 1]",
        "NumericalDerivative[List[0x1p+0, 0x1p+0], List[3, 4], 1]",
        "NDeriv[Power[x, 2], 2, 1/10]",
        "NDeriv[Power[x, 2], x, 0]",
        "NDeriv[Power[x, 2]]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
