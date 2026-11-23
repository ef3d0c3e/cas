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

/// giac's `bspline` (CAS-596): knots on `[0, 1]`, `[a, b]` and given,
/// complex, symbolic and matrix control points, the collapse of one column,
/// `piecewise`.
#[test]
fn b_splines() {
    check_exact(&[
        (
            "BSpline[List[1, 2, 3, 4, 5], Equation[x, Range[0, 1]], 2]",
            "List[List[Plus[1, Times[6, x], Times[-9/2, Power[x, 2]]], 0], \
             List[Plus[3/2, Times[3, x]], 0], \
             List[Plus[7/2, Times[-3, x], Times[9/2, Power[x, 2]]], 0]]",
        ),
        (
            "BSpline[List[0, 1, 0], x, 1]",
            "List[List[Times[2, x], 0], List[Plus[2, Times[-2, x]], 0]]",
        ),
        (
            "BSpline[List[List[1, 2], List[3, 4], List[5, 6]], x, 1]",
            "List[List[Plus[1, Times[4, x]], Plus[2, Times[4, x]]], \
             List[Plus[1, Times[4, x]], Plus[2, Times[4, x]]]]",
        ),
        (
            "BSpline[List[1, 2, 3], Equation[x, List[0, 1/2, 1]], 1]",
            "List[List[Plus[1, Times[2, x]], 0], List[Plus[1, Times[2, x]], 0]]",
        ),
        ("BSpline[List[List[1], List[2]], x, 1]", "Plus[1, x]"),
        (
            "BSpline[List[List[1], List[2], List[3]], x, 1]",
            "List[Plus[1, Times[2, x]], Plus[1, Times[2, x]]]",
        ),
        (
            "BSpline[List[1, 2, 3, 4], x, 3, Piecewise]",
            "List[Piecewise[And[GreaterEqual[x, 0], Less[x, 1]], Plus[1, Times[3, x]]], \
             Piecewise[And[GreaterEqual[x, 0], Less[x, 1]], 0]]",
        ),
    ]);
    check(&[
        (
            "BSpline[List[a, b, c], x, 1]",
            "List[List[Plus[a, Times[-2, a, x], Times[2, b, x]], 0], \
             List[Plus[Times[2, b], Times[-1, c], Times[-2, b, x], Times[2, c, x]], 0]]",
        ),
        (
            "BSpline[List[Complex[0, 2], Complex[1, -2], Complex[5, -1], Complex[6, -4], \
             Complex[8, 2], Complex[5, 1]], Equation[x, Range[0, 15]]]",
            "List[List[Plus[Times[3/5, x], Times[3/25, Power[x, 2]], Times[-11/750, Power[x, 3]]], \
             Plus[2, Times[-12/5, x], Times[27/50, Power[x, 2]], Times[-21/500, Power[x, 3]]]], \
             List[Plus[-3, Times[12/5, x], Times[-6/25, Power[x, 2]], Times[7/750, Power[x, 3]]], \
             Plus[-6, Times[12/5, x], Times[-21/50, Power[x, 2]], Times[11/500, Power[x, 3]]]], \
             List[Plus[41, Times[-54/5, x], Times[27/25, Power[x, 2]], Times[-13/375, Power[x, 3]]], \
             Plus[64, Times[-93/5, x], Times[42/25, Power[x, 2]], Times[-6/125, Power[x, 3]]]]]",
        ),
    ]);
    // Clamped: the first piece starts at the first control point and the
    // last ends at the last.
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let pieces = ctx
        .eval(&read(
            "BSpline[List[3, 1, 4, 1, 5, 9], Equation[x, Range[2, 7]], 3]",
        ))
        .unwrap();
    let at = |piece: &Expr, v: i64| ctx.eval(&piece.substitute(&[(x.clone(), Expr::from(v))]));
    let (first, last) = (&pieces.args()[0], pieces.args().last().unwrap());
    assert_eq!(at(first, 2), Ok(read("List[3, 0]")), "{pieces:?}");
    assert_eq!(at(last, 7), Ok(read("List[9, 0]")), "{pieces:?}");
}

/// giac's `fitspline` (CAS-596): exact least squares at points, the options,
/// several series, the variable form; decimals answer balls.
#[test]
fn spline_fits() {
    let six = "List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]]";
    let seven = "List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1], \
                 List[6, 0]]";
    let half = "List[1/2, 5/2]";
    check_exact(&[
        (&format!("SplineFit[{six}, {half}]"), "List[61/126, 1/2]"),
        (&format!("SplineFit[{seven}, List[3]]"), "List[169/553]"),
        (
            &format!("SplineFit[{six}, {half}, Equation[breakpoint, List[0, 2, 5]]]"),
            "List[2053/3276, 61/156]",
        ),
        (
            &format!("SplineFit[{six}, {half}, Equation[Exponent, 2]]"),
            "List[135/334, 1/2]",
        ),
        (
            &format!("SplineFit[{six}, {half}, 3]"),
            "List[229/170, 1/2]",
        ),
        (
            &format!("SplineFit[{six}, {half}, Equation[derive, 1]]"),
            "List[227/378, -109/378]",
        ),
        (
            "SplineFit[List[List[0, 0, 1], List[1, 1, 2], List[2, 0, 3], List[3, 1, 4], \
             List[4, 0, 5], List[5, 1, 6], List[6, 0, 8]], List[1/2, 5/2]]",
            "List[List[1461/2212, 106861/70784], List[783/2212, 247775/70784]]",
        ),
        (
            &format!("SplineFit[{six}, x, Equation[breakpoint, List[0, 1, 5]]]"),
            "List[Plus[Times[138/35, x], Times[-163/35, Power[x, 2]], Times[8/5, Power[x, 3]]], \
             Plus[8/5, Times[-6/7, x], Times[1/7, Power[x, 2]]]]",
        ),
    ]);
    // A decimal: balls near the exact fit.
    let ctx = Context::new();
    let found = ctx
        .eval(&read(&format!("SplineFit[{six}, List[0x1p-1, 5/2]]")))
        .unwrap();
    let gap = Expr::call(Symbol::ABS, [found.args()[0].clone() - read("61/126")]);
    let close = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-12)]);
    assert!(found.args()[0].to_rational().is_none(), "{found:?}");
    assert_eq!(ctx.eval(&close), Ok(Expr::from(Symbol::TRUE)), "{found:?}");
}

/// A fit's control points solve its normal equations `(AᵀA) Q = AᵀY`.
#[test]
fn spline_fit_normal_equations() {
    use super::splines::{basis, least_squares, padded};
    use cas_num::Rational;
    let r = |n: i64| Rational::from(n);
    let s: Vec<Rational> = (0..9).map(r).collect();
    let ys: Vec<Vec<Rational>> = [2, -1, 0, 3, 5, 4, 4, 7, 1]
        .into_iter()
        .map(|y| vec![r(y), r(y * y)])
        .collect();
    let t = padded(&[r(0), r(3), r(5), r(8)], 3);
    let n = 5;
    let q = least_squares(&s, &ys, &t, n, 3).unwrap();
    let a: Vec<Vec<Rational>> = s.iter().map(|u| basis(n, 3, u, &t).unwrap()).collect();
    for i in 0..=n {
        for k in 0..2 {
            let mut left = r(0);
            let mut right = r(0);
            for (row, y) in a.iter().zip(&ys) {
                for (j, qj) in q.iter().enumerate() {
                    left = &left + &(&(&row[i] * &row[j]) * &qj[k]);
                }
                right = &right + &(&row[i] * &y[k]);
            }
            assert_eq!(left, right, "row {i}, series {k}");
        }
    }
}

/// giac's `fitpoly` rows exactly: exact data answer exactly, the
/// coefficients without a variable highest first, a decimal outside the
/// range leaves the answer exact.
#[test]
fn polynomial_fits() {
    let squares = "List[0, 1, 4, 9, 16, 25, 36, 49, 64, 163/2]";
    let tens = "List[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]";
    let eight = "List[0, 1, 2, 3, 4, 5, 6, 7], List[1, 3, 2, 5, 4, 7, 6, 9]";
    let quadratic = "Plus[7/110, Times[-3/40, x], Times[89/88, Power[x, 2]]]";
    check_exact(&[
        (
            "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10], x]",
            "Plus[1, Power[x, 2]]",
        ),
        (
            "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10]]",
            "List[1, 0, 1]",
        ),
        (
            "PolynomialFit[List[List[0, 1], List[1, 2], List[2, 5], List[3, 10]], x]",
            "Plus[1, Power[x, 2]]",
        ),
        (
            "PolynomialFit[List[0, 1, 2, 3, 4], List[1, 2, 5, 10, 35/2], x]",
            "Plus[1, Times[-1/8, x], Times[59/48, Power[x, 2]], Times[-1/8, Power[x, 3]], \
             Times[1/48, Power[x, 4]]]",
        ),
        (
            &format!("PolynomialFit[{eight}, x, Equation[Exponent, 2]]"),
            "Plus[11/8, Times[121/168, x], Times[1/24, Power[x, 2]]]",
        ),
        (&format!("PolynomialFit[{tens}, {squares}, x]"), quadratic),
        (
            &format!("PolynomialFit[{tens}, {squares}, x, Equation[len, 3]]"),
            quadratic,
        ),
        (
            &format!(
                "PolynomialFit[{tens}, List[0, 1, 4, 9, 16, 25, 36, 49, 64, 0x1.46p+6], \
                 Equation[x, Range[2, 7]]]"
            ),
            "Power[x, 2]",
        ),
    ]);
    // Eight points, eight steps: the interpolating polynomial.
    check(&[(
        &format!("PolynomialFit[{eight}, x]"),
        &format!("Expand[InterpolatingPolynomial[{eight}, x]]"),
    )]);
    // A decimal: balls near the exact fit.
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "PolynomialFit[List[0, 1, 2, 3, 4], List[1, 2, 5, 10, 0x1.18p+4], x]",
        ))
        .unwrap();
    let exact = read(
        "Plus[1, Times[-1/8, x], Times[59/48, Power[x, 2]], Times[-1/8, Power[x, 3]], \
         Times[1/48, Power[x, 4]]]",
    );
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [Expr::call(Symbol::EXPAND, [found.clone() - exact])
                .substitute(&[(Expr::symbol("x"), Expr::from(3))])],
        )],
    );
    let close = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-12)]);
    assert!(found.to_rational().is_none(), "{found:?}");
    assert_eq!(ctx.eval(&close), Ok(Expr::from(Symbol::TRUE)), "{found:?}");
}

/// A data fit's coefficients solve the normal equations `(VᵀV) c = Vᵀy` of
/// its degree.
#[test]
fn polynomial_fit_normal_equations() {
    use super::fits::fit_data;
    use cas_num::Rational;
    let r = |n: i64| Rational::from(n);
    let xs: Vec<Rational> = [-3, -1, 0, 2, 3, 5, 8, 9, 12].into_iter().map(r).collect();
    let ys: Vec<Rational> = [2, -1, 0, 3, 5, 4, 4, 7, 1].into_iter().map(r).collect();
    for degree in 0..=4 {
        let c = fit_data(&xs, &ys, &r(0), -degree, 15);
        let n = c.len();
        assert_eq!(n, usize::try_from(degree).unwrap() + 1);
        let power = |x: &Rational, k: usize| (0..k).fold(r(1), |p, _| &p * x);
        for i in 0..n {
            let mut left = r(0);
            let mut right = r(0);
            for (x, y) in xs.iter().zip(&ys) {
                for (j, cj) in c.iter().enumerate() {
                    left = &left + &(&(&power(x, i) * &power(x, j)) * cj);
                }
                right = &right + &(&power(x, i) * y);
            }
            assert_eq!(left, right, "degree {degree}, row {i}");
        }
    }
}

/// A function's fit: certified balls that hold the prototype's
/// coefficients of `exp` on `[0, 1]` (Chebyshev's sum on 64 exact nodes,
/// degree 3, to 28 decimals).
#[test]
fn polynomial_fit_function() {
    use crate::expr::View;
    use cas_num::{Integer, Rational};
    let ctx = Context::new();
    let found = ctx
        .eval(&read("PolynomialFit[Exp[x], Equation[x, Range[0, 1]]]"))
        .unwrap();
    let list = ctx
        .eval(&Expr::call(
            Symbol::COEFFICIENT_LIST,
            [found.clone(), Expr::symbol("x")],
        ))
        .unwrap();
    let expected = [
        "2791073514660980515716119686",
        "4230085218483483250662041596",
        "10161116437138066846558795930",
        "9994825894939007903780512587",
    ];
    assert_eq!(list.args().len(), expected.len(), "{found:?}");
    let ten = Integer::from(10).pow(28);
    for (c, digits) in list.args().iter().zip(expected) {
        let q = Rational::new(digits.parse::<Integer>().unwrap(), ten.clone());
        let View::Ball(b) = c.view() else {
            panic!("{c:?} is not a ball");
        };
        assert!(b.contains_rational(&q), "{c:?} misses {digits}e-28");
    }
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
        // B-splines and fits (CAS-596): too few knots, knots not ascending,
        // a degree 0, decimal knots, a list in a list of points; a point
        // outside the data, derive ≥ p, n + p − 1 ≥ m, breakpoints not
        // covering the data, an unknown option, data not numbers.
        "BSpline[List[0, 1, 0, 1, 0], Equation[x, List[0, 1/2, 1]], 2]",
        "BSpline[List[1, 2, 3], Equation[x, List[0, 1, 1/2]], 1]",
        "BSpline[List[1, 2, 3], x, 0]",
        "BSpline[List[1, 2, 3], Equation[x, Range[0, 0x1p+0]], 1]",
        "BSpline[List[1, List[2], 3], x, 1]",
        "SplineFit[List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[7]]",
        "SplineFit[List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[1], Equation[derive, 3]]",
        "SplineFit[List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[1], 4]",
        "SplineFit[List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[1], Equation[breakpoint, List[1, 2, 5]]]",
        "SplineFit[List[List[0, 0], List[1, 1], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[1], foo]",
        "SplineFit[List[List[0, 0], List[1, y], List[2, 0], List[3, 1], List[4, 0], List[5, 1]], \
         List[1]]",
        // Polynomial fits (CAS-597): a range a > b, a variable without a
        // range, an unknown option, a limit 0, unequal lists, one point,
        // data not numbers; a function without a range, of two symbols,
        // of another variable, a bare variable, not real at the nodes.
        "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10], Equation[x, Range[3, 1]]]",
        "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10], Equation[x, 3]]",
        "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10], x, Equation[foo, 2]]",
        "PolynomialFit[List[0, 1, 2, 3], List[1, 2, 5, 10], x, Equation[Limit, 0]]",
        "PolynomialFit[List[0, 1, 2], List[1, 2, 5, 10], x]",
        "PolynomialFit[List[0], List[1], x]",
        "PolynomialFit[List[0, 1, 2, 3], List[1, y, 5, 10], x]",
        "PolynomialFit[Exp[x], x]",
        "PolynomialFit[Exp[Times[x, y]], Equation[x, Range[0, 1]]]",
        "PolynomialFit[Exp[x], Equation[y, Range[0, 1]]]",
        "PolynomialFit[x, Equation[x, Range[0, 1]]]",
        "PolynomialFit[Log[x], Equation[x, Range[-1, 1]]]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
