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
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
