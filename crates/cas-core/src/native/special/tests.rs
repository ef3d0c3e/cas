//! Tests of the special values (CAS-200): giac's table, the shifts and
//! signs, what stays.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
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

/// giac 1.9's `sin`, `cos` and `tan` at multiples of `pi/6` and `pi/4`.
#[test]
fn sines_cosines_and_tangents() {
    check(&[
        ("Sin[0]", "0"),
        ("Sin[Times[1/6, Pi]]", "1/2"),
        ("Sin[Times[1/4, Pi]]", "Times[1/2, Power[2, 1/2]]"),
        ("Sin[Times[1/3, Pi]]", "Times[1/2, Power[3, 1/2]]"),
        ("Sin[Times[1/2, Pi]]", "1"),
        ("Sin[Times[2/3, Pi]]", "Times[1/2, Power[3, 1/2]]"),
        ("Sin[Times[3/4, Pi]]", "Times[1/2, Power[2, 1/2]]"),
        ("Sin[Times[5/6, Pi]]", "1/2"),
        ("Sin[Pi]", "0"),
        ("Sin[Times[7/6, Pi]]", "-1/2"),
        ("Sin[Times[-1/3, Pi]]", "Times[-1/2, Power[3, 1/2]]"),
        ("Sin[Times[13/6, Pi]]", "1/2"),
        ("Cos[0]", "1"),
        ("Cos[Times[1/6, Pi]]", "Times[1/2, Power[3, 1/2]]"),
        ("Cos[Times[1/4, Pi]]", "Times[1/2, Power[2, 1/2]]"),
        ("Cos[Times[1/3, Pi]]", "1/2"),
        ("Cos[Times[1/2, Pi]]", "0"),
        ("Cos[Times[2/3, Pi]]", "-1/2"),
        ("Cos[Pi]", "-1"),
        ("Cos[Times[-1/4, Pi]]", "Times[1/2, Power[2, 1/2]]"),
        ("Cos[Times[5/3, Pi]]", "1/2"),
        ("Tan[0]", "0"),
        ("Tan[Times[1/6, Pi]]", "Times[1/3, Power[3, 1/2]]"),
        ("Tan[Times[1/4, Pi]]", "1"),
        ("Tan[Times[1/3, Pi]]", "Power[3, 1/2]"),
        ("Tan[Times[1/2, Pi]]", "ComplexInfinity"),
        ("Tan[Times[-1/4, Pi]]", "-1"),
        ("Tan[Times[2/3, Pi]]", "Times[-1, Power[3, 1/2]]"),
        ("Tan[Pi]", "0"),
    ]);
}

/// giac's `sin(x+pi)`, `cos(x+2*pi)`, `sin(pi-x)`, `sin(x+pi/2)`,
/// `cos(pi/2-x)`, and the signs of `sin(-x)`, `cos(-x)`, `tan(-2*x)`,
/// `sin(-1)`, `asin(-x)` and `atan(-x)`.
#[test]
fn shifts_and_signs() {
    check(&[
        ("Sin[Plus[x, Pi]]", "Times[-1, Sin[x]]"),
        ("Cos[Plus[x, Times[2, Pi]]]", "Cos[x]"),
        ("Sin[Plus[Pi, Times[-1, x]]]", "Sin[x]"),
        ("Sin[Plus[x, Times[1/2, Pi]]]", "Cos[x]"),
        ("Cos[Plus[Times[1/2, Pi], Times[-1, x]]]", "Sin[x]"),
        ("Tan[Plus[x, Pi]]", "Tan[x]"),
        ("Sin[Times[-1, x]]", "Times[-1, Sin[x]]"),
        ("Cos[Times[-1, x]]", "Cos[x]"),
        ("Tan[Times[-2, x]]", "Times[-1, Tan[Times[2, x]]]"),
        ("Sin[-1]", "Times[-1, Sin[1]]"),
        ("ArcSin[Times[-1, x]]", "Times[-1, ArcSin[x]]"),
        ("ArcTan[Times[-1, x]]", "Times[-1, ArcTan[x]]"),
    ]);
}

/// giac's `exp(0)`, `ln(1)`, `ln(e)`, `ln(e^3)`, `ln(e^(1/2))`,
/// `ln(exp(2))`, `exp(ln(5))` and `ln(1/e)`.
#[test]
fn exponentials_and_logarithms() {
    check(&[
        ("Exp[0]", "1"),
        ("Log[1]", "0"),
        ("Log[E]", "1"),
        ("Log[Power[E, 3]]", "3"),
        ("Log[Power[E, 1/2]]", "1/2"),
        ("Log[Exp[2]]", "2"),
        ("Exp[Log[5]]", "5"),
        ("Log[Power[E, -1]]", "-1"),
    ]);
}

/// giac's inverse functions at the values the others take.
#[test]
fn inverse_functions() {
    check(&[
        ("ArcSin[0]", "0"),
        ("ArcSin[1/2]", "Times[1/6, Pi]"),
        ("ArcSin[Times[1/2, Power[2, 1/2]]]", "Times[1/4, Pi]"),
        ("ArcSin[Times[1/2, Power[3, 1/2]]]", "Times[1/3, Pi]"),
        ("ArcSin[1]", "Times[1/2, Pi]"),
        ("ArcSin[-1/2]", "Times[-1/6, Pi]"),
        ("ArcSin[-1]", "Times[-1/2, Pi]"),
        ("ArcCos[0]", "Times[1/2, Pi]"),
        ("ArcCos[1/2]", "Times[1/3, Pi]"),
        ("ArcCos[-1]", "Pi"),
        ("ArcCos[Times[-1/2, Power[3, 1/2]]]", "Times[5/6, Pi]"),
        ("ArcCos[1]", "0"),
        ("ArcCos[Times[1/2, Power[2, 1/2]]]", "Times[1/4, Pi]"),
        ("ArcTan[0]", "0"),
        ("ArcTan[1]", "Times[1/4, Pi]"),
        ("ArcTan[Power[3, 1/2]]", "Times[1/3, Pi]"),
        ("ArcTan[Times[1/3, Power[3, 1/2]]]", "Times[1/6, Pi]"),
        ("ArcTan[-1]", "Times[-1/4, Pi]"),
        ("ArcTan[Power[3, -1/2]]", "Times[1/6, Pi]"),
    ]);
}

/// What stays: other multiples of `Pi`, other numbers, shifts that are not
/// quarter turns, `ArcCos` of a negative (giac's `pi/2+asin(x)`), sums
/// whose terms are negative (giac's `sin(-x-1)`), and `ArcSin[2]`.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Sin[Times[1/5, Pi]]",
        "Sin[Times[1/12, Pi]]",
        "Sin[1]",
        "Cos[x]",
        "Sin[Plus[x, Times[1/3, Pi]]]",
        "Tan[Plus[x, Times[1/2, Pi]]]",
        "ArcCos[Times[-1, x]]",
        "Sin[Plus[-1, Times[-1, x]]]",
        "ArcSin[2]",
        "Log[2]",
        "Exp[1]",
    ] {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e.head(), read(input).head(), "{input}: {e:?}");
    }
}

/// Compositions with the inverse functions (CAS-231), as giac evaluates
/// them.
#[test]
fn compositions() {
    check(&[
        ("Sin[ArcSin[x]]", "x"),
        (
            "Cos[ArcSin[x]]",
            "Power[Plus[1, Times[-1, Power[x, 2]]], 1/2]",
        ),
        (
            "Tan[ArcSin[x]]",
            "Times[x, Power[Plus[1, Times[-1, Power[x, 2]]], -1/2]]",
        ),
        (
            "Sin[ArcCos[x]]",
            "Power[Plus[1, Times[-1, Power[x, 2]]], 1/2]",
        ),
        ("Cos[ArcCos[x]]", "x"),
        (
            "Tan[ArcCos[x]]",
            "Times[Power[x, -1], Power[Plus[1, Times[-1, Power[x, 2]]], 1/2]]",
        ),
        (
            "Sin[ArcTan[x]]",
            "Times[x, Power[Plus[1, Power[x, 2]], -1/2]]",
        ),
        ("Cos[ArcTan[x]]", "Power[Plus[1, Power[x, 2]], -1/2]"),
        ("Tan[ArcTan[x]]", "x"),
        ("Sin[ArcSin[1/3]]", "1/3"),
        // Log of Exp for a real argument: x is real in the Real domain.
        ("Log[Exp[x]]", "x"),
        ("Log[Power[E, x]]", "x"),
        ("Log[Exp[Times[I, x]]]", "Log[Exp[Times[I, x]]]"),
        // The other order only in the principal range.
        ("ArcSin[Sin[x]]", "ArcSin[Sin[x]]"),
    ]);
    let mut ctx = Context::new();
    ctx.set_domain(crate::Domain::Complex);
    assert_eq!(ctx.eval(&read("Log[Exp[x]]")).unwrap(), read("Log[Exp[x]]"));
    let inside = Context::new().assuming(read("And[Greater[x, -1], Less[x, 1]]"));
    let positive = Context::new().assuming(read("And[Greater[x, 0], Less[x, 3]]"));
    for (ctx, input, expected) in [
        (&inside, "ArcSin[Sin[x]]", "x"),
        (&inside, "ArcTan[Tan[x]]", "x"),
        (&inside, "ArcCos[Cos[x]]", "ArcCos[Cos[x]]"),
        (&positive, "ArcCos[Cos[x]]", "x"),
        (&positive, "ArcSin[Sin[x]]", "ArcSin[Sin[x]]"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
}

/// The hyperbolic functions' values (CAS-238).
#[test]
fn hyperbolic_values() {
    check(&[
        ("Sinh[0]", "0"),
        ("Cosh[0]", "1"),
        ("Tanh[0]", "0"),
        ("ArcSinh[0]", "0"),
        ("ArcCosh[1]", "0"),
        ("ArcTanh[0]", "0"),
        ("ArcTanh[1]", "Infinity"),
        ("ArcTanh[-1]", "Times[-1, Infinity]"),
        ("ArcCosh[0]", "Times[I, Times[1/2, Pi]]"),
        // Odd and even.
        ("Sinh[Times[-1, x]]", "Times[-1, Sinh[x]]"),
        ("Cosh[Times[-2, x]]", "Cosh[Times[2, x]]"),
        ("Tanh[-3]", "Times[-1, Tanh[3]]"),
        ("ArcSinh[Times[-1, x]]", "Times[-1, ArcSinh[x]]"),
        // Imaginary arguments.
        ("Sinh[Times[I, x]]", "Times[I, Sin[x]]"),
        ("Cosh[Times[I, x]]", "Cos[x]"),
        ("Tanh[Times[I, x]]", "Times[I, Tan[x]]"),
        ("ArcSinh[Times[I, x]]", "Times[I, ArcSin[x]]"),
        // Of logarithms.
        ("Sinh[Log[2]]", "3/4"),
        ("Cosh[Log[2]]", "5/4"),
        ("Tanh[Log[2]]", "3/5"),
        // Compositions.
        ("Sinh[ArcSinh[x]]", "x"),
        ("Cosh[ArcCosh[x]]", "x"),
        ("Tanh[ArcTanh[x]]", "x"),
    ]);
}

/// The hyperbolic functions in `D`, `Limit` and `Integrate` (CAS-238).
#[test]
fn hyperbolic_calculus() {
    check(&[
        ("D[Sinh[x], x]", "Cosh[x]"),
        ("D[Cosh[x], x]", "Sinh[x]"),
        ("D[Tanh[x], x]", "Power[Cosh[x], -2]"),
        ("D[ArcSinh[x], x]", "Power[Plus[1, Power[x, 2]], -1/2]"),
        (
            "D[ArcCosh[x], x]",
            "Times[Power[Plus[-1, x], -1/2], Power[Plus[1, x], -1/2]]",
        ),
        (
            "D[ArcTanh[x], x]",
            "Power[Plus[1, Times[-1, Power[x, 2]]], -1]",
        ),
        ("Limit[Times[Sinh[x], Power[x, -1]], x, 0]", "1"),
        ("Limit[Times[Plus[-1, Cosh[x]], Power[x, -2]], x, 0]", "1/2"),
        ("Limit[Times[ArcTanh[x], Power[x, -1]], x, 0]", "1"),
        ("Limit[Times[ArcSinh[x], Power[x, -1]], x, 0]", "1"),
        ("Limit[Tanh[x], x, Infinity]", "1"),
        ("Integrate[Sinh[x], x]", "Cosh[x]"),
        (
            "Integrate[Cosh[Times[2, x]], x]",
            "Times[1/2, Sinh[Times[2, x]]]",
        ),
        ("Integrate[Tanh[x], x]", "Log[Cosh[x]]"),
    ]);
}
