//! Tests of the series engine (CAS-203): giac 1.9's expansions, Laurent
//! and Puiseux terms, cancellation, logarithms, and what has no series.

use cas_num::{Integer, Rational};

use super::{Series, expansion, log_t};
use crate::eval::{Eval, run};
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The series in `t` of `f`, in full form, to `order`.
fn expand(f: &str, order: i64) -> Option<Series> {
    let ctx = Context::new();
    let f = ctx.eval(&read(f)).unwrap();
    let t = Expr::symbol("t");
    let order = Rational::from(Integer::from(order));
    run(ctx.limits(), |recursion| {
        let mut ev = Eval::new(&ctx, recursion);
        expansion(&mut ev, &f, &t, &order).unwrap()
    })
}

/// `(f, order, valuation, ramification, coefficients)`: the expansion.
fn check(cases: &[(&str, i64, i64, i64, &[&str])]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for &(f, order, valuation, ramification, coeffs) in cases {
        let expected: Vec<Expr> = coeffs.iter().map(|c| ctx.eval(&read(c)).unwrap()).collect();
        match expand(f, order) {
            Some(s)
                if s.valuation == valuation
                    && s.ramification == ramification
                    && s.coeffs == expected => {}
            found => wrong.push(format!("{f}: {found:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `series` at 0 of exp, sin(x)/x, 1/(1-x), ln(1+x), sqrt(1+x),
/// tan(x), x/(exp(x)-1), atan, asin, cos(x)^2, exp(x)*sin(x),
/// sin(sin(x)), exp(a*x).
#[test]
fn taylor_series() {
    check(&[
        (
            "Exp[t]",
            6,
            0,
            1,
            &["1", "1", "1/2", "1/6", "1/24", "1/120"],
        ),
        (
            "Times[Sin[t], Power[t, -1]]",
            7,
            0,
            1,
            &["1", "0", "-1/6", "0", "1/120", "0", "-1/5040"],
        ),
        (
            "Power[Plus[1, Times[-1, t]], -1]",
            5,
            0,
            1,
            &["1", "1", "1", "1", "1"],
        ),
        ("Log[Plus[1, t]]", 5, 1, 1, &["1", "-1/2", "1/3", "-1/4"]),
        (
            "Power[Plus[1, t], 1/2]",
            4,
            0,
            1,
            &["1", "1/2", "-1/8", "1/16"],
        ),
        (
            "Tan[t]",
            8,
            1,
            1,
            &["1", "0", "1/3", "0", "2/15", "0", "17/315"],
        ),
        (
            "Times[t, Power[Plus[-1, Exp[t]], -1]]",
            5,
            0,
            1,
            &["1", "-1/2", "1/12", "0", "-1/720"],
        ),
        ("ArcTan[t]", 6, 1, 1, &["1", "0", "-1/3", "0", "1/5"]),
        ("ArcSin[t]", 6, 1, 1, &["1", "0", "1/6", "0", "3/40"]),
        ("Power[Cos[t], 2]", 5, 0, 1, &["1", "0", "-1", "0", "1/3"]),
        ("Times[Exp[t], Sin[t]]", 5, 1, 1, &["1", "1", "1/3", "0"]),
        ("Sin[Sin[t]]", 6, 1, 1, &["1", "0", "-1/3", "0", "1/10"]),
        (
            "Exp[Times[a, t]]",
            4,
            0,
            1,
            &[
                "1",
                "a",
                "Times[1/2, Power[a, 2]]",
                "Times[1/6, Power[a, 3]]",
            ],
        ),
    ]);
}

/// Laurent terms, cancellation that needs more of the children, and
/// Puiseux exponents: giac's 1/sin(x), 1/(x^2*(1-x)), (1-cos(x))/x^2 and
/// sqrt(x+x^2).
#[test]
fn laurent_and_puiseux_series() {
    check(&[
        (
            "Power[Sin[t], -1]",
            4,
            -1,
            1,
            &["1", "0", "1/6", "0", "7/360"],
        ),
        (
            "Times[Power[t, -2], Power[Plus[1, Times[-1, t]], -1]]",
            3,
            -2,
            1,
            &["1", "1", "1", "1", "1"],
        ),
        (
            "Times[Plus[1, Times[-1, Cos[t]]], Power[t, -2]]",
            3,
            0,
            1,
            &["1/2", "0", "-1/24"],
        ),
        (
            "Power[Plus[t, Power[t, 2]], 1/2]",
            3,
            1,
            2,
            &["1", "0", "1/2", "0", "-1/8"],
        ),
    ]);
}

/// `log t` stays a coefficient: `Log[t]` and `t Log[t]`.
#[test]
fn logarithms_of_vanishing_arguments() {
    let log = expand("Log[t]", 2).unwrap();
    assert_eq!((log.valuation, log.coeffs[0].clone()), (0, log_t()));
    let x_log = expand("Times[t, Log[t]]", 2).unwrap();
    assert_eq!((x_log.valuation, x_log.coeffs[0].clone()), (1, log_t()));
}

/// No series: an essential singularity, an unknown function of `t`.
#[test]
fn what_has_no_series() {
    assert_eq!(expand("Exp[Power[t, -1]]", 3), None);
    assert_eq!(expand("f[t]", 3), None);
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check_calls(cases: &[(&str, &str)]) {
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

/// giac's `series` and `taylor` (CAS-204): the forms of the call, at 0, at
/// another point, at infinity; Laurent and Puiseux terms; a series that
/// ends; a logarithm kept.
#[test]
fn the_builtin() {
    let exp5 = "Plus[1, x, Times[1/2, Power[x, 2]], Times[1/6, Power[x, 3]], Times[1/24, Power[x, 4]], Times[1/120, Power[x, 5]], Times[Power[x, 6], OrderSize[x]]]";
    check_calls(&[
        ("Series[Exp[x], List[x, 0, 5]]", exp5),
        ("Series[Exp[x], Equation[x, 0], 5]", exp5),
        ("Series[Exp[x], x, 0, 5]", exp5),
        ("Series[Exp[x], x]", exp5),
        // 1-x+1+(x-1)^2-(x-1)^3+(x-1)^4*order_size(x-1).
        (
            "Series[Power[x, -1], Equation[x, 1], 3]",
            "Plus[1, Times[-1, Plus[-1, x]], Power[Plus[-1, x], 2], Times[-1, Power[Plus[-1, x], 3]], Times[Power[Plus[-1, x], 4], OrderSize[Plus[-1, x]]]]",
        ),
        // 1/x+1/6*x+7/360*x^3+x^4*order_size(x).
        (
            "Series[Power[Sin[x], -1], x, 0, 3]",
            "Plus[Power[x, -1], Times[1/6, x], Times[7/360, Power[x, 3]], Times[Power[x, 4], OrderSize[x]]]",
        ),
        // 1/x+x: the series ends.
        (
            "Series[Plus[Power[x, -1], x], x, 0, 3]",
            "Plus[Power[x, -1], x]",
        ),
        // 1+1/x+1/2*(1/x)^2+1/6*(1/x)^3+(1/x)^4*order_size(1/x).
        (
            "Series[Exp[Power[x, -1]], Equation[x, Infinity], 3]",
            "Plus[1, Power[x, -1], Times[1/2, Power[x, -2]], Times[1/6, Power[x, -3]], Times[Power[x, -4], OrderSize[Power[x, -1]]]]",
        ),
        // sqrt(x)+1/2*sqrt(x)*x+sqrt(x)*x^2*order_size(x).
        (
            "Series[Power[Plus[x, Power[x, 2]], 1/2], x, 0, 2]",
            "Plus[Power[x, 1/2], Times[1/2, Power[x, 3/2]], Times[Power[x, 5/2], OrderSize[x]]]",
        ),
        // ln(x)*x.
        ("Series[Times[x, Log[x]], x, 0, 3]", "Times[x, Log[x]]"),
    ]);
    // What stays: no series, a variable that is no symbol, a bad order.
    let ctx = Context::new();
    for input in [
        "Series[Exp[Power[x, -1]], x, 0, 3]",
        "Series[f[x], x, 0, 3]",
        "Series[Exp[x], 2, 0, 3]",
        "Series[Exp[x], x, 0, -1]",
    ] {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e.head(), read(input).head(), "{input}: {e:?}");
    }
}
