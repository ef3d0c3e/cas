//! Tests of the incomplete gamma and beta functions (CAS-469): values at
//! 30 digits against 40-digit references (Python's `decimal`, by
//! Legendre's continued fraction and `betacf`; giac's doubles agree), the
//! ends, closed forms, derivatives, and what stays.

use cas_num::{Integer, Rational};

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
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

/// A decimal `d.ddd…e−k` as an exact rational.
fn decimal(text: &str) -> Rational {
    let (mantissa, exponent) = text.split_once('e').unwrap();
    let (whole, fraction) = mantissa.split_once('.').unwrap();
    let digits: Integer = format!("{whole}{fraction}").parse().unwrap();
    let shift = i64::try_from(fraction.len()).unwrap() - exponent.parse::<i64>().unwrap();
    let ten = Integer::from(10);
    if shift >= 0 {
        Rational::new(digits, ten.pow(u64::try_from(shift).unwrap()))
    } else {
        Rational::from(&digits * &ten.pow(u64::try_from(-shift).unwrap()))
    }
}

/// `N[input, 30]` is within `10⁻²⁸` (relatively) of the reference.
fn check_values(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, reference) in cases {
        let found = ctx
            .eval(&Expr::call(Symbol::N, [read(input), Expr::from(30)]))
            .unwrap();
        let r = Expr::from(decimal(reference));
        let gap = Expr::call(
            Symbol::ABS,
            [Expr::call(
                Symbol::PLUS,
                [found.clone(), Expr::from(-1) * r.clone()],
            )],
        );
        let bound = Expr::call(Symbol::TIMES, [r, read("Power[10, -28]")]);
        let small = ctx.eval(&Expr::call(Symbol::LESS, [gap, bound])).unwrap();
        if small.as_symbol() != Some(Symbol::TRUE) {
            wrong.push(format!("{input}: {found:?}, expected {reference}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn values() {
    check_values(&[
        // giac's 0.398820945392, 0.930519442787, 0.017464059206 and
        // 0.0889437231707.
        (
            "LowerGamma[5/2, 3/2]",
            "3.9882094539234462960638425321251850342984e-1",
        ),
        (
            "Gamma[5/2, 3/2]",
            "9.3051944278679239086724135929334038366832e-1",
        ),
        (
            "LowerGamma[5/2, 3/2, 1]",
            "3.00014164121372490900198448375331778846520914019585e-1",
        ),
        (
            "Gamma[5/2, 3/2, 1]",
            "6.99985835878627509099801551624668221153479085980415e-1",
        ),
        (
            "Beta[5/2, 3/2, 3/10]",
            "1.7464059205992957411354920801664887959867e-2",
        ),
        (
            "Beta[5/2, 3/2, 3/10, 1]",
            "8.8943723170665599354496267775358493405288e-2",
        ),
        (
            "Gamma[1/3, 100, 1]",
            "6.40321467929021216959579737374031238609132697157907e-46",
        ),
        (
            "Beta[1/2, 1/3, 99/100, 1]",
            "8.46158403063435878082629336362436842319912538748412e-1",
        ),
    ]);
}

#[test]
fn closed_forms_and_ends() {
    check(&[
        // giac's 1 − e^(−x) and 1 − e^(−x)(x + 1).
        ("LowerGamma[1, x]", "Plus[1, Times[-1, Exp[Times[-1, x]]]]"),
        (
            "LowerGamma[1, x, 1]",
            "Plus[1, Times[-1, Exp[Times[-1, x]]]]",
        ),
        ("Gamma[1, x]", "Exp[Times[-1, x]]"),
        ("LowerGamma[a, 0]", "0"),
        ("Gamma[a, 0]", "Gamma[a]"),
        ("Gamma[a, 0, 1]", "1"),
        ("Beta[a, b, 0]", "0"),
        ("Beta[a, b, 1]", "Beta[a, b]"),
        ("Beta[a, b, 1, 1]", "1"),
        // I_x(1, 1) = x, and I_x(1/2, 1/2) = (2/π) asin √x.
        ("Beta[1, 1, x, 1]", "x"),
        (
            "Beta[1/2, 1/2, x, 1]",
            "Times[2, Power[Pi, -1], ArcSin[Power[x, 1/2]]]",
        ),
        // What stays.
        ("LowerGamma[a, x]", "LowerGamma[a, x]"),
        ("Gamma[a, x, 2]", "Gamma[a, x, 2]"),
        ("Beta[a, b, x]", "Beta[a, b, x]"),
    ]);
}

#[test]
fn derivatives() {
    check(&[
        (
            "D[LowerGamma[a, x], x]",
            "Times[Power[x, Plus[a, -1]], Power[E, Times[-1, x]]]",
        ),
        (
            "D[Gamma[a, x], x]",
            "Times[-1, Power[x, Plus[a, -1]], Power[E, Times[-1, x]]]",
        ),
        (
            "D[Beta[a, b, x], x]",
            "Times[Power[x, Plus[a, -1]], Power[Plus[1, Times[-1, x]], Plus[b, -1]]]",
        ),
        (
            "D[Gamma[a, x, 1], x]",
            "Times[-1, Power[x, Plus[a, -1]], Power[E, Times[-1, x]], Power[Gamma[a], -1]]",
        ),
    ]);
}
