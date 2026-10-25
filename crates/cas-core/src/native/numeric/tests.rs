//! Tests of `N` and of the elementary functions of balls (CAS-84,
//! CAS-DOC-5 §9.2, §9.4), against MPFR.

use cas_num::Rational;
use cas_num::ball::{Ball, ComplexBall};
use proptest::prelude::*;
use rug::Float as Mpfr;
use rug::float::Constant;
use rug::ops::Pow;

use crate::{Context, Domain, EvalOptions, Expr, Limits, View, Warning, ball_decimal};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The value of `e` in `ctx`, with its warnings.
fn value(ctx: &Context, e: &Expr) -> (Expr, Vec<Warning>) {
    let (result, report) = ctx.eval_with(e, &EvalOptions::default());
    (
        result.unwrap_or_else(|err| panic!("{e:?}: {err}")),
        report.warnings,
    )
}

fn eval_in(ctx: &Context, text: &str) -> (Expr, Vec<Warning>) {
    value(ctx, &read(text))
}

fn eval(text: &str) -> (Expr, Vec<Warning>) {
    eval_in(&Context::new(), text)
}

fn ball_of(e: &Expr) -> Ball {
    match e.view() {
        View::Ball(b) => b.clone(),
        _ => panic!("{e:?} is not a ball"),
    }
}

fn complex_ball_of(e: &Expr) -> ComplexBall {
    match e.view() {
        View::ComplexBall(z) => z.clone(),
        _ => panic!("{e:?} is not a complex ball"),
    }
}

/// The significant digits of a decimal, without sign, point or exponent.
fn significant(s: &str) -> String {
    let mantissa = s.split(['e', 'E']).next().unwrap_or(s);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    digits.trim_start_matches('0').to_owned()
}

/// MPFR's `x` to `digits` significant digits, rounded to nearest.
fn mpfr_digits(x: &Mpfr, digits: usize) -> String {
    significant(&x.to_string_radix(10, Some(digits)))
}

/// Whether `b` contains the value MPFR brackets at `x` and its neighbours:
/// `x` is within an ulp of the value, so the ball must hold both ends.
fn contains_mpfr(b: &Ball, x: &Mpfr) -> bool {
    let (mut lo, mut hi) = (x.clone(), x.clone());
    lo.next_down();
    hi.next_up();
    [lo, hi].iter().all(|end| {
        let q = end.to_rational().expect("finite");
        let integer =
            |n: &rug::Integer| cas_num::Integer::from_str_radix(&n.to_string(), 10).unwrap();
        b.contains_rational(&Rational::new(integer(q.numer()), integer(q.denom())))
    })
}

#[test]
fn pi_to_a_thousand_digits_matches_mpfr() {
    let (e, warnings) = eval("N[Pi, 1000]");
    assert_eq!(warnings, []);
    let b = ball_of(&e);
    let pi = Mpfr::with_val(4000, Constant::Pi);
    assert!(contains_mpfr(&b, &pi));
    let shown = ball_decimal(&b, 53);
    assert_eq!(significant(&shown), mpfr_digits(&pi, 1000));
    // Fewer digits, within those of the working precision: the nearest
    // decimal.
    for (d, text) in [
        (1, "3.0"),
        (2, "3.1"),
        (5, "3.1416"),
        (15, "3.14159265358979"),
    ] {
        let (e, _) = eval(&format!("N[Pi, {d}]"));
        assert_eq!(ball_decimal(&ball_of(&e), 53), text, "{d} digits");
    }
    // `N[Pi]` at the working precision, as decimals show.
    let (e, _) = eval("N[Pi]");
    assert_eq!(ball_decimal(&ball_of(&e), 53), "3.14159265358979");
}

#[test]
fn cancellation_is_overcome_by_restarting() {
    // e^(π√163) is within 7.5e-13 of an integer: at the first precision,
    // t + 32 bits, the difference has no digit right; the loop doubles.
    let (e, warnings) =
        eval("N[Plus[Exp[Times[Pi, Sqrt[163]]], Times[-1, Power[640320, 3]], -744], 30]");
    assert_eq!(warnings, []);
    let b = ball_of(&e);
    let p = 1000;
    let mut x = Mpfr::with_val(p, 163);
    x.sqrt_mut();
    x *= Mpfr::with_val(p, Constant::Pi);
    x.exp_mut();
    x -= Mpfr::with_val(p, Mpfr::with_val(p, 640_320).pow(3u32));
    x -= 744;
    assert!(contains_mpfr(&b, &x));
    assert_eq!(significant(&ball_decimal(&b, 53)), mpfr_digits(&x, 30));
    assert!(b.rel_accuracy_bits() >= 100);
}

#[test]
fn zero_cannot_be_certified() {
    let (e, warnings) = eval("N[Sin[Pi]]");
    let b = ball_of(&e);
    assert!(b.contains_rational(&Rational::from(0)));
    // The cap: 8t + 1024 bits for t = 53.
    assert_eq!(
        warnings,
        [Warning::Inaccurate {
            expr: read("Sin[Pi]"),
            digits: 0,
            bits: 1448,
        }]
    );
    // The radius shows the precision reached.
    assert!(b.rad().exponent() < &cas_num::Integer::from(-1400));
    assert_eq!(
        warnings[0].to_string(),
        "the value of Sin[Pi] could not be told from 0 within 1448 bits"
    );
}

#[test]
fn the_limit_on_numbers_bounds_the_precision() {
    let mut ctx = Context::new();
    ctx.set_limits(Limits {
        number_bits: 256,
        ..Limits::default()
    });
    let (_, warnings) = eval_in(&ctx, "N[Sin[Pi]]");
    assert_eq!(
        warnings,
        [Warning::Inaccurate {
            expr: read("Sin[Pi]"),
            digits: 0,
            bits: 256,
        }]
    );
    // A target past the limit: unevaluated, with what is allowed.
    let (e, warnings) = eval_in(&ctx, "N[Pi, 100]");
    assert_eq!(e, read("N[Pi, 100]"));
    assert_eq!(
        warnings[0].to_string(),
        "the number of digits must be a positive integer up to 76: N[Pi, 100]"
    );
    let (e, _) = eval_in(&ctx, "N[Pi, 76]");
    assert_eq!(ball_of(&e).rel_accuracy_bits(), 254);
    // Not a count of digits.
    for text in [
        "N[Pi, 0]",
        "N[Pi, -3]",
        "N[Pi, 1/2]",
        "N[Pi, Decimal[2, 0]]",
    ] {
        let (e, warnings) = eval(text);
        assert_eq!(e.head(), Some(&read("N")), "{text}");
        assert_eq!(warnings.len(), 1, "{text}");
    }
    // Digits still to come.
    let (e, warnings) = eval("N[Pi, k]");
    assert_eq!((e, warnings), (read("N[Pi, k]"), vec![]));
}

/// `f[x]`: 0 below 1, 1/2 at 1, 1 above: the discontinuities of
/// CAS-DOC-5 §9.4.
fn piecewise() -> Context {
    let mut ctx = Context::new();
    ctx.define_delayed(
        &read("f[Pattern[x, Blank[]]]"),
        &read("If[Less[x, 1], 0, If[Equal[x, 1], 1/2, 1]]"),
    )
    .unwrap();
    ctx
}

#[test]
fn discontinuities_are_separated_or_left() {
    let ctx = piecewise();
    // Exact input decides exactly, without balls.
    assert_eq!(eval_in(&ctx, "f[1]"), (read("1/2"), vec![]));
    assert_eq!(eval_in(&ctx, "N[f[1]]").0, eval("N[1/2]").0);
    // A constant equal to 1 that the simplifier does not see: no precision
    // tells it from 1, and `N` stays, saying which comparison it could not
    // decide.
    let one = "Times[2, Sin[Times[1/6, Pi]]]";
    let (e, warnings) = eval_in(&ctx, &format!("N[f[{one}]]"));
    assert_eq!(e.head(), Some(&read("N")));
    assert!(
        matches!(&warnings[..], [Warning::Undecided { bits: 1448, .. }]),
        "{warnings:?}"
    );
    // Just above 1: more precision tells, within the cap of a larger target.
    let above = format!("Plus[{one}, Power[10, -1000]]");
    let (e, warnings) = eval_in(&ctx, &format!("N[f[{above}]]"));
    assert_eq!(e.head(), Some(&read("N")));
    assert!(matches!(&warnings[..], [Warning::Undecided { .. }]));
    let (e, warnings) = eval_in(&ctx, &format!("N[f[{above}], 400]"));
    assert_eq!((ball_of(&e).is_exact(), warnings), (true, vec![]));
    assert_eq!(ball_decimal(&ball_of(&e), 53), "1.0");
    // A measured input, 1 ± 10^-5: more precision gains nothing, and `N`
    // stops early, with the note, still undecided.
    let (e, warnings) = eval_in(&ctx, "N[f[Ball[1, 1/100000]]]");
    assert_eq!(e.head(), Some(&read("N")));
    assert!(
        matches!(
            &warnings[..],
            [
                Warning::InputLimited { .. },
                Warning::Undecided { bits: 170, .. }
            ]
        ),
        "{warnings:?}"
    );
    assert!(warnings[0].is_note());
}

#[test]
fn the_input_limits_the_accuracy() {
    // 1 ± 10^-5 knows 5 digits of its sine, whatever the precision.
    let (e, warnings) = eval("N[Sin[Ball[1, 1/100000]], 20]");
    let b = ball_of(&e);
    let sin1 = Mpfr::with_val(200, 1).sin();
    assert!(contains_mpfr(&b, &sin1));
    assert!(
        matches!(&warnings[..], [Warning::InputLimited { digits: 4, .. }]),
        "{warnings:?}"
    );
    // A decimal is a ball of the working precision (CAS-148).
    let (e, warnings) = eval("N[Decimal[1, -1], 30]");
    assert_eq!(ball_decimal(&ball_of(&e), 53), "0.1");
    assert!(matches!(
        &warnings[..],
        [Warning::InputLimited { digits: 15, .. }]
    ));
}

#[test]
fn n_of_its_result_changes_nothing() {
    for text in [
        "N[Pi]",
        "N[Pi, 5]",
        "N[Pi, 15]",
        "N[Pi, 16]",
        "N[Pi, 30]",
        "N[Times[1/3, Pi], 7]",
        "N[Power[-2, 1/2], 10]",
        "N[Exp[Complex[0, 1]], 10]",
        "N[Plus[Complex[1, 0], Times[Complex[0, 1], Power[10, -20]]], 5]",
    ] {
        let (once, warnings) = eval(text);
        assert_eq!(warnings, [], "{text}");
        let input = read(text);
        let again = Expr::apply(
            &read("N"),
            &[std::slice::from_ref(&once), &input.args()[1..]].concat(),
        );
        assert_eq!(value(&Context::new(), &again), (once, vec![]), "{text}");
    }
}

#[test]
fn n_approximates_values_only() {
    let ctx = Context::new();
    let n = |text: &str| eval_in(&ctx, &format!("N[{text}]")).0;
    // Symbols remain; numbers next to them become balls.
    let partial = n("Plus[x, Pi]");
    assert_eq!(partial.args().len(), 2);
    assert!(contains_mpfr(
        &ball_of(&partial.args()[0]),
        &Mpfr::with_val(100, Constant::Pi)
    ));
    // Kept exact: integer exponents, the units of a product, `E` under a
    // symbolic exponent, and the arguments of other functions.
    for text in [
        "Power[x, 2]",
        "Power[x, -1]",
        "Plus[x, Times[-1, y]]",
        "Times[Complex[0, 1], x]",
        "Power[E, x]",
        "f[1/3]",
        "f[Pi]",
        "Derivative[1][g][x]",
        "Hold[1/3]",
        "x",
        "True",
    ] {
        assert_eq!(n(text), ctx.eval(&read(text)).unwrap(), "{text}");
    }
    // Approximated: coefficients, rational exponents, comparisons, lists,
    // equations, the branches of `If`.
    for (text, balls) in [
        ("Times[2, x]", 1),
        ("Power[x, 1/3]", 1),
        ("Power[E, 2]", 1),
        ("Less[x, 1/3]", 1),
        ("List[1/3, Pi, 2]", 3),
        ("Equation[x, 1/3]", 1),
        ("If[Less[x, 0], 1/3, 2/3]", 3),
    ] {
        let found = n(text)
            .distinct()
            .filter(|e| matches!(e.view(), View::Ball(_)))
            .count();
        assert_eq!(found, balls, "{text}");
    }
}

/// A number's rational exponent stays exact in `N`, so that the odd root of
/// a negative is the real one in the Real domain (O7), as exact evaluation
/// reads it, and the principal one in the Complex domain (CAS-177). It was
/// a ball, which took the principal root in either.
#[test]
fn odd_roots_of_negatives_follow_the_domain() {
    let pi = || Mpfr::with_val(100, Constant::Pi);
    let mut real = Context::new();
    real.set_domain(Domain::Real);
    let (e, _) = eval_in(&real, "N[Power[Times[-1, Pi], 1/3]]");
    assert!(contains_mpfr(&ball_of(&e), &-pi().cbrt()), "{e:?}");
    let (e, _) = eval_in(&real, "N[Power[Sin[4], 1/3]]");
    assert!(
        contains_mpfr(&ball_of(&e), &Mpfr::with_val(100, 4).sin().cbrt()),
        "{e:?}"
    );
    let mut complex = Context::new();
    complex.set_domain(Domain::Complex);
    let (e, _) = eval_in(&complex, "N[Power[Times[-1, Pi], 1/3]]");
    // π^(1/3) e^(iπ/3).
    let z = complex_ball_of(&e);
    let third = Mpfr::with_val(100, pi() / 3);
    assert!(contains_mpfr(z.re(), &(pi().cbrt() * third.clone().cos())));
    assert!(contains_mpfr(z.im(), &(pi().cbrt() * third.sin())));
}

#[test]
fn complex_values() {
    // A negative real to an odd number of halves: its real part is exactly
    // 0, in either domain.
    for domain in [Domain::Real, Domain::Complex] {
        let mut ctx = Context::new();
        ctx.set_domain(domain);
        let z = complex_ball_of(&eval_in(&ctx, "N[Sqrt[-2]]").0);
        assert!(z.re().is_exact() && z.re().mid().is_zero(), "{domain:?}");
        assert!(contains_mpfr(z.im(), &Mpfr::with_val(100, 2).sqrt()));
        let z = complex_ball_of(&eval_in(&ctx, "Power[Decimal[-2, 0], -3/2]").0);
        assert!(z.re().is_exact() && z.re().mid().is_zero());
        // 1/(i 2^1.5) = i / 2^1.5... times i^-3 = i.
        let expected = Mpfr::with_val(100, 2).pow(Mpfr::with_val(100, -1.5));
        assert!(contains_mpfr(z.im(), &expected));
        // The principal logarithm of -1.
        let z = complex_ball_of(&eval_in(&ctx, "N[Log[-1]]").0);
        assert!(z.re().is_exact() && z.re().mid().is_zero());
        assert!(contains_mpfr(z.im(), &Mpfr::with_val(100, Constant::Pi)));
    }
    // Complex accuracy is that of the larger part: an inexact 0 beside 1
    // is no reason to raise the precision.
    let (e, warnings) = eval("N[Exp[Times[Complex[0, 1/2], Pi]], 20]");
    let z = complex_ball_of(&e);
    assert!(z.re().contains_rational(&Rational::from(0)));
    assert!(z.im().contains_rational(&Rational::from(1)));
    assert_eq!(warnings, []);
}

#[test]
fn elementary_functions_of_balls() {
    let half = "Decimal[5, -1]";
    let checks: [(&str, Mpfr); 8] = [
        ("Exp", Mpfr::with_val(200, 0.5).exp()),
        ("Log", Mpfr::with_val(200, 0.5).ln()),
        ("Sin", Mpfr::with_val(200, 0.5).sin()),
        ("Cos", Mpfr::with_val(200, 0.5).cos()),
        ("Tan", Mpfr::with_val(200, 0.5).tan()),
        ("ArcSin", Mpfr::with_val(200, 0.5).asin()),
        ("ArcCos", Mpfr::with_val(200, 0.5).acos()),
        ("ArcTan", Mpfr::with_val(200, 0.5).atan()),
    ];
    for (f, expected) in checks {
        let (e, warnings) = eval(&format!("{f}[{half}]"));
        let b = ball_of(&e);
        assert!(contains_mpfr(&b, &expected), "{f}");
        assert!(
            b.rel_accuracy_bits() >= 52,
            "{f}: {}",
            b.rel_accuracy_bits()
        );
        assert_eq!(warnings, [], "{f}");
    }
    // Exact arguments stay; lists thread.
    assert_eq!(eval("Sin[1/2]").0, read("Sin[1/2]"));
    assert_eq!(eval("Sin[x]").0, read("Sin[x]"));
    assert_eq!(
        eval(&format!("Sin[List[{half}, 1/2]]")).0.args()[1],
        read("Sin[1/2]")
    );
    // No real value: the principal complex one.
    let z = complex_ball_of(&eval("Log[Decimal[-2, 0]]").0);
    assert!(contains_mpfr(z.im(), &Mpfr::with_val(100, Constant::Pi)));
    let z = complex_ball_of(&eval("ArcSin[Decimal[2, 0]]").0);
    assert!(contains_mpfr(
        z.re(),
        &(Mpfr::with_val(100, Constant::Pi) / 2u32)
    ));
    // No value at all.
    let (e, warnings) = eval("Log[Decimal[0, 0]]");
    assert_eq!(e, read("DirectedInfinity[-1]"));
    assert!(matches!(&warnings[..], [Warning::Infinite(_)]));
    let (e, warnings) = eval("Log[Ball[0, 1/1000]]");
    assert_eq!(e, read("Indeterminate"));
    assert!(matches!(&warnings[..], [Warning::Indeterminate(_)]));
    // Complex balls.
    let z = complex_ball_of(&eval("Exp[Times[Complex[0, 1], Decimal[1, 0]]]").0);
    assert!(contains_mpfr(z.re(), &Mpfr::with_val(100, 1).cos()));
    assert!(contains_mpfr(z.im(), &Mpfr::with_val(100, 1).sin()));
}

#[test]
fn comparisons_inside_n_decide_only_what_they_prove() {
    let ctx = Context::new().numeric(100);
    let third = "Ball[1/3, 0]";
    // Overlapping inexact balls: undecided, with a warning for `N`.
    for text in [
        format!("Equal[{third}, 1/3]"),
        format!("Less[{third}, 1/3]"),
        format!("Unequal[{third}, 1/3]"),
        format!("Equal[{third}, {third}]"),
    ] {
        let (e, warnings) = eval_in(&ctx, &text);
        assert_eq!(e.head(), read(&text).head(), "{text}");
        assert!(
            matches!(&warnings[..], [Warning::Undecided { bits: 100, .. }]),
            "{text}"
        );
    }
    // Outside `N`, accurate balls that overlap are equal (D3).
    assert_eq!(eval(&format!("Equal[{third}, 1/3]")).0, read("True"));
    // Disjoint balls and exact values decide.
    for (text, value) in [
        (format!("Less[{third}, 1/2]"), "True"),
        (format!("Equal[{third}, 1/2]"), "False"),
        ("Equal[Decimal[5, -1], 1/2]".to_owned(), "True"),
        ("Less[x, 1]".to_owned(), "Less[x, 1]"),
    ] {
        assert_eq!(eval_in(&ctx, &text), (read(value), vec![]), "{text}");
    }
}

/// The digits a ball shows, `d` of them, zeros after the last as padding.
fn shown_digits_of(b: &Ball, d: usize) -> String {
    let mut shown = significant(&ball_decimal(b, 53))
        .trim_end_matches('0')
        .to_owned();
    while shown.len() < d {
        shown.push('0');
    }
    shown
}

/// Values at a tie of the `d`-th digit, and next to one, of both signs and
/// in several decades, none of them dyadic: `N[x, d]` shows `d` digits, within one unit of the
/// nearest, never a digit fewer (CAS-153: `N[-1422/8000, 4]`, -0.17775,
/// showed `-0.178`).
#[test]
fn digits_next_to_ties() {
    let cases: &[(i64, i64, u64)] = &[
        (-1422, 8000, 4),
        (1422, 8000, 4),
        (17765, 100_000, 4),
        (12345, 100_000, 4),
        (-12345, 100_000, 4),
        (12345, 100, 4),
        (-1_234_565, 100, 6),
        (2001, 2000, 4),
        (1_234_567_890_123_455, 10_000_000_000_000_000, 15),
        (-1_234_567_890_123_455, 10_000_000_000_000_000, 15),
        (25, 1000, 1),
        (35, 100, 1),
        (123_456_785, 1_000_000_000, 8),
        (17_769, 100_000, 4),
    ];
    for &(p, q, d) in cases {
        let text = format!("N[Times[{p}, Power[{q}, -1]], {d}]");
        let (e, warnings) = eval(&text);
        assert_eq!(warnings, vec![], "{text}");
        let b = ball_of(&e);
        assert!(
            b.contains_rational(&Rational::new(p.into(), q.into())),
            "{text}"
        );
        // Dyadic values are exact balls, which show their own digits.
        assert!(!b.is_exact(), "{text}");
        let digits = usize::try_from(d).unwrap();
        let shown = shown_digits_of(&b, digits);
        assert_eq!(shown.len(), digits, "{text}: {e}");
        let mut x = Mpfr::with_val(256, p);
        x /= q;
        let expected = mpfr_digits(&x, digits);
        let gap = shown
            .parse::<u128>()
            .unwrap()
            .abs_diff(expected.parse::<u128>().unwrap());
        assert!(gap <= 1, "{text}: {e} shows {shown}, against {expected}");
    }
    // 1234.5 + 10^-27 is not dyadic, but its midpoint is the tie 1234.5,
    // and half a unit of its 4th digit, 1/2, is dyadic: a radius of exactly
    // 1/2 would hold 1234 and 1235.
    let text =
        "N[Times[1234500000000000000000000000001, Power[1000000000000000000000000000, -1]], 4]";
    let (e, warnings) = eval(text);
    assert_eq!(warnings, vec![], "{text}");
    let b = ball_of(&e);
    assert!(!b.is_exact(), "{text}");
    let shown = shown_digits_of(&b, 4);
    assert!(
        shown == "1234" || shown == "1235",
        "{text}: {e} shows {shown}"
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `N[p/q, d]` holds `p/q`, shows its nearest decimal of `d` digits
    /// (within one unit, the printer's promise), and is its own `N`.
    #[test]
    fn n_of_a_fraction(p in -10_000i64..10_000, q in 1i64..10_000, d in 1u64..39) {
        prop_assume!(p != 0);
        let text = format!("N[Times[{p}, Power[{q}, -1]], {d}]");
        let (e, warnings) = eval(&text);
        prop_assert_eq!(&warnings, &vec![]);
        let b = match e.view() {
            View::Ball(b) => b.clone(),
            _ => return Err(TestCaseError::fail(format!("{text}: {e:?}"))),
        };
        prop_assert!(b.contains_rational(&Rational::new(p.into(), q.into())));
        let mut x = Mpfr::with_val(256, p);
        x /= q;
        let digits = usize::try_from(d).unwrap();
        // Exact values show their own digits; others `d` of them, or fewer
        // when the last are zeros.
        if !b.is_exact() {
            // `2.0e+3` shows one digit: zeros after the last are padding.
            let mut shown = significant(&ball_decimal(&b, 53)).trim_end_matches('0').to_owned();
            prop_assert!(shown.len() <= digits, "{}: {}", text, shown);
            while shown.len() < digits {
                shown.push('0');
            }
            let expected = mpfr_digits(&x, digits);
            let close = shown == expected
                || shown.len() == expected.len()
                    && (shown.parse::<u128>().ok().zip(expected.parse::<u128>().ok()))
                        .is_some_and(|(a, b)| a.abs_diff(b) <= 1);
            prop_assert!(close, "{}: {} against {}", text, shown, expected);
        }
        let again = Expr::apply(&read("N"), &[e.clone(), Expr::from(i64::try_from(d).unwrap())]);
        prop_assert_eq!(value(&Context::new(), &again), (e, vec![]));
    }
}
