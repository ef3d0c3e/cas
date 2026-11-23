//! Tests of `Area` (CAS-588, CAS-589): each rule's ball holds the exact
//! integral of polynomials of its degree and the exact rule sums, its
//! value on `exp` is giac's, a parameter gives an exact sum, and the forms
//! that stay or give `0` or `Indeterminate`.

use cas_num::ball::ComplexBall;
use cas_num::{Integer, Rational};

use crate::expr::View;
use crate::trace::EvalOptions;
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `Area[f, x == a..b, n, method]` in `FullForm`, `Equation` for `x = a..b`.
fn area(f: &str, a: &str, b: &str, n: &str, method: &str) -> String {
    format!("Area[{f}, Equation[x, Range[{a}, {b}]], {n}, {method}]")
}

/// The ball `Area` gives for `input`, without a warning.
fn ball(ctx: &Context, input: &str) -> ComplexBall {
    let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
    let found = found.unwrap();
    assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    match found.view() {
        View::Ball(b) => ComplexBall::from_real(b.clone()),
        View::ComplexBall(z) => z.clone(),
        _ => panic!("{input}: {found:?}"),
    }
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// Each rule is exact on polynomials of its degree: the ball holds the
/// integral.
#[test]
fn exact_on_its_degree() {
    let ctx = Context::new();
    let line = "Plus[1, Times[2, x]]";
    let cubic = "Plus[Power[x, 3], Times[-1, x]]";
    for (input, value) in [
        (area("3", "0", "2", "3", "LeftRectangle"), q(6, 1)),
        (area("3", "0", "2", "3", "RightRectangle"), q(6, 1)),
        (area(line, "0", "3", "4", "MiddlePoint"), q(12, 1)),
        (area(line, "0", "3", "4", "Trapezoid"), q(12, 1)),
        (area(cubic, "-1", "2", "2", "Simpson"), q(9, 4)),
        (
            area("Power[x, 3]", "0", "1", "1", "RombergTrapezoid"),
            q(1, 4),
        ),
        (
            area("Power[x, 3]", "0", "1", "1", "RombergMidpoint"),
            q(1, 4),
        ),
        (
            area("Power[x, 5]", "0", "1", "2", "RombergTrapezoid"),
            q(1, 6),
        ),
        (
            area("Power[x, 5]", "0", "1", "2", "RombergMidpoint"),
            q(1, 6),
        ),
        (
            area("Power[x, 7]", "0", "1", "3", "RombergTrapezoid"),
            q(1, 8),
        ),
    ] {
        let found = ball(&ctx, &input);
        assert!(found.is_real(), "{input}: {found:?}");
        assert!(found.re().contains_rational(&value), "{input}: {found:?}");
    }
}

/// The ball holds the rule's exact sum, not the integral.
#[test]
fn exact_sums() {
    let ctx = Context::new();
    let square = "Power[x, 2]";
    for (input, value) in [
        (area(square, "0", "1", "5", "Trapezoid"), q(17, 50)),
        (area(square, "0", "1", "5", "MiddlePoint"), q(33, 100)),
        (area(square, "1", "0", "5", "Trapezoid"), q(-17, 50)),
        (area(square, "0", "1", "2", "LeftRectangle"), q(1, 8)),
        (area(square, "0", "1", "2", "RightRectangle"), q(5, 8)),
        (area(square, "1", "1", "3", "Simpson"), q(0, 1)),
    ] {
        let found = ball(&ctx, &input);
        assert!(found.re().contains_rational(&value), "{input}: {found:?}");
    }
}

/// giac's values of `area(exp(x), x = 0..1, n, method)`.
#[test]
fn giac_values() {
    let ctx = Context::new();
    let exp = "Exp[x]";
    for (input, value) in [
        (area(exp, "0", "1", "4", "Trapezoid"), 1.727_221_904_56),
        (area(exp, "0", "1", "4", "MiddlePoint"), 1.713_815_279_77),
        (area(exp, "0", "1", "4", "LeftRectangle"), 1.512_436_676),
        (area(exp, "0", "1", "4", "RightRectangle"), 1.942_007_133_11),
        (area(exp, "0", "1", "4", "Simpson"), 1.718_284_154_7),
        (
            area(exp, "0", "1", "1", "RombergTrapezoid"),
            1.718_861_151_88,
        ),
        (
            area(exp, "0", "1", "2", "RombergTrapezoid"),
            1.718_282_687_92,
        ),
        (
            area(exp, "0", "1", "3", "RombergTrapezoid"),
            1.718_281_828_79,
        ),
        (
            area(exp, "0", "1", "1", "RombergMidpoint"),
            1.717_776_531_97,
        ),
        (
            area(exp, "0", "1", "2", "RombergMidpoint"),
            1.718_280_996_51,
        ),
        (
            area(exp, "0", "1", "3", "RombergMidpoint"),
            1.718_281_828_13,
        ),
    ] {
        let found = ball(&ctx, &input);
        let mid = found.re().mid().to_f64();
        assert!((mid - value).abs() < 1e-10, "{input}: {found:?}");
    }
    // gauss15 is NIntegrate's certified integral, E − 1.
    let found = ball(&ctx, &area(exp, "0", "1", "7", "Gauss15"));
    assert!((found.re().mid().to_f64() - 1.718_281_828_459_045).abs() < 1e-12);
}

/// Exact answers: a parameter's sum, the integral without a method, and
/// `n = 0`.
#[test]
fn exact_answers() {
    let ctx = Context::new();
    for (input, value) in [
        (
            area("Times[a, x]", "0", "1", "5", "Trapezoid"),
            "Times[a, Power[2, -1]]",
        ),
        (
            area("Power[x, 2]", "0", "c", "1", "Trapezoid"),
            "Times[Power[c, 3], Power[2, -1]]",
        ),
        (
            "Area[Exp[x], Equation[x, Range[0, 1]]]".into(),
            "Plus[-1, Exp[1]]",
        ),
        (
            "Area[Power[x, 2], Equal[x, Range[0, 1]]]".into(),
            "Power[3, -1]",
        ),
        (area("Exp[x]", "0", "1", "0", "LeftRectangle"), "0"),
        (area("Exp[x]", "0", "1", "0", "MiddlePoint"), "0"),
        (area("Exp[x]", "0", "1", "0", "Trapezoid"), "Indeterminate"),
        (area("Exp[x]", "0", "1", "0", "Simpson"), "Indeterminate"),
    ] {
        let found = ctx.eval(&read(&input)).unwrap();
        let expected = ctx.eval(&read(value)).unwrap();
        assert_eq!(found, expected, "{input}");
    }
    // Romberg with no level is one trapezoid: (1 + E)/2.
    let found = ball(&ctx, &area("Exp[x]", "0", "1", "0", "RombergTrapezoid"));
    assert!((found.re().mid().to_f64() - 1.859_140_914_229_522_5).abs() < 1e-12);
    // A rational or decimal n is floored: 5/2 and 2.7 are 2.
    for n in ["Times[5, Power[2, -1]]", "2.7"] {
        let found = ball(&ctx, &area("Power[x, 2]", "0", "1", n, "LeftRectangle"));
        assert!(found.re().contains_rational(&q(1, 8)), "{n}: {found:?}");
    }
}

/// What stays: a huge or negative `n`, `n` without a method, an unknown
/// method or integrand, another spec.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        area(
            "Exp[x]",
            "0",
            "1",
            "1000000000000000000000000000000",
            "Trapezoid",
        ),
        area("Exp[x]", "0", "1", "1000000", "RombergTrapezoid"),
        area("Exp[x]", "0", "1", "1000000", "RombergMidpoint"),
        area("Exp[x]", "0", "1", "-1", "Trapezoid"),
        area("Exp[x]", "0", "1", "4", "foo"),
        area("f[x]", "0", "1", "4", "Trapezoid"),
        "Area[Exp[x], Equation[x, Range[0, 1]], 4]".into(),
        "Area[Exp[x], List[x, 0, 1]]".into(),
        "Area[f[x], Equation[x, Range[0, 1]]]".into(),
    ] {
        let e = read(&input);
        assert_eq!(ctx.eval(&e).unwrap().head(), e.head(), "{input}");
    }
}
