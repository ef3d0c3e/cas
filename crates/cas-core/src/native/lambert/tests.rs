//! Tests of `LambertW` (CAS-356): giac's values `W(0)`, `W(e)`, `W(−1/e)`,
//! `W(z e^z)`, the branches, balls against giac's 50 digits, `E^W`, the
//! derivative, and what stays.

use cas_num::ball::Ball;

use crate::expr::View;
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

#[test]
fn values() {
    check(&[
        ("LambertW[0]", "0"),
        ("LambertW[E]", "1"),
        ("LambertW[Times[-1, Power[E, -1]]]", "-1"),
        ("LambertW[Times[-1, Power[E, -1]], -1]", "-1"),
        ("LambertW[Times[2, Power[E, 2]]]", "2"),
        ("LambertW[Times[1/2, Power[E, 1/2]]]", "1/2"),
        ("LambertW[Times[-2, Power[E, -2]], -1]", "-2"),
        // On the other branch, z e^z is not W's.
        (
            "LambertW[Times[-2, Power[E, -2]]]",
            "LambertW[Times[-2, Power[E, -2]]]",
        ),
        ("LambertW[E, -1]", "LambertW[E, -1]"),
        ("LambertW[x, 0]", "LambertW[x]"),
        ("LambertW[x, 1]", "LambertW[x, 1]"),
        ("LambertW[2]", "LambertW[2]"),
        // E^W = z/W.
        ("Power[E, LambertW[x]]", "Times[x, Power[LambertW[x], -1]]"),
        (
            "Power[E, Times[-2, LambertW[x, -1]]]",
            "Times[Power[x, -2], Power[LambertW[x, -1], 2]]",
        ),
        ("Exp[LambertW[x]]", "Times[x, Power[LambertW[x], -1]]"),
    ]);
}

/// Balls on both branches against giac's 50 digits; outside the domain the
/// call stays.
#[test]
fn numbers() {
    let mut ctx = Context::new();
    ctx.set_working_precision(150);
    for (input, digits) in [
        (
            "LambertW[1]",
            "0.56714329040978387299996866221035554975381578718651",
        ),
        (
            "LambertW[10]",
            "1.7455280027406993830743012648753899115352881290809",
        ),
        (
            "LambertW[-1/5]",
            "-0.25917110181907374505665195021540670571358833970089",
        ),
        (
            "LambertW[-1/5, -1]",
            "-2.5426413577735264242938061566618482901614749075294",
        ),
    ] {
        let found = ctx.eval(&read(&format!("N[{input}]"))).unwrap();
        let View::Ball(b) = found.view() else {
            panic!("{input}: {found:?}");
        };
        let want = Ball::from_decimal(digits, 200).unwrap();
        assert!(b.overlaps(&want), "{input}: {b:?} against {digits}");
        assert!(b.rel_accuracy_bits() >= 140, "{input}: {b:?}");
    }
    let outside = ctx.eval(&read("LambertW[N[-1/2]]")).unwrap();
    assert_eq!(outside.head(), Some(&read("LambertW")), "{outside:?}");
    let above = ctx.eval(&read("LambertW[N[1/2], -1]")).unwrap();
    assert_eq!(above.head(), Some(&read("LambertW")), "{above:?}");
}

/// giac's `diff(LambertW(x),x)`: `LambertW(x)/x/(1+LambertW(x))`.
#[test]
fn derivatives() {
    check(&[
        (
            "D[LambertW[x], x]",
            "Times[LambertW[x], Power[x, -1], Power[Plus[1, LambertW[x]], -1]]",
        ),
        (
            "D[LambertW[x, -1], x]",
            "Times[LambertW[x, -1], Power[x, -1], Power[Plus[1, LambertW[x, -1]], -1]]",
        ),
    ]);
}
