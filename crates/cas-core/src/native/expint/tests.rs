//! Tests of `ExpIntegralEi`, `SinIntegral`, `CosIntegral` and
//! `LogIntegral` (CAS-354): the values at 0 and the infinities, `Si`'s sign,
//! balls against giac's 50 digits, `Ci` of a negative ball, the
//! derivatives, and what stays.

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
        ("ExpIntegralEi[0]", "Times[-1, Infinity]"),
        ("CosIntegral[0]", "Times[-1, Infinity]"),
        ("SinIntegral[0]", "0"),
        ("ExpIntegralEi[Infinity]", "Infinity"),
        ("ExpIntegralEi[Times[-1, Infinity]]", "0"),
        ("SinIntegral[Infinity]", "Times[1/2, Pi]"),
        ("SinIntegral[Times[-1, Infinity]]", "Times[-1/2, Pi]"),
        ("CosIntegral[Infinity]", "0"),
        // Si is odd; the others keep their sign.
        ("SinIntegral[-2]", "Times[-1, SinIntegral[2]]"),
        (
            "SinIntegral[Times[-3, x]]",
            "Times[-1, SinIntegral[Times[3, x]]]",
        ),
        ("ExpIntegralEi[-2]", "ExpIntegralEi[-2]"),
        ("CosIntegral[-2]", "CosIntegral[-2]"),
        ("ExpIntegralEi[x]", "ExpIntegralEi[x]"),
        // Li is Ei of the logarithm, as giac writes it.
        ("LogIntegral[x]", "ExpIntegralEi[Log[x]]"),
        ("LogIntegral[1]", "Times[-1, Infinity]"),
        // Listable.
        ("SinIntegral[List[0, x]]", "List[0, SinIntegral[x]]"),
    ]);
}

/// `e` as a real ball, or the real and imaginary parts of a complex one.
fn ball_of(e: &Expr) -> (Ball, Option<Ball>) {
    match e.view() {
        View::Ball(b) => (b.clone(), None),
        View::ComplexBall(z) => (z.re().clone(), Some(z.im().clone())),
        _ => panic!("{e:?} is not a ball"),
    }
}

/// `N` of each function at 150 bits meets giac's 50 digits.
#[test]
fn numbers() {
    let mut ctx = Context::new();
    ctx.set_working_precision(150);
    for (input, digits) in [
        (
            "ExpIntegralEi[1]",
            "1.8951178163559367554665209343316342690170605817327",
        ),
        (
            "ExpIntegralEi[-30]",
            "-0.3021552010688812544815825045153697921167e-14",
        ),
        (
            "SinIntegral[1]",
            "0.94608307036718301494135331382317965781233795473811",
        ),
        (
            "SinIntegral[30]",
            "1.566756540030351110983731309006798166523",
        ),
        (
            "CosIntegral[1]",
            "0.33740392290096813466264620388915076999757803258573",
        ),
        (
            "CosIntegral[100000]",
            "0.35758791572935135690333184729128548684762053177217e-6",
        ),
    ] {
        let found = ctx.eval(&read(&format!("N[{input}]"))).unwrap();
        let (b, im) = ball_of(&found);
        assert!(im.is_none(), "{input}: {found:?}");
        let want = Ball::from_decimal(digits, 200).unwrap();
        assert!(b.overlaps(&want), "{input}: {b:?} against {digits}");
        assert!(b.rel_accuracy_bits() >= 140, "{input}: {b:?}");
    }
    // Ci(-1) = Ci(1) + iπ.
    let found = ctx.eval(&read("N[CosIntegral[-1]]")).unwrap();
    let (re, im) = ball_of(&found);
    let want =
        Ball::from_decimal("0.33740392290096813466264620388915076999757803258573", 200).unwrap();
    assert!(re.overlaps(&want), "{found:?}");
    assert!(im.unwrap().overlaps(&Ball::const_pi(200)), "{found:?}");
    // A ball meeting the pole stays indeterminate, with a warning.
    let found = ctx.eval(&read("ExpIntegralEi[Ball[0, 1/1024]]")).unwrap();
    assert_eq!(found, read("Indeterminate"));
}

/// The derivatives, giac's: `exp(x)/x`, `sin(x)/x`, `cos(x)/x`, and
/// `1/ln(x)` for `Li` through `Ei(ln x)`.
#[test]
fn derivatives() {
    check(&[
        ("D[ExpIntegralEi[x], x]", "Times[Power[E, x], Power[x, -1]]"),
        (
            "D[SinIntegral[Times[2, x]], x]",
            "Times[Power[x, -1], Sin[Times[2, x]]]",
        ),
        ("D[CosIntegral[x], x]", "Times[Power[x, -1], Cos[x]]"),
        ("D[LogIntegral[x], x]", "Power[Log[x], -1]"),
    ]);
}
