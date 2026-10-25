//! Tests of balls in expressions (CAS-83, CAS-DOC-5 §9.1): the builtin,
//! contagion, printing and comparisons.

use std::time::{Duration, Instant};

use cas_num::ball::{Ball, Float, Mag};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use crate::{Context, EvalError, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(text: &str) -> Expr {
    Context::new()
        .eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// The ball of a decimal at 53 bits, as the language will read it
/// (CAS-148).
fn decimal(s: &str) -> Ball {
    Ball::from_decimal(s, 53).unwrap()
}

fn ball_of(e: &Expr) -> Ball {
    match e.view() {
        View::Ball(b) => b.clone(),
        _ => panic!("{e:?} is not a ball"),
    }
}

#[test]
fn the_builtin_makes_balls() {
    // An exact midpoint rounds to the working precision, the radius grows.
    let b = ball_of(&eval("Ball[1/10, 0]"));
    assert_eq!(
        b,
        Ball::from_rational(&Rational::new(1.into(), 10.into()), 53)
    );
    let b = ball_of(&eval("Ball[3, 1/1000]"));
    assert!(b.contains_rational(&Rational::new(3001.into(), 1000.into())));
    assert!(!b.contains_rational(&Rational::new(3002.into(), 1000.into())));
    assert!(ball_of(&eval("Ball[2, 0]")).is_exact());
    // A ball midpoint gets more radius; a ball radius counts its bound.
    let wider = ball_of(&eval("Ball[Ball[1, 1/8], 1/8]"));
    assert_eq!(
        wider.rad().to_float().unwrap(),
        Float::from(1i64).mul_2exp(&Integer::from(-2))
    );
    // Entry by entry over a matrix: a scalar radius for every entry, a list
    // of radii paired with the rows.
    let m = eval("Ball[List[List[1, 2], List[3, 4]], 1/1000]");
    assert_eq!(m.args().len(), 2);
    for row in m.args() {
        for entry in row.args() {
            assert!(!ball_of(entry).is_exact());
        }
    }
    let paired = eval("Ball[List[1, 2], List[0, 1/2]]");
    assert!(ball_of(&paired.args()[0]).is_exact());
    assert!(!ball_of(&paired.args()[1]).is_exact());
    let ctx = Context::new();
    for text in [
        "Ball[List[1, 2], List[0]]",
        "Ball[List[List[1, 2]], List[List[0, 0, 0]]]",
        "Ball[1, List[0, 1]]",
    ] {
        assert_eq!(ctx.eval(&read(text)), Err(EvalError::Shape), "{text}");
    }
    // Symbols, a negative radius and complex midpoints stay.
    for text in [
        "Ball[x, 0]",
        "Ball[1, r]",
        "Ball[1, -1]",
        "Ball[Complex[1, 1], 0]",
        "Ball[1]",
    ] {
        assert_eq!(eval(text), read(text), "{text}");
    }
    // The working precision decides the rounding.
    let mut ctx = Context::new();
    ctx.set_working_precision(200);
    let b = ball_of(&ctx.eval(&read("Ball[1/3, 0]")).unwrap());
    assert!(b.rel_accuracy_bits() >= 199, "{b:?}");
}

#[test]
fn contagion() {
    let tenth = Expr::from(decimal("0.1"));
    let ctx = Context::new();
    let plus = |a: &Expr, b: &Expr| Expr::apply(&Symbol::PLUS.into(), &[a.clone(), b.clone()]);
    let times = |a: &Expr, b: &Expr| Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()]);
    // A ball absorbs exact numbers.
    let sum = ball_of(&ctx.eval(&plus(&tenth, &Expr::from(3))).unwrap());
    assert!(sum.contains_rational(&Rational::new(31.into(), 10.into())));
    let product = ball_of(
        &ctx.eval(&times(
            &tenth,
            &Expr::from(Rational::new(1.into(), 3.into())),
        ))
        .unwrap(),
    );
    assert!(product.contains_rational(&Rational::new(1.into(), 30.into())));
    // Coefficients and like terms: `x + 2.0 x` is `3.0 x`, and a ball
    // coefficient stays, even 1.0.
    let x = Expr::from(Symbol::new("x"));
    let two = Expr::from(Ball::from_integer(&Integer::from(2)));
    let three = Expr::from(Ball::from_integer(&Integer::from(3)));
    assert_eq!(
        ctx.eval(&plus(&x, &times(&two, &x))).unwrap(),
        times(&three, &x)
    );
    let one = Expr::from(Ball::from_integer(&Integer::from(1)));
    assert_eq!(ctx.eval(&times(&one, &x)).unwrap(), times(&one, &x));
    // Zeros of balls keep their mark: `0.0 x` and `x - 1.0 x` are `0.0`.
    let zero = Expr::from(Ball::ZERO);
    assert_eq!(ctx.eval(&times(&zero, &x)).unwrap(), zero);
    let minus_one = Expr::from(Ball::from_integer(&Integer::from(-1)));
    assert_eq!(ctx.eval(&plus(&x, &times(&minus_one, &x))).unwrap(), zero);
    assert_eq!(ctx.eval(&plus(&zero, &x)).unwrap(), x);
    // With a complex number, a complex ball.
    let z = ctx.eval(&plus(&tenth, &read("Complex[0, 1]"))).unwrap();
    let View::ComplexBall(z) = z.view() else {
        panic!("{z:?}");
    };
    assert!(z.re().overlaps(&decimal("0.1")) && z.im().is_exact());
    // Powers: integer, root, and a ball exponent; the real cube root of a
    // negative in the Real domain, a complex power in the Complex domain.
    let power = |b: &Expr, n: &str| Expr::apply(&Symbol::POWER.into(), &[b.clone(), read(n)]);
    let squared = ball_of(&ctx.eval(&power(&tenth, "2")).unwrap());
    assert!(squared.contains_rational(&Rational::new(1.into(), 100.into())));
    let root = ball_of(&ctx.eval(&power(&two, "1/2")).unwrap());
    assert!(root.mul(&root, 80).contains_rational(&Rational::from(2)));
    let cube = Expr::from(Ball::from_integer(&Integer::from(-8)));
    let real_root = ball_of(&ctx.eval(&power(&cube, "1/3")).unwrap());
    assert!(real_root.contains_rational(&Rational::from(-2)));
    let mut complex = Context::new();
    complex.set_domain(crate::Domain::Complex);
    let principal = complex.eval(&power(&cube, "1/3")).unwrap();
    assert!(
        matches!(principal.view(), View::ComplexBall(_)),
        "{principal:?}"
    );
    let two_to_tenth = ball_of(
        &ctx.eval(&Expr::apply(
            &Symbol::POWER.into(),
            &[Expr::from(2), tenth.clone()],
        ))
        .unwrap(),
    );
    assert!(two_to_tenth.overlaps(&Ball::from_decimal("1.0717734625362931", 53).unwrap()));
    // 1/0.0 as 1/0; a ball around 0 to a negative power has no value.
    assert_eq!(
        ctx.eval(&power(&zero, "-1")).unwrap(),
        read("ComplexInfinity")
    );
    let around_zero = Expr::from(Ball::new(Float::from(0i64), Mag::ONE));
    assert_eq!(
        ctx.eval(&power(&around_zero, "-1")).unwrap(),
        read("Indeterminate")
    );
    // Abs.
    let abs = |e: &Expr| Expr::apply(&Symbol::ABS.into(), std::slice::from_ref(e));
    assert_eq!(ctx.eval(&abs(&minus_one)).unwrap(), one);
}

#[test]
fn printing() {
    let show = |b: Ball| Expr::from(b).to_string();
    // The sum of the balls of 1/10 and 2/10 shows as 0.3 (CAS-DOC-5 §6.1).
    let sum = eval(&format!(
        "Plus[{:?}, {:?}]",
        Expr::from(decimal("0.1")),
        Expr::from(decimal("0.2"))
    ));
    assert_eq!(sum.to_string(), "0.3");
    assert_eq!(show(decimal("0.1")), "0.1");
    assert_eq!(show(Ball::from_integer(&Integer::from(2))), "2.0");
    let third = eval(&format!(
        "Times[1, Power[{:?}, -1]]",
        Expr::from(decimal("3.0"))
    ));
    assert_eq!(third.to_string(), "0.333333333333333");
    // A ball that has lost accuracy shows only its correct digits.
    let blurred = decimal("3.14159265358979").add_error(&Mag::from_float(
        &Float::from(1i64).mul_2exp(&Integer::from(-20)),
    ));
    assert_eq!(show(blurred), "3.14159");
    // Negative balls, in sums and products.
    assert_eq!(
        eval(&format!("Plus[x, {:?}]", Expr::from(decimal("-2.5")))).to_string(),
        "-2.5 + x"
    );
    assert_eq!(
        eval(&format!("Times[x, {:?}]", Expr::from(decimal("-2.5")))).to_string(),
        "-2.5*x"
    );
    // A complex ball.
    let z = eval(&format!(
        "Plus[{:?}, Complex[0, -2]]",
        Expr::from(decimal("0.5"))
    ));
    assert_eq!(z.to_string(), "0.5 - 2.0*I");
}

#[test]
fn comparisons() {
    let ctx = Context::new();
    let compare = |head: &str, a: &Expr, b: &Expr| {
        ctx.eval(&Expr::apply(
            &Symbol::new(head).into(),
            &[a.clone(), b.clone()],
        ))
        .unwrap()
    };
    let truth = |b: bool| Expr::from(if b { Symbol::TRUE } else { Symbol::FALSE });
    let sum = eval(&format!(
        "Plus[{:?}, {:?}]",
        Expr::from(decimal("0.1")),
        Expr::from(decimal("0.2"))
    ));
    let point_three = Expr::from(decimal("0.3"));
    // Accurate balls that overlap are equal; disjoint ones are not.
    assert_eq!(compare("Equal", &sum, &point_three), truth(true));
    assert_eq!(compare("Unequal", &sum, &point_three), truth(false));
    assert_eq!(compare("Less", &sum, &point_three), truth(false));
    assert_eq!(compare("LessEqual", &sum, &point_three), truth(true));
    let point_four = Expr::from(decimal("0.4"));
    assert_eq!(compare("Equal", &sum, &point_four), truth(false));
    assert_eq!(compare("Less", &sum, &point_four), truth(true));
    assert_eq!(compare("Greater", &point_four, &sum), truth(true));
    // Against exact numbers.
    let two = Expr::from(Ball::from_integer(&Integer::from(2)));
    assert_eq!(compare("Equal", &two, &Expr::from(2)), truth(true));
    assert_eq!(compare("Less", &Expr::from(1), &two), truth(true));
    assert_eq!(
        compare("Equal", &Expr::from(decimal("0.1")), &read("1/10")),
        truth(true)
    );
    // Within 8 bits of the working precision is accurate enough; 13 bits
    // short is not.
    let blur = |bits: i64| {
        let r = Mag::from_float(&Float::from(1i64).mul_2exp(&Integer::from(-bits)));
        Expr::from(decimal("0.3").add_error(&r))
    };
    assert_eq!(compare("Equal", &blur(50), &point_three), truth(true));
    let e = Expr::apply(&Symbol::EQUAL.into(), &[blur(42), point_three.clone()]);
    assert_eq!(ctx.eval(&e).unwrap(), e);
    // Overlapping balls that have lost accuracy stay unevaluated.
    let blurred = Expr::from(decimal("0.3").add_error(&Mag::from_float(
        &Float::from(1i64).mul_2exp(&Integer::from(-30)),
    )));
    let e = Expr::apply(
        &Symbol::EQUAL.into(),
        &[blurred.clone(), point_three.clone()],
    );
    assert_eq!(ctx.eval(&e).unwrap(), e);
    let e = Expr::apply(&Symbol::LESS.into(), &[blurred, point_three]);
    assert_eq!(ctx.eval(&e).unwrap(), e);
    // A ball around 0 is not certified 0.
    let around_zero = Expr::from(Ball::new(
        Float::from(0i64),
        Mag::ONE.mul_2exp(&Integer::from(-60)),
    ));
    let e = Expr::apply(&Symbol::EQUAL.into(), &[around_zero, Expr::from(0)]);
    assert_eq!(ctx.eval(&e).unwrap(), e);
    // A real ball is never a non-real complex number.
    assert_eq!(compare("Equal", &two, &read("Complex[2, 1]")), truth(false));
    // At a higher working precision, the 53-bit balls are not accurate.
    let mut fine = Context::new();
    fine.set_working_precision(100);
    let e = Expr::apply(&Symbol::EQUAL.into(), &[sum, Expr::from(decimal("0.3"))]);
    assert_eq!(fine.eval(&e).unwrap(), e);
}

#[test]
fn predicates() {
    let tenth = Expr::from(decimal("0.1"));
    let test = |head: &str, e: &Expr| {
        Context::new()
            .eval(&Expr::apply(
                &Symbol::new(head).into(),
                std::slice::from_ref(e),
            ))
            .unwrap()
    };
    let t = Expr::from(Symbol::TRUE);
    let f = Expr::from(Symbol::FALSE);
    assert_eq!(test("NumberQ", &tenth), t);
    assert_eq!(test("RealQ", &tenth), t);
    assert_eq!(
        test(
            "IntegerQ",
            &Expr::from(Ball::from_integer(&Integer::from(2)))
        ),
        f
    );
    assert_eq!(test("Positive", &tenth), t);
    assert_eq!(test("Negative", &tenth), f);
    let around_zero = Expr::from(Ball::new(Float::from(0i64), Mag::ONE));
    let e = Expr::apply(&Symbol::POSITIVE.into(), std::slice::from_ref(&around_zero));
    assert_eq!(Context::new().eval(&e).unwrap(), e);
}

/// Small exact and inexact balls.
fn ball() -> impl Strategy<Value = Expr> {
    (-40i64..40, -6i64..6, prop::option::of(-70i64..-20)).prop_map(|(n, e, r)| {
        let mid = Float::from(n).mul_2exp(&Integer::from(e));
        let rad = r.map_or(Mag::ZERO, |r| Mag::power_of_two(&Integer::from(r)));
        Expr::from(Ball::new(mid, rad))
    })
}

/// Expressions of balls, small integers, `x` and `I`.
fn expr() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        3 => ball(),
        2 => (-3i64..4).prop_map(Expr::from),
        1 => Just(Expr::from(Symbol::new("x"))),
        1 => Just(Expr::from(Symbol::I)),
    ];
    leaf.prop_recursive(3, 16, 3, |inner| {
        (
            prop::sample::select(vec![
                Symbol::PLUS,
                Symbol::TIMES,
                Symbol::POWER,
                Symbol::LIST,
            ]),
            prop::collection::vec(inner, 1..3),
        )
            .prop_map(|(head, args)| Expr::apply(&head.into(), &args))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    /// Evaluating a result with balls gives it back, and caches change
    /// nothing (CAS-53).
    #[test]
    fn ball_expressions_are_idempotent(e in expr()) {
        let ctx = Context::new();
        if let Ok(r) = ctx.eval(&e) {
            prop_assert!(ctx.eval(&r) == Ok(r.clone()), "{:?} gave {:?}", e, r);
            prop_assert!(ctx.eval_uncached(&e) == Ok(r.clone()), "{:?}", e);
            // Full form reads back.
            prop_assert_eq!(Expr::parse_full_form(&format!("{r:?}")).unwrap(), r);
        }
    }
}

/// A result keeps the bits of its least accurate ball, rounded back after
/// the guard bits (CAS-DOC-5 §9.1); exact balls do not lower it.
#[test]
fn results_keep_the_least_accuracy() {
    let ctx = Context::new();
    let plus = |a: &Ball, b: &Ball| {
        let e = Expr::apply(&Symbol::PLUS.into(), &[Expr::from(a), Expr::from(b)]);
        ball_of(&ctx.eval(&e).unwrap())
    };
    let coarse = decimal("0.1");
    let fine = Ball::from_rational(&Rational::new(1.into(), 3.into()), 200);
    let sum = plus(&coarse, &fine);
    assert!(sum.mid().bit_length() <= 53, "{sum:?}");
    // An exact 2.0 next to the 200-bit ball keeps its 200 bits.
    let exact = Ball::from_integer(&Integer::from(2));
    let sum = plus(&exact, &fine);
    assert!(sum.rel_accuracy_bits() >= 190, "{sum:?}");
    assert!(sum.mid().bit_length() > 53, "{sum:?}");
    // Two exact balls give an exact one.
    assert!(plus(&exact, &exact).is_exact());
}

/// Negation is exact: the builder negates balls as numbers.
#[test]
fn negation() {
    let b = decimal("0.1");
    assert_eq!(-Expr::from(b.clone()), Expr::from(b.neg()));
    let x = Expr::from(Symbol::new("x"));
    assert_eq!(
        format!("{:?}", -Expr::from(b) * x),
        format!("Times[{:?}, x]", Expr::from(decimal("-0.1")))
    );
}

/// Exact identities leave balls as they are: `Times[ball, x]` folds its
/// coefficient as `1 · ball`, which ball arithmetic would round, so the
/// product would not be a fixpoint (found by the `eval` fuzz target).
#[test]
fn identities_leave_balls_alone() {
    let ctx = Context::new();
    // Midpoints of 100 bits, longer than the balls' 60 bits of accuracy,
    // which arithmetic would round them to.
    let long = |n: i64| {
        let m = Integer::from(2).pow(99) + Integer::from(n);
        let mid = Float::from(m).mul_2exp(&Integer::from(-100));
        let ball = Ball::new(mid, Mag::power_of_two(&Integer::from(-61)));
        assert!(ball.mid().bit_length() > 90 && ball.rel_accuracy_bits() <= 62);
        ball
    };
    let z = cas_num::ball::ComplexBall::new(long(3), long(-5));
    for b in [Expr::from(long(7)), Expr::from(z)] {
        for head in [Symbol::TIMES, Symbol::PLUS] {
            let e = Expr::apply(&head.into(), &[b.clone(), Expr::from(Symbol::new("x"))]);
            assert_eq!(ctx.eval(&e), Ok(e.clone()), "{e:?}");
        }
    }
}

/// `Decimal[n, e]` is the ball of the literal at the working precision, or
/// at more bits for more digits than that reads back (CAS-148).
#[test]
fn decimals_make_balls() {
    assert_eq!(ball_of(&eval("Decimal[1, -1]")), decimal("0.1"));
    assert!(ball_of(&eval("Decimal[5, -1]")).is_exact());
    assert!(ball_of(&eval("Decimal[2, 10]")).is_exact());
    assert_eq!(eval("Decimal[0, 0]"), Expr::from(Ball::ZERO));
    // 28 digits are read at 95 bits, so they all show.
    let long = ball_of(&eval("Decimal[1000000000000000000000000001, -28]"));
    assert_eq!(
        long.to_shortest_decimal(28).as_deref(),
        Some("0.1000000000000000000000000001")
    );
    // The working precision is read when the literal is made.
    let mut ctx = Context::new();
    ctx.set_working_precision(200);
    let b = ball_of(&ctx.eval(&read("Decimal[1, -1]")).unwrap());
    assert_eq!(b, Ball::from_decimal("0.1", 200).unwrap());
    // What is not a literal stays.
    for input in ["Decimal[x, 1]", "Decimal[1/2, 1]", "Decimal[1]"] {
        assert_eq!(eval(input), read(input), "{input}");
    }
}

/// `exact` gives the simplest fraction that matches the digits a ball shows
/// (D4, CAS-148), and giac's answers where it has the same rule.
#[test]
fn exact_gives_the_simplest_fraction_shown() {
    let cases = [
        ("Exact[Decimal[1, -1]]", "1/10"),
        ("Exact[Decimal[333333333333, -12]]", "1/3"),
        ("Exact[Decimal[3333, -4]]", "3333/10000"),
        ("Exact[Decimal[7, -1]]", "7/10"),
        ("Exact[Decimal[14, -2]]", "7/50"),
        ("Exact[Decimal[314159, -5]]", "314159/100000"),
        ("Exact[Decimal[142857, -6]]", "142857/1000000"),
        ("Exact[Decimal[1428571428571, -13]]", "1/7"),
        ("Exact[Decimal[15, -6]]", "3/200000"),
        (
            "Exact[Decimal[25, -31]]",
            "1/400000000000000000000000000000",
        ),
        ("Exact[Decimal[123456, -3]]", "15432/125"),
        ("Exact[Decimal[-25, -1]]", "-5/2"),
        // A short decimal that a simpler fraction matches to 14 digits,
        // far from a coincidence by the rule (CAS-178).
        ("Exact[Decimal[-151515, -12]]", "-5/33000033"),
        // Exact balls are their value, however many digits it has.
        ("Exact[Decimal[5, -1]]", "1/2"),
        ("Exact[Ball[0x1p-60, 0x0p+0]]", "1/1152921504606846976"),
        ("Exact[Decimal[3, 0]]", "3"),
        // Results of arithmetic, by what they show.
        ("Exact[Plus[Decimal[1, -1], Decimal[2, -1]]]", "3/10"),
        ("Exact[Power[Decimal[3, 0], -1]]", "1/3"),
        // A ball that shows no digit holds 0 here.
        (
            "Exact[Plus[Decimal[1, -1], Decimal[2, -1], Times[-1, Decimal[3, -1]]]]",
            "0",
        ),
        // Inside expressions; exact parts stay.
        ("Exact[Times[Decimal[5, -1], x]]", "Times[1/2, x]"),
        ("Exact[List[Decimal[25, -2], 3, x]]", "List[1/4, 3, x]"),
        ("Exact[x]", "x"),
        (
            "Exact[Plus[Decimal[1, -1], Times[Decimal[2, -1], I]]]",
            "Complex[1/10, 1/5]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = format!("{:?}", eval(input));
        if found != expected {
            wrong.push(format!("{input}: {found}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `exact` makes fractions of a ball's midpoint and radius only within
/// `Limits::number_bits` bits, which their exponents tell first, and a ball
/// that shows no digit and holds 0 is 0 without them: each of these made
/// fractions of about 10^9 bits (CAS-170).
#[test]
fn exact_measures_far_balls_first() {
    let start = Instant::now();
    for ball in [
        "Ball[0x1p-1000000000, 0x0p+0]",
        "Ball[0x3p-1000000000, 0x1p-1000000000]",
    ] {
        let input = format!("Exact[{ball}]");
        assert_eq!(eval(&input), eval(ball), "{input}");
    }
    let input = "Exact[Ball[0x1p-1000000000, 0x1p-999999999]]";
    assert_eq!(format!("{:?}", eval(input)), "0", "{input}");
    let elapsed = start.elapsed();
    assert!(elapsed < Duration::from_secs(5), "{elapsed:?}");
}

proptest! {
    /// A decimal of at most 6 digits is its own exact value, at any
    /// magnitude, unless a simpler fraction matches it to 12 digits far from
    /// a coincidence, `q² w <= 10^-3` for its denominator `q` and half a unit
    /// `w` of the 12th digit: about one decimal in a thousand, as
    /// `-1.51515e-7` is `-5/33000033` (CAS-178).
    #[test]
    fn exact_keeps_short_decimals(n in -999_999i64..=999_999, e in -40i64..40) {
        let input = Expr::apply(&Symbol::EXACT.into(), &[Expr::apply(&Symbol::DECIMAL.into(), &[Expr::from(n), Expr::from(e)])]);
        let ten = |k: u64| Integer::from(10u64).pow(k);
        let value = if e >= 0 {
            Rational::from(Integer::from(n) * ten(e.unsigned_abs()))
        } else {
            Rational::new(Integer::from(n), ten(e.unsigned_abs()))
        };
        let found = Context::new().eval(&input).unwrap().to_rational().expect("a fraction");
        if found != value {
            // The 12th digit's power of 10, from the first digit's.
            let digits = i64::try_from(n.unsigned_abs().to_string().len()).unwrap();
            let k = digits - 1 + e - 11;
            let unit = if k >= 0 { Rational::from(ten(k.unsigned_abs())) } else { Rational::new(Integer::from(1), ten(k.unsigned_abs())) };
            let half = &unit / &Rational::from(2);
            let gap = if found >= value { &found - &value } else { &value - &found };
            prop_assert!(gap <= half, "{} gave {}", value, found);
            let q = Rational::from(found.denom().clone());
            let chance = &(&q * &q) * &half;
            prop_assert!(chance <= Rational::new(Integer::from(1), Integer::from(1000)), "{} gave {}", value, found);
            prop_assert!(found.denom() < value.denom(), "{} gave {}", value, found);
        }
    }

    /// For up to 15 digits, the fraction matches the decimal at the digits
    /// it shows, without the zeros that end it, and at least 12; and it is
    /// no more complex than the decimal.
    #[test]
    fn exact_matches_long_decimals(n in 1i64..999_999_999_999_999, e in -30i64..10) {
        let input = Expr::apply(&Symbol::EXACT.into(), &[Expr::apply(&Symbol::DECIMAL.into(), &[Expr::from(n), Expr::from(e)])]);
        let ten = |k: u64| Integer::from(10u64).pow(k);
        let value = if e >= 0 { Rational::from(Integer::from(n) * ten(e.unsigned_abs())) } else { Rational::new(Integer::from(n), ten(e.unsigned_abs())) };
        let found = Context::new().eval(&input).unwrap().to_rational().expect("a fraction");
        // Within half a unit of its last digit shown, or of the 12th.
        let digits = i64::try_from(n.to_string().len()).unwrap();
        let shown = i64::try_from(n.to_string().trim_end_matches('0').len()).unwrap();
        let unit = {
            let k = digits - 1 + e - shown.max(12) + 1;
            if k >= 0 { Rational::from(ten(k.unsigned_abs())) } else { Rational::new(Integer::from(1), ten(k.unsigned_abs())) }
        };
        let gap = if found >= value { &found - &value } else { &value - &found };
        prop_assert!(gap * Rational::from(2) <= unit, "{} gave {}", value, found);
        prop_assert!(found.denom() <= value.denom());
    }
}

/// `Digits` is the decimal digits the working precision shows (CAS-151).
#[test]
fn digits_reads_the_working_precision() {
    let mut ctx = Context::new();
    assert_eq!(ctx.eval(&read("Digits")), Ok(Expr::from(15)));
    ctx.set_working_precision(crate::digits_bits(30));
    assert_eq!(ctx.eval(&read("Digits")), Ok(Expr::from(30)));
    assert_eq!(crate::digits_bits(15), 53);
    assert_eq!(crate::digits_bits(30), 102);
}

proptest! {
    /// `digits_bits(n)` is the most bits that show `n` digits.
    #[test]
    fn digits_bits_is_the_most_bits_for_the_digits(n in 1u64..10_000_000) {
        let bits = crate::digits_bits(n);
        prop_assert_eq!(crate::shown_digits(bits), n);
        prop_assert_eq!(crate::shown_digits(bits + 1), n + 1);
    }
}
