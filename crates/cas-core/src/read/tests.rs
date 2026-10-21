use std::str::FromStr;

use cas_num::ball::{Ball, ComplexBall, Float, Mag};
use cas_num::{Complex, Integer, Rational};
use proptest::prelude::*;

use super::{FullFormErrorKind as K, reads_back};
use crate::build;
use crate::context::Context;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(text, the full form it reads as)`.
fn check(cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (text, expected) in cases {
        let found = format!("{:?}", read(text));
        if found != *expected {
            wrong.push(format!("{text:?}: read as {found}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn symbols() {
    check(&[
        ("x", "x"),
        ("x_1", "x_1"),
        ("_t", "_t"),
        ("θ", "θ"),
        // The language's escaped names (CAS-119).
        ("$E", "$E"),
        ("a$b", "a$b"),
        ("Plus", "Plus"),
    ]);
    assert_eq!(read("Plus"), Expr::from(Symbol::PLUS));
    assert_eq!(read("$E"), Expr::from(Symbol::new("$E")));
}

#[test]
fn numbers() {
    let big = "123456789012345678901234567890";
    check(&[
        ("42", "42"),
        ("-42", "-42"),
        ("-0", "0"),
        ("007", "7"),
        (big, big),
        ("-9223372036854775808", "-9223372036854775808"),
        ("9223372036854775808", "9223372036854775808"),
        ("3/2", "3/2"),
        ("-3/2", "-3/2"),
        ("6/4", "3/2"),
        ("4/2", "2"),
        ("0/5", "0"),
    ]);
    // One form per number: inline when it fits, a rational only when it is
    // not an integer.
    assert_eq!(read("-9223372036854775808").as_i64(), Some(i64::MIN));
    assert!(matches!(
        read("9223372036854775808").view(),
        View::Integer(_)
    ));
    assert_eq!(read("4/2").as_i64(), Some(2));
    assert_eq!(read(big), Expr::from(Integer::from_str(big).unwrap()));
    assert_eq!(
        read("-3/2"),
        Expr::from(Rational::from_str("-3/2").unwrap())
    );
}

#[test]
fn complex_numbers() {
    let q = |s: &str| Rational::from_str(s).unwrap();
    check(&[
        ("Complex[0, 1]", "Complex[0, 1]"),
        ("Complex[2, -3/4]", "Complex[2, -3/4]"),
        ("Complex[-1/2, 1]", "Complex[-1/2, 1]"),
        // A real value is a real number (CAS-51).
        ("Complex[2, 0]", "2"),
        ("Complex[1/2, 0]", "1/2"),
        // Anything else is an application.
        ("Complex[x, 1]", "Complex[x, 1]"),
        ("Complex[1]", "Complex[1]"),
        ("Complex[1, 2, 3]", "Complex[1, 2, 3]"),
        ("Complex[Complex[0, 1], 1]", "Complex[Complex[0, 1], 1]"),
        ("Complex[]", "Complex[]"),
        // A complex number as a head.
        ("Complex[0, 1][x]", "Complex[0, 1][x]"),
    ]);
    assert_eq!(read("Complex[0, 1]"), build::i());
    assert_eq!(
        read("Complex[2, -3/4]"),
        Expr::from(Complex::new(q("2"), q("-3/4")))
    );
    assert!(!read("Complex[x, 1]").is_number());
}

#[test]
fn strings() {
    check(&[
        (r#""s""#, r#""s""#),
        (r#""""#, r#""""#),
        (r#""a\"b\\c""#, r#""a\"b\\c""#),
        (r#""\n\t\r""#, r#""\n\t\r""#),
        (r#""\u{1}\u{7f}""#, r#""\u{1}\u{7f}""#),
        // Escapes of printable characters read as the characters.
        (r#""\u{3c0}\u{1F600}""#, "\"π😀\""),
        // Raw characters, line breaks among them.
        ("\"π\nx\"", r#""π\nx""#),
    ]);
    assert_eq!(read(r#""a\"b""#).as_str(), Some("a\"b"));
}

#[test]
fn applications() {
    check(&[
        ("f[]", "f[]"),
        ("f[x]", "f[x]"),
        ("f[x, y]", "f[x, y]"),
        ("f[a][b]", "f[a][b]"),
        ("f[][]", "f[][]"),
        ("f[g[h[x]], List[]]", "f[g[h[x]], List[]]"),
        // Heads that are not symbols.
        ("2[x]", "2[x]"),
        ("-2[x]", "-2[x]"),
        ("1/2[x]", "1/2[x]"),
        (r#""s"[x]"#, r#""s"[x]"#),
        // Nothing is evaluated or normalised.
        ("Plus[1, 2]", "Plus[1, 2]"),
        ("Times[y, x]", "Times[y, x]"),
        ("Plus[x, Plus[y, z]]", "Plus[x, Plus[y, z]]"),
        // Whitespace between tokens.
        (" f [ x ,y ] ", "f[x, y]"),
        ("f[\n  x,\n\ty\r\n]", "f[x, y]"),
        ("Complex [ 0 , 1 ]", "Complex[0, 1]"),
    ]);
    let f = Symbol::new("f");
    let [x, y] = ["x", "y"].map(|n| Expr::from(Symbol::new(n)));
    assert_eq!(read("f[x, y]"), Expr::apply(&f.into(), &[x.clone(), y]));
    let fa = Expr::apply(&f.into(), &[x]);
    assert_eq!(read("f[x][]"), Expr::apply(&fa, &[]));
}

#[test]
fn errors() {
    let cases = [
        ("", 0, K::End),
        ("   ", 3, K::End),
        ("f[", 2, K::End),
        ("f[x", 3, K::End),
        ("f[x,", 4, K::End),
        ("f[x,]", 4, K::Unexpected(']')),
        ("f[,x]", 2, K::Unexpected(',')),
        ("f[x y]", 4, K::Unexpected('y')),
        ("f]", 1, K::Unexpected(']')),
        ("f x", 2, K::Unexpected('x')),
        ("f[x]]", 4, K::Unexpected(']')),
        ("f[x] g", 5, K::Unexpected('g')),
        ("[x]", 0, K::Unexpected('[')),
        ("·", 0, K::Unexpected('·')),
        ("x+y", 1, K::Unexpected('+')),
        ("2x", 1, K::Unexpected('x')),
        ("2.5", 1, K::Unexpected('.')),
        ("-", 0, K::MissingDigits),
        ("-x", 0, K::MissingDigits),
        ("f[-]", 2, K::MissingDigits),
        ("1/", 2, K::End),
        ("1/x", 2, K::Unexpected('x')),
        ("1/-2", 2, K::Unexpected('-')),
        ("1/0", 0, K::ZeroDenominator),
        ("f[-3/0]", 2, K::ZeroDenominator),
        (r#""abc"#, 0, K::UnterminatedString),
        (r#"f["a]"#, 2, K::UnterminatedString),
        (r#""a\q""#, 2, K::InvalidEscape),
        (r#""\"#, 1, K::InvalidEscape),
        (r#""\u{110000}""#, 1, K::InvalidEscape),
        (r#""\u{d800}""#, 1, K::InvalidEscape),
        (r#""\u{}""#, 1, K::InvalidEscape),
        (r#""\u{1234567}""#, 1, K::InvalidEscape),
        (r#""\u{12g}""#, 1, K::InvalidEscape),
        (r#""\u1234""#, 1, K::InvalidEscape),
        (r#""\u{41""#, 1, K::InvalidEscape),
    ];
    let mut wrong = Vec::new();
    for (text, position, kind) in cases {
        match Expr::parse_full_form(text) {
            Ok(e) => wrong.push(format!("{text:?}: read as {e:?}")),
            Err(e) if (e.position, e.kind) != (position, kind) => {
                wrong.push(format!("{text:?}: {e:?}, expected {kind:?} at {position}"));
            }
            Err(e) => assert!(text.is_char_boundary(e.position)),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    let e = Expr::parse_full_form("f[x").unwrap_err();
    assert_eq!(e.to_string(), "at byte 3: unexpected end of text");
}

#[test]
fn deep_text_reads_in_constant_stack() {
    let n = 100_000;
    let nested = format!("{}x{}", "f[".repeat(n), "]".repeat(n));
    // Compared without `assert_eq!`, which would print 300 kB on failure.
    let same = format!("{:?}", read(&nested)) == nested;
    assert!(same, "nested applications read back");
    let heads = format!("x{}", "[]".repeat(n));
    let same = format!("{:?}", read(&heads)) == heads;
    assert!(same, "compound heads read back");
}

#[test]
fn assert_eval_reads_both_sides() {
    assert_eval!("f[x, Complex[0, 1]]", "f[x, Complex[0, 1]]");
    assert_eval!(&Context::new(), "Plus[1, 2]", "3",);
}

#[test]
#[should_panic(expected = "evaluated to f[x]")]
fn assert_eval_reports_the_result() {
    assert_eval!("f[x]", "f[y]");
}

#[test]
#[should_panic(expected = "\"f[\": at byte 2")]
fn assert_eval_reports_unreadable_text() {
    assert_eval!("f[", "f[]");
}

/// Balls with midpoints of up to 128 bits and exponents of any size, exact
/// or with a radius of 1 to 30 bits.
fn ball() -> impl Strategy<Value = Ball> {
    (
        any::<i128>(),
        prop_oneof![-200i64..200, Just(i64::MAX), Just(-(1i64 << 62))],
        prop::option::of((1u32..(1 << 30), -300i64..100)),
    )
        .prop_map(|(m, e, r)| {
            let mid = Float::from(Integer::from(m)).mul_2exp(&Integer::from(e));
            let rad = r.map_or(Mag::ZERO, |(n, k)| {
                Mag::from_float(&Float::from(i64::from(n)).mul_2exp(&Integer::from(k)))
            });
            Ball::new(mid, rad)
        })
}

/// Random expressions whose full form reads back: symbols with identifier
/// names, numbers of every kind, any strings, and applications with any
/// heads.
fn expr() -> impl Strategy<Value = Expr> {
    let q = (any::<i64>(), 1u64..1_000_000)
        .prop_map(|(n, d)| Rational::new(n.into(), Integer::from(d)));
    let leaf = prop_oneof![
        prop::sample::select(vec![
            "x", "y", "f", "$E", "x_1", "_t", "θ", "a$b", "Plus", "Sin", "List", "Complex", "Ball",
        ])
        .prop_map(|n| Expr::from(Symbol::new(n))),
        any::<i64>().prop_map(Expr::from),
        (any::<u64>(), any::<u64>(), any::<bool>()).prop_map(|(a, b, negative)| {
            let digits = format!("{}{a}{b:020}", if negative { "-" } else { "" });
            Expr::from(Integer::from_str(&digits).unwrap())
        }),
        q.clone().prop_map(Expr::from),
        (q.clone(), q.prop_filter("not real", |im| !im.is_zero()))
            .prop_map(|(re, im)| Expr::from(Complex::new(re, im))),
        ball().prop_map(Expr::from),
        (ball(), ball()).prop_map(|(re, im)| Expr::from(ComplexBall::new(re, im))),
        any::<String>().prop_map(|s| Expr::string(&s)),
        // Mostly what the printer escapes.
        "[\\x00-\\x1f\\x7f\"\\\\aπ]{0,6}".prop_map(|s| Expr::string(&s)),
    ];
    leaf.prop_recursive(4, 64, 5, |inner| {
        prop_oneof![
            (
                prop::sample::select(vec!["f", "g", "Plus", "Complex", "Ball", "$F"]),
                prop::collection::vec(inner.clone(), 0..5),
            )
                .prop_map(|(head, args)| Expr::apply(&Symbol::new(head).into(), &args)),
            (inner.clone(), prop::collection::vec(inner, 0..3))
                .prop_map(|(head, args)| Expr::apply(&head, &args)),
        ]
    })
    .prop_filter("reads back", reads_back)
}

/// `text` with whitespace from `spaces` after each `[` and `,` and before
/// each `]` outside strings.
fn spaced(text: &str, spaces: &[&str]) -> String {
    let mut out = String::new();
    let mut next = spaces.iter().cycle();
    let (mut in_string, mut escaped) = (false, false);
    for c in text.chars() {
        if in_string {
            in_string = escaped || c != '"';
            escaped = !escaped && c == '\\';
            out.push(c);
            continue;
        }
        if c == ']' {
            out.push_str(next.next().expect("a cycle"));
        }
        out.push(c);
        match c {
            '"' => in_string = true,
            '[' | ',' => out.push_str(next.next().expect("a cycle")),
            _ => {}
        }
    }
    out
}

proptest! {
    /// Full form reads back as the expression it was printed from.
    #[test]
    fn full_form_reads_back(e in expr()) {
        let text = format!("{e:?}");
        let back = read(&text);
        prop_assert!(back == e, "{} read as {:?}", text, back);
    }

    /// Whitespace between tokens changes nothing.
    #[test]
    fn whitespace_is_ignored(
        e in expr(),
        spaces in prop::collection::vec(
            prop::sample::select(vec!["", " ", "\n", "\t  ", "\r\n "]),
            1..6,
        ),
    ) {
        let text = spaced(&format!("{e:?}"), &spaces);
        let back = read(&text);
        prop_assert!(back == e, "{:?} read as {:?}", text, back);
    }

    /// Reading never panics, and what it reads prints as text that reads
    /// back the same.
    #[test]
    fn any_text_reads_or_fails_cleanly(text in "[a-z$_ \\[\\],0-9/\"\\\\-]{0,40}") {
        if let Ok(e) = Expr::parse_full_form(&text) {
            let printed = format!("{e:?}");
            prop_assert!(read(&printed) == e, "{:?} printed as {}", text, printed);
        }
    }
}

/// Balls read and print exactly, in hexadecimal (CAS-DOC-5 §6.2).
#[test]
fn balls() {
    let tenth = Ball::from_decimal("0.1", 53).unwrap();
    let e = Expr::from(tenth.clone());
    assert_eq!(format!("{e:?}"), "Ball[0x1.999999999999ap-4, 0x1p-57]");
    assert_eq!(read("Ball[0x1.999999999999ap-4, 0x1p-57]"), e);
    // Hexadecimal floats are exact balls, of any exponent; a dyadic
    // midpoint and a radius a `Mag` holds make the ball.
    let exact =
        |n: i64, e: i64| Expr::from(Ball::exact(Float::from(n).mul_2exp(&Integer::from(e))));
    assert_eq!(read("0x1.8p+1"), exact(3, 0));
    assert_eq!(read("-0x1p-2"), exact(-1, -2));
    assert_eq!(
        read("0x1p+100000000000000000000"),
        Expr::from(Ball::exact(
            Float::from(1i64).mul_2exp(&Integer::from_str("100000000000000000000").unwrap())
        ))
    );
    assert_eq!(
        read("Ball[3, 0]"),
        Expr::from(Ball::from_integer(&Integer::from(3)))
    );
    assert_eq!(
        read("Ball[1/4, 1/1024]").view(),
        View::Ball(&Ball::new(
            Float::from(1i64).mul_2exp(&Integer::from(-2)),
            Mag::power_of_two(&Integer::from(-10))
        ))
    );
    // Otherwise an application, which the builtin evaluates.
    for text in ["Ball[1/3, 0]", "Ball[1, -1]", "Ball[x, 0]", "Ball[1, 2, 3]"] {
        assert!(matches!(read(text).view(), View::Apply { .. }), "{text}");
    }
    // A radius of more than 30 bits is not a `Mag`'s.
    let wide = format!("Ball[0, {}]", (1u64 << 31) + 1);
    assert!(matches!(read(&wide).view(), View::Apply { .. }));
    // Complex balls, and a complex ball with an exact zero imaginary part.
    let z = Expr::from(ComplexBall::new(
        tenth.clone(),
        Ball::from_integer(&Integer::from(2)),
    ));
    assert_eq!(
        format!("{z:?}"),
        "Complex[Ball[0x1.999999999999ap-4, 0x1p-57], Ball[0x1p+1, 0x0p+0]]"
    );
    assert_eq!(read(&format!("{z:?}")), z);
    assert_eq!(read("Complex[Ball[0x1.999999999999ap-4, 0x1p-57], 0]"), e);
    assert_eq!(
        read("Complex[Ball[0x1p+0, 0x0p+0], 1/2]"),
        Expr::from(ComplexBall::new(
            Ball::from_integer(&Integer::from(1)),
            Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(-1)))
        ))
    );
    let errors = [("0x", 0), ("0xg", 0), ("Ball[0x1p, 0]", 5), ("0x1p+", 0)];
    for (text, position) in errors {
        let error = Expr::parse_full_form(text).unwrap_err();
        assert_eq!(
            (error.kind, error.position),
            (K::HexFloat, position),
            "{text}"
        );
    }
}
