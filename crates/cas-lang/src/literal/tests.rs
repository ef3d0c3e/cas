use std::borrow::Cow;
use std::fmt::Write as _;

use cas_num::{Integer, Rational};
use proptest::prelude::*;

use super::{
    MAX_DECIMAL_EXPONENT, Number, NumberError, decimal_balls, name, number_value, string_value,
    superscript_value,
};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(n.into(), d.into())
}

fn integer(n: impl Into<Integer>) -> Number {
    Number::Integer(n.into())
}

fn decimal(digits: impl Into<Integer>, exponent: i64) -> Number {
    Number::Decimal {
        digits: digits.into(),
        exponent,
    }
}

#[test]
fn numbers() {
    let cases = [
        ("0", integer(0)),
        ("42", integer(42)),
        ("007", integer(7)),
        ("1_000_000", integer(1_000_000)),
        ("0x1f", integer(31)),
        ("0xFF", integer(255)),
        ("0x1e5", integer(485)),
        ("0o17", integer(15)),
        ("0b101", integer(5)),
        ("0b1_0000", integer(16)),
        // Decimals keep their digits as written, trailing zeros too.
        ("0.1", decimal(1, -1)),
        ("0.50", decimal(50, -2)),
        ("00.0", decimal(0, -1)),
        ("1.5e-3", decimal(15, -4)),
        ("1.5E3", decimal(15, 2)),
        ("2E10", decimal(2, 10)),
        ("1e+5", decimal(1, 5)),
        ("1e-0", decimal(1, 0)),
        ("12.5e-1", decimal(125, -2)),
        ("1_0.0_1e1", decimal(1001, -1)),
        ("6.02e23", decimal(602, 21)),
    ];
    for (text, value) in cases {
        assert_eq!(number_value(text), Ok(value), "{text}");
    }
    assert_eq!(
        number_value("18446744073709551616"),
        Ok(integer(Integer::from(u64::MAX) + Integer::from(1)))
    );
}

#[test]
fn number_errors() {
    for text in [
        "", ".5", "1_", "1__0", "0b102", "x", "1 2", "1.", "-1", "1e", " 1", "²",
    ] {
        assert_eq!(number_value(text), Err(NumberError::Malformed), "{text:?}");
    }
    for text in [
        "1e10001",
        "1e-10001",
        "1e99999999999999999999",
        "5.5e-99999999",
    ] {
        assert_eq!(
            number_value(text),
            Err(NumberError::ExponentTooLarge),
            "{text:?}"
        );
    }
    // The bound itself is accepted.
    let bound = i64::try_from(MAX_DECIMAL_EXPONENT).unwrap();
    assert_eq!(
        number_value(&format!("1e{MAX_DECIMAL_EXPONENT}")),
        Ok(decimal(1, bound))
    );
    assert_eq!(
        number_value(&format!("1e-{MAX_DECIMAL_EXPONENT}")),
        Ok(decimal(1, -bound))
    );
    assert!(NumberError::ExponentTooLarge.to_string().contains("10^n"));
}

#[test]
fn strings() {
    let cases = [
        (r#""""#, ""),
        (r#""abc""#, "abc"),
        (r#""a\"b""#, "a\"b"),
        (r#""\\""#, "\\"),
        (r#""\n\t""#, "\n\t"),
        (r#""\u{3c0}\u{1F600}\u{0}""#, "π😀\0"),
        ("\"two\nlines\"", "two\nlines"),
        (r#""π \u{20} x""#, "π   x"),
    ];
    for (text, value) in cases {
        assert_eq!(string_value(text).as_deref(), Some(value), "{text}");
    }
    for text in [
        r#""\q""#,
        r#""abc"#,
        "abc",
        r#""a" "b""#,
        r#""\u{d800}""#,
        "",
    ] {
        assert_eq!(string_value(text), None, "{text:?}");
    }
}

#[test]
fn superscripts() {
    let cases = [
        ("²", 2),
        ("³", 3),
        ("¹⁰", 10),
        ("⁻¹", -1),
        ("⁺⁴", 4),
        ("⁰", 0),
        ("⁻⁰", 0),
        ("⁵⁶⁷⁸⁹", 56789),
    ];
    for (text, value) in cases {
        assert_eq!(
            superscript_value(text),
            Some(Integer::from(value)),
            "{text}"
        );
    }
    assert_eq!(
        superscript_value("¹⁸⁴⁴⁶⁷⁴⁴⁰⁷³⁷⁰⁹⁵⁵¹⁶¹⁶"),
        Some(Integer::from(u64::MAX) + Integer::from(1))
    );
    for text in ["⁻", "", "2", "²x", "⁻⁻¹"] {
        assert_eq!(superscript_value(text), None, "{text:?}");
    }
}

#[test]
fn names() {
    assert!(matches!(name("x_1"), Cow::Borrowed("x_1")));
    assert!(matches!(name("θ"), Cow::Borrowed("θ")));
    assert_eq!(name("π"), "pi");
    assert_eq!(name("ℯ"), "e");
    assert_eq!(name("∞"), "inf");
    // A decomposed accent composes; a precomposed one is borrowed.
    assert!(matches!(name("e\u{301}t\u{e9}"), Cow::Owned(s) if s == "\u{e9}t\u{e9}"));
    assert!(matches!(name("\u{e9}"), Cow::Borrowed("\u{e9}")));
}

/// `digits` with a `_` after each digit that `separators` marks, except the
/// last.
fn with_separators(digits: &str, separators: &[bool]) -> String {
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        out.push(c);
        if i + 1 < digits.len() && separators.get(i).copied().unwrap_or(false) {
            out.push('_');
        }
    }
    out
}

/// The value of decimal digits, a digit at a time.
fn naive_digits(digits: &str) -> Rational {
    digits.chars().fold(q(0, 1), |acc, c| {
        acc * q(10, 1) + q(i64::from(c.to_digit(10).unwrap()), 1)
    })
}

proptest! {
    #[test]
    fn decimals_keep_their_digits(
        int in "[0-9]{1,24}",
        frac in proptest::option::of("[0-9]{1,24}"),
        exponent in proptest::option::of((-60i64..60, any::<bool>(), any::<bool>())),
        separators in prop::collection::vec(any::<bool>(), 24),
    ) {
        let mut text = with_separators(&int, &separators);
        let mut digits = int.clone();
        if let Some(frac) = &frac {
            text.push('.');
            text.push_str(&with_separators(frac, &separators));
            digits.push_str(frac);
        }
        let mut shift = -i64::try_from(frac.as_ref().map_or(0, String::len)).unwrap();
        if let Some((e, upper, plus)) = exponent {
            text.push(if upper { 'E' } else { 'e' });
            if plus && e >= 0 {
                text.push('+');
            }
            text.push_str(&e.to_string());
            shift += e;
        }
        let expected = if frac.is_none() && exponent.is_none() {
            Number::Integer(naive_digits(&digits).numer().clone())
        } else {
            decimal(naive_digits(&digits).numer().clone(), shift)
        };
        prop_assert_eq!(number_value(&text), Ok(expected), "{}", text);
    }

    #[test]
    fn prefixed_numbers(n in any::<u64>()) {
        let value = Ok(Number::Integer(Integer::from(n)));
        prop_assert_eq!(number_value(&format!("0x{n:x}")), value.clone());
        prop_assert_eq!(number_value(&format!("0x{n:X}")), value.clone());
        prop_assert_eq!(number_value(&format!("0o{n:o}")), value.clone());
        prop_assert_eq!(number_value(&format!("0b{n:b}")), value);
    }

    /// Any string, escaped in every way the language allows, reads back.
    #[test]
    fn strings_round_trip(chars in prop::collection::vec((any::<char>(), 0u8..3), 0..32)) {
        let mut text = String::from("\"");
        let mut value = String::new();
        for (c, style) in chars {
            value.push(c);
            match (c, style) {
                ('"', _) => text.push_str("\\\""),
                ('\\', _) => text.push_str("\\\\"),
                ('\n', 0) => text.push_str("\\n"),
                ('\t', 0) => text.push_str("\\t"),
                (c, 1) => write!(text, "\\u{{{:x}}}", u32::from(c)).unwrap(),
                (c, _) => text.push(c),
            }
        }
        text.push('"');
        prop_assert_eq!(string_value(&text), Some(value));
    }

    #[test]
    fn superscripts_read_back(n in any::<i64>(), plus in any::<bool>()) {
        const DIGITS: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
        let mut text = String::new();
        if n < 0 {
            text.push('⁻');
        } else if plus {
            text.push('⁺');
        }
        for d in n.unsigned_abs().to_string().bytes() {
            text.push(DIGITS[usize::from(d - b'0')]);
        }
        prop_assert_eq!(superscript_value(&text), Some(Integer::from(n)));
    }
}

/// The session makes each `Decimal[digits, exponent]` a ball, held parts
/// included, and leaves the rest (CAS-148).
#[test]
fn decimals_become_balls() {
    use cas_core::{Expr, Symbol};
    use cas_num::ball::Ball;

    let lit = |digits: i64, exponent: i64| {
        Expr::apply(
            &Symbol::DECIMAL.into(),
            &[Expr::from(digits), Expr::from(exponent)],
        )
    };
    let held = Expr::apply(&Symbol::HOLD.into(), &[lit(1, -1)]);
    let e = Expr::apply(&Symbol::LIST.into(), &[lit(5, -1), held, Expr::symbol("x")]);
    let made = decimal_balls(&e, 53);
    let ball = |s: &str| Expr::from(Ball::from_decimal(s, 53).unwrap());
    let expected = Expr::apply(
        &Symbol::LIST.into(),
        &[
            ball("0.5"),
            Expr::apply(&Symbol::HOLD.into(), &[ball("0.1")]),
            Expr::symbol("x"),
        ],
    );
    assert_eq!(made, expected);
    // At another precision, another ball.
    assert_eq!(
        decimal_balls(&lit(1, -1), 200),
        Expr::from(Ball::from_decimal("0.1", 200).unwrap())
    );
}
