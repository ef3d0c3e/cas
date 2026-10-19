use std::borrow::Cow;
use std::fmt::Write as _;

use cas_num::{Integer, Natural, Rational};
use proptest::prelude::*;

use super::{
    MAX_DECIMAL_EXPONENT, NumberError, name, number_value, string_value, superscript_value,
};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(n.into(), d.into())
}

fn ten_pow(e: u64) -> Integer {
    Integer::from(Natural::from(10u32).pow(e))
}

#[test]
fn numbers() {
    let cases = [
        ("0", q(0, 1)),
        ("42", q(42, 1)),
        ("007", q(7, 1)),
        ("1_000_000", q(1_000_000, 1)),
        ("0x1f", q(31, 1)),
        ("0xFF", q(255, 1)),
        ("0x1e5", q(485, 1)),
        ("0o17", q(15, 1)),
        ("0b101", q(5, 1)),
        ("0b1_0000", q(16, 1)),
        ("0.1", q(1, 10)),
        ("0.50", q(1, 2)),
        ("00.0", q(0, 1)),
        ("1.5e-3", q(3, 2000)),
        ("1.5E3", q(1500, 1)),
        ("2E10", q(20_000_000_000, 1)),
        ("1e+5", q(100_000, 1)),
        ("1e-0", q(1, 1)),
        ("12.5e-1", q(5, 4)),
        ("1_0.0_1e1", q(1001, 10)),
    ];
    for (text, value) in cases {
        assert_eq!(number_value(text), Ok(value), "{text}");
    }
    assert_eq!(
        number_value("6.02e23"),
        Ok(Rational::from(Integer::from(602) * ten_pow(21)))
    );
    assert_eq!(
        number_value("18446744073709551616"),
        Ok(Rational::from(Integer::from(u64::MAX) + Integer::from(1)))
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
        "1e1000001",
        "1e-1000001",
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
    assert_eq!(
        number_value(&format!("1e{MAX_DECIMAL_EXPONENT}")),
        Ok(Rational::from(ten_pow(MAX_DECIMAL_EXPONENT)))
    );
    assert_eq!(
        number_value(&format!("1e-{MAX_DECIMAL_EXPONENT}")),
        Ok(Rational::new(
            Integer::from(1),
            ten_pow(MAX_DECIMAL_EXPONENT)
        ))
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
    fn decimals_are_exact(
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
        let mut expected = naive_digits(&digits);
        for _ in 0..shift.unsigned_abs() {
            expected = if shift > 0 { expected * q(10, 1) } else { expected / q(10, 1) };
        }
        prop_assert_eq!(number_value(&text), Ok(expected), "{}", text);
    }

    #[test]
    fn prefixed_numbers(n in any::<u64>()) {
        let value = Ok(Rational::from(Integer::from(n)));
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
