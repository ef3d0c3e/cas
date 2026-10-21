//! Decimal printing against the exact endpoints of the ball (every point
//! within one unit of the last printed digit; radii bound what they claim),
//! parsing against correctly rounded `Float::from_rational`, round trips,
//! and the required cases (`1 ± 10^-30`, fallbacks, huge exponents).
#![allow(
    clippy::many_single_char_names,
    reason = "balls b, digits n, strings s, values v, radii r"
)]

use proptest::prelude::*;

use super::{DecimalStyle, format};
use crate::ball::{Ball, Float, Mag, Round};
use crate::error::ParseErrorKind;
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses, dyadic_with};
use crate::{Integer, Rational};

/// `10^e` as a rational, any sign of `e`.
fn ten(e: i64) -> Rational {
    let p = Integer::from(10u64).pow(e.unsigned_abs());
    if e >= 0 {
        Rational::from(p)
    } else {
        Rational::new(Integer::ONE, p)
    }
}

/// A printed decimal (`[-]d[.ddd][e±x]`, or `inf`) as `(value, unit of the
/// last digit, significant digits)`, parsed independently of the code under
/// test.
fn read(s: &str) -> (Rational, Rational, usize) {
    let (mant, exp) = match s.split_once('e') {
        Some((m, x)) => (m, x.parse::<i64>().unwrap()),
        None => (s, 0),
    };
    let (neg, mant) = mant.strip_prefix('-').map_or((false, mant), |m| (true, m));
    let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
    assert!(
        int.bytes().chain(frac.bytes()).all(|c| c.is_ascii_digit()),
        "{s}"
    );
    let all = format!("{int}{frac}");
    let sig = all.trim_start_matches('0').len();
    let n = Integer::from_str_radix(&all, 10).unwrap();
    let n = if neg { -n } else { n };
    let last = exp - i64::try_from(frac.len()).unwrap();
    (&Rational::from(n) * &ten(last), ten(last), sig)
}

fn endpoints(b: &Ball) -> (Rational, Rational) {
    let m = b.mid().to_rational();
    let r = b.rad().to_float().unwrap().to_rational();
    (&m - &r, &m + &r)
}

fn abs(q: &Rational) -> Rational {
    if *q < Rational::from(0i64) {
        -q.clone()
    } else {
        q.clone()
    }
}

/// Checks a printed ball; returns the digits shown (0 for `[± R]`).
fn check_printed(b: &Ball, n: u64, style: DecimalStyle, s: &str) -> Result<usize, TestCaseError> {
    let (lo, hi) = endpoints(b);
    if let Some(r) = s.strip_prefix("[± ").and_then(|t| t.strip_suffix(']')) {
        let (r, _, _) = read(r);
        prop_assert!(abs(&lo) <= r && abs(&hi) <= r, "{s} does not bound {b:?}");
        return Ok(0);
    }
    let (digits, radius) = match s.strip_prefix('[').and_then(|t| t.strip_suffix(']')) {
        Some(t) => {
            prop_assert_eq!(style, DecimalStyle::Radius, "{}", s);
            let (d, r) = t.split_once(" ± ").unwrap();
            (d, Some(read(r).0))
        }
        None => (s, None),
    };
    let (v, unit, sig) = read(digits);
    prop_assert!(
        sig <= usize::try_from(n).unwrap(),
        "{s}: more than {n} digits"
    );
    for x in [&lo, &hi] {
        let d = abs(&(x - &v));
        prop_assert!(d <= unit, "{s}: {x} is more than one unit from {v}");
        if let Some(r) = &radius {
            prop_assert!(d <= *r, "{s}: radius misses {x}");
        }
    }
    if style == DecimalStyle::Radius && radius.is_none() {
        prop_assert!(lo == hi && lo == v, "{s}: no radius on an inexact decimal");
    }
    Ok(sig)
}

/// Balls with exact or relative radii, midpoints of up to 3 limbs.
fn ball() -> BoxedStrategy<Ball> {
    (
        dyadic_with(3, -400..=400),
        prop_oneof![Just(None), (0i64..220).prop_map(Some)],
    )
        .prop_map(|(m, j)| {
            let m = m.to_ours();
            let rad = match j {
                Some(j) if !m.is_zero() => Mag::power_of_two(&(m.exponent() - &Integer::from(j))),
                _ => Mag::ZERO,
            };
            Ball::new(m, rad)
        })
        .boxed()
}

proptest! {
    #[test]
    fn printed_digits_are_correct(b in ball(), n in 1u64..45) {
        for style in [DecimalStyle::Digits, DecimalStyle::Radius] {
            let s = b.to_decimal(n, style);
            let sig = check_printed(&b, n, style, &s)?;
            // Not much fewer digits than the ball knows.
            let known = b.rel_accuracy_bits().saturating_mul(30_103) / 100_000;
            let want = usize::try_from(known.saturating_sub(2).max(0)).unwrap().min(usize::try_from(n).unwrap());
            prop_assert!(sig >= want || b.is_exact(), "{s}: {sig} digits, the ball knows ~{known}");
        }
    }

    #[test]
    fn parse_rounds_to_nearest(n in any::<i64>(), e in -40i64..40, prec in 1u64..200, form in 0u8..3) {
        let s = match form {
            0 => format!("{n}e{e}"),
            1 => format!("{n}E+{e}").replace("+-", "-"),
            _ => {
                // The same value with a decimal point.
                let digits = n.unsigned_abs().to_string();
                let sign = if n < 0 { "-" } else { "" };
                format!("{sign}{}.{}e{}", &digits[..1], &digits[1..], e + i64::try_from(digits.len()).unwrap() - 1)
            }
        };
        let q = &Rational::from(n) * &ten(e);
        let b = Ball::from_decimal(&s, prec).unwrap();
        let (m, o) = Float::from_rational(&q, prec, Round::Nearest);
        prop_assert_eq!(b.mid(), &m, "{}", s);
        prop_assert_eq!(b.is_exact(), o == std::cmp::Ordering::Equal, "{}", s);
        prop_assert!(b.contains_rational(&q));
    }

    #[test]
    fn parse_print_round_trip(d in 1u64..u64::MAX, e10 in -40i64..40) {
        // A canonical decimal of the digits of d without trailing zeros.
        let mut ds = d.to_string();
        while ds.ends_with('0') { ds.pop(); }
        let k = u64::try_from(ds.len()).unwrap();
        let s = format(false, ds, e10, false, k);
        let b = Ball::from_decimal(&s, 4 * k + 20).unwrap();
        prop_assert_eq!(b.to_decimal(k, DecimalStyle::Digits), s);
    }
}

#[test]
fn one_plus_minus_ten_to_minus_30() {
    let r = Mag::from_float(&Float::from_rational(&ten(-30), 128, Round::Ceil).0);
    let digits = |n: usize| format!("1.{}", "0".repeat(n - 1));
    let one = Ball::from_integer(&Integer::ONE).add_error(&r);
    assert_eq!(one.to_decimal(31, DecimalStyle::Digits), digits(30));
    assert_eq!(one.to_decimal(10, DecimalStyle::Digits), digits(10));
    // Just below 1: still 1.000…, the digit-boundary case.
    let below = Float::from(1i64).sub(&Float::from(1i64).mul_2exp(&Integer::from(-200)));
    let below = Ball::new(below, r.clone());
    let s = below.to_decimal(40, DecimalStyle::Digits);
    assert!(s.starts_with("1.000000000000000000000000000"), "{s}");
    check_printed(&below, 40, DecimalStyle::Digits, &s).unwrap();
    let s = one.to_decimal(40, DecimalStyle::Radius);
    assert!(
        s.starts_with("[1.00000000000000000000000000000 ± 1.0"),
        "{s}"
    );
}

#[test]
fn known_forms() {
    let exact = |n: i64| Ball::from_integer(&Integer::from(n));
    let dec = |s: &str| Ball::from_decimal(s, 64).unwrap();
    let digits = DecimalStyle::Digits;
    assert_eq!(exact(1000).to_decimal(20, digits), "1000");
    assert_eq!(exact(-25).to_decimal(20, digits), "-25");
    assert_eq!(exact(0).to_decimal(5, digits), "0");
    assert_eq!(dec("0.5").to_decimal(20, digits), "0.5");
    assert_eq!(dec("-2.5").to_decimal(20, digits), "-2.5");
    assert_eq!(exact(10).pow(30, 128).to_decimal(5, digits), "1e+30");
    assert_eq!(dec("0.1").to_decimal(5, digits), "0.10000");
    assert_eq!(dec("0.000123").to_decimal(3, digits), "0.000123");
    assert_eq!(dec("1.5e-7").to_decimal(2, digits), "1.5e-7");
    assert_eq!(dec("123456789").to_decimal(4, digits), "1.235e+8");
    assert_eq!(dec("9.9999").to_decimal(3, digits), "10.0");
    let pi = Ball::const_pi(64);
    let s = pi.to_decimal(30, digits);
    assert!(s.starts_with("3.14159265358979323"), "{s}");
    check_printed(&pi, 30, digits, &s).unwrap();
    let s = pi.to_decimal(30, DecimalStyle::Radius);
    assert!(
        s.starts_with("[3.14159265358979323") && s.contains(" ± "),
        "{s}"
    );
    // No digit known.
    let wide = Ball::new(
        Float::from(1000i64),
        Mag::from_integer(&Integer::from(2000)),
    );
    assert_eq!(wide.to_decimal(5, digits), "[± 3e+3]");
    assert_eq!(
        Ball::new(Float::ZERO, Mag::ONE).to_decimal(5, digits),
        "[± 1]"
    );
    assert_eq!(Ball::INDETERMINATE.to_decimal(5, digits), "[± inf]");
    // Display.
    assert_eq!(format!("{}", exact(1000)), "1000");
    assert_eq!(format!("{pi:.3}"), "3.14");
    assert!(format!("{pi:#.5}").starts_with("[3.1416 ± "));
}

#[test]
fn huge_exponents() {
    // Beyond the exact path's cap: scaled by ball arithmetic.
    for s in ["1e400000", "-3.7e-400000", "1.25e+500001"] {
        let b = Ball::from_decimal(s, 64).unwrap();
        let (r, _, _) = read(&s.replace("e+", "e"));
        let r = r.to_rug();
        let br = bracket(64, |p, rnd| rug::Float::with_val_round(p, &r, rnd).0);
        let mid = b.mid().to_rational().to_rug();
        let rad = b.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, &br).unwrap();
        assert!(b.rel_accuracy_bits() >= 63, "{s}: {b:?}");
    }
    let b = Ball::from_decimal("1e1000000000", 64).unwrap();
    assert_eq!(b.to_decimal(5, DecimalStyle::Digits), "1.0000e+1000000000");
    let tiny = Ball::from_decimal("-4.5e-999999999999", 64).unwrap();
    assert_eq!(
        tiny.to_decimal(2, DecimalStyle::Digits),
        "-4.5e-999999999999"
    );
    // An exponent beyond i64 prints in hex.
    let far =
        Ball::exact(Float::from(3i64).mul_2exp(&(Integer::from(i64::MAX) * Integer::from(4))));
    assert!(
        far.to_decimal(5, DecimalStyle::Digits)
            .starts_with("[0x1.8p+")
    );
}

#[test]
fn parse_rejections() {
    for s in [
        "", "-", ".", "+.", "1e", "1e+", "1.2.3", "abc", " 1", "1 ", "e5", "1e5.5", "+-1", "--1",
        "1_000", "0x10", "inf", "1,5",
    ] {
        assert!(Ball::from_decimal(s, 64).is_err(), "{s:?} accepted");
    }
    assert_eq!(
        Ball::from_decimal("", 64).unwrap_err().kind(),
        ParseErrorKind::Empty
    );
    let far = Ball::from_decimal("1e99999999999999999999", 64).unwrap_err();
    assert_eq!(far.kind(), ParseErrorKind::ExponentRange);
    for (s, v) in [
        ("5.", 5i64),
        (".5e1", 5),
        ("+7", 7),
        ("-0", 0),
        ("00012e-1", 12),
    ] {
        let want = if s == "00012e-1" {
            Ball::from_decimal("1.2", 64).unwrap()
        } else {
            Ball::from_integer(&Integer::from(v))
        };
        assert_eq!(Ball::from_decimal(s, 64).unwrap(), want, "{s}");
    }
}

/// The significant digits of a decimal, without the zeros that end its
/// fraction: those of `2.0` and `1.0e-30` are 1.
fn significant(s: &str) -> usize {
    let mant = s.split_once('e').map_or(s, |(m, _)| m);
    let mant = mant.strip_prefix('-').unwrap_or(mant);
    let mant = match mant.split_once('.') {
        Some((int, frac)) => format!("{int}{}", frac.trim_end_matches('0')),
        None => mant.to_owned(),
    };
    mant.trim_start_matches('0')
        .trim_end_matches('0')
        .len()
        .max(1)
}

proptest! {
    /// The shortest decimal lies in the ball, has a point and at most `max`
    /// digits, and none with one digit fewer does: neither neighbour of the
    /// result on the coarser grid is in the ball, and a farther decimal
    /// would put one of them in it.
    #[test]
    fn shortest_decimals_are_inside_and_shortest(b in ball(), max in 1u64..45) {
        let Some(s) = b.to_shortest_decimal(max) else {
            return Ok(());
        };
        prop_assert!(s.contains('.'), "{}", s);
        let (v, _, _) = read(&s);
        let sig = significant(&s);
        let (lo, hi) = endpoints(&b);
        prop_assert!(lo <= v && v <= hi, "{s} is outside {b:?}");
        prop_assert!(sig <= usize::try_from(max).unwrap(), "{s}: more than {max} digits");
        if sig >= 2 {
            // The unit of the last significant digit, and the grid one
            // digit coarser.
            let (_, unit, all) = read(&s);
            let unit = &unit * &ten(i64::try_from(all - sig).unwrap());
            let coarse = &unit * &Rational::from(10i64);
            let k = &v / &coarse;
            for n in [k.floor(), k.ceil()] {
                let d = &Rational::from(n) * &coarse;
                prop_assert!(!(lo <= d && d <= hi), "{s}: {d} is shorter and inside {b:?}");
            }
        }
    }

    /// A decimal read at 53 bits shows with no more digits than it was
    /// written with.
    #[test]
    fn typed_decimals_show_short(d in 1u64..100_000_000_000_000_000, e10 in -30i64..30) {
        let mut ds = d.to_string();
        while ds.ends_with('0') { ds.pop(); }
        let k = u64::try_from(ds.len()).unwrap();
        let typed = format(false, ds, e10, false, k);
        let b = Ball::from_decimal(&typed, 53).unwrap();
        let s = b.to_shortest_decimal(17).unwrap();
        prop_assert!(significant(&s) <= usize::try_from(k).unwrap(), "{typed} shows as {s}");
    }
}

#[test]
fn shortest_decimals() {
    let d = |s: &str| Ball::from_decimal(s, 53).unwrap();
    let exact = |n: i64| Ball::from_integer(&Integer::from(n));
    let third = exact(1).div(&d("3.0"), 69).round_to_accuracy(53);
    let sum = d("0.1").add(&d("0.2"), 69).round_to_accuracy(53);
    let cases: [(Ball, u64, Option<&str>); 14] = [
        (d("0.1"), 15, Some("0.1")),
        (sum, 15, Some("0.3")),
        (exact(2), 15, Some("2.0")),
        (exact(-25).mul_2exp(&Integer::from(-1)), 15, Some("-12.5")),
        // 1/3 needs more digits than 53 bits carry.
        (third.clone(), 15, None),
        (third, 17, Some("0.3333333333333333")),
        (d("1e-30"), 15, Some("1.0e-30")),
        (d("0.001"), 15, Some("0.001")),
        (d("0.00001"), 15, Some("1.0e-5")),
        (exact(10).pow(20, 128), 15, Some("1.0e+20")),
        (d("123456789.5"), 15, Some("123456789.5")),
        // A ball around 0, and one that reaches the next power of 10.
        (Ball::new(Float::from(0i64), Mag::ONE), 15, Some("0.0")),
        (
            d("9.99").add_error(&Mag::from_float(
                &Float::from(1i64).mul_2exp(&Integer::from(-5)),
            )),
            15,
            Some("10.0"),
        ),
        // A tie between 0.2 and 0.3 goes to the even digit.
        (
            Ball::new(
                Float::from(1i64).mul_2exp(&Integer::from(-2)),
                Mag::from_float(&Float::from(1i64).mul_2exp(&Integer::from(-3))),
            ),
            15,
            Some("0.2"),
        ),
    ];
    for (b, max, expected) in cases {
        assert_eq!(b.to_shortest_decimal(max).as_deref(), expected, "{b:?}");
    }
    assert_eq!(Ball::INDETERMINATE.to_shortest_decimal(15), None);
}

/// The digits and the exponent of a decimal as written: `2.50e-3` is
/// `(250, -5)`.
fn literal_parts(s: &str) -> (Integer, i64) {
    let (mant, exp) = s
        .split_once(['e', 'E'])
        .map_or((s, 0), |(m, e)| (m, e.parse::<i64>().expect("an exponent")));
    let (neg, mant) = mant.strip_prefix('-').map_or((false, mant), |m| (true, m));
    let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
    let digits = crate::Natural::from_str_radix(&format!("{int}{frac}"), 10).expect("digits");
    let shift = i64::try_from(frac.len()).expect("short");
    (Integer::from_parts(neg, digits), exp - shift)
}

#[test]
fn literal_bits_keep_the_working_precision_to_the_round_trip() {
    for (digits, working, bits) in [
        (1, 53, 53),
        (17, 53, 53),
        // ⌈18 log2 10⌉ + 1 = 60 + 1.
        (18, 53, 61),
        (28, 53, 95),
        (32, 100, 100),
        (33, 100, 111),
        (1, 1, 1),
        (2, 1, 1),
        (3, 1, 11),
    ] {
        assert_eq!(
            Ball::literal_bits(digits, working),
            bits,
            "{digits} digits at {working}"
        );
    }
}

#[test]
fn literals() {
    let read = |s: &str, working: u64| {
        let (n, e) = literal_parts(s);
        Ball::from_literal(&n, e, working)
    };
    assert_eq!(read("0.1", 53), Ball::from_decimal("0.1", 53).unwrap());
    assert!(read("0.5", 53).is_exact() && read("2.0", 53).is_exact());
    assert!(read("2E10", 53).is_exact());
    // A long literal keeps its digits.
    let long = "0.1000000000000000000000000001";
    assert_eq!(
        read(long, 53).to_shortest_decimal(28).as_deref(),
        Some(long)
    );
    // The forms that read back, as an f64's shortest form would be.
    for (b, literal) in [
        (Ball::from_decimal("0.1", 53).unwrap(), "0.1"),
        // Ball sums truncate their midpoint: this one is the f64 nearest 0.3.
        (read("0.1", 53).add(&read("0.2", 53), 53), "0.3"),
        (read("0.30000000000000004", 53), "0.30000000000000004"),
        (
            Ball::from_rational(&Rational::new(1.into(), 3.into()), 53),
            "0.3333333333333333",
        ),
        (
            Ball::from_rational(&Rational::new(2.into(), 3.into()), 53),
            "0.6666666666666666",
        ),
        (
            Ball::exact(Float::from(&Integer::from(1)).mul_2exp(&Integer::from(-60))),
            "8.673617379884035e-19",
        ),
        (read("1e23", 53), "1.0e+23"),
        (read("-2.5", 53), "-2.5"),
        (Ball::ZERO, "0.0"),
        (read(long, 53), long),
    ] {
        assert_eq!(b.to_literal(53).as_deref(), Some(literal), "{b}");
    }
}

proptest! {
    /// A ball's literal reads back to its midpoint, and has no more digits
    /// than a float of the working precision needs when the midpoint fits
    /// it.
    #[test]
    fn literals_read_back(b in ball(), working in 1u64..300) {
        let s = b.to_literal(working).expect("finite");
        let (n, e) = literal_parts(&s);
        prop_assert!(s.contains('.'), "{s}");
        let back = Ball::from_literal(&n, e, working);
        prop_assert_eq!(back.mid(), b.mid(), "{}", s);
        if b.mid().bit_length() <= working {
            let round_trip = working * 30_103 / 100_000 + 2;
            prop_assert!(significant(&s) as u64 <= round_trip, "{} at {}", s, working);
        }
    }
}

/// Past `2^16` bits of exponent, the exact digits are not looked for: the
/// powers of 10 would take seconds. `to_decimal` still prints (CAS-148).
#[test]
fn exact_digits_stop_at_huge_exponents() {
    let near = Ball::from_literal(&Integer::from(499), 10_000, 53);
    assert_eq!(near.to_literal(53).as_deref(), Some("4.99e+10002"));
    let far = Ball::from_literal(&Integer::from(499), 100_000, 53);
    assert_eq!(far.to_literal(53), None);
    assert_eq!(far.to_shortest_decimal(15), None);
    assert!(far.to_decimal(15, DecimalStyle::Digits).starts_with("4.99"));
}
