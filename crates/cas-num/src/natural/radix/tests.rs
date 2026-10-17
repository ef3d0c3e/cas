//! `Natural` radix conversion checked against GMP's strings, `std`'s parsing
//! rules and `rug`'s formatting.

use proptest::prelude::*;
use rug::Integer;

use super::{FROM_STRING_DC_THRESHOLD, TO_STRING_DC_THRESHOLD};
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, ToRug, differential};
use crate::{Natural, ParseErrorKind};

fn radix() -> impl Strategy<Value = u32> {
    prop_oneof![2 => Just(10u32), 1 => Just(16u32), 1 => Just(2u32), 4 => 2..=36u32]
}

/// Values on both sides of the divide-and-conquer threshold and several
/// recursion levels above it.
fn large() -> BoxedStrategy<Mag> {
    let t = crate::testkit::tested(TO_STRING_DC_THRESHOLD);
    strategy::natural_with_len(prop_oneof![
        len_around(&[t, 2 * t, 4 * t]),
        (t..=5 * t).boxed()
    ])
    .boxed()
}

/// `s` with each letter's case chosen by `mask`.
fn mixed_case(s: &str, mask: u64) -> String {
    s.chars()
        .enumerate()
        .map(|(i, c)| {
            if mask >> (i % 64) & 1 == 1 {
                c.to_ascii_uppercase()
            } else {
                c
            }
        })
        .collect()
}

/// Short strings over digits, letters, signs and junk, for checking the
/// parsing rules against `std`.
fn junk() -> impl Strategy<Value = String> {
    proptest::string::string_regex("[0-9a-zA-Z+\\- _.]{0,6}|[+-]?[0-9a-f]{0,30}").unwrap()
}

differential! {
    fn to_string_radix(a in strategy::natural(), r in radix()) {
        ours: a.to_natural().to_string_radix(r),
        gmp: a.to_rug().to_string_radix(i32::try_from(r).unwrap()),
    }

    fn to_string_radix_large(a in large(), r in radix()) {
        ours: a.to_natural().to_string_radix(r),
        gmp: a.to_rug().to_string_radix(i32::try_from(r).unwrap()),
    }

    /// GMP's digits in mixed case, with leading zeros (so that divide-and-
    /// conquer parts start with zeros) and an optional `+`.
    fn from_str_radix(
        a in prop_oneof![strategy::natural(), large()],
        r in radix(),
        zeros in prop_oneof![Just(0usize), 0..2 * crate::testkit::tested(FROM_STRING_DC_THRESHOLD)],
        plus in any::<bool>(),
        mask in any::<u64>(),
    ) {
        ours: {
            let digits = a.to_rug().to_string_radix(i32::try_from(r).unwrap());
            let s = format!("{}{}{}", if plus { "+" } else { "" }, "0".repeat(zeros), mixed_case(&digits, mask));
            Natural::from_str_radix(&s, r).unwrap()
        },
        gmp: a.to_rug(),
    }

    fn display_and_debug(a in strategy::natural()) {
        ours: {
            let x = a.to_natural();
            prop_assert_eq!(format!("{x:?}"), x.to_string());
            x.to_string()
        },
        gmp: a.to_rug().to_string(),
    }
}

proptest! {
    /// Accepted strings, values and error kinds agree with `u128`, except
    /// where `u128` overflows; then the value agrees with GMP.
    #[test]
    fn parse_rules_match_std(s in junk(), r in radix()) {
        let ours = Natural::from_str_radix(&s, r);
        match u128::from_str_radix(&s, r) {
            Ok(v) => prop_assert_eq!(ours.unwrap(), Natural::from(v)),
            Err(e) => match e.kind() {
                std::num::IntErrorKind::PosOverflow => {
                    let g = Integer::from_str_radix(&s, i32::try_from(r).unwrap()).unwrap();
                    prop_assert!(ours.unwrap().to_rug() == g);
                }
                std::num::IntErrorKind::Empty => {
                    prop_assert_eq!(ours.unwrap_err().kind(), ParseErrorKind::Empty);
                }
                std::num::IntErrorKind::InvalidDigit => {
                    prop_assert_eq!(ours.unwrap_err().kind(), ParseErrorKind::InvalidDigit);
                }
                k => prop_assert!(false, "unexpected std error {:?}", k),
            },
        }
    }

    /// Width, fill, alignment, `+`, `#` and zero padding as `rug` formats.
    #[test]
    fn format_flags(a in strategy::natural_up_to(3)) {
        let (x, g) = (a.to_natural(), a.to_rug());
        macro_rules! same {
            ($($fmt:literal),*) => {$(
                prop_assert_eq!(format!($fmt, x), format!($fmt, g), "format {}", $fmt);
            )*};
        }
        same!("{}", "{:+}", "{:>50}", "{:<50}", "{:^50}", "{:*^50}", "{:050}", "{:+050}",
              "{:b}", "{:#b}", "{:o}", "{:#o}", "{:#030o}", "{:x}", "{:#X}", "{:#040x}");
    }

    #[test]
    fn round_trip(a in prop_oneof![strategy::natural(), large()], r in radix()) {
        let x = a.to_natural();
        prop_assert!(Natural::from_str_radix(&x.to_string_radix(r), r).unwrap() == x);
    }
}

#[test]
fn edge_cases() {
    let parse = |s| Natural::from_str_radix(s, 10).map_err(|e| e.kind());
    assert_eq!(parse(""), Err(ParseErrorKind::Empty));
    assert_eq!(parse("+"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("-0"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse(" 1"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1_000"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("+007"), Ok(Natural::from(7u64)));
    assert_eq!(parse("0000"), Ok(Natural::ZERO));
    assert_eq!(
        Natural::from_str_radix("zZ", 36),
        Ok(Natural::from(36 * 35 + 35u64))
    );
    assert_eq!(
        Natural::from_str_radix("2", 2).map_err(|e| e.kind()),
        Err(ParseErrorKind::InvalidDigit)
    );
    assert_eq!(Natural::ZERO.to_string_radix(7), "0");
    assert_eq!(
        "18446744073709551616".parse::<Natural>(),
        Ok(Natural::power_of_two(64))
    );
    assert_eq!(
        format!("{}", Natural::power_of_two(64)),
        "18446744073709551616"
    );
    assert_eq!(
        crate::ParseNumberError::new(ParseErrorKind::InvalidDigit).to_string(),
        "invalid digit found in string"
    );
}

/// Powers of the radix and their neighbours, from one digit to several
/// divide-and-conquer levels: `radix^e + 1` has a long run of zero chunks
/// between non-zero ones, which only padding emits.
#[test]
fn to_string_radix_powers() {
    for radix in [3u32, 7, 10, 36] {
        // Every exponent up to a few levels; past 64 limbs only the basecase runs.
        let max_digits = TO_STRING_DC_THRESHOLD.min(64) * 4 * 64;
        let mut e = 0u32;
        while (e as usize) < max_digits {
            let p = Integer::from(Integer::u_pow_u(radix, e));
            for v in [p.clone() - 1u32, p.clone(), p + 1u32] {
                if v < 0 {
                    continue;
                }
                let ours = Mag::from_rug(&v).to_natural().to_string_radix(radix);
                assert_eq!(ours, v.to_string_radix(i32::try_from(radix).unwrap()));
            }
            e = if e < 64 { e + 1 } else { e + e / 3 };
        }
    }
}

#[test]
#[should_panic(expected = "radix must be in 2..=36")]
fn radix_37_panics() {
    let _ = Natural::ZERO.to_string_radix(37);
}

#[test]
#[should_panic(expected = "radix must be in 2..=36")]
fn radix_1_panics() {
    let _ = Natural::from_str_radix("0", 1);
}

/// Times the decimal basecase against one divide-and-conquer split (a
/// division by the nearest power `base^(2^i)` and the basecase on both
/// halves), to place `TO_STRING_DC_THRESHOLD`. Run with
/// `cargo test --release -p cas-num tune_to_string -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
fn tune_to_string() {
    use std::time::Instant;
    let (base, k) = super::big_base(10);
    let digits = super::Digits { radix: 10, base, k };
    for n in [8, 12, 16, 20, 24, 28, 32, 40, 48, 64, 96] {
        let x = Natural::from_limbs(&vec![0x9e37_79b9_7f4a_7c15; n]);
        let mut power = Natural::from(base);
        while power.square().limb_len() <= n.div_ceil(2) {
            power = power.square();
        }
        let time = |f: &dyn Fn() -> Vec<u8>| {
            let reps = 2_000_000 / (n * n);
            let start = Instant::now();
            for _ in 0..reps {
                std::hint::black_box(f());
            }
            start.elapsed().as_secs_f64() * 1e9 / f64::from(u32::try_from(reps).unwrap())
        };
        let basecase = time(&|| {
            let mut out = Vec::new();
            digits.basecase(x.limbs(), None, &mut out);
            out
        });
        let split = time(&|| {
            let (q, r) = x.div_rem(&power);
            let mut out = Vec::new();
            digits.basecase(q.limbs(), None, &mut out);
            digits.basecase(r.limbs(), None, &mut out);
            out
        });
        println!("n = {n:3}: basecase {basecase:8.0} ns, split {split:8.0} ns");
    }
}
