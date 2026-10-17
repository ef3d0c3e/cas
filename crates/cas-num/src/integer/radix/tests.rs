//! `Integer` radix conversion checked against GMP's strings, `std`'s parsing
//! rules and `rug`'s formatting.

use proptest::prelude::*;
use rug::Integer as Gmp;

use super::Integer;
use crate::ParseErrorKind;
use crate::testkit::{ToRug, differential, strategy};

fn radix() -> impl Strategy<Value = u32> {
    prop_oneof![2 => Just(10u32), 4 => 2..=36u32]
}

differential! {
    fn to_string_radix(a in strategy::integer(), r in radix()) {
        ours: a.to_integer().to_string_radix(r),
        gmp: a.to_rug().to_string_radix(i32::try_from(r).unwrap()),
    }

    fn from_str_radix(a in strategy::integer(), r in radix(), plus in any::<bool>(), upper in any::<bool>()) {
        ours: {
            let mut s = a.to_rug().to_string_radix(i32::try_from(r).unwrap());
            if upper {
                s = s.to_uppercase();
            }
            if plus && !s.starts_with('-') {
                s.insert(0, '+');
            }
            Integer::from_str_radix(&s, r).unwrap()
        },
        gmp: a.to_rug(),
    }

    fn display_and_debug(a in strategy::integer()) {
        ours: {
            let x = a.to_integer();
            prop_assert_eq!(format!("{x:?}"), x.to_string());
            x.to_string()
        },
        gmp: a.to_rug().to_string(),
    }
}

proptest! {
    /// Accepted strings, values and error kinds agree with `i128`, except
    /// where `i128` overflows; then the value agrees with GMP.
    #[test]
    fn parse_rules_match_std(
        s in proptest::string::string_regex("[0-9a-zA-Z+\\- _]{0,6}|[+-]?[0-9a-f]{0,40}").unwrap(),
        r in radix(),
    ) {
        let ours = Integer::from_str_radix(&s, r);
        match i128::from_str_radix(&s, r) {
            Ok(v) => prop_assert_eq!(ours.unwrap(), Integer::from(v)),
            Err(e) => match e.kind() {
                std::num::IntErrorKind::PosOverflow | std::num::IntErrorKind::NegOverflow => {
                    let g = Gmp::from_str_radix(&s, i32::try_from(r).unwrap()).unwrap();
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

    /// Width, fill, alignment, `+`, `#` and zero padding as `rug` formats,
    /// with the sign in front of the prefix for negative values.
    #[test]
    fn format_flags(a in strategy::signed_up_to(3)) {
        let (x, g) = (a.to_integer(), a.to_rug());
        macro_rules! same {
            ($($fmt:literal),*) => {$(
                prop_assert_eq!(format!($fmt, x), format!($fmt, g), "format {}", $fmt);
            )*};
        }
        same!("{}", "{:+}", "{:>50}", "{:<50}", "{:^50}", "{:*^50}", "{:050}", "{:+050}",
              "{:b}", "{:#b}", "{:o}", "{:#o}", "{:#030o}", "{:x}", "{:#X}", "{:#040x}");
    }
}

#[test]
fn edge_cases() {
    let parse = |s| Integer::from_str_radix(s, 10).map_err(|e| e.kind());
    assert_eq!(parse(""), Err(ParseErrorKind::Empty));
    assert_eq!(parse("-"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("+-1"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("--1"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("-0"), Ok(Integer::ZERO));
    assert_eq!(parse("-9223372036854775808"), Ok(Integer::from(i64::MIN)));
    assert!(parse("9223372036854775808").unwrap().to_i64().is_none());
    assert_eq!(Integer::from(i64::MIN).to_string(), "-9223372036854775808");
    assert_eq!(Integer::from(-255).to_string_radix(16), "-ff");
    assert_eq!(format!("{:?}", Integer::from(-255)), "-255");
}
