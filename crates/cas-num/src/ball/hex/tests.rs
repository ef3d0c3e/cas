//! Hexadecimal floats: exact round trips (including exponents beyond `i64`),
//! the canonical spelling, the other accepted spellings, and rejections.

use proptest::prelude::*;

use crate::Integer;
use crate::ball::{Float, Mag};
use crate::testkit::float::{dyadic, huge_shift};

proptest! {
    #[test]
    fn round_trip(x in dyadic(), k in huge_shift()) {
        for v in [x.to_ours(), x.mul_2exp(&k).to_ours()] {
            let s = v.to_hex();
            prop_assert_eq!(Float::from_hex(&s).unwrap(), v.clone(), "{}", s);
            // Canonical: a leading 1 and no trailing zero digit.
            let body = s.trim_start_matches('-');
            if !v.is_zero() {
                prop_assert!(body.starts_with("0x1"));
                let mant = &body[..body.find('p').unwrap()];
                prop_assert!(!mant.contains('.') || !mant.ends_with('0'), "{}", s);
            }
        }
    }

    #[test]
    fn mag_round_trip(m in 1u32 << 29..1 << 30, e in -2000i64..2000) {
        let mag = Mag::from_float(&Float::from_integer_exp(&Integer::from(m), &Integer::from(e)));
        prop_assert_eq!(Mag::from_hex(&mag.to_hex()).unwrap(), mag);
    }
}

fn hex(s: &str) -> Float {
    Float::from_hex(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

#[test]
fn known_values() {
    let f = |x: f64| Float::from_f64(x).unwrap();
    assert_eq!(f(std::f64::consts::PI).to_hex(), "0x1.921fb54442d18p+1");
    assert_eq!(f(1.0).to_hex(), "0x1p+0");
    assert_eq!(f(0.5).to_hex(), "0x1p-1");
    assert_eq!(f(-3.0).to_hex(), "-0x1.8p+1");
    assert_eq!(f(0.1).to_hex(), "0x1.999999999999ap-4");
    assert_eq!(f(f64::MIN_POSITIVE).to_hex(), "0x1p-1022");
    assert_eq!(Float::ZERO.to_hex(), "0x0p+0");
    let huge = Float::from(1i64).mul_2exp(&(Integer::from(1i64 << 62) * Integer::from(1000i64)));
    assert_eq!(huge.to_hex(), "0x1p+4611686018427387904000");
    assert_eq!(Mag::INF.to_hex(), "inf");
    assert_eq!(Mag::ZERO.to_hex(), "0x0p+0");
}

#[test]
fn other_spellings() {
    assert_eq!(hex("0X1.8P1"), hex("0x1.8p+1"));
    assert_eq!(hex("0x18p-3"), hex("0x1.8p+1"));
    assert_eq!(hex("0x.8p0"), hex("0x1p-1"));
    assert_eq!(hex("0x8."), hex("0x1p+3"));
    assert_eq!(hex("+0x1"), hex("0x1p+0"));
    assert_eq!(hex("0x00010p-4"), hex("0x1p+0"));
    assert_eq!(hex("-0x0p+0"), Float::ZERO);
    assert_eq!(hex("0xABCdefp0"), Float::from(0xab_cdef_i64));
    assert_eq!(Mag::from_hex("inf").unwrap(), Mag::INF);
    // Rounded up to 30 bits.
    let m = Mag::from_hex("0x1.0000001p+0").unwrap();
    assert!(m.to_float().unwrap() > hex("0x1p+0"));
}

#[test]
fn rejections() {
    for s in [
        "",
        "0x",
        "0xp1",
        "1.5",
        "0x1.8p",
        "0x1g",
        "0x1.8p+",
        "0x1.8p1.5",
        "--0x1",
        "0x1p 1",
        " 0x1",
        "0x1..8",
        "0x1p++1",
    ] {
        assert!(Float::from_hex(s).is_err(), "{s:?} accepted");
    }
    assert!(Mag::from_hex("-0x1p+0").is_err());
    assert!(Mag::from_hex("Inf").is_err());
}
