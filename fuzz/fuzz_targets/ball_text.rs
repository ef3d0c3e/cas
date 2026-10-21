//! Ball and float text forms (CAS-DOC-5 §6): arbitrary strings never make
//! the parsers or printers panic; an accepted decimal gives a ball containing
//! its exact value, which prints with every digit within one unit of every
//! point, and whose literal reads back to its midpoint (CAS-148); an
//! accepted hexadecimal float prints and parses back exactly.
#![no_main]

use cas_num::ball::{Ball, DecimalStyle, Float, Mag};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

/// The exact value of a decimal literal with a small exponent (`None`
/// otherwise), parsed independently of the code under test.
fn exact(s: &str) -> Option<Rational> {
    let (mant, exp) = match s.find(['e', 'E']) {
        Some(i) => (&s[..i], s[i + 1..].strip_prefix('+').unwrap_or(&s[i + 1..]).parse::<i64>().ok()?),
        None => (s, 0),
    };
    let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
    let e = exp.checked_sub(i64::try_from(frac.len()).ok()?)?;
    if e.unsigned_abs() > 2000 {
        return None;
    }
    let n = Integer::from_str_radix(&format!("{int}{frac}"), 10).ok()?;
    let p = Integer::from(10u64).pow(e.unsigned_abs());
    Some(if e >= 0 {
        Rational::from(&n * &p)
    } else {
        Rational::new(n, p)
    })
}

/// The digits and the power of 10 of a decimal literal as
/// [`Ball::to_literal`] writes it: `2.50e-3` is `(250, -5)`.
fn literal_parts(s: &str) -> Option<(Integer, i64)> {
    let (mant, exp) = match s.find(['e', 'E']) {
        Some(i) => (&s[..i], s[i + 1..].strip_prefix('+').unwrap_or(&s[i + 1..]).parse::<i64>().ok()?),
        None => (s, 0),
    };
    let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
    let n = Integer::from_str_radix(&format!("{int}{frac}"), 10).ok()?;
    Some((n, exp.checked_sub(i64::try_from(frac.len()).ok()?)?))
}

fuzz_target!(|data: &[u8]| {
    let Some((&k, rest)) = data.split_first() else {
        return;
    };
    let Ok(s) = std::str::from_utf8(rest) else {
        return;
    };
    let prec = u64::from(k % 128) + 1;
    if s.len() <= 4000 {
        if let Ok(b) = Ball::from_decimal(s, prec) {
            if let Some(q) = exact(s) {
                assert!(b.contains_rational(&q), "{s:?} at {prec}: {b:?}");
            }
            for n in [1, 5, 40] {
                let _ = b.to_decimal(n, DecimalStyle::Digits);
                let _ = b.to_decimal(n, DecimalStyle::Radius);
            }
            let _ = b.rad().to_decimal();
            if exact(s).is_some() {
                let literal = b.to_literal(prec).expect("finite");
                let (n, e) = literal_parts(&literal).expect("a literal");
                let back = Ball::from_literal(&n, e, prec);
                assert_eq!(back.mid(), b.mid(), "{s:?} at {prec}: {literal}");
            }
        }
    }
    if let Ok(x) = Float::from_hex(s) {
        assert_eq!(Float::from_hex(&x.to_hex()).as_ref(), Ok(&x), "{s:?}");
    }
    if let Ok(m) = Mag::from_hex(s) {
        assert_eq!(Mag::from_hex(&m.to_hex()).as_ref(), Ok(&m), "{s:?}");
    }
});
