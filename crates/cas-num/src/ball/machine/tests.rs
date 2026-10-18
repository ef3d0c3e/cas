//! Machine balls against exact `Rational` evaluation: every result must
//! contain the exact result at the corners and midpoints of the inputs,
//! one exact operation costs at most half an ulp, the predicates agree with
//! exact comparisons, and the range checks send results back to `Ball`.
#![allow(
    clippy::many_single_char_names,
    clippy::float_cmp,
    reason = "operands a, b, x, y, r, s; exact float comparisons are the point"
)]

use std::cmp::Ordering;

use proptest::prelude::*;

use super::{EFT_MIN, MachineBall, half_ulp, next_down, next_up};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Rational};

fn rat(x: f64) -> Rational {
    Float::from_f64(x).unwrap().to_rational()
}

/// The exact interval of a machine ball.
fn interval(b: &MachineBall) -> (Rational, Rational) {
    let (m, r) = (rat(b.mid), rat(b.rad));
    (&m - &r, &m + &r)
}

fn holds(b: &MachineBall, v: &Rational) -> bool {
    let (lo, hi) = interval(b);
    lo <= *v && *v <= hi
}

/// Lower endpoint, midpoint and upper endpoint, exactly.
fn points(b: &MachineBall) -> [Rational; 3] {
    let (lo, hi) = interval(b);
    [lo, rat(b.mid), hi]
}

fn from_parts(sign: bool, biased: u64, frac: u64) -> f64 {
    f64::from_bits(u64::from(sign) << 63 | biased << 52 | (frac & ((1 << 52) - 1)))
}

/// Regular doubles: small integers, all-ones and power-of-two mantissas, and
/// random ones with exponents around 1, near overflow and near the subnormal
/// range.
fn double() -> BoxedStrategy<f64> {
    let frac = prop_oneof![Just(0u64), Just(u64::MAX), Just(1u64), any::<u64>()];
    let around = |lo: i64, hi: i64| {
        (any::<bool>(), lo..=hi, frac.clone())
            .prop_map(|(s, e, f)| from_parts(s, u64::try_from(e + 1023).unwrap(), f))
    };
    prop_oneof![
        2 => (-1000i32..=1000).prop_map(f64::from),
        4 => around(-60, 60),
        1 => around(990, 1023),
        1 => around(-1022, -960),
        1 => Just(0.0),
    ]
    .boxed()
}

/// A radius for midpoint `m`: zero, relative (`|m| 2^-k`), or unrelated.
fn radius(m: f64) -> BoxedStrategy<f64> {
    prop_oneof![
        2 => Just(0.0),
        3 => (0i32..60).prop_map(move |k| m.abs() * 2f64.powi(-k)),
        1 => double().prop_map(f64::abs),
    ]
    .prop_filter_map("regular", move |r| MachineBall::new(m, r).map(|b| b.rad))
    .boxed()
}

fn ball() -> BoxedStrategy<MachineBall> {
    double()
        .prop_flat_map(|m| (Just(m), radius(m)))
        .prop_map(|(m, r)| MachineBall::new(m, r).unwrap())
        .boxed()
}

/// Exact balls and radius at most `|m| / 2`, with exponents around 1: no
/// operation of these may fall back.
fn ordinary() -> BoxedStrategy<MachineBall> {
    let frac = prop_oneof![Just(0u64), Just(u64::MAX), any::<u64>()];
    (
        any::<bool>(),
        1023 - 60..=1023 + 60u64,
        frac,
        prop_oneof![Just(None), (1i32..60).prop_map(Some)],
    )
        .prop_map(|(s, e, f, k)| {
            let m = from_parts(s, e, f);
            MachineBall::new(m, k.map_or(0.0, |k| m.abs() * 2f64.powi(-k))).unwrap()
        })
        .boxed()
}

type Op = (
    &'static str,
    fn(&MachineBall, &MachineBall) -> Option<MachineBall>,
    fn(&Rational, &Rational) -> Option<Rational>,
);

const OPS: [Op; 4] = [
    ("add", MachineBall::add, |x, y| Some(x + y)),
    ("sub", MachineBall::sub, |x, y| Some(x - y)),
    ("mul", MachineBall::mul, |x, y| Some(x * y)),
    ("div", MachineBall::div, |x, y| x.checked_div(y)),
];

fn check_binary(a: &MachineBall, b: &MachineBall) -> Result<(), TestCaseError> {
    for (name, op, exact) in OPS {
        let Some(c) = op(a, b) else { continue };
        prop_assert!(
            MachineBall::new(c.mid, c.rad).is_some(),
            "{name}: irregular {c:?}"
        );
        for x in &points(a) {
            for y in &points(b) {
                let v = exact(x, y);
                prop_assert!(v.is_some(), "{name}: {b:?} contains 0 but gave {c:?}");
                let v = v.unwrap();
                prop_assert!(holds(&c, &v), "{name}({a:?}, {b:?}) = {c:?} misses {v}");
            }
        }
        if a.is_exact() && b.is_exact() && c.mid != 0.0 {
            // One rounding: at most half an ulp (or the smallest normal).
            prop_assert!(
                c.rad <= half_ulp(c.mid).max(f64::MIN_POSITIVE),
                "{name}: loose {c:?}"
            );
        }
    }
    Ok(())
}

fn check_sqrt(a: &MachineBall) -> Result<(), TestCaseError> {
    let Some(c) = a.sqrt() else {
        prop_assert!(a.mid < a.rad, "sqrt({a:?}) fell back");
        return Ok(());
    };
    let (lo, hi) = interval(&c);
    for x in &points(a) {
        prop_assert!(*x >= Rational::from(0i64), "sqrt of {a:?} reaches below 0");
        let below = lo <= Rational::from(0i64) || &lo * &lo <= *x;
        prop_assert!(
            below && *x <= &hi * &hi,
            "sqrt({a:?}) = {c:?} misses sqrt({x})"
        );
    }
    if a.is_exact() && c.mid != 0.0 {
        prop_assert!(
            c.rad <= half_ulp(c.mid).max(f64::MIN_POSITIVE),
            "sqrt: loose {c:?}"
        );
    }
    Ok(())
}

/// `x^n` for a rational `x != 0` or `n >= 0`.
fn rat_pow(x: &Rational, n: i32) -> Rational {
    let mut p = Rational::from(1i64);
    for _ in 0..n.unsigned_abs() {
        p = &p * x;
    }
    if n < 0 { Rational::from(1i64) / p } else { p }
}

proptest! {
    #[test]
    fn arithmetic_contains_exact((a, b) in (ball(), ball())) {
        check_binary(&a, &b)?;
        check_binary(&b, &a)?;
        check_sqrt(&a)?;
        check_sqrt(&a.abs())?;
    }

    #[test]
    fn ordinary_never_falls_back((a, b) in (ordinary(), ordinary())) {
        prop_assert!(a.add(&b).is_some() && a.sub(&b).is_some());
        prop_assert!(a.mul(&b).is_some() && a.div(&b).is_some());
        prop_assert!(a.abs().sqrt().is_some());
        check_binary(&a, &b)?;
    }

    #[test]
    fn powers((a, n) in (ball(), -6i32..=9)) {
        if let Some(c) = a.pow(n) {
            for x in &points(&a) {
                prop_assert!(n >= 0 || *x != Rational::from(0i64), "{a:?}^{n} through 0");
                let v = rat_pow(x, n);
                prop_assert!(holds(&c, &v), "{a:?}^{n} = {c:?} misses {v}");
            }
        }
    }

    #[test]
    fn contains_is_exact((a, t) in (ball(), any::<u8>())) {
        let edge = a.mid + a.rad;
        let cand = [a.mid, edge, next_up(edge), next_down(edge), a.mid - a.rad,
            next_down(a.mid - a.rad), f64::from(t) - 128.0];
        for x in cand {
            if x.is_finite() {
                prop_assert_eq!(a.contains(x), holds(&a, &rat(x)), "{:?} ∋ {}", a, x);
            }
        }
    }

    #[test]
    fn predicates_are_exact((a, b, k) in (ball(), ball(), 0u8..4)) {
        // Touching and nearly touching pairs as well as random ones.
        let b = match k {
            0 => b,
            1 => MachineBall::new(a.mid + a.rad + b.rad, b.rad).unwrap_or(b),
            2 => MachineBall::new(next_up(a.mid + a.rad + b.rad), b.rad).unwrap_or(b),
            _ => MachineBall::new(a.mid + (a.rad - b.rad), b.rad).unwrap_or(b),
        };
        let (alo, ahi) = interval(&a);
        let (blo, bhi) = interval(&b);
        prop_assert_eq!(a.overlaps(&b), alo <= bhi && blo <= ahi, "{:?} {:?}", a, b);
        prop_assert_eq!(a.contains_ball(&b), alo <= blo && bhi <= ahi, "{:?} ⊇ {:?}", a, b);
        let zero = Rational::from(0i64);
        let sign = if alo > zero {
            Some(Ordering::Greater)
        } else if ahi < zero {
            Some(Ordering::Less)
        } else if a.mid == 0.0 && a.rad == 0.0 {
            Some(Ordering::Equal)
        } else {
            None
        };
        prop_assert_eq!(a.sign(), sign);
        let cmp = if ahi < blo {
            Some(Ordering::Less)
        } else if bhi < alo {
            Some(Ordering::Greater)
        } else if a.is_exact() && b.is_exact() && a.mid == b.mid {
            Some(Ordering::Equal)
        } else {
            None
        };
        prop_assert_eq!(a.cmp(&b), cmp, "{:?} vs {:?}", a, b);
    }

    #[test]
    fn conversions(a in ball(), (num, den, prec, k) in (any::<i64>(), 1i64..=i64::MAX, 1u64..200, 0i64..200)) {
        let b = a.to_ball();
        for x in &points(&a) {
            prop_assert!(b.contains_rational(x));
        }
        // Back again, unless the Mag rounded the radius up past f64::MAX.
        match MachineBall::from_ball(&b) {
            Some(c) => prop_assert!(c.mid == a.mid && c.rad >= a.rad),
            None => prop_assert!(a.rad > f64::MAX / 2.0),
        }
        // A Ball at up to 200 bits, with a radius, into a machine ball.
        let q = Rational::new(Integer::from(num), Integer::from(den));
        let rad = Mag::power_of_two(&Integer::from(-k));
        let wide = Ball::from_rational(&q, prec).add_error(&rad);
        if let Some(c) = MachineBall::from_ball(&wide) {
            let r = rad_value(wide.rad());
            let m = wide.mid().to_rational();
            prop_assert!(holds(&c, &(&m - &r)) && holds(&c, &(&m + &r)), "{wide:?} → {c:?}");
        }
    }
}

fn rad_value(m: &Mag) -> Rational {
    m.to_float().unwrap().to_rational()
}

#[test]
fn exact_results_stay_exact() {
    let x = |v: f64| MachineBall::exact(v).unwrap();
    assert_eq!(x(3.0).add(&x(4.0)), Some(x(7.0)));
    assert_eq!(x(3.0).mul(&x(7.0)), Some(x(21.0)));
    assert_eq!(x(21.0).div(&x(7.0)), Some(x(3.0)));
    assert_eq!(x(1.0).div(&x(4.0)), Some(x(0.25)));
    assert_eq!(x(49.0).sqrt(), Some(x(7.0)));
    assert_eq!(x(2.0).pow(10), Some(x(1024.0)));
    assert_eq!(x(2.0).pow(-3), Some(x(0.125)));
    assert_eq!(x(0.0).sqrt(), Some(x(0.0)));
    let third = x(1.0).div(&x(3.0)).unwrap();
    assert!(!third.is_exact() && third.rad <= half_ulp(third.mid));
    // 0.1 + 0.2: the exact error of the sum, not half an ulp.
    let s = x(0.1).add(&x(0.2)).unwrap();
    assert!(s.rad > 0.0 && holds(&s, &(&rat(0.1) + &rat(0.2))));
}

#[test]
fn falls_back_outside_the_range() {
    let x = |v: f64| MachineBall::exact(v).unwrap();
    // Overflow and underflow of midpoints.
    assert_eq!(x(f64::MAX).add(&x(f64::MAX)), None);
    assert_eq!(x(f64::MAX).mul(&x(2.0)), None);
    assert_eq!(x(1e-300).mul(&x(1e-300)), None);
    assert_eq!(x(1e-300).div(&x(1e300)), None);
    assert_eq!(x(f64::MIN_POSITIVE).mul(&x(0.5)), None);
    // Radius overflow.
    let wide = MachineBall::new(1.0, f64::MAX).unwrap();
    assert_eq!(wide.mul(&MachineBall::new(4.0, 2.0).unwrap()), None);
    // Subnormal, infinite and NaN parts are refused.
    assert_eq!(MachineBall::new(f64::from_bits(1), 0.0), None);
    assert_eq!(MachineBall::new(1.0, f64::from_bits(1)), None);
    assert_eq!(MachineBall::new(f64::INFINITY, 0.0), None);
    assert_eq!(MachineBall::new(1.0, f64::NAN), None);
    assert_eq!(MachineBall::new(1.0, -1.0), None);
    // Division by a ball containing 0, the root of one reaching below 0.
    let around1 = MachineBall::new(1.0, 2.0).unwrap();
    assert_eq!(x(1.0).div(&around1), None);
    assert_eq!(around1.inv(), None);
    assert_eq!(around1.sqrt(), None);
    assert_eq!(x(-1.0).sqrt(), None);
    assert_eq!(MachineBall::new(0.0, 1.0).unwrap().sqrt(), None);
    // A Ball out of range, or indeterminate.
    let tiny = Ball::exact(Float::from_integer_exp(
        &Integer::ONE,
        &Integer::from(-2000),
    ));
    assert_eq!(MachineBall::from_ball(&tiny), None);
    assert_eq!(MachineBall::from_ball(&Ball::INDETERMINATE), None);
    let small_rad =
        Ball::exact(Float::from(1i64)).add_error(&Mag::power_of_two(&Integer::from(-1100)));
    // A radius below the normal range is rounded up instead.
    assert_eq!(
        MachineBall::from_ball(&small_rad),
        MachineBall::new(1.0, f64::MIN_POSITIVE)
    );
    // sqrt(1e16 ± 4 MIN_POSITIVE): radius about 2 MIN_POSITIVE / 1e8.
    let exact_root = MachineBall::new(1e16, 4.0 * f64::MIN_POSITIVE)
        .unwrap()
        .sqrt()
        .unwrap();
    assert!(exact_root.mid == 1e8 && exact_root.rad == f64::MIN_POSITIVE);
    // The caller's fallback: the same op through Ball contains the result.
    let (a, b) = (
        Float::from_f64(1e-300).unwrap(),
        Float::from_f64(1e-300).unwrap(),
    );
    let full = Ball::exact(a).mul(&Ball::exact(b), 53);
    assert!(full.contains_rational(&(&rat(1e-300) * &rat(1e-300))));
}

#[test]
fn half_ulp_values() {
    assert_eq!(half_ulp(1.0), 2f64.powi(-53));
    assert_eq!(half_ulp(1.5), 2f64.powi(-53));
    assert_eq!(half_ulp(-2.0), 2f64.powi(-52));
    assert_eq!(half_ulp(f64::MAX), 2f64.powi(970));
    assert_eq!(half_ulp(f64::MIN_POSITIVE), f64::from_bits(1));
    assert_eq!(half_ulp(2f64.powi(-970)), 2f64.powi(-1023));
    assert_eq!(EFT_MIN, 2f64.powi(-900));
    assert_eq!(next_up(0.0), f64::from_bits(1));
    assert_eq!(next_up(-f64::from_bits(1)), -0.0);
    assert_eq!(next_down(1.0), 1.0 - 2f64.powi(-53));
}

#[test]
fn via_ball_fallback() {
    let a = MachineBall::new(2.0, 1e-10).unwrap();
    let viab = a.via_ball(Ball::sqrt).unwrap();
    let native = a.sqrt().unwrap();
    assert!(viab.overlaps(&native));
    assert!(viab.contains(2f64.sqrt()));
    let pi = MachineBall::exact(1.0)
        .unwrap()
        .via_ball(|_, prec| Ball::const_pi(prec))
        .unwrap();
    assert!(pi.contains(std::f64::consts::PI));
    assert!(pi.rad <= 2f64.powi(-50));
}
