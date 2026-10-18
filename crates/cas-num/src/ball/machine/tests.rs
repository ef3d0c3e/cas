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
fn overflowing_check_product() {
    // q = MAX / b rounds so that q b overflows: not exact.
    let a = MachineBall::exact(f64::MAX).unwrap();
    let b = MachineBall::exact(5.604_110_778_259_519_4e17).unwrap();
    let q = a.div(&b).unwrap();
    assert!(!q.is_exact());
    assert!(holds(
        &q,
        &(&rat(f64::MAX) / &rat(5.604_110_778_259_519_4e17))
    ));
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

#[test]
fn functions_through_ball() {
    let x = MachineBall::new(0.75, 1e-12).unwrap();
    assert!(x.exp().unwrap().contains(0.75f64.exp()));
    assert!(x.expm1().unwrap().contains(0.75f64.exp_m1()));
    assert!(x.log().unwrap().contains(0.75f64.ln()));
    assert!(x.log1p().unwrap().contains(0.75f64.ln_1p()));
    let tiny = MachineBall::exact(1e-300).unwrap();
    assert!(tiny.expm1().unwrap().contains(1e-300));
    assert_eq!(MachineBall::exact(-1.0).unwrap().log(), None);
    assert_eq!(MachineBall::exact(1000.0).unwrap().exp(), None);
}

#[test]
fn trig_through_ball() {
    let x = MachineBall::new(0.75, 1e-12).unwrap();
    assert!(x.sin().unwrap().contains(0.75f64.sin()));
    assert!(x.cos().unwrap().contains(0.75f64.cos()));
    assert!(x.tan().unwrap().contains(0.75f64.tan()));
    assert!(x.atan().unwrap().contains(0.75f64.atan()));
    assert!(x.asin().unwrap().contains(0.75f64.asin()));
    assert!(x.acos().unwrap().contains(0.75f64.acos()));
    assert_eq!(MachineBall::exact(2.0).unwrap().asin(), None);
}

// ---- Double-double kernels (CAS-88).

mod kernels_against_mpfr {
    use proptest::prelude::*;
    use rug::Float as MpFloat;

    use super::super::MachineBall;
    use super::super::dd::Dd;
    use super::super::kernels::{self, KERNEL_ERR};
    use crate::ball::Ball;
    use crate::testkit::ToRug;

    fn value(d: Dd) -> MpFloat {
        MpFloat::with_val(400, d.hi) + d.lo
    }

    /// `|got - want| / |want|`, with `want` from MPFR at 400 bits.
    fn rel(got: Dd, want: &MpFloat) -> f64 {
        ((value(got) - want).abs() / want.clone().abs()).to_f64()
    }

    /// The kernels' budgets are below `2^-94`; check `2^-92` (the claim is
    /// `2^-90`).
    const CHECK: f64 = KERNEL_ERR / 4.0;

    fn mp(x: f64) -> MpFloat {
        MpFloat::with_val(400, x)
    }

    type F = (
        &'static str,
        fn(&MachineBall) -> Option<MachineBall>,
        fn(MpFloat) -> MpFloat,
    );

    proptest! {
        #[test]
        fn exp_within_budget(x in prop_oneof![-665.0f64..708.0, -700.0f64..-640.0, -1e-20f64..1e-20]) {
            if x < -665.0 { prop_assert!(kernels::exp(x).is_none()); return Ok(()); }
            let v = kernels::exp(x).unwrap();
            prop_assert!(rel(v, &mp(x).exp()) <= CHECK, "exp({}) off by {:e}", x, rel(v, &mp(x).exp()));
        }

        #[test]
        fn log_within_budget(e in -1022i32..1023, f in 1.0f64..2.0) {
            let x = f * 2f64.powi(e);
            if !x.is_normal() { return Ok(()); }
            let v = kernels::log(x).unwrap();
            let want = mp(x).ln();
            if want.is_zero() { prop_assert_eq!(v.to_f64(), 0.0); return Ok(()); }
            prop_assert!(rel(v, &want) <= CHECK, "log({}) off by {:e}", x, rel(v, &want));
        }

        #[test]
        fn sin_cos_within_budget(x in prop_oneof![-4.0f64..4.0, -1e6f64..1e6, -1e-30f64..1e-30, -1e-300f64..1e-300]) {
            if let Some((s, c)) = kernels::sin_cos(x) {
                let (ws, wc) = (mp(x).sin(), mp(x).cos());
                if !ws.is_zero() { prop_assert!(rel(s, &ws) <= CHECK, "sin({}) off by {:e}", x, rel(s, &ws)); }
                prop_assert!(rel(c, &wc) <= CHECK, "cos({}) off by {:e}", x, rel(c, &wc));
            }
        }

        #[test]
        fn atan_within_budget(x in prop_oneof![-2.0f64..2.0, -1e100f64..1e100, -1e-100f64..1e-100]) {
            let v = kernels::atan(x).unwrap();
            let want = mp(x).atan();
            if want.is_zero() { return Ok(()); }
            prop_assert!(rel(v, &want) <= CHECK, "atan({}) off by {:e}", x, rel(v, &want));
        }

        #[test]
        fn machine_functions_contain(m in prop_oneof![-30.0f64..30.0, -1e-3f64..1e-3], k in 0i32..60, z in any::<bool>()) {
            let r = if z { 0.0 } else { m.abs() * 2f64.powi(-k - 20) };
            let Some(b) = MachineBall::new(m, r) else { return Ok(()) };
            let fs: [F; 6] = [
                ("exp", MachineBall::exp, MpFloat::exp),
                ("log", MachineBall::log, MpFloat::ln),
                ("sin", MachineBall::sin, MpFloat::sin),
                ("cos", MachineBall::cos, MpFloat::cos),
                ("tan", MachineBall::tan, MpFloat::tan),
                ("atan", MachineBall::atan, MpFloat::atan),
            ];
            for (name, f, oracle) in fs {
                let Some(out) = f(&b) else {
                    prop_assert!(name == "log" && m - r <= 0.0 || name == "tan", "{} of {:?} fell through", name, b);
                    continue;
                };
                for end in [MpFloat::with_val(400, m) - r, MpFloat::with_val(400, m) + r] {
                    if name == "log" && end <= 0 { continue; }
                    let want = oracle(end);
                    let lo = MpFloat::with_val(400, out.mid()) - out.rad();
                    let hi = MpFloat::with_val(400, out.mid()) + out.rad();
                    prop_assert!(lo <= want && want <= hi, "{}({:?}) = {:?} misses {}", name, b, out, want);
                }
            }
        }
    }

    #[test]
    fn hard_inputs_and_fallbacks() {
        let contains = |b: &MachineBall, want: &MpFloat| {
            let lo = MpFloat::with_val(400, b.mid()) - b.rad();
            let hi = MpFloat::with_val(400, b.mid()) + b.rad();
            lo <= *want && *want <= hi
        };
        // Doubles next to multiples of π/2, up to the reduction's range and
        // beyond it: the kernel declines the closest ones, the Ball fallback
        // answers.
        for k in [
            1i64,
            2,
            3,
            7,
            1000,
            355,
            103_993,
            1 << 19,
            (1 << 20) + 1,
            1 << 40,
        ] {
            #[allow(clippy::cast_precision_loss, reason = "test points")]
            let x = (MpFloat::with_val(200, rug::float::Constant::Pi) * k / 2u32).to_f64();
            let b = MachineBall::exact(x).unwrap();
            let (s, c) = b.sin_cos().unwrap();
            assert!(
                contains(&s, &mp(x).sin()) && contains(&c, &mp(x).cos()),
                "sin_cos({x}): {s:?} {c:?}"
            );
        }
        // exp at the ends of its range and beyond (fallback, then None).
        for x in [-665.0, -690.4, -707.0, 708.0, -745.0, 709.7] {
            let b = MachineBall::exact(x).unwrap();
            if let Some(e) = b.exp() {
                assert!(contains(&e, &mp(x).exp()), "exp({x}) = {e:?}");
            }
        }
        assert_eq!(MachineBall::exact(800.0).unwrap().exp(), None);
        // log near 1, of the smallest normal and of the largest double.
        for x in [
            1.0 + f64::EPSILON,
            1.0 - f64::EPSILON / 2.0,
            f64::MIN_POSITIVE,
            f64::MAX,
            0.75,
        ] {
            let l = MachineBall::exact(x).unwrap().log().unwrap();
            assert!(contains(&l, &mp(x).ln()), "log({x}) = {l:?}");
        }
        // atan of tiny and huge arguments (beyond the kernel: Ball).
        for x in [
            1e-200,
            1e-310_f64.max(f64::MIN_POSITIVE),
            1e300,
            -1e150,
            0.5,
            1.0,
            -1.0,
        ] {
            let a = MachineBall::exact(x).unwrap().atan().unwrap();
            assert!(contains(&a, &mp(x).atan()), "atan({x}) = {a:?}");
        }
        // Wide balls take the Ball path.
        let wide = MachineBall::new(1.0, 0.5).unwrap();
        assert!(
            wide.exp().unwrap().contains(1.5f64.exp())
                && wide.exp().unwrap().contains(0.5f64.exp())
        );
        // Tight: one kernel evaluation costs about half an ulp.
        let e = MachineBall::exact(1.0).unwrap().exp().unwrap();
        assert!(e.rad() <= 1.0001 * f64::EPSILON * 1.36, "{e:?}");
    }

    #[test]
    fn constants_match_cas78() {
        let (ln2, half_pi) = kernels::constants_dd();
        let check = |d: Dd, b: Ball| {
            let want = b.mid().to_rational().to_rug();
            let got = MpFloat::with_val(400, d.hi) + d.lo;
            let err = (got - MpFloat::with_val(400, &want)).abs() / MpFloat::with_val(400, &want);
            assert!(err.to_f64() < 2f64.powi(-104), "{d:?}");
        };
        check(ln2, Ball::const_log2(300));
        check(half_pi, Ball::const_pi(300).mul_2exp(&-crate::Integer::ONE));
        let (e1, e2) = kernels::reduction_errors();
        assert!(
            e1 < 2f64.powi(-130) && e2 < 2f64.powi(-170),
            "{e1:e} {e2:e}"
        );
    }
}
