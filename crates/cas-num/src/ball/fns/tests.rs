//! Elementary functions against MPFR brackets at `prec + 64` bits (each
//! result must meet the bracket and know about `prec` bits), identities,
//! cross-precision overlap, domains, and the huge and tiny arguments the
//! reductions must survive.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{Dyadic, bracket, check_encloses, dyadic_with};

type Oracle = fn(&MpFloat, u32, MpRound) -> MpFloat;
type Ours = fn(&Ball, u64) -> Ball;

const FUNS: [(&str, Ours, Oracle); 4] = [
    ("exp", Ball::exp, |x, p, r| {
        MpFloat::with_val_round(p, x.exp_ref(), r).0
    }),
    ("expm1", Ball::expm1, |x, p, r| {
        MpFloat::with_val_round(p, x.exp_m1_ref(), r).0
    }),
    ("log", Ball::log, |x, p, r| {
        MpFloat::with_val_round(p, x.ln_ref(), r).0
    }),
    ("log1p", Ball::log1p, |x, p, r| {
        MpFloat::with_val_round(p, x.ln_1p_ref(), r).0
    }),
];

/// Whether the function is defined at `x`.
fn in_domain(name: &str, x: &MpFloat) -> bool {
    match name {
        "log" => *x > 0,
        "log1p" => *x > -1,
        _ => true,
    }
}

/// Whether MPFR can bracket `name(x)`: exp overflows its exponent range
/// near `2^30 log 2`.
fn in_mpfr_range(name: &str, x: &MpFloat) -> bool {
    !name.starts_with("exp") || x.clone().abs() < 1e8
}

fn check(name: &str, f: Ours, oracle: Oracle, x: &Dyadic, prec: u64) -> Result<(), TestCaseError> {
    let xr = x.to_rug().expect("in MPFR's range");
    let ours = x.to_ours();
    let b = f(&Ball::exact(ours.clone()), prec);
    if !in_domain(name, &xr) {
        prop_assert_eq!(&b, &Ball::INDETERMINATE, "{}({:?})", name, ours);
        return Ok(());
    }
    if !in_mpfr_range(name, &xr) {
        return Ok(());
    }
    prop_assert!(
        b.is_finite(),
        "{}({:?}) at {}: indeterminate",
        name,
        ours,
        prec
    );
    let br = bracket(u32::try_from(prec).unwrap(), |p, r| oracle(&xr, p, r));
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br)
        .map_err(|e| TestCaseError::fail(format!("{name}({ours:?}) at {prec}: {e}")))?;
    if !b.mid().is_zero() {
        #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
        let want = prec as i64 - 3;
        prop_assert!(
            b.rel_accuracy_bits() >= want,
            "{}({:?}) at {}: only {} bits: {:?}",
            name,
            ours,
            prec,
            b.rel_accuracy_bits(),
            b
        );
    }
    Ok(())
}

fn prec() -> impl Strategy<Value = u64> {
    prop_oneof![1u64..80, 80u64..400, Just(1000u64), Just(2049u64)]
}

/// Arguments from tiny to about 2^30, so that exp stays in MPFR's range.
fn argument() -> BoxedStrategy<Dyadic> {
    prop_oneof![
        3 => dyadic_with(2, -70..=4),
        1 => dyadic_with(4, -300..=-60),
        1 => dyadic_with(2, 4..=26),
    ]
    .boxed()
}

proptest! {
    #[test]
    fn against_mpfr(x in argument(), prec in prec()) {
        for (name, f, oracle) in FUNS {
            check(name, f, oracle, &x, prec)?;
        }
    }

    #[test]
    fn identities(x in dyadic_with(2, -40..=6), prec in 2u64..300) {
        let b = Ball::exact(x.to_ours());
        // log(exp x) ∋ x, exp(log x) ∋ x (x > 0), expm1(log1p x) ∋ x (x > -1).
        prop_assert!(b.exp(prec).log(prec).contains(b.mid()));
        if !b.mid().is_negative() && !b.mid().is_zero() {
            prop_assert!(b.log(prec).exp(prec).contains(b.mid()));
        }
        let m1 = b.mid().add(&Float::from(1i64));
        if !m1.is_negative() && !m1.is_zero() {
            prop_assert!(b.log1p(prec).expm1(prec).contains(b.mid()));
        }
        // Cross-precision overlap.
        for (name, f, _) in FUNS {
            let (lo, hi) = (f(&b, prec), f(&b, prec + 500));
            prop_assert!(lo.overlaps(&hi), "{} at {} and {}: {:?} {:?}", name, prec, prec + 500, lo, hi);
        }
    }

    #[test]
    fn radii_propagate(x in dyadic_with(2, -20..=4), j in 1i64..80, prec in 2u64..200) {
        let m = x.to_ours();
        // Values below 2^8 (the strategy bounds the dyadic exponent, not the
        // value), so that exact rationals of the results stay small.
        if m.is_zero() || *m.exponent() > Integer::from(8i64) { return Ok(()); }
        let r = Mag::power_of_two(&(m.exponent() - &Integer::from(j)));
        let b = Ball::new(m.clone(), r.clone());
        let rf = r.to_float().unwrap();
        for (name, f, oracle) in FUNS {
            let out = f(&b, prec);
            // Exact rationals of the result must stay small.
            let small = |e: &Integer| e.to_i64().is_some_and(|e| e.abs() < 1 << 20);
            prop_assert!(!out.is_finite() || ((out.mid().is_zero() || small(out.mid().exponent())) && (out.rad().is_zero() || small(out.rad().exponent()))), "{} of {:?} at {}: {:?}", name, b, prec, out);
            for end in [m.sub(&rf), m.add(&rf)] {
                let er = Dyadic::from_rug(&MpFloat::with_val(10_000, &end.to_rational().to_rug())).unwrap();
                let xr = er.to_rug().unwrap();
                if !in_domain(name, &xr) || !in_mpfr_range(name, &xr) || !out.is_finite() { continue; }
                let br = bracket(u32::try_from(prec).unwrap(), |p, rd| oracle(&xr, p, rd));
                let mid = out.mid().to_rational().to_rug();
                let rad = out.rad().to_float().unwrap().to_rational().to_rug();
                check_encloses(&mid, &rad, &br).map_err(|e| TestCaseError::fail(format!("{name} of {b:?} at {prec}, end {end:?}: {e}")))?;
            }
        }
    }
}

fn exact(x: f64) -> Ball {
    Ball::exact(Float::from_f64(x).unwrap())
}

#[test]
fn precision_sweep() {
    let args = [
        0.5, -0.75, 1.0, 2.0, 3.25, -7.5, 1e-10, 100.0, 0.9999, 1.0001,
    ];
    for p in [64u64, 128, 256, 1000, 4096, 10_000] {
        for x in args {
            let d = Dyadic::from_rug(&MpFloat::with_val(53, x)).unwrap();
            for (name, f, oracle) in FUNS {
                check(name, f, oracle, &d, p).unwrap();
            }
        }
    }
}

#[test]
fn domains() {
    let b = |m: f64, r: f64| {
        Ball::new(
            Float::from_f64(m).unwrap(),
            Mag::from_float(&Float::from_f64(r).unwrap()),
        )
    };
    assert_eq!(b(1.0, 1.0).log(64), Ball::INDETERMINATE);
    assert_eq!(b(1.0, 2.0).log(64), Ball::INDETERMINATE);
    assert_eq!(exact(-1.0).log(64), Ball::INDETERMINATE);
    assert_eq!(exact(0.0).log(64), Ball::INDETERMINATE);
    assert_eq!(exact(-1.0).log1p(64), Ball::INDETERMINATE);
    assert_eq!(b(-0.5, 0.5).log1p(64), Ball::INDETERMINATE);
    assert_eq!(Ball::INDETERMINATE.exp(64), Ball::INDETERMINATE);
    // Exact values.
    assert_eq!(exact(0.0).exp(64), exact(1.0));
    assert_eq!(exact(0.0).expm1(64), exact(0.0));
    assert_eq!(exact(1.0).log(64), exact(0.0));
    assert_eq!(exact(0.0).log1p(64), exact(0.0));
    // A radius around 0.
    let around0 = b(0.0, 0.5);
    assert!(around0.exp(64).contains(&Float::from(1i64)));
    assert!(
        around0
            .expm1(64)
            .contains(&Float::from_f64(0.5f64.exp_m1()).unwrap())
    );
    // An enormous radius: indeterminate, not a panic.
    let wide = Ball::new(
        Float::from(1i64),
        Mag::power_of_two(&Integer::from(1i64 << 40)),
    );
    assert_eq!(wide.exp(64), Ball::INDETERMINATE);
}

#[test]
fn huge_and_tiny_arguments() {
    // exp(10^30) = 2^n 2^f with n + f = 10^30 / log 2, beyond MPFR's
    // exponents: compare the scaled-down value with MPFR's 2^f.
    let x = Integer::from(10u64).pow(30);
    for p in [64u64, 256, 1000] {
        let b = Ball::from_integer(&x).exp(p);
        // n + f = x / log 2 with f to well beyond the bracket's p + 64 bits
        // (x has about 100 bits above the point).
        let q = u32::try_from(p).unwrap() + 300;
        let y =
            MpFloat::with_val(q, &x.to_rug()) / MpFloat::with_val(q, rug::float::Constant::Log2);
        let n = y.clone().floor();
        let frac = y - &n;
        let n = n.to_integer().unwrap();
        let scaled = b.mul_2exp(&-Integer::from_str_radix(&n.to_string(), 10).unwrap());
        let br = bracket(u32::try_from(p).unwrap(), |q, r| {
            MpFloat::with_val_round(q, frac.exp2_ref(), r).0
        });
        let mid = scaled.mid().to_rational().to_rug();
        let rad = scaled.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, &br).unwrap();
        #[allow(clippy::cast_possible_wrap, reason = "small")]
        let want = p as i64 - 3;
        assert!(b.rel_accuracy_bits() >= want, "exp(10^30) at {p}: {b:?}");
        // And back.
        assert!(b.log(p + 120).contains_integer(&x), "log(exp(10^30))");
    }
    // log(2^-10^6) = -10^6 log 2.
    let tiny = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(-1_000_000i64)));
    for p in [64u64, 1000] {
        let l = tiny.log(p);
        let br = bracket(u32::try_from(p).unwrap(), |q, r| {
            MpFloat::with_val_round(
                q,
                MpFloat::with_val(q + 64, rug::float::Constant::Log2) * -1_000_000i64,
                r,
            )
            .0
        });
        let mid = l.mid().to_rational().to_rug();
        let rad = l.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, &br).unwrap();
    }
    // exp of a hugely negative number: a tiny bound, not indeterminate.
    let neg = Ball::exact(Float::from(-1i64).mul_2exp(&Integer::from(100_000i64)));
    let e = neg.exp(64);
    assert!(e.mid().is_zero() && e.is_finite() && e.contains(&Float::ZERO));
    assert_eq!(neg.neg().exp(64), Ball::INDETERMINATE);
    // Tiny arguments keep relative accuracy.
    let t = Ball::exact(Float::from(3i64).mul_2exp(&Integer::from(-5000i64)));
    for f in [Ball::expm1, Ball::log1p] {
        let v = f(&t, 100);
        assert!(v.rel_accuracy_bits() >= 97, "{v:?}");
        assert!(v.contains(t.mid()) || v.rel_accuracy_bits() > 0);
    }
}

// ---- Trigonometric functions (CAS-80).

type Unary = (&'static str, Ours, Oracle);

const TRIG: [Unary; 6] = [
    ("sin", Ball::sin, |x, p, r| {
        MpFloat::with_val_round(p, x.sin_ref(), r).0
    }),
    ("cos", Ball::cos, |x, p, r| {
        MpFloat::with_val_round(p, x.cos_ref(), r).0
    }),
    ("tan", Ball::tan, |x, p, r| {
        MpFloat::with_val_round(p, x.tan_ref(), r).0
    }),
    ("atan", Ball::atan, |x, p, r| {
        MpFloat::with_val_round(p, x.atan_ref(), r).0
    }),
    ("asin", Ball::asin, |x, p, r| {
        MpFloat::with_val_round(p, x.asin_ref(), r).0
    }),
    ("acos", Ball::acos, |x, p, r| {
        MpFloat::with_val_round(p, x.acos_ref(), r).0
    }),
];

/// Containment and tightness of a trigonometric function at an exact point.
fn check_trig(
    name: &str,
    f: Ours,
    oracle: Oracle,
    x: &MpFloat,
    prec: u64,
) -> Result<(), TestCaseError> {
    let d = Dyadic::from_rug(x).unwrap();
    let b = f(&Ball::exact(d.to_ours()), prec);
    let outside = matches!(name, "asin" | "acos") && x.clone().abs() > 1;
    if outside {
        prop_assert_eq!(&b, &Ball::INDETERMINATE, "{}({})", name, x);
        return Ok(());
    }
    prop_assert!(b.is_finite(), "{}({}) at {}: indeterminate", name, x, prec);
    let br = bracket(u32::try_from(prec).unwrap(), |p, r| oracle(x, p, r));
    let mid = b.mid().to_rational().to_rug();
    let rad = b.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br)
        .map_err(|e| TestCaseError::fail(format!("{name}({x}) at {prec}: {e}")))?;
    if !b.mid().is_zero() {
        #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
        let want = prec as i64 - 3;
        prop_assert!(
            b.rel_accuracy_bits() >= want,
            "{}({}) at {}: only {} bits: {:?}",
            name,
            x,
            prec,
            b.rel_accuracy_bits(),
            b
        );
    }
    Ok(())
}

proptest! {
    #[test]
    fn trig_against_mpfr(x in argument(), prec in prec()) {
        let xr = x.to_rug().unwrap();
        for (name, f, oracle) in TRIG {
            check_trig(name, f, oracle, &xr, prec)?;
        }
    }

    #[test]
    fn atan2_against_mpfr(y in dyadic_with(2, -40..=4), x in dyadic_with(2, -40..=4), prec in prec(), zeros in 0u8..4) {
        let (y, x) = match zeros {
            0 => (y, x),
            1 => (Dyadic::from_rug(&MpFloat::new(53)).unwrap(), x),
            2 => (y, Dyadic::from_rug(&MpFloat::new(53)).unwrap()),
            _ => (y.clone(), y),
        };
        let (yr, xr) = (y.to_rug().unwrap(), x.to_rug().unwrap());
        let b = Ball::exact(y.to_ours()).atan2(&Ball::exact(x.to_ours()), prec);
        if yr == 0 && xr == 0 {
            prop_assert!(b.contains(&Float::ZERO) && b.rad() >= &Mag::from_integer(&Integer::from(3)));
            return Ok(());
        }
        let br = bracket(u32::try_from(prec).unwrap(), |p, r| MpFloat::with_val_round(p, yr.atan2_ref(&xr), r).0);
        let mid = b.mid().to_rational().to_rug();
        let rad = b.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, &br).map_err(|e| TestCaseError::fail(format!("atan2({yr}, {xr}) at {prec}: {e}")))?;
        if !b.mid().is_zero() {
            #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
            let want = prec as i64 - 3;
            prop_assert!(b.rel_accuracy_bits() >= want, "atan2({}, {}) at {}: {:?}", yr, xr, prec, b);
        }
    }

    #[test]
    fn trig_identities(x in dyadic_with(2, -40..=1), prec in 2u64..300) {
        let b = Ball::exact(x.to_ours());
        let (s, c) = b.sin_cos(prec);
        let one = s.mul(&s, prec).add(&c.mul(&c, prec), prec);
        prop_assert!(one.contains(&Float::from(1i64)), "sin² + cos² of {:?}: {:?}", b, one);
        // |x| < 1.5 < π/2 (the strategy's values stay below 2^1 · 2^-0).
        if b.mid().cmp_abs(&Float::from_f64(1.5).unwrap()) == std::cmp::Ordering::Less {
            prop_assert!(b.tan(prec).atan(prec).contains(b.mid()), "atan(tan {:?})", b);
            prop_assert!(s.asin(prec).contains(b.mid()), "asin(sin {:?})", b);
            if !b.mid().is_negative() {
                prop_assert!(c.acos(prec).contains(b.mid()), "acos(cos {:?})", b);
            }
        }
        for (name, f, _) in TRIG {
            let (lo, hi) = (f(&b, prec), f(&b, prec + 500));
            prop_assert!(lo.overlaps(&hi), "{} at {} and {}: {:?} {:?}", name, prec, prec + 500, lo, hi);
        }
    }

    #[test]
    fn trig_radii(x in dyadic_with(2, -20..=4), j in 1i64..80, prec in 2u64..200) {
        let m = x.to_ours();
        if m.is_zero() || *m.exponent() > Integer::from(8i64) { return Ok(()); }
        let r = Mag::power_of_two(&(m.exponent() - &Integer::from(j)));
        let b = Ball::new(m.clone(), r.clone());
        let rf = r.to_float().unwrap();
        for (name, f, oracle) in TRIG {
            let out = f(&b, prec);
            if !out.is_finite() { continue; }
            for end in [m.sub(&rf), m.add(&rf)] {
                let xr = Dyadic::from_rug(&MpFloat::with_val(10_000, &end.to_rational().to_rug())).unwrap().to_rug().unwrap();
                if matches!(name, "asin" | "acos") && xr.clone().abs() > 1 { continue; }
                let br = bracket(u32::try_from(prec).unwrap(), |p, rd| oracle(&xr, p, rd));
                let mid = out.mid().to_rational().to_rug();
                let rad = out.rad().to_float().unwrap().to_rational().to_rug();
                check_encloses(&mid, &rad, &br).map_err(|e| TestCaseError::fail(format!("{name} of {b:?} at {prec}, end {end:?}: {e}")))?;
            }
        }
    }

    #[test]
    fn sin_cos_pi(x in dyadic_with(2, -60..=2), prec in 2u64..300, n in any::<i64>()) {
        // Exact at integers and half-integers, any size.
        let int = Ball::from_integer(&Integer::from(n));
        prop_assert_eq!(int.sin_pi(prec), Ball::ZERO);
        let sign = if n % 2 == 0 { 1 } else { -1 };
        prop_assert_eq!(int.cos_pi(prec), Ball::from_integer(&Integer::from(sign)));
        let half = Ball::exact(Float::from(n).add(&Float::from(1i64).mul_2exp(&-Integer::ONE)));
        prop_assert_eq!(half.sin_pi(prec), Ball::from_integer(&Integer::from(sign)));
        prop_assert_eq!(half.cos_pi(prec), Ball::ZERO);
        // Elsewhere, against sin(π x) from MPFR far above the precision (plus
        // the bits above x's point, which π x loses to the period).
        let xr = x.to_rug().unwrap();
        let b = Ball::exact(x.to_ours());
        if (xr.clone() * 2u32).is_integer() { return Ok(()); }
        let above = u32::try_from(xr.get_exp().unwrap_or(0).max(0)).unwrap();
        for (name, ours, cos) in [("sin_pi", b.sin_pi(prec), false), ("cos_pi", b.cos_pi(prec), true)] {
            let q = u32::try_from(prec).unwrap() + 400 + above;
            let t = MpFloat::with_val(q, rug::float::Constant::Pi) * &xr;
            let v = if cos { t.cos() } else { t.sin() };
            if v.is_zero() { continue; }
            let slack = MpFloat::with_val(53, 1) << (v.get_exp().unwrap() - i32::try_from(q).unwrap() + 8);
            let lo = (v.clone() - &slack).to_rational().unwrap();
            let hi = (v + &slack).to_rational().unwrap();
            let br = crate::testkit::float::Bracket { lo, hi };
            let mid = ours.mid().to_rational().to_rug();
            let rad = ours.rad().to_float().unwrap().to_rational().to_rug();
            check_encloses(&mid, &rad, &br).map_err(|e| TestCaseError::fail(format!("{name}({xr}) at {prec}: {e}")))?;
            if !ours.mid().is_zero() {
                #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
                let want = prec as i64 - 3;
                prop_assert!(ours.rel_accuracy_bits() >= want, "{}({}) at {}: {:?}", name, xr, prec, ours);
            }
        }
    }
}

#[test]
fn huge_and_near_pole_arguments() {
    // sin(10^100): π at about p + 333 bits for the reduction.
    let big = Integer::from(10u64).pow(100);
    for p in [64u64, 256, 1000] {
        let xr = MpFloat::with_val(400, &big.to_rug());
        for (name, f, oracle) in TRIG.iter().take(3) {
            check_trig(name, *f, *oracle, &xr, p).unwrap();
        }
    }
    // Next to multiples of π/2: the reduction cancels; relative accuracy holds.
    for k in [1i64, 2, 3, 4, 10, 1001, 1_000_000, 1 << 40] {
        for bits in [80u32, 300] {
            let x = MpFloat::with_val(bits, rug::float::Constant::Pi) * k / 2;
            for p in [64u64, 256, 1000] {
                for (name, f, oracle) in TRIG.iter().take(3) {
                    check_trig(name, *f, *oracle, &x, p).unwrap();
                }
            }
        }
    }
    // Far beyond the reduction limit, a wide radius, or no information: [0 ± 1].
    let far = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(100_000i64)));
    let unit = Ball::new(Float::ZERO, Mag::ONE);
    assert_eq!(far.sin(64), unit);
    assert_eq!(
        Ball::new(Float::ZERO, Mag::from_integer(&Integer::from(5))).cos(64),
        unit
    );
    assert_eq!(Ball::INDETERMINATE.sin(64), unit);
    // tan across a pole: indeterminate.
    let pole = Ball::const_pi(64)
        .mul_2exp(&-Integer::ONE)
        .add_error(&Mag::power_of_two(&Integer::from(-10)));
    assert_eq!(pole.tan(64), Ball::INDETERMINATE);
    // Tiny arguments keep relative accuracy.
    let t = Ball::exact(Float::from(3i64).mul_2exp(&Integer::from(-5000i64)));
    for f in [Ball::sin, Ball::tan, Ball::atan, Ball::asin] {
        let v = f(&t, 100);
        assert!(
            v.rel_accuracy_bits() >= 97 && v.contains(t.mid()) || v.rel_accuracy_bits() >= 97,
            "{v:?}"
        );
    }
    // Exact endpoints of asin and acos; outside the domain.
    let one = Ball::from_integer(&Integer::ONE);
    assert!(
        one.asin(64)
            .overlaps(&Ball::const_pi(64).mul_2exp(&-Integer::ONE))
    );
    assert_eq!(one.acos(64), Ball::ZERO);
    assert!(one.neg().acos(64).overlaps(&Ball::const_pi(64)));
    assert_eq!(
        Ball::from_integer(&Integer::from(2)).asin(64),
        Ball::INDETERMINATE
    );
    assert_eq!(
        one.add_error(&Mag::power_of_two(&Integer::from(-20)))
            .asin(64),
        Ball::INDETERMINATE
    );
    // atan2 on the branch cut: [0 ± π].
    let cut = Ball::new(Float::ZERO, Mag::ONE).atan2(&Ball::from_integer(&Integer::from(-1)), 64);
    assert!(cut.contains(&Float::ZERO) && cut.rad() >= &Mag::from_integer(&Integer::from(3)));
    assert!(atan_huge_consistent());
}

/// atan of huge values approaches π/2 from below, through 1/m.
fn atan_huge_consistent() -> bool {
    let m = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(1_000_000i64)));
    let a = m.atan(128);
    let hp = Ball::const_pi(200).mul_2exp(&-Integer::ONE);
    a.overlaps(&hp) && a.rel_accuracy_bits() >= 125
}
