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
