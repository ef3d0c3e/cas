//! The cores behind each function (CAS-86): the baseline, the tables and
//! bit-burst, called directly on their reduced domains at any precision
//! (bit-burst is cheap at low precision, so its chunking and series get
//! many cases there), against MPFR brackets at `wp + 64` bits and against
//! each other (the baseline is the reference). The public functions at the
//! thresholds and above `BIT_BURST_MIN_PREC` close the loop.

use proptest::prelude::*;
use rug::Float as MpFloat;

use super::atan::atan_halving;
use super::exp::exp_small;
use super::trig::sin_cos_small;
use super::{BIT_BURST_MIN_PREC, TABLE_MAX_PREC, bitburst, tables};
use crate::Integer;
use crate::ball::{Ball, Float};
use crate::testkit::ToRug;
use crate::testkit::float::{Dyadic, bracket, check_encloses, dyadic_with};

type Oracle = fn(&MpFloat, u32, rug::float::Round) -> MpFloat;

/// `v` encloses `oracle(x)` and knows `wp - 4` bits of it.
fn check(name: &str, v: &Ball, x: &MpFloat, oracle: Oracle, wp: u64) -> Result<(), TestCaseError> {
    prop_assert!(v.is_finite(), "{} {}: indeterminate", name, x);
    let br = bracket(u32::try_from(wp).unwrap(), |p, r| oracle(x, p, r));
    let mid = v.mid().to_rational().to_rug();
    let rad = v.rad().to_float().unwrap().to_rational().to_rug();
    check_encloses(&mid, &rad, &br)
        .map_err(|e| TestCaseError::fail(format!("{name}({x}) at {wp}: {e}")))?;
    if !v.mid().is_zero() {
        #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
        let want = wp as i64 - 4;
        prop_assert!(
            v.rel_accuracy_bits() >= want,
            "{}({}) at {}: only {} bits",
            name,
            x,
            wp,
            v.rel_accuracy_bits()
        );
    }
    Ok(())
}

fn exp_oracle(x: &MpFloat, p: u32, r: rug::float::Round) -> MpFloat {
    MpFloat::with_val_round(p, x.exp_ref(), r).0
}

fn sin_oracle(x: &MpFloat, p: u32, r: rug::float::Round) -> MpFloat {
    MpFloat::with_val_round(p, x.sin_ref(), r).0
}

fn cos_oracle(x: &MpFloat, p: u32, r: rug::float::Round) -> MpFloat {
    MpFloat::with_val_round(p, x.cos_ref(), r).0
}

fn log1p_oracle(x: &MpFloat, p: u32, r: rug::float::Round) -> MpFloat {
    MpFloat::with_val_round(p, x.ln_1p_ref(), r).0
}

fn atan_oracle(x: &MpFloat, p: u32, r: rug::float::Round) -> MpFloat {
    MpFloat::with_val_round(p, x.atan_ref(), r).0
}

/// Every core at the exact `t` (each on the part of its domain `t` is in),
/// checked against MPFR and the baseline.
fn check_cores(t: &Float, wp: u64, tables_too: bool) -> Result<(), TestCaseError> {
    let x = t.to_rug().to_rug().expect("in MPFR's range");
    let overlap = |name: &str, a: &Ball, b: &Ball| {
        prop_assert!(
            a.overlaps(b),
            "{} at {}: {:?} vs baseline {:?}",
            name,
            wp,
            a,
            b
        );
        Ok(())
    };
    let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
    let quarter_pi = 0.785f64;
    // exp on |t| < 1/2.
    if t.cmp_abs(&half) == std::cmp::Ordering::Less {
        let base = exp_small(t, wp);
        check("exp baseline", &base, &x, exp_oracle, wp)?;
        let bb = bitburst::exp(t, wp);
        check("exp bit-burst", &bb, &x, exp_oracle, wp)?;
        overlap("exp bit-burst", &bb, &base)?;
        if tables_too {
            let tb = tables::exp(t, wp);
            check("exp table", &tb, &x, exp_oracle, wp)?;
            overlap("exp table", &tb, &base)?;
        }
    }
    // sin, cos on |t| < π/4.
    if t.to_f64().abs() < quarter_pi {
        let (bs, bc) = sin_cos_small(t, wp);
        check("sin baseline", &bs, &x, sin_oracle, wp)?;
        check("cos baseline", &bc, &x, cos_oracle, wp)?;
        let (s, c) = bitburst::sin_cos(t, wp);
        check("sin bit-burst", &s, &x, sin_oracle, wp)?;
        check("cos bit-burst", &c, &x, cos_oracle, wp)?;
        overlap("sin bit-burst", &s, &bs)?;
        overlap("cos bit-burst", &c, &bc)?;
        if tables_too {
            let (s, c) = tables::sin_cos(t, wp).expect("within the table");
            check("sin table", &s, &x, sin_oracle, wp)?;
            check("cos table", &c, &x, cos_oracle, wp)?;
            overlap("sin table", &s, &bs)?;
            overlap("cos table", &c, &bc)?;
        }
    }
    let a = t.abs();
    let xa = x.clone().abs();
    if a.cmp_abs(&Float::from(1i64)) != std::cmp::Ordering::Greater {
        // atan on [0, 1].
        let base = atan_halving(&Ball::exact(a.clone()), wp);
        check("atan baseline", &base, &xa, atan_oracle, wp)?;
        let bb = bitburst::atan(&a, wp);
        check("atan bit-burst", &bb, &xa, atan_oracle, wp)?;
        overlap("atan bit-burst", &bb, &base)?;
        if tables_too {
            let tb = tables::atan(&a, wp);
            check("atan table", &tb, &xa, atan_oracle, wp)?;
            overlap("atan table", &tb, &base)?;
        }
    }
    // log1p on [0, 1): only the tables have a core of their own.
    if tables_too && a.cmp_abs(&Float::from(1i64)) == std::cmp::Ordering::Less {
        check("log1p table", &tables::log1p(&a, wp), &xa, log1p_oracle, wp)?;
    }
    Ok(())
}

/// Reduced arguments: below 1 in magnitude, down to tiny, with long
/// mantissas.
fn reduced() -> BoxedStrategy<Dyadic> {
    prop_oneof![
        4 => dyadic_with(2, -20..=0),
        1 => dyadic_with(4, -200..=-20),
    ]
    .boxed()
}

proptest! {
    #[test]
    fn cores_against_mpfr(t in reduced(), wp in prop_oneof![2u64..200, 200u64..1200]) {
        check_cores(&t.to_ours(), wp, true)?;
    }
}

/// Points of each kind: table steps exactly, just around them, tiny, the
/// domain edges.
fn points() -> Vec<Float> {
    let f = |x: f64| Float::from_f64(x).unwrap();
    let step = |i: i64, k: i64| Float::from(i).mul_2exp(&Integer::from(-k));
    let mut v = vec![
        f(0.0),
        f(0.1),
        f(-0.3),
        f(0.499_999),
        f(-0.499_999),
        f(0.785),
        f(-0.785),
        f(0.75),
        f(0.999_999_999),
        f(1.0),
        f(1e-30),
        f(-3e-6),
        step(1, 8),
        step(201, 8),
        step(255, 8),
        step(1, 16),
        step(255, 16),
        step(65_535, 16),
        step(-127, 8),
    ];
    for x in v.clone() {
        if !x.is_zero() {
            let ulp = Float::from(1i64).mul_2exp(&(x.exponent() - Integer::from(70)));
            v.push(x.add(&ulp));
            v.push(x.sub(&ulp));
        }
    }
    v
}

#[test]
fn cores_at_the_edges() {
    for wp in [2u64, 64, 333, 2000] {
        for t in points() {
            check_cores(&t, wp, true).unwrap();
        }
    }
}

#[test]
fn cores_at_the_thresholds() {
    for t in [0.3f64, -0.4, 0.78, 0.9, 1e-12] {
        let t = Float::from_f64(t).unwrap();
        check_cores(&t, TABLE_MAX_PREC, true).unwrap();
        check_cores(&t, BIT_BURST_MIN_PREC, false).unwrap();
    }
}

/// The public functions at precisions whose working precision straddles
/// each threshold, and above bit-burst's.
#[test]
fn public_functions_across_the_thresholds() {
    type Fun = (&'static str, fn(&Ball, u64) -> Ball, Oracle);
    let funs: [Fun; 7] = [
        ("exp", Ball::exp, exp_oracle),
        ("expm1", Ball::expm1, |x, p, r| {
            MpFloat::with_val_round(p, x.exp_m1_ref(), r).0
        }),
        ("log", Ball::log, |x, p, r| {
            MpFloat::with_val_round(p, x.ln_ref(), r).0
        }),
        ("log1p", Ball::log1p, log1p_oracle),
        ("sin", Ball::sin, sin_oracle),
        ("cos", Ball::cos, cos_oracle),
        ("atan", Ball::atan, atan_oracle),
    ];
    let tw = TABLE_MAX_PREC - 33;
    let bw = BIT_BURST_MIN_PREC - 34;
    for prec in [tw - 1, tw + 1, bw - 1, bw + 1, 20_000] {
        for x in [2.5f64, -0.375, 0.0078, 1.0625, 13.0, 1e-9] {
            let xr = MpFloat::with_val(53, x);
            let b = Ball::exact(Float::from_f64(x).unwrap());
            for (name, f, oracle) in funs {
                if (name == "log" && x <= 0.0) || (name == "log1p" && x <= -1.0) {
                    continue;
                }
                check(name, &f(&b, prec), &xr, oracle, prec).unwrap();
            }
        }
    }
}
