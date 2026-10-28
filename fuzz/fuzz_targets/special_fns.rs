//! The special functions MPFR lacks (CAS-DOC-31 §8, CAS-353):
//!
//! - `si` and `ci` against their series summed in MPFR far above the
//!   precision, at the points of balls within 256 of 0, and at two
//!   precisions;
//! - `polygamma` by its recurrence `ψ^(m)(x + 1) − ψ^(m)(x) = (−1)^m
//!   m!/x^(m+1)` and at two precisions;
//! - `lambert_w` by `W e^W`, which must hold every point of the input, its
//!   branch and its domain.
//!
//! An exact input inside the domain gives a finite result with at least
//! `prec - 3` bits, but for `ci` past `sin`'s range, where `sin x` is
//! `[0 ± 1]`.
#![no_main]

use std::cmp::Ordering;

use arbitrary::{Result, Unstructured};
use cas_num::ball::{Ball, Float, Mag, Round};
use cas_num::{Integer, Natural};
use cas_num_fuzz::ball::{Range, float, from_mp, meets, points, precision, radius, to_mp, warm_up};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;
use rug::float::Constant;

/// The most bits fuzzed, as `ball_fns`'s for `zeta`: the series in MPFR
/// runs at about twice that.
const MAX_PRECISION: u64 = 2048;

/// Arguments of exponent above this have `sin` and `cos` in `[−1, 1]`.
const TRIG_HUGE_EXP: i64 = 1 << 16;

fuzz_target!(|data: &[u8]| {
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let op = u.int_in_range(0u8..=3)?;
    let prec = precision(u)?.min(MAX_PRECISION);
    let m = float(u, Range::Wide)?;
    let r = radius(u, &m, prec, Range::Wide)?;
    let x = Ball::new(m, r);
    match op {
        0 | 1 => sin_cos(u, op == 1, &x, prec),
        2 => polygamma(u, &x, prec),
        _ => lambert(u, &x, prec),
    }
}

fn exponent(x: &Float) -> i64 {
    let e = x.exponent();
    e.to_i64()
        .unwrap_or(if e.is_negative() { i64::MIN } else { i64::MAX })
}

/// An exact input's result has `prec - 3` bits, or is 0.
fn accurate(y: &Ball, prec: u64, what: &dyn Fn() -> String) {
    assert!(y.is_finite(), "{}: indeterminate", what());
    let want = i64::try_from(prec).expect("small") - 3;
    assert!(
        y.mid().is_zero() || y.rel_accuracy_bits() >= want,
        "{}: {} bits: {y:?}",
        what(),
        y.rel_accuracy_bits()
    );
}

/// `Si(x)` or `Ci(x)` by the series in MPFR at enough bits for an error far
/// below `2^-(prec + 64)`, bracketed by that much; `None` beyond 256.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    reason = "small sizes"
)]
fn series(x: &Float, prec: u64, cosine: bool) -> Option<(Float, Float)> {
    let xr = to_mp(x)?;
    let size = xr.to_f64().abs();
    if size > 256.0 || xr.is_zero() || (cosine && size <= 0.0) {
        return None;
    }
    let p = (prec as f64 + 2.0 * size * std::f64::consts::LOG2_E + 256.0) as u32;
    let xr = MpFloat::with_val(p, xr);
    let square = MpFloat::with_val(p, &xr * &xr);
    let (mut power, mut sum, mut j) = if cosine {
        (MpFloat::with_val(p, 1), MpFloat::new(p), 0u32)
    } else {
        (xr.clone(), xr.clone(), 1u32)
    };
    let tiny = MpFloat::with_val(p, MpFloat::i_exp(1, -(p as i32)));
    loop {
        j += 2;
        power *= &square;
        power /= (j - 1) * j;
        power = -power;
        let term = MpFloat::with_val(p, &power / j);
        sum += &term;
        if f64::from(j) > 2.0 * size && term.abs() < tiny {
            break;
        }
    }
    if cosine {
        sum += MpFloat::with_val(p, Constant::Euler);
        sum += MpFloat::with_val(p, xr.ln_ref());
    }
    let slack = MpFloat::with_val(p, MpFloat::i_exp(1, -(prec as i32) - 64));
    let lo = MpFloat::with_val(p, &sum - &slack);
    let hi = MpFloat::with_val(p, &sum + &slack);
    Some((from_mp(&lo)?, from_mp(&hi)?))
}

fn sin_cos(u: &mut Unstructured<'_>, cosine: bool, x: &Ball, prec: u64) -> Result<()> {
    let name = if cosine { "ci" } else { "si" };
    let f = |b: &Ball, p: u64| if cosine { b.ci(p) } else { b.si(p) };
    let y = f(x, prec);
    let what = || format!("{name}({x:?}) at {prec} bits");
    let low = Mag::from_float_lower(&x.mid().abs()).sub_lower(x.rad());
    if cosine && (x.mid().is_negative() || low.is_zero()) {
        assert_eq!(y, Ball::INDETERMINATE, "{}: outside the domain", what());
        return Ok(());
    }
    if !x.rad().is_finite() {
        return Ok(());
    }
    assert!(y.is_finite(), "{}: indeterminate", what());
    let z = f(x, prec + 64);
    assert!(
        y.overlaps(&z),
        "{}: {y:?} against {z:?} at {} bits",
        what(),
        prec + 64
    );
    for p in points(u, x, prec)? {
        if let Some((lo, hi)) = series(&p, prec, cosine) {
            assert!(
                meets(&y, &lo, &hi),
                "{}: {y:?} misses [{lo:?}, {hi:?}] at {p:?}",
                what()
            );
        }
    }
    if x.is_exact() && !(cosine && exponent(x.mid()) > TRIG_HUGE_EXP) {
        accurate(&y, prec, &what);
    }
    Ok(())
}

/// Whether the ball meets 0 or a negative integer, from its ends rounded
/// outward.
fn meets_pole(x: &Ball) -> bool {
    let Some(r) = x.rad().to_float() else {
        return true;
    };
    let bits = x.mid().bit_length().max(62) + 2;
    let lo = x.mid().sub_round(&r, bits, Round::Floor).0;
    let hi = x.mid().add_round(&r, bits, Round::Ceil).0;
    if lo.sign() == Ordering::Greater {
        return false;
    }
    // lo ≤ 0: 0 is in the ball, or an integer lies in [lo, hi] below it.
    if hi.sign() != Ordering::Less {
        return true;
    }
    ceil(&lo) <= floor(&hi)
}

/// `⌊x⌋` for `x < 0`, without a rational below 1 in magnitude.
fn floor(x: &Float) -> Integer {
    if exponent(x) <= 0 {
        return Integer::from(-1);
    }
    x.to_rational().floor()
}

/// `⌈x⌉` for `x ≤ 0`, likewise.
fn ceil(x: &Float) -> Integer {
    if x.is_zero() || exponent(x) <= 0 {
        return Integer::ZERO;
    }
    -(-x.to_rational()).floor()
}

fn polygamma(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let m: u64 = match u.int_in_range(0u8..=15)? {
        0 => u.int_in_range(13..=257)?,
        _ => u.int_in_range(0..=12)?,
    };
    let y = x.polygamma(m, prec);
    let what = || format!("ψ^({m})({x:?}) at {prec} bits");
    if m > 256 {
        assert_eq!(y, Ball::INDETERMINATE, "{}: past the order", what());
        return Ok(());
    }
    // Ends past 2^20 have rational floors too long to take.
    let far = |e: Option<Float>| e.is_none_or(|e| exponent(&e) > 1 << 20);
    let r = x.rad().to_float();
    if far(Some(x.mid().clone())) || far(r.clone()) {
        return Ok(());
    }
    if meets_pole(x) {
        assert_eq!(y, Ball::INDETERMINATE, "{}: a pole", what());
        return Ok(());
    }
    if !y.is_finite() {
        // A ball whose radius rounding reaches a pole may be indeterminate;
        // an exact one may not.
        assert!(!x.is_exact(), "{}: indeterminate", what());
        return Ok(());
    }
    let z = x.polygamma(m, prec + 64);
    assert!(y.overlaps(&z), "{}: {y:?} against {z:?}", what());
    if x.is_exact() {
        accurate(&y, prec, &what);
        // The recurrence, where x + 1 is short.
        let reach = 2 * i64::try_from(prec).expect("small") + 64;
        if (-reach..reach).contains(&exponent(x.mid())) {
            let next = Ball::exact(x.mid().add(&Float::from(1i64))).polygamma(m, prec + 64);
            if next.is_finite() {
                let q = prec + 64;
                let i = i64::try_from(m + 1).expect("small");
                let factorial = Ball::from_integer(&Integer::from(Natural::factorial(m)));
                let mut step = factorial.div(&x.pow(i, q), q);
                if m % 2 == 1 {
                    step = step.neg();
                }
                let difference = next.sub(&z, q);
                assert!(
                    difference.overlaps(&step),
                    "{}: ψ(x + 1) − ψ(x) = {difference:?}, not {step:?}",
                    what()
                );
            }
        }
    }
    Ok(())
}

fn lambert(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let lower: bool = u.arbitrary()?;
    // On branch −1 the ball is moved into [−1/e, 0) most of the time.
    let x = if lower && u.ratio(3u8, 4)? && !x.mid().is_zero() {
        let k = u.int_in_range(2i64..=200)?;
        let m = x.mid().abs().neg();
        let scale = Integer::from(-k) - m.exponent();
        Ball::new(m.mul_2exp(&scale), x.rad().mul_2exp(&scale))
    } else {
        x.clone()
    };
    let w = x.lambert_w(lower, prec);
    let what = || format!("W_{}({x:?}) at {prec} bits", if lower { -1 } else { 0 });
    let Some(r) = x.rad().to_float() else {
        return Ok(());
    };
    let bits = x.mid().bit_length().max(62) + 2;
    let hi = x.mid().add_round(&r, bits, Round::Ceil).0;
    // Below −6028/2^14 < −1/e, or on branch −1 from 0 up: outside the
    // domain.
    let below = Float::from(-6028i64).mul_2exp(&Integer::from(-14i64));
    if hi < below || (lower && hi.sign() != Ordering::Less) {
        assert_eq!(w, Ball::INDETERMINATE, "{}: outside the domain", what());
        return Ok(());
    }
    if !w.is_finite() {
        return Ok(());
    }
    let q = prec + 64;
    let back = w.mul(&w.exp(q), q);
    for p in points(u, &x, prec)? {
        assert!(
            back.contains(&p),
            "{}: W e^W = {back:?} misses {p:?}",
            what()
        );
    }
    let side = w.cmp(&Ball::from_integer(&Integer::from(-1)));
    let wrong = if lower {
        Ordering::Greater
    } else {
        Ordering::Less
    };
    assert!(side != Some(wrong), "{}: {w:?} off its branch", what());
    if x.is_exact() {
        accurate(&w, prec, &what);
    }
    Ok(())
}
