//! The Bessel and Airy functions of balls (CAS-DOC-34 §7, CAS-381), at
//! arguments below `2^20` in size:
//!
//! - `J_n` and `Y_n` against MPFR's `jn` and `yn`, where they stay quick
//!   (`|x| ≤ 2^12` or `n²/2 < |x|`), and `Ai` against its `ai` (for `|x| ≤
//!   64`), at the points of the ball;
//! - `I_n` and `K_n` by their Wronskian `I_n K_(n+1) + I_(n+1) K_n = 1/x`,
//!   and `Ai′`, `Bi` and `Bi′` by Airy's, `Ai Bi′ − Ai′ Bi = 1/π`, over
//!   the whole ball;
//! - the recurrences of `J`, `Y`, `I` and `K` in their order, over the
//!   ball;
//! - two precisions agree; `Y` and `K` are indeterminate unless the ball
//!   lies above 0, and every function past the largest order.
//!
//! An exact input inside the domain gives a finite result with at least
//! `prec - 3` bits.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::Integer;
use cas_num::ball::{Airy, Ball, Mag, Round};
use cas_num_fuzz::ball::{
    Oracle, Range, bracket, float, meets, points, precision, radius, to_mp, warm_up,
};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;

/// The most bits fuzzed.
const MAX_PRECISION: u64 = 2048;

/// The largest order the functions take, as cas-num's.
const MAX_ORDER: i64 = 1024;

fuzz_target!(|data: &[u8]| {
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let op = u.int_in_range(0u8..=5)?;
    let prec = precision(u)?.min(MAX_PRECISION);
    let m = float(u, Range::Narrow)?;
    let r = radius(u, &m, prec, Range::Narrow)?;
    let x = Ball::new(m, r);
    match op {
        0 => against_mpfr(u, &x, prec, false),
        1 => against_mpfr(u, &x, prec, true),
        2 => modified(u, &x, prec),
        3 => airy(u, &x, prec),
        _ => recurrences(u, &x, prec),
    }
}

/// An order: mostly small, sometimes up to the largest, rarely past it.
fn order(u: &mut Unstructured<'_>) -> Result<i64> {
    Ok(match u.int_in_range(0u8..=15)? {
        0 => u.int_in_range(-MAX_ORDER..=MAX_ORDER)?,
        1 => MAX_ORDER + u.int_in_range(1..=4)?,
        _ => u.int_in_range(-12..=12)?,
    })
}

/// Whether the ball lies above 0, as cas-num tells it: its lower end
/// rounded down at 64 bits.
fn positive(x: &Ball) -> bool {
    let Some(r) = x.rad().to_float() else {
        return false;
    };
    let lo = x.mid().sub_round(&r, 64, Round::Floor).0;
    !lo.is_negative() && !lo.is_zero()
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

/// The ball and the same at 64 more bits agree.
fn agrees(y: &Ball, z: &Ball, what: &dyn Fn() -> String) {
    assert!(
        y.overlaps(z),
        "{}: {y:?} against {z:?} at 64 more bits",
        what()
    );
}

fn against_mpfr(u: &mut Unstructured<'_>, x: &Ball, prec: u64, second: bool) -> Result<()> {
    let n = order(u)?;
    let f = |b: &Ball, p: u64| {
        if second {
            b.bessel_y(n, p)
        } else {
            b.bessel_j(n, p)
        }
    };
    let y = f(x, prec);
    let name = if second { "Y" } else { "J" };
    let what = || format!("{name}_{n}({x:?}) at {prec} bits");
    if n.abs() > MAX_ORDER || (second && !positive(x)) {
        assert_eq!(y, Ball::INDETERMINATE, "{}: outside the domain", what());
        return Ok(());
    }
    if !x.rad().is_finite() {
        return Ok(());
    }
    assert!(y.is_finite(), "{}: indeterminate", what());
    agrees(&y, &f(x, prec + 64), &what);
    let order = i32::try_from(n).expect("a small order");
    for p in points(u, x, prec)? {
        let Some(point) = to_mp(&p) else {
            continue;
        };
        // MPFR's expansion needs n²/2 < |x| (its jyn_asympt.c); short of
        // that it takes the series, about |x| terms at 1.44 |x| more bits,
        // too slow past 2^12 for the fuzzer's minute: J_(−1024) at about
        // −4.1·10⁵ timed out (CAS-381).
        let size = point.to_f64().abs();
        let square = (n * n) as f64;
        if size > 4096.0 && square / 2.0 >= size {
            continue;
        }
        let oracle = bracket(prec, |bits, rnd| {
            if second {
                MpFloat::with_val_round(bits, point.yn_ref(order), rnd).0
            } else {
                MpFloat::with_val_round(bits, point.jn_ref(order), rnd).0
            }
        });
        if let Oracle::Bracket(lo, hi) = oracle {
            assert!(
                meets(&y, &lo, &hi),
                "{}: {y:?} misses [{lo:?}, {hi:?}] at {p:?}",
                what()
            );
        }
    }
    if x.is_exact() {
        accurate(&y, prec, &what);
    }
    Ok(())
}

fn modified(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let n = order(u)?.clamp(-MAX_ORDER, MAX_ORDER - 1);
    let what = || format!("I_{n}, K_{n}({x:?}) at {prec} bits");
    let (i0, i1) = (x.bessel_i(n, prec), x.bessel_i(n + 1, prec));
    let (k0, k1) = (x.bessel_k(n, prec), x.bessel_k(n + 1, prec));
    if !positive(x) {
        assert_eq!(k0, Ball::INDETERMINATE, "{}: outside the domain", what());
    }
    if !x.rad().is_finite() {
        return Ok(());
    }
    assert!(i0.is_finite(), "{}: I indeterminate", what());
    agrees(&i0, &x.bessel_i(n, prec + 64), &what);
    if x.is_exact() {
        accurate(&i0, prec, &what);
    }
    if positive(x) {
        assert!(k0.is_finite(), "{}: K indeterminate", what());
        agrees(&k0, &x.bessel_k(n, prec + 64), &what);
        if x.is_exact() {
            accurate(&k0, prec, &what);
        }
        let q = prec + 64;
        let w = i0.mul(&k1, q).add(&i1.mul(&k0, q), q);
        let expected = Ball::from_integer(&Integer::ONE).div(x, q);
        assert!(
            w.overlaps(&expected),
            "{}: the Wronskian is {w:?}, not {expected:?}",
            what()
        );
    }
    Ok(())
}

fn airy(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let which = [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime][u.int_in_range(0..=3usize)?];
    let what = || format!("{which:?}({x:?}) at {prec} bits");
    let y = x.airy(which, prec);
    if !x.rad().is_finite() {
        return Ok(());
    }
    assert!(y.is_finite(), "{}: indeterminate", what());
    agrees(&y, &x.airy(which, prec + 64), &what);
    if x.is_exact() {
        accurate(&y, prec, &what);
    }
    let small = Mag::from_float(x.mid()).add(x.rad()) <= Mag::from_integer(&Integer::from(64));
    if which == Airy::Ai && small {
        for p in points(u, x, prec)? {
            let Some(point) = to_mp(&p) else {
                continue;
            };
            let oracle = bracket(prec, |bits, rnd| {
                MpFloat::with_val_round(bits, point.ai_ref(), rnd).0
            });
            if let Oracle::Bracket(lo, hi) = oracle {
                assert!(
                    meets(&y, &lo, &hi),
                    "{}: {y:?} misses [{lo:?}, {hi:?}] at {p:?}",
                    what()
                );
            }
        }
    }
    let q = prec + 64;
    let [ai, aip, bi, bip] =
        [Airy::Ai, Airy::AiPrime, Airy::Bi, Airy::BiPrime].map(|w| x.airy(w, prec));
    let w = ai.mul(&bip, q).sub(&aip.mul(&bi, q), q);
    let expected = Ball::from_integer(&Integer::ONE).div(&Ball::const_pi(q), q);
    assert!(w.overlaps(&expected), "{}: the Wronskian is {w:?}", what());
    Ok(())
}

fn recurrences(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let n = u.int_in_range(-40i64..=40)?;
    if !x.rad().is_finite() {
        return Ok(());
    }
    let q = prec + 64;
    let scale = Ball::from_integer(&Integer::from(2 * n)).div(x, q);
    let kinds: [(&str, fn(&Ball, i64, u64) -> Ball, i8); 4] = [
        ("J", Ball::bessel_j, 1),
        ("Y", Ball::bessel_y, 1),
        ("I", Ball::bessel_i, -1),
        ("K", Ball::bessel_k, 0),
    ];
    for (name, f, sign) in kinds {
        let (before, now, after) = (f(x, n - 1, prec), f(x, n, prec), f(x, n + 1, prec));
        if !(before.is_finite() && now.is_finite() && after.is_finite() && scale.is_finite()) {
            continue;
        }
        let left = match sign {
            1 => before.add(&after, q),
            -1 => before.sub(&after, q),
            _ => after.sub(&before, q),
        };
        let right = scale.mul(&now, q);
        assert!(
            left.overlaps(&right),
            "{name}: C_{} and C_{} at {x:?}: {left:?} against (2n/x) C_{n} = {right:?}",
            n - 1,
            n + 1
        );
    }
    Ok(())
}
