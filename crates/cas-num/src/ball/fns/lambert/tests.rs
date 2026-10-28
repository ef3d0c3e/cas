//! Tests of `Ball::lambert_w` (CAS-352): values on both branches against
//! giac's to 50 digits, next to the branch point and far out; `W e^W`
//! holds `x`; domains, and radii.

use proptest::prelude::*;
use rug::Float as MpFloat;

use crate::Integer;
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::dyadic_with;

fn float(x: f64) -> Float {
    Float::from_f64(x).unwrap()
}

/// `n/2^k` exactly.
fn dyadic(n: i64, k: i64) -> Float {
    Float::from(n).mul_2exp(&Integer::from(-k))
}

/// `b` holds the decimal `digits` (giac's 50) and has `prec - 3` bits.
fn check(what: &str, b: &Ball, prec: u64, digits: &str) -> Result<(), String> {
    if !b.is_finite() {
        return Err(format!("{what} at {prec}: indeterminate"));
    }
    let p = 400;
    let want = MpFloat::with_val(p, MpFloat::parse(digits).unwrap());
    let mid = MpFloat::with_val(p, b.mid().to_rational().to_rug());
    let rad = MpFloat::with_val(p, b.rad().to_float().unwrap().to_rational().to_rug());
    let error = MpFloat::with_val(p, &mid - &want).abs();
    let slack = MpFloat::with_val(p, &want * 1e-48).abs();
    if error > rad + slack {
        return Err(format!("{what} at {prec}: {b:?} misses {digits}"));
    }
    #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
    let bits = prec as i64 - 3;
    if b.rel_accuracy_bits() < bits {
        return Err(format!(
            "{what} at {prec}: {} bits: {b:?}",
            b.rel_accuracy_bits()
        ));
    }
    Ok(())
}

/// Both branches at moderate arguments, next to the branch point, far out
/// and next to 0.
#[test]
fn values() {
    let cases: [(Float, bool, &str); 14] = [
        (
            float(1.0),
            false,
            "0.567143290409783872999968662210355549753815787186512508135131",
        ),
        (
            float(10.0),
            false,
            "1.7455280027406993830743012648753899115352881290809",
        ),
        (
            dyadic(-1, 2),
            false,
            "-0.35740295618138890306881110405590475331659055507601",
        ),
        (
            dyadic(-1, 2),
            true,
            "-2.1532923641103496491690991500929813755362064853195",
        ),
        (
            dyadic(-47, 7),
            false,
            "-0.93988639805243454196954641321072348096296895907156",
        ),
        (
            dyadic(-47, 7),
            true,
            "-1.0626237132134927642013024316496868661248411292258",
        ),
        (
            dyadic(-6027, 14),
            false,
            "-0.98946609035690924185797336545490797353125329476037",
        ),
        (
            dyadic(-6027, 14),
            true,
            "-1.0106084086921758557607824104536393683453013815822",
        ),
        (
            dyadic(1, -100),
            false,
            "65.138206785153646139535697932661445210192114632140",
        ),
        (
            dyadic(1, 10),
            false,
            "0.97561022024675304998187487109890989177176540629522e-3",
        ),
        (
            dyadic(-1, 70),
            true,
            "-52.480749056918375519920405266209821025763326351797",
        ),
        (
            dyadic(-1, 70),
            false,
            "-0.84703294725430033906903996549337626839170546268365e-21",
        ),
        (
            dyadic(1, 100),
            false,
            "0.78886090522101180541172856528216392814542032093831e-30",
        ),
        (
            float(2.0_f64.powi(60)),
            false,
            "37.952495572432745910977525672145334368657632014413",
        ),
    ];
    let mut wrong = Vec::new();
    for (x, lower, digits) in &cases {
        for prec in [24, 53, 150] {
            let b = Ball::exact(x.clone()).lambert_w(*lower, prec);
            let branch = if *lower { -1 } else { 0 };
            if let Err(e) = check(&format!("W_{branch}({x:?})"), &b, prec, digits) {
                wrong.push(e);
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // W(e) = 1 and W(−1/e) = −1 are not exact inputs; W(0) = 0 is.
    assert_eq!(Ball::exact(Float::ZERO).lambert_w(false, 64), Ball::ZERO);
}

/// Outside each branch's domain, or not shown inside it.
#[test]
fn domains() {
    for x in [-0.37, -0.367_88, -1.0, -1e9] {
        for lower in [false, true] {
            assert_eq!(
                Ball::exact(float(x)).lambert_w(lower, 64),
                Ball::INDETERMINATE,
                "W({x}), lower {lower}"
            );
        }
    }
    for x in [0.0, 1.0, 1e9] {
        assert_eq!(
            Ball::exact(float(x)).lambert_w(true, 64),
            Ball::INDETERMINATE
        );
    }
    // A ball reaching below −1/e, or above 0 on branch −1.
    let around = Ball::new(float(-0.367_8), Mag::power_of_two(&Integer::from(-10)));
    assert_eq!(around.lambert_w(false, 64), Ball::INDETERMINATE);
    let straddle = Ball::new(float(-0.001), Mag::power_of_two(&Integer::from(-8)));
    assert_eq!(straddle.lambert_w(true, 64), Ball::INDETERMINATE);
    assert!(straddle.lambert_w(false, 64).is_finite());
    // Far out: an exponent past f64's.
    let huge = Float::from(1i64).mul_2exp(&Integer::from(1i64 << 20));
    let w = Ball::exact(huge).lambert_w(false, 100);
    assert!(w.is_finite() && w.rel_accuracy_bits() >= 97, "{w:?}");
    let tiny = Float::from(-1i64).mul_2exp(&Integer::from(-2000i64));
    let w = Ball::exact(tiny.clone()).lambert_w(true, 100);
    assert!(w.is_finite() && w.rel_accuracy_bits() >= 97, "{w:?}");
    let w = Ball::exact(tiny.clone()).lambert_w(false, 100);
    assert!(w.contains(&tiny) || w.rel_accuracy_bits() >= 97, "{w:?}");
}

/// A ball's value holds the values at its ends.
#[test]
fn radii() {
    for (m, e, lower) in [
        (1.0, -10i64, false),
        (-0.25, -12, false),
        (-0.25, -12, true),
        (-0.3, -6, true),
        (1e6, -20, false),
    ] {
        let r = Mag::power_of_two(&Integer::from(e));
        let rf = r.to_float().unwrap();
        let b = Ball::new(float(m), r).lambert_w(lower, 100);
        assert!(b.is_finite(), "W({m} ± 2^{e})");
        for end in [float(m).sub(&rf), float(m).add(&rf)] {
            let at = Ball::exact(end).lambert_w(lower, 200);
            assert!(
                b.overlaps(&at) && b.contains(at.mid()),
                "W({m} ± 2^{e}): {b:?} against {at:?}"
            );
        }
    }
}

proptest! {
    /// `W(x) e^W(x)` holds `x`, with `prec - 3` bits, on both branches.
    #[test]
    fn inverse(
        x in dyadic_with(2, 0..=0),
        top in -30i64..=40,
        lower in any::<bool>(),
        prec in prop_oneof![2u64..120, 120u64..500],
    ) {
        let x = x.to_ours();
        if x.is_zero() {
            return Ok(());
        }
        // Branch 0 above −1/e, branch −1 between −1/e and 0.
        let x = if lower {
            x.abs().neg().mul_2exp(&(&Integer::from(top.min(-2)) - x.exponent()))
        } else {
            x.mul_2exp(&(&Integer::from(if x.is_negative() { top.min(-2) } else { top }) - x.exponent()))
        };
        let w = Ball::exact(x.clone()).lambert_w(lower, prec);
        prop_assume!(w.is_finite());
        let q = prec + 64;
        let back = w.mul(&w.exp(q), q);
        prop_assert!(back.contains(&x), "W({:?}) = {:?}: {:?}", x, w, back);
        #[allow(clippy::cast_possible_wrap, reason = "small precisions")]
        let bits = prec as i64 - 3;
        prop_assert!(w.rel_accuracy_bits() >= bits, "W({:?}) at {}: {:?}", x, prec, w);
        let wrong_side = if lower {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Less
        };
        let side = w.cmp(&Ball::from_integer(&Integer::from(-1)));
        prop_assert!(side != Some(wrong_side), "W({:?}) = {:?} off its branch", x, w);
    }
}

/// Arguments of exponents near and beyond `i64`'s (CAS-353): past `2^61`
/// the estimates saturated and Newton's steps ran away, past `Float`'s exact
/// bits (the `special_fns` fuzz target): a huge `x`, and a tiny one on branch
/// −1, are indeterminate; a tiny one on branch 0 is `x` within `2x²`.
#[test]
fn huge_exponents() {
    for far in [
        Integer::from(1i64 << 62) * Integer::from(4),
        Integer::from(1i64 << 62),
    ] {
        let huge = Float::from(1i64).mul_2exp(&far);
        assert_eq!(Ball::exact(huge).lambert_w(false, 64), Ball::INDETERMINATE);
        let tiny = Float::from(-1i64).mul_2exp(&-far.clone());
        assert_eq!(
            Ball::exact(tiny.clone()).lambert_w(true, 64),
            Ball::INDETERMINATE
        );
        let principal = Ball::exact(tiny.clone()).lambert_w(false, 64);
        assert!(
            principal.is_finite() && principal.contains(&tiny),
            "{principal:?}"
        );
        assert!(principal.rel_accuracy_bits() >= 61, "{principal:?}");
    }
    let near = Integer::from(1i64 << 40);
    let big = Ball::exact(Float::from(1i64).mul_2exp(&near)).lambert_w(false, 64);
    assert!(big.is_finite() && big.rel_accuracy_bits() >= 61, "{big:?}");
    let tiny = Float::from(-1i64).mul_2exp(&-near);
    let lower = Ball::exact(tiny.clone()).lambert_w(true, 64);
    assert!(lower.is_finite() && lower.mid().is_negative(), "{lower:?}");
}
