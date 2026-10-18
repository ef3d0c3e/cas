//! Complex balls against MPC (`rug::Complex`) brackets at `prec + 192` bits,
//! rounded down and up in each part: both parts must meet their brackets,
//! with a radius within a few ulps of the result's magnitude; branch cuts
//! from both sides; identities; Gaussian rationals.

use proptest::prelude::*;
use rug::float::Round as MpRound;
use rug::{Complex as MpComplex, Float as MpFloat};

use super::{ComplexBall, MUL3_PREC};
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{Bracket, Dyadic, check_encloses, dyadic_with};
use crate::{Complex, Integer, Rational};

fn ball(d: &Dyadic) -> Ball {
    Ball::exact(d.to_ours())
}

fn mp(re: &Dyadic, im: &Dyadic, p: u32) -> MpComplex {
    MpComplex::with_val(p, (re.to_rug().unwrap(), im.to_rug().unwrap()))
}

/// Both parts' brackets of an MPC computation at `prec + 192` bits: MPC's
/// `log` at low precisions misses `log|z|` for `|z|` within about `2^-126`
/// of 1 (at 66 bits, `(1 - 2^-127) + 4.5e-28 i` gives `+1.03e-55` for
/// `-5.88e-39`), so the oracle runs well above where that was seen.
fn brackets(prec: u64, f: impl Fn(u32, (MpRound, MpRound)) -> MpComplex) -> (Bracket, Bracket) {
    let p = u32::try_from(prec).unwrap() + 192;
    let lo = f(p, (MpRound::Down, MpRound::Down));
    let hi = f(p, (MpRound::Up, MpRound::Up));
    let part = |l: &MpFloat, h: &MpFloat| Bracket {
        lo: l.to_rational().expect("finite"),
        hi: h.to_rational().expect("finite"),
    };
    (part(lo.real(), hi.real()), part(lo.imag(), hi.imag()))
}

fn check(
    what: &str,
    z: &ComplexBall,
    prec: u64,
    br: &(Bracket, Bracket),
) -> Result<(), TestCaseError> {
    prop_assert!(z.is_finite(), "{}: indeterminate {:?}", what, z);
    for (part, b) in [(z.re(), &br.0), (z.im(), &br.1)] {
        let mid = part.mid().to_rational().to_rug();
        let rad = part.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, b).map_err(|e| TestCaseError::fail(format!("{what}: {e}")))?;
    }
    // Radii within a few ulps of the larger part.
    let size = Mag::from_float(z.re().mid()).max(Mag::from_float(z.im().mid()));
    if !size.is_zero() {
        let tol = size.mul_2exp(&(Integer::from(6i64) - Integer::from(prec)));
        prop_assert!(
            *z.re().rad() <= tol && *z.im().rad() <= tol,
            "{}: loose {:?}",
            what,
            z
        );
    }
    Ok(())
}

fn prec() -> impl Strategy<Value = u64> {
    prop_oneof![2u64..80, 80u64..300, Just(1000u64)]
}

fn part() -> BoxedStrategy<Dyadic> {
    // Values below 2^4 (mantissas of one limb at 2^-60, two at 2^-124).
    prop_oneof![3 => dyadic_with(2, -130..=-124), 1 => dyadic_with(1, -80..=-60)].boxed()
}

proptest! {
    #[test]
    fn against_mpc(a in part(), b in part(), c in part(), d in part(), prec in prec()) {
        let z = ComplexBall::new(ball(&a), ball(&b));
        let w = ComplexBall::new(ball(&c), ball(&d));
        let (zr, wr) = (mp(&a, &b, 200), mp(&c, &d, 200));
        let zero = zr == 0;
        check("mul", &z.mul(&w, prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, &zr * &wr, r).0))?;
        if wr != 0 {
            check("div", &z.div(&w, prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, &zr / &wr, r).0))?;
        }
        check("exp", &z.exp(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.exp_ref(), r).0))?;
        check("sin", &z.sin(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.sin_ref(), r).0))?;
        check("cos", &z.cos(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.cos_ref(), r).0))?;
        check("sqrt", &z.sqrt(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.sqrt_ref(), r).0))?;
        if !zero {
            check("log", &z.log(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.ln_ref(), r).0))?;
            check("pow", &z.pow(&w, prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, rug::ops::Pow::pow(&zr, &wr), r).0))?;
            // atan has poles at ±i, not hit by dyadic points of this size.
            check("atan", &z.atan(prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, zr.atan_ref(), r).0))?;
        }
        // abs and arg are real balls.
        let real = |what: &str, b: &Ball, f: &dyn Fn(u32, MpRound) -> MpFloat| -> Result<(), TestCaseError> {
            let br = crate::testkit::float::bracket(u32::try_from(prec).unwrap(), f);
            let mid = b.mid().to_rational().to_rug();
            let rad = b.rad().to_float().unwrap().to_rational().to_rug();
            check_encloses(&mid, &rad, &br).map_err(|e| TestCaseError::fail(format!("{what}: {e}")))
        };
        real("abs", &z.abs(prec), &|p, r| MpFloat::with_val_round(p, zr.abs_ref(), r).0)?;
        if !zero {
            real("arg", &z.arg(prec), &|p, r| MpFloat::with_val_round(p, zr.arg_ref(), r).0)?;
        }
    }

    #[test]
    fn identities(a in part(), b in part(), prec in 2u64..300, n in -6i64..8) {
        let z = ComplexBall::new(ball(&a), ball(&b));
        let exact = Complex::new(a.to_ours().to_rational(), b.to_ours().to_rational());
        prop_assert!(z.exp(prec).log(prec).contains_complex(&exact) || b.to_ours().cmp_abs(&Float::from(3i64)) == std::cmp::Ordering::Greater);
        if !(z.is_exact() && z.re().mid().is_zero() && z.im().mid().is_zero()) {
            prop_assert!(z.log(prec).exp(prec).contains_complex(&exact), "exp(log {:?})", z);
            let s = z.sqrt(prec);
            prop_assert!(s.mul(&s, prec).contains_complex(&exact), "sqrt² of {:?}", z);
            prop_assert!(s.re().sign() != Some(std::cmp::Ordering::Less));
        }
        let (s, c) = (z.sin(prec), z.cos(prec));
        let one = s.mul(&s, prec).add(&c.mul(&c, prec), prec);
        prop_assert!(one.contains_complex(&Complex::from(Rational::from(1i64))));
        // Integer powers by squaring against repeated products.
        let mut p = ComplexBall::from_real(Ball::from_integer(&Integer::ONE));
        for _ in 0..n.unsigned_abs() { p = p.mul(&z, prec + 64); }
        if n < 0 { p = p.inv(prec + 64); }
        prop_assert!(z.pow_int(n, prec).overlaps(&p) || !p.is_finite());
        // Cross-precision overlap.
        for (lo, hi) in [(z.exp(prec), z.exp(prec + 400)), (z.sqrt(prec), z.sqrt(prec + 400)), (z.atan(prec), z.atan(prec + 400))] {
            prop_assert!(lo.overlaps(&hi), "{:?} vs {:?}", lo, hi);
        }
    }

    #[test]
    fn three_product_mul(a in dyadic_with(80, -5000..=-5000), b in dyadic_with(80, -5000..=-5000), c in dyadic_with(80, -5000..=-5000), d in dyadic_with(80, -5000..=-5000)) {
        let prec = MUL3_PREC + 100;
        let (z, w) = (ComplexBall::new(ball(&a), ball(&b)), ComplexBall::new(ball(&c), ball(&d)));
        let (zr, wr) = (mp(&a, &b, 6000), mp(&c, &d, 6000));
        check("mul3", &z.mul(&w, prec), prec, &brackets(prec, |p, r| MpComplex::with_val_round(p, &zr * &wr, r).0))?;
    }
}

#[test]
fn branch_cuts() {
    let pi = Ball::const_pi(128);
    let at = |im: Ball| ComplexBall::new(Ball::from_integer(&Integer::from(-4)), im);
    let tiny = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(-60)));
    // Above the cut: arg ≈ +π, sqrt ≈ +2i; below: -π, -2i; on it: the
    // principal +π and +2i.
    for (im, sign) in [(tiny.clone(), 1i64), (tiny.neg(), -1), (Ball::ZERO, 1)] {
        let z = at(im);
        let l = z.log(100);
        // arg = ±(π - atan(2^-62)), within 2^-50 of ±π.
        let off = l
            .im()
            .sub(&pi.mul(&Ball::from_integer(&Integer::from(sign)), 128), 128);
        assert!(
            off.abs().mid() < &Float::from(1i64).mul_2exp(&Integer::from(-50)),
            "{l:?}"
        );
        assert!(l.im().rel_accuracy_bits() >= 95);
        let s = z.sqrt(100);
        assert!(
            s.im()
                .overlaps(&Ball::from_integer(&Integer::from(2 * sign))),
            "{s:?}"
        );
    }
    // Straddling it: both sides.
    let z = at(Ball::new(
        Float::ZERO,
        Mag::power_of_two(&Integer::from(-60)),
    ));
    let l = z.log(100);
    assert!(
        l.im().contains(&Float::from_f64(3.1).unwrap())
            && l.im().contains(&Float::from_f64(-3.1).unwrap())
    );
    let s = z.sqrt(100);
    assert!(s.im().contains(&Float::from(2i64)) && s.im().contains(&Float::from(-2i64)));
    // log of 0, or of a ball around it: indeterminate.
    assert_eq!(ComplexBall::ZERO.log(64), ComplexBall::INDETERMINATE);
    let around0 = ComplexBall::new(Ball::new(Float::ZERO, Mag::ONE), Ball::ZERO);
    assert_eq!(around0.log(64), ComplexBall::INDETERMINATE);
    // atan's cuts on the imaginary axis beyond ±i: right and left of 2i.
    let right = ComplexBall::new(tiny.clone(), Ball::from_integer(&Integer::from(2)));
    let left = ComplexBall::new(tiny.neg(), Ball::from_integer(&Integer::from(2)));
    let hp = Ball::const_pi(128).mul_2exp(&-Integer::ONE);
    // Real parts within 2^-50 of +π/2 and -π/2 (they differ by O(2^-60)).
    let near = |v: &Ball, w: &Ball| {
        v.sub(w, 128).abs().mid() < &Float::from(1i64).mul_2exp(&Integer::from(-50))
    };
    assert!(near(right.atan(100).re(), &hp), "{:?}", right.atan(100));
    assert!(near(left.atan(100).re(), &hp.neg()), "{:?}", left.atan(100));
    // Exactly on the cuts: the side continuous with Re z > 0 (MPC's +0).
    for y in [2i64, -2] {
        let on = ComplexBall::new(Ball::ZERO, Ball::from_integer(&Integer::from(y)));
        let v = on.atan(100);
        assert!(near(v.re(), &hp), "atan({y}i) = {v:?}");
    }
    // A ball straddling the cut below -i: both sides.
    let across = ComplexBall::new(
        Ball::new(Float::ZERO, Mag::power_of_two(&Integer::from(-60))),
        Ball::from_integer(&Integer::from(-2)),
    );
    let v = across.atan(100);
    assert!(
        v.re().contains(&Float::from_f64(1.5).unwrap())
            && v.re().contains(&Float::from_f64(-1.5).unwrap()),
        "{v:?}"
    );
}

#[test]
fn gaussian_rationals_and_reals() {
    let q = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
    let z = Complex::new(q(1, 3), q(-5, 7));
    let b = ComplexBall::from_complex(&z, 100);
    assert!(b.contains_complex(&z) && !b.is_exact());
    let dy = Complex::new(q(3, 4), q(-1, 8));
    assert!(ComplexBall::from_complex(&dy, 100).is_exact());
    // Real inputs stay real.
    let x = ComplexBall::from_real(Ball::from_integer(&Integer::from(2)));
    for f in [
        ComplexBall::exp,
        ComplexBall::sin,
        ComplexBall::cos,
        ComplexBall::sqrt,
        ComplexBall::atan,
        ComplexBall::log,
    ] {
        assert!(f(&x, 64).is_real(), "{:?}", f(&x, 64));
    }
    assert_eq!(
        ComplexBall::from_real(Ball::from_integer(&Integer::from(-9)))
            .sqrt(64)
            .im(),
        &Ball::from_integer(&Integer::from(3))
    );
    // Exact Gaussian integer arithmetic stays exact.
    let i3 = ComplexBall::new(
        Ball::from_integer(&Integer::from(2)),
        Ball::from_integer(&Integer::from(3)),
    );
    let sq = i3.mul(&i3, 64);
    assert_eq!(
        sq,
        ComplexBall::new(
            Ball::from_integer(&Integer::from(-5)),
            Ball::from_integer(&Integer::from(12))
        )
    );
    assert_eq!(sq.div(&i3, 64), i3);
    assert_eq!(i3.pow_int(2, 64), sq);
    assert_eq!(format!("{i3:.3}"), "2 + 3*I");
}
