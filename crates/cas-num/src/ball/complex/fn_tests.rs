//! CAS-102's functions against MPC: exact inputs at random points, on and
//! next to the branch cuts, near 0 and for large imaginary parts; inexact
//! balls must contain the value at their corners (across a cut, both
//! sides); identities, real-axis agreement and the exact values.
//!
//! The brackets are MPC's directed roundings at `prec + 192` bits where the
//! function is one MPC operation; `expm1` is `exp` rounded the same way
//! minus 1, `log1p` the logarithm of the exact `1 + z`. `root` and
//! `sin_pi`/`cos_pi` are MPC's `pow(z, 1/n)` and `sin(π z)` at `prec +
//! 300` bits, rounded inputs, with an explicit margin above their error.

#![allow(
    clippy::many_single_char_names,
    reason = "z, w, p, v, n and the parts a, b as in the oracles"
)]

use proptest::prelude::*;
use rug::float::{Constant, Round as MpRound};
use rug::ops::Pow;
use rug::{Complex as MpComplex, Float as MpFloat};

use super::ComplexBall;
use super::tests::{ball, brackets, check, mp, part, prec};
use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{Bracket, Dyadic, check_encloses};
use crate::{Complex, Integer, Rational};

type Ours = fn(&ComplexBall, u64) -> ComplexBall;
type Oracle = fn(&MpComplex, u64) -> (Bracket, Bracket);

/// A bracket from one correctly rounded MPC operation.
macro_rules! mpc {
    ($method:ident) => {
        |z: &MpComplex, prec: u64| {
            brackets(prec, |p, r| MpComplex::with_val_round(p, z.$method(), r).0)
        }
    };
}

/// `[v - eps, v + eps]` in each part.
fn around(v: &MpComplex, eps: &MpFloat) -> (Bracket, Bracket) {
    let part = |x: &MpFloat| Bracket {
        lo: MpFloat::with_val_round(x.prec(), x - eps, MpRound::Down)
            .0
            .to_rational()
            .expect("finite"),
        hi: MpFloat::with_val_round(x.prec(), x + eps, MpRound::Up)
            .0
            .to_rational()
            .expect("finite"),
    };
    (part(v.real()), part(v.imag()))
}

fn expm1_oracle(z: &MpComplex, prec: u64) -> (Bracket, Bracket) {
    // 140 more bits: |z| >= 2^-130 cancels against 1 by at most that much.
    brackets(prec + 140, |p, r| {
        let e = MpComplex::with_val_round(p, z.exp_ref(), r).0;
        let re = MpFloat::with_val_round(p, e.real() - 1u32, r.0).0;
        MpComplex::with_val(p, (re, e.imag()))
    })
}

fn log1p_oracle(z: &MpComplex, prec: u64) -> (Bracket, Bracket) {
    // The parts are below 2^5 with bits down to 2^-200: 1 + z is exact.
    let w = MpComplex::with_val(512, z + 1u32);
    brackets(prec, |p, r| MpComplex::with_val_round(p, w.ln_ref(), r).0)
}

/// `sin(π z)` or `cos(π z)` at `prec + 300` bits: the rounded `π z` is
/// within `2 |π z| 2^-p` and `|sin'|, |cos'| <= cosh(π b)`, so `(1 + |z|)
/// cosh(π b) 2^-(p - 8)` covers the error.
fn pi_oracle(z: &MpComplex, prec: u64, cos: bool) -> (Bracket, Bracket) {
    let p = u32::try_from(prec).unwrap() + 300;
    let pi = MpFloat::with_val(p, Constant::Pi);
    let arg = MpComplex::with_val(p, z * &pi);
    let v = if cos {
        MpComplex::with_val(p, arg.cos_ref())
    } else {
        MpComplex::with_val(p, arg.sin_ref())
    };
    let size = MpFloat::with_val(64, z.abs_ref()) + 1u32;
    let grow = MpFloat::with_val(64, arg.imag().clone().abs().cosh());
    let eps = (size * grow) >> (p - 8);
    around(&v, &eps)
}

/// MPC's `pow(z, 1/n)` at `prec + 300` bits, `1/n` rounded: the rounding
/// moves the result by `|v| |log z| 2^-p`, and `|log z| < 2^7` for these
/// inputs.
fn root_oracle(z: &MpComplex, prec: u64, n: u32) -> (Bracket, Bracket) {
    let p = u32::try_from(prec).unwrap() + 300;
    let w = MpComplex::with_val(p, MpFloat::with_val(p, 1u32) / n);
    let v = MpComplex::with_val(p, z.pow(&w));
    let size = MpFloat::with_val(64, v.abs_ref()) + 1u32;
    let eps = size >> (p - 16);
    around(&v, &eps)
}

const FUNCTIONS: [(&str, Ours, Oracle); 12] = [
    ("tan", ComplexBall::tan, mpc!(tan_ref)),
    ("tanh", ComplexBall::tanh, mpc!(tanh_ref)),
    ("sinh", ComplexBall::sinh, mpc!(sinh_ref)),
    ("cosh", ComplexBall::cosh, mpc!(cosh_ref)),
    ("asin", ComplexBall::asin, mpc!(asin_ref)),
    ("acos", ComplexBall::acos, mpc!(acos_ref)),
    ("asinh", ComplexBall::asinh, mpc!(asinh_ref)),
    ("acosh", ComplexBall::acosh, mpc!(acosh_ref)),
    ("atanh", ComplexBall::atanh, mpc!(atanh_ref)),
    ("expm1", ComplexBall::expm1, expm1_oracle),
    ("log1p", ComplexBall::log1p, log1p_oracle),
    ("sin_pi", ComplexBall::sin_pi, |z, p| pi_oracle(z, p, false)),
];

/// Whether `f` has a pole or a logarithmic singularity at `z`.
fn singular(name: &str, z: &MpComplex) -> bool {
    let at = |re: i32, im: i32| *z.real() == re && *z.imag() == im;
    match name {
        "atanh" => at(1, 0) || at(-1, 0),
        "log1p" => at(-1, 0),
        _ => false,
    }
}

/// Checks every function at the exact point `re + im i`.
fn check_point(re: &Dyadic, im: &Dyadic, prec: u64) -> Result<(), TestCaseError> {
    let z = ComplexBall::new(ball(re), ball(im));
    let zr = mp(re, im, 2000);
    for (name, ours, oracle) in FUNCTIONS {
        if !singular(name, &zr) {
            check(name, &ours(&z, prec), prec, &oracle(&zr, prec))?;
        }
    }
    check("cos_pi", &z.cos_pi(prec), prec, &pi_oracle(&zr, prec, true))?;
    for n in [3u32, 4, 7] {
        let name = format!("root({n})");
        check(
            &name,
            &z.root(u64::from(n), prec),
            prec,
            &root_oracle(&zr, prec, n),
        )?;
    }
    Ok(())
}

/// Exact 0, ±1 and ±2^-k as well as [`part`]'s values: points on and next
/// to the cuts and the branch points.
fn edge_part() -> BoxedStrategy<Dyadic> {
    let pow2 = |k: i64, neg: bool| {
        Dyadic::new(
            &rug::Integer::from(if neg { -1 } else { 1 }),
            &rug::Integer::from(k),
        )
    };
    prop_oneof![
        4 => part(),
        1 => Just(Dyadic::default()),
        1 => (any::<bool>()).prop_map(move |n| pow2(0, n)),
        1 => (-200i64..-40, any::<bool>()).prop_map(move |(k, n)| pow2(k, n)),
    ]
    .boxed()
}

/// A part below `2^20` in magnitude, for large imaginary parts (MPC's
/// `tan` of `a + 2^20 i` still stays in MPFR's exponent range).
fn large_part() -> BoxedStrategy<Dyadic> {
    (crate::testkit::strategy::signed_up_to(1), 0i64..=20)
        .prop_map(|(m, e)| Dyadic::new(&m.to_rug(), &rug::Integer::from(e - 64)))
        .boxed()
}

/// `x` as an MPC-ready dyadic.
fn dyadic(x: &Float) -> Dyadic {
    x.to_rug()
}

proptest! {
    #[test]
    fn against_mpc(a in part(), b in part(), prec in prec()) {
        check_point(&a, &b, prec)?;
    }

    #[test]
    fn on_and_near_the_cuts(a in edge_part(), b in edge_part(), prec in 2u64..300) {
        check_point(&a, &b, prec)?;
        check_point(&b, &a, prec)?;
    }

    #[test]
    fn large_imaginary_parts(a in part(), b in large_part(), prec in 2u64..300) {
        let z = ComplexBall::new(ball(&a), ball(&b));
        let zr = mp(&a, &b, 2000);
        for (name, ours, oracle) in &FUNCTIONS[..2] {
            check(name, &ours(&z, prec), prec, &oracle(&zr, prec))?;
            // i z exactly (an owned operand would set the product's precision).
            let izr = MpComplex::with_val(2000, (-zr.imag().clone(), zr.real()));
            check(name, &ours(&z.mul_i(), prec), prec, &oracle(&izr, prec))?;
        }
    }

    /// An inexact ball contains the values at its corners and centre,
    /// which for a ball across a cut lie on both sides.
    #[test]
    fn balls_contain_their_points(
        a in edge_part(), b in edge_part(),
        ra in prop_oneof![Just(None), (-140i64..-50).prop_map(Some)],
        rb in prop_oneof![Just(None), (-140i64..-50).prop_map(Some)],
        prec in 2u64..200,
    ) {
        let rad = |k: Option<i64>| k.map_or(Mag::ZERO, |k| Mag::power_of_two(&Integer::from(k)));
        let (ra, rb) = (rad(ra), rad(rb));
        let z = ComplexBall::new(Ball::new(a.to_ours(), ra.clone()), Ball::new(b.to_ours(), rb.clone()));
        let (fa, fb) = (ra.to_float().unwrap(), rb.to_float().unwrap());
        let corners = [(0, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)];
        let mut results: Vec<(&str, ComplexBall, Oracle)> = FUNCTIONS
            .iter()
            .map(|(name, ours, oracle)| (*name, ours(&z, prec), *oracle))
            .collect();
        results.push(("cos_pi", z.cos_pi(prec), |z, p| pi_oracle(z, p, true)));
        results.push(("root(3)", z.root(3, prec), |z, p| root_oracle(z, p, 3)));
        for (sa, sb) in corners {
            let shift = |c: Float, side: i32, by: &Float| match side.cmp(&0) {
                std::cmp::Ordering::Equal => c,
                std::cmp::Ordering::Greater => c.add(by),
                std::cmp::Ordering::Less => c.sub(by),
            };
            let (pa, pb) = (shift(a.to_ours(), sa, &fa), shift(b.to_ours(), sb, &fb));
            let zr = mp(&dyadic(&pa), &dyadic(&pb), 2000);
            for (name, w, oracle) in &results {
                if !w.is_finite() || singular(name, &zr) {
                    continue;
                }
                let (bre, bim) = oracle(&zr, prec);
                for (part, br) in [(w.re(), bre), (w.im(), bim)] {
                    let mid = part.mid().to_rational().to_rug();
                    let r = part.rad().to_float().unwrap().to_rational().to_rug();
                    check_encloses(&mid, &r, &br).map_err(|e| TestCaseError::fail(format!("{name}({z:?}) at {zr}: {e}")))?;
                }
            }
        }
    }

    #[test]
    fn identities(a in part(), b in part(), prec in 2u64..300) {
        let z = ComplexBall::new(ball(&a), ball(&b));
        let exact = Complex::new(a.to_ours().to_rational(), b.to_ours().to_rational());
        let unit = Complex::from(Rational::from(1i64));
        let wp = prec + 32;
        // cosh^2 - sinh^2 = 1, tan = sin / cos, tanh = sinh / cosh.
        let (ch, sh) = (z.cosh(prec), z.sinh(prec));
        prop_assert!(ch.mul(&ch, wp).sub(&sh.mul(&sh, wp), wp).contains_complex(&unit));
        let quotient = z.sin(wp).div(&z.cos(wp), wp);
        prop_assert!(z.tan(prec).overlaps(&quotient) || !quotient.is_finite(), "tan {:?}", z);
        prop_assert!(z.tanh(prec).overlaps(&sh.div(&ch, wp)) || !ch.is_finite(), "tanh {:?}", z);
        // f(f^-1(z)) = z wherever the inverse is defined.
        for (f, g) in [
            (ComplexBall::sin as Ours, ComplexBall::asin as Ours),
            (ComplexBall::cos, ComplexBall::acos),
            (ComplexBall::sinh, ComplexBall::asinh),
            (ComplexBall::cosh, ComplexBall::acosh),
            (ComplexBall::tanh, ComplexBall::atanh),
            (ComplexBall::expm1, ComplexBall::log1p),
        ] {
            let inv = g(&z, prec);
            if inv.is_finite() {
                prop_assert!(f(&inv, prec).contains_complex(&exact), "f(f^-1({:?})) = {:?}", z, f(&inv, prec));
            }
        }
        let r = z.root(5, prec);
        prop_assert!(r.pow_int(5, prec).contains_complex(&exact), "root(5)^5 of {:?}", z);
        // sin_pi(z) = sin(π z), cos_pi likewise.
        let pz = z.mul(&ComplexBall::from_real(Ball::const_pi(wp)), wp);
        prop_assert!(z.sin_pi(prec).overlaps(&pz.sin(wp)) && z.cos_pi(prec).overlaps(&pz.cos(wp)));
    }

    #[test]
    fn real_axis(a in part(), prec in 2u64..300) {
        let x = ball(&a);
        let z = ComplexBall::from_real(x.clone());
        let one = Ball::from_integer(&Integer::ONE);
        let mut cases: Vec<(&str, ComplexBall, Ball)> = vec![
            ("tan", z.tan(prec), x.tan(prec)),
            ("tanh", z.tanh(prec), x.tanh(prec)),
            ("sinh", z.sinh(prec), x.sinh(prec)),
            ("cosh", z.cosh(prec), x.cosh(prec)),
            ("asinh", z.asinh(prec), x.asinh(prec)),
            ("expm1", z.expm1(prec), x.expm1(prec)),
            ("sin_pi", z.sin_pi(prec), x.sin_pi(prec)),
            ("cos_pi", z.cos_pi(prec), x.cos_pi(prec)),
        ];
        let inside = x.abs().cmp(&one) == Some(std::cmp::Ordering::Less);
        if inside || x.abs() == one {
            cases.push(("asin", z.asin(prec), x.asin(prec)));
            cases.push(("acos", z.acos(prec), x.acos(prec)));
        }
        if inside {
            cases.push(("atanh", z.atanh(prec), x.atanh(prec)));
        }
        if x.cmp(&one).is_some_and(|o| o != std::cmp::Ordering::Less) {
            cases.push(("acosh", z.acosh(prec), x.acosh(prec)));
        }
        if x.sign() != Some(std::cmp::Ordering::Less) {
            cases.push(("root", z.root(3, prec), x.root(3, prec)));
            cases.push(("log1p", z.log1p(prec), x.log1p(prec)));
        }
        for (name, c, r) in cases {
            prop_assert!(c.is_real() && *c.re() == r, "{}({:?}) = {:?}, real {:?}", name, x, c, r);
        }
    }
}

#[test]
fn cut_sides() {
    let int = |n: i64| Ball::from_integer(&Integer::from(n));
    let tiny = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(-60)));
    let across = Ball::new(Float::ZERO, Mag::power_of_two(&Integer::from(-60)));
    let hp = Ball::const_pi(128).mul_2exp(&-Integer::ONE);
    let pi = Ball::const_pi(128);
    let near = |v: &Ball, w: &Ball| {
        v.sub(w, 128).abs().mid() < &Float::from(1i64).mul_2exp(&Integer::from(-40))
    };
    // asin(2 + 0i) = π/2 + i acosh 2 (the +0 side), conjugate below,
    // both across.
    let acosh2 = int(2).acosh(128);
    for (im, sign) in [(Ball::ZERO, 1), (tiny.clone(), 1), (tiny.neg(), -1)] {
        let v = ComplexBall::new(int(2), im.clone()).asin(100);
        assert!(
            near(v.re(), &hp) && near(v.im(), &acosh2.mul(&int(sign), 128)),
            "asin(2 + {im:?} i) = {v:?}"
        );
        let w = ComplexBall::new(int(2), im).acos(100);
        assert!(near(w.im(), &acosh2.mul(&int(-sign), 128)), "acos: {w:?}");
    }
    let v = ComplexBall::new(int(2), across.clone()).asin(100);
    assert!(
        v.im().contains(&Float::from_f64(1.3).unwrap())
            && v.im().contains(&Float::from_f64(-1.3).unwrap())
    );
    // acosh(-2 + 0i) = acosh 2 + iπ, -iπ below, both across.
    for (im, sign) in [(Ball::ZERO, 1), (tiny.neg(), -1)] {
        let v = ComplexBall::new(int(-2), im).acosh(100);
        assert!(
            near(v.re(), &acosh2) && near(v.im(), &pi.mul(&int(sign), 128)),
            "acosh: {v:?}"
        );
    }
    let v = ComplexBall::new(int(-2), across.clone()).acosh(100);
    assert!(
        v.im().contains(&Float::from_f64(3.1).unwrap())
            && v.im().contains(&Float::from_f64(-3.1).unwrap())
    );
    // asinh(0 + 2i): real part > 0 (the side of Re z > 0), left of the cut < 0.
    for (re, sign) in [(Ball::ZERO, 1), (tiny.neg(), -1)] {
        let v = ComplexBall::new(re, int(2)).asinh(100);
        assert!(
            near(v.re(), &acosh2.mul(&int(sign), 128)) && near(v.im(), &hp),
            "asinh: {v:?}"
        );
    }
    let v = ComplexBall::new(across.clone(), int(2)).asinh(100);
    assert!(
        v.re().contains(&Float::from_f64(1.3).unwrap())
            && v.re().contains(&Float::from_f64(-1.3).unwrap())
    );
    // atanh(±2 + 0i): imaginary part +π/2, -π/2 below.
    for x in [2, -2] {
        for (im, sign) in [(Ball::ZERO, 1), (tiny.neg(), -1)] {
            let v = ComplexBall::new(int(x), im).atanh(100);
            assert!(near(v.im(), &hp.mul(&int(sign), 128)), "atanh({x}): {v:?}");
        }
        let v = ComplexBall::new(int(x), across.clone()).atanh(100);
        assert!(
            v.im().contains(&Float::from_f64(1.5).unwrap())
                && v.im().contains(&Float::from_f64(-1.5).unwrap())
        );
    }
    // log1p(-3 + 0i) = log 2 + iπ; root(-8, 3) = 1 + i√3, the principal root.
    let v = ComplexBall::new(int(-3), Ball::ZERO).log1p(100);
    assert!(near(v.re(), &Ball::const_log2(128)) && near(v.im(), &pi));
    let v = ComplexBall::from_real(int(-8)).root(3, 100);
    assert!(
        near(v.re(), &int(1)) && near(v.im(), &int(3).sqrt(128)),
        "{v:?}"
    );
    let v = ComplexBall::new(int(-8), across).root(3, 100);
    assert!(
        v.im().contains(&Float::from_f64(1.7).unwrap())
            && v.im().contains(&Float::from_f64(-1.7).unwrap())
    );
    // Poles and singular points: indeterminate.
    assert!(!ComplexBall::from_real(int(1)).atanh(64).is_finite());
    assert!(!ComplexBall::from_real(int(-1)).log1p(64).is_finite());
}

#[test]
fn exact_values() {
    let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
    for k in -6i64..=6 {
        let n = ComplexBall::from_real(Ball::from_integer(&Integer::from(k)));
        let h = ComplexBall::from_real(Ball::exact(Float::from(k).add(&half)));
        let sign = Ball::from_integer(&Integer::from(if k % 2 == 0 { 1 } else { -1 }));
        assert_eq!(n.sin_pi(64), ComplexBall::ZERO);
        assert_eq!(n.cos_pi(64), ComplexBall::from_real(sign.clone()));
        assert_eq!(h.sin_pi(64), ComplexBall::from_real(sign));
        assert_eq!(h.cos_pi(64), ComplexBall::ZERO);
        // sin(π(k + bi)) is purely imaginary.
        let w = ComplexBall::new(
            Ball::from_integer(&Integer::from(k)),
            Ball::exact(half.clone()),
        );
        assert_eq!(w.sin_pi(64).re(), &Ball::ZERO);
        assert_eq!(w.cos_pi(64).im(), &Ball::ZERO);
    }
    assert_eq!(ComplexBall::ZERO.root(5, 64), ComplexBall::ZERO);
    assert_eq!(
        ComplexBall::from_real(Ball::from_integer(&Integer::from(32))).root(5, 64),
        ComplexBall::from_real(Ball::from_integer(&Integer::from(2)))
    );
    assert_eq!(ComplexBall::ZERO.asin(64), ComplexBall::ZERO);
    assert_eq!(ComplexBall::ZERO.atanh(64), ComplexBall::ZERO);
    assert_eq!(ComplexBall::ZERO.expm1(64), ComplexBall::ZERO);
}

#[test]
fn tan_far_from_the_real_axis() {
    // |Im z| = 2^70000, beyond exp's range: ±i within 2^(4 - 2^31).
    let huge = Ball::exact(Float::from(1i64).mul_2exp(&Integer::from(70_000)));
    let tiny = Mag::power_of_two(&Integer::from(-1_000_000));
    for (im, sign) in [(huge.clone(), 1i64), (huge.neg(), -1)] {
        let t = ComplexBall::new(Ball::from_integer(&Integer::from(3)), im).tan(200);
        assert!(
            *t.re().rad() < tiny && t.re().contains(&Float::ZERO),
            "{t:?}"
        );
        assert!(
            t.im().contains(&Float::from(sign)) && *t.im().rad() < tiny,
            "{t:?}"
        );
    }
    // A ball around the pole π/2: indeterminate; next to it, finite.
    let around = |m: f64, k: i64| {
        Ball::new(
            Float::from_f64(m).unwrap(),
            Mag::power_of_two(&Integer::from(k)),
        )
    };
    let pole = ComplexBall::new(around(1.5, -3), around(0.0, -4));
    assert!(!pole.tan(64).is_finite() && !pole.mul_i().tanh(64).is_finite());
    let beside = ComplexBall::new(around(1.25, -4), around(0.0, -4));
    assert!(beside.tan(64).is_finite() && beside.mul_i().tanh(64).is_finite());
}

#[test]
fn narrow_balls_are_tight() {
    type Deriv = fn(&ComplexBall, &ComplexBall, &ComplexBall, u64) -> ComplexBall;
    // An input radius ρ = 2^-60 in each part moves the value by about |f'|
    // √2 ρ; the radius must stay within twice that, with 2^-50 for |f'| over
    // the rectangle beyond its value at the midpoint (|f''| ρ √2, the whole
    // term at a critical point such as sin_pi(3/2)). Composing the formulas
    // instead gave 4000 |f'| ρ for tan(-0.3 + 4i) and 47000 |f'| ρ for cosh
    // near 0.
    let prec = 200;
    let rho = Mag::power_of_two(&Integer::from(-60));
    let ex = |v: f64| Ball::exact(Float::from_f64(v).unwrap());
    let pi = ComplexBall::from_real(Ball::const_pi(prec));
    let one = ComplexBall::from_real(Ball::from_integer(&Integer::ONE));
    let fs: [(&str, Ours, Deriv); 6] = [
        ("tan", ComplexBall::tan, |z, one, _, p| {
            one.add(&z.tan(p).pow_int(2, p), p)
        }),
        ("tanh", ComplexBall::tanh, |z, one, _, p| {
            one.sub(&z.tanh(p).pow_int(2, p), p)
        }),
        ("sinh", ComplexBall::sinh, |z, _, _, p| z.cosh(p)),
        ("cosh", ComplexBall::cosh, |z, _, _, p| z.sinh(p)),
        ("sin_pi", ComplexBall::sin_pi, |z, _, pi, p| {
            pi.mul(&z.cos_pi(p), p)
        }),
        ("cos_pi", ComplexBall::cos_pi, |z, _, pi, p| {
            pi.mul(&z.sin_pi(p), p)
        }),
    ];
    for (x, y) in [
        (0.7, 0.3),
        (2.5, 1e-3),
        (-0.3, 4.0),
        (1e-5, 2e-5),
        (0.999, 1e-4),
        (1.5, 0.0),
    ] {
        let z = ComplexBall::new(ex(x), ex(y));
        let ball = ComplexBall::new(
            Ball::new(z.re().mid().clone(), rho.clone()),
            Ball::new(z.im().mid().clone(), rho.clone()),
        );
        for (name, f, deriv) in fs {
            let d = deriv(&z, &one, &pi, prec).abs(64);
            let slope = Mag::from_float(d.mid())
                .add(d.rad())
                .add(&Mag::power_of_two(&Integer::from(-50)));
            let bound = slope.mul(&rho).mul_2exp(&Integer::from(2));
            let w = f(&ball, prec);
            assert!(
                *w.re().rad() <= bound && *w.im().rad() <= bound,
                "{name}({x} + {y}i): {w:?}, |f'| = {d:?}"
            );
        }
    }
}

#[test]
fn straddling_balls_stay_finite() {
    // A real part across 0 (atanh folds it onto [0, 1/2], not [-1/2, 1/2],
    // whose 4X / |1 - z|^2 would reach below -1) and across 1 (|1 - z|^2
    // squares a ball across 0 as [0, 1/4], not [-1/4, 1/4]): both finite
    // and containing the values at their edges.
    let half = Mag::power_of_two(&-Integer::ONE);
    let quarter = Float::from(1i64).mul_2exp(&Integer::from(-2));
    for centre in [0i64, 1] {
        let z = ComplexBall::new(
            Ball::new(Float::from(centre), half.clone()),
            Ball::exact(quarter.clone()),
        );
        let w = z.atanh(64);
        assert!(w.is_finite(), "atanh({z:?}) = {w:?}");
        // Re atanh(±1/2 + i/4) = ±log1p(6.4)/4 = ±0.5004 at the edges of
        // the first ball. Folding x onto [0, 1/2] gives the quotient [0,
        // 10.7] and a radius of 1.8 (log1p of that ball bounds its slope at
        // the lower end); without the fold the quotient is [0, 32] and the
        // radius 4.7.
        if centre == 0 {
            assert!(
                *w.re().rad() <= Mag::from_integer(&Integer::from(2)),
                "{w:?}"
            );
        }
        for edge in [-1i64, 0, 1] {
            // centre + edge/2 + i/4, with the pole 1 + i/4 excluded by the
            // imaginary part.
            let re = Dyadic::new(
                &rug::Integer::from(2 * centre + edge),
                &rug::Integer::from(-1),
            );
            let im = Dyadic::new(&rug::Integer::from(1), &rug::Integer::from(-2));
            let zr = mp(&re, &im, 200);
            let br = brackets(64, |p, r| MpComplex::with_val_round(p, zr.atanh_ref(), r).0);
            for (part, b) in [(w.re(), &br.0), (w.im(), &br.1)] {
                let mid = part.mid().to_rational().to_rug();
                let rad = part.rad().to_float().unwrap().to_rational().to_rug();
                check_encloses(&mid, &rad, b).unwrap();
            }
        }
    }
}

#[test]
fn found_by_proptest() {
    // Inputs the proptests found, kept as fixed cases: atanh next to |z| =
    // 1, whose imaginary part had 2^44 radii before atan2 took x / y for
    // |y| >= |x|, and cos_pi next to a half-integer, before its exact fold.
    let d = |m: &str, e: i64| {
        Dyadic::new(
            &rug::Integer::from_str_radix(m, 16).unwrap(),
            &rug::Integer::from(e),
        )
    };
    for (a, b, prec) in [
        (
            d("139d1", -31),
            d("80000000000000000000000000000001", -127),
            15,
        ),
        (d("7ffffffffff", -61), d("1", 0), 21),
        (d("0", 0), d("7fffffffffffffffffffffffffffffff", -128), 6),
    ] {
        check_point(&a, &b, prec).unwrap();
        check_point(&b, &a, prec).unwrap();
    }
}
