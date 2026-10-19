//! Machine balls against MPFR (CAS-113). `None` is always allowed: the
//! caller redoes the operation with `Ball`. A `Some` must meet MPFR's
//! bracket at 192 bits at the midpoint, both ends and an interior point of
//! every input (every pair of them for the binary operations), and a point
//! where MPFR has no finite value (outside the domain, at a pole, or past
//! MPFR's exponents, far beyond `f64`'s) must give `None`.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::Integer;
use cas_num::ball::{Float, MachineBall};
use cas_num_fuzz::ball::{Oracle, bracket, to_mp, warm_up};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;
use rug::float::{Constant, Round, Special};
use rug::ops::Pow;

type Unary = fn(&MachineBall) -> Option<MachineBall>;
type Binary = fn(&MachineBall, &MachineBall) -> Option<MachineBall>;
type Mpfr = fn(&MpFloat, u32, Round) -> MpFloat;
type Mpfr2 = fn(&MpFloat, &MpFloat, u32, Round) -> MpFloat;

/// MPFR's correctly rounded `x.$f()`.
macro_rules! mpfr {
    ($f:ident) => {
        |x: &MpFloat, p: u32, r: Round| MpFloat::with_val_round(p, x.$f(), r).0
    };
}

const UNARY: [(&str, Unary, Mpfr); 20] = [
    ("exp", MachineBall::exp, mpfr!(exp_ref)),
    ("expm1", MachineBall::expm1, mpfr!(exp_m1_ref)),
    ("log", MachineBall::log, mpfr!(ln_ref)),
    ("log1p", MachineBall::log1p, mpfr!(ln_1p_ref)),
    ("sin", MachineBall::sin, mpfr!(sin_ref)),
    ("cos", MachineBall::cos, mpfr!(cos_ref)),
    ("tan", MachineBall::tan, mpfr!(tan_ref)),
    ("atan", MachineBall::atan, mpfr!(atan_ref)),
    ("asin", MachineBall::asin, mpfr!(asin_ref)),
    ("acos", MachineBall::acos, mpfr!(acos_ref)),
    ("sinh", MachineBall::sinh, mpfr!(sinh_ref)),
    ("cosh", MachineBall::cosh, mpfr!(cosh_ref)),
    ("tanh", MachineBall::tanh, mpfr!(tanh_ref)),
    ("asinh", MachineBall::asinh, mpfr!(asinh_ref)),
    ("acosh", MachineBall::acosh, mpfr!(acosh_ref)),
    ("atanh", MachineBall::atanh, mpfr!(atanh_ref)),
    ("sin_pi", MachineBall::sin_pi, mpfr!(sin_pi_ref)),
    ("cos_pi", MachineBall::cos_pi, mpfr!(cos_pi_ref)),
    ("sqrt", MachineBall::sqrt, mpfr!(sqrt_ref)),
    ("inv", MachineBall::inv, |x, p, r| {
        quotient(&MpFloat::with_val(1, 1u32), x, p, r)
    }),
];

const BINARY: [(&str, Binary, Mpfr2); 6] = [
    ("add", MachineBall::add, |a, b, p, r| {
        MpFloat::with_val_round(p, a + b, r).0
    }),
    ("sub", MachineBall::sub, |a, b, p, r| {
        MpFloat::with_val_round(p, a - b, r).0
    }),
    ("mul", MachineBall::mul, |a, b, p, r| {
        MpFloat::with_val_round(p, a * b, r).0
    }),
    ("div", MachineBall::div, quotient),
    ("atan2", MachineBall::atan2, |y, x, p, r| {
        MpFloat::with_val_round(p, y.atan2_ref(x), r).0
    }),
    ("pow", MachineBall::pow_ball, |a, b, p, r| {
        MpFloat::with_val_round(p, a.pow(b), r).0
    }),
];

/// `a / b` rounded once from the exact rational quotient, not by
/// `mpfr_div`, which misrounds some halfway quotients (see the test kit's
/// `div_round`); an infinity or NaN for `b = 0`.
fn quotient(a: &MpFloat, b: &MpFloat, p: u32, r: Round) -> MpFloat {
    if b.is_zero() {
        let s = if a.is_zero() {
            Special::Nan
        } else {
            Special::Infinity
        };
        return MpFloat::with_val(p, s);
    }
    let q = a.to_rational().expect("finite") / b.to_rational().expect("finite");
    MpFloat::with_val_round(p, &q, r).0
}

fuzz_target!(|data: &[u8]| {
    // The fallbacks to `Ball` read the caches too.
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let op = u.int_in_range(0..=UNARY.len() + BINARY.len() + 1)?;
    let Some(x) = ball(u)? else { return Ok(()) };
    if let Some(&(name, f, mp)) = UNARY.get(op) {
        let y = f(&x);
        for p in points(u, &x)? {
            let answer = at(&p, mp);
            check(y, &answer, || format!("{name}({x:?}) = {y:?}, at {p:?}"));
        }
        return Ok(());
    }
    match op - UNARY.len() {
        k if k < BINARY.len() => {
            let (name, f, mp) = BINARY[k];
            let Some(w) = ball(u)? else { return Ok(()) };
            let y = f(&x, &w);
            let (px, pw) = (points(u, &x)?, points(u, &w)?);
            for a in &px {
                for b in &pw {
                    let (Some(am), Some(bm)) = (to_mp(a), to_mp(b)) else {
                        continue;
                    };
                    let answer = bracket(192, |q, r| mp(&am, &bm, q, r));
                    check(y, &answer, || {
                        format!("{name}({x:?}, {w:?}) = {y:?}, at {a:?}, {b:?}")
                    });
                }
            }
        }
        k if k == BINARY.len() => {
            let n: i32 = if u.ratio(1u8, 4)? {
                u.arbitrary()?
            } else {
                u.int_in_range(-40..=40)?
            };
            let y = x.pow(n);
            for p in points(u, &x)? {
                let answer = at(&p, |pm, q, r| MpFloat::with_val_round(q, pm.pow(n), r).0);
                check(y, &answer, || format!("pow({x:?}, {n}) = {y:?}, at {p:?}"));
            }
        }
        _ => {
            let n: u32 = if u.ratio(1u8, 4)? {
                u.int_in_range(17..=1 << 20)?
            } else {
                u.int_in_range(1..=16)?
            };
            let y = x.root(u64::from(n));
            for p in points(u, &x)? {
                let answer = at(&p, |pm, q, r| {
                    MpFloat::with_val_round(q, pm.root_ref(n), r).0
                });
                check(y, &answer, || format!("root({x:?}, {n}) = {y:?}, at {p:?}"));
            }
        }
    }
    Ok(())
}

/// MPFR's bracket of `f` at the point `p`, at 192 bits.
fn at(p: &Float, f: impl Fn(&MpFloat, u32, Round) -> MpFloat) -> Oracle {
    let pm = to_mp(p).expect("f64 points are in MPFR's range");
    bracket(192, |q, r| f(&pm, q, r))
}

/// `None` is always allowed; a `Some` must meet MPFR's bracket, and MPFR
/// must have a finite value.
fn check(y: Option<MachineBall>, answer: &Oracle, what: impl Fn() -> String) {
    let Some(b) = y else { return };
    match answer {
        Oracle::Bracket(lo, hi) => {
            // The machine ball's own ends, exactly.
            let (m, r) = (exact(b.mid()), exact(b.rad()));
            assert!(
                *lo <= m.add(&r) && m.sub(&r) <= *hi,
                "{}: misses [{lo:?}, {hi:?}]",
                what()
            );
        }
        Oracle::Nan | Oracle::Infinite => panic!("{}: MPFR has no finite value", what()),
    }
}

fn exact(x: f64) -> Float {
    Float::from_f64(x).expect("finite")
}

/// A machine ball, or `None` where `MachineBall::new` refuses the parts
/// (subnormal, infinite or NaN).
fn ball(u: &mut Unstructured<'_>) -> Result<Option<MachineBall>> {
    let m = double(u)?;
    let r = radius(u, m)?;
    Ok(MachineBall::new(m, r))
}

/// A double: small integers and halves, powers of two, values next to 1
/// and to multiples of π/2, random mantissas near 1, or raw bits.
fn double(u: &mut Unstructured<'_>) -> Result<f64> {
    let sign = if u.arbitrary()? { -1.0 } else { 1.0 };
    Ok(match u.int_in_range(0u8..=7)? {
        0 => f64::from(u.int_in_range(-2000i32..=2000)?) / 2.0,
        1 => sign * 2f64.powi(u.int_in_range(-1022..=1023)?),
        2 => sign * (1.0 + f64::from(u.int_in_range(-8i32..=8)?) * f64::EPSILON),
        3 => {
            let k = u.int_in_range(-(1i64 << 22)..=1 << 22)?;
            let x = (MpFloat::with_val(200, Constant::Pi) * k / 2u32).to_f64();
            let step: i64 = u.int_in_range(-2..=2)?;
            f64::from_bits(x.to_bits().wrapping_add_signed(step))
        }
        4 => f64::from_bits(u.arbitrary()?),
        _ => {
            let e = u64::try_from(1023 + u.int_in_range(-60i64..=60)?).expect("biased");
            sign * f64::from_bits(e << 52 | (u.arbitrary::<u64>()? >> 12))
        }
    })
}

/// A radius for the midpoint `m`: zero, relative, a few ulps, absolute, or
/// raw bits.
fn radius(u: &mut Unstructured<'_>, m: f64) -> Result<f64> {
    Ok(match u.int_in_range(0u8..=7)? {
        0..=2 => 0.0,
        3 | 4 => m.abs() * 2f64.powi(-u.int_in_range(1..=60)?),
        5 => m.abs() * f64::EPSILON * f64::from(u.int_in_range(1u8..=8)?),
        6 => 2f64.powi(u.int_in_range(-1022..=10)?),
        _ => f64::from_bits(u.arbitrary()?).abs(),
    })
}

/// The midpoint, both ends and `m + t r` for a fuzzer-chosen `t` in
/// `[-1, 1]`, exactly (from the `f64` radius itself: `to_ball` would round
/// it up, past the ball).
fn points(u: &mut Unstructured<'_>, b: &MachineBall) -> Result<Vec<Float>> {
    let (m, r) = (exact(b.mid()), exact(b.rad()));
    if r.is_zero() {
        return Ok(vec![m]);
    }
    let t = Float::from(u.int_in_range(-(1i64 << 16)..=1 << 16)?).mul_2exp(&Integer::from(-16i64));
    let inner = m.add(&r.mul(&t));
    Ok(vec![m.clone(), m.sub(&r), m.add(&r), inner])
}
