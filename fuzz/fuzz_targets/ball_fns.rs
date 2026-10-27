//! Real ball functions against MPFR (CAS-113). At the midpoint, both ends
//! and an interior point of every input ball (every pair of them for
//! `atan2` and `pow_ball`), the result must meet MPFR's bracket at
//! `prec + 64` bits, and a point outside the domain or at a pole must make
//! it the indeterminate ball. An exact input inside the domain gives a
//! finite result with at least `prec - 3` bits, the property tests' bound.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::Integer;
use cas_num::ball::{Ball, Float};
use cas_num_fuzz::ball::{
    Oracle, Range, bracket, float, meets, points, precision, radius, to_mp, warm_up,
};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;
use rug::float::Round;
use rug::ops::Pow;

type Ours = fn(&Ball, u64) -> Ball;
type Mpfr = fn(&MpFloat, u32, Round) -> MpFloat;

/// MPFR's correctly rounded `x.$f()`.
macro_rules! mpfr {
    ($f:ident) => {
        |x: &MpFloat, p: u32, r: Round| MpFloat::with_val_round(p, x.$f(), r).0
    };
}

const UNARY: [(&str, Ours, Mpfr); 20] = [
    ("exp", Ball::exp, mpfr!(exp_ref)),
    ("expm1", Ball::expm1, mpfr!(exp_m1_ref)),
    ("log", Ball::log, mpfr!(ln_ref)),
    ("log1p", Ball::log1p, mpfr!(ln_1p_ref)),
    ("sin", Ball::sin, mpfr!(sin_ref)),
    ("cos", Ball::cos, mpfr!(cos_ref)),
    ("tan", Ball::tan, mpfr!(tan_ref)),
    ("sin_pi", Ball::sin_pi, mpfr!(sin_pi_ref)),
    ("cos_pi", Ball::cos_pi, mpfr!(cos_pi_ref)),
    ("atan", Ball::atan, mpfr!(atan_ref)),
    ("asin", Ball::asin, mpfr!(asin_ref)),
    ("acos", Ball::acos, mpfr!(acos_ref)),
    ("sinh", Ball::sinh, mpfr!(sinh_ref)),
    ("cosh", Ball::cosh, mpfr!(cosh_ref)),
    ("tanh", Ball::tanh, mpfr!(tanh_ref)),
    ("asinh", Ball::asinh, mpfr!(asinh_ref)),
    ("acosh", Ball::acosh, mpfr!(acosh_ref)),
    ("atanh", Ball::atanh, mpfr!(atanh_ref)),
    ("sqrt", Ball::sqrt, mpfr!(sqrt_ref)),
    ("gamma", Ball::gamma, mpfr!(gamma_ref)),
];

/// `sin` and `cos` give `[0 ± 1]` for arguments of exponent above this
/// (`fns/trig.rs`), and `tan` the indeterminate ball.
const HUGE_EXP: i64 = 1 << 16;

fuzz_target!(|data: &[u8]| {
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let op = u.int_in_range(0..=UNARY.len() + 3)?;
    let prec = precision(u)?;
    let x = ball(u, prec)?;
    match op.checked_sub(UNARY.len()) {
        None => unary(u, UNARY[op], &x, prec),
        Some(0) => root(u, &x, prec),
        Some(1) => pow_int(u, &x, prec),
        Some(k) => {
            let y = ball(u, prec)?;
            binary(u, k == 2, &x, &y, prec)
        }
    }
}

fn ball(u: &mut Unstructured<'_>, prec: u64) -> Result<Ball> {
    let m = float(u, Range::Wide)?;
    let r = radius(u, &m, prec, Range::Wide)?;
    Ok(Ball::new(m, r))
}

/// Checks `y` against MPFR's answer at one point of the input; `pole` tells
/// a pole's infinity from an overflow past MPFR's exponent range.
fn check(y: &Ball, answer: &Oracle, pole: bool, what: impl Fn() -> String) {
    match answer {
        Oracle::Bracket(lo, hi) => {
            assert!(
                meets(y, lo, hi),
                "{}: {y:?} misses [{lo:?}, {hi:?}]",
                what()
            );
        }
        Oracle::Nan => assert_eq!(*y, Ball::INDETERMINATE, "{}: outside the domain", what()),
        Oracle::Infinite if pole => assert_eq!(*y, Ball::INDETERMINATE, "{}: a pole", what()),
        Oracle::Infinite => {}
    }
}

/// For an exact input where MPFR has a value: a finite result and, when
/// `accurate`, at least `prec - 3` bits.
fn check_exact(y: &Ball, answer: &Oracle, prec: u64, accurate: bool, what: impl Fn() -> String) {
    if let Oracle::Bracket(..) = answer {
        assert!(y.is_finite(), "{}: indeterminate", what());
        let want = i64::try_from(prec).expect("small precision") - 3;
        if accurate && !y.mid().is_zero() {
            let bits = y.rel_accuracy_bits();
            assert!(bits >= want, "{}: {y:?} has {bits} bits", what());
        }
    }
}

/// Poles at exact points, where MPFR's infinities are not overflows.
fn pole(name: &str, x: &Float) -> bool {
    let at = |v: i64| *x == Float::from(v);
    match name {
        "log" => x.is_zero(),
        "log1p" => at(-1),
        "atanh" => at(1) || at(-1),
        // MPFR's gamma is ±∞ at ±0, NaN at the negative integers.
        "gamma" => x.is_zero(),
        _ => false,
    }
}

fn unary(
    u: &mut Unstructured<'_>,
    (name, f, mp): (&str, Ours, Mpfr),
    x: &Ball,
    prec: u64,
) -> Result<()> {
    let y = f(x, prec);
    for p in points(u, x, prec)? {
        let Some(pm) = to_mp(&p) else { continue };
        let answer = bracket(prec + 64, |q, r| mp(&pm, q, r));
        let what = || format!("{name}({x:?}) at {prec} bits, at {p:?}");
        check(&y, &answer, pole(name, &p), what);
        let huge = matches!(name, "sin" | "cos" | "tan") && *p.exponent() > Integer::from(HUGE_EXP);
        if x.is_exact() && !huge {
            check_exact(&y, &answer, prec, true, what);
        }
    }
    Ok(())
}

fn root(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let n: u32 = if u.ratio(1u8, 4)? {
        u.int_in_range(17..=1 << 20)?
    } else {
        u.int_in_range(1..=16)?
    };
    let y = x.root(u64::from(n), prec);
    for p in points(u, x, prec)? {
        let Some(pm) = to_mp(&p) else { continue };
        let answer = bracket(prec + 64, |q, r| {
            MpFloat::with_val_round(q, pm.root_ref(n), r).0
        });
        let what = || format!("root({x:?}, {n}) at {prec} bits, at {p:?}");
        check(&y, &answer, false, what);
        if x.is_exact() {
            check_exact(&y, &answer, prec, true, what);
        }
    }
    Ok(())
}

/// `Ball::pow`, by squaring: its rounding errors grow with `log2 n`, so
/// only containment and finiteness are checked.
fn pow_int(u: &mut Unstructured<'_>, x: &Ball, prec: u64) -> Result<()> {
    let n: i64 = match u.int_in_range(0u8..=7)? {
        0 => u.arbitrary()?,
        1 => u.int_in_range(-(1 << 20)..=1 << 20)?,
        _ => u.int_in_range(-64..=64)?,
    };
    let y = x.pow(n, prec);
    for p in points(u, x, prec)? {
        let Some(pm) = to_mp(&p) else { continue };
        let answer = bracket(prec + 64, |q, r| {
            MpFloat::with_val_round(q, (&pm).pow(n), r).0
        });
        let what = || format!("pow({x:?}, {n}) at {prec} bits, at {p:?}");
        check(&y, &answer, p.is_zero() && n < 0, what);
        if x.is_exact() {
            check_exact(&y, &answer, prec, false, what);
        }
    }
    Ok(())
}

/// `atan2(x, y)` (`x` the ordinate) or `pow_ball(x, y)`, at every pair of
/// points of the two balls.
fn binary(u: &mut Unstructured<'_>, atan2: bool, x: &Ball, y: &Ball, prec: u64) -> Result<()> {
    let (name, z) = if atan2 {
        ("atan2", x.atan2(y, prec))
    } else {
        ("pow", x.pow_ball(y, prec))
    };
    let (px, py) = (points(u, x, prec)?, points(u, y, prec)?);
    for a in &px {
        for b in &py {
            let (Some(am), Some(bm)) = (to_mp(a), to_mp(b)) else {
                continue;
            };
            let answer = bracket(prec + 64, |q, r| {
                if atan2 {
                    MpFloat::with_val_round(q, am.atan2_ref(&bm), r).0
                } else {
                    MpFloat::with_val_round(q, (&am).pow(&bm), r).0
                }
            });
            let what = || format!("{name}({x:?}, {y:?}) at {prec} bits, at {a:?}, {b:?}");
            check(&z, &answer, !atan2 && a.is_zero() && b.is_negative(), what);
            if x.is_exact() && y.is_exact() {
                check_exact(&z, &answer, prec, true, what);
            }
        }
    }
    Ok(())
}
