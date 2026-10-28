//! Complex ball functions against MPC (CAS-113). At the centre, the four
//! corners and an interior point of every input (pairs of them for `mul`,
//! `div` and `pow`), each part of the result must meet the oracle's bracket
//! at `prec + 192` bits, and a pole must make the result non-finite. An
//! exact input where the oracle has a value gives a finite result whose
//! radii are within `2^(6 - prec)` of its larger part, the property tests'
//! bound, for parts below `2^8`.
//!
//! The oracle is MPC's function where there is one and it can be trusted,
//! and quick: `pow_int`'s is the exact power by Gaussian integers (CAS-362).
//! MPC's `log` misses `log|z|` next to `|z| = 1` (see [`mpfr_log`]), so
//! `log`, `log1p`, `pow` and `root` take the logarithm from MPFR instead:
//! `log1p` of the exact `1 + z`, and `pow` and `root` as `exp(w log z)` far
//! above the precision with an explicit margin, as `sin_pi` and `cos_pi`.
//! `expm1` is `exp` minus 1 with the bits the subtraction cancels.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::Integer;
use cas_num::ball::{Ball, ComplexBall, Float, Mag};
use cas_num_fuzz::ball::{
    Oracle, Range, bracket, float, from_mp, meets, points, precision, radius, to_mp, warm_up,
};
use libfuzzer_sys::fuzz_target;
use rug::float::{Constant, Round};
use rug::ops::Pow;
use rug::{Complex as MpComplex, Float as MpFloat};

type Ours = fn(&ComplexBall, u64) -> ComplexBall;
type Mpc = fn(&MpComplex, u32, (Round, Round)) -> MpComplex;

/// MPC's correctly rounded `z.$f()`.
macro_rules! mpc {
    ($f:ident) => {
        |z: &MpComplex, p: u32, r: (Round, Round)| MpComplex::with_val_round(p, z.$f(), r).0
    };
}

const UNARY: [(&str, Ours, Mpc); 16] = [
    ("exp", ComplexBall::exp, mpc!(exp_ref)),
    ("log", ComplexBall::log, mpfr_log),
    ("sqrt", ComplexBall::sqrt, mpc!(sqrt_ref)),
    ("sin", ComplexBall::sin, mpc!(sin_ref)),
    ("cos", ComplexBall::cos, mpc!(cos_ref)),
    ("tan", ComplexBall::tan, mpc!(tan_ref)),
    ("sinh", ComplexBall::sinh, mpc!(sinh_ref)),
    ("cosh", ComplexBall::cosh, mpc!(cosh_ref)),
    ("tanh", ComplexBall::tanh, mpc!(tanh_ref)),
    ("asin", ComplexBall::asin, mpc!(asin_ref)),
    ("acos", ComplexBall::acos, mpc!(acos_ref)),
    ("atan", ComplexBall::atan, mpc!(atan_ref)),
    ("asinh", ComplexBall::asinh, mpc!(asinh_ref)),
    ("acosh", ComplexBall::acosh, mpc!(acosh_ref)),
    ("atanh", ComplexBall::atanh, mpc!(atanh_ref)),
    ("inv", ComplexBall::inv, mpc!(recip_ref)),
];

/// The other operations, after [`UNARY`].
const OTHERS: [&str; 11] = [
    "expm1", "log1p", "sin_pi", "cos_pi", "root", "pow_int", "abs", "arg", "mul", "div", "pow",
];

/// MPC's answer at one point.
#[derive(Clone, Debug)]
#[allow(
    clippy::large_enum_variant,
    reason = "one answer per oracle call, never stored"
)]
enum Answer {
    /// Brackets of the real and imaginary parts.
    Bracket([(Float, Float); 2]),
    /// An infinity or NaN: with the inputs' magnitudes below `2^20`, a pole.
    NonFinite,
    /// No answer to check against.
    Unknown,
}

fuzz_target!(|data: &[u8]| {
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let op = u.int_in_range(0..=UNARY.len() + OTHERS.len() - 1)?;
    let prec = precision(u)?;
    let z = ball(u, prec)?;
    let name = UNARY
        .get(op)
        .map_or_else(|| OTHERS[op - UNARY.len()], |f| f.0);
    match name {
        "abs" | "arg" => real(u, name, &z, prec),
        "mul" | "div" | "pow" => {
            let w = ball(u, prec)?;
            binary(u, name, &z, &w, prec)
        }
        "root" | "pow_int" => {
            let n: i64 = match u.int_in_range(0u8..=3)? {
                0 => u.int_in_range(-(1 << 16)..=1 << 16)?,
                _ => u.int_in_range(-16..=16)?,
            };
            if name == "root" {
                let n = n.unsigned_abs().max(1);
                let w = z.root(n, prec);
                let n = u32::try_from(n).expect("below 2^16");
                unary(u, "root", &z, &w, prec, true, &|zm, prec| {
                    root_oracle(zm, prec, n)
                })
            } else {
                let w = z.pow_int(n, prec);
                let f = move |zm: &MpComplex, prec: u64| power_oracle(zm, n, prec);
                unary(u, "pow_int", &z, &w, prec, false, &f)
            }
        }
        "expm1" => unary(u, name, &z, &z.expm1(prec), prec, true, &expm1_oracle),
        "log1p" => unary(u, name, &z, &z.log1p(prec), prec, true, &log1p_oracle),
        "sin_pi" => unary(u, name, &z, &z.sin_pi(prec), prec, true, &|zm, p| {
            pi_oracle(zm, p, false)
        }),
        "cos_pi" => unary(u, name, &z, &z.cos_pi(prec), prec, true, &|zm, p| {
            pi_oracle(zm, p, true)
        }),
        _ => {
            let (_, f, mp) = UNARY[op];
            let oracle = move |zm: &MpComplex, prec: u64| mpc_bracket(prec, |p, r| mp(zm, p, r));
            unary(u, name, &z, &f(&z, prec), prec, true, &oracle)
        }
    }
}

/// A complex ball with parts below `2^20` and radii up to `2^20`.
fn ball(u: &mut Unstructured<'_>, prec: u64) -> Result<ComplexBall> {
    let mut part = || -> Result<Ball> {
        let m = float(u, Range::Narrow)?;
        let r = radius(u, &m, prec, Range::Narrow)?;
        Ok(Ball::new(m, r))
    };
    let re = part()?;
    Ok(ComplexBall::new(re, part()?))
}

/// The centre, the corners and an interior point of `z`, exactly (each
/// part's midpoint, ends and `m + t r`), without repeats.
fn cpoints(u: &mut Unstructured<'_>, z: &ComplexBall, prec: u64) -> Result<Vec<MpComplex>> {
    let (a, b) = (points(u, z.re(), prec)?, points(u, z.im(), prec)?);
    let mut out: Vec<MpComplex> = Vec::new();
    if a.is_empty() || b.is_empty() {
        return Ok(out);
    }
    for (i, j) in [(0, 0), (1, 1), (1, 2), (2, 1), (2, 2), (3, 3)] {
        let (x, y) = (&a[i.min(a.len() - 1)], &b[j.min(b.len() - 1)]);
        let (Some(x), Some(y)) = (to_mp(x), to_mp(y)) else {
            continue;
        };
        let p = MpComplex::with_val((x.prec(), y.prec()), (x, y));
        if !out.contains(&p) {
            out.push(p);
        }
    }
    Ok(out)
}

/// MPC's brackets from `f` at `prec + 192` bits.
fn mpc_bracket(prec: u64, f: impl Fn(u32, (Round, Round)) -> MpComplex) -> Answer {
    let p = u32::try_from(prec + 192).expect("small precision");
    let lo = f(p, (Round::Down, Round::Down));
    let hi = f(p, (Round::Up, Round::Up));
    let part = |l: &MpFloat, h: &MpFloat| {
        let (l, h) = (from_mp(l)?, from_mp(h)?);
        assert!(l <= h, "oracle bracket reversed: {l:?} > {h:?}");
        Some((l, h))
    };
    match (part(lo.real(), hi.real()), part(lo.imag(), hi.imag())) {
        (Some(re), Some(im)) => Answer::Bracket([re, im]),
        _ => Answer::NonFinite,
    }
}

/// `[v - eps, v + eps]` in each part, rounded outward.
fn around(v: &MpComplex, eps: &MpFloat) -> Answer {
    let part = |x: &MpFloat| {
        let lo = MpFloat::with_val_round(x.prec(), x - eps, Round::Down).0;
        let hi = MpFloat::with_val_round(x.prec(), x + eps, Round::Up).0;
        Some((from_mp(&lo)?, from_mp(&hi)?))
    };
    match (part(v.real()), part(v.imag())) {
        (Some(re), Some(im)) => Answer::Bracket([re, im]),
        _ => Answer::NonFinite,
    }
}

/// The exponent of the smaller non-zero part, or 0.
fn small_exp(z: &MpComplex) -> i32 {
    [z.real(), z.imag()]
        .into_iter()
        .filter_map(MpFloat::get_exp)
        .min()
        .unwrap_or(0)
}

/// `exp(z) - 1`: `Re e^z - 1` cancels about twice the bits of the smaller
/// part below the point (`cos b - 1 ≈ -b^2/2`), so `exp` runs that much
/// higher; each step rounds in the bracket's direction.
fn expm1_oracle(z: &MpComplex, prec: u64) -> Answer {
    let below = u64::from(small_exp(z).min(0).unsigned_abs());
    mpc_bracket(prec + 2 * below + 64, |p, r| {
        let e = MpComplex::with_val_round(p, z.exp_ref(), r).0;
        let re = MpFloat::with_val_round(p, e.real() - 1u32, r.0).0;
        MpComplex::with_val(p, (re, e.imag()))
    })
}

/// `log z` from MPFR, each part rounded in its direction: `log|z| =
/// log1p(|z|^2 - 1) / 2` with `|z|^2 - 1` exact, and `arg z = atan2(b, a)`.
/// MPC's `log` misses `log|z|` next to `|z| = 1`: at 400 bits it put the
/// real part of `log(1 + z)` for `z ≈ 2^-280 - 2^-217 i` 2^-134 off,
/// relative (CAS-113; the complex tests saw it at 66 bits).
fn mpfr_log(z: &MpComplex, p: u32, r: (Round, Round)) -> MpComplex {
    let (a, b) = (z.real(), z.imag());
    // |z|^2 - 1 exactly: the bits from the top of the largest term to the
    // lowest bit of the squares.
    let (mut top, mut low) = (1, 0);
    for x in [a, b] {
        if let Some(e) = x.get_exp() {
            top = top.max(2 * e);
            low = low.min(2 * (e - i32::try_from(x.prec()).expect("small")));
        }
    }
    let q = u32::try_from(top - low + 4).expect("positive");
    let n = MpFloat::with_val(q, a.square_ref()) + MpFloat::with_val(q, b.square_ref());
    let n = MpFloat::with_val(q, n - 1u32);
    let re = MpFloat::with_val_round(p, n.ln_1p_ref(), r.0).0 >> 1u32;
    let im = MpFloat::with_val_round(p, b.atan2_ref(a), r.1).0;
    MpComplex::with_val(p, (re, im))
}

/// `log(1 + z)` with `1 + z` exact.
fn log1p_oracle(z: &MpComplex, prec: u64) -> Answer {
    let re = z.real();
    let bits = re.get_exp().map_or(2, |e| {
        let low = e - i32::try_from(re.prec()).expect("small");
        u32::try_from(e.max(1) - low.min(0) + 2).expect("positive")
    });
    let w = MpComplex::with_val((bits, z.imag().prec()), z + 1u32);
    mpc_bracket(prec, |p, r| mpfr_log(&w, p, r))
}

/// `|x|` rounded up, at 64 bits.
fn abs_up(x: &MpComplex) -> MpFloat {
    MpFloat::with_val_round(64, x.abs_ref(), Round::Up).0
}

/// `exp(t)` for `t = w log z` (`w` given, or `1/n` when `w` is `None`) at
/// `prec + 300` bits, with [`mpfr_log`] rounded to nearest. Rounding
/// `log z` and `t` moves `t` by at most `(|w| |log z| + |t|) 2^(2-p)`,
/// which `exp` turns into `|v|` times as much, plus its own rounding of
/// `|v| 2^(1-p)`: the margin is `|v| (|w| |log z| + |t| + 1) 2^-(p-4)`.
fn exp_log_oracle(z: &MpComplex, w: Option<&MpComplex>, n: u32, prec: u64) -> Answer {
    let p = u32::try_from(prec + 300).expect("small precision");
    let l = mpfr_log(z, p, (Round::Nearest, Round::Nearest));
    let (t, wl) = match w {
        Some(w) => (MpComplex::with_val(p, w * &l), abs_up(w) * abs_up(&l)),
        None => (MpComplex::with_val(p, &l / n), abs_up(&l)),
    };
    let v = MpComplex::with_val(p, t.exp_ref());
    let grow = MpFloat::with_val_round(64, wl + abs_up(&t), Round::Up).0 + 1u32;
    around(
        &v,
        &(MpFloat::with_val_round(64, abs_up(&v) * grow, Round::Up).0 >> (p - 4)),
    )
}

/// `sin(π z)` or `cos(π z)` at `prec + 300` bits: the rounded `π z` is
/// within `2 |π z| 2^-p`, and `|sin'|, |cos'| <= cosh(π b)`, so
/// `(1 + |z|) cosh(π b) 2^-(p - 8)` covers the error (as the complex tests).
fn pi_oracle(z: &MpComplex, prec: u64, cos: bool) -> Answer {
    let p = u32::try_from(prec + 300).expect("small precision");
    let pi = MpFloat::with_val(p, Constant::Pi);
    let arg = MpComplex::with_val(p, z * &pi);
    let v = if cos {
        MpComplex::with_val(p, arg.cos_ref())
    } else {
        MpComplex::with_val(p, arg.sin_ref())
    };
    let size = MpFloat::with_val(64, z.abs_ref()) + 1u32;
    let grow = MpFloat::with_val(64, arg.imag().clone().abs().cosh());
    around(&v, &((size * grow) >> (p - 8)))
}

/// The most bits of an exact power the oracle computes.
const MAX_EXACT_BITS: u64 = 1 << 22;

/// `x 2^s` as a rational.
fn scaled(x: rug::Integer, s: i64) -> rug::Rational {
    let shift = u32::try_from(s.unsigned_abs()).expect("bounded");
    if s >= 0 {
        rug::Rational::from(x << shift)
    } else {
        rug::Rational::from((x, rug::Integer::from(1) << shift))
    }
}

/// `z^n` for the exact point `z` by Gaussian integers (CAS-362): every power
/// of a dyadic point is an exact dyadic, which MPC's `pow` rounds by raising
/// its precision, 71 s for `(−16 − 16i)^29435` at 9940 bits. `Unknown` past
/// [`MAX_EXACT_BITS`].
fn power_oracle(zm: &MpComplex, n: i64, prec: u64) -> Answer {
    let parts = |x: &MpFloat| -> Option<(rug::Integer, i64)> {
        if x.is_zero() {
            return Some((rug::Integer::new(), 0));
        }
        let (m, e) = x.to_integer_exp()?;
        Some((m, i64::from(e)))
    };
    let (Some((a, ea)), Some((b, eb))) = (parts(zm.real()), parts(zm.imag())) else {
        return Answer::Unknown;
    };
    let k = n.unsigned_abs();
    if a == 0 && b == 0 {
        return match n.cmp(&0) {
            std::cmp::Ordering::Less => Answer::NonFinite,
            std::cmp::Ordering::Equal => {
                around(&MpComplex::with_val(64, (1, 0)), &MpFloat::new(64))
            }
            std::cmp::Ordering::Greater => around(&MpComplex::new(64), &MpFloat::new(64)),
        };
    }
    // a + b i over the common exponent e.
    let e = if a == 0 {
        eb
    } else if b == 0 {
        ea
    } else {
        ea.min(eb)
    };
    let lift = |m: rug::Integer, at: i64| m << u32::try_from(at - e).expect("small");
    let (a, b) = (lift(a, ea), lift(b, eb));
    let bits = u64::from(a.significant_bits().max(b.significant_bits())) + 1;
    if bits.saturating_mul(k) > MAX_EXACT_BITS
        || e.unsigned_abs().saturating_mul(k) > MAX_EXACT_BITS
    {
        return Answer::Unknown;
    }
    // (a + b i)^k by squaring.
    let (mut x, mut y) = (rug::Integer::from(1), rug::Integer::new());
    let (mut p, mut q) = (a, b);
    let mut rest = k;
    while rest > 0 {
        if rest & 1 == 1 {
            let nx = rug::Integer::from(&x * &p) - rug::Integer::from(&y * &q);
            let ny = rug::Integer::from(&x * &q) + rug::Integer::from(&y * &p);
            (x, y) = (nx, ny);
        }
        rest >>= 1;
        if rest > 0 {
            let np = rug::Integer::from(&p * &p) - rug::Integer::from(&q * &q);
            let nq = rug::Integer::from(&p * &q) << 1u32;
            (p, q) = (np, nq);
        }
    }
    let shift = e * i64::try_from(k).expect("small");
    let (re, im) = if n >= 0 {
        (scaled(x, shift), scaled(y, shift))
    } else {
        // 1/((x + y i) 2^s) = (x − y i)/((x² + y²) 2^s).
        let d = rug::Integer::from(&x * &x) + rug::Integer::from(&y * &y);
        let re = rug::Rational::from((x, d.clone())) * scaled(rug::Integer::from(1), -shift);
        let im = rug::Rational::from((-y, d)) * scaled(rug::Integer::from(1), -shift);
        (re, im)
    };
    let p = u32::try_from(prec + 192).expect("small precision");
    let round = |v: &rug::Rational, r: Round| from_mp(&MpFloat::with_val_round(p, v, r).0);
    match (
        round(&re, Round::Down),
        round(&re, Round::Up),
        round(&im, Round::Down),
        round(&im, Round::Up),
    ) {
        (Some(rl), Some(rh), Some(il), Some(ih)) => Answer::Bracket([(rl, rh), (il, ih)]),
        _ => Answer::Unknown,
    }
}

/// The principal `z^(1/n) = exp(log(z) / n)`, and 0 at 0.
fn root_oracle(z: &MpComplex, prec: u64, n: u32) -> Answer {
    if z.real().is_zero() && z.imag().is_zero() {
        let zero = MpComplex::new(64);
        return around(&zero, &MpFloat::new(64));
    }
    exp_log_oracle(z, None, n, prec)
}

/// The principal `z^w`: MPC's special cases at `z = 0`, which take no
/// logarithm, and [`exp_log_oracle`] elsewhere.
fn pow_oracle(z: &MpComplex, w: &MpComplex, prec: u64) -> Answer {
    if z.real().is_zero() && z.imag().is_zero() {
        return mpc_bracket(prec, |p, r| MpComplex::with_val_round(p, z.pow(w), r).0);
    }
    exp_log_oracle(z, Some(w), 1, prec)
}

/// Checks `w` against MPC's answer at one point.
fn check(w: &ComplexBall, answer: &Answer, what: impl Fn() -> String) {
    match answer {
        Answer::Bracket([(rl, rh), (il, ih)]) => assert!(
            meets(w.re(), rl, rh) && meets(w.im(), il, ih),
            "{}: {w:?} misses [{rl:?}, {rh:?}] + [{il:?}, {ih:?}] i",
            what()
        ),
        Answer::NonFinite => assert!(!w.is_finite(), "{}: a pole, but {w:?}", what()),
        Answer::Unknown => {}
    }
}

/// For an exact input where MPC has a value: a finite result, and when
/// `tight` and the input's parts are below `2^8` (and zero or above
/// `2^-250`), radii within `2^(6 - prec)` of the larger part.
fn check_exact(
    z: &ComplexBall,
    w: &ComplexBall,
    answer: &Answer,
    prec: u64,
    tight: bool,
    what: impl Fn() -> String,
) {
    if !matches!(answer, Answer::Bracket(_)) {
        return;
    }
    assert!(w.is_finite(), "{}: not finite", what());
    let ordinary = |b: &Ball| {
        let e = b.mid().exponent();
        b.mid().is_zero() || (*e <= Integer::from(8i64) && *e >= Integer::from(-250i64))
    };
    if !tight || !ordinary(z.re()) || !ordinary(z.im()) {
        return;
    }
    let size = Mag::from_float(w.re().mid()).max(Mag::from_float(w.im().mid()));
    if !size.is_zero() {
        let tol = size.mul_2exp(&(Integer::from(6i64) - Integer::from(prec)));
        assert!(
            *w.re().rad() <= tol && *w.im().rad() <= tol,
            "{}: loose {w:?}",
            what()
        );
    }
}

fn unary(
    u: &mut Unstructured<'_>,
    name: &str,
    z: &ComplexBall,
    w: &ComplexBall,
    prec: u64,
    tight: bool,
    oracle: &dyn Fn(&MpComplex, u64) -> Answer,
) -> Result<()> {
    for zm in cpoints(u, z, prec)? {
        let answer = oracle(&zm, prec);
        let what = || format!("{name}({z:?}) at {prec} bits, at {zm}");
        check(w, &answer, what);
        if z.is_exact() {
            check_exact(z, w, &answer, prec, tight, what);
        }
    }
    Ok(())
}

/// `abs` and `arg`, real balls, against MPFR's brackets of MPC's values.
fn real(u: &mut Unstructured<'_>, name: &str, z: &ComplexBall, prec: u64) -> Result<()> {
    let y = if name == "abs" {
        z.abs(prec)
    } else {
        z.arg(prec)
    };
    for zm in cpoints(u, z, prec)? {
        let answer = bracket(prec + 64, |p, r| {
            if name == "abs" {
                MpFloat::with_val_round(p, zm.abs_ref(), r).0
            } else {
                MpFloat::with_val_round(p, zm.arg_ref(), r).0
            }
        });
        let what = || format!("{name}({z:?}) at {prec} bits, at {zm}");
        match &answer {
            Oracle::Bracket(lo, hi) => assert!(
                meets(&y, lo, hi),
                "{}: {y:?} misses [{lo:?}, {hi:?}]",
                what()
            ),
            Oracle::Nan | Oracle::Infinite => panic!("{}: MPC has no value", what()),
        }
        if z.is_exact() {
            assert!(y.is_finite(), "{}: not finite", what());
        }
    }
    Ok(())
}

/// `mul`, `div` and `pow` at pairs of points: each point of one ball with
/// the other's centre, and a few fuzzer-chosen pairs.
fn binary(
    u: &mut Unstructured<'_>,
    name: &str,
    z: &ComplexBall,
    w: &ComplexBall,
    prec: u64,
) -> Result<()> {
    let v = match name {
        "mul" => z.mul(w, prec),
        "div" => z.div(w, prec),
        _ => z.pow(w, prec),
    };
    let (pz, pw) = (cpoints(u, z, prec)?, cpoints(u, w, prec)?);
    if pz.is_empty() || pw.is_empty() {
        return Ok(());
    }
    let mut pairs: Vec<(usize, usize)> = (0..pz.len())
        .map(|i| (i, 0))
        .chain((1..pw.len()).map(|j| (0, j)))
        .collect();
    for _ in 0..3 {
        pairs.push((u.choose_index(pz.len())?, u.choose_index(pw.len())?));
    }
    for (i, j) in pairs {
        let (a, b) = (&pz[i], &pw[j]);
        let answer = match name {
            "mul" => mpc_bracket(prec, |p, r| MpComplex::with_val_round(p, a * b, r).0),
            "div" => mpc_bracket(prec, |p, r| MpComplex::with_val_round(p, a / b, r).0),
            _ => pow_oracle(a, b, prec),
        };
        let what = || format!("{name}({z:?}, {w:?}) at {prec} bits, at {a}, {b}");
        check(&v, &answer, what);
        if z.is_exact() && w.is_exact() {
            check_exact(z, &v, &answer, prec, true, what);
        }
    }
    Ok(())
}
