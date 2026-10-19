//! Machine-ball elementary functions (CAS-DOC-5 §5.4, CAS-88): argument
//! reduction and a fixed polynomial in double-double ([`Dd`]), with an
//! a-priori relative error bound per kernel, so the `f64` midpoint costs
//! half an ulp plus a radius term `KERNEL_ERR |v|`. No platform libm.
//!
//! The reduction constants (π/2 and log 2 split Cody–Waite style
//! \[CW80\], the table `atan(k/8)`) and the Taylor coefficients are built
//! once, at first use, from the `Ball` constants and exact rationals:
//! correct by construction, no hand-typed digits. Each kernel's error budget
//! is written next to it; every one stays below `2^-94` and is claimed as
//! [`KERNEL_ERR`] `= 2^-90`. Horner's rule's errors are the `γ_n` bounds of
//! \[Hig02\].
//!
//! The CAS-103 kernels build on these pieces ([`exp_dd`] and [`log_dd`] of
//! double-double arguments, [`log1p_dd`] on `[1/√2 - 1, √2 - 1]`, the
//! `sin`/`cos` polynomials): `sinh` and `cosh` from `e^|m|` (the sine
//! series at `-m^2` below 1/2), `tanh` as their quotient, `asinh`, `acosh`
//! and `atanh` through `log` where its argument is at least 1.41 and
//! `log1p` (or the `atanh` series) below, with no cancellation in either,
//! and `sin(π m)`, `cos(π m)` after an exact reduction modulo 2. Their
//! budgets stay below `2^-92`. `pow` and `root` go through `exp(t)`, `t =
//! y log x` or `log(x) / n`, and claim `(1 + |t|) 2^-93`, as `log x`'s
//! relative error becomes `t`'s absolute one.
//!
//! Arguments outside a kernel's range return `None` from the kernel; the
//! `MachineBall` method then falls back to `Ball` at 64 bits.

use std::sync::OnceLock;

use super::dd::{Dd, horner};
use crate::ball::{Ball, Float, Round};
use crate::{Integer, Rational};

/// The relative error bound claimed for every kernel's double-double value.
pub(super) const KERNEL_ERR: f64 = f64::from_bits((1023 - 90) << 52);

/// Degrees: `exp` in `t` (24), `sin(t)/t` and `cos t` in `t^2` (15), the
/// `atanh` and `atan` series in the square (22 and 14).
const EXP_DEG: usize = 24;
const TRIG_DEG: usize = 15;
const ATANH_DEG: usize = 22;
const ATAN_DEG: usize = 14;

/// Constants built at first use.
struct Consts {
    /// log 2 = l[0] + l[1] + l[2] + δ, l[0], l[1] of 42 bits (so `k l[i]`
    /// is exact for `|k| < 2^11`), |δ| <= `ln2_err`.
    ln2: [f64; 3],
    ln2_err: f64,
    ln2_dd: Dd,
    /// π/2 = p[0] + … + p[3] + δ, p[0], p[1] of 33 bits (`k p[i]` exact for
    /// `|k| < 2^20`), |δ| <= `half_pi_err`.
    half_pi: [f64; 4],
    half_pi_err: f64,
    half_pi_dd: Dd,
    /// `1/k!` for `k <= EXP_DEG`.
    exp_c: Vec<Dd>,
    /// (-1)^j / (2j+1)! and (-1)^j / (2j)!.
    sin_c: Vec<Dd>,
    cos_c: Vec<Dd>,
    /// 1/(2j+1) and (-1)^j/(2j+1).
    atanh_c: Vec<Dd>,
    atan_c: Vec<Dd>,
    /// atan(k/8), k = 0..=8.
    atan_table: Vec<Dd>,
}

/// The double-double nearest a float (`|x - dd| <= 2^-106 |x|` roughly, as
/// both parts round to nearest), and the exact rest.
fn to_dd(x: &Float) -> (Dd, Float) {
    let hi = x.round(53, Round::Nearest).0;
    let rest = x.sub(&hi);
    let lo = rest.round(53, Round::Nearest).0;
    (
        Dd {
            hi: hi.to_f64(),
            lo: lo.to_f64(),
        },
        rest.sub(&lo),
    )
}

fn rat_dd(q: &Rational) -> Dd {
    to_dd(&Float::from_rational(q, 300, Round::Nearest).0).0
}

/// `x` split into parts of `bits[i]` bits each and the bound on what is
/// left (`|x - Σ parts| <= err`), from a ball `x` at 300 bits.
fn split(x: &Ball, bits: &[u64]) -> (Vec<f64>, f64) {
    let mut rest = x.mid().clone();
    let mut parts = Vec::with_capacity(bits.len());
    for &b in bits {
        let p = rest.round(b, Round::Nearest).0;
        rest = rest.sub(&p);
        parts.push(p.to_f64());
    }
    let left = crate::ball::Mag::from_float(&rest).add(x.rad());
    let err = left
        .to_float()
        .map_or(f64::INFINITY, |f| f.round(53, Round::Up).0.to_f64());
    (parts, err)
}

fn consts() -> &'static Consts {
    static CONSTS: OnceLock<Consts> = OnceLock::new();
    CONSTS.get_or_init(|| {
        let fact = |n: u64| (1..=n).fold(Integer::ONE, |a, k| &a * &Integer::from(k));
        let q = |n: Integer, d: Integer| Rational::new(n, d);
        let sign = |j: u64| {
            if j % 2 == 0 {
                Integer::ONE
            } else {
                -Integer::ONE
            }
        };
        let deg = |d: usize| u64::try_from(d).expect("small");
        let ln2 = Ball::const_log2(300);
        let (l, ln2_err) = split(&ln2, &[42, 42, 53]);
        let half_pi = Ball::const_pi(300).mul_2exp(&-Integer::ONE);
        let (p, half_pi_err) = split(&half_pi, &[33, 33, 53, 53]);
        let atan_table = (0..=8i64)
            .map(|k| {
                let b = Ball::from_rational(&q(Integer::from(k), Integer::from(8)), 300).atan(300);
                to_dd(b.mid()).0
            })
            .collect();
        Consts {
            ln2: [l[0], l[1], l[2]],
            ln2_err,
            ln2_dd: to_dd(ln2.mid()).0,
            half_pi: [p[0], p[1], p[2], p[3]],
            half_pi_err,
            half_pi_dd: to_dd(half_pi.mid()).0,
            exp_c: (0..=deg(EXP_DEG))
                .map(|k| rat_dd(&q(Integer::ONE, fact(k))))
                .collect(),
            sin_c: (0..=deg(TRIG_DEG))
                .map(|j| rat_dd(&q(sign(j), fact(2 * j + 1))))
                .collect(),
            cos_c: (0..=deg(TRIG_DEG))
                .map(|j| rat_dd(&q(sign(j), fact(2 * j))))
                .collect(),
            atanh_c: (0..=deg(ATANH_DEG))
                .map(|j| rat_dd(&q(Integer::ONE, Integer::from(2 * j + 1))))
                .collect(),
            atan_c: (0..=deg(ATAN_DEG))
                .map(|j| rat_dd(&q(sign(j), Integer::from(2 * j + 1))))
                .collect(),
            atan_table,
        }
    })
}

/// `e^m` for `-665 <= m <= 708`: below, the scaled low part `lo 2^k`
/// (about `2^-53 e^m`) would fall into the subnormal range and lose its
/// bits (`e^-690` came out `2^-79` off), so those go to `Ball`.
pub(super) fn exp(m: f64) -> Option<Dd> {
    exp_dd(Dd::from_f64(m))
}

/// `e^m` for a double-double `m` with `-665 <= m.hi <= 708` (see [`exp`]).
///
/// Budget, relative:
///
/// ```text
/// m = k log 2 + t, |k| <= 1022 < 2^11: k l0, k l1 exact, m.hi - k l0
///   exact (Sterbenz, |k| >= 1), three double-double sums (with m.lo, at
///   most 2^-44), Cody–Waite rest
///   |t - t_dd| <= 3 OP_ERR |t| + 1022 ln2_err < 2^-100    -> e^t: 2^-99
/// Horner, degree 24, |t| <= 0.35: γ_48 e^|t| / e^t
///   <= 49 · 2^-102 · e^0.7                                  < 2^-95
/// coefficients 2^-106 e^0.7                                 < 2^-104
/// truncation 0.35^25 / 25!                                  < 2^-121
/// total                                                     < 2^-94
/// ```
///
/// The scaling by `2^k` is exact (the result and its low part stay normal).
pub(super) fn exp_dd(m: Dd) -> Option<Dd> {
    if !(-665.0..=708.0).contains(&m.hi) {
        return None;
    }
    let c = consts();
    #[allow(clippy::cast_possible_truncation, reason = "|m / log 2| <= 1022")]
    let k = (m.hi * std::f64::consts::LOG2_E).round() as i32;
    let kf = f64::from(k);
    let r1 = m.hi - kf * c.ln2[0];
    let mut t = Dd::from_f64(r1)
        .sub(Dd::prod(kf, c.ln2[1]))
        .sub(Dd::prod(kf, c.ln2[2]));
    if m.lo != 0.0 {
        t = t.add(Dd::from_f64(m.lo));
    }
    debug_assert!(c.ln2_err < 1e-36 && t.hi.abs() < 0.35);
    Some(horner(&c.exp_c, t).mul_pow2(k))
}

/// `log x` for a positive normal `x`.
///
/// Budget, relative:
///
/// ```text
/// x = f 2^e, f in [1/√2, √2); log f = log1p(f - 1), f - 1 exact:
///   within 56 OP_ERR (log1p_dd)                             < 2^-96
/// e log 2: |e| 2^-106 log 2 absolute, against |log x| >= |e| log 2 - 0.35
///   (or e = 0); the final sum OP_ERR
/// total                                                     < 2^-94
/// ```
pub(super) fn log(x: f64) -> Option<Dd> {
    if !(x.is_normal() && x > 0.0) {
        return None;
    }
    let c = consts();
    let bits = x.to_bits();
    let raw = i32::try_from((bits >> 52) & 0x7ff).expect("11 bits") - 1023;
    // f = x 2^-raw in [1, 2); fold into [1/√2, √2).
    let mut e = raw;
    let mut f = f64::from_bits((bits & ((1 << 52) - 1)) | (1023 << 52));
    if f > std::f64::consts::SQRT_2 {
        f /= 2.0;
        e += 1;
    }
    let lf = log1p_dd(Dd::from_f64(f - 1.0));
    Some(c.ln2_dd.mul_f64(f64::from(e)).add(lf))
}

/// `log(1 + u)` for a double-double `u` in `[1/√2 - 1, √2 - 1]`, as `2
/// atanh(s)`, `s = u / (2 + u)` in `[-0.1716, 0.1716]`.
///
/// Budget, relative:
///
/// ```text
/// 2 + u >= 1.7 within OP_ERR (exact for a double u), s within 2 OP_ERR,
///   y = s^2 within 5 OP_ERR; Horner on positive terms γ_44 ≈ 45 OP_ERR
///   (y's error moves the series by at most y/3 of it), coefficients
///   2^-106, truncation y^23 / 47 < 2^-120; the product by 2s 2 OP_ERR
///   and s's own 2 OP_ERR: within 56 OP_ERR                   < 2^-96
/// ```
pub(super) fn log1p_dd(u: Dd) -> Dd {
    let c = consts();
    let mut den = Dd::sum(2.0, u.hi);
    if u.lo != 0.0 {
        den = den.add(Dd::from_f64(u.lo));
    }
    let s = u.div(den);
    let y = s.mul(s);
    horner(&c.atanh_c, y).mul(s).mul_f64(2.0)
}

/// `log v` for a double-double `v` with a positive normal `v.hi`, as
/// `log(v.hi) + δ`, `δ = v.lo / v.hi`.
///
/// Budget, absolute: `log(v.hi)` within `2^-94 |log v.hi|` ([`log`]);
/// `log1p(δ) = δ` within `δ^2 / 2 <= 2^-107` for `|δ| <= 2^-53`, and `δ`
/// rounds by `2^-106`; the sum `OP_ERR` relative. Against `|log v| >= 1/3`
/// (every caller) that is below `2^-94 + 2^-103 + 2^-102 < 2^-93.9`
/// relative.
pub(super) fn log_dd(v: Dd) -> Option<Dd> {
    Some(log(v.hi)?.add(Dd::from_f64(v.lo / v.hi)))
}

/// `(sin m, cos m)` for `|m| <= 2^20`, `None` when `m` is within `2^-40` of
/// a non-zero multiple of π/2 (relative accuracy of the reduced argument).
///
/// Budget, relative to each result:
///
/// ```text
/// m = k π/2 + t, |k| < 2^20: Cody–Waite rest 2^20 half_pi_err < 2^-150,
///   three subtractions 3 OP_ERR |t|; against |t| >= 2^-40 t is within 2^-100
/// Horner, degree 15, y = t^2 <= 0.62: the alternating series against its
///   value (sinh|t| / sin|t|, cosh t / cos t <= 1.9 for |t| <= π/4 + 2^-30):
///   1.9 γ_30                                                < 2^-96
/// coefficients 2^-106 · 1.9; truncation y^16 / 31!          < 2^-120
/// t × S                                                     OP_ERR
/// total                                                     < 2^-94
/// ```
pub(super) fn sin_cos(m: f64) -> Option<(Dd, Dd)> {
    // m is finite (a machine-ball midpoint).
    if m.abs() > 1_048_576.0 {
        return None;
    }
    // Below 2^-300, sin m = m and cos m = 1 within m^2/6 < 2^-600 relative
    // (and the squares below would underflow).
    if m.abs() < f64::from_bits((1023 - 300) << 52) {
        return Some((Dd::from_f64(m), Dd::ONE));
    }
    let c = consts();
    #[allow(clippy::cast_possible_truncation, reason = "|m| <= 2^20")]
    let k = (m * std::f64::consts::FRAC_2_PI).round() as i64;
    #[allow(clippy::cast_precision_loss, reason = "|k| < 2^20")]
    let kf = k as f64;
    let t = if k == 0 {
        Dd::from_f64(m)
    } else {
        let r1 = m - kf * c.half_pi[0];
        Dd::from_f64(r1)
            .sub(Dd::prod(kf, c.half_pi[1]))
            .sub(Dd::prod(kf, c.half_pi[2]))
            .sub(Dd::prod(kf, c.half_pi[3]))
    };
    debug_assert!(c.half_pi_err < 1e-50);
    if t.hi.abs() < f64::from_bits((1023 - 40) << 52) && k != 0 {
        return None;
    }
    if t.hi == 0.0 {
        return Some((Dd::ZERO, Dd::ONE));
    }
    let (s, co) = sin_cos_poly(c, t);
    Some(match k.rem_euclid(4) {
        0 => (s, co),
        1 => (co, s.neg()),
        2 => (s.neg(), co.neg()),
        _ => (co.neg(), s),
    })
}

/// `(sin t, cos t)` for a double-double `|t| <= π/4 + 2^-30`, `|t| >=
/// 2^-300`: the Horner part of [`sin_cos`]'s budget, within `2^-96` of each
/// at the given `t`.
#[inline]
fn sin_cos_poly(c: &Consts, t: Dd) -> (Dd, Dd) {
    let y = t.mul(t);
    (horner(&c.sin_c, y).mul(t), horner(&c.cos_c, y))
}

/// `atan x` for `2^-500 <= |x| <= 2^500` (or 0).
///
/// Budget, relative:
///
/// ```text
/// |x| > 1 folds to π/2 - atan(1/x): 1/x within OP_ERR, the difference
///   at least π/4, no cancellation
/// c = k/8 nearest y, y' = (y - c)/(1 + y c), |y'| <= 1/16: within 4 OP_ERR
///   of exact relative to y (against a result >= atan(1/16) for k >= 1,
///   or relative to y' = y for k = 0)
/// Horner, degree 14 in y'^2 <= 2^-8, ratio below 1.01: γ_28  < 2^-97
/// coefficients 2^-106, truncation 2^-120, the table 2^-106, two sums
/// total                                                     < 2^-94
/// ```
pub(super) fn atan(x: f64) -> Option<Dd> {
    let a = x.abs();
    if x == 0.0 {
        return Some(Dd::ZERO);
    }
    let lim = f64::from_bits((1023 + 500) << 52);
    if !(a >= 1.0 / lim && a <= lim) {
        return None;
    }
    let c = consts();
    let (y, fold) = if a > 1.0 {
        (Dd::ONE.div(Dd::from_f64(a)), true)
    } else {
        (Dd::from_f64(a), false)
    };
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "0 <= 8y <= 8"
    )]
    let k = (y.hi * 8.0).round() as usize;
    #[allow(clippy::cast_precision_loss, reason = "k <= 8")]
    let cf = k as f64 / 8.0;
    let yr = if k == 0 {
        y
    } else {
        y.sub(Dd::from_f64(cf)).div(y.mul_f64(cf).add(Dd::ONE))
    };
    let at = horner(&c.atan_c, yr.mul(yr)).mul(yr);
    let mut v = c.atan_table[k].add(at);
    if fold {
        v = c.half_pi_dd.sub(v);
    }
    Some(if x < 0.0 { v.neg() } else { v })
}

/// Below this, `sinh m`, `asinh m` and `atanh m` are `m` within `m^2 / 3 <=
/// 2^-600` relative, and `cosh m` is 1 within `m^2` (where `m^2` would
/// leave the normal range of the double-doubles below).
const TINY: f64 = f64::from_bits((1023 - 300) << 52);

/// Above this, `asinh m = log(2m)` and `acosh m = log(2m)` within `1/(4
/// m^2) < 2^-1000` absolute (where `m^2` would overflow).
const LARGE: f64 = f64::from_bits((1023 + 500) << 52);

/// `(sinh m, cosh m)` for `|m| <= 708`.
///
/// Budget, relative to each:
///
/// ```text
/// |m| < 2^-300: m and 1
/// E = e^|m| within 2^-94 (exp)
/// |m| >= 40: e^-|m| < 2^-57, both E/2 within 2^-113 more        < 2^-93.9
/// cosh = (E + 1/E)/2: 1/E within 2^-94 + OP_ERR, the sum of
///   positives OP_ERR more                                      < 2^-93.9
/// 1/2 <= |m| < 40: sinh = (E - 1/E)/2, the difference magnifying
///   by coth(1/2) < 2.17                                        < 2^-92.8
/// |m| < 1/2: sinh = m Σ m^2j/(2j+1)!, the sine series at -m^2
///   (positive terms): m^2 exact, Horner γ_30 < 2^-97, truncation
///   (1/4)^16/33! < 2^-150, coefficients 2^-106, the product by m
///   OP_ERR                                                     < 2^-96
/// ```
pub(super) fn sinh_cosh(m: f64) -> Option<(Dd, Dd)> {
    let a = m.abs();
    if a > 708.0 {
        return None;
    }
    if a < TINY {
        return Some((Dd::from_f64(m), Dd::ONE));
    }
    let e = exp(a)?;
    // Halving by scaling: a two-product with 0.5 would overflow in
    // Veltkamp's split for e above 2^996 (|m| from about 690).
    let (s, c) = if a >= 40.0 {
        let h = e.mul_pow2(-1);
        (h, h)
    } else {
        let inv = Dd::ONE.div(e);
        let s = if a >= 0.5 {
            e.sub(inv).mul_pow2(-1)
        } else {
            horner(&consts().sin_c, Dd::prod(a, a).neg()).mul_f64(a)
        };
        (s, e.add(inv).mul_pow2(-1))
    };
    Some((if m < 0.0 { s.neg() } else { s }, c))
}

/// `tanh m` for every finite `m`.
///
/// Budget, relative: from `|m| >= 40`, `±1` within `2 e^-80 < 2^-113`;
/// below, `sinh / cosh` from [`sinh_cosh`], `2^-92.8 + 2^-93.9 + OP_ERR <
/// 2^-92`.
pub(super) fn tanh(m: f64) -> Option<Dd> {
    if m.abs() >= 40.0 {
        return m.is_finite().then(|| Dd::from_f64(1f64.copysign(m)));
    }
    let (s, c) = sinh_cosh(m)?;
    Some(s.div(c))
}

/// `asinh m` for every finite `m`.
///
/// Budget, relative, with `a = |m|`:
///
/// ```text
/// a < 2^-300: m
/// a > 2^500: log a + log 2, the sum of positives OP_ERR        < 2^-93.9
/// 0.35 <= a <= 2^500: log v, v = a + √(a^2 + 1) >= 1.409: a^2
///   exact, v within 3 OP_ERR, which moves log v by as much
///   absolute; log_dd against log v >= 0.343                    < 2^-93.8
/// a < 0.35: log1p(u), u = a + a^2 / (1 + √(1 + a^2)) <= 0.41,
///   positive terms within 5 OP_ERR, which log1p passes on at
///   most as much (u / ((1 + u) log1p u) <= 1); log1p_dd 56 OP_ERR < 2^-96
/// ```
pub(super) fn asinh(m: f64) -> Option<Dd> {
    let a = m.abs();
    if !m.is_finite() {
        return None;
    }
    if a < TINY {
        return Some(Dd::from_f64(m));
    }
    let v = if a > LARGE {
        log(a)?.add(consts().ln2_dd)
    } else {
        let a2 = Dd::prod(a, a);
        let root = a2.add(Dd::ONE).sqrt();
        if a >= 0.35 {
            log_dd(root.add(Dd::from_f64(a)))?
        } else {
            log1p_dd(a2.div(root.add(Dd::ONE)).add(Dd::from_f64(a)))
        }
    };
    Some(if m < 0.0 { v.neg() } else { v })
}

/// `acosh m` for `m >= 1`.
///
/// Budget, relative, with `t = m - 1`:
///
/// ```text
/// m = 1: 0
/// m > 2^500: log m + log 2                                     < 2^-93.9
/// 1.06 < m <= 2^500: log v, v = m + √((m - 1)(m + 1)) >= 1.41,
///   m ± 1 exact, v within 3 OP_ERR; log_dd against log v >= 0.344 < 2^-93.8
/// m <= 1.06: log1p(u), u = t + √(t (t + 2)) <= 0.412, t exact
///   (Sterbenz), u within 3 OP_ERR; log1p_dd 56 OP_ERR            < 2^-96
/// ```
pub(super) fn acosh(m: f64) -> Option<Dd> {
    if !(m >= 1.0 && m.is_finite()) {
        return None;
    }
    if m == 1.0 {
        return Some(Dd::ZERO);
    }
    if m > LARGE {
        return Some(log(m)?.add(consts().ln2_dd));
    }
    let t = Dd::sum(m, -1.0);
    if m <= 1.06 {
        let u = t.mul(t.add(Dd::from_f64(2.0))).sqrt().add(t);
        Some(log1p_dd(u))
    } else {
        log_dd(t.mul(Dd::sum(m, 1.0)).sqrt().add(Dd::from_f64(m)))
    }
}

/// `atanh m` for `|m| < 1`.
///
/// Budget, relative, with `a = |m|`:
///
/// ```text
/// a < 2^-300: m
/// a <= 0.1716: m Σ m^2j / (2j+1), log1p_dd's series at s = m: m^2
///   exact, Horner γ_44 ≈ 45 OP_ERR, truncation < 2^-120, the product
///   OP_ERR                                                     < 2^-96
/// a > 0.1716: log(v) / 2, v = (1 + a) / (1 - a) >= 1.414: 1 ± a
///   exact, v within OP_ERR; log_dd against log v >= 0.346       < 2^-93.8
/// ```
pub(super) fn atanh(m: f64) -> Option<Dd> {
    let a = m.abs();
    if a.is_nan() || a >= 1.0 {
        return None;
    }
    if a < TINY {
        return Some(Dd::from_f64(m));
    }
    let v = if a <= 0.1716 {
        horner(&consts().atanh_c, Dd::prod(a, a)).mul_f64(a)
    } else {
        log_dd(Dd::sum(1.0, a).div(Dd::sum(1.0, -a)))?.mul_pow2(-1)
    };
    Some(if m < 0.0 { v.neg() } else { v })
}

/// `(sin(π m), cos(π m), exact)` for every finite `m`: `exact` at integers
/// and half-integers, whose values (0 and ±1) are exact; `None` when the
/// reduced argument is below `2^-900` but not 0 (its double-double would
/// leave the normal range).
///
/// Budget, relative to each result:
///
/// ```text
/// r = m mod 2 in [-1, 1), the folds onto b in [0, 1/2] and onto
///   1/2 - b for b > 1/4: exact (fmod, Sterbenz)
/// t = π b' (b' <= 1/4): π/2 within 2^-106, the product by 2b' OP_ERR:
///   t within 2^-101.9, which moves sin t by as much relative and cos t
///   by at most tan(t) t <= 0.79 times as much                  < 2^-101
/// the polynomials (sin_cos_poly)                               < 2^-96
/// total                                                       < 2^-95
/// ```
pub(super) fn sin_cos_pi(m: f64) -> Option<(Dd, Dd, bool)> {
    if !m.is_finite() {
        return None;
    }
    let mut r = m % 2.0;
    if r >= 1.0 {
        r -= 2.0;
    } else if r < -1.0 {
        r += 2.0;
    }
    // sin(π r) = sgn(r) sin(π |r|), cos(π r) = cos(π |r|); a > 1/2 folds
    // onto 1 - a, which keeps the sine and negates the cosine.
    let mut a = r.abs();
    let neg_c = a > 0.5;
    if neg_c {
        a = 1.0 - a;
    }
    let (s, c, exact) = if a == 0.0 {
        (Dd::ZERO, Dd::ONE, true)
    } else if a == 0.5 {
        (Dd::ONE, Dd::ZERO, true)
    } else {
        // sin(π a) = cos(π (1/2 - a)) and the other way.
        let (b, swap) = if a <= 0.25 {
            (a, false)
        } else {
            (0.5 - a, true)
        };
        if b < f64::from_bits((1023 - 900) << 52) {
            return None;
        }
        let k = consts();
        let (s, c) = sin_cos_poly(k, k.half_pi_dd.mul_f64(2.0 * b));
        if swap { (c, s, false) } else { (s, c, false) }
    };
    Some((
        if r < 0.0 { s.neg() } else { s },
        if neg_c { c.neg() } else { c },
        exact,
    ))
}

/// `(x^y, err, log x)` for a positive normal `x`: `exp(t)`, `t = y log x`
/// as a double-double, within `err = (1 + |t|) 2^-93` relative, for
/// `-665 <= t <= 708`.
///
/// Budget, relative: `log x` within `2^-94` and the product by `y`
/// `OP_ERR`, so `t` is within `|t| (2^-94 + 2^-102)` absolute, which `exp`
/// turns into as much relative; `exp_dd`'s own `2^-94`: in all below `(1 +
/// 1.01 |t|) 2^-94`.
pub(super) fn pow(x: f64, y: f64) -> Option<(Dd, f64, Dd)> {
    let l = log(x)?;
    let t = l.mul_f64(y);
    let v = exp_dd(t)?;
    Some((v, (1.0 + t.hi.abs()) * f64::from_bits((1023 - 93) << 52), l))
}

/// `(x^(1/n), err)` for a positive normal `x` and `3 <= n < 2^53`:
/// `exp(log(x) / n)` within `err = (1 + |t|) 2^-93`, `t = log(x) / n`, as
/// for [`pow`] (the division is one more `OP_ERR`).
pub(super) fn root(x: f64, n: u64) -> Option<(Dd, f64)> {
    if !(3..1 << 53).contains(&n) {
        return None;
    }
    #[allow(clippy::cast_precision_loss, reason = "n < 2^53")]
    let t = log(x)?.div(Dd::from_f64(n as f64));
    Some((
        exp_dd(t)?,
        (1.0 + t.hi.abs()) * f64::from_bits((1023 - 93) << 52),
    ))
}

#[cfg(test)]
pub(super) fn reduction_errors() -> (f64, f64) {
    let c = consts();
    (c.ln2_err, c.half_pi_err)
}

#[cfg(test)]
pub(super) fn constants_dd() -> (Dd, Dd) {
    let c = consts();
    (c.ln2_dd, c.half_pi_dd)
}
