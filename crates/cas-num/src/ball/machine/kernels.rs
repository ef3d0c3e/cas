//! Machine-ball elementary functions (CAS-DOC-5 §5.4, CAS-88): argument
//! reduction and a fixed polynomial in double-double ([`Dd`]), with an
//! a-priori relative error bound per kernel, so the `f64` midpoint costs
//! half an ulp plus a radius term `KERNEL_ERR |v|`. No platform libm.
//!
//! The reduction constants (π/2 and log 2 split Cody–Waite style, the table
//! `atan(k/8)`) and the Taylor coefficients are built once, at first use,
//! from the `Ball` constants and exact rationals: correct by construction,
//! no hand-typed digits. Each kernel's error budget is written next to it;
//! every one stays below `2^-94` and is claimed as [`KERNEL_ERR`] `= 2^-90`.
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
///
/// Budget, relative:
///
/// ```text
/// m = k log 2 + t, |k| <= 1022 < 2^11: k l0, k l1 exact, m - k l0 exact
///   (Sterbenz, |k| >= 1), two double-double subtractions, Cody–Waite rest
///   |t - t_dd| <= 2 OP_ERR |t| + 1022 ln2_err < 2^-100    -> e^t: 2^-99
/// Horner, degree 24, |t| <= 0.35: γ_48 e^|t| / e^t
///   <= 49 · 2^-102 · e^0.7                                  < 2^-95
/// coefficients 2^-106 e^0.7                                 < 2^-104
/// truncation 0.35^25 / 25!                                  < 2^-121
/// total                                                     < 2^-94
/// ```
///
/// The scaling by `2^k` is exact (the result and its low part stay normal).
pub(super) fn exp(m: f64) -> Option<Dd> {
    if !(-665.0..=708.0).contains(&m) {
        return None;
    }
    let c = consts();
    #[allow(clippy::cast_possible_truncation, reason = "|m / log 2| <= 1022")]
    let k = (m * std::f64::consts::LOG2_E).round() as i32;
    let kf = f64::from(k);
    let r1 = m - kf * c.ln2[0];
    let t = Dd::from_f64(r1)
        .sub(Dd::prod(kf, c.ln2[1]))
        .sub(Dd::prod(kf, c.ln2[2]));
    debug_assert!(c.ln2_err < 1e-36 && t.hi.abs() < 0.35);
    Some(horner(&c.exp_c, t).mul_pow2(k))
}

/// `log x` for a positive normal `x`.
///
/// Budget, relative:
///
/// ```text
/// x = f 2^e, f in [1/√2, √2); log f = 2 atanh(s), s = (f - 1)/(f + 1),
///   |s| <= 0.1716; f - 1 and f + 1 exact, s within OP_ERR, y = s^2
///   within 3 OP_ERR; Horner on positive terms γ_44 ≈ 45 OP_ERR,
///   coefficients 2^-106, truncation y^23 / 47 < 2^-120, the product
///   by 2s 2 OP_ERR: log f within 52 OP_ERR                  < 2^-96
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
    let s = Dd::from_f64(f - 1.0).div(Dd::sum(f, 1.0));
    let y = s.mul(s);
    let atanh = horner(&c.atanh_c, y).mul(s).mul_f64(2.0);
    Some(c.ln2_dd.mul_f64(f64::from(e)).add(atanh))
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
    let y = t.mul(t);
    let s = horner(&c.sin_c, y).mul(t);
    let co = horner(&c.cos_c, y);
    Some(match k.rem_euclid(4) {
        0 => (s, co),
        1 => (co, s.neg()),
        2 => (s.neg(), co.neg()),
        _ => (co.neg(), s),
    })
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
