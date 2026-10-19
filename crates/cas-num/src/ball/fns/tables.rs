//! Table-based argument reduction for medium precision (CAS-DOC-5 §7.2,
//! after Johansson, "Efficient implementation of elementary functions in
//! the medium-precision range", ARITH 2015, \[Joh15\]).
//!
//! Two lookups shrink an argument from below 1 to below `2^-16`, so the
//! baseline series that follows needs no halvings (at the precisions the
//! tables serve, the baselines' `r ≈ p^(1/3)` is at most 17) and fewer
//! terms. With `a = trunc(|t| 2^16) = a1 2^8 + a2`, `0 <= a2 < 2^8`:
//!
//! ```text
//! exp:      e^t = e^(a1/2^8) e^(a2/2^16) e^x,          x = t - a / 2^16
//! sin, cos: angle addition of a1/2^8, a2/2^16 and x
//! log:      log(1 + u) = log(1 + a1/2^8) + log(1 + a2/2^16) + log(1 + u2),
//!           u1 = (2^8 u - a1) / (2^8 + a1),  u2 = (2^16 u1 - a2) / (2^16 + a2)
//! atan:     atan x = atan(a1/2^8) + atan(a2/2^16) + atan(x2),
//!           x1 = (2^8 x - a1) / (2^8 + a1 x),  x2 = (2^16 x1 - a2) / (2^16 + a2 x1)
//! ```
//!
//! (for `log` and `atan`, `a2` is read off `u1` and `x1`). The reductions
//! run on balls and the entries are balls, so there is nothing to prove
//! beyond the identities: every error is in a radius.
//!
//! The tables are built at first use by the baseline code: `exp` and
//! `sin`/`cos` entries as powers (rotations) of the baseline's value at one
//! step, `log` entries as `2 atanh(c / (2 + c))` and `atan` entries by the
//! baseline `atan`. A table records its precision; a request above it
//! rebuilds at `max(p, 2 · held)`, as the constants do (CAS-78), up to
//! [`TABLE_PREC_MAX`]. The step sizes `2^-8`, `2^-16` are untuned
//! placeholders (CAS-86); at 4600 bits the four tables take about 900 KB.

use std::sync::{Arc, PoisonError, RwLock};

use super::atan::{atan_halving, atan_series, atanh_series};
use super::exp::{exp_small, to_fixed};
use super::trig::sin_cos_small;
use super::{TABLE_MAX_PREC, one};
use crate::Integer;
use crate::ball::{Ball, Float};

/// The coarse step `2^-R1` and the fine step `2^-(R1 + R2)`.
const R1: i64 = 8;
const R2: i64 = 8;

/// The highest precision a table is built at: the working precision of
/// the fast paths plus the guard bits they ask for.
const TABLE_PREC_MAX: u64 = TABLE_MAX_PREC + 64;

/// The lowest one, so that small precisions do not rebuild in small steps.
const TABLE_PREC_MIN: u64 = 128;

/// Entries `coarse[i]` and `fine[j]`, valid at `prec` bits.
pub(super) struct Table<T> {
    prec: u64,
    coarse: Vec<T>,
    fine: Vec<T>,
}

/// A table built at first use.
struct Lazy<T>(RwLock<Option<Arc<Table<T>>>>);

impl<T> Lazy<T> {
    const fn new() -> Lazy<T> {
        Lazy(RwLock::new(None))
    }

    /// The table at `prec` bits or more, built by `build` at the target
    /// precision when the one held is short.
    fn get(&self, prec: u64, build: fn(u64) -> Table<T>) -> Arc<Table<T>> {
        let have = {
            let c = self.0.read().unwrap_or_else(PoisonError::into_inner);
            match &*c {
                Some(t) if t.prec >= prec => return Arc::clone(t),
                Some(t) => t.prec,
                None => 0,
            }
        };
        // Built without the lock held, as for the constants.
        let target = prec
            .max(have.saturating_mul(2))
            .clamp(TABLE_PREC_MIN, TABLE_PREC_MAX)
            .max(prec);
        let t = Arc::new(build(target));
        debug_assert_eq!(t.prec, target);
        let mut c = self.0.write().unwrap_or_else(PoisonError::into_inner);
        if c.as_ref().is_none_or(|c| c.prec < target) {
            *c = Some(Arc::clone(&t));
        }
        t
    }

    fn clear(&self) {
        *self.0.write().unwrap_or_else(PoisonError::into_inner) = None;
    }
}

static EXP: Lazy<Ball> = Lazy::new();
static SIN_COS: Lazy<(Ball, Ball)> = Lazy::new();
static LOG: Lazy<Ball> = Lazy::new();
static ATAN: Lazy<Ball> = Lazy::new();

/// Frees the tables; the next use rebuilds them.
pub(in crate::ball) fn clear() {
    EXP.clear();
    SIN_COS.clear();
    LOG.clear();
    ATAN.clear();
}

/// `2^k` as an exact float.
fn pow2(k: i64) -> Float {
    Float::from(1i64).mul_2exp(&Integer::from(k))
}

/// `i / 2^k` as an exact float.
fn at(i: u64, k: i64) -> Float {
    Float::from_integer_exp(&Integer::from(i), &Integer::from(-k))
}

/// `trunc(|x| 2^k)` as a `u64`, for `|x| 2^k < 2^63`.
fn digits(x: &Float, k: i64) -> u64 {
    to_fixed(&x.abs(), k.unsigned_abs())
        .to_i64()
        .and_then(|d| u64::try_from(d).ok())
        .expect("a small argument")
}

/// `x^0, x^1, …, x^n` as balls at `prec` bits, each from the last: the
/// relative error grows by about `2^-prec` a step.
fn powers<T: Clone>(first: T, step: &T, n: usize, mul: impl Fn(&T, &T) -> T) -> Vec<T> {
    let mut v = Vec::with_capacity(n + 1);
    v.push(first);
    for _ in 0..n {
        let next = mul(v.last().expect("non-empty"), step);
        v.push(next);
    }
    v
}

/// `e^(i/2^8)` for `|i| <= 2^7` (`coarse[i + 2^7]`) and `e^(j/2^16)` for
/// `|j| < 2^8` (`fine[j + 2^8 - 1]`).
fn build_exp(prec: u64) -> Table<Ball> {
    // Up to 2^8 products, 8 bits of drift.
    let q = prec + 16;
    let mul = |a: &Ball, b: &Ball| a.mul(b, q);
    let side = |k: i64, n: usize| {
        let up = powers(one(), &exp_small(&pow2(-k), q), n, mul);
        let down = powers(one(), &exp_small(&pow2(-k).neg(), q), n, mul);
        down.into_iter().skip(1).rev().chain(up).collect::<Vec<_>>()
    };
    Table {
        prec,
        coarse: side(R1, 1 << (R1 - 1)),
        fine: side(R1 + R2, (1 << R2) - 1),
    }
}

/// `(sin, cos)(i/2^8)` for `0 <= i <= 202` (`π/4 < 202/2^8`) and
/// `(sin, cos)(j/2^16)` for `0 <= j < 2^8`.
fn build_sin_cos(prec: u64) -> Table<(Ball, Ball)> {
    let q = prec + 16;
    let rotate = |(s, c): &(Ball, Ball), (s1, c1): &(Ball, Ball)| {
        (
            s.mul(c1, q).add(&c.mul(s1, q), q),
            c.mul(c1, q).sub(&s.mul(s1, q), q),
        )
    };
    let side =
        |k: i64, n: usize| powers((Ball::ZERO, one()), &sin_cos_small(&pow2(-k), q), n, rotate);
    Table {
        prec,
        coarse: side(R1, 202),
        fine: side(R1 + R2, (1 << R2) - 1),
    }
}

/// `log(1 + i/2^8)` for `0 <= i <= 2^8` and `log(1 + j/2^16)` for `0 <= j
/// <= 2^8`, as `2 atanh(z)`, `z = c / (2 + c) <= 1/3`. `atanh' = 1 / (1 -
/// z^2) <= 9/8` there, so `z`'s radius adds at most twice itself.
fn build_log(prec: u64) -> Table<Ball> {
    let q = prec + 16;
    let entry = |i: u64, k: i64| {
        let c = Ball::exact(at(i, k));
        let z = c.div(&c.add(&Ball::from_integer(&Integer::from(2)), q), q);
        atanh_series(z.mid(), q)
            .add_error(&z.rad().mul_2exp(&Integer::ONE))
            .mul_2exp(&Integer::ONE)
    };
    Table {
        prec,
        coarse: (0..=1 << R1).map(|i| entry(i, R1)).collect(),
        fine: (0..=1 << R2).map(|j| entry(j, R1 + R2)).collect(),
    }
}

/// `atan(i/2^8)` for `0 <= i <= 2^8` and `atan(j/2^16)` for `0 <= j <= 2^8`.
fn build_atan(prec: u64) -> Table<Ball> {
    let q = prec + 16;
    Table {
        prec,
        coarse: (0..=1 << R1)
            .map(|i| atan_halving(&Ball::exact(at(i, R1)), q))
            .collect(),
        fine: (0..=1 << R2)
            .map(|j| atan_series(&at(j, R1 + R2), q))
            .collect(),
    }
}

/// `usize` of a table index.
fn index(i: u64) -> usize {
    usize::try_from(i).expect("a table index")
}

/// A ball containing `exp(t)` for an exact `|t| < 1/2`, at about `wp <=
/// TABLE_MAX_PREC` bits.
pub(super) fn exp(t: &Float, wp: u64) -> Ball {
    let tab = EXP.get(wp + 16, build_exp);
    let a = digits(t, R1 + R2);
    let (a1, a2) = (index(a >> R2), index(a & ((1 << R2) - 1)));
    let half = 1usize << (R1 - 1);
    let last = (1usize << R2) - 1;
    let (c, f) = if t.is_negative() {
        (&tab.coarse[half - a1], &tab.fine[last - a2])
    } else {
        (&tab.coarse[half + a1], &tab.fine[last + a2])
    };
    // x = t - (±a) / 2^16, exact, |x| < 2^-16 with t's sign.
    let shift = at(a, R1 + R2);
    let x = t.sub(&if t.is_negative() { shift.neg() } else { shift });
    let q = wp + 8;
    c.mul(f, q).mul(&exp_small(&x, q), q)
}

/// Balls containing `sin(t)` and `cos(t)` for an exact `|t| < 0.789` (just
/// above `π/4`), at about `wp <= TABLE_MAX_PREC` bits relative to each;
/// `None` beyond the table.
pub(super) fn sin_cos(t: &Float, wp: u64) -> Option<(Ball, Ball)> {
    let tab = SIN_COS.get(wp + 16, build_sin_cos);
    let a = digits(t, R1 + R2);
    let (a1, a2) = (index(a >> R2), index(a & ((1 << R2) - 1)));
    let (s1, c1) = tab.coarse.get(a1)?;
    let (s2, c2) = &tab.fine[a2];
    let x = t.abs().sub(&at(a, R1 + R2));
    let q = wp + 16;
    let (sx, cx) = sin_cos_small(&x, q);
    // Angles a1/2^8 + a2/2^16 + x, all >= 0 and summing to at most 0.79:
    // the sines add positive terms, the cosines (>= 0.7) subtract at most
    // 2^-8 of themselves.
    let s12 = s1.mul(c2, q).add(&c1.mul(s2, q), q);
    let c12 = c1.mul(c2, q).sub(&s1.mul(s2, q), q);
    let s = s12.mul(&cx, q).add(&c12.mul(&sx, q), q);
    let c = c12.mul(&cx, q).sub(&s12.mul(&sx, q), q);
    Some((if t.is_negative() { s.neg() } else { s }, c))
}

/// A ball containing `atanh(z)` for every `z` in a ball with `|z| <= 1/2`,
/// at about `wp` bits: `atanh' <= 4/3` there.
fn atanh_ball(z: &Ball, wp: u64) -> Ball {
    atanh_series(z.mid(), wp).add_error(&z.rad().mul_2exp(&Integer::ONE))
}

/// A ball containing `log(1 + u)` for an exact `0 <= u < 1`, at about `wp
/// <= TABLE_MAX_PREC` bits relative to the result. All three terms are
/// non-negative but for `u2`, which is at worst a rounding below 0.
pub(super) fn log1p(u: &Float, wp: u64) -> Ball {
    let tab = LOG.get(wp + 16, build_log);
    let q = wp + 24;
    let a1 = digits(u, R1);
    let n1 = u.mul_2exp(&Integer::from(R1)).sub(&at(a1, 0));
    let u1 = Ball::exact(n1).div(&Ball::from_integer(&Integer::from((1 << R1) + a1)), q);
    // u1 < 2^-8, so a2 <= 2^8 (at 2^8 only by rounding).
    let a2 = if u1.mid().is_negative() {
        0
    } else {
        digits(u1.mid(), R1 + R2).min(1 << R2)
    };
    let u2 = u1
        .mul_2exp(&Integer::from(R1 + R2))
        .sub(&Ball::from_integer(&Integer::from(a2)), q)
        .div(
            &Ball::from_integer(&Integer::from((1u64 << (R1 + R2)) + a2)),
            q,
        );
    // log(1 + u2) = 2 atanh(u2 / (2 + u2)).
    let two = Ball::from_integer(&Integer::from(2));
    let z = u2.div(&two.add(&u2, q), q);
    let l = atanh_ball(&z, q).mul_2exp(&Integer::ONE);
    tab.coarse[index(a1)]
        .add(&tab.fine[index(a2)], q)
        .add(&l, q)
}

/// A ball containing `atan(x)` for an exact `0 <= x <= 1`, at about `wp <=
/// TABLE_MAX_PREC` bits relative to the result.
pub(super) fn atan(x: &Float, wp: u64) -> Ball {
    let tab = ATAN.get(wp + 16, build_atan);
    let q = wp + 24;
    let a1 = digits(x, R1);
    let k1 = at(a1, 0);
    // x1 = (2^8 x - a1) / (2^8 + a1 x), numerator and denominator exact.
    let num = x.mul_2exp(&Integer::from(R1)).sub(&k1);
    let den = Float::from(1i64 << R1).add(&k1.mul(x));
    let x1 = Ball::exact(num).div(&Ball::exact(den), q);
    let a2 = if x1.mid().is_negative() {
        0
    } else {
        digits(x1.mid(), R1 + R2).min(1 << R2)
    };
    let k2 = Ball::from_integer(&Integer::from(a2));
    let x2 = x1.mul_2exp(&Integer::from(R1 + R2)).sub(&k2, q).div(
        &Ball::from_integer(&Integer::from(1i64 << (R1 + R2))).add(&k2.mul(&x1, q), q),
        q,
    );
    if !x2.is_finite() {
        return Ball::INDETERMINATE;
    }
    // atan is 1-Lipschitz.
    let r = atan_series(x2.mid(), q).add_error(x2.rad());
    tab.coarse[index(a1)]
        .add(&tab.fine[index(a2)], q)
        .add(&r, q)
}
