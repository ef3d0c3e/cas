//! The bit-burst algorithm for high precision (CAS-DOC-5 §7.2; Brent 1976,
//! as in Arb's `arb_exp_arf_bb` and `arb_atan_arf_bb`).
//!
//! An exact `|t| < 1` splits into chunks of doubling width,
//!
//! ```text
//! t = Σ_j p_j / 2^(b_(j+1)),   b_0 = 0, b_1 = 8, b_(j+1) = 2 b_j,
//! ```
//!
//! with `|p_j| < 2^(b_(j+1) - b_j)`, so `x_j = p_j / 2^(b_(j+1))` is below
//! `2^-b_j` and its series needs about `p / b_j` terms of `b_j`-bit
//! integers: binary splitting ([`bsplit`]) sums each in `O(M(p) log p)`, and
//! there are `log2(p / 8)` chunks. Chunks stop a few bits past the working
//! precision; the rest goes to the baseline's tiny-argument bound.
//!
//! - `exp(t) = Π exp(x_j)`, each from `Σ x^k / k!` (`p(k) = p_j`, `q(k) = k
//!   2^b`), tail `<= 2 |x|^n / n!` for `|x| <= 1/2`.
//! - `sin`, `cos` by angle addition from each chunk's `cos x = Σ (-x^2)^k /
//!   (2k)!` and `sin x = x Σ (-x^2)^k / (2k+1)!`; alternating series of
//!   decreasing terms for `|x| < 1`, so each tail is at most its first term.
//! - `atan x = atan(c) + atan((x - c) / (1 + x c))` with `c` the next chunk
//!   of the current argument (a ball, so the residual's rounding is in its
//!   radius) and `atan c = Σ (-1)^k c^(2k+1) / (2k+1)` (`p(0) = p_j`, `q(0) =
//!   2^b`, `p(k) = -p_j^2`, `q(k) = 2^(2b)`, `b(k) = 2k + 1`), tail at most
//!   `|c|^(2n+1)`. Arguments above 1/2 first take `atan x = π/4 - atan((1 -
//!   x) / (1 + x))`.
//!
//! `log` needs no series of its own: its Newton iteration and certificate
//! run on `exp` ([`super::log`]), which is bit-burst at these precisions.

use super::atan::atan_series;
use super::exp::{exp_small, to_fixed};
use super::factorial_lower;
use super::trig::sin_cos_small;
use crate::ball::bsplit::{self, Series};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Natural};

/// The first chunk boundary `b_1`.
const FIRST: u64 = 8;

/// `p / 2^b` as an exact float.
fn at(p: &Integer, b: u64) -> Float {
    Float::from_integer_exp(p, &-Integer::from(b))
}

/// `2^k` as an integer.
fn pow2(k: u64) -> Integer {
    Integer::from_parts(false, Natural::from(1u64) << k)
}

/// An upper bound of `|p| / 2^b`.
fn bound(p: &Integer, b: u64) -> Mag {
    Mag::from_integer(p).mul_2exp(&-Integer::from(b))
}

/// The chunks `(p_j, b_(j+1))` of an exact `|t| < 1`, zero chunks left out,
/// up to bit `limit` and the exact rest below `2^-limit`.
fn chunks(t: &Float, limit: u64) -> (Vec<(Integer, u64)>, Float) {
    let mut out = Vec::new();
    let mut rest = t.clone();
    let mut b = FIRST;
    while !rest.is_zero() {
        let p = to_fixed(&rest, b);
        if !p.is_zero() {
            rest = rest.sub(&at(&p, b));
            out.push((p, b));
        }
        if b >= limit {
            break;
        }
        b = (2 * b).min(limit);
    }
    (out, rest)
}

/// The fewest terms `n` with `n · per + log2(den(n)) >= target`, `den(n)`
/// the factorial-like growth `Π step(k)` (an estimate: it only sets the
/// accuracy, the tails are bounded exactly).
fn terms(per: u64, target: u64, step: impl Fn(u64) -> f64) -> u64 {
    #[allow(clippy::cast_precision_loss, reason = "an estimate")]
    let (per, target) = (per as f64, target as f64);
    let (mut n, mut bits) = (1u64, step(1).log2());
    #[allow(clippy::cast_precision_loss, reason = "an estimate")]
    while (n as f64) * per + bits < target {
        n += 1;
        bits += step(n).log2();
    }
    n
}

/// `-log2 |x|` from below for `x = p / 2^b`: `|p| < 2^len`.
fn per_term(p: &Integer, b: u64) -> u64 {
    b.saturating_sub(p.unsigned_abs().bit_length())
}

/// `exp(p / 2^b)`: `Σ x^k / k!`.
struct ExpChunk<'a> {
    p: &'a Integer,
    b: u64,
}

impl Series for ExpChunk<'_> {
    fn p(&self, k: u64) -> Option<Integer> {
        Some(if k == 0 { Integer::ONE } else { self.p.clone() })
    }

    fn q(&self, k: u64) -> Integer {
        if k == 0 {
            Integer::ONE
        } else {
            Integer::from_parts(false, Natural::from(k) << self.b)
        }
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }
}

/// `exp(p / 2^b)` for `|p / 2^b| <= 1/2` at `q` bits.
#[allow(clippy::cast_precision_loss, reason = "term count estimates")]
fn exp_chunk(p: &Integer, b: u64, q: u64) -> Ball {
    let n = terms(per_term(p, b), q + 2, |k| k as f64);
    let (t, d, _) = bsplit::sum(&ExpChunk { p, b }, n);
    // Tail Σ_{k >= n} |x|^k / k! <= 2 |x|^n / n!.
    let tail = bound(p, b)
        .pow(n)
        .mul_2exp(&Integer::ONE)
        .div(&factorial_lower(n));
    bsplit::quotient(&t, &d, q).add_error(&tail)
}

/// A ball containing `exp(t)` for an exact `|t| < 1/2`, at about `wp` bits.
pub(super) fn exp(t: &Float, wp: u64) -> Ball {
    let q = wp + 32;
    let (parts, rest) = chunks(t, q + 16);
    // |rest| < 2^-(q + 16): the baseline's tiny bound 1 + rest ± rest^2.
    let mut acc = exp_small(&rest, q);
    for (p, b) in &parts {
        acc = acc.mul(&exp_chunk(p, *b, q), q);
    }
    acc
}

/// `Σ (-x^2)^k / ((2k + s - 1)(2k + s))…`: the cosine (`odd = false`) or
/// `sin(x) / x` (`odd = true`) series of `x = p / 2^b`.
struct TrigChunk<'a> {
    p2: &'a Integer,
    b: u64,
    odd: bool,
}

impl Series for TrigChunk<'_> {
    fn p(&self, k: u64) -> Option<Integer> {
        Some(if k == 0 {
            Integer::ONE
        } else {
            -self.p2.clone()
        })
    }

    fn q(&self, k: u64) -> Integer {
        if k == 0 {
            return Integer::ONE;
        }
        let k = u128::from(k);
        let f = if self.odd {
            2 * k * (2 * k + 1)
        } else {
            (2 * k - 1) * 2 * k
        };
        Integer::from_parts(false, Natural::from(f) << (2 * self.b))
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }
}

/// `(sin, cos)(p / 2^b)` for `|p / 2^b| < 1` at `q` bits.
#[allow(clippy::cast_precision_loss, reason = "term count estimates")]
fn sin_cos_chunk(p: &Integer, b: u64, q: u64) -> (Ball, Ball) {
    let n = terms(2 * per_term(p, b), q + 2, |k| ((2 * k - 1) * 2 * k) as f64);
    let p2 = p * p;
    let x = bound(p, b);
    let series = |odd: bool| {
        let (t, d, _) = bsplit::sum(&TrigChunk { p2: &p2, b, odd }, n);
        // The first omitted term: |x|^(2n) / (2n)! or / (2n + 1)!.
        let tail = x.pow(2 * n).div(&factorial_lower(2 * n + u64::from(odd)));
        bsplit::quotient(&t, &d, q).add_error(&tail)
    };
    let s = Ball::exact(at(p, b)).mul(&series(true), q);
    (s, series(false))
}

/// Balls containing `sin(t)` and `cos(t)` for an exact `|t| < 1`, at about
/// `wp` bits relative to each.
pub(super) fn sin_cos(t: &Float, wp: u64) -> (Ball, Ball) {
    let q = wp + 32;
    let (parts, rest) = chunks(t, q + 16);
    // |rest| < 2^-(q + 16): the baseline's tiny bounds.
    let (mut s, mut c) = sin_cos_small(&rest, q);
    for (p, b) in &parts {
        let (s1, c1) = sin_cos_chunk(p, *b, q);
        (s, c) = (
            s.mul(&c1, q).add(&c.mul(&s1, q), q),
            c.mul(&c1, q).sub(&s.mul(&s1, q), q),
        );
    }
    (s, c)
}

/// `atan(p / 2^b)`: `Σ (-1)^k x^(2k+1) / (2k + 1)`.
struct AtanChunk<'a> {
    p: &'a Integer,
    p2: &'a Integer,
    b: u64,
}

impl Series for AtanChunk<'_> {
    fn p(&self, k: u64) -> Option<Integer> {
        Some(if k == 0 {
            self.p.clone()
        } else {
            -self.p2.clone()
        })
    }

    fn q(&self, k: u64) -> Integer {
        pow2(if k == 0 { self.b } else { 2 * self.b })
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, k: u64) -> Option<Integer> {
        Some(Integer::from(2 * u128::from(k) + 1))
    }
}

/// `atan(p / 2^b)` for `|p / 2^b| <= 1/2` at `q` bits.
fn atan_chunk(p: &Integer, b: u64, q: u64) -> Ball {
    // (2n + 1) per >= q + 2, per >= 1 as |x| <= 1/2.
    let n = (q + 2).div_ceil(2 * per_term(p, b).max(1)).max(1);
    let p2 = p * p;
    let (t, d, _) = bsplit::sum(&AtanChunk { p, p2: &p2, b }, n);
    let tail = bound(p, b).pow(2 * n + 1);
    bsplit::quotient(&t, &d, q).add_error(&tail)
}

/// A ball containing `atan(x)` for an exact `0 <= x <= 1`, at about `wp`
/// bits relative to the result.
pub(super) fn atan(x: &Float, wp: u64) -> Ball {
    let q = wp + 32;
    let one = Ball::from_integer(&Integer::ONE);
    let half = Float::from(1i64).mul_2exp(&-Integer::ONE);
    // atan x = π/4 - atan((1 - x) / (1 + x)) above 1/2, where the residual
    // is below 1/3: no cancellation, the result is above atan(1/2) > 0.46.
    let (base, mut y, flip) = if *x > half {
        let bx = Ball::exact(x.clone());
        let r = one.sub(&bx, q).div(&one.add(&bx, q), q);
        (Ball::const_pi(q).mul_2exp(&-Integer::from(2)), r, true)
    } else {
        (Ball::ZERO, Ball::exact(x.clone()), false)
    };
    let mut sum = Ball::ZERO;
    let mut b = FIRST;
    loop {
        // c = trunc(mid · 2^b) / 2^b, so the residual is below 2^-b; its
        // own rounding is in the ball.
        let p = to_fixed(y.mid(), b);
        if !p.is_zero() {
            let c = Ball::exact(at(&p, b));
            sum = sum.add(&atan_chunk(&p, b, q), q);
            y = y.sub(&c, q).div(&one.add(&y.mul(&c, q), q), q);
            if !y.is_finite() {
                return Ball::INDETERMINATE;
            }
        }
        if b >= q {
            break;
        }
        b = (2 * b).min(q);
    }
    // |y| < 2^-q: atan_series' tiny bound; atan is 1-Lipschitz.
    let sum = sum.add(&atan_series(y.mid(), q).add_error(y.rad()), q);
    if flip { base.sub(&sum, q) } else { sum }
}
