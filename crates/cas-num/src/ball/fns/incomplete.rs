//! The regularized incomplete gamma and beta functions of real balls
//! (CAS-DOC-46 §2, CAS-468): `P(a, x) = γ(a, x)/Γ(a)`, `Q(a, x) = Γ(a,
//! x)/Γ(a) = 1 − P(a, x)` and `I_x(a, b) = B_x(a, b)/B(a, b)`.
//!
//! At exact points, `a, b > 0` and `x ≥ 0` (`x ≤ 1` for `I`):
//!
//! - **The series** `P(a, x) = xᵃ e^(−x)/Γ(a + 1) Σ_k xᵏ/((a + 1) ⋯ (a +
//!   k))` (DLMF 8.7.1), of positive terms: once `2x ≤ a + k + 1` the ratio
//!   of terms stays below 1/2, so the rest is at most the last term; summed
//!   until that is below `2^-(wp + 12)` of the sum.
//! - **The asymptotic expansion** `Γ(a, x) = x^(a−1) e^(−x) Σ_{n<N} u_n`,
//!   `u_n = (a − 1)(a − 2) ⋯ (a − n)/xⁿ`, whose remainder is at most the
//!   first term left out once `N ≥ a − 1` (DLMF 8.11.2 and the text after
//!   it): for `x` past `(wp + 16) log 2 + 2 abs(a) + 4`, while the terms
//!   still fall.
//! - `P` is the series, or `1 − Q` past that bound; `Q` the expansion past
//!   it, or `1 − P` at more bits, doubled until the difference keeps the
//!   precision (it cancels where `P` is near 1).
//! - **`I_x(a, b) = xᵃ (1 − x)ᵇ/(a B(a, b)) Σ_k (a + b)_k/(a + 1)_k xᵏ`**
//!   (DLMF 8.17.8), of positive terms, for `x ≤ (a + 1)/(a + b + 2)`, and
//!   `1 − I_(1−x)(b, a)` above (DLMF 8.17.4). The term ratio `x (a + b +
//!   k)/(a + k + 1)` is monotonic in `k` and tends to `x`, so past `k` it is
//!   at most `r = max(ratio_k, x) < 1` and the rest at most `t_k r/(1 −
//!   r)`.
//!
//! Past [`MAX_TERMS`] terms a value is the indeterminate ball. For balls,
//! the monotonicity (`P` falls in `a` and rises in `x`; `I_x(a, b)` falls
//! in `a` and rises in `b` and `x`) bounds the values over the box by
//! those at two corners, each an exact point. A box that meets `a ≤ 0`,
//! `b ≤ 0`, `x < 0` or, for `I`, `x > 1`, is indeterminate.

use super::gamma::log_gamma_positive;
use super::{one, upper, working};
use crate::Integer;
use crate::ball::{Ball, Float, Mag, Round};

/// The most terms of a series or an expansion.
const MAX_TERMS: u64 = 1 << 20;

/// `2^-(wp + 12)`, the relative size of a last term.
fn small(wp: u64) -> Mag {
    Mag::power_of_two(&-Integer::from(wp + 12))
}

fn integer(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

/// Whether `t` is below `small` of `sum`.
fn negligible(t: &Ball, sum: &Ball, small: &Mag) -> bool {
    upper(t) <= Mag::from_float_lower(sum.mid()).mul_lower(small)
}

/// `P(a, x)` for exact `a > 0` and `x > 0` by the series, at about `wp`
/// bits; `None` past [`MAX_TERMS`].
fn p_series(a: &Float, x: &Float, wp: u64) -> Option<Ball> {
    let q = wp + 32;
    let (ab, xb) = (Ball::exact(a.clone()), Ball::exact(x.clone()));
    let twice = x.mul_2exp(&Integer::ONE);
    let small = small(wp);
    let mut term = one();
    let mut sum = one();
    let mut k = 0;
    loop {
        k += 1;
        if k > MAX_TERMS {
            return None;
        }
        term = term.mul(&xb, q).div(&ab.add(&integer(k), q), q);
        sum = sum.add(&term, q);
        // x/(a + k + 1) ≤ 1/2, and so for every later ratio.
        let halving = twice <= a.add(&Float::from(Integer::from(k + 1)));
        if halving && negligible(&term, &sum, &small) {
            break;
        }
    }
    let sum = sum.add_error(&upper(&term));
    let log_scale = ab
        .mul(&xb.log(q), q)
        .sub(&xb, q)
        .sub(&log_gamma_positive(&ab.add(&one(), q), q), q);
    Some(log_scale.exp(q).mul(&sum, q))
}

/// Whether `x` is past the bound of the asymptotic expansion.
#[allow(clippy::cast_precision_loss, reason = "a threshold")]
fn large(a: &Float, x: &Float, wp: u64) -> bool {
    let (a, x) = (a.to_f64(), x.to_f64());
    x.is_finite() && x >= (wp as f64 + 16.0) * std::f64::consts::LN_2 + 2.0 * a.abs() + 4.0
}

/// `Q(a, x)` for exact `a > 0` and `x` past [`large`] by the asymptotic
/// expansion, at about `wp` bits; `None` when its terms grow first.
fn q_asymptotic(a: &Float, x: &Float, wp: u64) -> Option<Ball> {
    let q = wp + 32;
    let (ab, xb) = (Ball::exact(a.clone()), Ball::exact(x.clone()));
    let small = small(wp);
    let mut term = one();
    let mut sum = one();
    let mut n: u64 = 0;
    let rest = loop {
        n += 1;
        if n > MAX_TERMS {
            return None;
        }
        // (a − n)/x from one term to the next.
        term = term.mul(&ab.sub(&integer(n), q), q).div(&xb, q);
        let past = Float::from(Integer::from(n)) >= a.sub(&Float::from(Integer::ONE));
        if past && negligible(&term, &sum, &small) {
            break upper(&term);
        }
        // Growing terms: |a − n − 1| ≥ x.
        let next = a.sub(&Float::from(Integer::from(n + 1))).abs();
        if next >= *x {
            return None;
        }
        sum = sum.add(&term, q);
    };
    let sum = sum.add_error(&rest);
    let log_scale = ab
        .sub(&one(), q)
        .mul(&xb.log(q), q)
        .sub(&xb, q)
        .sub(&log_gamma_positive(&ab, q), q);
    Some(log_scale.exp(q).mul(&sum, q))
}

/// `1 − f(w)` at `w = wp` bits and more, doubled until it keeps `wp` bits,
/// up to `4 wp + 4096` more.
fn complement(wp: u64, f: impl Fn(u64) -> Option<Ball>) -> Ball {
    let want = i64::try_from(wp).unwrap_or(i64::MAX);
    let mut extra = 0;
    loop {
        let Some(v) = f(wp + extra) else {
            return Ball::INDETERMINATE;
        };
        let c = one().sub(&v, wp + extra + 32);
        if c.rel_accuracy_bits() >= want || extra >= 4 * wp + 4096 || !c.is_finite() {
            return c;
        }
        extra = (2 * extra).max(64).min(4 * wp + 4096);
    }
}

/// `P(a, x)` at an exact point, `a > 0`, `x ≥ 0`.
fn p_exact(a: &Float, x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    if large(a, x, wp)
        && let Some(q) = q_asymptotic(a, x, wp)
    {
        return one().sub(&q, wp + 32);
    }
    p_series(a, x, wp).unwrap_or(Ball::INDETERMINATE)
}

/// `Q(a, x)` at an exact point, `a > 0`, `x ≥ 0`.
fn q_exact(a: &Float, x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return one();
    }
    if large(a, x, wp)
        && let Some(q) = q_asymptotic(a, x, wp)
    {
        return q;
    }
    complement(wp, |w| p_series(a, x, w))
}

/// `I_x(a, b)` for exact `a, b > 0` and `0 < x ≤ (a + 1)/(a + b + 2)` by
/// the series; `None` past [`MAX_TERMS`].
fn i_series(a: &Float, b: &Float, x: &Float, wp: u64) -> Option<Ball> {
    let q = wp + 32;
    let (ab, bb, xb) = (
        Ball::exact(a.clone()),
        Ball::exact(b.clone()),
        Ball::exact(x.clone()),
    );
    let sum_ab = ab.add(&bb, q);
    let small = small(wp);
    let mut term = one();
    let mut sum = one();
    let mut k = 0;
    let rest = loop {
        // The ratio x (a + b + k)/(a + k + 1) to the next term.
        let ratio = xb
            .mul(&sum_ab.add(&integer(k), q), q)
            .div(&ab.add(&integer(k + 1), q), q);
        // Past k the ratios are at most r = max(ratio_k, x) < 1.
        let r = upper(&ratio).max(Mag::from_float(x));
        if r < Mag::ONE && negligible(&term, &sum, &small) {
            let rest = Ball::exact(r.to_float()?);
            let factor = rest.div(&one().sub(&rest, q), q);
            break upper(&term).mul(&upper(&factor));
        }
        k += 1;
        if k > MAX_TERMS {
            return None;
        }
        term = term.mul(&ratio, q);
        sum = sum.add(&term, q);
    };
    let sum = sum.add_error(&rest);
    let log_beta = log_gamma_positive(&ab, q)
        .add(&log_gamma_positive(&bb, q), q)
        .sub(&log_gamma_positive(&sum_ab, q), q);
    let one_minus = Ball::exact(Float::from(Integer::ONE).sub(x));
    let log_scale = ab
        .mul(&xb.log(q), q)
        .add(&bb.mul(&one_minus.log(q), q), q)
        .sub(&ab.log(q), q)
        .sub(&log_beta, q);
    Some(log_scale.exp(q).mul(&sum, q))
}

/// Whether `x ≤ (a + 1)/(a + b + 2)`, where the series of `I_x(a, b)` is
/// taken.
fn below_mean(a: &Float, b: &Float, x: &Float) -> bool {
    let one = Float::from(Integer::ONE);
    let two = Float::from(Integer::from(2));
    x.mul(&a.add(b).add(&two)) <= a.add(&one)
}

/// `I_x(a, b)` at an exact point, `a, b > 0`, `0 ≤ x ≤ 1`.
fn i_exact(a: &Float, b: &Float, x: &Float, wp: u64) -> Ball {
    let one_f = Float::from(Integer::ONE);
    if x.is_zero() {
        return Ball::ZERO;
    }
    if *x == one_f {
        return one();
    }
    if below_mean(a, b, x) {
        return i_series(a, b, x, wp).unwrap_or(Ball::INDETERMINATE);
    }
    let y = one_f.sub(x);
    complement(wp, |w| i_series(b, a, &y, w))
}

/// The ends of a finite ball, rounded outward at `prec` bits; an exact
/// ball's are its point, as it is.
fn ends(b: &Ball, prec: u64) -> Option<(Float, Float)> {
    if b.is_exact() {
        return Some((b.mid().clone(), b.mid().clone()));
    }
    let r = b.rad().to_float()?;
    Some((
        b.mid().sub_round(&r, prec, Round::Floor).0,
        b.mid().add_round(&r, prec, Round::Ceil).0,
    ))
}

/// Whether every point of the ball is above 0, and its ends.
fn positive(b: &Ball, prec: u64) -> Option<(Float, Float)> {
    let (lo, hi) = ends(b, prec)?;
    (!lo.is_negative() && !lo.is_zero()).then_some((lo, hi))
}

/// Whether every point of the ball is at least 0, and its ends.
fn nonnegative(b: &Ball, prec: u64) -> Option<(Float, Float)> {
    let (lo, hi) = ends(b, prec)?;
    (!lo.is_negative()).then_some((lo, hi))
}

/// The value at one point, or the union of the values at the two corners
/// that bound it.
fn over(prec: u64, low: Ball, high: &Ball, exact: bool) -> Ball {
    let value = if exact {
        low
    } else {
        low.union(high, prec + 32)
    };
    value.round_to_accuracy(prec)
}

impl Ball {
    /// `P(a, x) = γ(a, x)/Γ(a)`, the regularized lower incomplete gamma
    /// function, for every `a > 0` and `x ≥ 0` in the balls, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn gamma_p(a: &Ball, x: &Ball, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::gamma_p: precision 0");
        let (Some((a_lo, a_hi)), Some((x_lo, x_hi))) =
            (positive(a, prec + 64), nonnegative(x, prec + 64))
        else {
            return Ball::INDETERMINATE;
        };
        let wp = working(prec);
        let exact = a.is_exact() && x.is_exact();
        let low = p_exact(&a_hi, &x_lo, wp);
        let high = if exact {
            Ball::ZERO
        } else {
            p_exact(&a_lo, &x_hi, wp)
        };
        over(prec, low, &high, exact)
    }

    /// `Q(a, x) = Γ(a, x)/Γ(a) = 1 − P(a, x)`, the regularized upper
    /// incomplete gamma function, for every `a > 0` and `x ≥ 0` in the
    /// balls, at `prec` bits, relative to its value however small.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn gamma_q(a: &Ball, x: &Ball, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::gamma_q: precision 0");
        let (Some((a_lo, a_hi)), Some((x_lo, x_hi))) =
            (positive(a, prec + 64), nonnegative(x, prec + 64))
        else {
            return Ball::INDETERMINATE;
        };
        let wp = working(prec);
        let exact = a.is_exact() && x.is_exact();
        let low = q_exact(&a_lo, &x_hi, wp);
        let high = if exact {
            Ball::ZERO
        } else {
            q_exact(&a_hi, &x_lo, wp)
        };
        over(prec, low, &high, exact)
    }

    /// `I_x(a, b) = B_x(a, b)/B(a, b)`, the regularized incomplete beta
    /// function, for every `a, b > 0` and `0 ≤ x ≤ 1` in the balls, at
    /// `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn beta_i(a: &Ball, b: &Ball, x: &Ball, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::beta_i: precision 0");
        let (Some((a_lo, a_hi)), Some((b_lo, b_hi)), Some((x_lo, x_hi))) = (
            positive(a, prec + 64),
            positive(b, prec + 64),
            nonnegative(x, prec + 64),
        ) else {
            return Ball::INDETERMINATE;
        };
        if x_hi > Float::from(Integer::ONE) {
            return Ball::INDETERMINATE;
        }
        let wp = working(prec);
        let exact = a.is_exact() && b.is_exact() && x.is_exact();
        let low = i_exact(&a_hi, &b_lo, &x_lo, wp);
        let high = if exact {
            Ball::ZERO
        } else {
            i_exact(&a_lo, &b_hi, &x_hi, wp)
        };
        over(prec, low, &high, exact)
    }
}

#[cfg(test)]
mod tests;
