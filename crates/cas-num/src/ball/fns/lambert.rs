//! Lambert's W of real balls (CAS-DOC-31 §2.3, CAS-352): the inverse of
//! `f(w) = w e^w` on its branch 0, `W ≥ −1` for `x ≥ −1/e`, and its branch
//! −1, `W ≤ −1` for `−1/e ≤ x < 0`.
//!
//! At an exact `x`:
//!
//! 1. A first approximation: `x` itself for a tiny `x` on branch 0; `L₁ −
//!    L₂ + L₂/L₁` with `L₁ = log |x|` and `L₂ = log |L₁|` far out, from the
//!    exponent where `f64` has no room; `−1 ± p − p²/3 + 11p³/72`, `p =
//!    √(2 (1 + e x))`, next to the branch point; polished by Newton's steps
//!    in `f64` where it has room.
//! 2. Newton's steps `w ← w − (w e^w − x)/(e^w (w + 1))` on the midpoint,
//!    the precision doubling up to the working one, then until the step is
//!    below it.
//! 3. A certificate: `f` is increasing on `[−1, ∞)` and decreasing on `(−∞,
//!    −1]`, so `W(x)` lies between `w ± ε` where the balls `f(w − ε)` and
//!    `f(w + ε)` lie on either side of `x`, an end past −1 put at −1, where
//!    `f(−1) = −1/e` is on the right side as `x > −1/e` was shown. `ε`
//!    starts from the residual over the slope, `2 |f(w) − x|/|f′(w)|`, and
//!    quadruples until the balls decide, at most 64 times.
//!
//! A ball whose lower end is not shown above `−1/e`, or on branch −1 whose
//! upper end is not shown below 0, is indeterminate. A ball with a radius
//! has the hull of its ends' values, as `W` is monotone on each branch. An
//! exact input is computed again at more bits while its relative accuracy
//! falls short, which next to the branch point, where `f′ → 0`, it does.

use std::cmp::Ordering;
use std::f64::consts::{E, LN_2};

use super::{exponent, one, refined, upper, working};
use crate::Integer;
use crate::ball::{Ball, Float, Mag, Round};

/// How `x` compares with `−1/e`, at more bits until the balls decide.
fn against_branch_point(x: &Float) -> Option<Ordering> {
    let mut p = 64 + x.bit_length();
    let xb = Ball::exact(x.clone());
    for _ in 0..8 {
        let point = Ball::exact(Float::from(-1i64)).exp(p).neg();
        if let Some(order) = xb.cmp(&point) {
            return Some(order);
        }
        p *= 2;
    }
    None
}

/// `W` next to the branch point from `p = ±√(2 (1 + e x))`, in `f64`.
fn near_branch_point(x: f64, lower: bool) -> f64 {
    let p = (2.0 * (1.0 + E * x)).max(0.0).sqrt();
    let p = if lower { -p } else { p };
    -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p
}

/// `L₁ − L₂ + L₂/L₁` for `L₁ = log |x|`, the asymptotic start far out.
fn far(ln_abs: f64) -> f64 {
    let l2 = ln_abs.abs().ln();
    ln_abs - l2 + l2 / ln_abs
}

/// Newton's steps in `f64` from `w`, while they are finite.
fn polish(x: f64, mut w: f64) -> f64 {
    for _ in 0..32 {
        let e = w.exp();
        let step = (w * e - x) / (e * (w + 1.0));
        if !step.is_finite() {
            break;
        }
        w -= step;
        if step.abs() <= 1e-15 * w.abs() {
            break;
        }
    }
    w
}

/// A first approximation of `W(x)` on the branch.
#[allow(clippy::cast_precision_loss, reason = "an estimate")]
fn start(x: &Float, lower: bool) -> Float {
    let e = exponent(x);
    if !lower && e < -60 {
        // W(x) = x − x² + …
        return x.clone();
    }
    if e.abs() > 1000 {
        let ln_abs = e as f64 * LN_2;
        return Float::from_f64(far(ln_abs)).expect("finite");
    }
    let xf = x.to_f64();
    let w = if xf < -0.25 {
        near_branch_point(xf, lower)
    } else if lower {
        far(xf.abs().ln())
    } else if xf > 3.0 {
        far(xf.ln())
    } else {
        xf.ln_1p()
    };
    let w = polish(xf, w);
    // Keep the start on its branch.
    let w = if lower { w.min(-1.0) } else { w.max(-1.0) };
    Float::from_f64(w).unwrap_or_else(|| Float::from(-1i64))
}

/// `(w e^w, e^w)` for an exact `w` at `q` bits.
fn f_and_exp(w: &Float, q: u64) -> (Ball, Ball) {
    let wb = Ball::exact(w.clone());
    let ew = wb.exp(q);
    (wb.mul(&ew, q), ew)
}

/// The midpoint `W(x)` at about `wp` bits by Newton's steps, kept on the
/// branch: a step past −1 goes halfway to it instead, and at −1 itself,
/// where `f′ = 0`, the midpoint moves off to the branch's side.
fn newton(x: &Float, lower: bool, wp: u64) -> Float {
    let xb = Ball::exact(x.clone());
    let minus_one = Float::from(-1i64);
    let mut w = start(x, lower);
    let target = wp + 32;
    let mut p = 64u64;
    let mut at_target = 0;
    for _ in 0..256 {
        p = (2 * p).min(target);
        let q = p + 32;
        let (f, ew) = f_and_exp(&w, q);
        let slope = ew.mul(&Ball::exact(w.clone()).add(&one(), q), q);
        let step = f.sub(&xb, q).div(&slope, q);
        if !step.is_finite() {
            let nudge = Float::from(1i64).mul_2exp(&-Integer::from(p / 2));
            w = if lower { w.sub(&nudge) } else { w.add(&nudge) };
            continue;
        }
        let next = w.sub_round(step.mid(), q, Round::Nearest).0;
        let crossed = if lower {
            next > minus_one
        } else {
            next < minus_one
        };
        w = if crossed {
            w.add(&minus_one).mul_2exp(&-Integer::ONE)
        } else {
            next
        };
        if p == target {
            at_target += 1;
            let small = Mag::from_float(&w).mul(&Mag::power_of_two(&-Integer::from(wp + 24)));
            if (!crossed && upper(&step) <= small) || at_target >= 12 {
                break;
            }
        }
    }
    w
}

/// `W(x)` for an exact `x ≠ 0` on the branch, its domain shown, at about
/// `wp` bits: Newton's midpoint and its certificate.
fn certified(x: &Float, lower: bool, wp: u64) -> Ball {
    if x.is_zero() {
        return Ball::ZERO;
    }
    let w = newton(x, lower, wp);
    let q = wp + 32;
    let xb = Ball::exact(x.clone());
    // ε from the residual over the slope, and at least 2^-wp of w.
    let (f, ew) = f_and_exp(&w, q);
    let residual = upper(&f.sub(&xb, q));
    let slope = ew.mul(&Ball::exact(w.clone()).add(&one(), q), q);
    let slope_low = Mag::from_float_lower(&slope.mid().abs()).sub_lower(slope.rad());
    let floor = Mag::from_float(&w).mul(&Mag::power_of_two(&-Integer::from(wp)));
    let mut eps = if slope_low.is_zero() {
        floor.clone()
    } else {
        residual.mul_2exp(&Integer::ONE).div(&slope_low).add(&floor)
    };
    let minus_one = Float::from(-1i64);
    for _ in 0..64 {
        let Some(e) = eps.to_float() else {
            return Ball::INDETERMINATE;
        };
        let mut lo = w.sub_round(&e, q, Round::Floor).0;
        let mut hi = w.add_round(&e, q, Round::Ceil).0;
        // An end past −1 goes to −1, where f(−1) = −1/e is on its side.
        let (mut lo_clamped, mut hi_clamped) = (false, false);
        if !lower && lo < minus_one {
            lo = minus_one.clone();
            lo_clamped = true;
        }
        if lower && hi > minus_one {
            hi = minus_one.clone();
            hi_clamped = true;
        }
        let side = |t: &Float, clamped: bool, want: Ordering| {
            clamped || f_and_exp(t, q).0.cmp(&xb) == Some(want)
        };
        let (below, above) = if lower {
            (Ordering::Greater, Ordering::Less)
        } else {
            (Ordering::Less, Ordering::Greater)
        };
        if side(&lo, lo_clamped, below) && side(&hi, hi_clamped, above) {
            let mid = lo.add(&hi).mul_2exp(&-Integer::ONE);
            let rad = Mag::from_float(&hi.sub(&mid));
            return Ball::new(mid, rad);
        }
        eps = eps.mul_2exp(&Integer::from(2));
    }
    Ball::INDETERMINATE
}

impl Ball {
    /// `W(x)` on branch 0 (`W ≥ −1`), or on branch −1 (`W ≤ −1`) for
    /// `lower`, for every `x` in the ball, at `prec` bits; indeterminate
    /// where the ball is not shown inside the branch's domain, `x ≥ −1/e`,
    /// and `x < 0` for branch −1.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn lambert_w(&self, lower: bool, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::lambert_w: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let Some(r) = self.rad().to_float() else {
            return Ball::INDETERMINATE;
        };
        let m = self.mid();
        let bits = m.bit_length().max(62) + 2;
        let lo = m.sub_round(&r, bits, Round::Floor).0;
        let hi = m.add_round(&r, bits, Round::Ceil).0;
        if against_branch_point(&lo) != Some(Ordering::Greater)
            || (lower && hi.sign() != Ordering::Less)
        {
            return Ball::INDETERMINATE;
        }
        if self.is_exact() {
            return refined(self, prec, |wp| certified(m, lower, wp)).round_to_accuracy(prec);
        }
        let wp = working(prec);
        let (a, b) = (certified(&lo, lower, wp), certified(&hi, lower, wp));
        a.union(&b, wp).round_to_accuracy(prec)
    }
}

#[cfg(test)]
mod tests;
