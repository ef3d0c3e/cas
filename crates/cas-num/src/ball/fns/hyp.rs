//! Hyperbolic functions and their inverses (CAS-DOC-5 §7.3), written so that
//! nothing cancels near 0:
//!
//! ```text
//! sinh x  = e (e + 2) / (2 (e + 1)),      e = expm1(|x|), odd
//! cosh x  = (E + 1/E) / 2,                E = exp(|x|)
//! tanh x  = e2 / (e2 + 2),                e2 = expm1(2|x|), odd
//! asinh x = log1p(|x| + x^2 / (1 + √(1 + x^2))), odd
//! acosh x = log1p(t + √(t (t + 2))),      t = x - 1
//! atanh x = log1p(2|x| / (1 - |x|)) / 2,  odd
//! ```
//!
//! At the exact midpoint, then widened by a derivative bound on the ball
//! (`z = |m| + ρ`): `ρ cosh z <= ρ (e^z + 1) / 2` for `sinh`, `ρ sinh z <= ρ
//! e^z / 2` for `cosh`, `ρ` for `tanh` and `asinh`, `ρ / (1 - z)` for
//! `atanh`. `acosh`, whose derivative is unbounded at 1, is the ball
//! composition of the formula. Outside a domain, or for `sinh`/`cosh` of
//! arguments past `exp`'s range, the indeterminate ball.

use super::exp::{HUGE_EXP, exp_exact, expm1_exact};
use super::{one, working};
use crate::Integer;
use crate::ball::{Ball, Float, Mag};

/// `|m| + ρ` from above.
fn reach(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// Whether `|m|`'s exponent is beyond `exp`'s reduction.
fn huge(m: &Float) -> bool {
    !m.is_zero() && *m.exponent() > Integer::from(HUGE_EXP)
}

/// `v` negated when `neg`.
fn signed(v: Ball, neg: bool) -> Ball {
    if neg { v.neg() } else { v }
}

impl Ball {
    /// `sinh x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn sinh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::sinh: precision 0");
        let m = self.mid();
        if !self.is_finite() || huge(m) {
            return Ball::INDETERMINATE;
        }
        let core = if m.is_zero() {
            Ball::ZERO
        } else {
            let wp = working(prec);
            let e = expm1_exact(&m.abs(), wp);
            let num = e.mul(&e.add(&Ball::from_integer(&Integer::from(2)), wp), wp);
            let v = num.div(&e.add(&one(), wp), wp).mul_2exp(&-Integer::ONE);
            signed(v, m.is_negative())
        };
        let spread = self
            .rad()
            .mul(&reach(self).exp().add(&Mag::ONE).mul_2exp(&-Integer::ONE));
        core.add_error(&spread).round_to_accuracy(prec)
    }

    /// `cosh x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn cosh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::cosh: precision 0");
        let m = self.mid();
        if !self.is_finite() || huge(m) {
            return Ball::INDETERMINATE;
        }
        let wp = working(prec);
        let big = exp_exact(&m.abs(), wp);
        let core = big.add(&one().div(&big, wp), wp).mul_2exp(&-Integer::ONE);
        let spread = self.rad().mul(&reach(self).exp().mul_2exp(&-Integer::ONE));
        core.add_error(&spread).round_to_accuracy(prec)
    }

    /// `tanh x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn tanh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::tanh: precision 0");
        let m = self.mid();
        if !self.is_finite() {
            return Ball::new(Float::ZERO, Mag::ONE);
        }
        let wp = working(prec);
        let core = if m.is_zero() {
            Ball::ZERO
        } else if huge(&m.mul_2exp(&Integer::ONE)) {
            // 1 - tanh x <= 2 e^(-2x) < 2^-x, far below 2^-(wp + 10).
            let one_ish = one().add_error(&Mag::power_of_two(&-Integer::from(wp + 10)));
            signed(one_ish, m.is_negative())
        } else {
            let e2 = expm1_exact(&m.abs().mul_2exp(&Integer::ONE), wp);
            let v = e2.div(&e2.add(&Ball::from_integer(&Integer::from(2)), wp), wp);
            signed(v, m.is_negative())
        };
        core.add_error(self.rad()).round_to_accuracy(prec)
    }

    /// `asinh x` for every `x` in the ball, at `prec` bits.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn asinh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::asinh: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        let core = if m.is_zero() {
            Ball::ZERO
        } else {
            let wp = working(prec);
            let a = Ball::exact(m.abs());
            let a2 = a.mul(&a, wp);
            let s = one().add(&a2, wp).sqrt(wp);
            let t = a.add(&a2.div(&one().add(&s, wp), wp), wp);
            signed(t.log1p(wp), m.is_negative())
        };
        core.add_error(self.rad()).round_to_accuracy(prec)
    }

    /// `acosh x` for every `x` in the ball, at `prec` bits; indeterminate
    /// unless the ball lies in `x >= 1`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn acosh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::acosh: precision 0");
        let wp = working(prec);
        let t = self.sub(&one(), wp);
        let root = t
            .mul(&t.add(&Ball::from_integer(&Integer::from(2)), wp), wp)
            .sqrt(wp);
        if !root.is_finite() {
            return Ball::INDETERMINATE;
        }
        t.add(&root, wp).log1p(prec)
    }

    /// `atanh x` for every `x` in the ball, at `prec` bits; indeterminate
    /// unless the ball lies in `(-1, 1)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn atanh(&self, prec: u64) -> Ball {
        assert!(prec >= 1, "Ball::atanh: precision 0");
        if !self.is_finite() {
            return Ball::INDETERMINATE;
        }
        let m = self.mid();
        if m.cmp_abs(&Float::from(1i64)) != std::cmp::Ordering::Less {
            return Ball::INDETERMINATE;
        }
        // 1 - z from below, z = |m| + ρ (1 - |m| rounded down, not exact: for
        // a tiny |m| that would be as long as its exponent); 0 when the ball
        // reaches ±1. atanh' = 1 / ((1 - z)(1 + z)) <= 1 / (1 - z).
        let one_minus = Float::from(1i64)
            .sub_round(&m.abs(), 64, crate::ball::Round::Floor)
            .0;
        let gap = Mag::from_float_lower(&one_minus).sub_lower(self.rad());
        if gap.is_zero() {
            return Ball::INDETERMINATE;
        }
        let core = if m.is_zero() {
            Ball::ZERO
        } else {
            let wp = working(prec);
            let a = Ball::exact(m.abs());
            let q = a.mul_2exp(&Integer::ONE).div(&one().sub(&a, wp), wp);
            signed(q.log1p(wp).mul_2exp(&-Integer::ONE), m.is_negative())
        };
        core.add_error(&self.rad().div(&gap))
            .round_to_accuracy(prec)
    }
}
