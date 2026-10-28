//! The asymptotic expansions of the Bessel functions for large arguments
//! (CAS-DOC-34 §2.1, §2.2, CAS-377), for a rational order `ν ≥ 0` and an
//! exact `x > 0`, each with the remainder DLMF proves. With `a_0(ν) = 1` and
//! `a_(k+1)(ν) = a_k(ν) (4ν² − (2k + 1)²)/(8(k + 1))`:
//!
//! - **Hankel's** (DLMF 10.17.3, 10.17.4): `J_ν(x) = √(2/(πx)) (P cos ω −
//!   Q sin ω)` and `Y_ν(x) = √(2/(πx)) (P sin ω + Q cos ω)`, `ω = x − νπ/2
//!   − π/4`, with `P = Σ_{k<ℓ} (−1)^k a_{2k}(ν)/x^(2k)` and `Q = Σ_{k<ℓ′}
//!   (−1)^k a_{2k+1}(ν)/x^(2k+1)`. Each is within its first term left out
//!   for `ℓ ≥ max(ν/2 − 1/4, 1)` and `ℓ′ ≥ max(ν/2 − 3/4, 1)` (DLMF
//!   10.17(iii)).
//! - **`K`'s** (DLMF 10.40.2): `K_ν(x) = √(π/(2x)) e^(−x) Σ_{k<ℓ}
//!   a_k(ν)/x^k`, within the first term left out for `ℓ ≥ max(ν − 1/2, 1)`
//!   (DLMF 10.40(ii)).
//! - **`I`'s** (DLMF 10.40.1): `I_ν(x) = e^x (2πx)^(−1/2) Σ_{k<ℓ} (−1)^k
//!   a_k(ν)/x^k − (sin νπ/π) K_ν(x)`. DLMF bounds no remainder for `I`, but
//!   by DLMF 10.34.2 `I_ν(x) = (i/π) (K_ν(x e^(πi)) − e^(−νπi) K_ν(x))`.
//!   The expansion of `(i/π) K_ν(x e^(πi))`, `K`'s at the phase π, is the
//!   sum, and its remainder is within `2 χ(ℓ) |a_ℓ(ν)| x^(−ℓ) exp(|ν² −
//!   1/4| π/(2x))` (DLMF 10.40.11, 10.40.12, `χ(1) = π/2`), where `χ(ℓ) =
//!   √π Γ(ℓ/2 + 1)/Γ(ℓ/2 + 1/2) < √(π (ℓ/2 + 1))` by Gautschi's
//!   inequality. The real part of the rest, `(i/π) e^(−νπi) K_ν(x)`, is
//!   `(sin νπ/π) K_ν(x)`, which vanishes for an integer `ν`.
//!
//! The terms fall until `k` is about `2x`, to about `e^(−2x)`, so the
//! expansions serve `wp` bits from `x` about `wp ln 2 / 2`. Where the terms
//! grow again before they reach `2^-wp` of the sum, the expansion gives
//! `None`. For `ν² > 2x` the first terms grow, to about `e^(ν²/(2x))`,
//! before they fall: their sums take that many more bits.

use std::f64::consts::{LN_2, LOG2_E};

use super::super::{exponent, one, upper};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Rational};

/// The order `ν ≥ 0` of an expansion.
#[derive(Clone, Debug)]
pub(in crate::ball::fns) struct Order {
    nu: Rational,
    /// `4ν²` as `four / den²`: `four = 4 p²` for `ν = p/den`.
    four: Integer,
    den_square: Integer,
}

impl Order {
    /// The order `ν`.
    ///
    /// # Panics
    ///
    /// If `ν` is negative.
    pub(in crate::ball::fns) fn new(nu: Rational) -> Order {
        assert!(!nu.is_negative(), "Order::new: a negative order");
        let p = nu.numer().clone();
        let den = nu.denom().clone();
        Order {
            four: Integer::from(4) * &p * &p,
            den_square: &den * &den,
            nu,
        }
    }

    pub(in crate::ball::fns) fn nu(&self) -> &Rational {
        &self.nu
    }

    /// `ν` in `f64`, for the counts.
    fn approx(&self) -> f64 {
        Ball::from_rational(&self.nu, 64).mid().to_f64()
    }

    fn is_integer(&self) -> bool {
        self.nu.is_integer()
    }

    /// `a_(k+1)/a_k = (4ν² − (2k + 1)²)/(8(k + 1))` as an integer quotient.
    fn ratio(&self, k: u64) -> (Integer, Integer) {
        let odd = Integer::from(2 * k + 1);
        let numerator = &self.four - &(&(&odd * &odd) * &self.den_square);
        let denominator = Integer::from(8 * (k + 1)) * &self.den_square;
        (numerator, denominator)
    }
}

/// The terms `a_k(ν)/x^k` from `k = 0`, at `q` bits.
struct Terms<'o> {
    order: &'o Order,
    x: Ball,
    q: u64,
    k: u64,
    term: Ball,
}

impl<'o> Terms<'o> {
    fn new(order: &'o Order, x: &Float, q: u64) -> Terms<'o> {
        Terms {
            order,
            x: Ball::exact(x.clone()),
            q,
            k: 0,
            term: one(),
        }
    }

    /// The current term, `a_k/x^k`, and its index.
    fn current(&self) -> (u64, &Ball) {
        (self.k, &self.term)
    }

    /// Moves to the next term.
    fn advance(&mut self) {
        let (numerator, denominator) = self.order.ratio(self.k);
        let factor = Ball::from_integer(&numerator).div(&Ball::from_integer(&denominator), self.q);
        self.term = self.term.mul(&factor, self.q).div(&self.x, self.q);
        self.k += 1;
    }
}

/// `|x|` in `f64`, saturated past its range.
fn magnitude(x: &Float) -> f64 {
    if exponent(x) > 1000 {
        return f64::INFINITY;
    }
    x.to_f64().abs()
}

/// The least `ℓ ≥ 1` with `ℓ ≥ bound`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a small positive count"
)]
fn at_least(bound: f64) -> u64 {
    bound.ceil().max(1.0) as u64
}

/// The most terms an expansion takes for `x`, past which they grow:
/// about `2x` for a small order, more for a large one.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a positive count, capped"
)]
fn most_terms(order: &Order, x: f64) -> u64 {
    (4.0 * x + 4.0 * order.approx() + 16.0).min(1e6) as u64
}

/// The extra bits for the growth of the first terms, `ν²/(2x) log₂ e`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a small positive estimate"
)]
fn growth(order: &Order, x: f64) -> u64 {
    let nu = order.approx();
    (nu * nu / (2.0 * x) * LOG2_E).clamp(0.0, 1e6).ceil() as u64
}

/// `2^-bits`.
fn tiny(bits: u64) -> Mag {
    Mag::power_of_two(&-Integer::from(bits))
}

/// Hankel's `P` and `Q` at `x`, each with its remainder; `None` where the
/// terms grow before they reach `2^-wp`.
fn hankel_sums(order: &Order, x: &Float, wp: u64) -> Option<(Ball, Ball)> {
    let ax = magnitude(x);
    let q = wp + growth(order, ax) + 32;
    let small = tiny(wp + 16);
    let nu = order.approx();
    let (p_least, q_least) = (at_least(nu / 2.0 - 0.25), at_least(nu / 2.0 - 0.75));
    let mut terms = Terms::new(order, x, q);
    let (mut p, mut p_rest) = (Ball::ZERO, None);
    let (mut s, mut s_rest) = (Ball::ZERO, None);
    for _ in 0..most_terms(order, ax) {
        let (k, term) = terms.current();
        let size = upper(term);
        // P takes the even k = 2j, Q the odd k = 2j + 1, signs (−1)^j.
        let j = k / 2;
        let signed = if j % 2 == 0 { term.clone() } else { term.neg() };
        if k % 2 == 0 {
            if p_rest.is_none() {
                if j >= p_least && size <= small {
                    p_rest = Some(size);
                } else {
                    p = p.add(&signed, q);
                }
            }
        } else if s_rest.is_none() {
            if j >= q_least && size <= small {
                s_rest = Some(size);
            } else {
                s = s.add(&signed, q);
            }
        }
        if let (Some(p_rest), Some(s_rest)) = (&p_rest, &s_rest) {
            return Some((p.add_error(p_rest), s.add_error(s_rest)));
        }
        terms.advance();
    }
    None
}

/// `(J_ν(x), Y_ν(x))` by Hankel's expansion, at about `wp` bits; `None`
/// where its terms grow before they reach them.
pub(in crate::ball::fns) fn hankel(order: &Order, x: &Float, wp: u64) -> Option<(Ball, Ball)> {
    let (p, s) = hankel_sums(order, x, wp)?;
    let q = wp + 32;
    let xb = Ball::exact(x.clone());
    let pi = Ball::const_pi(q);
    // ω = x − θ, θ = (ν/2 + 1/4) π: cos ω and sin ω from those of x and θ,
    // so that a large x is reduced once, exactly.
    let quarter = Rational::new(Integer::ONE, Integer::from(4));
    let theta = Ball::from_rational(
        &(&(order.nu() / &Rational::from(Integer::from(2))) + &quarter),
        q,
    )
    .mul(&pi, q);
    let (sin_x, cos_x) = xb.sin_cos(q);
    let (sin_t, cos_t) = theta.sin_cos(q);
    let cos_w = cos_x.mul(&cos_t, q).add(&sin_x.mul(&sin_t, q), q);
    let sin_w = sin_x.mul(&cos_t, q).sub(&cos_x.mul(&sin_t, q), q);
    let scale = Ball::from_integer(&Integer::from(2))
        .div(&pi.mul(&xb, q), q)
        .sqrt(q);
    let j = scale.mul(&p.mul(&cos_w, q).sub(&s.mul(&sin_w, q), q), q);
    let y = scale.mul(&p.mul(&sin_w, q).add(&s.mul(&cos_w, q), q), q);
    Some((j, y))
}

/// `K`'s sum `Σ_{k<ℓ} a_k/x^k` with its remainder; `None` where its terms
/// grow before they reach `2^-wp`.
fn k_sum(order: &Order, x: &Float, wp: u64, q: u64) -> Option<Ball> {
    let ax = magnitude(x);
    let small = tiny(wp + 16);
    let least = at_least(order.approx() - 0.5);
    let mut terms = Terms::new(order, x, q);
    let mut sum = Ball::ZERO;
    for _ in 0..most_terms(order, ax) {
        let (k, term) = terms.current();
        let size = upper(term);
        if k >= least && size <= small {
            return Some(sum.add_error(&size));
        }
        sum = sum.add(term, q);
        terms.advance();
    }
    None
}

/// `K_ν(x)` by its expansion, at about `wp` bits; `None` where its terms
/// grow before they reach them.
pub(in crate::ball::fns) fn k_expansion(order: &Order, x: &Float, wp: u64) -> Option<Ball> {
    let q = wp + growth(order, magnitude(x)) + 32;
    let sum = k_sum(order, x, wp, q)?;
    let xb = Ball::exact(x.clone());
    // √(π/(2x)) e^(−x).
    let scale = Ball::const_pi(q)
        .div(&xb.mul_2exp(&Integer::ONE), q)
        .sqrt(q)
        .mul(&xb.neg().exp(q), q);
    Some(scale.mul(&sum, q))
}

/// `2 χ(ℓ) exp(|ν² − 1/4| π/(2x))` from above, with
/// `χ(ℓ) < √(π (ℓ/2 + 1)) ≤ ℓ + 2`.
fn i_factor(order: &Order, x: &Float, l: u64) -> Mag {
    let chi = Mag::from_integer(&Integer::from(l + 2));
    // |ν² − 1/4| π/(2x) ≤ 2 |ν² − 1/4| / x.
    let nu = order.nu();
    let offset = &(nu * nu) - &Rational::new(Integer::ONE, Integer::from(4));
    let exponent_bound = Mag::from_integer(&offset.numer().abs())
        .div(&Mag::from_integer_lower(offset.denom()))
        .mul_2exp(&Integer::ONE)
        .div(&Mag::from_float_lower(x));
    chi.mul_2exp(&Integer::ONE).mul(&exponent_bound.exp())
}

/// `I_ν(x)` by its expansion, at about `wp` bits; `None` where its terms
/// grow before they reach them.
pub(in crate::ball::fns) fn i_expansion(order: &Order, x: &Float, wp: u64) -> Option<Ball> {
    let ax = magnitude(x);
    let q = wp + growth(order, ax) + 32;
    let small = tiny(wp + 16);
    let mut terms = Terms::new(order, x, q);
    let mut sum = Ball::ZERO;
    let mut found = None;
    for _ in 0..most_terms(order, ax) {
        let (k, term) = terms.current();
        if k >= 1 {
            let bound = upper(term).mul(&i_factor(order, x, k));
            if bound <= small {
                found = Some(sum.add_error(&bound));
                break;
            }
        }
        let signed = if k % 2 == 0 { term.clone() } else { term.neg() };
        sum = sum.add(&signed, q);
        terms.advance();
    }
    let sum = found?;
    let xb = Ball::exact(x.clone());
    // e^x (2πx)^(−1/2).
    let scale = xb.exp(q).div(
        &Ball::const_pi(q)
            .mul(&xb, q)
            .mul_2exp(&Integer::ONE)
            .sqrt(q),
        q,
    );
    let value = scale.mul(&sum, q);
    if order.is_integer() {
        return Some(value);
    }
    // − (sin νπ/π) K_ν(x).
    let k = k_expansion(order, x, wp)?;
    let pi = Ball::const_pi(q);
    let sine = Ball::from_rational(order.nu(), q).mul(&pi, q).sin(q);
    Some(value.sub(&sine.div(&pi, q).mul(&k, q), q))
}

/// The arguments from which the expansions serve `wp` bits: their least
/// term, about `e^(−2x)`, is then below `2^-(wp + 16)`.
#[allow(clippy::cast_precision_loss, reason = "a threshold")]
pub(in crate::ball::fns) fn threshold(wp: u64) -> f64 {
    (wp as f64 + 24.0) * LN_2 / 2.0
}

#[cfg(test)]
mod tests;
