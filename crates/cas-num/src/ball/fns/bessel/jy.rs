//! `J_n` and `Y_n` of integer order (CAS-DOC-34 §2.1, CAS-378), at an
//! exact `x`:
//!
//! - **Series**, below [`threshold`] or for `n > |x|`: `J_n(x) = (x/2)^n
//!   Σ_k (−x²/4)^k / (k! (n + k)!)`. Once `(k + 1)(n + k + 1) ≥ x²/2` each
//!   term is at most half the one before, so what follows a term is at most
//!   it. The terms reach about `e^|x|` where `J_n` is about `|x|^(−1/2)`, so
//!   the sum takes the bits of its largest term above that more.
//! - `Y_0` and `Y_1` by DLMF 10.8.1, with `ψ(k + 1) = H_k − γ`: `Y_0(x) =
//!   (2/π) (ln(x/2) + γ) J_0(x) + (2/π) Σ_{k≥1} (−1)^(k+1) H_k (x²/4)^k /
//!   (k!)²` and `Y_1(x) = −2/(πx) + (2/π) (ln(x/2) + γ) J_1(x) − (1/π)
//!   Σ_{k≥0} (H_k + H_(k+1)) (−1)^k (x/2)^(2k+1) / (k! (k + 1)!)`.
//! - **Hankel's expansion** ([`hankel`]) from [`threshold`] on, for `J_0`,
//!   `J_1`, `Y_0` and `Y_1`.
//! - **The forward recurrence** `C_(k+1) = (2k/x) C_k − C_(k−1)` for the
//!   other orders: of `Y` always, as `Y` is the dominant solution, and of
//!   `J` while `k < |x|`, where it is stable too.
//!
//! `J_n(−x) = (−1)^n J_n(x)`; `Y_n` takes `x > 0`.

use super::super::{exponent, one, upper};
use super::expansions::{Order, hankel, threshold};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Natural};

/// `|x|` in `f64`, saturated past its range.
pub(super) fn magnitude(x: &Float) -> f64 {
    if exponent(x) > 1000 {
        return f64::INFINITY;
    }
    x.to_f64().abs()
}

fn small_integer(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

/// The bits of an order, for the recurrences' rounding.
pub(super) fn order_bits(n: u64) -> u64 {
    u64::from(u64::BITS - n.leading_zeros())
}

/// Extra bits for a sum whose terms reach `2^bits` above its result.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a small positive estimate"
)]
pub(super) fn lift(bits: f64) -> u64 {
    bits.clamp(0.0, 1e7).ceil() as u64
}

/// `log₂` of the largest term of `Σ_k (x/2)^(2k+n) / (k! (n + k)!)` over its
/// first, and that first term's `log₂`, in `f64`.
#[allow(clippy::cast_precision_loss, reason = "estimates")]
pub(super) fn growth(n: u64, x: f64) -> (f64, f64) {
    let half = (x / 2.0).log2();
    // log₂ of the first term, (x/2)^n / n!.
    let mut first = n as f64 * half;
    for j in 2..=n {
        first -= (j as f64).log2();
    }
    let (mut now, mut most) = (0.0f64, 0.0f64);
    let mut k = 0u64;
    loop {
        let step = 2.0 * half - ((k + 1) as f64).log2() - ((n + k + 1) as f64).log2();
        if step <= 0.0 {
            break;
        }
        now += step;
        most = most.max(now);
        k += 1;
    }
    (most, first)
}

/// The series of `J_n(x)` for an exact `x ≥ 0`, at about `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn j_series(n: u64, x: &Float, wp: u64) -> Ball {
    let ax = magnitude(x);
    let (most, first) = growth(n, ax);
    // The result is at least about min(first term, |x|^(−1/2)).
    let floor = first.min(-0.5 * ax.max(1.0).log2() - 1.0);
    let q = wp + lift(first + most - floor) + 32;
    let half = Ball::exact(x.mul_2exp(&-Integer::ONE));
    let u = half.mul(&half, q).neg();
    let factorial = Ball::from_integer(&Integer::from(Natural::factorial(n)));
    let n_i64 = i64::try_from(n).expect("a small order");
    let mut term = half.pow(n_i64, q).div(&factorial, q);
    let mut sum = term.clone();
    let mut largest = upper(&term);
    let square_half = ax * ax / 2.0;
    let mut k = 0u64;
    let rest = loop {
        k += 1;
        term = term.mul(&u, q).div(&small_integer(k * (n + k)), q);
        sum = sum.add(&term, q);
        let size = upper(&term);
        if size > largest {
            largest = size.clone();
        }
        let halving = ((k + 1) * (n + k + 1)) as f64 >= square_half;
        if halving && size <= largest.mul(&Mag::power_of_two(&-Integer::from(q))) {
            break size;
        }
    };
    sum.add_error(&rest)
}

/// `Y_0(x)` and `Y_1(x)` for an exact `x > 0` by their series, at about
/// `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn y_series(x: &Float, wp: u64) -> (Ball, Ball) {
    let ax = magnitude(x);
    let (most, first) = growth(0, ax);
    let floor = first.min(-0.5 * ax.max(1.0).log2() - 1.0);
    // H_k grows as log k: a few bits more.
    let q = wp + lift(first + most - floor + ax.max(2.0).log2().log2()) + 40;
    let j0 = j_series(0, x, q);
    let j1 = j_series(1, x, q);
    let xb = Ball::exact(x.clone());
    let half = Ball::exact(x.mul_2exp(&-Integer::ONE));
    let quarter = half.mul(&half, q);
    // c_k = (x²/4)^k / (k!)² for Y_0, d_k = (−x²/4)^k (x/2) / (k! (k+1)!)
    // for Y_1, and H_k.
    let mut c = one();
    let mut d = half.clone();
    let mut harmonic = Ball::ZERO;
    let mut sum0 = Ball::ZERO;
    // k = 0: H_0 + H_1 = 1.
    let mut sum1 = d.clone();
    let mut largest = upper(&d);
    let square_half = ax * ax / 2.0;
    let mut k = 0u64;
    let rest = loop {
        k += 1;
        let kb = small_integer(k);
        harmonic = harmonic.add(&one().div(&kb, q), q);
        let next_harmonic = harmonic.add(&one().div(&small_integer(k + 1), q), q);
        c = c.mul(&quarter, q).div(&small_integer(k * k), q);
        d = d.mul(&quarter, q).neg().div(&small_integer(k * (k + 1)), q);
        // (−1)^(k+1) H_k c_k and (H_k + H_(k+1)) d_k.
        let t0 = harmonic.mul(&c, q);
        let t0 = if k % 2 == 1 { t0 } else { t0.neg() };
        let t1 = harmonic.add(&next_harmonic, q).mul(&d, q);
        sum0 = sum0.add(&t0, q);
        sum1 = sum1.add(&t1, q);
        let size = upper(&t0).add(&upper(&t1));
        if size > largest {
            largest = size.clone();
        }
        // From (k + 1)² ≥ x²/2 the terms of c and d at least halve, and
        // H_(k+2)/H_k ≤ 3/2 then (k ≥ 1): each next term is at most 3/4
        // of this one, and what follows at most 3 times it.
        let halving = ((k + 1) * (k + 1)) as f64 >= square_half;
        if halving && size <= largest.mul(&Mag::power_of_two(&-Integer::from(q))) {
            break size.mul(&Mag::from_integer(&Integer::from(3)));
        }
    };
    let pi = Ball::const_pi(q);
    let log = xb
        .mul_2exp(&-Integer::ONE)
        .log(q)
        .add(&Ball::const_euler(q), q);
    let two_over_pi = small_integer(2).div(&pi, q);
    let y0 = two_over_pi.mul(&log.mul(&j0, q).add(&sum0.add_error(&rest), q), q);
    let y1 = two_over_pi
        .mul(&log.mul(&j1, q), q)
        .sub(&two_over_pi.div(&xb, q), q)
        .sub(&sum1.add_error(&rest).div(&pi, q), q);
    (y0, y1)
}

/// `C_n` from `C_0 = c0` and `C_1 = c1` by `C_(k+1) = (2k/x) C_k ∓ C_(k−1)`,
/// `−` for `J` and `Y`, `+` for `K`.
pub(super) fn forward(c0: Ball, c1: Ball, n: u64, x: &Float, q: u64, plus: bool) -> Ball {
    if n == 0 {
        return c0;
    }
    let xb = Ball::exact(x.clone());
    let (mut before, mut now) = (c0, c1);
    for k in 1..n {
        let scaled = small_integer(2 * k).div(&xb, q).mul(&now, q);
        let next = if plus {
            scaled.add(&before, q)
        } else {
            scaled.sub(&before, q)
        };
        before = now;
        now = next;
    }
    now
}

/// `(J_0(x), J_1(x), Y_0(x), Y_1(x))` by Hankel's expansion for an exact
/// `x > 0`, at about `wp` bits; `None` below its reach.
fn hankel_pairs(x: &Float, wp: u64) -> Option<(Ball, Ball, Ball, Ball)> {
    let zero = Order::new(crate::Rational::from(Integer::ZERO));
    let one = Order::new(crate::Rational::from(Integer::ONE));
    let x = Ball::exact(x.clone());
    let (j0, y0) = hankel(&zero, &x, wp)?;
    let (j1, y1) = hankel(&one, &x, wp)?;
    Some((j0, j1, y0, y1))
}

/// `J_n(x)` for an exact `x`, at about `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "an order against an argument")]
pub(super) fn j_exact(n: u64, x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return if n == 0 { one() } else { Ball::ZERO };
    }
    let ax = x.abs();
    let size = magnitude(x);
    let q = wp + order_bits(n) + 16;
    let value = if size >= threshold(q) && n as f64 <= size {
        match hankel_pairs(&ax, q) {
            Some((j0, j1, _, _)) => forward(j0, j1, n, &ax, q, false),
            None => j_series(n, &ax, wp),
        }
    } else {
        j_series(n, &ax, wp)
    };
    if x.is_negative() && n % 2 == 1 {
        value.neg()
    } else {
        value
    }
}

/// `Y_n(x)` for an exact `x > 0`, at about `wp` bits.
pub(super) fn y_exact(n: u64, x: &Float, wp: u64) -> Ball {
    let q = wp + order_bits(n) + 16;
    let (y0, y1) = match (magnitude(x) >= threshold(q))
        .then(|| hankel_pairs(x, q))
        .flatten()
    {
        Some((_, _, y0, y1)) => (y0, y1),
        None => y_series(x, q),
    };
    forward(y0, y1, n, x, q, false)
}

/// `√(J_k(x)² + Y_k(x)²)` from above at an exact `x > 0`: the modulus,
/// which falls as `x` grows (Nicholson's formula, DLMF 10.9.30, has
/// `K_0(2x sinh t)` under the integral).
pub(super) fn modulus(k: u64, x: &Float) -> Mag {
    let (j, y) = (j_exact(k, x, 64), y_exact(k, x, 64));
    let (j, y) = (upper(&j), upper(&y));
    j.mul(&j).add(&y.mul(&y)).sqrt()
}
