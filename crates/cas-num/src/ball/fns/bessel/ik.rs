//! `I_n` and `K_n` of integer order (CAS-DOC-34 §2.2, CAS-379), at an
//! exact `x`:
//!
//! - **`I_n`**: the series `(x/2)^n Σ_k (x²/4)^k / (k! (n + k)!)` below
//!   [`threshold`], whose terms have one sign: once `(k + 1)(n + k + 1) ≥
//!   x²/2` each is at most half the one before, and what follows a term is
//!   at most it. From [`threshold`] on, the expansion ([`i_expansion`]).
//!   `I_n(−x) = (−1)^n I_n(x)`.
//! - **`K_0` and `K_1`** by DLMF 10.31.1 below [`threshold`], with
//!   `ψ(k + 1) = H_k − γ`: `K_0(x) = −(ln(x/2) + γ) I_0(x) + Σ_{k≥1} H_k
//!   (x²/4)^k / (k!)²` and `K_1(x) = 1/x + (ln(x/2) + γ) I_1(x) − (x/4)
//!   Σ_{k≥0} (H_k + H_(k+1)) (x²/4)^k / (k! (k + 1)!)`. `K` is about
//!   `e^(−x)` where the terms reach `e^x`: `2x log₂ e` bits more. From
//!   [`threshold`] on, the expansion ([`k_expansion`]).
//! - **`K_n`** by the forward recurrence `K_(k+1) = (2k/x) K_k + K_(k−1)`,
//!   whose terms have one sign. `K_n` takes `x > 0`.

use std::f64::consts::LOG2_E;

use super::super::{one, upper};
use super::expansions::{Order, i_expansion, k_expansion, threshold};
use super::jy::{forward, growth, lift, magnitude, order_bits};
use crate::ball::{Ball, Float, Mag};
use crate::{Integer, Natural, Rational};

fn small_integer(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

/// The series of `I_n(x)` for an exact `x ≥ 0`, at about `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn i_series(n: u64, x: &Float, wp: u64) -> Ball {
    let q = wp + 32;
    let ax = magnitude(x);
    let half = Ball::exact(x.mul_2exp(&-Integer::ONE));
    let u = half.mul(&half, q);
    let factorial = Ball::from_integer(&Integer::from(Natural::factorial(n)));
    let n_i64 = i64::try_from(n).expect("a small order");
    let mut term = half.pow(n_i64, q).div(&factorial, q);
    let mut sum = term.clone();
    let square_half = ax * ax / 2.0;
    let small = Mag::power_of_two(&-Integer::from(q));
    let mut k = 0u64;
    let rest = loop {
        k += 1;
        term = term.mul(&u, q).div(&small_integer(k * (n + k)), q);
        sum = sum.add(&term, q);
        let size = upper(&term);
        let halving = ((k + 1) * (n + k + 1)) as f64 >= square_half;
        if halving && size <= upper(&sum).mul(&small) {
            break size;
        }
    };
    sum.add_error(&rest)
}

/// `I_n(x)` for an exact `x`, at about `wp` bits.
pub(super) fn i_exact(n: u64, x: &Float, wp: u64) -> Ball {
    if x.is_zero() {
        return if n == 0 { one() } else { Ball::ZERO };
    }
    let ax = x.abs();
    let order = Order::new(Rational::from(Integer::from(n)));
    let value = if magnitude(x) >= threshold(wp) {
        i_expansion(&order, &Ball::exact(ax.clone()), wp).unwrap_or_else(|| i_series(n, &ax, wp))
    } else {
        i_series(n, &ax, wp)
    };
    if x.is_negative() && n % 2 == 1 {
        value.neg()
    } else {
        value
    }
}

/// `K_0(x)` and `K_1(x)` for an exact `x > 0` by their series, at about
/// `wp` bits.
#[allow(clippy::cast_precision_loss, reason = "term counts")]
fn k_series(x: &Float, wp: u64) -> (Ball, Ball) {
    let ax = magnitude(x);
    let (most, first) = growth(0, ax);
    // The terms reach 2^(first + most); K is about e^(−x).
    let q = wp + lift(first + most + ax * LOG2_E + ax.max(2.0).log2().log2()) + 40;
    let i0 = i_series(0, x, q);
    let i1 = i_series(1, x, q);
    let xb = Ball::exact(x.clone());
    let half = Ball::exact(x.mul_2exp(&-Integer::ONE));
    let quarter = half.mul(&half, q);
    // c_k = (x²/4)^k / (k!)², e_k = (x²/4)^k / (k! (k + 1)!), and H_k.
    let mut c = one();
    let mut e = one();
    let mut harmonic = Ball::ZERO;
    let mut sum0 = Ball::ZERO;
    // k = 0: (H_0 + H_1) e_0 = 1.
    let mut sum1 = one();
    let mut largest = Mag::from_integer(&Integer::ONE);
    let square_half = ax * ax / 2.0;
    let mut k = 0u64;
    let rest = loop {
        k += 1;
        harmonic = harmonic.add(&one().div(&small_integer(k), q), q);
        let next_harmonic = harmonic.add(&one().div(&small_integer(k + 1), q), q);
        c = c.mul(&quarter, q).div(&small_integer(k * k), q);
        e = e.mul(&quarter, q).div(&small_integer(k * (k + 1)), q);
        let t0 = harmonic.mul(&c, q);
        let t1 = harmonic.add(&next_harmonic, q).mul(&e, q);
        sum0 = sum0.add(&t0, q);
        sum1 = sum1.add(&t1, q);
        let size = upper(&t0).add(&upper(&t1));
        if size > largest {
            largest = size.clone();
        }
        // As for Y's series: from (k + 1)² ≥ x²/2 each next term is at most
        // 3/4 of this one, and what follows at most 3 times it.
        let halving = ((k + 1) * (k + 1)) as f64 >= square_half;
        if halving && size <= largest.mul(&Mag::power_of_two(&-Integer::from(q))) {
            break size.mul(&Mag::from_integer(&Integer::from(3)));
        }
    };
    let log = xb
        .mul_2exp(&-Integer::ONE)
        .log(q)
        .add(&Ball::const_euler(q), q);
    let k0 = sum0.add_error(&rest).sub(&log.mul(&i0, q), q);
    let k1 = one().div(&xb, q).add(&log.mul(&i1, q), q).sub(
        &xb.mul_2exp(&-Integer::from(2))
            .mul(&sum1.add_error(&rest), q),
        q,
    );
    (k0, k1)
}

/// `K_n(x)` for an exact `x > 0`, at about `wp` bits.
pub(super) fn k_exact(n: u64, x: &Float, wp: u64) -> Ball {
    let q = wp + order_bits(n) + 16;
    let large = magnitude(x) >= threshold(q);
    let expanded = large
        .then(|| {
            let zero = Order::new(Rational::from(Integer::ZERO));
            let one = Order::new(Rational::from(Integer::ONE));
            let xb = Ball::exact(x.clone());
            Some((k_expansion(&zero, &xb, q)?, k_expansion(&one, &xb, q)?))
        })
        .flatten();
    let (k0, k1) = expanded.unwrap_or_else(|| k_series(x, q));
    forward(k0, k1, n, x, q, true)
}
