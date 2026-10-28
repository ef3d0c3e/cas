//! Gauss–Legendre rules (CAS-DOC-26 §2, CAS-322): the nodes `x_{n,k}`,
//! the zeros of the Legendre polynomial `P_n`, and the weights
//! `w_{n,k} = 2/((1 − x²) P_n′(x)²)`, each in a ball that contains it, as
//! Arb's `arb_hypgeom_legendre_p_ui_root` gives them \[Joh17\]:
//!
//! - `P_n` and `P_n′` by the three-term recurrence
//!   `(m + 1) P_{m+1} = (2m + 1) x P_m − m P_{m−1}`, stable on `[−1, 1]`,
//!   and `P_n′ = n (x P_n − P_{n−1})/(x² − 1)`. In ball arithmetic, which
//!   sees no correlation between `P_m` and `P_{m−1}`, a radius grows by
//!   `1 + √2` at each step, `2^(1.27n)` in all, for the rounding errors and
//!   for the radius of `X` alike: the work runs at [`EXTRA`]`·n` bits more,
//!   `2.6n`, so that `P_n′(X)` stays away from 0;
//! - Newton's iteration on exact midpoints, from Tricomi's approximation
//!   `(1 − (n − 1)/(8n³)) cos(π(4k − 1)/(4n + 2))` of the `k`-th largest
//!   zero, at 64 bits and then twice as many each round;
//! - the interval Newton test \[Moo66\]: every zero of `P_n` in `X = x̃ ± r`
//!   is in `N = x̃ − P_n(x̃)/P_n′(X)`, so that `N ⊆ X` proves a zero in `N`,
//!   the only one in `X` when `0 ∉ P_n′(X)`. The radius starts at twice
//!   the Newton step's ball `P_n(x̃)/P_n′(x̃)`, and grows by `2²` up to
//!   [`ATTEMPTS`] times, then the precision doubles;
//! - the zeros in `(0, 1)` only, each apart from the next, `⌊n/2⌋` of them,
//!   as many as `P_n` has: their negatives are the others, and 0 is one for
//!   an odd `n`.
//!
//! The rules are cached process-wide per `(n, prec)`, at most
//! [`MAX_CACHED`] of them, cleared with the constants
//! ([`super::clear_constant_cache`]).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the method: P_n at x, its zeros k, the slope d, the weights w"
)]

use std::collections::HashMap;
use std::sync::{Arc, LazyLock, PoisonError, RwLock};

use super::{Ball, Float, Mag};
use crate::{Integer, Rational};

/// A rule: the nodes in descending order, each with its weight.
pub type Rule = Arc<[(Ball, Ball)]>;

/// The rules kept at most; past it the cache starts afresh.
pub const MAX_CACHED: usize = 64;

/// The widenings of the radius at one precision before it doubles.
const ATTEMPTS: u32 = 8;

/// The extra bits per point of the rule, in tenths: `2.6n` bits pay for the
/// radii's growth through the recurrence, `(1 + √2)^n` twice, at a node
/// and over its ball.
const EXTRA: u64 = 26;

/// The precisions a zero is tried at, doubling: past them a rule is
/// computed at twice the precision, which then always succeeds in practice.
const ROUNDS: u32 = 6;

static RULES: LazyLock<RwLock<HashMap<(usize, u64), Rule>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

/// Frees the cached rules.
pub(super) fn clear_rules() {
    RULES
        .write()
        .unwrap_or_else(PoisonError::into_inner)
        .clear();
}

/// The `n`-point Gauss–Legendre rule on `[−1, 1]` at `prec` bits: the nodes
/// in descending order with their weights, balls that contain them, of
/// radii about `2^(−prec)`. `∫_{−1}^{1} f ≈ Σ w_k f(x_k)`.
///
/// # Panics
///
/// If `n` is 0 or `prec` is below 2.
#[must_use]
pub fn gauss_legendre(n: usize, prec: u64) -> Rule {
    assert!(n >= 1 && prec >= 2, "a rule of {n} points at {prec} bits");
    let key = (n, prec);
    if let Some(rule) = RULES
        .read()
        .unwrap_or_else(PoisonError::into_inner)
        .get(&key)
    {
        return rule.clone();
    }
    let rule = compute(n, prec);
    let mut rules = RULES.write().unwrap_or_else(PoisonError::into_inner);
    if rules.len() >= MAX_CACHED {
        rules.clear();
    }
    rules.insert(key, rule.clone());
    rule
}

fn integer(n: usize) -> Ball {
    Ball::from_integer(&Integer::from(i64::try_from(n).expect("a count in memory")))
}

/// `(P_n(x), P_n′(x))` at `prec` bits, `n ≥ 1`.
fn legendre(n: usize, x: &Ball, prec: u64) -> (Ball, Ball) {
    let one = integer(1);
    let (mut previous, mut current) = (one.clone(), x.clone());
    for m in 1..n {
        let a = integer(2 * m + 1).mul(x, prec).mul(&current, prec);
        let b = integer(m).mul(&previous, prec);
        let next = a.sub(&b, prec).div(&integer(m + 1), prec);
        previous = current;
        current = next;
    }
    let numerator = x.mul(&current, prec).sub(&previous, prec);
    let denominator = x.mul(x, prec).sub(&one, prec);
    let derivative = integer(n).mul(&numerator, prec).div(&denominator, prec);
    (current, derivative)
}

/// Tricomi's approximation of the `k`-th largest zero of `P_n`, `1 ≤ k`.
fn initial(n: usize, k: usize) -> Float {
    #[allow(
        clippy::cast_precision_loss,
        reason = "an initial approximation, which Newton's iteration refines"
    )]
    let (n, k) = (n as f64, k as f64);
    let theta = std::f64::consts::PI * (4.0 * k - 1.0) / (4.0 * n + 2.0);
    let x = (1.0 - (n - 1.0) / (8.0 * n * n * n)) * theta.cos();
    Float::from_f64(x).expect("a finite approximation")
}

/// The working precision for `n` points at `prec` bits (the module's
/// [`EXTRA`] bits).
fn working(n: usize, prec: u64) -> u64 {
    let n = u64::try_from(n).expect("a count in memory");
    prec + 32 + (EXTRA * n).div_ceil(10)
}

/// An upper bound of `|b|`.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(&b.mid().abs()).add(b.rad())
}

/// The `k`-th largest zero of `P_n` in a ball of radius about
/// `2^(−prec)`, `k ≤ n/2`; `None` when the test fails at every precision
/// tried.
fn zero(n: usize, k: usize, prec: u64) -> Option<Ball> {
    let mut x = initial(n, k);
    let mut target = working(n, prec);
    for _ in 0..ROUNDS {
        let mut wp = 64;
        loop {
            for _ in 0..3 {
                let point = Ball::exact(x.clone());
                let (p, d) = legendre(n, &point, wp);
                x = point.sub(&p.div(&d, wp), wp).mid().clone();
            }
            if wp >= target {
                break;
            }
            wp = (2 * wp).min(target);
        }
        let point = Ball::exact(x.clone());
        let (p, d) = legendre(n, &point, target);
        let floor = Mag::power_of_two(&-Integer::from(i64::try_from(target).ok()?));
        let mut radius = upper(&p.div(&d, target))
            .mul_2exp(&Integer::from(1))
            .add(&floor);
        for _ in 0..ATTEMPTS {
            let around = Ball::new(x.clone(), radius.clone());
            let (_, slope) = legendre(n, &around, target);
            if !slope.contains_rational(&Rational::ZERO) {
                let newton = point.sub(&p.div(&slope, target), target);
                if around.contains_ball(&newton) {
                    return Some(newton);
                }
            }
            radius = radius.mul_2exp(&Integer::from(2));
        }
        target *= 2;
    }
    None
}

fn compute(n: usize, prec: u64) -> Rule {
    let wp = working(n, prec);
    let mut positive: Vec<Ball> = Vec::with_capacity(n / 2);
    for k in 1..=n / 2 {
        let x = zero(n, k, prec)
            .unwrap_or_else(|| panic!("the zero {k} of P_{n} at {prec} bits is not validated"));
        assert!(
            x.sign() == Some(std::cmp::Ordering::Greater)
                && positive.last().is_none_or(|last| !last.overlaps(&x)),
            "the zero {k} of P_{n} at {prec} bits is not apart"
        );
        positive.push(x);
    }
    let two = integer(2);
    let one = integer(1);
    let weight = |x: &Ball| {
        let (_, d) = legendre(n, x, wp);
        let w = two.div(&one.sub(&x.mul(x, wp), wp).mul(&d.mul(&d, wp), wp), wp);
        w.round_to_accuracy(prec)
    };
    let mut rule: Vec<(Ball, Ball)> = Vec::with_capacity(n);
    for x in &positive {
        rule.push((x.round_to_accuracy(prec), weight(x)));
    }
    if n % 2 == 1 {
        let origin = Ball::ZERO;
        let w = weight(&origin);
        rule.push((origin, w));
    }
    for (x, w) in positive.iter().rev().map(|x| (x.neg(), weight(x))) {
        rule.push((x.round_to_accuracy(prec), w));
    }
    rule.into()
}

#[cfg(test)]
mod tests;
