//! π, log 2, e, Euler's γ and Catalan's G by binary splitting, with a global
//! cache (CAS-DOC-5 §7.1).
//!
//! Each constant is computed at `prec + GUARD` bits by [`bsplit`] sums and a
//! few ball operations, with the series tails bounded rigorously and added
//! to the radius, so the result is a ball like any other: correct whatever
//! the term count, and tight because the count is chosen for the precision.
//! γ is in [`euler`].
//!
//! The cache keeps, per constant, the ball at the highest precision computed
//! so far. A request at or below it rounds the cached ball (half an ulp more
//! radius); a request above it recomputes at `max(prec, 2·cached)`, so a
//! sequence of growing requests costs a geometric series. It memoises a pure
//! function, so callers stay pure; the memory it holds is bounded by the
//! largest precision ever requested, and [`clear_constant_cache`] frees it.
#![allow(
    clippy::many_single_char_names,
    reason = "the series' p, q, a, b, t, d and the precisions w, n of the formulas"
)]

use std::sync::{PoisonError, RwLock};

use super::bsplit::{self, Series};
use super::{Ball, Mag};
use crate::Integer;

mod euler;
#[cfg(test)]
mod tests;

/// Bits computed beyond the requested precision. The ball operations after
/// the sum cost a few ulps each, well inside this.
const GUARD: u64 = 32;

/// Chudnovsky's `1/π = 12 Σ (-1)^k (6k)! (A + B k) / ((3k)! (k!)^3 C^(3k+3/2))`:
/// `p(k) = -(6k-5)(2k-1)(6k-1)`, `q(k) = k^3 C^3 / 24`, `a(k) = A + B k`
/// (`p(0) = q(0) = 1`), so that `S = Σ a(k) Π p/q` and `π = 426880 √10005 / S`
/// (`426880 √10005 = C^(3/2) / 12`).
struct Chudnovsky;

const CHUD_A: u64 = 13_591_409;
const CHUD_B: u64 = 545_140_134;
/// `C^3` for `C = 640320`.
const CHUD_C3: u64 = 262_537_412_640_768_000;
/// About 47.11 bits per term: `log2(C^3 / 1728)`, rounded down.
const CHUD_BITS_PER_TERM: u64 = 47;

impl Series for Chudnovsky {
    fn p(&self, k: u64) -> Option<Integer> {
        if k == 0 {
            return Some(Integer::ONE);
        }
        // 72 k^3 fits an i128 below k = 2^40, far beyond any precision.
        let k = i128::from(k);
        let (x, y, z) = (6 * k - 5, 2 * k - 1, 6 * k - 1);
        Some(match x.checked_mul(y).and_then(|xy| xy.checked_mul(z)) {
            Some(p) => Integer::from(-p),
            None => -(&(&Integer::from(x) * &Integer::from(y)) * &Integer::from(z)),
        })
    }

    fn q(&self, k: u64) -> Integer {
        if k == 0 {
            return Integer::ONE;
        }
        let c = u128::from(CHUD_C3 / 24);
        match u128::from(k)
            .checked_pow(3)
            .and_then(|k3| k3.checked_mul(c))
        {
            Some(q) => Integer::from(q),
            None => &Integer::from(k).pow(3) * &Integer::from(c),
        }
    }

    fn a(&self, k: u64) -> Option<Integer> {
        Some(Integer::from(
            u128::from(CHUD_A) + u128::from(CHUD_B) * u128::from(k),
        ))
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }
}

/// `Σ_{k >= n} |t_k| <= 2 a(n) (1728 / C^3)^n`: `|t_k| = a(k) (6k)! / ((3k)!
/// (k!)^3 C^(3k))`, and `(6k)! / ((3k)! (k!)^3) = C(6k, 3k) · (3k)! / (k!)^3
/// <= 2^(6k) 3^(3k) = 1728^k` (a binomial and a multinomial coefficient,
/// each at most the sum of its row). With `r = 1728 / C^3 < 2^-47` and
/// `a(n + j) <= (1 + j) a(n)`, the tail is at most `a(n) r^n / (1 - r)^2`.
fn chudnovsky_tail(n: u64) -> Mag {
    let r = Mag::from_integer(&Integer::from(1728u64))
        .div(&Mag::from_integer_lower(&Integer::from(CHUD_C3)));
    let a = Chudnovsky.a(n).expect("a(k) is never 1");
    r.pow(n).mul(&Mag::from_integer(&a)).mul_2exp(&Integer::ONE)
}

fn pi_uncached(prec: u64) -> Ball {
    let w = prec + GUARD;
    let n = w / CHUD_BITS_PER_TERM + 2;
    let (t, d, _) = bsplit::sum(&Chudnovsky, n);
    let s = bsplit::quotient(&t, &d, w).add_error(&chudnovsky_tail(n));
    let c = Ball::from_integer(&Integer::from(10_005u64)).sqrt(w);
    let c = c.mul(&Ball::from_integer(&Integer::from(426_880u64)), w);
    c.div(&s, w)
}

/// `atanh(1/m) = Σ 1 / ((2k+1) m^(2k+1))`: `q(0) = m`, `q(k) = m^2`,
/// `b(k) = 2k + 1`, `p = a = 1`.
struct Atanh(u64);

impl Series for Atanh {
    fn p(&self, _: u64) -> Option<Integer> {
        None
    }

    fn q(&self, k: u64) -> Integer {
        let m = u128::from(self.0);
        Integer::from(if k == 0 { m } else { m * m })
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, k: u64) -> Option<Integer> {
        Some(Integer::from(2 * u128::from(k) + 1))
    }
}

/// `atanh(1/m)` at `w` bits from `n` terms: the tail is at most
/// `m^-(2n+1) / (1 - m^-2) <= 2 m^-(2n+1)`.
fn atanh_inv(m: u64, w: u64) -> Ball {
    // 2 log2(m) bits per term (the count only sets the accuracy).
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "a term count estimate, positive and far below 2^53"
    )]
    let n = (w as f64 / (2.0 * (m as f64).log2())).ceil() as u64 + 1;
    let (t, d, _) = bsplit::sum(&Atanh(m), n);
    bsplit::quotient(&t, &d, w).add_error(&atanh_tail(m, n))
}

/// `2 m^-(2n+1)`, see [`atanh_inv`].
fn atanh_tail(m: u64, n: u64) -> Mag {
    let inv_m = Mag::ONE.div(&Mag::from_integer_lower(&Integer::from(m)));
    inv_m.pow(2 * n + 1).mul_2exp(&Integer::ONE)
}

/// `log 2 = 18 atanh(1/26) - 2 atanh(1/4801) + 8 atanh(1/8749)`.
fn log2_uncached(prec: u64) -> Ball {
    let w = prec + GUARD;
    let term = |c: i64, m: u64| atanh_inv(m, w).mul(&Ball::from_integer(&Integer::from(c)), w);
    term(18, 26).add(&term(-2, 4801), w).add(&term(8, 8749), w)
}

/// `e = Σ 1/k!`: `q(0) = 1`, `q(k) = k`, `p = a = b = 1`.
struct E;

impl Series for E {
    fn p(&self, _: u64) -> Option<Integer> {
        None
    }

    fn q(&self, k: u64) -> Integer {
        Integer::from(k.max(1))
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }
}

/// `e` from `n` terms: `Q = (n-1)!` and the tail is at most
/// `(1/n!) (1 + 1/(n+1) + 1/(n+1)^2 + …) <= 2 / n! = 2 / (n Q)`.
fn e_uncached(prec: u64) -> Ball {
    let w = prec + GUARD;
    // The smallest n with log2(n!) > w + 1 (the count only sets the accuracy).
    #[allow(
        clippy::cast_precision_loss,
        reason = "a term count estimate, far below 2^53"
    )]
    let n = {
        let (mut n, mut bits) = (1u64, 0.0f64);
        while bits <= (w + 1) as f64 {
            n += 1;
            bits += (n as f64).log2();
        }
        n
    };
    let (t, d, q) = bsplit::sum(&E, n);
    bsplit::quotient(&t, &d, w).add_error(&e_tail(n, &q))
}

/// `2 / (n Q)` with `Q = (n-1)!`, see [`e_uncached`].
fn e_tail(n: u64, q: &Integer) -> Mag {
    let nf = Mag::from_integer_lower(q).mul_lower(&Mag::from_integer_lower(&Integer::from(n)));
    Mag::ONE.mul_2exp(&Integer::ONE).div(&nf)
}

/// `2G = Σ (-8)^k (3k + 2) / ((2k+1)^3 C(2k, k)^3)` (Lima, "A rapidly
/// converging Ramanujan-type series for Catalan's constant", 2012, arXiv
/// 1207.3139, Theorem 1; Guillera's series shifted by one index). The ratio
/// of consecutive `h_k = (-8)^k / ((2k+1)^3 C(2k, k)^3)` is `-k^3 / (2k+1)^3`:
/// `p(k) = -k^3`, `q(k) = (2k+1)^3`, `a(k) = 3k + 2` (`p(0) = q(0) = 1`).
/// Three bits per term, and no π or logarithm, unlike Ramanujan's
/// `G = π/8 log(2 + √3) + 3/8 Σ 1 / ((2k+1)^2 C(2k, k))` (two bits).
struct Catalan;

impl Series for Catalan {
    fn p(&self, k: u64) -> Option<Integer> {
        let k = u128::from(k);
        Some(if k == 0 {
            Integer::ONE
        } else {
            -Integer::from(k * k * k)
        })
    }

    fn q(&self, k: u64) -> Integer {
        // (2k+1)^3 fits a u128 below k = 2^41, far beyond any precision.
        let x = 2 * u128::from(k) + 1;
        Integer::from(x * x * x)
    }

    fn a(&self, k: u64) -> Option<Integer> {
        Some(Integer::from(3 * u128::from(k) + 2))
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }
}

/// `Σ_{k >= n} |t_k| <= 5/4 (3n + 2) 8^-n` for `n >= 1`: `|h_k| <= 8^-k`, as
/// `C(2k, k) >= 4^k / (2k + 1)` (the largest of the `2k + 1` binomials that
/// sum to `4^k`), and from `k = n` on the ratio of `|t_k| = (3k + 2) |h_k|`
/// is at most `(3n + 5) / (3n + 2) / 8 <= 1/5`.
fn catalan_tail(n: u64) -> Mag {
    let first = Integer::from(5 * (3 * u128::from(n) + 2));
    Mag::from_integer(&first).mul_2exp(&-Integer::from(3 * u128::from(n) + 2))
}

fn catalan_uncached(prec: u64) -> Ball {
    let w = prec + GUARD;
    // 3n >= w + bits(w) + 8 puts the tail below 2^-w, its factor 3n + 2
    // included (the count only sets the accuracy).
    catalan_terms((w + u64::from(w.ilog2()) + 6).div_ceil(3) + 1, w)
}

/// G at `w` bits from `n >= 1` terms and the tail bound.
fn catalan_terms(n: u64, w: u64) -> Ball {
    let (t, d, _) = bsplit::sum(&Catalan, n);
    bsplit::quotient(&t, &d, w)
        .add_error(&catalan_tail(n))
        .mul_2exp(&-Integer::ONE)
}

/// A constant's cached ball and the precision it was computed for (0: none).
struct Cached {
    prec: u64,
    value: Ball,
}

type Cache = RwLock<Cached>;

const fn empty() -> Cached {
    Cached {
        prec: 0,
        value: Ball::ZERO,
    }
}

static PI: Cache = RwLock::new(empty());
static LOG2: Cache = RwLock::new(empty());
static E_CACHE: Cache = RwLock::new(empty());
static EULER: Cache = RwLock::new(empty());
static CATALAN: Cache = RwLock::new(empty());

/// The constant at `prec` bits through `cache`; see the module comment.
fn cached(cache: &Cache, prec: u64, compute: fn(u64) -> Ball) -> Ball {
    assert!(prec >= 1, "constant at precision 0");
    let have = {
        let c = cache.read().unwrap_or_else(PoisonError::into_inner);
        if c.prec >= prec {
            return c.value.round_to_accuracy(prec);
        }
        c.prec
    };
    // Computed without the lock held: readers at lower precisions are not
    // kept waiting, and a racing thread at most repeats the work.
    let target = prec.max(have.saturating_mul(2));
    let value = compute(target);
    let out = value.round_to_accuracy(prec);
    let mut c = cache.write().unwrap_or_else(PoisonError::into_inner);
    if c.prec < target {
        *c = Cached {
            prec: target,
            value,
        };
    }
    out
}

impl Ball {
    /// π at `prec` bits: the midpoint within about half an ulp, from the
    /// cache when it holds `prec` bits or more.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn const_pi(prec: u64) -> Ball {
        cached(&PI, prec, pi_uncached)
    }

    /// `log 2` at `prec` bits, as [`Ball::const_pi`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn const_log2(prec: u64) -> Ball {
        cached(&LOG2, prec, log2_uncached)
    }

    /// `e` at `prec` bits, as [`Ball::const_pi`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn const_e(prec: u64) -> Ball {
        cached(&E_CACHE, prec, e_uncached)
    }

    /// Euler's constant `γ = 0.5772…` at `prec` bits, as [`Ball::const_pi`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn const_euler(prec: u64) -> Ball {
        cached(&EULER, prec, euler::euler_uncached)
    }

    /// Catalan's constant `G = 1 - 1/3^2 + 1/5^2 - … = 0.9159…` at `prec`
    /// bits, as [`Ball::const_pi`].
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn const_catalan(prec: u64) -> Ball {
        cached(&CATALAN, prec, catalan_uncached)
    }
}

/// Frees the cached constants and the elementary functions' tables; the
/// next request recomputes from scratch.
pub fn clear_constant_cache() {
    for cache in [&PI, &LOG2, &E_CACHE, &EULER, &CATALAN] {
        *cache.write().unwrap_or_else(PoisonError::into_inner) = empty();
    }
    super::fns::clear_tables();
}
