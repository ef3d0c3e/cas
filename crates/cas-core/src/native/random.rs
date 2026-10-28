//! Random numbers (CAS-DOC-33): Xcas's `srand`, `rand`, `randint` and
//! the rest, drawn from the context's generator, xoshiro256\*\*
//! ([`Generator`]), through [`Eval::draw`], so that a call that draws is
//! evaluated again each time it is met (§2).
//!
//! - **The seed** (§1, CAS-370): `SeedRandom[n]` (`srand`, `randseed`)
//!   seeds the generator with the integer `n` and gives `n`;
//!   `SeedRandom[]` seeds it from the clock and gives the seed it took, so
//!   that a run can be repeated. A new context starts from the seed 0.
//! - **Uniform numbers** (§3, CAS-371): `Random[]` (`rand`, `random`), a
//!   real in `[0, 1)`: an exact ball `m 2^-wp`, `m` drawn below `2^wp` at
//!   the working precision; `Random[n]`, an integer in `[0, n)`, or `(n,
//!   0]` for `n < 0`; `Random[a, b]`, the real `a + (b − a) u` between two
//!   real numbers, in balls; `Random[k, a, b]`, `k` distinct integers
//!   between `a` and `b`, in a random order; `Random[l]`, an element of the
//!   list `l`, and `Random[k, l]`, `k` of its elements at distinct places.
//!   `RandomInteger[a, b]` (`randint`) is an integer between `a` and `b`,
//!   both included, and `RandomInteger[n]` one between 1 and `n`.
//!
//! - **Variates** of the distributions (§4, CAS-372): [`variates`].
//! - **Vectors, matrices and polynomials, permutations and samples** (§5,
//!   §6, CAS-373): [`collections`].
//!
//! Integers are drawn by rejection below the width of their range, exact
//! at any size; distinct ones by the first `k` swaps of a Fisher–Yates
//! shuffle of the range, which a map keeps. A count past what the range or
//! the list holds, or past `Limits::terms`, and an integer form given a
//! number that is no integer, leave the call with a warning, where giac
//! errs; a symbol leaves it without one.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use cas_num::ball::{Ball, Float};
use cas_num::random::Generator;
use cas_num::{Integer, Natural};
use foldhash::fast::FixedState;
use hashbrown::HashMap;

use super::numeric::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// A seed from the clock: its nanoseconds, plus a count of the seeds taken
/// so that two in one tick differ.
fn clock_seed() -> Integer {
    static TAKEN: AtomicU64 = AtomicU64::new(0);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let taken = u128::from(TAKEN.fetch_add(1, Ordering::Relaxed));
    Integer::from(nanos.wrapping_add(taken << 64))
}

/// `SeedRandom[n]`, `SeedRandom[]`.
pub(super) fn seed_random(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let seed = match e.args() {
        [] => clock_seed(),
        [n] => match n.to_integer() {
            Some(n) => n,
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    ev.draw(|g| *g = Generator::seeded(&seed));
    Ok(Some(Expr::from(seed)))
}

/// A real in `[0, 1)` at the working precision: `wp` random bits, an exact
/// ball.
pub(super) fn unit(ev: &mut Eval<'_>) -> Ball {
    let wp = ev.context().working_precision();
    let m = ev.draw(|g| g.bits(wp));
    Ball::exact(Float::from_integer_exp(
        &Integer::from(m),
        &-Integer::from(wp),
    ))
}

/// An integer between `a` and `b`, both included.
fn between(ev: &mut Eval<'_>, a: &Integer, b: &Integer) -> Integer {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let width = (hi - lo).unsigned_abs() + Natural::from(1u64);
    lo + &Integer::from(ev.draw(|g| g.below(&width)))
}

/// `k` distinct naturals below `n`, `k ≤ n`, in a random order: the first
/// `k` swaps of a Fisher–Yates shuffle of `0 … n − 1`, the places moved
/// kept in a map.
fn distinct(ev: &mut Eval<'_>, k: u64, n: &Natural) -> Vec<Natural> {
    let mut moved: HashMap<Natural, Natural, FixedState> = HashMap::default();
    let mut out = Vec::new();
    for i in 0..k {
        let i = Natural::from(i);
        let rest = n - &i;
        let j = &i + &ev.draw(|g| g.below(&rest));
        let at_j = moved.get(&j).cloned().unwrap_or_else(|| j.clone());
        let at_i = moved.get(&i).cloned().unwrap_or(i);
        moved.insert(j, at_i);
        out.push(at_j);
    }
    out
}

/// A count of draws: a natural number within `Limits::terms`; `Err` with
/// the warning's message when it is a number that does not fit, `Ok(None)`
/// for a symbol.
fn count(ev: &Eval<'_>, k: &Expr) -> Result<Option<u64>, &'static str> {
    if !k.is_number() {
        return Ok(None);
    }
    let Some(k) = k.to_integer().filter(|k| !k.is_negative()) else {
        return Err("a count is a natural number");
    };
    let most = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    match k.to_i64().and_then(|k| u64::try_from(k).ok()) {
        Some(k) if k <= most => Ok(Some(k)),
        _ => Err("too many draws"),
    }
}

/// `e` as a real ball at `wp` bits: a rational or a ball as it is, a
/// numeric expression approximated; `None` for another.
fn real(ev: &Eval<'_>, e: &Expr, wp: u64) -> Option<Ball> {
    match e.view() {
        View::Ball(b) => return Some(b.clone()),
        _ => {
            if let Some(q) = e.to_rational() {
                return Some(Ball::from_rational(&q, wp));
            }
        }
    }
    let value = approximation(ev.context(), e, wp).ok()?;
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `Random[…]`.
pub(super) fn random(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [] => Ok(Some(Expr::from(unit(ev)))),
        [l] if is_list(l) => element(ev, l.args()),
        [n] => {
            let Some(n) = n.to_integer() else {
                return Ok(None);
            };
            if n.is_zero() {
                return Ok(None);
            }
            // [0, n) or (n, 0].
            let toward = if n.is_negative() {
                &n + &Integer::ONE
            } else {
                &n - &Integer::ONE
            };
            Ok(Some(Expr::from(between(ev, &Integer::ZERO, &toward))))
        }
        [k, l] if is_list(l) => sample(ev, e, k, l.args()),
        [a, b] => uniform(ev, a, b),
        [k, a, b] => lotto(ev, e, k, a, b),
        _ => Ok(None),
    }
}

/// `Random[a, b]`: `a + (b − a) u` in balls at the working precision.
fn uniform(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Native {
    let wp = ev.context().working_precision();
    let (Some(low), Some(high)) = (real(ev, a, wp), real(ev, b, wp)) else {
        return Ok(None);
    };
    let fraction = unit(ev);
    let drawn = low.add(&high.sub(&low, wp).mul(&fraction, wp), wp);
    Ok(Some(Expr::from(drawn)))
}

/// `Random[l]`: an element of `l`.
fn element(ev: &mut Eval<'_>, items: &[Expr]) -> Native {
    if items.is_empty() {
        return Ok(None);
    }
    let n = Natural::from(items.len() as u64);
    let i = ev.draw(|g| g.below(&n));
    let i = usize::try_from(&i).expect("below the length");
    Ok(Some(items[i].clone()))
}

/// `Random[k, l]`: `k` elements of `l` at distinct places, in a random
/// order (`RandomSample` with its arguments the other way).
pub(super) fn sample(ev: &mut Eval<'_>, e: &Expr, k: &Expr, items: &[Expr]) -> Native {
    let k = match count(ev, k) {
        Ok(Some(k)) => k,
        Ok(None) => return Ok(None),
        Err(why) => return refuse(ev, e, why),
    };
    let n = items.len() as u64;
    if k > n {
        return refuse(ev, e, "more draws than elements");
    }
    let places = distinct(ev, k, &Natural::from(n));
    let picked = places.iter().map(|i| {
        let i = usize::try_from(i).expect("below the length");
        items[i].clone()
    });
    Ok(Some(list(picked.collect())))
}

/// `Random[k, a, b]`: `k` distinct integers between `a` and `b`.
fn lotto(ev: &mut Eval<'_>, e: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Native {
    if !(a.is_number() && b.is_number()) {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (a.to_integer(), b.to_integer()) else {
        return refuse(ev, e, "the ends are integers");
    };
    let k = match count(ev, k) {
        Ok(Some(k)) => k,
        Ok(None) => return Ok(None),
        Err(why) => return refuse(ev, e, why),
    };
    let lo = if a <= b { &a } else { &b };
    let width = (&a - &b).unsigned_abs() + Natural::from(1u64);
    if Natural::from(k) > width {
        return refuse(ev, e, "more draws than integers in the range");
    }
    let picked = distinct(ev, k, &width)
        .into_iter()
        .map(|i| Expr::from(lo + &Integer::from(i)));
    Ok(Some(list(picked.collect())))
}

/// `RandomInteger[n]`, `RandomInteger[a, b]`.
pub(super) fn random_integer(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (a, b) = match e.args() {
        [n] => (Expr::from(1), n.clone()),
        [a, b] => (a.clone(), b.clone()),
        _ => return Ok(None),
    };
    if !(a.is_number() && b.is_number()) {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (a.to_integer(), b.to_integer()) else {
        return refuse(ev, e, "the ends are integers");
    };
    Ok(Some(Expr::from(between(ev, &a, &b))))
}

pub(super) mod collections;
pub(super) mod variates;

#[cfg(test)]
mod tests;
