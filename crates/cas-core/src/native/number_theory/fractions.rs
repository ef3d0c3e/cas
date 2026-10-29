//! Continued fractions (CAS-DOC-39 §5, CAS-424): `ContinuedFraction`,
//! Xcas's `dfc`, and `FromContinuedFraction`, its `dfc2f`.
//!
//! - `ContinuedFraction[x]` of a rational: its quotients by Euclid's
//!   algorithm, floors, so that `−7/3` is `{−3, 1, 2}`; giac's
//!   `dfc(355/113)` is `[3,7,16]`.
//! - Of a quadratic irrational `a + b √D`, `a` and `b` rational: written
//!   `(P + √N)/Q` with `Q | N − P²`, its quotients `⌊(P + √N)/Q⌋` exactly,
//!   `⌊(P + s)/Q⌋` for `Q > 0` and `⌊(P + s + 1)/Q⌋` for `Q < 0`, `s =
//!   ⌊√N⌋`, and the state `(P, Q)` carried as `P′ = a Q − P`, `Q′ = (N −
//!   P′²)/Q`; the period starts at the first reduced state, above 1 with
//!   its conjugate in `(−1, 0)` (Galois), and ends when that state comes
//!   back, the last entry as a list: `√7` is `{2, {1, 1, 1, 4}}` and the
//!   golden ratio `{{1}}`. giac writes `dfc(sqrt(2),5)` as `[1,2,[2]]`, a
//!   period started late. Past [`MAX_PERIODIC`] quotients the call stays.
//! - `ContinuedFraction[x, n]` of a real number: the first `n` quotients,
//!   each floor proven by `N`'s ball lying between two integers, the
//!   precision doubled up to [`MAX_BITS`] where it does not; the quotients
//!   proven, as many as they are. A rational's has its quotients up to `n`.
//! - `FromContinuedFraction[{a₀, …}]`: the rational number; with a period
//!   last, the quadratic irrational: the period's value `y` solves `y =
//!   (A y + B)/(C y + D)`, its matrix `Π (pᵢ 1; 1 0)`, the positive root,
//!   then the quotients before it apply, the result rationalized as `r + s
//!   √Δ`.

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The highest precision of the numeric quotients.
const MAX_BITS: u64 = 1 << 14;

/// The most quotients asked.
const MAX_TERMS: i64 = 1 << 16;

/// The most quotients of a periodic expansion, before and in its period.
const MAX_PERIODIC: usize = 1 << 16;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// A rational's quotients, at most `limit`.
fn rational_quotients(x: &Rational, limit: usize) -> Vec<Integer> {
    let (mut p, mut q) = (x.numer().clone(), x.denom().clone());
    let mut out = Vec::new();
    while !q.is_zero() && out.len() < limit {
        let (a, r) = p.div_rem_euclid(&q);
        out.push(a);
        p = q;
        q = r;
    }
    out
}

/// `x` as `a + b √D`, `b ≠ 0`, `D > 1` an integer: the terms of its sum.
fn quadratic(x: &Expr) -> Option<(Rational, Rational, Integer)> {
    let terms: &[Expr] = match head(x) {
        Some(Symbol::PLUS) => x.args(),
        _ => std::slice::from_ref(x),
    };
    let mut a = Rational::ZERO;
    let mut root: Option<(Rational, Integer)> = None;
    for t in terms {
        if let Some(q) = t.to_rational() {
            a = &a + &q;
            continue;
        }
        let (c, power) = match (head(t), t.args()) {
            (Some(Symbol::TIMES), [c, p]) => (c.to_rational()?, p),
            _ => (Rational::ONE, t),
        };
        let (Some(Symbol::POWER), [d, half]) = (head(power), power.args()) else {
            return None;
        };
        let d = d.to_integer()?;
        if half.to_rational()? != Rational::new(Integer::ONE, Integer::from(2)) || d <= Integer::ONE
        {
            return None;
        }
        if root.is_some() {
            return None;
        }
        root = Some((c, d));
    }
    let (b, d) = root?;
    Some((a, b, d))
}

/// `⌊(P + √N)/Q⌋`, `N` not a square, `s = ⌊√N⌋`.
fn quotient(p: &Integer, s: &Integer, q: &Integer) -> Integer {
    let top = if q.is_negative() {
        p + &(s + &Integer::ONE)
    } else {
        p + s
    };
    Rational::new(top, q.clone()).floor()
}

/// Whether `(P + √N)/Q` is reduced, above 1 with its conjugate in `(−1,
/// 0)`, `s = ⌊√N⌋`: by Galois's theorem, its expansion is then purely
/// periodic, so the period starts at the first reduced state.
fn reduced(p: &Integer, q: &Integer, s: &Integer) -> bool {
    q.is_positive() && p <= s && &(p + s) >= q && s < &(p + q)
}

/// The quotients of `(a + b √D)` with the start of its period; `None` past
/// [`MAX_PERIODIC`] quotients.
fn periodic(a: &Rational, b: &Rational, d: &Integer) -> Option<(Vec<Integer>, usize)> {
    // a = a₁/m, b = b₁/m: (a₁ + b₁ √D)/m = (P + √N)/Q, N = b₁² D.
    let m = a.denom().lcm(b.denom());
    let a1 = (a * &Rational::from(m.clone())).numer().clone();
    let b1 = (b * &Rational::from(m.clone())).numer().clone();
    let (mut p, mut q) = if b1.is_negative() { (-a1, -m) } else { (a1, m) };
    let mut n = &(&b1 * &b1) * d;
    // Q | N − P², or both scaled by |Q|.
    let (_, rest) = (&n - &(&p * &p)).div_rem_euclid(&q);
    if !rest.is_zero() {
        let size = Integer::from(q.unsigned_abs());
        p = &p * &size;
        n = &(&n * &size) * &size;
        q = &q * &size;
    }
    let (s, _) = n.sqrt_rem();
    let mut start: Option<(usize, Integer, Integer)> = None;
    let mut out = Vec::new();
    while out.len() <= MAX_PERIODIC {
        if let Some((k, p0, q0)) = &start {
            if (&p, &q) == (p0, q0) {
                return Some((out, *k));
            }
        } else if reduced(&p, &q, &s) {
            start = Some((out.len(), p.clone(), q.clone()));
        }
        let c = quotient(&p, &s, &q);
        let next = &(&c * &q) - &p;
        q = (&n - &(&next * &next)).div_exact(&q);
        p = next;
        out.push(c);
    }
    None
}

/// The lower and upper ends of a ball, rational; `None` when infinite.
fn ends(b: &Ball) -> Option<(Rational, Rational)> {
    let mid = b.mid().to_rational();
    let rad = b.rad().to_float()?.to_rational();
    Some((&mid - &rad, &mid + &rad))
}

/// Up to `n` quotients of the real `x` proven at `bits`.
fn numeric(ev: &mut Eval<'_>, x: &Expr, n: usize, bits: u64) -> Result<Vec<Integer>, EvalError> {
    let digits = Expr::from(i64::try_from(bits * 3 / 10 + 2).expect("small"));
    let value = ev.eval(&Expr::call(Symbol::N, [x.clone(), digits]))?;
    let View::Ball(start) = value.view() else {
        return Ok(Vec::new());
    };
    let mut y = start.clone();
    let mut out = Vec::new();
    while out.len() < n {
        let Some((lo, hi)) = ends(&y) else {
            break;
        };
        let a = lo.floor();
        if hi.floor() != a {
            break;
        }
        out.push(a.clone());
        let rest = y.sub(&Ball::from_integer(&a), bits);
        y = Ball::from_integer(&Integer::ONE).div(&rest, bits);
    }
    Ok(out)
}

/// `ContinuedFraction[x]` and `ContinuedFraction[x, n]`.
pub(in crate::native) fn continued_fraction(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, n) = match e.args() {
        [x] => (x, None),
        [x, n] => {
            let Some(n) = n.as_i64().filter(|n| (1..=MAX_TERMS).contains(n)) else {
                return Ok(None);
            };
            (x, Some(usize::try_from(n).expect("small")))
        }
        _ => return Ok(None),
    };
    let ints = |v: Vec<Integer>| v.into_iter().map(Expr::from).collect::<Vec<_>>();
    if let Some(q) = x.to_rational() {
        return Ok(Some(list(ints(rational_quotients(
            &q,
            n.unwrap_or(usize::MAX),
        )))));
    }
    if n.is_none()
        && let Some((a, b, d)) = quadratic(x)
    {
        let Some((quotients, start)) = periodic(&a, &b, &d) else {
            return Ok(None);
        };
        let mut items = ints(quotients);
        let period = items.split_off(start);
        items.push(list(period));
        return Ok(Some(list(items)));
    }
    let Some(n) = n else {
        return Ok(None);
    };
    let mut bits = 64 * u64::try_from(n).unwrap_or(u64::MAX).min(64) + 64;
    loop {
        let found = numeric(ev, x, n, bits)?;
        if found.len() == n || bits >= MAX_BITS {
            return Ok((!found.is_empty()).then(|| list(ints(found))));
        }
        bits *= 2;
    }
}

/// The matrix of the quotients `Π (aᵢ 1; 1 0)`: `[p_k, p_(k−1), q_k,
/// q_(k−1)]`, the convergents' numerators and denominators.
fn matrix(quotients: &[Integer]) -> [Integer; 4] {
    let (mut p, mut p0, mut q, mut q0) = (Integer::ONE, Integer::ZERO, Integer::ZERO, Integer::ONE);
    for a in quotients {
        let next_p = &(a * &p) + &p0;
        let next_q = &(a * &q) + &q0;
        (p0, p, q0, q) = (p, next_p, q, next_q);
    }
    [p, p0, q, q0]
}

/// `FromContinuedFraction[{a₀, …}]`.
pub(in crate::native) fn from_continued_fraction(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [quotients] = e.args() else {
        return Ok(None);
    };
    if head(quotients) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let items = quotients.args();
    let (prefix, period) = match items.split_last() {
        Some((last, rest)) if head(last) == Some(Symbol::LIST) => (rest, Some(last.args())),
        _ => (items, None),
    };
    let integers = |v: &[Expr]| v.iter().map(Expr::to_integer).collect::<Option<Vec<_>>>();
    let Some(prefix) = integers(prefix) else {
        return Ok(None);
    };
    let [p, p0, q, q0] = matrix(&prefix);
    let Some(period) = period else {
        if prefix.is_empty() {
            return Ok(None);
        }
        return Ok(Some(Expr::from(Rational::new(p, q))));
    };
    let Some(period) = integers(period) else {
        return Ok(None);
    };
    if period.is_empty() || period.iter().any(|a| *a < Integer::ONE) {
        return Ok(None);
    }
    // y = (A y + B)/(C y + D): C y² + (D − A) y − B = 0.
    let [a, b, c, d] = matrix(&period);
    let (u, delta, v) = (
        &a - &d,
        &(&(&d - &a) * &(&d - &a)) + &(&(&c * &b) * &Integer::from(4)),
        &c * &Integer::from(2),
    );
    // y = (u + √Δ)/v; x = (p y + p0)/(q y + q0) = (α + β√Δ)/(γ + δ√Δ).
    let (alpha, beta) = (&(&p * &u) + &(&p0 * &v), p.clone());
    let (gamma, delta_coefficient) = (&(&q * &u) + &(&q0 * &v), q.clone());
    let denominator = &(&gamma * &gamma) - &(&(&delta_coefficient * &delta_coefficient) * &delta);
    let rational = &(&alpha * &gamma) - &(&(&beta * &delta_coefficient) * &delta);
    let irrational = &(&beta * &gamma) - &(&alpha * &delta_coefficient);
    let r = Rational::new(rational, denominator.clone());
    let s = Rational::new(irrational, denominator);
    let half = Rational::new(Integer::ONE, Integer::from(2));
    Ok(Some(
        Expr::from(r) + Expr::from(s) * Expr::from(delta).pow(Expr::from(half)),
    ))
}

#[cfg(test)]
mod tests;
