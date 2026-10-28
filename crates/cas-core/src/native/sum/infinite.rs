//! Infinite sums (CAS-DOC-13 §4, CAS-226): `Sum[f, k, a, Infinity]`, giac's
//! `sum(f,k,a,inf)`.
//!
//! - A term that does not tend to 0 (`Limit`, CAS-DOC-10) diverges, to the
//!   infinity of its limit's sign (giac's `+infinity` for `sum(1,k,1,inf)`);
//!   without a sign, or a limit, the search goes on: `(-1)^k/k` has none
//!   here, and `(-1)^k` stays where giac answers `undef`.
//! - A pole at an integer `k ≥ a`, a root of a factor of the denominator
//!   linear in `k` alone, leaves the sum, as giac's `undef` for
//!   `sum(1/k^2,k,0,inf)`.
//! - Known series by the ratio `r = t(k + 1)/t(k)` of the term (CAS-225)
//!   and its first term `t(a)`: geometric, `r = c`, `t(a)/(1 − c)` for
//!   `|c| < 1`; exponential, `r = c/(k + 1)`, `t(a) a!/c^a (e^c −
//!   Σ_{j<a} c^j/j!)`; logarithmic, `r = c k/(k + 1)` from `a ≥ 1`,
//!   `t(a) a/c^a (−ln(1 − c) − Σ_{1≤j<a} c^j/j)` for `−1 ≤ c < 1`, and the
//!   harmonic series' infinity for `c = 1`; `c/k^s` from `a ≥ 1`, `s ≥ 2`,
//!   `c` times `Zeta[s]` less the terms before `a`, and `c (−1)^k/k^s` with
//!   `−(1 − 2^(1−s)) Zeta[s]` (CAS-245): exact for an even `s`, `Zeta[3]`
//!   for an odd one, as giac's.
//! - A hypergeometric antidifference `T` (Gosper's, CAS-225):
//!   `lim_{n→∞} T(n) − T(a)`, when `Limit` finds it; without one, a
//!   convergent rational function by `PolyGamma` (CAS-DOC-31 §6).
//! - A sum term by term, each term's sum found or left: finite sums add,
//!   with the others under one `Sum`; infinities of one sign, with only
//!   finite sums beside them, are that infinity.
//! - `−Infinity` as the lower bound: `Σ_{k=−∞}^{b} f(k) =
//!   Σ_{k=−b}^{∞} f(−k)`; both bounds infinite, the halves from 0 and 1.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of series: t in k from a, its ratio r and c"
)]

use std::cmp::Ordering;

use cas_num::algebra::{Integers, Ring};
use cas_num::{Integer, Rational};
use cas_poly::SparseRing;

use super::super::approximation;
use super::super::polynomial::{convert, factored};
use super::gosper;
use crate::assume::{Sign, sign_of};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::Fraction;
use crate::symbol::Symbol;

/// The highest `2m` of `ζ(2m)` written out.
const MAX_ZETA: i64 = 64;

/// The largest integer start of a known series but the geometric one: the
/// terms before it are subtracted exactly, and those of `1/k^64` have
/// denominators of `64 log₂ lcm(1, …, a)` bits.
const MAX_START: i64 = 256;

/// The sum of a series, when found.
#[derive(Clone, Debug, PartialEq)]
enum Total {
    Finite(Expr),
    /// `Infinity` (`true`) or `-Infinity`.
    Infinite(bool),
}

impl Total {
    fn expr(self) -> Expr {
        match self {
            Total::Finite(e) => e,
            Total::Infinite(true) => Expr::from(Symbol::INFINITY),
            Total::Infinite(false) => -Expr::from(Symbol::INFINITY),
        }
    }

    /// The sum of two totals; `None` for infinities of opposite signs.
    fn add(self, other: Total) -> Option<Total> {
        match (self, other) {
            (Total::Finite(x), Total::Finite(y)) => Some(Total::Finite(x + y)),
            (Total::Infinite(s), Total::Finite(_)) | (Total::Finite(_), Total::Infinite(s)) => {
                Some(Total::Infinite(s))
            }
            (Total::Infinite(s), Total::Infinite(u)) => (s == u).then_some(Total::Infinite(s)),
        }
    }
}

/// Which infinity `x` is: `Some(true)` for `Infinity`, `Some(false)` for
/// `-Infinity`, evaluated.
pub(super) fn infinity(x: &Expr) -> Option<bool> {
    if *x == Expr::from(Symbol::INFINITY) {
        return Some(true);
    }
    (*x == Expr::call(Symbol::DIRECTED_INFINITY, [Expr::from(-1)])).then_some(false)
}

/// Whether `x` holds an infinity or `Indeterminate`.
fn undefined(x: &Expr) -> bool {
    [
        Symbol::INFINITY,
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
        Symbol::INDETERMINATE,
    ]
    .into_iter()
    .any(|s| x.contains(&Expr::from(s)))
}

/// `Sum[f, k, a, b]` with an infinite bound; `None` when it stays.
pub(super) fn sum(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let reflected = |ev: &mut Eval<'_>| ev.eval(&f.substitute(&[(k.clone(), -k.clone())]));
    let total = match (infinity(a), infinity(b)) {
        (None, Some(true)) => series(ev, e, f, k, a)?,
        (Some(false), None) => {
            let g = reflected(ev)?;
            let start = ev.eval(&-b.clone())?;
            series(ev, e, &g, k, &start)?
        }
        (Some(false), Some(true)) => {
            let g = reflected(ev)?;
            match (
                series(ev, e, f, k, &Expr::from(0))?,
                series(ev, e, &g, k, &Expr::from(1))?,
            ) {
                (Some(x), Some(y)) => x.add(y),
                _ => None,
            }
        }
        _ => None,
    };
    Ok(total.map(Total::expr))
}

/// `Σ_{k=a}^{∞} f`: the whole term, else term by term.
fn series(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    k: &Expr,
    a: &Expr,
) -> Result<Option<Total>, EvalError> {
    if let Some(total) = term(ev, e, f, k, a)? {
        return Ok(Some(total));
    }
    if f.head().and_then(Expr::as_symbol) != Some(Symbol::PLUS) {
        return Ok(None);
    }
    let mut total: Option<Total> = None;
    let mut open = Vec::new();
    for t in f.args() {
        match term(ev, e, t, k, a)? {
            Some(s) => {
                total = match total {
                    None => Some(s),
                    Some(so_far) => match so_far.add(s) {
                        Some(sum) => Some(sum),
                        None => return Ok(None),
                    },
                };
            }
            None => open.push(t.clone()),
        }
    }
    Ok(match total {
        Some(Total::Finite(s)) if !open.is_empty() => {
            let rest = Expr::call(
                Symbol::SUM,
                [
                    Expr::call(Symbol::PLUS, open),
                    k.clone(),
                    a.clone(),
                    Expr::from(Symbol::INFINITY),
                ],
            );
            Some(Total::Finite(s + rest))
        }
        Some(Total::Infinite(_)) if !open.is_empty() => None,
        total => total,
    })
}

/// `Σ_{k=a}^{∞} t` for one term `t`.
fn term(
    ev: &mut Eval<'_>,
    e: &Expr,
    t: &Expr,
    k: &Expr,
    a: &Expr,
) -> Result<Option<Total>, EvalError> {
    if !t.contains(k) {
        if *t == Expr::from(0) {
            return Ok(Some(Total::Finite(Expr::from(0))));
        }
        return Ok(sign(ev, t).map(Total::Infinite));
    }
    if let Some(tail) = limit(ev, t, k)?
        && tail != Expr::from(0)
    {
        return Ok(sign(ev, &tail).map(Total::Infinite));
    }
    if pole(ev, e, t, k, a)? {
        return Ok(None);
    }
    if let Some(total) = known(ev, t, k, a)? {
        return Ok(Some(total));
    }
    let Some(antidifference) = gosper::antidifference(ev, e, t, k)? else {
        // A rational function by `PolyGamma` (CAS-DOC-31 §6).
        return Ok(super::rational::infinite(ev, t, k, a)?.map(Total::Finite));
    };
    let Some(end) = limit(ev, &antidifference, k)? else {
        return Ok(None);
    };
    if let Some(s) = infinity(&end) {
        return Ok(Some(Total::Infinite(s)));
    }
    let start = at(ev, &antidifference, k, a)?;
    if undefined(&end) || undefined(&start) {
        return Ok(None);
    }
    Ok(Some(Total::Finite(ev.eval(&(end - start))?)))
}

/// `Limit[f, k, Infinity]`, evaluated; `None` when it stays.
fn limit(ev: &mut Eval<'_>, f: &Expr, k: &Expr) -> Result<Option<Expr>, EvalError> {
    let value = ev.eval(&Expr::call(
        Symbol::LIMIT,
        [f.clone(), k.clone(), Expr::from(Symbol::INFINITY)],
    ))?;
    Ok((value.head().and_then(Expr::as_symbol) != Some(Symbol::LIMIT)).then_some(value))
}

/// The sign of a real `x`, positive `true`: by the assumptions, else by its
/// ball; `None` when neither tells, as for `x` or `ComplexInfinity`.
fn sign(ev: &Eval<'_>, x: &Expr) -> Option<bool> {
    match sign_of(ev.context(), x) {
        Sign::Positive => return Some(true),
        Sign::Negative => return Some(false),
        _ => {}
    }
    if undefined(x) {
        return None;
    }
    let value = approximation(ev.context(), x, 128).ok()?;
    match value.view() {
        View::Ball(b) => match b.sign()? {
            Ordering::Greater => Some(true),
            Ordering::Less => Some(false),
            Ordering::Equal => None,
        },
        _ => None,
    }
}

/// `f` at `k = x`, evaluated.
fn at(ev: &mut Eval<'_>, f: &Expr, k: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&f.substitute(&[(k.clone(), x.clone())]))
}

/// Whether `f` has a pole at `k = a + j` for an integer `j ≥ 0`, `a`
/// rational: a root of a factor of its denominator linear in `k` alone.
fn pole(ev: &mut Eval<'_>, e: &Expr, f: &Expr, k: &Expr, a: &Expr) -> Result<bool, EvalError> {
    let Some(a) = a.to_rational() else {
        return Ok(false);
    };
    let Some(p) = convert(ev, e, &[f.clone(), k.clone()], true)? else {
        return Ok(false);
    };
    let (Some(kv), Fraction::Integers { den, .. }) = (p.variable(k), &p.fractions()[0]) else {
        return Ok(false);
    };
    let ring = SparseRing::new(Integers, p.generators().len());
    if ring.degree(den, kv).unwrap_or(0) == 0 {
        return Ok(false);
    }
    let Some(powers) = factored(ev, e, &ring, den, &ring.one())? else {
        return Ok(false);
    };
    Ok(powers.factors.iter().any(|(g, _)| {
        gosper::linear_root(&ring, g, kv).is_some_and(|root| {
            let j = &root - &a;
            j.is_integer() && !j.is_negative()
        })
    }))
}

/// `1 − x`, evaluated, positive.
fn below_one(ev: &mut Eval<'_>, x: &Expr) -> Result<bool, EvalError> {
    let gap = ev.eval(&(Expr::from(1) - x.clone()))?;
    Ok(sign(ev, &gap) == Some(true))
}

/// The known series of §4 for the term `t`; `None` for another.
fn known(ev: &mut Eval<'_>, t: &Expr, k: &Expr, a: &Expr) -> Result<Option<Total>, EvalError> {
    let start = a.as_i64().filter(|s| (0..=MAX_START).contains(s));
    if let Some(start) = start
        && let Some(total) = zeta(ev, t, k, start)?
    {
        return Ok(Some(Total::Finite(total)));
    }
    let Some(r) = gosper::ratio(ev, t, k)? else {
        return Ok(None);
    };
    let r = ev.eval(&Expr::call(Symbol::TOGETHER, [r]))?;
    // Geometric: r = c.
    if !r.contains(k) {
        let modulus = ev.eval(&Expr::call(Symbol::ABS, [r.clone()]))?;
        if !below_one(ev, &modulus)? {
            return Ok(None);
        }
        let first = at(ev, t, k, a)?;
        return Ok(Some(Total::Finite(first / (Expr::from(1) - r))));
    }
    let Some(start) = start else {
        return Ok(None);
    };
    // Exponential: r = c/(k + 1).
    let c = ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [r * (k.clone() + Expr::from(1))],
    ))?;
    if !c.contains(k) {
        let first = at(ev, t, k, a)?;
        let head: Vec<Expr> = (0..start)
            .map(|j| c.pow(j) / Expr::from(factorial(j)))
            .collect();
        let scale = first * Expr::from(factorial(start)) / c.pow(start);
        let rest = Expr::call(Symbol::EXP, [c]) - Expr::call(Symbol::PLUS, head);
        return Ok(Some(Total::Finite(scale * rest)));
    }
    // Logarithmic: r = c k/(k + 1), from k = 1.
    let c = ev.eval(&Expr::call(Symbol::TOGETHER, [c / k.clone()]))?;
    if c.contains(k) || start < 1 {
        return Ok(None);
    }
    let first = at(ev, t, k, a)?;
    if c == Expr::from(1) {
        return Ok(sign(ev, &first).map(Total::Infinite));
    }
    // −1 ≤ c < 1.
    let shifted = ev.eval(&(c.clone() + Expr::from(1)))?;
    if !below_one(ev, &c)? || (shifted != Expr::from(0) && sign(ev, &shifted) != Some(true)) {
        return Ok(None);
    }
    let head: Vec<Expr> = (1..start).map(|j| c.pow(j) / Expr::from(j)).collect();
    let scale = first * Expr::from(start) / c.pow(start);
    let rest = -Expr::call(Symbol::LOG, [Expr::from(1) - c]) - Expr::call(Symbol::PLUS, head);
    Ok(Some(Total::Finite(scale * rest)))
}

/// `n!`.
fn factorial(n: i64) -> Integer {
    (1..=n).fold(Integer::ONE, |f, i| &f * &Integer::from(i))
}

/// `t` as `c/k^s` or `c (−1)^k/k^s` with `c` free of `k`: `c`, `s` and
/// whether it alternates.
fn power_term(t: &Expr, k: &Expr) -> Option<(Expr, i64, bool)> {
    let factors: Vec<&Expr> = match t.head().and_then(Expr::as_symbol) {
        Some(Symbol::TIMES) => t.args().iter().collect(),
        _ => vec![t],
    };
    let mut c = Vec::new();
    let mut s = None;
    let mut alternating = false;
    for x in factors {
        if !x.contains(k) {
            c.push(x.clone());
            continue;
        }
        match (x.head().and_then(Expr::as_symbol), x.args()) {
            (Some(Symbol::POWER), [base, n]) if base == k && s.is_none() => {
                s = Some(n.as_i64()?.checked_neg()?);
            }
            (Some(Symbol::POWER), [base, n])
                if *base == Expr::from(-1) && n == k && !alternating =>
            {
                alternating = true;
            }
            _ => return None,
        }
    }
    Some((Expr::call(Symbol::TIMES, c), s?, alternating))
}

/// `Σ_{k≥start} c/k^s`, `s ≥ 2`, or with `(−1)^k`: `c` times `Zeta[s]` (or
/// `−(1 − 2^(1−s)) Zeta[s]`) less the terms before `start`, `Zeta[s]` exact
/// for an even `s` and staying for an odd one (CAS-245); `None` for
/// another term.
fn zeta(ev: &mut Eval<'_>, t: &Expr, k: &Expr, start: i64) -> Result<Option<Expr>, EvalError> {
    let Some((c, s, alternating)) = power_term(t, k) else {
        return Ok(None);
    };
    if !(2..=MAX_ZETA).contains(&s) || start < 1 {
        return Ok(None);
    }
    let exponent = u64::try_from(s).expect("a small power");
    let mut value = Expr::call(Symbol::ZETA, [Expr::from(s)]);
    if alternating {
        // −(1 − 2^(1−s)).
        let half = Rational::new(Integer::ONE, Integer::from(2).pow(exponent - 1));
        value = Expr::from(&half - &Rational::ONE) * value;
    }
    let mut head = Rational::ZERO;
    for j in 1..start {
        let sign = if alternating && j % 2 == 1 { -1 } else { 1 };
        let term = Rational::new(Integer::from(sign), Integer::from(j).pow(exponent));
        head = &head + &term;
    }
    Ok(Some(ev.eval(&(c * (value - Expr::from(head))))?))
}
