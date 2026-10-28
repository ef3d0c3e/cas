//! Piecewise integrands (CAS-DOC-25 §5, CAS-316): `Integrate[f, x, a, b]`
//! where `f` has `Abs[u]`, `Sign[u]` or `HeavisideTheta[u]` with `u`
//! depending on `x` splits at the zeros of the `u` strictly between the
//! bounds:
//!
//! - `−β/α` for a linear `u = αx + β`, `α` and `β` free of `x`;
//! - for a polynomial over ℚ of degree [`MAX_DEGREE`] at most, the roots of
//!   its linear factors and the real ones of its quadratic factors, in
//!   radicals; a factor of a higher degree must have no real root;
//! - for `sin v` or `cos v`, `v = αx + β` with real numbers `α` and `β`,
//!   between finite numeric bounds, the points where `v` is `kπ` or
//!   `(k + 1/2)π`, [`MAX_POINTS`] at most (CAS-317): `|sin x|` over a
//!   period.
//!
//! On each piece, the sign `σ` of `u` at its middle (`p ± 1` on an infinite
//! piece) makes `Abs[u]` `σu`, `Sign[u]` `σ` and `HeavisideTheta[u]` 1 or
//! 0, and the pieces integrate as before (CAS-198, CAS-207): the integral
//! is their sum. Two points are ordered by their infinities, else by the
//! sign of their difference under the assumptions, else by its ball, so
//! that 0 lies between `−T/2` and `T/2` for `T > 0`. Another argument, an
//! order or a sign that none of these tells, or a piece that stays leaves
//! the integral, the piece with its warning.

use std::cmp::Ordering;

use cas_num::algebra::Rationals;
use cas_num::{Integer, Rational};
use cas_poly::DenseRing;

use super::definite::{ball_ends, end, infinite};
use super::{polynomial, slope};
use crate::assume::{Sign, sign_of};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

/// The heads that switch at the zeros of their argument.
const SWITCHES: [Symbol; 3] = [Symbol::ABS, Symbol::SIGN, Symbol::HEAVISIDE_THETA];

/// Polynomial arguments of a higher degree are not factored.
pub(super) const MAX_DEGREE: usize = 8;

/// The points between the bounds, at most: as many pieces and one more.
const MAX_POINTS: usize = 64;

/// The bits of the balls that tell signs.
const BITS: u64 = 128;

/// The distinct `Abs[u]`, `Sign[u]` and `HeavisideTheta[u]` of `f` with
/// `u` depending on `x`.
pub(super) fn switches(f: &Expr, x: &Expr) -> Vec<Expr> {
    f.distinct()
        .filter(|g| match (g.head().and_then(Expr::as_symbol), g.args()) {
            (Some(h), [u]) => SWITCHES.contains(&h) && u.contains(x),
            _ => false,
        })
        .cloned()
        .collect()
}

/// The sign of a real `v`: by its value, the assumptions, or its ball;
/// `None` when none tells.
fn sign(ev: &Eval<'_>, v: &Expr) -> Option<Ordering> {
    if let Some(r) = v.to_rational() {
        return Some(r.cmp(&Rational::ZERO));
    }
    match sign_of(ev.context(), v) {
        Sign::Positive => return Some(Ordering::Greater),
        Sign::Negative => return Some(Ordering::Less),
        Sign::Zero => return Some(Ordering::Equal),
        _ => {}
    }
    let value = approximation(ev.context(), v, BITS).ok()?;
    match value.view() {
        View::Ball(b) => b.sign(),
        _ => None,
    }
}

/// How `p` compares with `q` on the extended real line: by their
/// infinities, else by the sign of `q − p`; `None` when it is not told.
fn order(ev: &mut Eval<'_>, p: &Expr, q: &Expr) -> Result<Option<Ordering>, EvalError> {
    match (end(p), end(q)) {
        (Some(s), Some(t)) => return Ok(Some(s.cmp(&t))),
        (Some(true), None) | (None, Some(false)) => return Ok(Some(Ordering::Greater)),
        (Some(false), None) | (None, Some(true)) => return Ok(Some(Ordering::Less)),
        (None, None) if infinite(p) || infinite(q) => return Ok(None),
        (None, None) => {}
    }
    let difference = ev.eval(&(q.clone() - p.clone()))?;
    Ok(sign(ev, &difference).map(Ordering::reverse))
}

/// The real zeros of `u` (the module's list), as exact points, between
/// the bounds `lo < hi` and maybe beyond; `None` for another `u`, or a
/// factor of a higher degree with a real root.
fn zeros(
    ev: &mut Eval<'_>,
    e: &Expr,
    u: &Expr,
    x: &Expr,
    bounds: [&Expr; 2],
) -> Result<Option<Vec<Expr>>, EvalError> {
    if let (Some(h @ (Symbol::SIN | Symbol::COS)), [v]) =
        (u.head().and_then(Expr::as_symbol), u.args())
    {
        return periodic_zeros(ev, h == Symbol::COS, v, x, bounds);
    }
    if let Some(alpha) = slope(u, x) {
        let alpha = ev.eval(&alpha)?;
        if alpha == Expr::from(0) {
            return Ok(None);
        }
        let beta = ev.eval(&u.substitute(&[(x.clone(), Expr::from(0))]))?;
        return Ok(Some(vec![ev.eval(&(-(beta / alpha)))?]));
    }
    let Some(p) = polynomial(ev, e, u, x)? else {
        return Ok(None);
    };
    if p.degree().is_none_or(|d| d > MAX_DEGREE) {
        return Ok(None);
    }
    let qx = DenseRing::new(Rationals);
    let (two, four) = (
        Rational::from(Integer::from(2)),
        Rational::from(Integer::from(4)),
    );
    let half = Expr::from(Rational::new(Integer::from(1), Integer::from(2)));
    let mut roots = Vec::new();
    for (factor, _) in qx.factor(&p).factors() {
        let c = factor.coeffs();
        match factor.degree() {
            Some(0) | None => {}
            Some(1) => roots.push(Expr::from(-(&c[0] / &c[1]))),
            Some(2) => {
                let discriminant = &(&c[1] * &c[1]) - &(&four * &(&c[0] * &c[2]));
                if discriminant.is_negative() {
                    continue;
                }
                let root = Expr::from(discriminant).pow(half.clone());
                for s in [-1, 1] {
                    let r = (Expr::from(-c[1].clone()) + Expr::from(s) * root.clone())
                        / Expr::from(&two * &c[2]);
                    roots.push(ev.eval(&r)?);
                }
            }
            Some(_) if qx.count_real_roots(factor, None, None) == 0 => {}
            Some(_) => return Ok(None),
        }
    }
    Ok(Some(roots))
}

/// The zeros `((k + s)π − β)/α` of `sin v` (`s = 0`), or of `cos v` (`s =
/// 1/2`) when `cosine`, for `v = αx + β` with real numbers `α ≠ 0` and `β`
/// and finite numeric bounds (CAS-317): for the integers `k` that may put
/// them between the bounds, by the balls of `v/π` at the bounds, one more
/// on each side, for [`points`] to sort out; `None` for another `v` or
/// bounds, or more than [`MAX_POINTS`].
fn periodic_zeros(
    ev: &mut Eval<'_>,
    cosine: bool,
    v: &Expr,
    x: &Expr,
    bounds: [&Expr; 2],
) -> Result<Option<Vec<Expr>>, EvalError> {
    let Some(alpha) = slope(v, x) else {
        return Ok(None);
    };
    let alpha = ev.eval(&alpha)?;
    let beta = ev.eval(&v.substitute(&[(x.clone(), Expr::from(0))]))?;
    let pi = Expr::from(Symbol::PI);
    let mut ends = Vec::with_capacity(4);
    for bound in bounds {
        let turns = ev.eval(&((alpha.clone() * bound.clone() + beta.clone()) / pi.clone()))?;
        let Some((low, high)) = ball_ends(ev, &turns) else {
            return Ok(None);
        };
        ends.extend([low, high]);
    }
    let (Some(low), Some(high)) = (ends.iter().min(), ends.iter().max()) else {
        return Ok(None);
    };
    let (Some(first), Some(last)) = (low.floor().to_i64(), high.ceil().to_i64()) else {
        return Ok(None);
    };
    let (first, last) = (first.saturating_sub(1), last.saturating_add(1));
    if alpha == Expr::from(0)
        || last.saturating_sub(first) > i64::try_from(MAX_POINTS).expect("small")
    {
        return Ok(None);
    }
    let shift = Expr::from(Rational::new(
        Integer::from(i64::from(cosine)),
        Integer::from(2),
    ));
    let mut zeros = Vec::new();
    for k in first..=last {
        let z = ((Expr::from(k) + shift.clone()) * pi.clone() - beta.clone()) / alpha.clone();
        zeros.push(ev.eval(&z)?);
    }
    Ok(Some(zeros))
}

/// The zeros of the switches' arguments strictly between `lo < hi`,
/// ascending and distinct; `None` when one is not found or not placed.
fn points(
    ev: &mut Eval<'_>,
    e: &Expr,
    switches: &[Expr],
    x: &Expr,
    [lo, hi]: [&Expr; 2],
) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut points: Vec<Expr> = Vec::new();
    for s in switches {
        let Some(zeros) = zeros(ev, e, &s.args()[0], x, [lo, hi])? else {
            return Ok(None);
        };
        for z in zeros {
            let (Some(above), Some(below)) = (order(ev, lo, &z)?, order(ev, &z, hi)?) else {
                return Ok(None);
            };
            if above != Ordering::Less || below != Ordering::Less {
                continue;
            }
            let mut at = points.len();
            for (i, p) in points.iter().enumerate() {
                match order(ev, &z, p)? {
                    Some(Ordering::Greater) => {}
                    Some(Ordering::Less) => {
                        at = i;
                        break;
                    }
                    Some(Ordering::Equal) => {
                        at = usize::MAX;
                        break;
                    }
                    None => return Ok(None),
                }
            }
            if at != usize::MAX {
                points.insert(at, z);
            }
            if points.len() > MAX_POINTS {
                return Ok(None);
            }
        }
    }
    Ok(Some(points))
}

/// `∫ₐᵇ f` split at the zeros of its switches' arguments (the module's
/// rule), for `f` with [`switches`]; `None` when it does not apply.
pub(super) fn integral(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
    [a, b]: [&Expr; 2],
) -> Found {
    let switches = switches(f, x);
    let Some(direction) = order(ev, a, b)? else {
        return Ok(None);
    };
    let (lo, hi) = match direction {
        Ordering::Equal => return Ok(Some(Expr::from(0))),
        Ordering::Less => (a, b),
        Ordering::Greater => (b, a),
    };
    let Some(inside) = points(ev, e, &switches, x, [lo, hi])? else {
        return Ok(None);
    };
    let mut ends = Vec::with_capacity(inside.len() + 2);
    ends.push(lo.clone());
    ends.extend(inside);
    ends.push(hi.clone());
    let mut pieces = Vec::with_capacity(ends.len() - 1);
    for pair in ends.windows(2) {
        let (p, q) = (&pair[0], &pair[1]);
        let middle = match (end(p), end(q)) {
            (Some(false), Some(true)) => Expr::from(0),
            (Some(false), _) => q.clone() - Expr::from(1),
            (_, Some(true)) => p.clone() + Expr::from(1),
            _ => (p.clone() + q.clone()) / Expr::from(2),
        };
        let middle = ev.eval(&middle)?;
        let mut rules = Vec::with_capacity(switches.len());
        for s in &switches {
            let u = &s.args()[0];
            let value = ev.eval(&u.substitute(&[(x.clone(), middle.clone())]))?;
            let positive = match sign(ev, &value) {
                Some(Ordering::Greater) => true,
                Some(Ordering::Less) => false,
                _ => return Ok(None),
            };
            let by = match s.head().and_then(Expr::as_symbol) {
                Some(Symbol::ABS) if positive => u.clone(),
                Some(Symbol::ABS) => -u.clone(),
                Some(Symbol::SIGN) => Expr::from(if positive { 1 } else { -1 }),
                _ => Expr::from(i64::from(positive)),
            };
            rules.push((s.clone(), by));
        }
        let piece = Expr::call(
            Symbol::INTEGRATE,
            [f.substitute(&rules), x.clone(), p.clone(), q.clone()],
        );
        let value = ev.eval(&piece)?;
        if value.contains(&Expr::from(Symbol::INTEGRATE)) {
            return Ok(None);
        }
        pieces.push(value);
    }
    let sum: Expr = pieces.into_iter().sum();
    Ok(Some(if direction == Ordering::Less {
        sum
    } else {
        -sum
    }))
}
