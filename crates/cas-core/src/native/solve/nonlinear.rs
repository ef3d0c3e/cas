//! Two polynomial equations in two unknowns (CAS-DOC-11 §6, CAS-215):
//! `Solve[{f = 0, g = 0}, {x, y}]`, by the resultant in `y`.
//!
//! - The equations converted together, polynomials in `x` and `y` and
//!   parameters free of both, over ℤ: their numerators.
//! - `R(x) = res_y(f, g)` eliminates `y`, and is solved for `x` as an
//!   equation in one unknown (§2, §3). `R = 0` when `f` and `g` share a
//!   factor that has `y`, and so infinitely many solutions: the call stays.
//! - At each root `a`, the equation of the lower degree in `y` that still
//!   has `y` there, of one degree the one whose leading coefficient in `y`
//!   is free of `x`, is solved for `y`, and each `(a, b)` kept when both
//!   equations hold at it, their values 0 or numbers whose balls hold 0
//!   (the check of CAS-213): `[x^2+y^2=1,x-y=0]` gives `±(√2/2, √2/2)`.
//! - Each pair once; when all are numbers, by `x` then `y`, ascending.

use cas_num::algebra::Integers;
use cas_num::ball::Float;
use cas_poly::{Cancel, SparseRing};

use super::super::approximation;
use super::super::polynomial::{Dense, convert, within_dense};
use super::transcendental::holds;
use super::{BITS, roots, sides};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::Fraction;
use crate::symbol::Symbol;

/// The solutions of the two `equations` in the two `unknowns`, as
/// `{{a, b}, …}`; `None` when they are not two polynomial equations, or
/// their solutions are not found here.
pub(super) fn system(
    ev: &mut Eval<'_>,
    e: &Expr,
    equations: &[Expr],
    unknowns: &[Expr],
    real: bool,
) -> Result<Option<Expr>, EvalError> {
    let ([f, g], [x, y]) = (equations, unknowns) else {
        return Ok(None);
    };
    let mut fs = Vec::with_capacity(2);
    for eq in [f, g] {
        fs.push(match sides(eq) {
            Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
            None => eq.clone(),
        });
    }
    let Some(p) = convert(ev, e, &fs, false)? else {
        return Ok(None);
    };
    let (Some(kx), Some(ky)) = (p.variable(x), p.variable(y)) else {
        return Ok(None);
    };
    let inside = p
        .generators()
        .iter()
        .enumerate()
        .any(|(j, h)| j != kx && j != ky && (h.contains(x) || h.contains(y)));
    if inside || within_dense(ev, e, &p, Dense::All).is_none() {
        return Ok(None);
    }
    let [
        Fraction::Integers { num: a, .. },
        Fraction::Integers { num: b, .. },
    ] = p.fractions()
    else {
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, p.generators().len());
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let resultant = ring
        .resultant_with(a, b, ky, cancel)
        .map_err(|_| EvalError::Interrupted)?;
    if resultant.is_zero() {
        return Ok(None);
    }
    let r = ev.eval(&p.polynomial_expr(&resultant))?;
    let Some(xs) = roots(ev, e, &r, x, real)? else {
        return Ok(None);
    };
    // The equations by their degree in y, the lower first, and of one
    // degree first the one whose leading coefficient in y is free of x,
    // which y is then not divided by: x + y = 3 before x y = 1.
    let mut order: Vec<((u64, bool), &Expr)> = [(a, &fs[0]), (b, &fs[1])]
        .into_iter()
        .map(|(n, f)| {
            let lead = ring.lead_coeff(n, ky);
            let by_x = ring.degree(&lead, kx).unwrap_or(0) > 0;
            ((ring.degree(n, ky).unwrap_or(0), by_x), f)
        })
        .collect();
    order.sort_by_key(|&(key, _)| key);
    let mut pairs: Vec<(Expr, Expr)> = Vec::new();
    for xa in xs {
        let Some(ys) = y_roots(ev, e, &order, x, &xa, y, real)? else {
            return Ok(None);
        };
        for yb in ys {
            let rules = [(x.clone(), xa.clone()), (y.clone(), yb.clone())];
            let ctx = ev.context();
            if fs.iter().all(|f| holds(ctx, f, &rules))
                && !pairs.contains(&(xa.clone(), yb.clone()))
            {
                pairs.push((xa.clone(), yb));
            }
        }
    }
    let pairs = ordered(ev.context(), pairs);
    Ok(Some(Expr::call(
        Symbol::LIST,
        pairs
            .into_iter()
            .map(|(a, b)| Expr::call(Symbol::LIST, [a, b])),
    )))
}

/// The roots in `y` at `x = xa` of the first equation of `order` that still
/// has `y` there; none when neither has.
fn y_roots(
    ev: &mut Eval<'_>,
    e: &Expr,
    order: &[((u64, bool), &Expr)],
    x: &Expr,
    xa: &Expr,
    y: &Expr,
    real: bool,
) -> Result<Option<Vec<Expr>>, EvalError> {
    for &(_, f) in order {
        let h = ev.eval(&f.substitute(&[(x.clone(), xa.clone())]))?;
        if !h.contains(y) {
            continue;
        }
        return roots(ev, e, &h, y, real);
    }
    // Both free of y at xa: every y, which these methods do not list.
    Ok(None)
}

/// The pairs by `x` then `y` ascending when all are real numbers;
/// otherwise as they are.
fn ordered(ctx: &Context, pairs: Vec<(Expr, Expr)>) -> Vec<(Expr, Expr)> {
    let mid = |e: &Expr| -> Option<Float> {
        match approximation(ctx, e, BITS).ok()?.view() {
            View::Ball(b) if b.is_finite() => Some(b.mid().clone()),
            _ => None,
        }
    };
    let keys: Option<Vec<(Float, Float)>> = pairs
        .iter()
        .map(|(a, b)| Some((mid(a)?, mid(b)?)))
        .collect();
    let Some(keys) = keys else {
        return pairs;
    };
    let mut keyed: Vec<((Float, Float), (Expr, Expr))> = keys.into_iter().zip(pairs).collect();
    keyed.sort_by(|(p, _), (q, _)| p.cmp(q));
    keyed.into_iter().map(|(_, pair)| pair).collect()
}
