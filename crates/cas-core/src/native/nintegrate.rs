//! `NIntegrate` (CAS-DOC-26 §4, CAS-324): Xcas's `nInt`, `romberg` and
//! `gaussquad`, a ball that contains `∫ₐᵇ f(x) dx` for real bounds, at the
//! working precision `p`, by adaptive Gauss–Legendre quadrature with a
//! rigorous error bound, as Arb's `acb_calc_integrate` \[Joh18\]:
//!
//! - The integrand is compiled once into ball operations (`analytic`,
//!   §3); a part that is no number, as a parameter, or a head without a
//!   ball function leaves the call.
//! - Bounds that are balls integrate between their midpoints, and their
//!   radii add `r |f|` over each end's ball to the result's radius.
//! - The segments start as `[a, b]`, each with its share of the tolerance
//!   `ε = 2^(−p) max(1, |I₀|)` by length, `I₀` an 8-point estimate. A
//!   segment of middle `m` and half-length `h` is tried in turn:
//!   1. for `ρ` in [`RHOS`], `f(m + h t)` on the box `[−A, A] + i[−B, B]`
//!      around the Bernstein ellipse `E_ρ`, `A = (ρ + 1/ρ)/2`,
//!      `B = (ρ − 1/ρ)/2`, in the analytic mode, its switches continued
//!      from the segment: its modulus bounds `M`, and the `n`-point rule
//!      errs by at most `|h| 64M/(15(1 − ρ⁻²)ρ^(2n))` \[Tre08\]; the least
//!      `n` of [`POINTS`] within the share, over the `ρ`s, sums
//!      `h Σ w_k f(m + h x_k)`, the bound added to its radius;
//!   2. the crude enclosure `2h f([m − h, m + h])` in the real mode, when
//!      its radius is within `ε √(2h/(b − a))`: a jump or a kink, isolated
//!      in a short segment. A jump's segment encloses `2h [−1, 1]` times
//!      the jump, which a share in proportion to `2h` would never accept;
//!      the shares of `√(2h)` along the halves towards a point add up to a
//!      few `ε`;
//!   3. otherwise the halves, split exactly, so that segments go below the
//!      working precision, within [`MAX_SEGMENTS`] segments and
//!      [`MAX_EVALUATIONS`] evaluations; past them the call stays, with a
//!      warning.
//! - The result is the sum of the segments' balls, real where each is.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of quadrature: f in x from a to b, m ± h, n points, weights w"
)]

mod analytic;

use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};

use cas_num::Integer;
use cas_num::Rational;
use cas_num::ball::{Ball, ComplexBall, Float, Mag, gauss_legendre};

use self::analytic::{Mode, Node};
use super::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The point counts of the rules tried.
const POINTS: [usize; 7] = [4, 8, 16, 32, 64, 128, 256];

/// The ellipses tried, `ρ = p/q`, from the widest.
const RHOS: [(i64, i64); 5] = [(8, 1), (4, 1), (2, 1), (3, 2), (5, 4)];

/// The segments looked at, at most.
pub(super) const MAX_SEGMENTS: usize = 1024;

/// The evaluations of the integrand, at most.
pub(super) const MAX_EVALUATIONS: usize = 1 << 17;

/// The guard bits of the evaluations.
const GUARD: u64 = 32;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// `f`, `x` and the bounds of `NIntegrate[f, x, a, b]` or
/// `NIntegrate[f, {x, a, b}]`, `x` a symbol.
fn parts(e: &Expr) -> Option<(&Expr, &Expr, [&Expr; 2])> {
    let (f, x, a, b) = match e.args() {
        [f, x, a, b] => (f, x, a, b),
        [f, list] if list.head() == Some(&Expr::from(Symbol::LIST)) => match list.args() {
            [x, a, b] => (f, x, a, b),
            _ => return None,
        },
        _ => return None,
    };
    x.as_symbol()?;
    Some((f, x, [a, b]))
}

/// A real bound as a ball at `prec` bits; `None` for another.
fn bound(ev: &Eval<'_>, a: &Expr, prec: u64) -> Option<Ball> {
    let value = match a.view() {
        View::Ball(b) => return Some(b.clone()),
        _ => approximation(ev.context(), a, prec).ok()?,
    };
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

/// An upper bound of `|b|`.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(&b.mid().abs()).add(b.rad())
}

/// The ball `m ± r`.
fn around(m: &Ball, r: &Ball) -> Ball {
    m.add_error(&upper(r))
}

/// `z` widened by `error`: in its real part alone when it is real, the
/// integral of a function real on the segment.
fn widen(z: &ComplexBall, error: &Mag) -> ComplexBall {
    if z.is_real() {
        ComplexBall::from_real(z.re().add_error(error))
    } else {
        ComplexBall::new(z.re().add_error(error), z.im().add_error(error))
    }
}

pub(super) fn nintegrate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x, [a, b])) = parts(e) else {
        return Ok(None);
    };
    let prec = ev.context().working_precision();
    let wp = prec + GUARD;
    let (Some(a), Some(b)) = (bound(ev, a, wp), bound(ev, b, wp)) else {
        return Ok(None);
    };
    let Some(tree) = Node::compile(ev.context(), f, x, wp) else {
        return Ok(None);
    };
    let Some(value) = integral(&tree, &a, &b, prec, ev.context().interrupt_flag())? else {
        ev.warn(Warning::Native {
            message: "the integral is not enclosed within the budget of segments".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    ev.step("NIntegrate", e, &Expr::from(&value));
    Ok(Some(if value.is_real() {
        Expr::from(value.re())
    } else {
        Expr::from(&value)
    }))
}

/// What a segment's integral cost.
struct Budget {
    segments: usize,
    evaluations: usize,
}

/// `∫ₐᵇ f` in a ball, the module's rule; `None` past the budget or where
/// an end's value is not finite.
fn integral(
    tree: &Node,
    a: &Ball,
    b: &Ball,
    prec: u64,
    interrupt: &AtomicBool,
) -> Result<Option<ComplexBall>, EvalError> {
    let wp = prec + GUARD;
    let (lo, hi) = (Ball::exact(a.mid().clone()), Ball::exact(b.mid().clone()));
    let mut total = ComplexBall::ZERO;
    // The ends' radii: |∫ between an end and its midpoint| ≤ r sup |f|.
    for end in [a, b] {
        if !end.rad().is_zero() {
            let Some(value) = tree.eval(&ComplexBall::from_real(end.clone()), Mode::Real, wp)
            else {
                return Ok(None);
            };
            let error = end.rad().mul(&upper(&value.abs(wp)));
            total = widen(&total, &error);
        }
    }
    let length = hi.sub(&lo, wp);
    if length.mid().is_zero() {
        return Ok(Some(total));
    }
    let coarse = gauss_sum(tree, &lo, &hi, 8, wp)
        .map_or(Mag::from_integer(&Integer::from(1)), |v| {
            upper(&v.abs(wp)).max(Mag::from_integer(&Integer::from(1)))
        });
    let Ok(bits) = i64::try_from(prec) else {
        return Ok(None);
    };
    let tolerance = coarse.mul_2exp(&-Integer::from(bits));
    let length = upper(&length);
    let mut budget = Budget {
        segments: 0,
        evaluations: 0,
    };
    let mut work = vec![(lo, hi)];
    while let Some((lo, hi)) = work.pop() {
        budget.segments += 1;
        if budget.segments > MAX_SEGMENTS || budget.evaluations > MAX_EVALUATIONS {
            return Ok(None);
        }
        if interrupt.load(AtomicOrdering::Relaxed) {
            return Err(EvalError::Interrupted);
        }
        let width = upper(&hi.sub(&lo, wp));
        let share = tolerance.mul(&width).div(&length);
        if let Some(v) = gauss(tree, &lo, &hi, &share, wp, &mut budget) {
            total = total.add(&v, wp);
            continue;
        }
        // A jump's segment encloses `width · [−1, 1]`: within a share of
        // `√width`, it is accepted once short enough.
        let crude_share = tolerance.mul(&width.div(&length).sqrt());
        if let Some(v) = crude(tree, &lo, &hi, &crude_share, wp) {
            total = total.add(&v, wp);
            continue;
        }
        // Exactly, so that segments go below the working precision.
        let Some(sum) = lo.mid().checked_add(hi.mid()) else {
            return Ok(None);
        };
        let middle = Ball::exact(sum.mul_2exp(&-Integer::from(1)));
        work.push((middle.clone(), hi));
        work.push((lo, middle));
    }
    Ok(Some(total))
}

/// The middle and the half-length of `[lo, hi]`.
fn middle_and_half(lo: &Ball, hi: &Ball, wp: u64) -> (Ball, Ball) {
    let half = Integer::from(-1);
    (
        lo.add(hi, wp).mul_2exp(&half),
        hi.sub(lo, wp).mul_2exp(&half),
    )
}

/// `h Σ w_k f(m + h x_k)` for the `n`-point rule in the analytic mode of
/// `tree`, which has no switch left; `None` where an evaluation fails.
fn gauss_sum(tree: &Node, lo: &Ball, hi: &Ball, n: usize, wp: u64) -> Option<ComplexBall> {
    let (m, h) = middle_and_half(lo, hi, wp);
    let rule = gauss_legendre(n, wp);
    let mut sum = ComplexBall::ZERO;
    for (x, w) in rule.iter() {
        let point = ComplexBall::from_real(m.add(&h.mul(x, wp), wp));
        let value = tree.eval(&point, Mode::Analytic, wp)?;
        sum = sum.add(&value.mul(&ComplexBall::from_real(w.clone()), wp), wp);
    }
    Some(sum.mul(&ComplexBall::from_real(h), wp))
}

/// The Gauss–Legendre value of a segment with its error bound, the least
/// rule within `share` over the ellipses; `None` when none is.
fn gauss(
    tree: &Node,
    lo: &Ball,
    hi: &Ball,
    share: &Mag,
    wp: u64,
    budget: &mut Budget,
) -> Option<ComplexBall> {
    let (m, h) = middle_and_half(lo, hi, wp);
    let segment = ComplexBall::from_real(around(&m, &h));
    let resolved;
    let tree = if tree.has_switch() {
        resolved = tree.resolved(&segment, wp)?;
        &resolved
    } else {
        tree
    };
    let one = Ball::from_integer(&Integer::from(1));
    let mut best: Option<(usize, Mag)> = None;
    for (p, d) in RHOS {
        let rho = Ball::from_rational(&q(p, d), wp);
        let inverse = rho.inv(wp);
        let a = rho.add(&inverse, wp).mul_2exp(&Integer::from(-1));
        let b = rho.sub(&inverse, wp).mul_2exp(&Integer::from(-1));
        let disk = ComplexBall::new(
            Ball::new(Float::ZERO, upper(&a)),
            Ball::new(Float::ZERO, upper(&b)),
        );
        let z = ComplexBall::from_real(m.clone())
            .add(&disk.mul(&ComplexBall::from_real(h.clone()), wp), wp);
        budget.evaluations += 1;
        let Some(value) = tree.eval(&z, Mode::Analytic, wp) else {
            continue;
        };
        // |h| 64 M / (15 (1 − ρ⁻²) ρ^(2n)).
        let factor = Ball::from_rational(&q(64, 15), wp)
            .mul(&value.abs(wp), wp)
            .mul(&h.abs(), wp)
            .div(&one.sub(&inverse.mul(&inverse, wp), wp), wp);
        for &n in &POINTS {
            if best.as_ref().is_some_and(|(k, _)| *k <= n) {
                break;
            }
            let decay = rho.pow(i64::try_from(2 * n).expect("a small count"), wp);
            let error = upper(&factor.div(&decay, wp));
            if error <= *share {
                best = Some((n, error));
                break;
            }
        }
    }
    let (n, error) = best?;
    budget.evaluations += n;
    gauss_sum(tree, lo, hi, n, wp).map(|v| widen(&v, &error))
}

/// `2h f([lo, hi])` in the real mode, when its radius is within `share`.
fn crude(tree: &Node, lo: &Ball, hi: &Ball, share: &Mag, wp: u64) -> Option<ComplexBall> {
    let (m, h) = middle_and_half(lo, hi, wp);
    let segment = ComplexBall::from_real(around(&m, &h));
    let value = tree.eval(&segment, Mode::Real, wp)?;
    let width = ComplexBall::from_real(hi.sub(lo, wp));
    let enclosed = value.mul(&width, wp);
    (*enclosed.re().rad() <= *share && *enclosed.im().rad() <= *share).then_some(enclosed)
}

#[cfg(test)]
mod tests;
