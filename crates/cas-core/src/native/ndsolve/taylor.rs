//! Taylor's method (CAS-DOC-35 §4, CAS-392) for `y′ = f(t, y)`, `y` a
//! vector, from `(t_k, Y_k)`, `Y_k` balls:
//!
//! 1. **The Taylor coefficients** at `(t_k, Y_k)` by automatic
//!    differentiation: `y_0 = Y_k` and `y_(j+1) = F_j/(j + 1)`, `F_j` the
//!    `j`-th coefficient of `f(t_k + s, y(t_k + s))`, which the jets of `y`
//!    to `j` give ([`coefficients`]).
//! 2. **The a priori enclosure** over `T = t_k + [0, h]`: a box `B` with
//!    `Y_k + [0, h] f(T, B) ⊆ B`. By Picard and Lindelöf the solution exists
//!    on `T` and stays in `B`; then `Y_k + [0, h] f(T, B)`, tighter, holds
//!    it too. `B` starts from `Y_k` and that box widened, at most
//!    [`ENCLOSURE_TRIES`] times.
//! 3. **The remainder**: the `N`-th coefficient of the solution over `(T,
//!    B)`, by the same recurrence on those balls. By Taylor's theorem with
//!    Lagrange's remainder, `y(t_k + h) ∈ Σ_{j<N} y_j h^j + h^N y_N(T, B)`.
//! 4. **The step**: `h` takes the last coefficient `y_(N−1)` to about
//!    `2^-p` of the solution, `|h| ≈ (2^-p max(1, |y|)/|y_(N−1)|)^(1/(N−1))`,
//!    a power of 2, at most what is left to `t1`; it halves where the
//!    enclosure fails or the remainder is wider than `2^(8−p)` of the
//!    solution.
//!
//! `t` and `h` are balls, so that `t1` need not be dyadic, and `h` may be
//! negative. Past [`MAX_STEPS`] steps, or a step shorter than `2^-p` of
//! `|t1 − t0|`, the integration fails: a solution that blows up, or a stiff
//! problem.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the method: y′ = f(t, y), steps h, order n, bits q"
)]

use cas_num::Integer;
use cas_num::ball::{Ball, Float, Jet, Mag};

use super::rhs::Node;

/// The points of the steps: `(t, y)`.
pub(super) type Curve = Vec<(Ball, Vec<Ball>)>;

/// The most steps of one integration.
pub(super) const MAX_STEPS: usize = 4096;

/// The tries of the a priori enclosure at one step length.
const ENCLOSURE_TRIES: usize = 8;

/// The halvings of one step.
const HALVINGS: usize = 40;

/// The order of the method at `prec` bits: about `p/3`, within `[8, 60]`.
pub(super) fn order(prec: u64) -> usize {
    usize::try_from((prec / 3).clamp(8, 60)).expect("small")
}

/// The Taylor coefficients of the solution at `(t0, y0)` to the order
/// `order`, one jet per component; `None` outside a domain.
pub(super) fn coefficients(
    f: &[Node],
    t0: &Ball,
    y0: &[Ball],
    order: usize,
    prec: u64,
) -> Option<Vec<Jet>> {
    let mut ys: Vec<Vec<Ball>> = y0.iter().map(|b| vec![b.clone()]).collect();
    for j in 0..order.saturating_sub(1) {
        let len = j + 1;
        let mut at = vec![Jet::variable(t0.clone(), len)];
        at.extend(ys.iter().map(|c| Jet::new(c.clone())));
        let divisor = Ball::from_integer(&Integer::from(u64::try_from(j + 1).expect("small")));
        let mut next = Vec::with_capacity(f.len());
        for fi in f {
            let jet = fi.jet(&at, prec)?;
            next.push(jet.coeffs()[j].div(&divisor, prec));
        }
        for (c, n) in ys.iter_mut().zip(next) {
            c.push(n);
        }
    }
    Some(ys.into_iter().map(Jet::new).collect())
}

/// The ball of all points between `a` and `b`, two balls.
fn hull(a: &Ball, b: &Ball, prec: u64) -> Ball {
    a.union(b, prec)
}

/// `b` widened: twice its radius plus `2^-p` of its size and of 1.
fn widened(b: &Ball, prec: u64) -> Ball {
    let size = Mag::from_float(b.mid()).add(&Mag::from_integer(&Integer::ONE));
    let extra = size.mul(&Mag::power_of_two(&-Integer::from(prec)));
    b.add_error(&b.rad().add(&extra))
}

/// The enclosure over `T = t + [0, h]` of the solution from `y` at `t`;
/// `None` when not found.
fn enclosure(f: &[Node], t: &Ball, y: &[Ball], h: &Ball, prec: u64) -> Option<(Ball, Vec<Ball>)> {
    let span = hull(&Ball::ZERO, h, prec);
    let times = t.add(&span, prec);
    let advance = |b: &[Ball]| -> Option<Vec<Ball>> {
        let mut at = vec![times.clone()];
        at.extend(b.iter().cloned());
        let mut out = Vec::with_capacity(y.len());
        for (fi, yi) in f.iter().zip(y) {
            out.push(yi.add(&span.mul(&fi.ball(&at, prec)?, prec), prec));
        }
        Some(out)
    };
    let mut b: Vec<Ball> = advance(y)?.iter().map(|x| widened(x, prec)).collect();
    for _ in 0..ENCLOSURE_TRIES {
        let next = advance(&b)?;
        if next.iter().zip(&b).all(|(n, x)| x.contains_ball(n)) {
            return Some((times, next));
        }
        b = next
            .iter()
            .zip(&b)
            .map(|(n, x)| widened(&n.union(x, prec), prec))
            .collect();
    }
    None
}

/// `Σ_j c_j h^j` for the coefficients `c`, by Horner's scheme.
fn horner(c: &[Ball], h: &Ball, prec: u64) -> Ball {
    c.iter()
        .rev()
        .fold(Ball::ZERO, |acc, cj| acc.mul(h, prec).add(cj, prec))
}

/// One step of length `h` from `y` at `t`, at `prec` bits; `None` when the
/// enclosure fails or the step widens the radius past twice the old one
/// and `2^(8−target)` of the size.
#[allow(clippy::too_many_arguments, reason = "the step's state and precisions")]
fn step(
    f: &[Node],
    t: &Ball,
    y: &[Ball],
    taylor: &[Jet],
    h: &Ball,
    n: usize,
    prec: u64,
    target: u64,
) -> Option<Vec<Ball>> {
    let (times, boxed) = enclosure(f, t, y, h, prec)?;
    let remainder = coefficients(f, &times, &boxed, n + 1, prec)?;
    let hn = h.pow(i64::try_from(n).expect("small"), prec);
    let mut out = Vec::with_capacity(y.len());
    for (jet, rem) in taylor.iter().zip(&remainder) {
        let body = horner(&jet.coeffs()[..n], h, prec);
        let value = body.add(&hn.mul(&rem.coeffs()[n], prec), prec);
        if !value.is_finite() {
            return None;
        }
        out.push(value);
    }
    // The step's own error, against the growth of the radius it carries.
    let grown = out.iter().zip(y).any(|(v, old)| {
        let size = Mag::from_float(v.mid()).add(&Mag::from_integer(&Integer::ONE));
        let allowance = size.mul(&Mag::power_of_two(&Integer::from(
            8 - i64::try_from(target).unwrap_or(i64::MAX),
        )));
        *v.rad() > old.rad().add(&allowance).mul_2exp(&Integer::ONE)
    });
    (!grown).then_some(out)
}

/// The step's length from the coefficients: a power of 2 whose `j`-th
/// power times the coefficient `c_j`, for the last two, `j = n − 1` and
/// `n`, is about `2^-p` of the solution; the whole way where both are 0.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "an estimate"
)]
fn length(taylor: &[Jet], n: usize, prec: u64) -> Float {
    let mut log2_h = f64::INFINITY;
    for jet in taylor {
        let size = jet.coeffs()[0].mid().to_f64().abs().max(1.0).log2();
        for (j, c) in jet.coeffs().iter().enumerate().skip(n - 1).take(2) {
            let m = c.mid().to_f64().abs() + c.rad().to_float().map_or(0.0, |r| r.to_f64());
            if m == 0.0 || !m.is_finite() {
                continue;
            }
            // |c| h^j ≤ 2^-p max(1, |y|).
            let bound = (size - prec as f64 - m.log2()) / j as f64;
            log2_h = log2_h.min(bound);
        }
    }
    let e = if log2_h.is_finite() {
        log2_h.floor().clamp(-1e6, 1e6) as i64
    } else {
        1 << 20
    };
    Float::from(1i64).mul_2exp(&Integer::from(e))
}

/// The solution of `y′ = f(t, y)`, `y(t0) = y0`, at `t1`, with the points
/// of the steps; `None` when the integration fails (§4).
pub(super) fn integrate(
    f: &[Node],
    t0: &Ball,
    y0: &[Ball],
    t1: &Ball,
    prec: u64,
) -> Option<(Vec<Ball>, Curve)> {
    let n = order(prec);
    let q = prec + 32;
    let total = t1.sub(t0, q);
    let backward = total.mid().is_negative();
    let shortest = Mag::from_float(total.mid()).mul(&Mag::power_of_two(&-Integer::from(prec)));
    let (mut t, mut y) = (t0.clone(), y0.to_vec());
    let mut points = vec![(t.clone(), y.clone())];
    for _ in 0..MAX_STEPS {
        let left = t1.sub(&t, q);
        if left.mid().is_zero()
            || left.contains_rational(&cas_num::Rational::ZERO) && left.rad() <= &shortest
        {
            return Some((y, points));
        }
        let taylor = coefficients(f, &t, &y, n + 1, q)?;
        let mut h = Ball::exact(length(&taylor, n, prec));
        if backward {
            h = h.neg();
        }
        let mut last = false;
        // Not past t1: the rest of the way when the step reaches it.
        if Mag::from_float(h.mid()) >= Mag::from_float(left.mid()) {
            h = left.clone();
            last = true;
        }
        let mut taken = None;
        for _ in 0..HALVINGS {
            if Mag::from_float(h.mid()) < shortest {
                return None;
            }
            if let Some(next) = step(f, &t, &y, &taylor, &h, n, q, prec) {
                taken = Some(next);
                break;
            }
            h = h.mul_2exp(&-Integer::ONE);
            last = false;
        }
        let next = taken?;
        t = if last { t1.clone() } else { t.add(&h, q) };
        y = next;
        points.push((t.clone(), y.clone()));
        if last {
            return Some((y, points));
        }
    }
    None
}

#[cfg(test)]
mod tests;
