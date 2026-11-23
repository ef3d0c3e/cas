//! Sampling and intervals (CAS-593): Xcas's `linear_interpolate`,
//! `interval2center` and `center2interval`, after giac's misc.cc.
//!
//! - **`LinearInterpolate[m, xmin, xmax, xstep]`** (`gen_interpolate` of
//!   degree 1): `m` the two rows `[xs, ys]` of points, sorted by abscissa;
//!   the polygonal line through them sampled at `xmin, xmin + |xstep|, …`
//!   within `xmax`, answered as `[[x…], [y…]]`. Exact on exact input, where
//!   giac answers doubles (`19/3` for its `6.33333333333`).
//!   - A range past the points (giac's "Values out of range"), a step 0
//!     (giac's timeout), `xmax < xmin`, a repeated abscissa, an entry no
//!     real number, or more than [`MAX_ELEMENTS`] samples leave the call.
//! - **`IntervalToCenter[r]`** (`_interval2center`): a range `a .. b` gives
//!   `(a + b)/2`, a decimal its midpoint, a list each element's centre, and
//!   anything else itself, as giac's `interval2center(3)` is `3`.
//! - **`CenterToInterval[cs]` and `CenterToInterval[cs, c]`**
//!   (`centres2intervalles`): the classes of the real centres `m_0, m_1, …`
//!   as a sequence of ranges. The first starts at `c`, by default
//!   `m_0 + (m_0 − m_1)/2`; each ends at the midpoint of its centre and the
//!   next, the last as far past its centre as it starts before it. Fewer
//!   than two centres, one not real, or a class that does not end past its
//!   start leave the call (giac's "Bad Argument Value").
//!
//! Decimals are compared by their midpoints in doubles, never converted to
//! rationals: a huge decimal would not fit.

use std::cmp::Ordering;

use cas_num::Rational;
use cas_num::ball::{Ball, Float, Round};

use super::super::lists::MAX_ELEMENTS;
use super::{Native, is_list};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// A real number: exact, or a decimal's midpoint in doubles.
enum Real {
    Exact(Rational),
    Double(f64),
}

impl Real {
    fn of(e: &Expr) -> Option<Real> {
        if let Some(q) = e.to_rational() {
            return Some(Real::Exact(q));
        }
        match e.view() {
            View::Ball(b) if b.is_finite() => Some(Real::Double(b.mid().to_f64() + 0.0)),
            _ => None,
        }
    }

    fn double(&self) -> f64 {
        match self {
            Real::Exact(q) => Float::from_rational(q, 64, Round::Nearest).0.to_f64() + 0.0,
            Real::Double(x) => *x,
        }
    }

    fn cmp(&self, other: &Real) -> Ordering {
        match (self, other) {
            (Real::Exact(a), Real::Exact(b)) => a.cmp(b),
            _ => self.double().total_cmp(&other.double()),
        }
    }
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// The ends of `Range[a, b]`, the language's `a .. b`.
fn range(r: &Expr) -> Option<(&Expr, &Expr)> {
    match (r.head().and_then(Expr::as_symbol), r.args()) {
        (Some(Symbol::RANGE), [a, b]) => Some((a, b)),
        _ => None,
    }
}

/// Whether `lo, lo + step, …` within `hi` has at most [`MAX_ELEMENTS`]
/// values, and at least one; `step > 0`.
fn few_samples(lo: &Real, hi: &Real, step: &Real) -> bool {
    if hi.cmp(lo) == Ordering::Less {
        return false;
    }
    if let (Real::Exact(lo), Real::Exact(hi), Real::Exact(step)) = (lo, hi, step) {
        return (&(hi - lo) / step)
            .floor()
            .to_i64()
            .is_some_and(|n| n < MAX_ELEMENTS);
    }
    let bound = i32::try_from(MAX_ELEMENTS).map_or(f64::MAX, f64::from);
    let n = (hi.double() - lo.double()) / step.double();
    n.is_finite() && n < bound
}

pub(in crate::native) fn linear_interpolate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, xmin, xmax, xstep] = e.args() else {
        return Ok(None);
    };
    let [xs, ys] = m.args() else {
        return Ok(None);
    };
    if !is_list(m) || !is_list(xs) || !is_list(ys) || xs.args().len() != ys.args().len() {
        return Ok(None);
    }
    let mut points = Vec::with_capacity(xs.args().len());
    for (x, y) in xs.args().iter().zip(ys.args()) {
        let (Some(at), Some(_)) = (Real::of(x), Real::of(y)) else {
            return Ok(None);
        };
        points.push((at, x.clone(), y.clone()));
    }
    points.sort_by(|a, b| a.0.cmp(&b.0));
    if points.len() < 2
        || points
            .windows(2)
            .any(|w| w[0].0.cmp(&w[1].0) == Ordering::Equal)
    {
        return Ok(None);
    }
    let step = ev.eval(&Expr::call(Symbol::ABS, [xstep.clone()]))?;
    let (Some(lo), Some(hi), Some(by)) = (Real::of(xmin), Real::of(xmax), Real::of(&step)) else {
        return Ok(None);
    };
    let last = points.len() - 1;
    if by.double() <= 0.0
        || lo.cmp(&points[0].0) == Ordering::Less
        || hi.cmp(&points[last].0) == Ordering::Greater
        || !few_samples(&lo, &hi, &by)
    {
        return Ok(None);
    }
    let (mut xs_out, mut ys_out) = (Vec::new(), Vec::new());
    let mut k = 0;
    for i in 0..MAX_ELEMENTS {
        let x = ev.eval(&(xmin.clone() + Expr::from(i) * step.clone()))?;
        let Some(at) = Real::of(&x) else {
            return Ok(None);
        };
        if at.cmp(&hi) == Ordering::Greater {
            break;
        }
        while k < last && points[k + 1].0.cmp(&at) != Ordering::Greater {
            k += 1;
        }
        let y = if k == last {
            points[last].2.clone()
        } else {
            let ((_, x0, y0), (_, x1, y1)) = (&points[k], &points[k + 1]);
            ev.eval(
                &(y0.clone()
                    + (x.clone() - x0.clone()) * (y1.clone() - y0.clone())
                        / (x1.clone() - x0.clone())),
            )?
        };
        xs_out.push(x);
        ys_out.push(y);
    }
    let result = list(vec![list(xs_out), list(ys_out)]);
    ev.step("LinearInterpolate", e, &result);
    Ok(Some(result))
}

/// The centre of a range, a decimal or each element of a list.
fn centre(ev: &mut Eval<'_>, r: &Expr) -> Result<Expr, EvalError> {
    if is_list(r) {
        let mut items = Vec::with_capacity(r.args().len());
        for x in r.args() {
            items.push(centre(ev, x)?);
        }
        return Ok(list(items));
    }
    if let Some((a, b)) = range(r) {
        return ev.eval(&((a.clone() + b.clone()) / Expr::from(2)));
    }
    if let View::Ball(b) = r.view() {
        return Ok(Expr::from(Ball::exact(b.mid().clone())));
    }
    Ok(r.clone())
}

pub(in crate::native) fn interval_to_center(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [r] = e.args() else {
        return Ok(None);
    };
    let result = centre(ev, r)?;
    ev.step("IntervalToCenter", e, &result);
    Ok(Some(result))
}

pub(in crate::native) fn center_to_interval(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (cs, c) = match e.args() {
        [cs] => (cs, None),
        [cs, c] => (cs, Some(c)),
        _ => return Ok(None),
    };
    let centres = cs.args();
    if !is_list(cs) || centres.len() < 2 || centres.iter().any(|m| Real::of(m).is_none()) {
        return Ok(None);
    }
    let mut start = if let Some(c) = c {
        c.clone()
    } else {
        let (m0, m1) = (centres[0].clone(), centres[1].clone());
        ev.eval(&(m0.clone() + (m0 - m1) / Expr::from(2)))?
    };
    let Some(mut from) = Real::of(&start) else {
        return Ok(None);
    };
    let mut classes = Vec::with_capacity(centres.len());
    for (k, m) in centres.iter().enumerate() {
        let end = match centres.get(k + 1) {
            Some(next) => ev.eval(&((m.clone() + next.clone()) / Expr::from(2)))?,
            None => ev.eval(&(m.clone() + (m.clone() - start.clone())))?,
        };
        let Some(to) = Real::of(&end) else {
            return Ok(None);
        };
        if to.cmp(&from) != Ordering::Greater {
            return Ok(None);
        }
        classes.push(Expr::call(Symbol::RANGE, [start, end.clone()]));
        (start, from) = (end, to);
    }
    let result = Expr::call(Symbol::SEQUENCE, classes);
    ev.step("CenterToInterval", e, &result);
    Ok(Some(result))
}
