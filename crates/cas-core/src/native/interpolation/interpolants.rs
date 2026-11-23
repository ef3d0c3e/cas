//! Rational and trigonometric interpolation (CAS-594): Xcas's `ratinterp`
//! and `triginterp`, after giac's optimization.cc (`_ratinterp`,
//! `_triginterp`).
//!
//! - **`RationalInterpolation[xs, ys, x0, d]`** (`_ratinterp`), Floater and
//!   Hormann's rational interpolant of `n + 1` points, in barycentric form
//!   `r = Σ w_k y_k/(x − x_k) ÷ Σ w_k/(x − x_k)` with
//!   `w_k = Σ_{i = max(0, k−d)}^{min(k, n−d)} (−1)^i Π_{j = i…i+d, j ≠ k} 1/(x_k − x_j)`,
//!   put together into one rational function.
//!   - The forms are giac's: `(xs, ys[, x0][, d])`, and `(m[, x0][, d])`
//!     with `m` a two-column matrix of points, a first argument that is a
//!     matrix making the matrix form. `x0` is `x` when left out; a list of
//!     points gives `r` at each, another expression `r` there. `d` is
//!     `⌊n/2 + 1⌋` when left out.
//!   - Fewer than two points, a `d` not an integer of `0 … n` (giac's
//!     "Expected an integer 0<=d<=|X|", and its `undef` for `d = n + 1`),
//!     and a repeated abscissa among those a weight takes (giac's `undef`)
//!     leave the call.
//! - **`TrigonometricInterpolation[ys, x == a .. b]` and
//!   `TrigonometricInterpolation[ys, a, b, x]`** (`triginterp`): the
//!   trigonometric polynomial through `n ≥ 2` values at equally spaced
//!   points of `a … b`, of period `T = (b − a) n/(n − 1)`:
//!   `mean(ys) + Σ_{j=1}^{⌊n/2⌋} c_j (A_j cos(2πjx/T) + B_j sin(2πjx/T))`,
//!   `A_j = Σ_k y_k cos(2πj(a/T + k/n))`, `B_j` with `sin`, `c_j = 1/n` for
//!   the last `j` of an even `n` and `2/n` otherwise. Exact on exact input;
//!   a decimal in `ys`, `a` or `b` gives the result numerically, as giac's
//!   `evalf`. One value (giac's "At least two data points are required")
//!   and `a = b` (giac's `undef`) leave the call.
//!
//! Decimals are compared with 0 by their midpoints in doubles, never
//! converted to rationals: a huge decimal would not fit. Each call that
//! answers takes one step, its builtin's name.

use cas_num::Rational;

use super::{MAX_POINTS, Native, is_list, vanishes};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The most `(n + 1)(d + 1)²` of a rational interpolation: the weights'
/// factors, and so the terms put together.
const MAX_RATIONAL_WORK: usize = 1 << 16;

/// The most `n ⌊n/2⌋` of a trigonometric interpolation: its sines and
/// cosines.
const MAX_TRIG_TERMS: usize = 1 << 13;

/// Whether a ball occurs in `x`.
pub(super) fn has_ball(x: &Expr) -> bool {
    matches!(x.view(), View::Ball(_) | View::ComplexBall(_)) || x.args().iter().any(has_ball)
}

/// The rows of a matrix: a nonempty list of lists of one nonzero length.
fn matrix_rows(m: &Expr) -> Option<&[Expr]> {
    let rows = m.args();
    let width = rows.first().filter(|r| is_list(r))?.args().len();
    (is_list(m) && width > 0 && rows.iter().all(|r| is_list(r) && r.args().len() == width))
        .then_some(rows)
}

/// The weights `w_k` with exact abscissas; `None` for a repeated one.
fn rational_weights(xs: &[Rational], d: usize) -> Option<Vec<Expr>> {
    let n = xs.len() - 1;
    let mut weights = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let mut w = Rational::from(0);
        for i in k.saturating_sub(d)..=k.min(n - d) {
            let mut p = Rational::from(if i % 2 == 0 { 1 } else { -1 });
            for j in (i..=i + d).filter(|&j| j != k) {
                let gap = &xs[k] - &xs[j];
                if gap.is_zero() {
                    return None;
                }
                p /= gap;
            }
            w += p;
        }
        weights.push(Expr::from(w));
    }
    Some(weights)
}

/// The weights `w_k`, evaluated; `None` for a repeated abscissa.
fn weights(ev: &mut Eval<'_>, xs: &[Expr], d: usize) -> Result<Option<Vec<Expr>>, EvalError> {
    if let Some(exact) = xs.iter().map(Expr::to_rational).collect::<Option<Vec<_>>>() {
        return Ok(rational_weights(&exact, d));
    }
    let n = xs.len() - 1;
    let mut weights = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let mut terms = Vec::new();
        for i in k.saturating_sub(d)..=k.min(n - d) {
            let mut factors = vec![Expr::from(if i % 2 == 0 { 1 } else { -1 })];
            for j in (i..=i + d).filter(|&j| j != k) {
                let gap = ev.eval(&(xs[k].clone() - xs[j].clone()))?;
                if vanishes(&gap) {
                    return Ok(None);
                }
                factors.push(gap.pow(-1));
            }
            terms.push(Expr::call(Symbol::TIMES, factors));
        }
        weights.push(ev.eval(&Expr::call(Symbol::PLUS, terms))?);
    }
    Ok(Some(weights))
}

pub(in crate::native) fn rational_interpolation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    let Some(first) = args.first().filter(|a| is_list(a)) else {
        return Ok(None);
    };
    let (xs, ys, options): (Vec<Expr>, Vec<Expr>, _) = if let Some(rows) = matrix_rows(first) {
        if rows[0].args().len() != 2 {
            return Ok(None);
        }
        let xs = rows.iter().map(|r| r.args()[0].clone()).collect();
        let ys = rows.iter().map(|r| r.args()[1].clone()).collect();
        (xs, ys, &args[1..])
    } else {
        match args.get(1) {
            Some(ys) if is_list(ys) && ys.args().len() == first.args().len() => {
                (first.args().to_vec(), ys.args().to_vec(), &args[2..])
            }
            _ => return Ok(None),
        }
    };
    if !(2..=MAX_POINTS).contains(&xs.len()) || options.len() > 2 {
        return Ok(None);
    }
    let n = xs.len() - 1;
    let d = match options.get(1) {
        None => n / 2 + 1,
        Some(d) => match d.as_i64().and_then(|d| usize::try_from(d).ok()) {
            // d = n + 1 is giac's undef: no weight has a term.
            Some(d) if d <= n => d,
            _ => return Ok(None),
        },
    };
    if (n + 1) * (d + 1) * (d + 1) > MAX_RATIONAL_WORK {
        return Ok(None);
    }
    let Some(w) = weights(ev, &xs, d)? else {
        return Ok(None);
    };
    let x0 = options
        .first()
        .cloned()
        .unwrap_or_else(|| Expr::symbol("x"));
    let t = if x0.as_symbol().is_some() {
        x0.clone()
    } else {
        Expr::symbol("x")
    };
    let mut numerator = Vec::with_capacity(n + 1);
    let mut denominator = Vec::with_capacity(n + 1);
    for ((xk, yk), wk) in xs.iter().zip(&ys).zip(w) {
        let pole = (t.clone() - xk.clone()).pow(-1);
        numerator.push(wk.clone() * yk.clone() * pole.clone());
        denominator.push(wk * pole);
    }
    let r = Expr::call(Symbol::PLUS, numerator) / Expr::call(Symbol::PLUS, denominator);
    let r = ev.eval(&Expr::call(Symbol::TOGETHER, [r]))?;
    let at = |p: &Expr| Expr::call(Symbol::TOGETHER, [r.substitute(&[(t.clone(), p.clone())])]);
    let result = if x0.as_symbol().is_some() {
        r.clone()
    } else if is_list(&x0) {
        let values = x0.args().iter().map(at).collect::<Vec<_>>();
        ev.eval(&Expr::call(Symbol::LIST, values))?
    } else {
        ev.eval(&at(&x0))?
    };
    ev.step("RationalInterpolation", e, &result);
    Ok(Some(result))
}

pub(in crate::native) fn trigonometric_interpolation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (ys, a, b, x) = match e.args() {
        [ys, eq] => match (eq.head().and_then(Expr::as_symbol), eq.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL), [x, r]) => {
                match (r.head().and_then(Expr::as_symbol), r.args()) {
                    (Some(Symbol::RANGE), [a, b]) => (ys, a, b, x),
                    _ => return Ok(None),
                }
            }
            _ => return Ok(None),
        },
        [ys, a, b, x] => (ys, a, b, x),
        _ => return Ok(None),
    };
    if !is_list(ys) || x.as_symbol().is_none() {
        return Ok(None);
    }
    let ys = ys.args();
    let n = ys.len();
    let half = n / 2;
    if !(2..=MAX_POINTS).contains(&n) || n * half > MAX_TRIG_TERMS {
        return Ok(None);
    }
    let count = |m: usize| Expr::from(i64::try_from(m).expect("a length in memory"));
    let width = ev.eval(&(b.clone() - a.clone()))?;
    if vanishes(&width) {
        return Ok(None);
    }
    let period = ev.eval(&(width * count(n) / count(n - 1)))?;
    let pi = Expr::from(Symbol::PI);
    let mut terms = vec![Expr::call(Symbol::PLUS, ys.to_vec()) / count(n)];
    for j in 1..=half {
        let c = Expr::from(if n % 2 == 0 && j == half { 1 } else { 2 }) / count(n);
        let mut cosines = Vec::with_capacity(n);
        let mut sines = Vec::with_capacity(n);
        for (k, y) in ys.iter().enumerate() {
            let angle = Expr::from(2)
                * pi.clone()
                * count(j)
                * (a.clone() / period.clone() + count(k) / count(n));
            let angle = ev.eval(&angle)?;
            cosines.push(y.clone() * Expr::call(Symbol::COS, [angle.clone()]));
            sines.push(y.clone() * Expr::call(Symbol::SIN, [angle]));
        }
        let along = Expr::from(2) * pi.clone() * count(j) * x.clone() / period.clone();
        let along = ev.eval(&along)?;
        let cosine = ev.eval(&(c.clone() * Expr::call(Symbol::PLUS, cosines)))?;
        let sine = ev.eval(&(c * Expr::call(Symbol::PLUS, sines)))?;
        terms.push(cosine * Expr::call(Symbol::COS, [along.clone()]));
        terms.push(sine * Expr::call(Symbol::SIN, [along]));
    }
    let mut result = ev.eval(&Expr::call(Symbol::PLUS, terms))?;
    if ys.iter().any(has_ball) || has_ball(a) || has_ball(b) {
        result = ev.eval(&Expr::call(Symbol::N, [result]))?;
    }
    ev.step("TrigonometricInterpolation", e, &result);
    Ok(Some(result))
}
