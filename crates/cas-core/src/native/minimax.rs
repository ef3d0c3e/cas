//! Minimax approximation (CAS-599), Xcas's `minimax`, after giac's
//! optimization.cc (`_minimax`, `chebyshev_nodes`, `skew_nodes`).
//!
//! - **`MinimaxApproximation[f, x == a..b, n]`**, `n ≥ 1`: the polynomial
//!   `p` of degree at most `n` that makes `max |f − p|` over `[a, b]`
//!   least; **`MinimaxApproximation[f, x == a..b, {n, m}]`**, `n ≥ 0`,
//!   `m ≥ 1`: the rational function `p/q`, degrees at most `n` and `m`.
//!   The answer is `{E, p}` or `{E, p/q}`, `E` a ball holding `max |f −
//!   p/q|` over `[a, b]` for the answer's own `p/q`.
//! - Remez's second algorithm (E. Remez, 1934; E. W. Cheney, *Introduction
//!   to Approximation Theory*, ch. 3; L. N. Trefethen, *Approximation Theory
//!   and Approximation Practice*, ch. 10), as giac's `_minimax`:
//!   `N = n + m + 2` reference points, first `a`, `b` and the `n + m`
//!   Chebyshev nodes;
//!   on them the levelled system `p(x_i) − (f_i − (−1)^i E')(q(x_i) − 1) +
//!   (−1)^i E = f_i`, `q`'s constant term 1, solved again with `E' = E`
//!   until `E` settles (giac's linearisation of `f − p/q = ±E`, at most 20
//!   rounds; once for a polynomial); then the exchange; the iterate with the
//!   least `max |e|`, `e = f − p/q`, is kept.
//! - Options, giac's: `rand` skews the first reference 25 times, `x ↦
//!   sign(x)|x|^s` or `x ↦ ((x + 1)/2)^s·2 − 1` on `[−1, 1]`, `s` uniform on
//!   `[3/4, 5/4]`, each time from the last, and starts from the one with the
//!   least `|E|`; `rand = k` does that `k` times from the Chebyshev nodes and
//!   keeps the best answer; `limit = k` iterations (100); `threshold = t >
//!   1` stops at `max |e| ≤ t·min |e|` over the reference, giac's rule
//!   (giac's default `t` 1.02).
//!
//! Departures from giac:
//!
//! - Doubles in the Chebyshev basis on `t ∈ [−1, 1]`, `x = c + h t` with
//!   `c`, `h` from the doubles of `a` and `b` (giac: doubles, monomials in
//!   `x`); the inner rounds stop at `|ΔE| ≤ 2^−30 |E|` (giac: `10^−5`).
//! - The exchange (giac: one zero of `e` between consecutive points by
//!   Brent, one extremum by golden section between consecutive zeros) takes
//!   `e` on a grid of [`SUB`] steps between consecutive points of `{a} ∪
//!   reference ∪ {b}`, the largest `|e|` of each run of one sign refined by
//!   golden section, the runs merged keeping the larger, padded with an end
//!   of the interval when fewer than `N`, and the `N` consecutive ones
//!   holding the largest `|e|` with the largest least `|e|`. giac's exchange
//!   loses the alternation on symmetric cases (`|x|`, degree 2, from
//!   Chebyshev's nodes, interpolates: `E = 0`) and on non-smooth `f` at high
//!   degree.
//! - Stops at `max |e| − min |e| ≤ max(2^−20 min |e|, 64 ε max |f|)` over the
//!   reference, then one more iteration; not stopped within the limit, a
//!   warning and the call stays (giac answers its best).
//! - An `e` left with too few runs of one sign, all `|e| ≤ 64 ε max |f|`
//!   (`f` is a `p/q`), has converged.
//! - The answer: the doubles of the coefficients exactly, Chebyshev to
//!   monomials in `x` over the rationals, `q`'s constant term 1 (its first
//!   term not below `2^−40` of `q`'s size on `[a, b]` when the constant is
//!   rounding, as for `1/x`); a term whose largest value over `[a, b]` is
//!   below `E·2^−40` dropped (`E` at least `64 ε max |f|`); each coefficient
//!   an exact ball at the least precision that moves `p` (and `q`) by at most
//!   `E·2^−40` over `[a, b]`.
//! - `E` certified by branch and bound over a cover of `[a, b]`, its ends the
//!   balls of `a` and `b`: a piece `m ± h` is bounded by the mean value form
//!   `e(m) + e'(I)·[−h, h]` when `D[f, x]` compiles, else by `f(I) −
//!   p(I)/q(I)`, `p`, `q` and their derivatives by Taylor shifts at `m` in
//!   balls; the lower bound `L` the largest `|e|` at a centre; pieces below
//!   `L` dropped, the others bisected while the largest upper bound `U >
//!   L(1 + 2^−20)`, at most `min(terms/N², 2^16)` pieces (`terms` of the
//!   context's limits, which also bound `N²·limit·reps`). `q` vanishing on `[a, b]` (giac
//!   answers with a mean error by area and a warning) and a bound not found
//!   leave the call, with a warning.
//! - `rand` draws from the context's generator.

#![allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    reason = "the names of the algorithm: f, p, q, x, t, e, E, n, m, N, a, b, c, h"
)]

use std::cmp::Ordering;

use cas_num::ball::{Ball, ComplexBall, Float, Mag, Round};
use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::native::nintegrate::analytic::{Mode, Node};
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// Bits beyond the working precision for the bound of the error.
const GUARD: u64 = 32;
/// The largest `n + m`.
const MAX_ORDER: usize = 100;
/// The rounds of the levelled system of a rational function.
const ITC_MAX: usize = 20;
/// The skewed references tried at the first iteration, giac's.
const TRIES: usize = 25;
/// Grid steps between consecutive reference points in the exchange.
const SUB: usize = 16;
/// The most pieces of the cover.
const MAX_PIECES: usize = 1 << 16;

/// The options of a call.
#[derive(Clone, Copy, Debug)]
pub(super) struct Options {
    skew: bool,
    reps: usize,
    limit: usize,
    threshold: Option<f64>,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            skew: false,
            reps: 1,
            limit: 100,
            threshold: None,
        }
    }
}

/// Remez's answer on `t ∈ [−1, 1]`.
#[derive(Clone, Debug)]
pub(super) struct Approximation {
    /// `p`'s Chebyshev coefficients.
    p: Vec<f64>,
    /// `q`'s, `q[0] = 1`; `[1]` for a polynomial.
    q: Vec<f64>,
    /// The reference of the last exchange, ascending.
    #[cfg_attr(not(test), expect(dead_code, reason = "read by the tests"))]
    reference: Vec<f64>,
    /// `max |e|` over it.
    emax: f64,
    /// `max |f|` over it.
    scale: f64,
    #[cfg_attr(not(test), expect(dead_code, reason = "read by the tests"))]
    iterations: usize,
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &'static str) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

fn small(k: usize) -> f64 {
    f64::from(u32::try_from(k).expect("a count below MAX_ORDER"))
}

/// `Σ c_k T_k(t)` by Clenshaw's recurrence; `c` not empty.
pub(super) fn chebyshev(c: &[f64], t: f64) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = 2.0 * t * b1 - b2 + ck;
        b2 = b1;
        b1 = b0;
    }
    t * b1 - b2 + c[0]
}

/// `T_0(t) … T_{k−1}(t)`.
fn basis(t: f64, k: usize) -> Vec<f64> {
    let mut out = vec![1.0, t];
    out.truncate(k);
    while out.len() < k {
        let j = out.len();
        out.push(2.0 * t * out[j - 1] - out[j - 2]);
    }
    out
}

/// `a x = b` by Gaussian elimination with partial pivoting; `None` when
/// singular or not finite.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        let p = a[pivot][col];
        if !p.is_finite() || p.abs().total_cmp(&0.0) != Ordering::Greater {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (top, bottom) = a.split_at_mut(col + 1);
        let row = &top[col];
        let (btop, bbottom) = b.split_at_mut(col + 1);
        for (r, bi) in bottom.iter_mut().zip(bbottom.iter_mut()) {
            let factor = r[col] / p;
            for (x, y) in r[col..].iter_mut().zip(&row[col..]) {
                *x -= factor * y;
            }
            *bi -= factor * btop[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = a[i][i + 1..]
            .iter()
            .zip(&x[i + 1..])
            .map(|(u, v)| u * v)
            .sum();
        x[i] = (b[i] - s) / a[i][i];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

fn alternating(i: usize) -> f64 {
    if i % 2 == 0 { 1.0 } else { -1.0 }
}

/// The levelled system on the reference `ts` with values `fv`, from `E' =
/// e`: the solution (`p`'s coefficients, `q`'s but the first, `E`) and `E`.
fn levelled(ts: &[f64], fv: &[f64], n: usize, m: usize, e: f64) -> Option<(Vec<f64>, f64)> {
    let big_n = n + m + 2;
    let mut e = e;
    let mut sol = None;
    for _ in 0..ITC_MAX {
        let rows = ts
            .iter()
            .zip(fv)
            .enumerate()
            .map(|(i, (&t, &fi))| {
                let s = alternating(i);
                let b = basis(t, n.max(m) + 1);
                let g = s * e - fi;
                let mut row = Vec::with_capacity(big_n);
                row.extend_from_slice(&b[..=n]);
                row.extend(b[1..=m].iter().map(|v| v * g));
                row.push(s);
                row
            })
            .collect();
        let x = solve(rows, fv.to_vec())?;
        let new = x[big_n - 1];
        let settled = m == 0 || (new - e).abs() <= new.abs() * 2f64.powi(-30);
        e = new;
        sol = Some(x);
        if settled {
            break;
        }
    }
    sol.map(|x| (x, e))
}

/// The argument of the largest `g` on `[lo, hi]` by golden section, the
/// ends included.
fn golden(g: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    let ratio = f64::midpoint(5f64.sqrt(), 1.0);
    let (mut a, mut b) = (lo, hi);
    let mut c = b - (b - a) / ratio;
    let mut d = a + (b - a) / ratio;
    for _ in 0..200 {
        if (c - d).abs() <= 2f64.powi(-40) * c.abs().max(1.0) {
            break;
        }
        if g(c).total_cmp(&g(d)) == Ordering::Greater {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / ratio;
        d = a + (b - a) / ratio;
    }
    let mut x = f64::midpoint(a, b);
    for y in [lo, hi] {
        if g(y).total_cmp(&g(x)) == Ordering::Greater {
            x = y;
        }
    }
    x
}

/// The exchange: `N` points where `e` alternates in sign, its value at each
/// (0 for a padded end); `Err` with the largest `|e|` on the grid when fewer
/// than `N` remain, `Err(NaN)` when `e` is not finite there.
fn exchange(
    e: &dyn Fn(f64) -> f64,
    reference: &[f64],
    big_n: usize,
) -> Result<Vec<(f64, f64)>, f64> {
    let mut points = vec![-1.0];
    points.extend_from_slice(reference);
    points.push(1.0);
    points.sort_by(f64::total_cmp);
    points.dedup_by(|u, v| u.total_cmp(v) == Ordering::Equal);
    let mut grid = vec![-1.0];
    for w in points.windows(2) {
        let (u, v) = (w[0], w[1]);
        for k in 1..SUB {
            grid.push(u + (v - u) * small(k) / small(SUB));
        }
        grid.push(v);
    }
    let values: Vec<f64> = grid.iter().map(|&t| e(t)).collect();
    if values.iter().any(|v| !v.is_finite()) {
        return Err(f64::NAN);
    }
    let positive = |v: f64| v.total_cmp(&0.0) != Ordering::Less;
    // One point per run of one sign: its largest |e|, refined.
    let mut alt: Vec<(f64, f64)> = Vec::new();
    let mut k = 0;
    while k < grid.len() {
        let sign = positive(values[k]);
        let (mut j, mut top) = (k, k);
        while j + 1 < grid.len() && positive(values[j + 1]) == sign {
            j += 1;
            if values[j].abs() > values[top].abs() {
                top = j;
            }
        }
        let s = if sign { 1.0 } else { -1.0 };
        let lo = grid[top.saturating_sub(1).max(k)];
        let hi = grid[(top + 1).min(j)];
        let x = golden(&|t| s * e(t), lo, hi);
        let v = e(x);
        // The refined point keeps the run's sign, else the grid's.
        alt.push(if positive(v) == sign {
            (x, v)
        } else {
            (grid[top], values[top])
        });
        k = j + 1;
    }
    if alt.len() < big_n && alt[0].0 > -1.0 {
        alt.insert(0, (-1.0, 0.0));
    }
    if alt.len() < big_n && alt[alt.len() - 1].0 < 1.0 {
        alt.push((1.0, 0.0));
    }
    if alt.len() < big_n {
        return Err(values.iter().map(|v| v.abs()).fold(0.0, f64::max));
    }
    let g = (0..alt.len())
        .max_by(|&i, &j| alt[i].1.abs().total_cmp(&alt[j].1.abs()).then(j.cmp(&i)))
        .expect("N points");
    let first = g.saturating_sub(big_n - 1);
    let last = g.min(alt.len() - big_n);
    let least = |s: usize| {
        alt[s..s + big_n]
            .iter()
            .map(|r| r.1.abs())
            .fold(f64::INFINITY, f64::min)
    };
    let start = (first..=last)
        .max_by(|&i, &j| least(i).total_cmp(&least(j)).then(j.cmp(&i)))
        .expect("a window");
    Ok(alt[start..start + big_n].to_vec())
}

/// `a`, `b` and the `k` Chebyshev nodes on `[−1, 1]`, ascending
/// (`chebyshev_nodes`).
fn chebyshev_nodes(k: usize) -> Vec<f64> {
    let mut nodes = vec![-1.0];
    for i in (1..=k).rev() {
        nodes.push((small(2 * i - 1) * std::f64::consts::PI / small(2 * k)).cos());
    }
    nodes.push(1.0);
    nodes
}

/// giac's `skew_nodes` on `[−1, 1]`, the ends kept.
fn skew(nodes: &mut [f64], draw: &mut dyn FnMut() -> u64) {
    let u = f64::from(u32::try_from(draw() >> 32).expect("32 bits")) * 2f64.powi(-32);
    let s = 0.75 + 0.5 * u;
    let symmetric = draw() & 1 == 1;
    let last = nodes.len() - 1;
    for t in &mut nodes[1..last] {
        *t = if symmetric {
            t.abs().powf(s).copysign(*t)
        } else {
            f64::midpoint(*t, 1.0).powf(s) * 2.0 - 1.0
        };
    }
}

/// Remez's algorithm for `f` on `[−1, 1]`; `None` when no repetition
/// converged within the limit.
pub(super) fn remez(
    f: &dyn Fn(f64) -> f64,
    n: usize,
    m: usize,
    o: &Options,
    draw: &mut dyn FnMut() -> u64,
) -> Option<Approximation> {
    let k = n + m;
    let big_n = k + 2;
    let mut answer: Option<Approximation> = None;
    for _ in 0..o.reps {
        let mut nodes = chebyshev_nodes(k);
        let mut level = 0.0;
        let mut best: Option<Approximation> = None;
        let (mut converged, mut iterations) = (false, 0);
        while iterations < o.limit {
            iterations += 1;
            let first = o.skew && iterations == 1;
            let mut chosen: Option<(Vec<f64>, Vec<f64>, f64)> = None;
            for _ in 0..=(if first { TRIES } else { 0 }) {
                if first {
                    skew(&mut nodes, draw);
                }
                let fv: Vec<f64> = nodes.iter().map(|&t| f(t)).collect();
                if fv.iter().any(|v| !v.is_finite()) {
                    return None;
                }
                let Some((sol, e)) = levelled(&nodes, &fv, n, m, level) else {
                    continue;
                };
                if chosen.as_ref().is_none_or(|c| e.abs() < c.2.abs()) {
                    chosen = Some((nodes.clone(), sol, e));
                }
            }
            let Some((reference, sol, e)) = chosen else {
                break;
            };
            level = e;
            let p = sol[..=n].to_vec();
            let mut q = vec![1.0];
            q.extend_from_slice(&sol[n + 1..big_n - 1]);
            let err = |t: f64| f(t) - chebyshev(&p, t) / chebyshev(&q, t);
            let new = match exchange(&err, &reference, big_n) {
                Ok(new) => new,
                Err(emax) => {
                    // No alternation left: e vanishes to rounding (f is
                    // p/q), else a failure.
                    let scale = reference.iter().map(|&t| f(t).abs()).fold(0.0, f64::max);
                    if emax <= 64.0 * f64::EPSILON * scale {
                        best = Some(Approximation {
                            p,
                            q,
                            reference,
                            emax,
                            scale,
                            iterations,
                        });
                        converged = true;
                    }
                    break;
                }
            };
            nodes = new.iter().map(|r| r.0).collect();
            let emin = new.iter().map(|r| r.1.abs()).fold(f64::INFINITY, f64::min);
            let emax = new.iter().map(|r| r.1.abs()).fold(0.0, f64::max);
            let scale = nodes.iter().map(|&t| f(t).abs()).fold(0.0, f64::max);
            if best.as_ref().is_none_or(|b| emax < b.emax) {
                best = Some(Approximation {
                    p,
                    q,
                    reference: nodes.clone(),
                    emax,
                    scale,
                    iterations,
                });
            }
            if converged {
                break;
            }
            converged = match o.threshold {
                Some(t) => emax <= t * emin,
                None => emax - emin <= (2f64.powi(-20) * emin).max(64.0 * f64::EPSILON * scale),
            };
        }
        if let Some(b) = best.filter(|_| converged)
            && answer.as_ref().is_none_or(|a| b.emax < a.emax)
        {
            answer = Some(b);
        }
    }
    answer
}

/// `x == a..b`: the variable and the ends.
fn domain(v: &Expr) -> Option<(&Expr, &Expr, &Expr)> {
    let head = |e: &Expr| e.head().and_then(Expr::as_symbol);
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, r]) = (head(v), v.args()) else {
        return None;
    };
    let (Some(Symbol::RANGE), [a, b]) = (head(r), r.args()) else {
        return None;
    };
    x.as_symbol()
        .is_some_and(|s| !s.is_builtin())
        .then_some((x, a, b))
}

/// `n ≥ 1`, or `{n, m}` with `n ≥ 0`, `m ≥ 1`.
fn order(e: &Expr) -> Option<(usize, usize)> {
    let size = |e: &Expr| e.as_i64().and_then(|k| usize::try_from(k).ok());
    if let Some(n) = size(e) {
        return (n >= 1).then_some((n, 0));
    }
    let (Some(Symbol::LIST), [n, m]) = (e.head().and_then(Expr::as_symbol), e.args()) else {
        return None;
    };
    let (n, m) = (size(n)?, size(m)?);
    (m >= 1).then_some((n, m))
}

/// The options; `None` for another (giac's errors).
fn options(given: &[Expr]) -> Option<Options> {
    let mut o = Options::default();
    for option in given {
        if option.as_symbol() == Some(Symbol::RANDOM) {
            o.skew = true;
            continue;
        }
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [lh, rh]) =
            (option.head().and_then(Expr::as_symbol), option.args())
        else {
            return None;
        };
        let positive = || {
            rh.as_i64()
                .filter(|&k| k > 0)
                .and_then(|k| usize::try_from(k).ok())
        };
        match lh.as_symbol() {
            Some(Symbol::RANDOM) => {
                o.reps = positive()?;
                o.skew = true;
            }
            Some(Symbol::LIMIT) => o.limit = positive()?,
            Some(Symbol::THRESHOLD) => {
                let t = match rh.view() {
                    View::Ball(b) if b.is_finite() => b.mid().to_f64(),
                    _ => to_f64(&rh.to_rational()?),
                };
                o.threshold = Some(t).filter(|t| t.is_finite() && *t > 1.0);
                o.threshold?;
            }
            _ => return None,
        }
    }
    Some(o)
}

fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// The finite ball of a real number.
fn ball(ev: &Eval<'_>, e: &Expr, prec: u64) -> Option<Ball> {
    if let Some(q) = e.to_rational() {
        return Some(Ball::from_rational(&q, prec));
    }
    let value = match e.view() {
        View::Ball(b) => return b.is_finite().then(|| b.clone()),
        _ => approximation(ev.context(), e, prec).ok()?,
    };
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

/// A real value of `tree` at the ball `z`.
fn value(tree: &Node, z: &Ball, prec: u64) -> Option<Ball> {
    tree.eval(&ComplexBall::from_real(z.clone()), Mode::Real, prec)
        .filter(|v| v.is_real() && v.re().is_finite())
        .map(|v| v.re().clone())
}

/// `Σ c_j T_j((x − c)/h)` as monomials in `x`, the constant first.
fn monomials(coefficients: &[f64], c: &Rational, h: &Rational) -> Vec<Rational> {
    let zero = Rational::from(0);
    let u = [-(c / h), h.recip()];
    let mut out = vec![zero.clone(); coefficients.len()];
    let (mut previous, mut current): (Vec<Rational>, Vec<Rational>) =
        (vec![], vec![Rational::from(1)]);
    for (j, &cj) in coefficients.iter().enumerate() {
        let cj = Float::from_f64(cj).map_or(zero.clone(), |v| v.to_rational());
        for (o, t) in out.iter_mut().zip(&current) {
            *o += &(&cj * t);
        }
        if j + 1 < coefficients.len() {
            // T_{j+1} = u T_j for j = 0, else 2u T_j − T_{j−1}.
            let mut next = vec![zero.clone(); current.len() + 1];
            let two = if j == 0 {
                Rational::from(1)
            } else {
                Rational::from(2)
            };
            for (i, t) in current.iter().enumerate() {
                let t2 = &two * t;
                next[i] += &(&t2 * &u[0]);
                next[i + 1] += &(&t2 * &u[1]);
            }
            for (i, t) in previous.iter().enumerate() {
                next[i] -= t;
            }
            previous = std::mem::replace(&mut current, next);
        }
    }
    out
}

/// `Σ |c_k| M^k`.
fn size_on(c: &[Rational], big_m: &Rational) -> Rational {
    let mut s = Rational::from(0);
    let mut power = Rational::from(1);
    for ck in c {
        s += &(&ck.abs() * &power);
        power = &power * big_m;
    }
    s
}

/// Drops the terms below `tau` over `|x| ≤ M` and rounds the others to
/// exact balls at the precision that moves the sum by at most `tau`.
fn round(c: &[Rational], big_m: &Rational, tau: &Rational) -> Vec<Ball> {
    let mut power = Rational::from(1);
    let mut kept = Vec::with_capacity(c.len());
    for ck in c {
        kept.push(if &ck.abs() * &power < *tau {
            Rational::from(0)
        } else {
            ck.clone()
        });
        power = &power * big_m;
    }
    let s = size_on(&kept, big_m);
    // |c − round(c)| ≤ 2^−prec |c|: the sum moves by at most 2^−prec s.
    let prec = if tau.is_zero() || s.is_zero() {
        53
    } else {
        // log2(s/tau) < bits(num) − bits(den) + 1.
        let ratio = &s / tau;
        let bits = (ratio.numer().bit_length() + 1).saturating_sub(ratio.denom().bit_length());
        bits.clamp(2, 1 << 16)
    };
    kept.iter()
        .map(|ck| Ball::exact(Float::from_rational(ck, prec, Round::Nearest).0))
        .collect()
}

/// The coefficients of `p(m + s)` in `s`.
fn shift(p: &[Ball], m: &Ball, prec: u64) -> Vec<Ball> {
    let mut c = p.to_vec();
    let d = c.len();
    for i in 0..d {
        for j in (i..d - 1).rev() {
            let t = m.mul(&c[j + 1], prec);
            c[j] = c[j].add(&t, prec);
        }
    }
    c
}

/// `Σ c_k s^k` and `Σ k c_k s^(k−1)` at the ball `s` by Horner.
fn horner(c: &[Ball], s: &Ball, prec: u64) -> (Ball, Ball) {
    let mut v = Ball::ZERO;
    let mut dv = Ball::ZERO;
    for ck in c.iter().rev() {
        dv = dv.mul(s, prec).add(&v, prec);
        v = v.mul(s, prec).add(ck, prec);
    }
    (v, dv)
}

/// An upper bound of `|b|`, `+∞` for a ball not finite.
fn upper(b: &Ball) -> Mag {
    if !b.is_finite() {
        return Mag::INF;
    }
    Mag::from_float(b.mid()).add(b.rad())
}

/// A lower bound of `|b|`.
fn lower(b: &Ball) -> Mag {
    if !b.is_finite() {
        return Mag::ZERO;
    }
    Mag::from_float_lower(b.mid()).sub_lower(b.rad())
}

/// What the bound of a piece found.
enum Piece {
    /// The lower bound at the centre (when the centre is inside), the upper.
    Bound(Option<Mag>, Mag),
    /// `q` changes sign or vanishes.
    Pole,
}

/// The bound of the error on the piece `[l, r]`.
struct Certifier<'a> {
    f: &'a Node,
    df: Option<&'a Node>,
    p: &'a [Ball],
    q: &'a [Ball],
    prec: u64,
    /// The sign of `q` at the first centre.
    q_sign: Option<Ordering>,
}

impl Certifier<'_> {
    fn piece(&mut self, l: &Ball, r: &Ball) -> Piece {
        let prec = self.prec;
        let mid = l.mid().add(r.mid()).mul_2exp(&Integer::from(-1));
        let half = Mag::from_float(&r.mid().sub(l.mid()).mul_2exp(&Integer::from(-1)));
        let radius = half.add(&std::cmp::max(l.rad().clone(), r.rad().clone()));
        let centre = Ball::exact(mid.clone());
        let s = Ball::new(Float::ZERO, radius.clone());
        let interval = Ball::new(mid, radius);
        let (pc, qc) = (shift(self.p, &centre, prec), shift(self.q, &centre, prec));
        let (p0, q0) = (&pc[0], &qc[0]);
        let (pi, dpi) = horner(&pc, &s, prec);
        let (qi, dqi) = horner(&qc, &s, prec);
        // q at the centre: its sign must stay.
        match q0.sign() {
            Some(Ordering::Equal) | None => return Piece::Pole,
            Some(sign) => {
                if *self.q_sign.get_or_insert(sign) != sign {
                    return Piece::Pole;
                }
            }
        }
        let inside =
            l.cmp(&centre) == Some(Ordering::Less) && centre.cmp(r) == Some(Ordering::Less);
        let at_centre = value(self.f, &centre, prec).map(|fc| fc.sub(&p0.div(q0, prec), prec));
        let low = at_centre.as_ref().filter(|_| inside).map(lower);
        if qi.sign().is_none_or(|o| o == Ordering::Equal) {
            return Piece::Bound(low, Mag::INF);
        }
        let natural = value(self.f, &interval, prec)
            .map_or(Mag::INF, |fi| upper(&fi.sub(&pi.div(&qi, prec), prec)));
        let mean = match (self.df, &at_centre) {
            (Some(df), Some(ec)) => value(df, &interval, prec).map_or(Mag::INF, |dfi| {
                let dr = dpi
                    .mul(&qi, prec)
                    .sub(&pi.mul(&dqi, prec), prec)
                    .div(&qi.mul(&qi, prec), prec);
                upper(&ec.add(&dfi.sub(&dr, prec).mul(&s, prec), prec))
            }),
            _ => Mag::INF,
        };
        Piece::Bound(low, std::cmp::min(natural, mean))
    }
}

/// The certified `[L, U]` of `max |f − p/q|` on `[lo, hi]`: `Err(true)` for
/// a pole, `Err(false)` for no finite bound within `cap` pieces.
fn certify(c: &mut Certifier<'_>, lo: &Ball, hi: &Ball, cap: usize) -> Result<(Mag, Mag), bool> {
    let mut live: Vec<(Ball, Ball, Mag)> = Vec::new();
    let mut best = Mag::ZERO;
    let add = |c: &mut Certifier<'_>,
               l: Ball,
               r: Ball,
               live: &mut Vec<(Ball, Ball, Mag)>,
               best: &mut Mag| match c.piece(&l, &r) {
        Piece::Pole => Err(true),
        Piece::Bound(low, up) => {
            if let Some(low) = low {
                *best = std::cmp::max(best.clone(), low);
            }
            live.push((l, r, up));
            Ok(())
        }
    };
    add(c, lo.clone(), hi.clone(), &mut live, &mut best)?;
    let mut pieces = 1;
    loop {
        live.retain(|p| p.2 >= best);
        let top = live.iter().map(|p| p.2.clone()).max().unwrap_or(Mag::ZERO);
        let target = best.add(&best.mul_2exp(&Integer::from(-20)));
        if top <= target || pieces >= cap {
            return if top.is_finite() {
                Ok((best, top))
            } else {
                Err(false)
            };
        }
        let (split, keep): (Vec<_>, Vec<_>) = std::mem::take(&mut live)
            .into_iter()
            .partition(|p| p.2 > target);
        live = keep;
        for (l, r, up) in split {
            if pieces >= cap {
                live.push((l, r, up));
                continue;
            }
            let m = Ball::exact(l.mid().add(r.mid()).mul_2exp(&Integer::from(-1)));
            add(c, l, m.clone(), &mut live, &mut best)?;
            add(c, m, r, &mut live, &mut best)?;
            pieces += 2;
        }
    }
}

/// The answer's coefficients from Remez's (the constant first): monomials
/// in `x` exactly, `q`'s constant term 1, small terms dropped and the others
/// rounded (the module's rules); `None` for `q` = 0.
fn coefficients(
    r: &Approximation,
    m: usize,
    (cq, hq): (&Rational, &Rational),
    big_m: &Rational,
) -> Option<(Vec<Ball>, Vec<Ball>)> {
    let (mut p, mut q) = (monomials(&r.p, cq, hq), monomials(&r.q, cq, hq));
    if m > 0 {
        // The first term of q not below 2^−40 of its size on [a, b] made 1,
        // those before it (rounding, as for f = 1/x) dropped.
        let floor = &size_on(&q, big_m) / &Rational::from(Integer::from(1_i64 << 40));
        let mut power = Rational::from(1);
        let mut lead = None;
        for v in &mut q {
            if &v.abs() * &power >= floor {
                lead = Some(v.clone());
                break;
            }
            *v = Rational::from(0);
            power = &power * big_m;
        }
        let lead = lead?;
        for v in p.iter_mut().chain(q.iter_mut()) {
            *v = &*v / &lead;
        }
    }
    let size = r.emax.max(64.0 * f64::EPSILON * r.scale);
    let tau = Float::from_f64(size * 2f64.powi(-40)).map_or(Rational::from(0), |v| v.to_rational());
    Some((round(&p, big_m, &tau), round(&q, big_m, &tau)))
}

/// `Σ c_k x^k`, the zero terms left out.
fn poly(x: &Expr, c: Vec<Ball>) -> Expr {
    let terms: Vec<Expr> = c
        .into_iter()
        .enumerate()
        .filter(|(_, v)| !v.mid().is_zero())
        .map(|(k, v)| Expr::from(v) * x.pow(i64::try_from(k).expect("a degree below MAX_ORDER")))
        .collect();
    Expr::call(Symbol::PLUS, terms)
}

/// giac's `minimax`: `{E, p}` or `{E, p/q}`; `None` where it stays.
pub(in crate::native) fn minimax_approximation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, range, degrees, given @ ..] = e.args() else {
        return Ok(None);
    };
    let (Some((x, a, b)), Some((n, m)), Some(o)) = (domain(range), order(degrees), options(given))
    else {
        return Ok(None);
    };
    // The work bounded before any, as a product (C0019).
    let big_n = n + m + 2;
    let terms = ev.context().limits().terms;
    let work = (big_n * big_n)
        .checked_mul(o.limit)
        .and_then(|w| w.checked_mul(o.reps));
    if n + m > MAX_ORDER || work.is_none_or(|w| w > terms) {
        warn(ev, e, "minimax: past the limits");
        return Ok(None);
    }
    let cap = (terms / (big_n * big_n)).min(MAX_PIECES);
    let wp = ev.context().working_precision() + GUARD;
    let (Some(lo), Some(hi)) = (ball(ev, a, wp), ball(ev, b, wp)) else {
        return Ok(None);
    };
    if lo.cmp(&hi) != Some(Ordering::Less) {
        return Ok(None);
    }
    let (af, bf) = (lo.mid().to_f64(), hi.mid().to_f64());
    let (Some(aq), Some(bq)) = (Float::from_f64(af), Float::from_f64(bf)) else {
        return Ok(None);
    };
    let (aq, bq) = (aq.to_rational(), bq.to_rational());
    let two = Rational::from(2);
    let (cq, hq) = (&(&aq + &bq) / &two, &(&bq - &aq) / &two);
    let (c, h) = (to_f64(&cq), to_f64(&hq));
    if !(c.is_finite() && h.is_finite() && hq.is_positive()) {
        return Ok(None);
    }
    let Some(tree) = Node::compile(ev.context(), f, x, wp) else {
        return Ok(None);
    };
    let fx = |t: f64| {
        Float::from_f64(c + h * t)
            .and_then(|z| value(&tree, &Ball::exact(z), 64))
            .map_or(f64::NAN, |v| v.mid().to_f64())
    };
    let Some(r) = remez(&fx, n, m, &o, &mut || {
        ev.draw(cas_num::random::Generator::next_u64)
    }) else {
        warn(ev, e, "minimax: no convergence within the iteration limit");
        return Ok(None);
    };
    // Monomials, q's constant term 1.
    let big_m = std::cmp::max(aq.abs(), bq.abs());
    let Some((pb, qb)) = coefficients(&r, m, (&cq, &hq), &big_m) else {
        return Ok(None);
    };
    let df = ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?;
    let dtree = Node::compile(ev.context(), &df, x, wp);
    let mut certifier = Certifier {
        f: &tree,
        df: dtree.as_ref(),
        p: &pb,
        q: &qb,
        prec: wp,
        q_sign: None,
    };
    let (low, up) = match certify(&mut certifier, &lo, &hi, cap) {
        Ok(bound) => bound,
        Err(pole) => {
            warn(
                ev,
                e,
                if pole {
                    "minimax: the approximation has a pole on the interval"
                } else {
                    "minimax: no bound of the error found"
                },
            );
            return Ok(None);
        }
    };
    let (Some(l), Some(u)) = (low.to_float(), up.to_float()) else {
        return Ok(None);
    };
    let error = Ball::new(
        l.add(&u).mul_2exp(&Integer::from(-1)),
        Mag::from_float(&u.sub(&l).mul_2exp(&Integer::from(-1))),
    );
    let answer = if m == 0 {
        poly(x, pb)
    } else {
        poly(x, pb) / poly(x, qb)
    };
    Ok(Some(ev.eval(&Expr::call(
        Symbol::LIST,
        [Expr::from(error), answer],
    ))?))
}
