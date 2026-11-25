//! Empirical mode decomposition (CAS-635), giac's `emd`
//! (signalprocessing.cc:5131-5281, `emd` and `_emd`).
//!
//! `EmpiricalModeDecomposition[data, options…]`: `data` a list of real
//! numbers (rationals or decimals; anything else leaves the call with
//! giac's warning "Data must be numeric"), the options giac's
//! `threshold = v` (0.1, `v ≥ 0`), `limit = m` (the most IMFs, `m ≥ 1`;
//! unbounded), `extrema = k` (the fewest extrema of an IMF, `k ≥ 2`; 2) and
//! `residue = true|false` (true); another option, or a value giac refuses,
//! leaves the call. The answer is `Sequence[imfs, residue]`, or `imfs` with
//! `residue = false`.
//!
//! The algorithm, giac's on doubles: `x` the samples; each round sifts `h =
//! x`: the extrema of `h` (index `i` in `1 … n−3` a minimum when `h_(i+1) −
//! h_i > 0 > h_i − h_(i−1)`, or on a plateau `h_(i+1) = h_i` with
//! `h_i − h_(i−1) < 0 < h_(i+2) − h_(i+1)`; maxima the same reversed),
//! both ends knots of both envelopes; `ne` the interior extrema; `ne <
//! extrema` ends the decomposition with `h` the residue. Else the upper and
//! lower envelopes are GSL's splines through the knots (natural cubic,
//! `gsl_interp_cspline`, interpolation/cspline.c `cspline_init` and
//! `cspline_eval`, the system solved by linalg/tridiag.c `solve_tridiag`;
//! linear through two knots, interpolation/linear.c `linear_eval`), `m` their
//! mean `0.5 lo + 0.5 up`, `var = Σ m_i²/(h_i² + 10⁻⁷) / n` and `h ← h − m`,
//! until `var ≤ threshold` and the zero crossings (pairs `0 … n−4` with
//! `h_(i−1) = 0` or a sign change) differ from `ne` by at most 1. `h` is then
//! an IMF and `x ← x − h`; `limit` IMFs end it with `x` the residue.
//!
//! Every decision (extrema, crossings, the stopping test) is taken on
//! doubles computed as giac's are, operation for operation, so the IMFs are
//! giac's. The values are balls (C0064): each IMF's last sift is recomputed
//! in balls at the working precision from the exact doubles before it, the
//! residue is the data less the IMFs' balls (or the last sift's ball).
//!
//! Bounds: the sifts times `n` at most `Limits::terms` (a threshold of 0
//! may never stop: giac loops), beyond which the call stays with a
//! warning; empty data stays.

#![allow(
    clippy::many_single_char_names,
    reason = "GSL's names: x, y, h, c, d, g, m"
)]

use cas_num::ball::{Ball, Float, Round};

use super::{is_list, list, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::statistics::random_variables::equation;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The arithmetic of a sift: doubles for giac's decisions, balls for the
/// values.
trait Arith {
    type T: Clone;
    /// An exact constant (a knot's abscissa, 0.5, 2, 3).
    fn of(&self, x: f64) -> Self::T;
    fn add(&self, a: &Self::T, b: &Self::T) -> Self::T;
    fn sub(&self, a: &Self::T, b: &Self::T) -> Self::T;
    fn mul(&self, a: &Self::T, b: &Self::T) -> Self::T;
    fn div(&self, a: &Self::T, b: &Self::T) -> Self::T;
}

struct Doubles;

impl Arith for Doubles {
    type T = f64;
    fn of(&self, x: f64) -> f64 {
        x
    }
    fn add(&self, a: &f64, b: &f64) -> f64 {
        a + b
    }
    fn sub(&self, a: &f64, b: &f64) -> f64 {
        a - b
    }
    fn mul(&self, a: &f64, b: &f64) -> f64 {
        a * b
    }
    fn div(&self, a: &f64, b: &f64) -> f64 {
        a / b
    }
}

/// Balls at a precision.
struct Balls(u64);

impl Arith for Balls {
    type T = Ball;
    fn of(&self, x: f64) -> Ball {
        Ball::exact(Float::from_f64(x).expect("a finite constant"))
    }
    fn add(&self, a: &Ball, b: &Ball) -> Ball {
        a.add(b, self.0)
    }
    fn sub(&self, a: &Ball, b: &Ball) -> Ball {
        a.sub(b, self.0)
    }
    fn mul(&self, a: &Ball, b: &Ball) -> Ball {
        a.mul(b, self.0)
    }
    fn div(&self, a: &Ball, b: &Ball) -> Ball {
        a.div(b, self.0)
    }
}

#[allow(
    clippy::cast_precision_loss,
    reason = "indices of at most Limits::terms samples, exact in a double"
)]
fn abscissa(i: usize) -> f64 {
    i as f64
}

/// GSL's `solve_tridiag`: the symmetric system of diagonal `d`, off
/// diagonal `o` (one shorter) and right side `b`, at least 2 rows.
fn tridiagonal<A: Arith>(a: &A, d: &[A::T], o: &[A::T], b: &[A::T]) -> Vec<A::T> {
    let n = d.len();
    let mut alpha = Vec::with_capacity(n);
    let mut gamma = Vec::with_capacity(n);
    alpha.push(d[0].clone());
    gamma.push(a.div(&o[0], &alpha[0]));
    for i in 1..n - 1 {
        alpha.push(a.sub(&d[i], &a.mul(&o[i - 1], &gamma[i - 1])));
        gamma.push(a.div(&o[i], &alpha[i]));
    }
    alpha.push(a.sub(&d[n - 1], &a.mul(&o[n - 2], &gamma[n - 2])));
    let mut z = Vec::with_capacity(n);
    z.push(b[0].clone());
    for i in 1..n {
        let next = a.sub(&b[i], &a.mul(&gamma[i - 1], &z[i - 1]));
        z.push(next);
    }
    let c: Vec<A::T> = z.iter().zip(&alpha).map(|(z, al)| a.div(z, al)).collect();
    let mut x = c.clone();
    for i in (0..n - 1).rev() {
        x[i] = a.sub(&c[i], &a.mul(&gamma[i], &x[i + 1]));
    }
    x
}

/// GSL's spline through the knots `(xs_j, ys_j)`, `xs` increasing indices
/// from 0 to `n − 1`, at `0 … n − 1`: linear through two knots, else the
/// natural cubic.
fn spline<A: Arith>(a: &A, xs: &[usize], ys: &[A::T], n: usize) -> Vec<A::T> {
    let size = xs.len();
    let x = |i: usize| a.of(abscissa(xs[i]));
    if size == 2 {
        let dx = a.sub(&x(1), &x(0));
        let dy = a.sub(&ys[1], &ys[0]);
        return (0..n)
            .map(|t| {
                let delx = a.sub(&a.of(abscissa(t)), &x(0));
                a.add(&ys[0], &a.mul(&a.div(&delx, &dx), &dy))
            })
            .collect();
    }
    let (one, two, three) = (a.of(1.0), a.of(2.0), a.of(3.0));
    let max_index = size - 1;
    let sys = max_index - 1;
    let mut off = Vec::with_capacity(sys);
    let mut diag = Vec::with_capacity(sys);
    let mut g = Vec::with_capacity(sys);
    for i in 0..sys {
        let h_i = a.sub(&x(i + 1), &x(i));
        let h_ip1 = a.sub(&x(i + 2), &x(i + 1));
        let ydiff_i = a.sub(&ys[i + 1], &ys[i]);
        let ydiff_ip1 = a.sub(&ys[i + 2], &ys[i + 1]);
        // The knots are distinct: no zero step.
        let g_i = a.div(&one, &h_i);
        let g_ip1 = a.div(&one, &h_ip1);
        diag.push(a.mul(&two, &a.add(&h_ip1, &h_i)));
        off.push(h_ip1);
        g.push(a.mul(
            &three,
            &a.sub(&a.mul(&ydiff_ip1, &g_ip1), &a.mul(&ydiff_i, &g_i)),
        ));
    }
    let mut c = vec![a.of(0.0); size];
    if sys == 1 {
        c[1] = a.div(&g[0], &diag[0]);
    } else {
        let solved = tridiagonal(a, &diag, &off[..sys - 1], &g);
        c[1..max_index].clone_from_slice(&solved);
    }
    let mut k = 0;
    (0..n)
        .map(|t| {
            while k < size - 2 && xs[k + 1] <= t {
                k += 1;
            }
            let dx = a.sub(&x(k + 1), &x(k));
            let dy = a.sub(&ys[k + 1], &ys[k]);
            let delx = a.sub(&a.of(abscissa(t)), &x(k));
            let (c_i, c_ip1) = (&c[k], &c[k + 1]);
            // b = dy/dx − dx (c_(i+1) + 2 c_i)/3, d = (c_(i+1) − c_i)/(3 dx).
            let b = a.sub(
                &a.div(&dy, &dx),
                &a.div(&a.mul(&dx, &a.add(c_ip1, &a.mul(&two, c_i))), &three),
            );
            let d = a.div(&a.sub(c_ip1, c_i), &a.mul(&three, &dx));
            let inner = a.add(c_i, &a.mul(&delx, &d));
            a.add(&ys[k], &a.mul(&delx, &a.add(&b, &a.mul(&delx, &inner))))
        })
        .collect()
}

/// The knots of the two envelopes (maxima, minima; both ends in each) and
/// the zero crossings, giac's scan.
fn extrema(h: &[f64]) -> (Vec<usize>, Vec<usize>, usize) {
    let n = h.len();
    let (mut maxima, mut minima, mut crossings) = (vec![0], vec![0], 0);
    for i in 1..n.saturating_sub(2) {
        let (h1, h2, h3, h4) = (h[i - 1], h[i], h[i + 1], h[i + 2]);
        if h1 == 0.0 || h1 * h2 < 0.0 {
            crossings += 1;
        }
        let (d, dl, dr) = (h3 - h2, h2 - h1, h4 - h3);
        if (d == 0.0 && dl < 0.0 && dr > 0.0) || (dl < 0.0 && d > 0.0) {
            minima.push(i);
        } else if (d == 0.0 && dl > 0.0 && dr < 0.0) || (dl > 0.0 && d < 0.0) {
            maxima.push(i);
        }
    }
    maxima.push(n - 1);
    minima.push(n - 1);
    (maxima, minima, crossings)
}

/// The mean of the envelopes of `h` through the knots.
fn mean<A: Arith>(a: &A, h: &[A::T], maxima: &[usize], minima: &[usize]) -> Vec<A::T> {
    let n = h.len();
    let values = |knots: &[usize]| knots.iter().map(|&i| h[i].clone()).collect::<Vec<_>>();
    let up = spline(a, maxima, &values(maxima), n);
    let lo = spline(a, minima, &values(minima), n);
    let half = a.of(0.5);
    lo.iter()
        .zip(&up)
        .map(|(l, u)| a.add(&a.mul(&half, l), &a.mul(&half, u)))
        .collect()
}

/// A sift `h − m` in balls from the exact doubles `h`; `None` when one is
/// not finite.
fn ball_sift(h: &[f64], maxima: &[usize], minima: &[usize], prec: u64) -> Option<Vec<Ball>> {
    let a = Balls(prec);
    let h = h
        .iter()
        .map(|&v| Float::from_f64(v).map(Ball::exact))
        .collect::<Option<Vec<_>>>()?;
    let m = mean(&a, &h, maxima, minima);
    Some(h.iter().zip(&m).map(|(h, m)| a.sub(h, m)).collect())
}

/// A real number as giac's double: a rational to nearest, a ball's
/// midpoint; `None` otherwise or when not finite.
fn double(v: &Expr) -> Option<f64> {
    let x = if let View::Ball(b) = v.view() {
        b.mid().to_f64()
    } else {
        let q = v.to_rational()?;
        Float::from_rational(&q, 53, Round::Nearest).0.to_f64()
    };
    x.is_finite().then_some(x)
}

/// A sample as a ball at `prec` bits.
fn sample(v: &Expr, prec: u64) -> Option<Ball> {
    match v.view() {
        View::Ball(b) => Some(b.clone()),
        _ => v.to_rational().map(|q| Ball::from_rational(&q, prec)),
    }
}

/// giac's options, or `None` for one giac refuses.
struct Options {
    threshold: f64,
    limit: usize,
    extrema: usize,
    residue: bool,
}

fn options(given: &[Expr]) -> Option<Options> {
    let mut o = Options {
        threshold: 0.1,
        limit: usize::MAX,
        extrema: 2,
        residue: true,
    };
    for option in given {
        let (k, v) = equation(option)?;
        match k.as_symbol()? {
            Symbol::THRESHOLD => o.threshold = double(v).filter(|t| *t >= 0.0)?,
            Symbol::LIMIT => o.limit = usize::try_from(v.as_i64()?).ok().filter(|m| *m >= 1)?,
            Symbol::EXTREMA => {
                o.extrema = usize::try_from(v.as_i64()?).ok().filter(|k| *k >= 2)?;
            }
            Symbol::RESIDUE => {
                o.residue = match v.as_symbol()? {
                    Symbol::TRUE => true,
                    Symbol::FALSE => false,
                    _ => return None,
                };
            }
            _ => return None,
        }
    }
    Some(o)
}

/// `EmpiricalModeDecomposition[data, options…]`, giac's `emd`.
pub(in crate::native) fn empirical_mode_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [data, given @ ..] = e.args() else {
        return Ok(None);
    };
    if !is_list(data) || data.args().is_empty() {
        return Ok(None);
    }
    let Some(o) = options(given) else {
        return Ok(None);
    };
    let n = data.args().len();
    let limit = ev.context().limits().terms;
    if n > limit {
        return Ok(None);
    }
    let prec = ev.context().working_precision();
    let (Some(mut x), Some(mut x_ball)) = (
        data.args().iter().map(double).collect::<Option<Vec<_>>>(),
        data.args()
            .iter()
            .map(|v| sample(v, prec))
            .collect::<Option<Vec<_>>>(),
    ) else {
        warn(ev, e, "Data must be numeric");
        return Ok(None);
    };
    let a = Doubles;
    let budget = limit / n;
    let mut sifts = 0;
    let mut imfs: Vec<Vec<Ball>> = Vec::new();
    let residue = loop {
        let mut h = x.clone();
        // The doubles before the last sift and its knots.
        let mut last: Option<(Vec<f64>, Vec<usize>, Vec<usize>)> = None;
        let mut is_residue = false;
        loop {
            sifts += 1;
            if sifts > budget {
                warn(
                    ev,
                    e,
                    format!("more than {budget} sifts of {n} samples (Limits::terms)"),
                );
                return Ok(None);
            }
            let (maxima, minima, crossings) = extrema(&h);
            let ne = maxima.len() + minima.len() - 4;
            if ne < o.extrema {
                is_residue = true;
                break;
            }
            let not_imf = crossings.abs_diff(ne) > 1;
            let m = mean(&a, &h, &maxima, &minima);
            let mut var = 0.0;
            for (m, h) in m.iter().zip(&h) {
                var += (m * m) / (h * h + 1e-7);
            }
            var /= abscissa(n);
            let next: Vec<f64> = h.iter().zip(&m).map(|(h, m)| h - m).collect();
            last = Some((h, maxima, minima));
            h = next;
            let again = var > o.threshold || not_imf;
            if !again {
                break;
            }
        }
        let sifted = match &last {
            Some((before, maxima, minima)) => {
                let Some(found) = ball_sift(before, maxima, minima, prec) else {
                    warn(ev, e, "a sift that is not finite");
                    return Ok(None);
                };
                Some(found)
            }
            None => None,
        };
        if is_residue {
            break sifted.unwrap_or(x_ball);
        }
        // The loop ran at least once: its first test always holds.
        let Some(imf) = sifted else {
            return Ok(None);
        };
        x = x.iter().zip(&h).map(|(x, h)| x - h).collect();
        x_ball = x_ball
            .iter()
            .zip(&imf)
            .map(|(x, h)| x.sub(h, prec))
            .collect();
        imfs.push(imf);
        if imfs.len() == o.limit {
            break x_ball;
        }
    };
    let as_list = |v: Vec<Ball>| list(v.into_iter().map(Expr::from).collect());
    let imfs = list(imfs.into_iter().map(as_list).collect());
    let result = if o.residue {
        Expr::call(Symbol::SEQUENCE, [imfs, as_list(residue)])
    } else {
        imfs
    };
    ev.step("emd", e, &result);
    Ok(Some(result))
}
