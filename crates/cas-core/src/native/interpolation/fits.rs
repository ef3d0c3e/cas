//! Polynomial fits (CAS-597): Xcas's `fitpoly`, after giac's
//! optimization.cc (`fitpoly`, `_fitpoly`).
//!
//! - **`PolynomialFit[X, Y, x | x == a..b | a..b, opts]`**, or a two-column
//!   matrix of points in place of `X, Y`: the least-squares polynomial by
//!   Forsythe's orthogonal polynomials on the points, `q_0 = 1`,
//!   `q_1 = x − Σx_i/m`, `q_n = (x − α_n) q_{n−1} − β_n q_{n−2}`, with
//!   `α_n = Σ x_i q_{n−1}(x_i)²/‖q_{n−1}‖²` and
//!   `β_n = Σ x_i q_{n−1}(x_i) q_{n−2}(x_i)/‖q_{n−2}‖²`; the fit is
//!   `Σ_{i ≤ n} q_i ⟨y, q_i⟩/‖q_i‖²`.
//!   - The degree, giac's rule: the residuals `ρ_0 = ‖y‖² − ⟨y, q_0⟩²/m`,
//!     `ρ_n = ρ_{n−1} − ⟨y, q_n⟩²/‖q_n‖²`; the loop ends at `ρ_n = 0`, at
//!     `n = degree`, or where the last `w` residuals' population standard
//!     deviation over the mean of all is at most `tol` (then `n − w + 1`);
//!     `n < m` and `n ≤ limit + w`, then at most `m − 1` and `limit`.
//!     `‖q_n‖ = 0` (repeated abscissas) ends it at `n − 1`.
//!   - Everything exact over the rationals, `tol` giac's double `0.01` when
//!     left out taken exactly. A decimal among the points is taken as its
//!     exact rational and the answer is balls at the working precision.
//!   - A range keeps the points with `a ≤ x_i ≤ b`. With a variable, the
//!     expanded polynomial; without one, its coefficients, the highest
//!     first (giac prints its internal `poly1[…]`).
//! - **`PolynomialFit[f, x == a..b | a..b, opts]`**, `f` with one free
//!   symbol: Chebyshev's interpolation at the `N + 1` zeros
//!   `z_j = a + (b − a)(1 + cos θ_j)/2`, `θ_j = (2j + 1)π/(2(N + 1))`,
//!   `N = max(63, limit + 1)`: `c_n = 2/(N + 1) Σ_j f(z_j) cos(n θ_j)`.
//!   - The degree, giac's rule: the first `n` with
//!     `|f₂ − Σ_{i ≤ n} c_i²/(1 + [i = 0])| ≤ tol f₂/f₂ₐ`, decided on the
//!     doubles of the midpoints, where `f₂ = (2/π)∫₀^π f(z(θ))² dθ` and
//!     `f₂ₐ = ∫_a^b f²` (both by `NIntegrate`), or `degree`, at most
//!     `limit`; `tol` is `10⁻⁵` when left out.
//!   - `f(z_j)` by `NIntegrate`'s ball tree, `cos(n θ_j)` by `cos_pi` of
//!     the exact `n(2j + 1)/(2(N + 1))`; `Σ c_i T_i(2(x − a)/(b − a) − 1)`
//!     (`c_0` halved) expanded by the recurrence of `T` in balls at the
//!     working precision plus [`GUARD`] bits: each coefficient a ball that
//!     holds giac's formula's value on exact nodes.
//! - Options: `limit = k` (`limit`, 15 by default), `threshold = t` or a
//!   bare positive real (`tol`), `len = w` or `length = w` (the window `w`,
//!   5 by default, at most the number of points given), `degree = n`.
//!
//! What stays: giac's errors (unequal lists, fewer than two points, data
//! not real numbers, an option of another form or out of range, a range
//! not real or with `a ≥ b`, a function of no or several symbols, a
//! function without a range), a range that keeps no point, a variable
//! other than the function's, a node where `f` is not a finite real, an
//! integral `NIntegrate` leaves, and any size past the caps: [`MAX_POINTS`]
//! points, [`MAX_DEGREE`] for `limit` and `degree`, [`MAX_FIT_WORK`] for
//! `m (n + 1)` of a data fit (`n` the most steps) and nodes `· (limit + 1)`.
//! Each call that answers takes one step, its builtin's name.

use cas_num::ball::{Ball, ComplexBall, Float};
use cas_num::{Integer, Rational};

use super::interpolants::matrix_rows;
use super::splines::{number, real};
use super::{MAX_POINTS, Native, is_list};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::native::nintegrate::analytic::{Mode, Node};
use crate::symbol::Symbol;

/// The highest `limit` and `degree`.
const MAX_DEGREE: usize = 64;

/// The most `m (n + 1)` of a data fit, and nodes `· (limit + 1)` of a
/// function's.
const MAX_FIT_WORK: usize = 1 << 13;

/// The guard bits of a function's balls.
const GUARD: u64 = 32;

/// A polynomial, its coefficients from the constant up.
type Poly = Vec<Rational>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// Whether `e` is a variable: a symbol not a builtin's, `len` aside (an
/// option's name in giac).
fn is_variable(e: &Expr) -> bool {
    e.as_symbol().is_some_and(|s| !s.is_builtin()) && *e != Expr::symbol("len")
}

/// The argument after the data or the function.
enum Domain {
    /// Not a variable or a range: an option.
    Not,
    /// giac's error: `x == v` with `v` not a range, or `v == a..b` with
    /// `v` not a variable.
    Bad,
    Is(Option<Expr>, Option<(Expr, Expr)>),
}

fn domain(v: &Expr) -> Domain {
    if is_variable(v) {
        return Domain::Is(Some(v.clone()), None);
    }
    let range = |r: &Expr| match (head(r), r.args()) {
        (Some(Symbol::RANGE), [a, b]) => Some((a.clone(), b.clone())),
        _ => None,
    };
    if let Some(r) = range(v) {
        return Domain::Is(None, Some(r));
    }
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [lh, rh]) = (head(v), v.args()) else {
        return Domain::Not;
    };
    match (is_variable(lh), range(rh)) {
        (true, Some(r)) => Domain::Is(Some(lh.clone()), Some(r)),
        (false, None) => Domain::Not,
        _ => Domain::Bad,
    }
}

/// The options of a fit.
struct Options {
    limit: usize,
    tol: Option<Rational>,
    window: usize,
    degree: Option<usize>,
}

/// A positive integer, or a nonnegative one with `zero`, as a size.
fn size(e: &Expr, zero: bool) -> Option<usize> {
    e.as_i64()
        .filter(|&k| k > 0 || (zero && k == 0))
        .and_then(|k| usize::try_from(k).ok())
}

/// A positive real number exactly.
fn positive(e: &Expr) -> Option<Rational> {
    real(e).map(|(q, _)| q).filter(|q| *q > Rational::from(0))
}

/// The options; `None` for another (giac's errors) or past the caps.
fn options(given: &[Expr]) -> Option<Options> {
    let mut o = Options {
        limit: 15,
        tol: None,
        window: 5,
        degree: None,
    };
    for option in given {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [lh, rh]) = (head(option), option.args())
        else {
            o.tol = Some(positive(option)?);
            continue;
        };
        match lh.as_symbol() {
            Some(Symbol::LIMIT) => o.limit = size(rh, false).filter(|&k| k <= MAX_DEGREE)?,
            Some(Symbol::THRESHOLD) => o.tol = Some(positive(rh)?),
            Some(Symbol::LENGTH) => o.window = size(rh, false)?,
            Some(_) if *lh == Expr::symbol("len") => o.window = size(rh, false)?,
            Some(Symbol::EXPONENT) => {
                o.degree = Some(size(rh, true).filter(|&k| k <= MAX_DEGREE)?);
            }
            _ => return None,
        }
    }
    Some(o)
}

fn count(n: usize) -> Rational {
    Rational::from(i64::try_from(n).expect("a size in memory"))
}

fn dot(u: &[Rational], v: &[Rational]) -> Rational {
    u.iter()
        .zip(v)
        .fold(Rational::from(0), |s, (a, b)| &s + &(a * b))
}

/// giac's `fitpoly` on the points `(xs, ys)`, `m ≥ 1`, with the tolerance,
/// `d` the window when positive and minus the degree otherwise, and the
/// limit: the coefficients, the constant first.
pub(super) fn fit_data(
    xs: &[Rational],
    ys: &[Rational],
    tol: &Rational,
    d: i64,
    limit: usize,
) -> Vec<Rational> {
    let zero = Rational::from(0);
    let m = xs.len();
    let bound = limit + usize::try_from(d.max(0)).expect("a window in memory");
    let mut q: Vec<Poly> = Vec::new();
    let mut qy: Vec<Vec<Rational>> = Vec::new();
    let (mut fq, mut qq, mut rho): (Vec<Rational>, Vec<Rational>, Vec<Rational>) =
        (Vec::new(), Vec::new(), Vec::new());
    let mut n = 0;
    while n < m && n <= bound {
        let (poly, values) = if n == 0 {
            (vec![Rational::from(1)], vec![Rational::from(1); m])
        } else if n == 1 {
            let mean = &xs.iter().fold(zero.clone(), |s, x| &s + x) / &qq[0];
            let values = xs.iter().map(|x| x - &mean).collect();
            (vec![-mean, Rational::from(1)], values)
        } else {
            let (last, before) = (&qy[n - 1], &qy[n - 2]);
            let xq: Vec<Rational> = xs.iter().zip(last).map(|(x, v)| x * v).collect();
            let alpha = &dot(&xq, last) / &qq[n - 1];
            let beta = &dot(&xq, before) / &qq[n - 2];
            let mut poly = vec![zero.clone(); n + 1];
            for (i, c) in q[n - 1].iter().enumerate() {
                poly[i + 1] = &poly[i + 1] + c;
                poly[i] = &poly[i] - &(&alpha * c);
            }
            for (i, c) in q[n - 2].iter().enumerate() {
                poly[i] = &poly[i] - &(&beta * c);
            }
            let values = xs
                .iter()
                .zip(last.iter().zip(before))
                .map(|(x, (u, v))| &(&(x - &alpha) * u) - &(&beta * v))
                .collect();
            (poly, values)
        };
        let norm = dot(&values, &values);
        if norm.is_zero() {
            // Repeated abscissas: q_n vanishes at the points.
            n -= 1;
            break;
        }
        let f = dot(ys, &values);
        let previous = if n == 0 {
            dot(ys, ys)
        } else {
            rho[n - 1].clone()
        };
        let r = &previous - &(&(&f * &f) / &norm);
        rho.push(if r < zero { zero.clone() } else { r });
        q.push(poly);
        qy.push(values);
        fq.push(f);
        qq.push(norm);
        let at = i64::try_from(n).expect("a degree in memory");
        if rho[n].is_zero() || (d <= 0 && at == -d) {
            break;
        }
        let w = usize::try_from(d).unwrap_or(0);
        if w > 0 && n + 1 >= w {
            let window = &rho[n + 1 - w..];
            let mean_w = &window.iter().fold(zero.clone(), |s, v| &s + v) / &count(w);
            let var = &window.iter().fold(zero.clone(), |s, v| {
                let e = v - &mean_w;
                &s + &(&e * &e)
            }) / &count(w);
            let mean = &rho.iter().fold(zero.clone(), |s, v| &s + v) / &count(n + 1);
            if &(&(tol * tol) * &mean) * &mean >= var {
                n -= w - 1;
                break;
            }
        }
        n += 1;
    }
    let n = n.min(m - 1).min(limit);
    let mut coefficients = vec![zero; n + 1];
    for i in 0..=n {
        let scale = &fq[i] / &qq[i];
        for (j, c) in q[i].iter().enumerate() {
            coefficients[j] = &coefficients[j] + &(&scale * c);
        }
    }
    coefficients
}

/// A real number as a finite ball at `prec` bits; `None` for another.
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

/// The double of a real integral `NIntegrate` gives; `None` when it stays.
fn integral(
    ev: &mut Eval<'_>,
    g: Expr,
    x: &Expr,
    a: Expr,
    b: Expr,
) -> Result<Option<f64>, EvalError> {
    let value = ev.eval(&Expr::call(Symbol::NINTEGRATE, [g, x.clone(), a, b]))?;
    Ok(match value.view() {
        View::Ball(v) if v.is_finite() => Some(v.mid().to_f64()),
        _ => None,
    })
}

/// giac's `fitpoly` of `f` in `x` on `[a, b]` (the module's function
/// path): the polynomial with ball coefficients; `None` where it stays.
fn fit_function(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    (a, b): (&Expr, &Expr),
    tol: f64,
    limit: usize,
    degree: Option<usize>,
) -> Result<Option<Expr>, EvalError> {
    let wp = ev.context().working_precision() + GUARD;
    let (Some(lo), Some(hi)) = (ball(ev, a, wp), ball(ev, b, wp)) else {
        return Ok(None);
    };
    if lo.cmp(&hi) != Some(std::cmp::Ordering::Less) {
        return Ok(None);
    }
    let Some(tree) = Node::compile(ev.context(), f, x, wp) else {
        return Ok(None);
    };
    let nodes = 64.max(limit + 2);
    let width = hi.sub(&lo, wp);
    let one = Ball::from_integer(&Integer::from(1));
    let half = |v: &Ball| v.mul_2exp(&Integer::from(-1));
    let denominator = Integer::from(2 * i64::try_from(nodes).expect("a node count"));
    let angle = |k: usize| {
        let k = Integer::from(i64::try_from(k).expect("a node index"));
        Ball::from_rational(&Rational::new(k, denominator.clone()), wp).cos_pi(wp)
    };
    let mut values = Vec::with_capacity(nodes);
    for j in 0..nodes {
        let z = lo.add(&half(&width.mul(&one.add(&angle(2 * j + 1), wp), wp)), wp);
        let Some(v) = tree
            .eval(&ComplexBall::from_real(z), Mode::Real, wp)
            .filter(|v| v.is_real() && v.re().is_finite())
        else {
            return Ok(None);
        };
        values.push(v.re().clone());
    }
    // The integrals after the nodes, the cheaper refusal first.
    let (mut f2, mut tol) = (0.0, tol);
    if degree.is_none() {
        let shifted = a.clone()
            + (b.clone() - a.clone()) * (Expr::from(1) + Expr::call(Symbol::COS, [x.clone()]))
                / Expr::from(2);
        let g = f.substitute(&[(x.clone(), shifted)]).pow(2);
        let pi = Expr::from(Symbol::PI);
        let Some(weighted) = integral(ev, g, x, Expr::from(0), pi)? else {
            return Ok(None);
        };
        let Some(plain) = integral(ev, f.clone().pow(2), x, a.clone(), b.clone())? else {
            return Ok(None);
        };
        f2 = 2.0 * weighted / std::f64::consts::PI;
        tol = tol * f2 / plain;
    }
    let scale = Ball::from_rational(
        &Rational::new(
            Integer::from(2),
            Integer::from(i64::try_from(nodes).expect("a node count")),
        ),
        wp,
    );
    let mut c: Vec<Ball> = Vec::new();
    let mut fp = 0.0;
    for n in 0..=limit {
        let mut sum = Ball::ZERO;
        for (j, v) in values.iter().enumerate() {
            sum = sum.add(&v.mul(&angle(n * (2 * j + 1)), wp), wp);
        }
        let cf = sum.mul(&scale, wp);
        let mid = cf.mid().to_f64();
        fp += if n == 0 { mid * mid / 2.0 } else { mid * mid };
        c.push(cf);
        if degree == Some(n) || (degree.is_none() && tol >= (f2 - fp).abs()) {
            break;
        }
    }
    // Σ c_i T_i(u), u = α x + β, by T_{k+1} = 2u T_k − T_{k−1}.
    let alpha = Ball::from_integer(&Integer::from(2)).div(&width, wp);
    let beta = lo.add(&hi, wp).div(&width, wp).neg();
    let top = c.len() - 1;
    let mut p = vec![Ball::ZERO; top + 1];
    let (mut previous, mut current) = (vec![one], vec![beta.clone(), alpha.clone()]);
    p[0] = half(&c[0]);
    for (k, ck) in c.iter().enumerate().skip(1) {
        for (pi, t) in p.iter_mut().zip(&current) {
            *pi = pi.add(&ck.mul(t, wp), wp);
        }
        if k < top {
            let mut next = vec![Ball::ZERO; current.len() + 1];
            for (i, t) in current.iter().enumerate() {
                let t2 = t.mul_2exp(&Integer::from(1));
                next[i] = next[i].add(&t2.mul(&beta, wp), wp);
                next[i + 1] = next[i + 1].add(&t2.mul(&alpha, wp), wp);
            }
            for (i, t) in previous.iter().enumerate() {
                next[i] = next[i].sub(t, wp);
            }
            previous = std::mem::replace(&mut current, next);
        }
    }
    let terms = p
        .into_iter()
        .enumerate()
        .map(|(k, v)| Expr::from(v) * x.pow(i64::try_from(k).expect("a degree below MAX_DEGREE")));
    Ok(Some(ev.eval(&Expr::call(Symbol::PLUS, terms))?))
}

/// The one variable of `e`, symbols not a builtin's; `None` for none or
/// several.
fn the_variable(e: &Expr) -> Option<Expr> {
    fn collect(e: &Expr, found: &mut Vec<Expr>) {
        if e.as_symbol().is_some_and(|s| !s.is_builtin()) {
            if !found.contains(e) {
                found.push(e.clone());
            }
            return;
        }
        for a in e.args() {
            collect(a, found);
        }
    }
    let mut found = Vec::new();
    collect(e, &mut found);
    match found.as_slice() {
        [x] => Some(x.clone()),
        _ => None,
    }
}

/// `v` as a double.
fn double(v: &Rational) -> f64 {
    Ball::from_rational(v, 64).mid().to_f64()
}

/// The module's data path on the points `(xs, ys)`, with the variable,
/// the range and the options; `None` where it stays.
fn fit_points(
    ev: &mut Eval<'_>,
    (xs, ys): (Vec<Expr>, Vec<Expr>),
    var: Option<Expr>,
    range: Option<&(Expr, Expr)>,
    o: &Options,
) -> Result<Option<Expr>, EvalError> {
    let m = xs.len();
    if m != ys.len() || !(2..=MAX_POINTS).contains(&m) {
        return Ok(None);
    }
    let bounds = match range {
        Some((a, b)) => {
            let (Some((a, _)), Some((b, _))) = (real(a), real(b)) else {
                return Ok(None);
            };
            if a >= b {
                return Ok(None);
            }
            Some((a, b))
        }
        None => None,
    };
    let (mut kept_x, mut kept_y, mut decimal) = (Vec::new(), Vec::new(), false);
    for (x, y) in xs.iter().zip(&ys) {
        let (Some((x, dx)), Some((y, dy))) = (real(x), real(y)) else {
            return Ok(None);
        };
        if bounds.as_ref().is_some_and(|(a, b)| x < *a || x > *b) {
            continue;
        }
        decimal |= dx || dy;
        kept_x.push(x);
        kept_y.push(y);
    }
    let kept = kept_x.len();
    if kept == 0 {
        return Ok(None);
    }
    let window = o.window.min(m);
    let d = match o.degree {
        Some(k) => -i64::try_from(k).expect("a degree below MAX_DEGREE"),
        None => i64::try_from(window).expect("a window below MAX_POINTS"),
    };
    let steps = kept.min(o.limit + usize::try_from(d.max(0)).expect("a window") + 1);
    if kept * steps > MAX_FIT_WORK {
        return Ok(None);
    }
    let tol = o.tol.clone().unwrap_or_else(|| {
        Float::from_f64(0.01)
            .expect("a finite double")
            .to_rational()
    });
    let coefficients = fit_data(&kept_x, &kept_y, &tol, d, o.limit);
    let prec = ev.context().working_precision();
    Ok(Some(match var {
        Some(x) => {
            let terms: Vec<Expr> = coefficients
                .into_iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(|(k, c)| {
                    number(c, decimal, prec)
                        * x.pow(i64::try_from(k).expect("a degree below MAX_DEGREE"))
                })
                .collect();
            ev.eval(&Expr::call(Symbol::PLUS, terms))?
        }
        None => Expr::call(
            Symbol::LIST,
            coefficients
                .into_iter()
                .rev()
                .map(|c| number(c, decimal, prec)),
        ),
    }))
}

pub(in crate::native) fn polynomial_fit(ev: &mut Eval<'_>, e: &Expr) -> Native {
    // The data (or the function) and the arguments after it.
    let (data, rest) = match e.args() {
        [first, rest @ ..] if is_list(first) => {
            if let Some(rows) = matrix_rows(first) {
                if rows[0].args().len() != 2 {
                    return Ok(None);
                }
                let column =
                    |k: usize| -> Vec<Expr> { rows.iter().map(|r| r.args()[k].clone()).collect() };
                (Some((column(0), column(1))), rest)
            } else {
                let [ys, rest @ ..] = rest else {
                    return Ok(None);
                };
                if !is_list(ys) {
                    return Ok(None);
                }
                (Some((first.args().to_vec(), ys.args().to_vec())), rest)
            }
        }
        [_, rest @ ..] => (None, rest),
        [] => return Ok(None),
    };
    let (var, range, given) = match rest.split_first().map(|(v, more)| (domain(v), more)) {
        None | Some((Domain::Not, _)) => (None, None, rest),
        Some((Domain::Bad, _)) => return Ok(None),
        Some((Domain::Is(var, range), more)) => (var, range, more),
    };
    let Some(o) = options(given) else {
        return Ok(None);
    };
    let result = if let Some(data) = data {
        let Some(result) = fit_points(ev, data, var, range.as_ref(), &o)? else {
            return Ok(None);
        };
        result
    } else {
        let f = &e.args()[0];
        if f.is_number() || is_variable(f) {
            return Ok(None);
        }
        let Some(x) = the_variable(f) else {
            return Ok(None);
        };
        let Some((a, b)) = &range else {
            return Ok(None);
        };
        if var.is_some_and(|v| v != x) {
            return Ok(None);
        }
        let nodes = 64.max(o.limit + 2);
        if nodes * (o.limit + 1) > MAX_FIT_WORK {
            return Ok(None);
        }
        let tol = o.tol.as_ref().map_or(1e-5, double);
        let Some(result) = fit_function(ev, f, &x, (a, b), tol, o.limit, o.degree)? else {
            return Ok(None);
        };
        result
    };
    ev.step("PolynomialFit", e, &result);
    Ok(Some(result))
}
