//! Interpolation (CAS-DOC-27 §2–§4): Xcas's `lagrange` and `interp`,
//! `spline` and `thiele`, as `InterpolatingPolynomial`, `Spline` and
//! `ThieleInterpolation`.
//!
//! - **`InterpolatingPolynomial`** (CAS-328), the polynomial of degree
//!   `< n` through `n` points, by Newton's divided differences
//!   `c_k = f[x_0, …, x_k]`, exactly, in Newton's nested form
//!   `c_0 + (t − x_0)(c_1 + (t − x_1)(c_2 + …))`, as giac's
//!   `(x-2+4)*(x-1)+4` and Mathematica's `4 + (−1 + x)(4 + (−2 + x))`.
//!   - The forms are giac's `(xs, ys)`, `(xs, ys, t)` and `([xs, ys])`, and
//!     Mathematica's `({{x_0, y_0}, …}, t)`; the variable is `x` when left
//!     out.
//! - **`Spline[xs, ys, t, d]`** (CAS-329), giac's natural spline of degree
//!   `d ≥ 1` through `n + 1` points of increasing abscissas: `n`
//!   polynomials `s_i = Σ_k a_{i,k} (t − x_i)^k`, one per interval.
//!   - Each meets its two points, and the derivatives of orders
//!     `1 … d − 1` agree at the inner points.
//!   - giac's end conditions, read from its answers for `d` from 2 to 5:
//!     the derivatives of orders `⌈d/2⌉ … d − 1` vanish at `x_0`, and those
//!     of orders `⌊d/2⌋ + 1 … d − 1` at `x_n`. For `d = 3` these are the
//!     natural cubic spline's `s″ = 0`.
//!   - The `n(d + 1)` coefficients solve one linear system (`LinSolve`).
//! - **`ThieleInterpolation[xs, ys, t]`** (CAS-330), Thiele's continued
//!   fraction `a_0 + (t − x_0)/(a_1 + (t − x_1)/(a_2 + …))`, its `a_k` the
//!   inverse differences `φ_k(x_k)`, where `φ_0 = y` and
//!   `φ_k(x_i) = (x_i − x_{k−1})/(φ_{k−1}(x_i) − φ_{k−1}(x_{k−1}))`; put
//!   together, as giac's `(-11*x+6)/(x-6)`. Where `φ_{k−1}` is constant on
//!   the points left, the fraction ends at `a_{k−1}`: `1/x` through four of
//!   its points, as giac answers, and a constant. Where it is so on some
//!   of them only, the call stays.
//!
//! Lists of unequal lengths, a repeated abscissa (a division by 0), and
//! for `Spline` abscissas not shown increasing, leave the call. Each call
//! that answers takes one step, its builtin's name.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of interpolation: points x and y, the variable t, the degree d, coefficients a and c"
)]

use super::solve::linear_system;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The points at most, so that the systems stay small.
const MAX_POINTS: usize = 256;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The abscissas, the values and the variable of a call's points
/// (the module's forms); `None` for another form.
fn points(args: &[Expr]) -> Option<(Vec<Expr>, Vec<Expr>, Expr)> {
    let x = || Expr::symbol("x");
    let (xs, ys, t) = match args {
        [xs, ys] if is_list(xs) && is_list(ys) => (xs.args().to_vec(), ys.args().to_vec(), x()),
        [xs, ys, t] if is_list(xs) && is_list(ys) => {
            (xs.args().to_vec(), ys.args().to_vec(), t.clone())
        }
        [pairs, t] if is_list(pairs) && t.as_symbol().is_some() => {
            let mut xs = Vec::with_capacity(pairs.args().len());
            let mut ys = Vec::with_capacity(pairs.args().len());
            for pair in pairs.args() {
                match pair.args() {
                    [a, b] if is_list(pair) => {
                        xs.push(a.clone());
                        ys.push(b.clone());
                    }
                    _ => return None,
                }
            }
            (xs, ys, t.clone())
        }
        [both] if is_list(both) => match both.args() {
            [xs, ys] if is_list(xs) && is_list(ys) => (xs.args().to_vec(), ys.args().to_vec(), x()),
            _ => return None,
        },
        _ => return None,
    };
    (xs.len() == ys.len() && !xs.is_empty() && xs.len() <= MAX_POINTS).then_some((xs, ys, t))
}

/// `a/b` evaluated; `None` when `b` is 0.
fn quotient(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<Expr>, EvalError> {
    let b = ev.eval(b)?;
    if b == Expr::from(0) {
        return Ok(None);
    }
    ev.eval(&(a.clone() / b)).map(Some)
}

pub(super) fn interpolating_polynomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys, t)) = points(e.args()) else {
        return Ok(None);
    };
    // Divided differences, in place: c_i = f[x_{i−j}, …, x_i].
    let mut c = ys;
    let n = c.len();
    for j in 1..n {
        for i in (j..n).rev() {
            let Some(d) = quotient(
                ev,
                &(c[i].clone() - c[i - 1].clone()),
                &(xs[i].clone() - xs[i - j].clone()),
            )?
            else {
                return Ok(None);
            };
            c[i] = d;
        }
    }
    let mut p = c[n - 1].clone();
    for i in (0..n - 1).rev() {
        p = c[i].clone() + (t.clone() - xs[i].clone()) * p;
    }
    let p = ev.eval(&p)?;
    ev.step("InterpolatingPolynomial", e, &p);
    Ok(Some(p))
}

/// `k!/(k − j)!`, the coefficient of the `j`-th derivative of `u^k` at
/// `u^(k−j)`.
fn falling(k: usize, j: usize) -> i64 {
    ((k - j + 1)..=k)
        .map(|m| i64::try_from(m).expect("a small degree"))
        .product()
}

pub(super) fn spline(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [xs, ys, t, d] = e.args() else {
        return Ok(None);
    };
    let Some((xs, ys, t)) = points(&[xs.clone(), ys.clone(), t.clone()]) else {
        return Ok(None);
    };
    let Some(d) = d
        .as_i64()
        .and_then(|d| usize::try_from(d).ok())
        .filter(|d| (1..=16).contains(d))
    else {
        return Ok(None);
    };
    let n = xs.len() - 1;
    if n == 0 || n * (d + 1) > 4 * MAX_POINTS {
        return Ok(None);
    }
    for pair in xs.windows(2) {
        let less = ev.eval(&Expr::call(
            Symbol::LESS,
            [pair[0].clone(), pair[1].clone()],
        ))?;
        if less.as_symbol() != Some(Symbol::TRUE) {
            return Ok(None);
        }
    }
    let unknown = |i: usize, k: usize| Expr::symbol(&format!("Spline$a{i}k{k}"));
    let mut widths = Vec::with_capacity(n);
    for i in 0..n {
        widths.push(ev.eval(&(xs[i + 1].clone() - xs[i].clone()))?);
    }
    // The j-th derivative of s_i at x_{i+1}.
    let derivative_at_end = |i: usize, j: usize| -> Expr {
        Expr::call(
            Symbol::PLUS,
            (j..=d)
                .map(|k| {
                    Expr::from(falling(k, j))
                        * unknown(i, k)
                        * widths[i].pow(i64::try_from(k - j).expect("a small degree"))
                })
                .collect::<Vec<_>>(),
        )
    };
    let equation = |l: Expr, r: Expr| Expr::call(Symbol::EQUATION, [l, r]);
    let mut equations = Vec::with_capacity(n * (d + 1));
    for i in 0..n {
        equations.push(equation(unknown(i, 0), ys[i].clone()));
        equations.push(equation(derivative_at_end(i, 0), ys[i + 1].clone()));
        if i + 1 < n {
            for j in 1..d {
                equations.push(equation(
                    derivative_at_end(i, j),
                    Expr::from(falling(j, j)) * unknown(i + 1, j),
                ));
            }
        }
    }
    for j in d.div_ceil(2)..d {
        equations.push(equation(unknown(0, j), Expr::from(0)));
    }
    for j in (d / 2 + 1)..d {
        equations.push(equation(derivative_at_end(n - 1, j), Expr::from(0)));
    }
    let unknowns: Vec<Expr> = (0..n)
        .flat_map(|i| (0..=d).map(move |k| (i, k)))
        .map(|(i, k)| unknown(i, k))
        .collect();
    let mut evaluated = Vec::with_capacity(equations.len());
    for eq in equations {
        evaluated.push(ev.eval(&eq)?);
    }
    let Some(found) = linear_system(ev, e, &evaluated, &unknowns)? else {
        return Ok(None);
    };
    let [values] = found.args() else {
        return Ok(None);
    };
    if values.args().len() != unknowns.len() {
        return Ok(None);
    }
    let mut pieces = Vec::with_capacity(n);
    for (x, a) in xs.iter().zip(values.args().chunks(d + 1)) {
        let u = t.clone() - x.clone();
        let terms: Vec<Expr> = a
            .iter()
            .enumerate()
            .map(|(k, c)| c.clone() * u.pow(i64::try_from(k).expect("a small degree")))
            .collect();
        pieces.push(Expr::call(Symbol::PLUS, terms));
    }
    let result = ev.eval(&Expr::call(Symbol::LIST, pieces))?;
    ev.step("Spline", e, &result);
    Ok(Some(result))
}

pub(super) fn thiele_interpolation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys, t)) = points(e.args()) else {
        return Ok(None);
    };
    let n = xs.len();
    let mut phi = ys;
    let mut a = vec![phi[0].clone()];
    for k in 1..n {
        let mut denominators = Vec::with_capacity(n - k);
        for value in &phi[k..] {
            denominators.push(ev.eval(&(value.clone() - phi[k - 1].clone()))?);
        }
        let zero = Expr::from(0);
        // All 0: φ_{k−1} is constant on the rest, and the fraction ends.
        if denominators.iter().all(|d| *d == zero) {
            break;
        }
        if denominators.contains(&zero) {
            return Ok(None);
        }
        let mut next = phi.clone();
        for (i, d) in (k..n).zip(denominators) {
            next[i] = ev.eval(&((xs[i].clone() - xs[k - 1].clone()) / d))?;
        }
        phi = next;
        a.push(phi[k].clone());
    }
    let mut r = a[a.len() - 1].clone();
    for k in (0..a.len() - 1).rev() {
        r = a[k].clone() + (t.clone() - xs[k].clone()) / r;
    }
    let result = ev.eval(&Expr::call(Symbol::TOGETHER, [r]))?;
    ev.step("ThieleInterpolation", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
