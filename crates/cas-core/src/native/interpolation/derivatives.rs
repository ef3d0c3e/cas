//! Numerical derivatives (CAS-595): Xcas's `numdiff` and `nDeriv`, after
//! giac's optimization.cc (`_numdiff`) and ti89.cc (`_nDeriv`).
//!
//! - **`NumericalDerivative[xs, ys, x0, m_1, m_2, …]`** (`_numdiff`): the
//!   `m`-th derivative at `x0` of the polynomial of degree `< n` through
//!   the `n` points `(x_k, y_k)`; the orders are `[1]` when left out, each a
//!   nonnegative integer. One order answers a value, several a list.
//!   - giac takes Sadiq and Viswanath's finite-difference weights `w_{k,m}`
//!     and answers `Σ_k w_{k,m} y_k`. Those weights are the `m`-th
//!     derivatives at `x0` of the Lagrange basis polynomials `ℓ_k`: the
//!     interpolant is `p = Σ_k y_k ℓ_k`, so `Σ_k w_{k,m} y_k = p^(m)(x0)`,
//!     which is computed here from Newton's form instead. The divided
//!     differences `c_k` give `p = c_0 + (t − x_0)(c_1 + (t − x_1)(…))`;
//!     its Taylor coefficients `d_j` at `x0`, up to the largest order `M`,
//!     come by multiplying out the nesting from the inside,
//!     `d_j ← (x0 − x_i) d_j + d_{j−1}`, `d_0 ← (x0 − x_i) d_0 + c_i`; then
//!     `p^(m)(x0) = m! d_m`. An order of `n` or more answers 0.
//!   - Input all numeric with a decimal is taken numerically and answers
//!     numerically, as giac's `evalf` branch; other input answers exactly.
//!   - Fewer than three arguments, lists of unequal lengths, an order not a
//!     nonnegative integer (giac's errors), a repeated abscissa (giac's
//!     `undef`), and more than 12 points with an abscissa or `x0` not a
//!     number (the unexpanded products would nest exponentially) leave
//!     the call.
//! - **`NDeriv[f, x, h]`** (`_nDeriv`): the central difference quotient
//!   `(f(x + h) − f(x − h))/(2h)`, `h = 1/1000` when left out (giac's
//!   `0.001`, here exact, as giac's answer otherwise). A variable not a
//!   symbol, and `h` 0, leave the call.
//!
//! Each call that answers takes one step, its builtin's name.

use super::interpolants::has_ball;
use super::{MAX_POINTS, Native, divided_differences, is_list, vanishes};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The most `n (M + 1)` of a numerical derivative: the products of the
/// Taylor coefficients.
const MAX_DERIVATIVE_WORK: usize = 1 << 13;

/// The most points when an abscissa or `x0` is not a number.
const MAX_SYMBOLIC_POINTS: usize = 12;

pub(in crate::native) fn numerical_derivative(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [xs, ys, x0, orders @ ..] = e.args() else {
        return Ok(None);
    };
    if !is_list(xs) || !is_list(ys) {
        return Ok(None);
    }
    let n = xs.args().len();
    if n != ys.args().len() || n == 0 || n > MAX_POINTS {
        return Ok(None);
    }
    let mut ms = Vec::with_capacity(orders.len().max(1));
    for m in orders {
        match m.as_i64().and_then(|m| usize::try_from(m).ok()) {
            Some(m) => ms.push(m),
            None => return Ok(None),
        }
    }
    if ms.is_empty() {
        ms.push(1);
    }
    // Orders of n or more answer 0: clamp before any work.
    let top = ms.iter().copied().max().unwrap_or(0).min(n - 1);
    if n * (top + 1) > MAX_DERIVATIVE_WORK {
        return Ok(None);
    }
    let all = || xs.args().iter().chain(ys.args()).chain([x0]);
    let numeric = all().all(Expr::is_numeric) && all().any(has_ball);
    let mut xs = xs.args().to_vec();
    let mut ys = ys.args().to_vec();
    let mut x0 = x0.clone();
    if numeric {
        for v in xs.iter_mut().chain(ys.iter_mut()).chain([&mut x0]) {
            *v = ev.eval(&Expr::call(Symbol::N, [v.clone()]))?;
        }
    }
    // Products of gaps that are not numbers stay unexpanded and nest: the
    // Taylor coefficients grow as binomials in n, so few points only.
    let number = |v: &Expr| v.to_rational().is_some() || matches!(v.view(), View::Ball(_));
    if n > MAX_SYMBOLIC_POINTS && !(number(&x0) && xs.iter().all(number)) {
        return Ok(None);
    }
    let Some(c) = divided_differences(ev, &xs, ys)? else {
        return Ok(None);
    };
    let mut d = vec![Expr::from(0); top + 1];
    for (ci, xi) in c.iter().zip(&xs).rev() {
        let gap = ev.eval(&(x0.clone() - xi.clone()))?;
        for k in (1..=top).rev() {
            d[k] = ev.eval(&(gap.clone() * d[k].clone() + d[k - 1].clone()))?;
        }
        d[0] = ev.eval(&(ci.clone() + gap * d[0].clone()))?;
    }
    let mut values = Vec::with_capacity(ms.len());
    for &m in &ms {
        let value = if m >= n {
            Expr::from(0)
        } else {
            let order = Expr::from(i64::try_from(m).expect("an order below MAX_POINTS"));
            Expr::call(Symbol::FACTORIAL, [order]) * d[m].clone()
        };
        values.push(value);
    }
    let result = match <[Expr; 1]>::try_from(values) {
        Ok([value]) => value,
        Err(values) => Expr::call(Symbol::LIST, values),
    };
    let result = if numeric {
        Expr::call(Symbol::N, [result])
    } else {
        result
    };
    let result = ev.eval(&result)?;
    ev.step("NumericalDerivative", e, &result);
    Ok(Some(result))
}

pub(in crate::native) fn n_deriv(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, x, h) = match e.args() {
        [f, x] => (f, x, Expr::from(1) / Expr::from(1000)),
        [f, x, h] => (f, x, h.clone()),
        _ => return Ok(None),
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let h = ev.eval(&h)?;
    if vanishes(&h) {
        return Ok(None);
    }
    let ahead = f.substitute(&[(x.clone(), x.clone() + h.clone())]);
    let behind = f.substitute(&[(x.clone(), x.clone() - h.clone())]);
    let result = ev.eval(&((ahead - behind) / (Expr::from(2) * h)))?;
    ev.step("NDeriv", e, &result);
    Ok(Some(result))
}
