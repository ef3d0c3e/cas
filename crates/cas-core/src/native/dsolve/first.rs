//! First-order equations that are not linear (CAS-DOC-16 §4.2–§4.5,
//! CAS-260), each written `M(x, y) + N(x, y) y′ = 0` with `y′ = F = −M/N`,
//! in this order:
//!
//! - **Separable**, `F = f(x) g(y)`, the factors of `F` apart once its
//!   exponentials of sums are products: `H(x) − G(y) = c_0` for `H = ∫f`,
//!   `G = ∫ dy/g`, giac's constant (`y′ = y²` is `1/(−x + c_0)`).
//! - **Bernoulli**, `F = −p y + q yᵏ`, `k ≠ 0, 1` free of `x`: `v = y^(1−k)`
//!   is linear, `v′ + (1 − k) p v = (1 − k) q`, and `y = v^(1/(1−k))`, both
//!   signs where that root is even.
//! - **Homogeneous**, `F(x, v x)` free of `x`: separable in `v = y/x`,
//!   `v′ = (F(1, v) − v)/x`.
//! - **Exact**, `∂M/∂y = ∂N/∂x`, or so once multiplied by `μ(x) =
//!   e^(∫ (M_y − N_x)/N dx)` or `μ(y) = e^(∫ (N_x − M_y)/M dy)`, the ratio
//!   free of the other variable: the potential `Φ = ∫ M dx + h(y)`, `h′ = N −
//!   ∂(∫M dx)/∂y`, and `Φ = c_0`.
//!
//! The relation found is solved for `y` by `Solve`, each root a branch of
//! the list answered; integrals that stay, or a relation `Solve` does not
//! invert, leave the method. Logarithms of absolute values lose them, which
//! the constant absorbs.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equations: M + N y′ = 0, y′ = F(x, y)"
)]

use super::{Ode, exponential, first_order_linear, without_absolute_values};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The branches of a first-order equation of degree 1 in `y′` that is not
/// linear, by the first method that applies; `None` when none does.
pub(super) fn nonlinear(
    ev: &mut Eval<'_>,
    ode: &Ode,
    x: &Expr,
    c: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let (y, dy) = (&ode.ys[0], &ode.ys[1]);
    let n = ev.eval(&Expr::call(Symbol::D, [ode.residual.clone(), dy.clone()]))?;
    if n.contains(dy) || n == Expr::from(0) {
        return Ok(None);
    }
    let m = ev.eval(&ode.residual.substitute(&[(dy.clone(), Expr::from(0))]))?;
    let f = ev.eval(&Expr::call(Symbol::TOGETHER, [-(m.clone() / n.clone())]))?;
    if let Some(branches) = separable(ev, &f, x, y, c)? {
        return Ok(Some(branches));
    }
    if let Some(branches) = bernoulli(ev, &f, x, y, c)? {
        return Ok(Some(branches));
    }
    if let Some(branches) = homogeneous(ev, &f, x, y, c)? {
        return Ok(Some(branches));
    }
    exact(ev, &m, &n, x, y, c)
}

/// `∫ f dv` without absolute values in its logarithms; `None` when an
/// integral stays.
fn integral(ev: &mut Eval<'_>, f: &Expr, v: &Expr) -> Result<Option<Expr>, EvalError> {
    let i = ev.eval(&Expr::call(Symbol::INTEGRATE, [f.clone(), v.clone()]))?;
    if i.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(None);
    }
    Ok(Some(without_absolute_values(&i)))
}

/// The roots in `y` of `lhs = rhs`, as `Solve` finds them; `None` when it
/// does not.
fn roots(
    ev: &mut Eval<'_>,
    lhs: Expr,
    rhs: Expr,
    y: &Expr,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let equation = Expr::call(Symbol::EQUATION, [lhs, rhs]);
    let found = ev.eval(&Expr::call(Symbol::SOLVE, [equation, y.clone()]))?;
    let solved = found.head() == Some(&Expr::from(Symbol::LIST))
        && !found.args().is_empty()
        && found
            .args()
            .iter()
            .all(|r| !r.contains(y) && !r.contains(&Expr::from(Symbol::SOLVE)));
    Ok(solved.then(|| found.args().to_vec()))
}

/// The list of the branches.
fn branches(roots: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, roots)
}

/// `e` with each exponential of a sum a product: `e^(x + y)` is `e^x e^y`,
/// so that its factors can be told apart.
fn split_exponentials(e: &Expr) -> Expr {
    e.map_bottom_up(|part| {
        let exponent = match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(Symbol::POWER), [b, s]) if b.as_symbol() == Some(Symbol::E) => s,
            (Some(Symbol::EXP), [s]) => s,
            _ => return part.clone(),
        };
        if exponent.head() != Some(&Expr::from(Symbol::PLUS)) {
            return part.clone();
        }
        let factors = exponent
            .args()
            .iter()
            .map(|t| Expr::from(Symbol::E).pow(t.clone()));
        Expr::call(Symbol::TIMES, factors)
    })
}

/// `F = f(x) g(y)`: `H(x) − G(y) = c` solved for `y`.
fn separable(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    y: &Expr,
    c: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some((fx, gy)) = separate(ev, f, x, y)? else {
        return Ok(None);
    };
    let Some(h) = integral(ev, &fx, x)? else {
        return Ok(None);
    };
    let Some(g) = integral(ev, &(Expr::from(1) / gy), y)? else {
        return Ok(None);
    };
    Ok(roots(ev, g, h - c.clone(), y)?.map(branches))
}

/// `(f, g)` with `F = f(x) g(y)`, `g` holding `y`; `None` when a factor of
/// `F` holds both, or none holds `y`.
fn separate(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    y: &Expr,
) -> Result<Option<(Expr, Expr)>, EvalError> {
    let split = ev.eval(&split_exponentials(f))?;
    let factored = ev.eval(&Expr::call(Symbol::FACTOR, [split]))?;
    let factors = if factored.head() == Some(&Expr::from(Symbol::TIMES)) {
        factored.args().to_vec()
    } else {
        vec![factored]
    };
    let (mut in_x, mut in_y) = (Vec::new(), Vec::new());
    for factor in factors {
        match (factor.contains(x), factor.contains(y)) {
            (true, true) => return Ok(None),
            (_, true) => in_y.push(factor),
            _ => in_x.push(factor),
        }
    }
    if in_y.is_empty() {
        return Ok(None);
    }
    let product = |v: Vec<Expr>| Expr::call(Symbol::TIMES, v);
    Ok(Some((ev.eval(&product(in_x))?, ev.eval(&product(in_y))?)))
}

/// `(c, m)` for a term `c y^m` of a sum, `c` free of `y` and `m` of `x`
/// and `y`.
fn power_term(term: &Expr, x: &Expr, y: &Expr) -> Option<(Expr, Expr)> {
    let factors = if term.head() == Some(&Expr::from(Symbol::TIMES)) {
        term.args().to_vec()
    } else {
        vec![term.clone()]
    };
    let (mut coefficient, mut exponent) = (Vec::new(), Vec::new());
    for factor in factors {
        if !factor.contains(y) {
            coefficient.push(factor);
        } else if factor == *y {
            exponent.push(Expr::from(1));
        } else {
            match (factor.head().and_then(Expr::as_symbol), factor.args()) {
                (Some(Symbol::POWER), [b, m]) if b == y && !m.contains(x) && !m.contains(y) => {
                    exponent.push(m.clone());
                }
                _ => return None,
            }
        }
    }
    Some((
        Expr::call(Symbol::TIMES, coefficient),
        Expr::call(Symbol::PLUS, exponent),
    ))
}

/// `F = −p y + q yᵏ`: `v = y^(1−k)` solves a linear equation.
fn bernoulli(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    y: &Expr,
    c: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone()]))?;
    let terms = if expanded.head() == Some(&Expr::from(Symbol::PLUS)) {
        expanded.args().to_vec()
    } else {
        vec![expanded]
    };
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    for term in &terms {
        let Some((coefficient, m)) = power_term(term, x, y) else {
            return Ok(None);
        };
        let m = ev.eval(&m)?;
        match groups.iter_mut().find(|(e, _)| *e == m) {
            Some((_, cs)) => cs.push(coefficient),
            None => groups.push((m, vec![coefficient])),
        }
    }
    let one = Expr::from(1);
    let [(m1, c1), (m2, c2)] = groups.as_slice() else {
        return Ok(None);
    };
    let ((_, linear), (k, power)) = if *m1 == one {
        ((m1, c1), (m2, c2))
    } else if *m2 == one {
        ((m2, c2), (m1, c1))
    } else {
        return Ok(None);
    };
    if *k == Expr::from(0) {
        return Ok(None);
    }
    let sum = |cs: &[Expr]| Expr::call(Symbol::PLUS, cs.to_vec());
    // y′ = −p y + q yᵏ.
    let p = ev.eval(&-sum(linear))?;
    let q = ev.eval(&sum(power))?;
    let one_minus_k = ev.eval(&(one.clone() - k.clone()))?;
    let a0 = ev.eval(&(one_minus_k.clone() * p))?;
    let rhs = ev.eval(&(one_minus_k.clone() * q))?;
    let v = first_order_linear(ev, x, &a0, &one, &rhs, c)?;
    let root = ev.eval(&(one / one_minus_k))?;
    let y = ev.eval(&v.pow(root.clone()))?;
    let even = root
        .to_rational()
        .and_then(|r| r.denom().to_i64())
        .is_some_and(|d| d % 2 == 0);
    Ok(Some(branches(if even {
        let negative = ev.eval(&-y.clone())?;
        vec![y, negative]
    } else {
        vec![y]
    })))
}

/// `F(x, v x)` free of `x`: `v = y/x` is separable.
fn homogeneous(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    y: &Expr,
    c: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let v = Expr::from(Symbol::new("DSolve$v"));
    let at = f.substitute(&[(y.clone(), v.clone() * x.clone())]);
    let g = ev.eval(&Expr::call(Symbol::TOGETHER, [at]))?;
    if g.contains(x) {
        return Ok(None);
    }
    // v′ = (g(v) − v)/x.
    let slope = ev.eval(&Expr::call(Symbol::TOGETHER, [g - v.clone()]))?;
    if slope == Expr::from(0) {
        return Ok(None);
    }
    let Some(h) = integral(ev, &(Expr::from(1) / x.clone()), x)? else {
        return Ok(None);
    };
    let Some(gv) = integral(ev, &(Expr::from(1) / slope), &v)? else {
        return Ok(None);
    };
    let Some(vs) = roots(ev, gv, h - c.clone(), &v)? else {
        return Ok(None);
    };
    let mut ys = Vec::with_capacity(vs.len());
    for root in vs {
        ys.push(ev.eval(&(x.clone() * root))?);
    }
    Ok(Some(branches(ys)))
}

/// `M + N y′ = 0` exact, or so with an integrating factor in one variable:
/// its potential `Φ = c` solved for `y`.
fn exact(
    ev: &mut Eval<'_>,
    m: &Expr,
    n: &Expr,
    x: &Expr,
    y: &Expr,
    c: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let d = |ev: &mut Eval<'_>, e: &Expr, v: &Expr| {
        ev.eval(&Expr::call(Symbol::D, [e.clone(), v.clone()]))
    };
    let together = |ev: &mut Eval<'_>, e: Expr| ev.eval(&Expr::call(Symbol::TOGETHER, [e]));
    let my = d(ev, m, y)?;
    let nx = d(ev, n, x)?;
    let gap = together(ev, my - nx)?;
    let (m, n) = if gap == Expr::from(0) {
        (m.clone(), n.clone())
    } else {
        let in_x = together(ev, gap.clone() / n.clone())?;
        let in_y = together(ev, -(gap / m.clone()))?;
        let factor = if !in_x.contains(y) {
            integral(ev, &in_x, x)?
        } else if !in_y.contains(x) {
            integral(ev, &in_y, y)?
        } else {
            None
        };
        let Some(log) = factor else {
            return Ok(None);
        };
        let mu = exponential(ev, &log)?;
        (
            ev.eval(&(mu.clone() * m.clone()))?,
            ev.eval(&(mu * n.clone()))?,
        )
    };
    let Some(p) = integral(ev, &m, x)? else {
        return Ok(None);
    };
    let py = d(ev, &p, y)?;
    let h_prime = together(ev, n - py)?;
    if h_prime.contains(x) {
        return Ok(None);
    }
    let Some(h) = integral(ev, &h_prime, y)? else {
        return Ok(None);
    };
    Ok(roots(ev, p + h, c.clone(), y)?.map(branches))
}
