//! Recurrences (CAS-DOC-24 §4, CAS-310): `RSolve`, Xcas's `rsolve`.
//!
//! - **Forms**: `RSolve[eq, u[n]]`, with conditions after it, or a list of
//!   equations in a list of unknowns, `RSolve[{eq₁, …}, {u[n], v[n], …},
//!   conditions…]`. The unknowns occur at shifts `n + k` with integer
//!   `k`; a condition is an equation in their values at numbers, `u[0] =
//!   1`.
//! - **Constant coefficients**: each equation is moved up by its least
//!   shift (`u(n) = 2u(n − 1) + 1` as `u(n + 1) = 2u(n) + 1`). Then it must
//!   be linear in the `uⱼ(n + k)` with coefficients free of `n`, its
//!   forcing `g(n)` what is left, with a nonsingular matrix of the highest
//!   shifts. The Z-transform makes `uⱼ(n + k)` into `zᵏ (Uⱼ − Σ_{i<k} uⱼ(i)
//!   z⁻ⁱ)`, the initial values `uⱼ(0), …` naming the constants `c_0, c_1,
//!   …`, and `g` into `ZTransform`'s. `LinSolve` solves for the `Uⱼ`, and
//!   each is split by the constants and inverted (`InverseZTransform`), as
//!   `DSolve` solves systems (CAS-DOC-20).
//! - **First order with coefficients in `n`**, one unknown: `u(n + 1) =
//!   a(n) u(n) + b(n)` is `P(n) (c_0 + Σ_{j<n} b(j)/P(j + 1))`, with `P(n) =
//!   Π_{j<n} a(j)`, by `Product` and `Sum`: `u(n + 1) = (n + 1) u(n)` is `c_0
//!   n!`. A product or sum without a closed form leaves the call.
//! - **Conditions** fix the constants (`LinSolve`), and none consistent
//!   gives `{}`.
//! - The answer is `{u(n)}`, or `{{u(n), v(n), …}}` for a list of unknowns,
//!   expanded, as giac's `[2*2^n-1]` and `[[n,1]]`. What these do not solve
//!   stays, with a warning.
//!
//! Each call that solves takes one step, `"RSolve"`.

#![allow(
    clippy::many_single_char_names,
    reason = "the recurrence's names: u(n + k), a(n), b(n), z and U"
)]

use super::dsolve::constants;
use super::solve::linear_system;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// A recurrence: its equations as residuals, the heads of its unknowns, the
/// variable, the conditions, and whether the unknowns came as a list.
struct Problem {
    residuals: Vec<Expr>,
    heads: Vec<Expr>,
    n: Expr,
    conditions: Vec<Expr>,
    listed: bool,
}

/// `l − r` of an equation, or the expression itself.
fn residual(eq: &Expr) -> Expr {
    match (eq.head().and_then(Expr::as_symbol), eq.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
        _ => eq.clone(),
    }
}

/// The recurrence of `RSolve[eqs, unknowns, conditions…]`.
fn problem(e: &Expr) -> Option<Problem> {
    let [eqs, unknowns, conditions @ ..] = e.args() else {
        return None;
    };
    let listed = is_list(unknowns);
    let unknowns: Vec<Expr> = if listed {
        unknowns.args().to_vec()
    } else {
        vec![unknowns.clone()]
    };
    let mut heads = Vec::with_capacity(unknowns.len());
    let mut variable: Option<Expr> = None;
    for u in &unknowns {
        let View::Apply { head, args: [n] } = u.view() else {
            return None;
        };
        if head.as_symbol().is_none_or(Symbol::is_builtin)
            || n.as_symbol().is_none_or(Symbol::is_builtin)
            || variable.as_ref().is_some_and(|v| v != n)
        {
            return None;
        }
        variable = Some(n.clone());
        heads.push(head.clone());
    }
    let equations: Vec<Expr> = if is_list(eqs) {
        eqs.args().to_vec()
    } else {
        vec![eqs.clone()]
    };
    Some(Problem {
        residuals: equations.iter().map(residual).collect(),
        heads,
        n: variable?,
        conditions: conditions.to_vec(),
        listed,
    })
}

/// The shift `k` of an argument `n + k`, an integer.
fn shift(ev: &mut Eval<'_>, a: &Expr, n: &Expr) -> Result<Option<i64>, EvalError> {
    let k = ev.eval(&(a - n))?;
    Ok(if k.contains(n) { None } else { k.as_i64() })
}

/// The shifts of the unknowns `heads` of `n` in `r`: `(j, k)` for each
/// `uⱼ(n + k)`; `None` when an unknown occurs otherwise.
fn shifts(
    ev: &mut Eval<'_>,
    r: &Expr,
    heads: &[Expr],
    n: &Expr,
) -> Result<Option<Vec<(usize, i64)>>, EvalError> {
    let mut found = Vec::new();
    for part in r.distinct() {
        let View::Apply { head, args } = part.view() else {
            continue;
        };
        let Some(j) = heads.iter().position(|h| h == head) else {
            continue;
        };
        let [a] = args else {
            return Ok(None);
        };
        let Some(k) = shift(ev, a, n)? else {
            return Ok(None);
        };
        found.push((j, k));
    }
    Ok(Some(found))
}

/// The symbol of `uⱼ(n + k)`.
fn shifted_symbol(j: usize, k: i64) -> Expr {
    Expr::from(Symbol::new(&format!("RSolve$u{j}s{k}")))
}

/// `RSolve[…]`.
pub(super) fn rsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(mut p) = problem(e) else {
        return Ok(None);
    };
    // Each equation moved up by its least shift.
    let mut orders = vec![0i64; p.heads.len()];
    for r in &mut p.residuals {
        let Some(found) = shifts(ev, r, &p.heads, &p.n)? else {
            return Ok(None);
        };
        let Some(least) = found.iter().map(|&(_, k)| k).min() else {
            return Ok(None);
        };
        if least != 0 {
            *r = ev.eval(&r.substitute(&[(p.n.clone(), &p.n - least)]))?;
        }
        for (j, k) in found {
            orders[j] = orders[j].max(k - least);
        }
    }
    let components = if let Some(c) = constant_coefficients(ev, e, &p, &orders)? {
        c
    } else if let Some(c) = first_order(ev, e, &p)? {
        c
    } else {
        ev.warn(Warning::Native {
            message: "no method solves this recurrence".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    let result = if p.conditions.is_empty() {
        Some(answer(&p, components))
    } else {
        fix(ev, e, &p, &components)?
    };
    if let Some(r) = &result {
        ev.step("RSolve", e, r);
    }
    Ok(result)
}

/// `{u}` or `{{u, v, …}}`.
fn answer(p: &Problem, components: Vec<Expr>) -> Expr {
    if p.listed {
        Expr::call(Symbol::LIST, [Expr::call(Symbol::LIST, components)])
    } else {
        Expr::call(Symbol::LIST, components)
    }
}

/// The residuals with each `uⱼ(n + k)` as its symbol.
fn symbolic(p: &Problem, ev: &mut Eval<'_>) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut out = Vec::with_capacity(p.residuals.len());
    for r in &p.residuals {
        let mut failed = false;
        let mut shifted = Vec::new();
        for part in r.distinct() {
            let View::Apply { head, args: [a] } = part.view() else {
                continue;
            };
            let Some(j) = p.heads.iter().position(|h| h == head) else {
                continue;
            };
            match shift(ev, a, &p.n)? {
                Some(k) => shifted.push((part.clone(), shifted_symbol(j, k))),
                None => failed = true,
            }
        }
        if failed {
            return Ok(None);
        }
        out.push(ev.eval(&r.substitute(&shifted))?);
    }
    Ok(Some(out))
}

/// A system with constant coefficients: equation `i` is `Σⱼₖ a[i][j][k]
/// uⱼ(n + k) = forcing[i]`.
struct Linear {
    a: Vec<Vec<Vec<Expr>>>,
    forcing: Vec<Expr>,
}

impl Linear {
    /// The system of `p`, its unknowns of the highest shifts `orders`;
    /// `None` when it is not linear with constant coefficients and a
    /// nonsingular matrix of the highest shifts.
    fn of(ev: &mut Eval<'_>, p: &Problem, orders: &[i64]) -> Result<Option<Linear>, EvalError> {
        let m = p.heads.len();
        let Some(residuals) = symbolic(p, ev)? else {
            return Ok(None);
        };
        let symbols: Vec<Expr> = (0..m)
            .flat_map(|j| (0..=orders[j]).map(move |k| shifted_symbol(j, k)))
            .collect();
        let mut a: Vec<Vec<Vec<Expr>>> = Vec::with_capacity(m);
        let mut forcing = Vec::with_capacity(m);
        for r in &residuals {
            let mut row = Vec::with_capacity(m);
            for (j, &order) in orders.iter().enumerate() {
                let mut coefficients = Vec::new();
                for k in 0..=order {
                    let c = ev.eval(&Expr::call(Symbol::D, [r.clone(), shifted_symbol(j, k)]))?;
                    if c.contains(&p.n) || symbols.iter().any(|s| c.contains(s)) {
                        return Ok(None);
                    }
                    coefficients.push(c);
                }
                row.push(coefficients);
            }
            a.push(row);
            let zero: Vec<(Expr, Expr)> =
                symbols.iter().map(|s| (s.clone(), Expr::from(0))).collect();
            forcing.push(ev.eval(&-r.substitute(&zero))?);
        }
        let leading: Vec<Expr> = a
            .iter()
            .map(|row| {
                Expr::call(
                    Symbol::LIST,
                    row.iter()
                        .zip(orders)
                        .map(|(c, &k)| c[usize::try_from(k).unwrap_or(0)].clone()),
                )
            })
            .collect();
        let det = ev.eval(&Expr::call(
            Symbol::DET,
            [Expr::call(Symbol::LIST, leading)],
        ))?;
        if det == Expr::from(0) || det.head() == Some(&Expr::from(Symbol::DET)) {
            return Ok(None);
        }
        Ok(Some(Linear { a, forcing }))
    }
}

/// The general solution of a system with constant coefficients by
/// Z-transforms, the initial values the constants; `None` when it is not
/// one, or a transform stays.
fn constant_coefficients(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Problem,
    orders: &[i64],
) -> Result<Option<Vec<Expr>>, EvalError> {
    if p.residuals.len() != p.heads.len() || orders.iter().all(|&k| k == 0) {
        return Ok(None);
    }
    let Some(system) = Linear::of(ev, p, orders)? else {
        return Ok(None);
    };
    let count = usize::try_from(orders.iter().sum::<i64>()).unwrap_or(0);
    let constants = constants(e, count);
    let constants = &constants[..count];
    let z = Expr::from(Symbol::new("RSolve$z"));
    let Some(transforms) = z_domain(ev, e, p, &system, orders, constants, &z)? else {
        return Ok(None);
    };
    let mut components = Vec::with_capacity(transforms.len());
    for v in &transforms {
        // Uⱼ = H + Σ c ∂Uⱼ/∂c, each part inverted.
        let zero: Vec<(Expr, Expr)> = constants
            .iter()
            .map(|c| (c.clone(), Expr::from(0)))
            .collect();
        let base = ev.eval(&v.substitute(&zero))?;
        let Some(h) = inverse(ev, &base, &z, &p.n)? else {
            return Ok(None);
        };
        let mut terms = vec![h];
        for c in constants {
            let part = ev.eval(&Expr::call(Symbol::D, [v.clone(), c.clone()]))?;
            let Some(g) = inverse(ev, &part, &z, &p.n)? else {
                return Ok(None);
            };
            terms.push(c * g);
        }
        components.push(ev.eval(&Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::PLUS, terms)],
        ))?);
    }
    Ok(Some(components))
}

/// The transforms `Uⱼ` of the unknowns: the system transformed, `uⱼ(n + k)`
/// as `zᵏ (Uⱼ − Σ_{i<k} uⱼ(i) z⁻ⁱ)` with the initial values `constants`,
/// solved by `LinSolve`; `None` when a forcing's transform stays.
fn z_domain(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Problem,
    system: &Linear,
    orders: &[i64],
    constants: &[Expr],
    z: &Expr,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let m = p.heads.len();
    // The constant of uⱼ(i) is initial[j][i].
    let mut rest = constants.iter();
    let initial: Vec<Vec<Expr>> = orders
        .iter()
        .map(|&k| {
            rest.by_ref()
                .take(usize::try_from(k).unwrap_or(0))
                .cloned()
                .collect()
        })
        .collect();
    let transforms: Vec<Expr> = (0..m)
        .map(|j| Expr::from(Symbol::new(&format!("RSolve$U{j}"))))
        .collect();
    let mut equations = Vec::with_capacity(m);
    for (row, g) in system.a.iter().zip(&system.forcing) {
        let mut left = Vec::new();
        let mut right = Vec::new();
        if *g != Expr::from(0) {
            let image = ev.eval(&Expr::call(
                Symbol::ZTRANSFORM,
                [g.clone(), p.n.clone(), z.clone()],
            ))?;
            if image.contains(&Expr::from(Symbol::ZTRANSFORM)) {
                return Ok(None);
            }
            right.push(image);
        }
        for (j, coefficients) in row.iter().enumerate() {
            for (k, c) in coefficients.iter().enumerate() {
                if *c == Expr::from(0) {
                    continue;
                }
                let k = i64::try_from(k).unwrap_or(0);
                left.push(c * z.pow(k) * &transforms[j]);
                for (i, value) in initial[j]
                    .iter()
                    .enumerate()
                    .take(usize::try_from(k).unwrap_or(0))
                {
                    let i = i64::try_from(i).unwrap_or(0);
                    right.push(c * z.pow(k - i) * value);
                }
            }
        }
        equations.push(Expr::call(
            Symbol::EQUATION,
            [
                Expr::call(Symbol::PLUS, left),
                Expr::call(Symbol::PLUS, right),
            ],
        ));
    }
    let Some(found) = linear_system(ev, e, &equations, &transforms)? else {
        return Ok(None);
    };
    let [values] = found.args() else {
        return Ok(None);
    };
    Ok((values.args().len() == m).then(|| values.args().to_vec()))
}

/// The inverse transform of `image`, expanded first; `None` when it stays.
fn inverse(ev: &mut Eval<'_>, image: &Expr, z: &Expr, n: &Expr) -> Result<Option<Expr>, EvalError> {
    if *image == Expr::from(0) {
        return Ok(Some(Expr::from(0)));
    }
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [image.clone()]))?;
    let f = ev.eval(&Expr::call(
        Symbol::INVERSE_ZTRANSFORM,
        [together, z.clone(), n.clone()],
    ))?;
    Ok((!f.contains(&Expr::from(Symbol::INVERSE_ZTRANSFORM))).then_some(f))
}

/// One unknown of first order with coefficients in `n`: `u(n + 1) = a(n)
/// u(n) + b(n)`, `P(n) (c_0 + Σ_{j<n} b(j)/P(j + 1))`; `None` when it is
/// not one, or a product or sum stays.
fn first_order(ev: &mut Eval<'_>, e: &Expr, p: &Problem) -> Result<Option<Vec<Expr>>, EvalError> {
    let ([_], [_]) = (p.heads.as_slice(), p.residuals.as_slice()) else {
        return Ok(None);
    };
    let Some(residuals) = symbolic(p, ev)? else {
        return Ok(None);
    };
    let r = &residuals[0];
    let (now, next) = (shifted_symbol(0, 0), shifted_symbol(0, 1));
    let alpha = ev.eval(&Expr::call(Symbol::D, [r.clone(), next.clone()]))?;
    let beta = ev.eval(&Expr::call(Symbol::D, [r.clone(), now.clone()]))?;
    let zero = [(now.clone(), Expr::from(0)), (next.clone(), Expr::from(0))];
    let rest = ev.eval(&r.substitute(&zero))?;
    let symbols_left = |x: &Expr| x.contains(&now) || x.contains(&next);
    // Linear in u(n), u(n + 1), and nothing else of u.
    if alpha == Expr::from(0)
        || [&alpha, &beta].iter().any(|c| symbols_left(c))
        || r.distinct().any(|part| part.head() == Some(&p.heads[0]))
    {
        return Ok(None);
    }
    let a = ev.eval(&(-(beta / &alpha)))?;
    let b = ev.eval(&(-(rest / &alpha)))?;
    let j = Expr::from(Symbol::new("RSolve$j"));
    let at = |x: &Expr, k: &Expr| x.substitute(&[(p.n.clone(), k.clone())]);
    // `Product[f, j, 0, top]`, as `Sum`'s forms: giac's `product(f, j, 0, top)`.
    let product = |ev: &mut Eval<'_>, top: Expr| {
        ev.eval(&Expr::call(
            Symbol::PRODUCT,
            [at(&a, &j), j.clone(), Expr::from(0), top],
        ))
    };
    let big_p = product(ev, &p.n - 1)?;
    if big_p.contains(&Expr::from(Symbol::PRODUCT)) {
        return Ok(None);
    }
    let constant = constants(e, 1).swap_remove(0);
    let sum = if b == Expr::from(0) {
        Expr::from(0)
    } else {
        let next_p = big_p.substitute(&[(p.n.clone(), &j + 1)]);
        let s = ev.eval(&Expr::call(
            Symbol::SUM,
            [at(&b, &j) / next_p, j.clone(), Expr::from(0), &p.n - 1],
        ))?;
        if s.contains(&Expr::from(Symbol::SUM)) {
            return Ok(None);
        }
        s
    };
    let u = ev.eval(&Expr::call(Symbol::EXPAND, [big_p * (constant + sum)]))?;
    Ok(Some(vec![u]))
}

/// The general solution's `components` with their constants fixed by the
/// conditions: `{}` when none meets them; `None` when a condition is not
/// an equation in values of the unknowns.
fn fix(ev: &mut Eval<'_>, e: &Expr, p: &Problem, components: &[Expr]) -> Native {
    let mut equations = Vec::with_capacity(p.conditions.len());
    for condition in &p.conditions {
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) =
            (condition.head().and_then(Expr::as_symbol), condition.args())
        else {
            return Ok(None);
        };
        let at = (l.clone() - r.clone()).replace(|part| {
            let View::Apply { head, args: [a] } = part.view() else {
                return None;
            };
            let j = p.heads.iter().position(|h| h == head)?;
            Some(components[j].substitute(&[(p.n.clone(), a.clone())]))
        });
        if at.contains(&p.n) || p.heads.iter().any(|h| at.contains(h)) {
            return Ok(None);
        }
        let value = ev.eval(&at)?;
        equations.push(Expr::call(Symbol::EQUATION, [value, Expr::from(0)]));
    }
    let candidates = constants(e, components.len().max(8) * 4);
    let used: Vec<Expr> = candidates
        .into_iter()
        .filter(|c| components.iter().any(|v| v.contains(c)))
        .collect();
    if used.is_empty() {
        return Ok(None);
    }
    let Some(found) = linear_system(ev, e, &equations, &used)? else {
        return Ok(None);
    };
    match found.args() {
        [] => Ok(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new()))),
        [values] if values.args().len() == used.len() => {
            let rules: Vec<(Expr, Expr)> = used
                .iter()
                .cloned()
                .zip(values.args().iter().cloned())
                .collect();
            let mut fixed = Vec::with_capacity(components.len());
            for v in components {
                fixed.push(ev.eval(&Expr::call(Symbol::EXPAND, [v.substitute(&rules)]))?);
            }
            Ok(Some(answer(p, fixed)))
        }
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests;
