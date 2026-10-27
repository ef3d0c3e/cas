//! Particular solutions of linear equations with a right side (CAS-DOC-16
//! §5, CAS-262), `Σ aₖ y⁽ᵏ⁾ = f` over the basis of the homogeneous one:
//!
//! - **Undetermined coefficients**: `f`, its products of sines and cosines
//!   reduced (`TrigReduce`), expanded into terms `c x^j e^(α x + γ)` times
//!   `cos(β x)`, `sin(β x)` or neither,
//!   grouped by `(α, β)`. A group's trial is `x^s e^(α x) (Q₁(x) cos(β x) +
//!   Q₂(x) sin(β x))`, `s` the multiplicity of `α + β i` as a
//!   characteristic root and `Q₁`, `Q₂` of the group's degree with unknown
//!   coefficients, which the linear system of `L[trial] = f` fixes, each
//!   power of `x` times the cosine and the sine apart.
//! - **Variation of parameters** for the other terms: `Σᵢ uᵢ yᵢ`, `uᵢ′` the
//!   solution of the Wronskian system `Σᵢ yᵢ⁽ᵏ⁾ uᵢ′ = 0` for `k < n − 1`
//!   and `f/aₙ` for `k = n − 1`, simplified, then integrated; an integral
//!   that stays stays in the answer.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equation: Σ aₖ y⁽ᵏ⁾ = f, trial coefficients"
)]

use cas_num::Complex;

use super::super::solve::linear_system;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// What multiplies `c x^j e^(α x)` in a term of the right side.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Trig {
    None,
    Cos,
    Sin,
}

/// A term `c x^j e^(α x) trig(β x)` of the right side.
struct Term {
    expr: Expr,
    power: u32,
    alpha: Expr,
    beta: Expr,
}

/// The `k`-th derivative of `e` in `x`.
fn derivative(ev: &mut Eval<'_>, e: &Expr, x: &Expr, k: usize) -> Result<Expr, EvalError> {
    if k == 0 {
        return Ok(e.clone());
    }
    let order = Expr::from(i64::try_from(k).unwrap_or(i64::MAX));
    let variable = Expr::call(Symbol::LIST, [x.clone(), order]);
    ev.eval(&Expr::call(Symbol::D, [e.clone(), variable]))
}

/// `Σ aₖ e⁽ᵏ⁾`.
fn operator(ev: &mut Eval<'_>, a: &[Expr], e: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let mut terms = Vec::with_capacity(a.len());
    for (k, ak) in a.iter().enumerate() {
        if *ak != Expr::from(0) {
            terms.push(ak.clone() * derivative(ev, e, x, k)?);
        }
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms))
}

/// A particular solution of `Σ aₖ y⁽ᵏ⁾ = f`, `basis` and `roots` (with
/// their multiplicities) those of the homogeneous equation; `None` when a
/// system has no solution.
pub(super) fn particular(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    f: &Expr,
    x: &Expr,
    basis: &[Expr],
    roots: &[(Expr, usize)],
) -> Result<Option<Expr>, EvalError> {
    let reduced = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [Expr::call(Symbol::TRIG_REDUCE, [f.clone()])],
    ))?;
    let terms = if reduced.head() == Some(&Expr::from(Symbol::PLUS)) {
        reduced.args().to_vec()
    } else {
        vec![reduced]
    };
    let mut groups: Vec<(Expr, Expr, Vec<Term>)> = Vec::new();
    let mut rest = Vec::new();
    for t in terms {
        match classify(ev, &t, x)? {
            Some(term) => match groups
                .iter_mut()
                .find(|(a, b, _)| *a == term.alpha && *b == term.beta)
            {
                Some((_, _, ts)) => ts.push(term),
                None => groups.push((term.alpha.clone(), term.beta.clone(), vec![term])),
            },
            None => rest.push(t),
        }
    }
    let mut solution = Vec::new();
    for (alpha, beta, group) in groups {
        match undetermined(ev, e, a, x, roots, &alpha, &beta, &group)? {
            Some(y) => solution.push(y),
            None => rest.extend(group.into_iter().map(|t| t.expr)),
        }
    }
    if !rest.is_empty() {
        let g = ev.eval(&Expr::call(Symbol::PLUS, rest))?;
        let Some(y) = variation(ev, e, a, x, basis, &g)? else {
            return Ok(None);
        };
        solution.push(y);
    }
    ev.eval(&Expr::call(
        Symbol::EXPAND,
        [Expr::call(Symbol::PLUS, solution)],
    ))
    .map(Some)
}

/// `β` for `u = β x + γ`, `β` and `γ` free of `x`, and `γ = 0` unless
/// `offset`: an exponential's constant factor `e^γ` is the term's, a sine's
/// or a cosine's phase is not.
fn slope(ev: &mut Eval<'_>, u: &Expr, x: &Expr, offset: bool) -> Result<Option<Expr>, EvalError> {
    let beta = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
    if beta.contains(x) {
        return Ok(None);
    }
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - beta.clone() * x.clone()],
    ))?;
    let fits = if offset {
        !gap.contains(x)
    } else {
        gap == Expr::from(0)
    };
    Ok(fits.then_some(beta))
}

/// A term of the right side as `c x^j e^(α x) trig(β x)`; `None` when it
/// is not one.
fn classify(ev: &mut Eval<'_>, t: &Expr, x: &Expr) -> Result<Option<Term>, EvalError> {
    let factors = if t.head() == Some(&Expr::from(Symbol::TIMES)) {
        t.args().to_vec()
    } else {
        vec![t.clone()]
    };
    let (mut power, mut alphas) = (0u32, Vec::new());
    let mut trig = (Trig::None, Expr::from(0));
    for factor in factors {
        if !factor.contains(x) {
            continue;
        }
        if factor == *x {
            power += 1;
            continue;
        }
        let (head, args) = (factor.head().and_then(Expr::as_symbol), factor.args());
        match (head, args) {
            (Some(Symbol::POWER), [b, k]) if b == x => {
                let Some(k) = k.as_i64().and_then(|k| u32::try_from(k).ok()) else {
                    return Ok(None);
                };
                power += k;
            }
            (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => {
                let Some(alpha) = slope(ev, u, x, true)? else {
                    return Ok(None);
                };
                alphas.push(alpha);
            }
            (Some(Symbol::EXP), [u]) => {
                let Some(alpha) = slope(ev, u, x, true)? else {
                    return Ok(None);
                };
                alphas.push(alpha);
            }
            (Some(s @ (Symbol::COS | Symbol::SIN)), [u]) if trig.0 == Trig::None => {
                let Some(beta) = slope(ev, u, x, false)? else {
                    return Ok(None);
                };
                let kind = if s == Symbol::COS {
                    Trig::Cos
                } else {
                    Trig::Sin
                };
                trig = (kind, beta);
            }
            _ => return Ok(None),
        }
    }
    let alpha = ev.eval(&Expr::call(Symbol::PLUS, alphas))?;
    Ok(Some(Term {
        expr: t.clone(),
        power,
        alpha,
        beta: trig.1,
    }))
}

/// The particular solution of the group `(α, β)` of terms by undetermined
/// coefficients; `None` when its system has no solution.
#[allow(
    clippy::too_many_arguments,
    reason = "the equation, its roots and the group"
)]
fn undetermined(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    x: &Expr,
    roots: &[(Expr, usize)],
    alpha: &Expr,
    beta: &Expr,
    group: &[Term],
) -> Result<Option<Expr>, EvalError> {
    let degree = group.iter().map(|t| t.power).max().unwrap_or(0);
    let d = usize::try_from(degree).unwrap_or(usize::MAX);
    let trig = *beta != Expr::from(0);
    // The resonance: α + β i's multiplicity as a root.
    let root = ev.eval(&(alpha.clone() + Expr::from(Complex::I) * beta.clone()))?;
    let s = roots
        .iter()
        .find(|(r, _)| *r == root)
        .map_or(0, |(_, m)| *m);
    let count = (d + 1) * if trig { 2 } else { 1 };
    let unknowns: Vec<Expr> = (0..count)
        .map(|i| Expr::from(Symbol::new(&format!("DSolve$u{i}"))))
        .collect();
    let polynomial = |offset: usize| {
        Expr::call(
            Symbol::PLUS,
            (0..=d).map(|j| unknowns[offset + j].clone() * x.pow(i64::try_from(j).unwrap_or(0))),
        )
    };
    let angle = beta.clone() * x.clone();
    let (cos, sin) = (
        ev.eval(&Expr::call(Symbol::COS, [angle.clone()]))?,
        ev.eval(&Expr::call(Symbol::SIN, [angle]))?,
    );
    let shape = if trig {
        polynomial(0) * cos.clone() + polynomial(d + 1) * sin.clone()
    } else {
        polynomial(0)
    };
    let s_power = x.pow(i64::try_from(s).unwrap_or(0));
    let exponential = Expr::from(Symbol::E).pow(alpha.clone() * x.clone());
    let trial = ev.eval(&(s_power * exponential.clone() * shape))?;
    let lhs = operator(ev, a, &trial, x)?;
    let rhs = Expr::call(Symbol::PLUS, group.iter().map(|t| t.expr.clone()));
    let scaled = (lhs - rhs) * Expr::from(Symbol::E).pow(-(alpha.clone() * x.clone()));
    let residual = ev.eval(&Expr::call(Symbol::EXPAND, [scaled]))?;
    // The cosine and the sine as symbols, then each power of x apart.
    let (c, s_symbol) = (
        Expr::from(Symbol::new("DSolve$cos")),
        Expr::from(Symbol::new("DSolve$sin")),
    );
    let residual = if trig {
        residual.substitute(&[(cos, c.clone()), (sin, s_symbol.clone())])
    } else {
        residual
    };
    let residual = ev.eval(&Expr::call(Symbol::EXPAND, [residual]))?;
    let top = ev.eval(&Expr::call(Symbol::EXPONENT, [residual.clone(), x.clone()]))?;
    let Some(top) = top.as_i64().and_then(|t| usize::try_from(t).ok()) else {
        return Ok(None);
    };
    let parts: Vec<Expr> = if trig {
        let mut v = Vec::with_capacity(2);
        for symbol in [&c, &s_symbol] {
            v.push(ev.eval(&Expr::call(
                Symbol::COEFFICIENT,
                [residual.clone(), symbol.clone(), Expr::from(1)],
            ))?);
        }
        let rest = ev.eval(&Expr::call(
            Symbol::EXPAND,
            [residual.substitute(&[
                (c.clone(), Expr::from(0)),
                (s_symbol.clone(), Expr::from(0)),
            ])],
        ))?;
        if rest != Expr::from(0) {
            return Ok(None);
        }
        v
    } else {
        vec![residual]
    };
    let mut equations = Vec::new();
    for part in &parts {
        if part.contains(&c) || part.contains(&s_symbol) {
            return Ok(None);
        }
        for j in 0..=top {
            let coefficient = ev.eval(&Expr::call(
                Symbol::COEFFICIENT,
                [
                    part.clone(),
                    x.clone(),
                    Expr::from(i64::try_from(j).unwrap_or(0)),
                ],
            ))?;
            if coefficient.contains(x) {
                return Ok(None);
            }
            equations.push(Expr::call(Symbol::EQUATION, [coefficient, Expr::from(0)]));
        }
    }
    let Some(values) = solved(ev, e, &equations, &unknowns)? else {
        return Ok(None);
    };
    let rules: Vec<(Expr, Expr)> = unknowns.into_iter().zip(values).collect();
    ev.eval(&trial.substitute(&rules)).map(Some)
}

/// The values of `unknowns` from the linear system `equations`, as
/// `LinSolve` finds the one solution; `None` when there is none.
fn solved(
    ev: &mut Eval<'_>,
    e: &Expr,
    equations: &[Expr],
    unknowns: &[Expr],
) -> Result<Option<Vec<Expr>>, EvalError> {
    let Some(found) = linear_system(ev, e, equations, unknowns)? else {
        return Ok(None);
    };
    let [solution] = found.args() else {
        return Ok(None);
    };
    let values = solution.args().to_vec();
    let free = values
        .iter()
        .all(|v| unknowns.iter().all(|u| !v.contains(u)));
    Ok((values.len() == unknowns.len() && free).then_some(values))
}

/// The particular solution for `g` by variation of parameters over
/// `basis`.
fn variation(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    x: &Expr,
    basis: &[Expr],
    g: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let n = basis.len();
    let Some(an) = a.get(n) else {
        return Ok(None);
    };
    let ws: Vec<Expr> = (0..n)
        .map(|i| Expr::from(Symbol::new(&format!("DSolve$w{i}"))))
        .collect();
    let mut equations = Vec::with_capacity(n);
    for k in 0..n {
        let mut row = Vec::with_capacity(n);
        for (b, w) in basis.iter().zip(&ws) {
            row.push(derivative(ev, b, x, k)? * w.clone());
        }
        let rhs = if k + 1 == n {
            g.clone() / an.clone()
        } else {
            Expr::from(0)
        };
        equations.push(Expr::call(
            Symbol::EQUATION,
            [Expr::call(Symbol::PLUS, row), rhs],
        ));
    }
    let Some(slopes) = solved(ev, e, &equations, &ws)? else {
        return Ok(None);
    };
    let mut terms = Vec::with_capacity(n);
    for (b, w) in basis.iter().zip(slopes) {
        let simplified = ev.eval(&Expr::call(Symbol::SIMPLIFY, [w]))?;
        let u = ev.eval(&Expr::call(Symbol::INTEGRATE, [simplified, x.clone()]))?;
        terms.push(u * b.clone());
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}
