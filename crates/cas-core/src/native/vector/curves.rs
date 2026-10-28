//! Curves (CAS-DOC-37 §4, CAS-410): `Curvature`, `Frenet` and `Evolute`,
//! giac's `curvature`, `frenet` and `evolute`, of a curve `C(t)`: a list of
//! two or three expressions in `t`, or an expression `f`, the graph `{t,
//! f}`.
//!
//! - `Curvature[C, t]`: in the plane `κ = (x′y″ − y′x″)/(x′² +
//!   y′²)^(3/2)`, signed as giac's; in space `|C′ × C″|/|C′|³`.
//! - `Frenet[C, t]`, as giac's lists: in the plane `{κ, M, T, N}`, the unit
//!   tangent `T = C′/|C′|`, the normal `N`, which is `T` turned by `+π/2`,
//!   and the centre of curvature `M = C + N/κ`; in space `{κ, M, τ, T, N,
//!   B}`, with the binormal `B = (C′ × C″)/|C′ × C″|`, `N = B × T` and the
//!   torsion `τ = ((C′ × C″)·C‴)/|C′ × C″|²`.
//! - `Evolute[C, t]`: in the plane, the centres of curvature, `C + (x′² +
//!   y′²)/(x′y″ − y′x″) (−y′, x′)`, made `Together`.
//! - A third argument `t0`: the answer at `t = t0`, evaluated.
//!
//! The curvature and the torsion are made `Simplify`: a circle's is a
//! number. `t` is a symbol; otherwise the call stays.

use super::entries;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn half() -> Expr {
    Expr::from(cas_num::Rational::new(
        cas_num::Integer::ONE,
        cas_num::Integer::from(2),
    ))
}

/// `Σ a_k b_k`, built.
fn dot(a: &[Expr], b: &[Expr]) -> Expr {
    Expr::call(
        Symbol::PLUS,
        a.iter().zip(b).map(|(x, y)| x.clone() * y.clone()),
    )
}

/// `a × b` in space, built.
fn cross(a: &[Expr], b: &[Expr]) -> Vec<Expr> {
    let (a1, a2, a3) = (&a[0], &a[1], &a[2]);
    let (b1, b2, b3) = (&b[0], &b[1], &b[2]);
    vec![
        a2.clone() * b3.clone() - a3.clone() * b2.clone(),
        a3.clone() * b1.clone() - a1.clone() * b3.clone(),
        a1.clone() * b2.clone() - a2.clone() * b1.clone(),
    ]
}

/// `v/s`, entry by entry, built.
fn scaled(v: &[Expr], s: &Expr) -> Vec<Expr> {
    v.iter().map(|x| x.clone() / s.clone()).collect()
}

/// The curve's points and their first three derivatives, evaluated; `None`
/// for a curve that is not one, or a `t` that is no symbol.
pub(super) fn curve(
    ev: &mut Eval<'_>,
    c: &Expr,
    t: &Expr,
) -> Result<Option<[Vec<Expr>; 4]>, EvalError> {
    if t.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    let points = match entries(c) {
        Some(items) if (2..=3).contains(&items.len()) => items.to_vec(),
        Some(_) => return Ok(None),
        None => vec![t.clone(), c.clone()],
    };
    let mut derivatives = vec![points];
    for k in 0..3 {
        let next = derivatives[k]
            .iter()
            .map(|x| ev.eval(&Expr::call(Symbol::D, [x.clone(), t.clone()])))
            .collect::<Result<Vec<_>, _>>()?;
        derivatives.push(next);
    }
    let [r, d1, d2, d3]: [Vec<Expr>; 4] = derivatives.try_into().expect("four");
    Ok(Some([r, d1, d2, d3]))
}

/// `e` at `t = t0`, or as it is without `t0`, evaluated.
fn at(ev: &mut Eval<'_>, e: &Expr, t: &Expr, t0: Option<&Expr>) -> Result<Expr, EvalError> {
    match t0 {
        Some(t0) => ev.eval(&e.substitute(&[(t.clone(), t0.clone())])),
        None => ev.eval(e),
    }
}

/// The curvature from the derivatives, built: signed in the plane.
fn kappa(d1: &[Expr], d2: &[Expr]) -> Expr {
    let speed_squared = dot(d1, d1);
    let three_halves = Expr::from(3) * half();
    if d1.len() == 2 {
        (d1[0].clone() * d2[1].clone() - d1[1].clone() * d2[0].clone())
            / speed_squared.pow(three_halves)
    } else {
        let c = cross(d1, d2);
        dot(&c, &c).pow(half()) / speed_squared.pow(three_halves)
    }
}

/// The parts of a call: the curve, `t` and `t0`.
fn arguments(e: &Expr) -> Option<(&Expr, &Expr, Option<&Expr>)> {
    match e.args() {
        [c, t] => Some((c, t, None)),
        [c, t, t0] => Some((c, t, Some(t0))),
        _ => None,
    }
}

/// `Curvature[C, t]` and `Curvature[C, t, t0]`.
pub(in crate::native) fn curvature(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((c, t, t0)) = arguments(e) else {
        return Ok(None);
    };
    let Some([_, d1, d2, _]) = curve(ev, c, t)? else {
        return Ok(None);
    };
    let k = Expr::call(Symbol::SIMPLIFY, [kappa(&d1, &d2)]);
    let result = at(ev, &k, t, t0)?;
    ev.step("curvature", e, &result);
    Ok(Some(result))
}

/// `Frenet[C, t]` and `Frenet[C, t, t0]`.
pub(in crate::native) fn frenet(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((c, t, t0)) = arguments(e) else {
        return Ok(None);
    };
    let Some([r, d1, d2, d3]) = curve(ev, c, t)? else {
        return Ok(None);
    };
    let k = Expr::call(Symbol::SIMPLIFY, [kappa(&d1, &d2)]);
    let speed = dot(&d1, &d1).pow(half());
    let tangent = scaled(&d1, &speed);
    let (normal, rest) = if r.len() == 2 {
        (vec![-tangent[1].clone(), tangent[0].clone()], Vec::new())
    } else {
        let c = cross(&d1, &d2);
        let size = dot(&c, &c);
        let binormal = scaled(&c, &size.pow(half()));
        let normal = cross(&binormal, &tangent);
        let torsion = Expr::call(Symbol::SIMPLIFY, [dot(&c, &d3) / size]);
        (normal, vec![torsion, list(binormal)])
    };
    let centre: Vec<Expr> = r
        .iter()
        .zip(&normal)
        .map(|(x, n)| x.clone() + n.clone() / k.clone())
        .collect();
    let items = if rest.is_empty() {
        vec![k, list(centre), list(tangent), list(normal)]
    } else {
        let [torsion, binormal]: [Expr; 2] = rest.try_into().expect("two");
        vec![
            k,
            list(centre),
            torsion,
            list(tangent),
            list(normal),
            binormal,
        ]
    };
    let result = at(ev, &list(items), t, t0)?;
    ev.step("frenet", e, &result);
    Ok(Some(result))
}

/// `Evolute[C, t]`, in the plane.
pub(in crate::native) fn evolute(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [c, t] = e.args() else {
        return Ok(None);
    };
    let Some([r, d1, d2, _]) = curve(ev, c, t)? else {
        return Ok(None);
    };
    if r.len() != 2 {
        return Ok(None);
    }
    let scale = dot(&d1, &d1) / (d1[0].clone() * d2[1].clone() - d1[1].clone() * d2[0].clone());
    let normal = [-d1[1].clone(), d1[0].clone()];
    let centre: Vec<Expr> = r
        .iter()
        .zip(normal)
        .map(|(x, n)| Expr::call(Symbol::TOGETHER, [x.clone() + scale.clone() * n]))
        .collect();
    let result = list(centre);
    ev.step("evolute", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
