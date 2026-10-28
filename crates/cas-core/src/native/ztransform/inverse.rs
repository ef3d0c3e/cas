//! The inverse Z-transform (CAS-DOC-24 §3, CAS-309):
//! `InverseZTransform[F, z, n]`, Xcas's `invztrans(F, z, n)`.
//!
//! - **Linearity**, factors free of `z` out.
//! - `e^(a/z)` is `aⁿ/n!`.
//! - **Rational functions** by residues: `f(n) = Σ_r Res(F(z) zⁿ⁻¹, r)`,
//!   with `R = F/z = N/(z D)`, `N` of a degree at most `D`'s (otherwise `F`
//!   is no one-sided transform, and the call stays).
//!   - A pole `r ≠ 0` of multiplicity `m` gives `rⁿ P(n)`, `P(n) =
//!     1/(m−1)! dᵐ⁻¹/dzᵐ⁻¹[(z − r)ᵐ R(z) (z/r)ⁿ]` at `z = r`, the factor
//!     `(z/r)ⁿ` being 1 there.
//!   - With real coefficients, a conjugate pair `r = ρe^(±iθ)` counts
//!     `2ρⁿ (Re P cos nθ − Im P sin nθ)`, so `z/(z² + 1)` is `sin(nπ/2)`.
//!   - The pole at 0 gives Kronecker deltas: the coefficient `c_j` of
//!     `z⁻ʲ` in `R`'s principal part is `c_j KroneckerDelta[n − j + 1]`.
//!   - The denominator's roots are `Solve`'s (`roots_counted`).
//!
//! The answer is expanded: `1/(z − 2)` is `2ⁿ/2 − KroneckerDelta[n]/2`,
//! where giac answers `2^n/2*Heaviside(n-1)`, equal on the naturals.

#![allow(
    clippy::many_single_char_names,
    reason = "the transform's names: F(z), f(n), N/D, roots r, ρ and θ"
)]

use cas_num::{Complex, Natural};

use super::super::parts::parts;
use super::super::solve::roots_counted;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `InverseZTransform[F, z, n]`.
pub(in crate::native) fn inverse_z_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [image, z, n] = e.args() else {
        return Ok(None);
    };
    if z.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    let fresh = Expr::from(Symbol::new("InverseZTransform$n"));
    let Some(original) = inverse(ev, e, image, z, &fresh)? else {
        return Ok(None);
    };
    let original = original.substitute(&[(fresh, n.clone())]);
    let result = ev.eval(&Expr::call(Symbol::EXPAND, [original]))?;
    ev.step("InverseZTransform", e, &result);
    Ok(Some(result))
}

/// The inverse transform of `f` in `z` into `n`; `None` when no rule
/// applies.
fn inverse(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    z: &Expr,
    n: &Expr,
) -> Result<Option<Expr>, EvalError> {
    if !f.contains(z) {
        return Ok(Some(
            f.clone() * Expr::call(Symbol::KRONECKER_DELTA, [n.clone()]),
        ));
    }
    let (head, args) = (f.head().and_then(Expr::as_symbol), f.args());
    if head == Some(Symbol::PLUS) {
        let mut originals = Vec::with_capacity(args.len());
        for term in args {
            let Some(g) = inverse(ev, e, term, z, n)? else {
                return Ok(None);
            };
            originals.push(g);
        }
        return Ok(Some(Expr::call(Symbol::PLUS, originals)));
    }
    let factors = if head == Some(Symbol::TIMES) {
        args.to_vec()
    } else {
        vec![f.clone()]
    };
    let (constant, varying): (Vec<Expr>, Vec<Expr>) =
        factors.into_iter().partition(|g| !g.contains(z));
    let times = |v: Vec<Expr>| match v.len() {
        0 => Expr::from(1),
        1 => v.into_iter().next().expect("one"),
        _ => Expr::call(Symbol::TIMES, v),
    };
    let scale = times(constant);
    let rest = times(varying);
    // e^(a/z) is aⁿ/n!.
    if let Some(a) = reciprocal_exponent(ev, &rest, z)? {
        let g = a.pow(n.clone()) / Expr::call(Symbol::FACTORIAL, [n.clone()]);
        return Ok(Some(scale * g));
    }
    Ok(rational(ev, e, &rest, z, n)?.map(|g| scale * g))
}

/// `a` with `f = e^(a/z)`, `a` free of `z`.
fn reciprocal_exponent(ev: &mut Eval<'_>, f: &Expr, z: &Expr) -> Result<Option<Expr>, EvalError> {
    let u = match (f.head().and_then(Expr::as_symbol), f.args()) {
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => u,
        (Some(Symbol::EXP), [u]) => u,
        _ => return Ok(None),
    };
    let a = ev.eval(&(u * z))?;
    Ok((!a.contains(z)).then_some(a))
}

/// The degree of a polynomial `p` in `z` with coefficients free of it.
fn degree(ev: &mut Eval<'_>, p: &Expr, z: &Expr) -> Result<Option<i64>, EvalError> {
    let d = ev.eval(&Expr::call(Symbol::EXPONENT, [p.clone(), z.clone()]))?;
    let Some(k) = d.as_i64() else {
        return Ok(None);
    };
    for j in 0..=k {
        let c = ev.eval(&Expr::call(
            Symbol::COEFFICIENT,
            [p.clone(), z.clone(), Expr::from(j)],
        ))?;
        if c.contains(z) {
            return Ok(None);
        }
    }
    Ok(Some(k))
}

/// Whether the coefficients of `p` in `z`, of degree `k`, are real.
fn real(ev: &mut Eval<'_>, p: &Expr, z: &Expr, k: i64) -> Result<bool, EvalError> {
    for j in 0..=k {
        let c = ev.eval(&Expr::call(
            Symbol::COEFFICIENT,
            [p.clone(), z.clone(), Expr::from(j)],
        ))?;
        if parts(ev.context(), &c).is_none_or(|(_, im)| im != Expr::from(0)) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The inverse transform of a rational function of `z` by residues; `None`
/// when it is not one, not a one-sided transform, or its denominator's
/// roots are not found.
fn rational(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    z: &Expr,
    n: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [f.clone()]))?;
    let num = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let den = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    let (Some(top), Some(bottom)) = (degree(ev, &num, z)?, degree(ev, &den, z)?) else {
        return Ok(None);
    };
    if top > bottom {
        return Ok(None);
    }
    let found = if bottom == 0 {
        Vec::new()
    } else {
        let Some(found) = roots_counted(ev, e, &den, z)? else {
            return Ok(None);
        };
        found
    };
    if i64::try_from(found.len()).ok() != Some(bottom) {
        return Ok(None);
    }
    // R = N/(z D): the roots of z D, 0 among them.
    let mut distinct: Vec<(Expr, usize)> = vec![(Expr::from(0), 1)];
    for root in found {
        match distinct.iter_mut().find(|(r, _)| *r == root) {
            Some((_, m)) => *m += 1,
            None => distinct.push((root, 1)),
        }
    }
    let lead = ev.eval(&Expr::call(
        Symbol::COEFFICIENT,
        [den.clone(), z.clone(), Expr::from(bottom)],
    ))?;
    let real = real(ev, &num, z, top)? && real(ev, &den, z, bottom)?;
    let mut used = vec![false; distinct.len()];
    let mut terms = Vec::with_capacity(distinct.len());
    for i in 0..distinct.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let (root, m) = distinct[i].clone();
        // (z − r)ᵐ R(z): N over c z^(…) Π_{r′≠r} (z − r′)^(m′).
        let others: Vec<Expr> = distinct
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, (r, k))| (z - r).pow(i64::try_from(*k).unwrap_or(1)))
            .collect();
        let local = num.clone() / (lead.clone() * Expr::call(Symbol::TIMES, others));
        if root == Expr::from(0) {
            terms.extend(principal_part(ev, &local, z, n, m)?);
            continue;
        }
        let p = polynomial_part(ev, local, z, n, &root, m)?;
        let pair = if real {
            conjugate(ev, &root, &distinct, &used)?
        } else {
            None
        };
        match pair {
            Some(k) => {
                used[k] = true;
                // 2ρⁿ (Re P cos nθ − Im P sin nθ).
                let Some((re, im)) = parts(ev.context(), &p) else {
                    return Ok(None);
                };
                let rho = Expr::call(Symbol::ABS, [root.clone()]);
                let theta = Expr::call(Symbol::ARG, [root.clone()]);
                let angle = n * theta;
                let wave = re * Expr::call(Symbol::COS, [angle.clone()])
                    - im * Expr::call(Symbol::SIN, [angle]);
                terms.push(2 * rho.pow(n.clone()) * wave);
            }
            None => terms.push(root.pow(n.clone()) * p),
        }
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// `k!` as an expression.
fn factorial(k: usize) -> Expr {
    Expr::from(cas_num::Integer::from(Natural::factorial(
        u64::try_from(k).unwrap_or(0),
    )))
}

/// The terms of the pole at 0 of multiplicity `m`: `c_j
/// KroneckerDelta[n − j + 1]` for the coefficients `c_j` of `z⁻ʲ` in `R`'s
/// principal part, from the derivatives of `local = zᵐ R`.
fn principal_part(
    ev: &mut Eval<'_>,
    local: &Expr,
    z: &Expr,
    n: &Expr,
    m: usize,
) -> Result<Vec<Expr>, EvalError> {
    let mut terms = Vec::new();
    for j in 1..=m {
        let order = m - j;
        let derivative = if order == 0 {
            local.clone()
        } else {
            let spec = Expr::call(
                Symbol::LIST,
                [z.clone(), Expr::from(i64::try_from(order).unwrap_or(0))],
            );
            ev.eval(&Expr::call(Symbol::D, [local.clone(), spec]))?
        };
        let at_zero = derivative.substitute(&[(z.clone(), Expr::from(0))]);
        let c = ev.eval(&(at_zero / factorial(order)))?;
        if c != Expr::from(0) {
            let delta = Expr::call(
                Symbol::KRONECKER_DELTA,
                [n - i64::try_from(j - 1).unwrap_or(0)],
            );
            terms.push(c * delta);
        }
    }
    Ok(terms)
}

/// `P(n) = 1/(m−1)! dᵐ⁻¹/dzᵐ⁻¹[local (z/r)ⁿ]` at `z = r`, `local = (z −
/// r)ᵐ R`, expanded.
fn polynomial_part(
    ev: &mut Eval<'_>,
    local: Expr,
    z: &Expr,
    n: &Expr,
    root: &Expr,
    m: usize,
) -> Result<Expr, EvalError> {
    let g = local * (z / root).pow(n.clone());
    let derivative = if m == 1 {
        g
    } else {
        let spec = Expr::call(
            Symbol::LIST,
            [z.clone(), Expr::from(i64::try_from(m - 1).unwrap_or(0))],
        );
        ev.eval(&Expr::call(Symbol::D, [g, spec]))?
    };
    let at_root = derivative.substitute(&[(z.clone(), root.clone())]);
    ev.eval(&Expr::call(Symbol::EXPAND, [at_root / factorial(m - 1)]))
}

/// The unused root `k` that is `root`'s conjugate, `root` not real.
fn conjugate(
    ev: &mut Eval<'_>,
    root: &Expr,
    distinct: &[(Expr, usize)],
    used: &[bool],
) -> Result<Option<usize>, EvalError> {
    let Some((re, im)) = parts(ev.context(), root) else {
        return Ok(None);
    };
    let (re, im) = (ev.eval(&re)?, ev.eval(&im)?);
    if im == Expr::from(0) {
        return Ok(None);
    }
    let other = ev.eval(&(re - Expr::from(Complex::I) * im))?;
    Ok(distinct
        .iter()
        .enumerate()
        .find(|(k, (r, _))| !used[*k] && *r == other)
        .map(|(k, _)| k))
}
