//! The inverse Laplace transform (CAS-DOC-18 §4, CAS-275):
//! `InverseLaplaceTransform[F, s, t]`, Xcas's `ilaplace(F, s, t)`, with
//! `ilaplace(F, s)` into `s` and `ilaplace(F)` in `x` as for the transform.
//!
//! - **Linearity**, factors free of `s` out; a constant `c` is `c δ(t)`.
//! - **Delays**: `e^(−a s + b) G(s)`, `a ≥ 0`, is `e^b θ(t − a) g(t − a)`.
//! - **Powers of `s`**: `s^(−ν)`, `ν > 0`, is `t^(ν−1)/Γ(ν)` (giac's
//!   `ilaplace(1/sqrt(s),s,t)` is wrong).
//! - **Rational functions** by residues: `F = N/D`, `D`'s roots with their
//!   multiplicities (`roots_counted`), `f(t) = Σ_r Res(F(s) e^(st), r)`,
//!   `Res = 1/(m−1)! dᵐ⁻¹/dsᵐ⁻¹[N/(c Π_{r′≠r} (s − r′)^(m′)) e^(st)]` at
//!   `s = r`, `c` the leading coefficient. With real coefficients a
//!   conjugate pair counts twice the real part of one residue (`parts`), so
//!   `e^(αt)` times cosines and sines. The polynomial part's `c sᵏ` is the
//!   `k`-th derivative of `c δ(t)`.
//!
//! The answer expanded, its exponentials `E^u`.

#![allow(
    clippy::many_single_char_names,
    reason = "the transform's names: F(s), f(t), N/D, roots r, e^(−a s)"
)]

use std::cmp::Ordering;

use cas_num::{Complex, Natural};

use super::super::parts::parts;
use super::super::rounding::real_sign;
use super::super::solve::roots_counted;
use super::linear;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `InverseLaplaceTransform[F, s, t]`, `[F, s]` and `[F]`.
pub(in crate::native) fn inverse_laplace_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let x = Expr::symbol("x");
    let (image, s, t) = match e.args() {
        [f] => (f, x.clone(), x),
        [f, s] => (f, s.clone(), s.clone()),
        [f, s, t] => (f, s.clone(), t.clone()),
        _ => return Ok(None),
    };
    if s.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    let fresh = Expr::from(Symbol::new("InverseLaplaceTransform$t"));
    let Some(original) = inverse(ev, e, image, &s, &fresh)? else {
        return Ok(None);
    };
    let original = original.substitute(&[(fresh, t)]);
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [original]))?;
    // Exp[u] as E^u, as the rest of the CAS writes exponentials.
    let result = ev.eval(&expanded.map_bottom_up(|part| {
        match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(Symbol::EXP), [u]) => Expr::from(Symbol::E).pow(u.clone()),
            _ => part.clone(),
        }
    }))?;
    ev.step("InverseLaplaceTransform", e, &result);
    Ok(Some(result))
}

/// The inverse transform of `f` in `s` into `t`; `None` when no rule
/// applies.
fn inverse(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    s: &Expr,
    t: &Expr,
) -> Result<Option<Expr>, EvalError> {
    if !f.contains(s) {
        return Ok(Some(
            f.clone() * Expr::call(Symbol::DIRAC_DELTA, [t.clone()]),
        ));
    }
    let (head, args) = (f.head().and_then(Expr::as_symbol), f.args());
    if head == Some(Symbol::PLUS) {
        let mut originals = Vec::with_capacity(args.len());
        for term in args {
            let Some(g) = inverse(ev, e, term, s, t)? else {
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
        factors.into_iter().partition(|g| !g.contains(s));
    let times = |v: Vec<Expr>| match v.len() {
        0 => Expr::from(1),
        1 => v.into_iter().next().expect("one"),
        _ => Expr::call(Symbol::TIMES, v),
    };
    let scale = times(constant);
    // A delay: e^(−a s + b) G(s) is e^b θ(t − a) g(t − a).
    if let Some(i) = varying.iter().position(|g| exponent(g).is_some()) {
        let u = exponent(&varying[i]).expect("an exponential").clone();
        let Some((slope, b)) = linear(ev, &u, s)? else {
            return Ok(None);
        };
        let a = ev.eval(&-slope)?;
        if !matches!(
            real_sign(ev.context(), &a),
            Some(Ordering::Greater | Ordering::Equal)
        ) {
            return Ok(None);
        }
        let mut rest = varying;
        rest.remove(i);
        let Some(g) = inverse(ev, e, &times(rest), s, t)? else {
            return Ok(None);
        };
        let delayed = g.substitute(&[(t.clone(), t.clone() - a.clone())]);
        let step = Expr::call(Symbol::HEAVISIDE_THETA, [t.clone() - a]);
        return Ok(Some(scale * Expr::from(Symbol::E).pow(b) * step * delayed));
    }
    let rest = times(varying);
    // A power of s: s^(−ν), ν > 0, is t^(ν−1)/Γ(ν).
    if let (Some(Symbol::POWER), [b, p]) = (rest.head().and_then(Expr::as_symbol), rest.args())
        && b == s
        && let Some(nu) = p
            .to_rational()
            .map(|p| -p)
            .filter(cas_num::Rational::is_positive)
        && !nu.is_integer()
    {
        let nu = Expr::from(nu);
        let g = t.pow(nu.clone() - Expr::from(1)) / Expr::call(Symbol::GAMMA, [nu]);
        return Ok(Some(scale * g));
    }
    Ok(rational(ev, e, &rest, s, t)?.map(|g| scale * g))
}

/// The exponent `u` of `e^u`.
fn exponent(e: &Expr) -> Option<&Expr> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => Some(u),
        (Some(Symbol::EXP), [u]) => Some(u),
        _ => None,
    }
}

/// The inverse transform of a rational function of `s` by residues; `None`
/// when it is not one, or its denominator's roots are not found.
fn rational(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    s: &Expr,
    t: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [f.clone()]))?;
    let n = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let d = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    let polynomial = |ev: &mut Eval<'_>, p: &Expr| -> Result<Option<i64>, EvalError> {
        let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [p.clone(), s.clone()]))?;
        let Some(k) = degree.as_i64() else {
            return Ok(None);
        };
        // A polynomial: its coefficients free of s.
        for j in 0..=k {
            let c = ev.eval(&Expr::call(
                Symbol::COEFFICIENT,
                [p.clone(), s.clone(), Expr::from(j)],
            ))?;
            if c.contains(s) {
                return Ok(None);
            }
        }
        Ok(Some(k))
    };
    let (Some(_), Some(degree)) = (polynomial(ev, &n)?, polynomial(ev, &d)?) else {
        return Ok(None);
    };
    if degree == 0 {
        return Ok(None);
    }
    // The polynomial part, deltas and their derivatives; the proper rest.
    let qr = ev.eval(&Expr::call(
        Symbol::POLYNOMIAL_QUOTIENT_REMAINDER,
        [n, d.clone(), s.clone()],
    ))?;
    let [q, r] = qr.args() else {
        return Ok(None);
    };
    let mut terms = Vec::new();
    if *q != Expr::from(0) {
        let Some(top) = polynomial(ev, q)? else {
            return Ok(None);
        };
        for k in 0..=top {
            let c = ev.eval(&Expr::call(
                Symbol::COEFFICIENT,
                [q.clone(), s.clone(), Expr::from(k)],
            ))?;
            if c == Expr::from(0) {
                continue;
            }
            let delta = Expr::call(Symbol::DIRAC_DELTA, [t.clone()]);
            let derivative = if k == 0 {
                delta
            } else {
                Expr::call(
                    Symbol::D,
                    [delta, Expr::call(Symbol::LIST, [t.clone(), Expr::from(k)])],
                )
            };
            terms.push(c * derivative);
        }
    }
    if *r != Expr::from(0) {
        let Some(sum) = residues(ev, e, r, &d, s, t, degree)? else {
            return Ok(None);
        };
        terms.push(sum);
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// `Σ_r Res(N/D e^(st), r)` over `D`'s roots, of degree `degree`, a
/// conjugate pair as twice one's real part where `N` and `D` are real.
#[allow(
    clippy::too_many_arguments,
    reason = "the fraction, its variables and degree"
)]
fn residues(
    ev: &mut Eval<'_>,
    e: &Expr,
    n: &Expr,
    d: &Expr,
    s: &Expr,
    t: &Expr,
    degree: i64,
) -> Result<Option<Expr>, EvalError> {
    let Some(found) = roots_counted(ev, e, d, s)? else {
        return Ok(None);
    };
    if i64::try_from(found.len()).ok() != Some(degree) {
        return Ok(None);
    }
    let mut distinct: Vec<(Expr, usize)> = Vec::new();
    for root in found {
        match distinct.iter_mut().find(|(r, _)| *r == root) {
            Some((_, m)) => *m += 1,
            None => distinct.push((root, 1)),
        }
    }
    let lead = ev.eval(&Expr::call(
        Symbol::COEFFICIENT,
        [d.clone(), s.clone(), Expr::from(degree)],
    ))?;
    let ctx = ev.context();
    let real = |p: &Expr| -> bool {
        (0..=degree).all(|k| {
            let c = Expr::call(Symbol::COEFFICIENT, [p.clone(), s.clone(), Expr::from(k)]);
            ctx.eval(&c)
                .ok()
                .and_then(|c| parts(ctx, &c))
                .is_some_and(|(_, im)| im == Expr::from(0))
        })
    };
    let real = real(n) && real(d);
    let mut used = vec![false; distinct.len()];
    let mut terms = Vec::with_capacity(distinct.len());
    for i in 0..distinct.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let (root, m) = distinct[i].clone();
        // N e^(st) over c Π_{r′≠r} (s − r′)^(m′).
        let others = distinct
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, (r, k))| (s.clone() - r.clone()).pow(i64::try_from(*k).unwrap_or(1)));
        let g = n.clone() * Expr::call(Symbol::EXP, [s.clone() * t.clone()])
            / (lead.clone() * Expr::call(Symbol::TIMES, others.collect::<Vec<_>>()));
        let order = Expr::call(
            Symbol::LIST,
            [s.clone(), Expr::from(i64::try_from(m - 1).unwrap_or(0))],
        );
        let derivative = if m == 1 {
            g
        } else {
            ev.eval(&Expr::call(Symbol::D, [g, order]))?
        };
        let factorial = Expr::from(cas_num::Integer::from(Natural::factorial(
            u64::try_from(m - 1).unwrap_or(0),
        )));
        let residue =
            ev.eval(&(derivative.substitute(&[(s.clone(), root.clone())]) / factorial))?;
        // A conjugate pair: twice the real part.
        let pair = if real {
            conjugate(ev, &root, &distinct, &used)?
        } else {
            None
        };
        match pair {
            Some(k) => {
                used[k] = true;
                let Some((re, _)) = parts(ev.context(), &residue) else {
                    return Ok(None);
                };
                terms.push(Expr::from(2) * re);
            }
            None => terms.push(residue),
        }
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
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
