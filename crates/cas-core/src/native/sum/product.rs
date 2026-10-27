//! Symbolic products (CAS-DOC-13 §5, CAS-227): `Product[f, k, a, b]` with
//! symbolic bounds, or too many terms to write out.
//!
//! - Factor by factor: the factors rational in `k` together, each other
//!   alone; those without a product here stay under one `Product`.
//! - `c^g(k)` and `e^g(k)`, `c` free of `k`: `c^Σg` by `Sum` (giac leaves
//!   `product(2^k,k,1,n)`); `h^n` for an integer `n`, `(Π h)^n`.
//! - A rational function of `k`, factored over ℚ with the parameters,
//!   `c ∏ (αᵢ k + βᵢ)^eᵢ` and factors free of `k`: those and `c` and the
//!   `αᵢ^eᵢ` to the power `b − a + 1`, and with `γ = β/α`
//!   `Π_{k=a}^{b} (k + γ) = (b + γ)!/(a − 1 + γ)!`, as giac's `(a+n)!/(a!)`
//!   for `k + a`. Factorials whose shifts differ by an integer are one,
//!   `(b + γ + j)! = (b + γ)! (b + γ + 1) ⋯ (b + γ + j)`, so `(k + 1)/k`
//!   is `n + 1`. A factor of the numerator that vanishes at an integer
//!   `k ≥ a` makes the product 0, as giac's `product(k,k,0,n)`; one of the
//!   denominator leaves it, and so does a factor of degree 2 or more.

use cas_num::algebra::Integers;
use cas_poly::SparseRing;

use super::super::polynomial::{convert, factored};
use super::gosper::rational_in;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::Fraction;
use crate::symbol::Symbol;

/// The largest integer difference of two shifts whose factorials are one.
const MAX_OFFSET: i64 = 64;

/// `b − a + 1`, built.
fn count(a: &Expr, b: &Expr) -> Expr {
    b.clone() - a.clone() + Expr::from(1)
}

/// `Π_{k=a}^{b} f` for an `f` with `k` in it; `None` when no factor has a
/// product here.
pub(super) fn product(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let factors: Vec<Expr> = match f.head().and_then(Expr::as_symbol) {
        Some(Symbol::TIMES) => f.args().to_vec(),
        _ => vec![f.clone()],
    };
    let (rational, others): (Vec<Expr>, Vec<Expr>) =
        factors.into_iter().partition(|x| rational_in(x, k));
    let mut closed = Vec::new();
    let mut open = Vec::new();
    if rational.iter().any(|x| x.contains(k)) {
        let r = Expr::call(Symbol::TIMES, rational.clone());
        match rational_product(ev, e, &r, k, a, b)? {
            Some(p) => closed.push(p),
            None => open.extend(rational),
        }
    } else {
        closed.extend(rational.into_iter().map(|x| x.pow(count(a, b))));
    }
    for x in others {
        match factor(ev, e, &x, k, a, b)? {
            Some(p) => closed.push(p),
            None => open.push(x),
        }
    }
    if closed.is_empty() {
        return Ok(None);
    }
    if !open.is_empty() {
        closed.push(Expr::call(
            Symbol::PRODUCT,
            [
                Expr::call(Symbol::TIMES, open),
                k.clone(),
                a.clone(),
                b.clone(),
            ],
        ));
    }
    Ok(Some(Expr::call(Symbol::TIMES, closed)))
}

/// The product of a factor `x` that is not rational in `k`.
fn factor(
    ev: &mut Eval<'_>,
    e: &Expr,
    x: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    Ok(match (x.head().and_then(Expr::as_symbol), x.args()) {
        (Some(Symbol::POWER), [c, g]) if !c.contains(k) => {
            exponent_sum(ev, g, k, a, b)?.map(|s| c.pow(s))
        }
        (Some(Symbol::EXP), [g]) => {
            exponent_sum(ev, g, k, a, b)?.map(|s| Expr::call(Symbol::EXP, [s]))
        }
        (Some(Symbol::POWER), [h, n]) if n.as_i64().is_some() => {
            product(ev, e, h, k, a, b)?.map(|p| p.pow(n.clone()))
        }
        _ => None,
    })
}

/// `Σ_{k=a}^{b} g`, when `Sum` finds it.
fn exponent_sum(
    ev: &mut Eval<'_>,
    g: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let s = ev.eval(&Expr::call(
        Symbol::SUM,
        [g.clone(), k.clone(), a.clone(), b.clone()],
    ))?;
    Ok((!s.contains(&Expr::from(Symbol::SUM))).then_some(s))
}

/// The product of `r`, a rational function of `k`; `None` when it has a
/// factor of degree 2 or more, or a pole at an integer `k ≥ a`.
fn rational_product(
    ev: &mut Eval<'_>,
    e: &Expr,
    r: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some(p) = convert(ev, e, &[r.clone(), k.clone()], false)? else {
        return Ok(None);
    };
    let Some(kv) = p.variable(k) else {
        return Ok(None);
    };
    if p.generators()
        .iter()
        .enumerate()
        .any(|(j, g)| j != kv && g.contains(k))
    {
        return Ok(None);
    }
    let Fraction::Integers { num, den } = &p.fractions()[0] else {
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, p.generators().len());
    let Some(powers) = factored(ev, e, &ring, num, den)? else {
        return Ok(None);
    };
    let n = count(a, b);
    let mut parts = vec![powers.constant().pow(n.clone())];
    // The shifts γ of the linear factors k + γ, with their exponents.
    let mut shifts: Vec<(Expr, i64)> = Vec::new();
    for (g, m) in &powers.factors {
        let power = Expr::from(*m) * n.clone();
        match ring.degree(g, kv).unwrap_or(0) {
            0 => parts.push(p.polynomial_expr(g).pow(power)),
            1 => {
                let (alpha, beta) = (ring.coeff(g, kv, 1), ring.coeff(g, kv, 0));
                parts.push(p.polynomial_expr(&alpha).pow(power));
                shifts.push((ev.eval(&p.ratio_expr(&beta, &alpha))?, *m));
            }
            _ => return Ok(None),
        }
    }
    // Zeros and poles at integers k = −γ ≥ a: a + γ an integer, not
    // positive. One past a numeric b, or a pole, leaves the product.
    if let Some(a) = a.to_rational() {
        let mut zero = false;
        for (gamma, m) in &shifts {
            let Some(gamma) = gamma.to_rational() else {
                continue;
            };
            let at = &a + &gamma;
            if at.is_integer() && !at.is_positive() {
                let past = b.to_rational().is_some_and(|b| b < -gamma);
                if *m < 0 || past {
                    return Ok(None);
                }
                zero = true;
            }
        }
        if zero {
            return Ok(Some(Expr::from(0)));
        }
    }
    parts.extend(factorials(ev, &shifts, a, b)?);
    Ok(Some(Expr::call(Symbol::TIMES, parts)))
}

/// `Π ((b + γ)!/(a − 1 + γ)!)^e` over the shifts `(γ, e)`, those that differ
/// by an integer through the factorial of the smallest, built.
fn factorials(
    ev: &mut Eval<'_>,
    shifts: &[(Expr, i64)],
    a: &Expr,
    b: &Expr,
) -> Result<Vec<Expr>, EvalError> {
    // Classes: a shift, and the offsets from it of the members with their
    // exponents.
    let mut classes: Vec<(Expr, Vec<(i64, i64)>)> = Vec::new();
    for (gamma, m) in shifts {
        let mut placed = false;
        for (base, members) in &mut classes {
            let offset = ev.eval(&Expr::call(
                Symbol::TOGETHER,
                [gamma.clone() - base.clone()],
            ))?;
            if let Some(d) = offset.as_i64().filter(|d| d.abs() <= MAX_OFFSET) {
                members.push((d, *m));
                placed = true;
                break;
            }
        }
        if !placed {
            classes.push((gamma.clone(), vec![(0, *m)]));
        }
    }
    let mut parts = Vec::new();
    for (base, members) in classes {
        let low = members.iter().map(|&(d, _)| d).min().unwrap_or(0);
        let base = base + Expr::from(low);
        let total: i64 = members.iter().map(|&(_, m)| m).sum();
        if total != 0 {
            let upper = Expr::call(Symbol::FACTORIAL, [b.clone() + base.clone()]);
            let lower = Expr::call(
                Symbol::FACTORIAL,
                [a.clone() + base.clone() - Expr::from(1)],
            );
            parts.push((upper / lower).pow(total));
        }
        // (b + γ + j)!/(a − 1 + γ + j)! over those of γ: (b + γ + 1) ⋯
        // (b + γ + j) over (a + γ) ⋯ (a + γ + j − 1).
        for (d, m) in members {
            for j in 0..d - low {
                let up = b.clone() + base.clone() + Expr::from(j + 1);
                let down = a.clone() + base.clone() + Expr::from(j);
                parts.push((up / down).pow(m));
            }
        }
    }
    Ok(parts)
}
