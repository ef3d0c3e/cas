//! Sums of rational functions through `PolyGamma` (CAS-DOC-31 §6, CAS-359):
//! what Gosper's algorithm leaves of a rational summand.
//!
//! - `Apart` (CAS-DOC-29 §2) splits it into a polynomial and terms
//!   `c (α k + β)^(−j)`, `c`, `α` and `β` free of `k`; a part of another
//!   shape, as an irreducible quadratic's, leaves the sum.
//! - From `a` to `b`: the polynomial by its antidifference, and with `β′ =
//!   β/α`, `Σ_{k=a}^{b} (k + β′)^(−j) = (−1)^(j−1)/(j − 1)! (ψ^(j−1)(b + 1 +
//!   β′) − ψ^(j−1)(a + β′))`: `sum(1/k,k,1,n)` is `Psi(n+1) + euler_gamma`,
//!   as giac's.
//! - To infinity, when the polynomial is 0 and the coefficients of the `j =
//!   1` terms add up to 0, so that it converges: `Σ_{k=a}^{∞} = −Σ c′
//!   (−1)^(j−1)/(j − 1)! ψ^(j−1)(a + β′)`, the `ψ(b + 1 + β′)` growing alike
//!   as `log b`; `PolyGamma`'s exact values close it: `Σ_{k≥1} 1/(k (2k +
//!   1)) = 2 − 2 log 2` and `Σ_{k≥1} 1/(2k − 1)² = π²/8`.
//! - A pole at `k = a + i` for an integer `i ≥ 0`, `a + β′` an integer
//!   `≤ 0`, leaves the sum.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most terms of `Apart`'s split summed.
const MAX_TERMS: usize = 64;

/// A term `c (k + shift)^(−j)`.
struct Term {
    c: Expr,
    shift: Expr,
    j: i64,
}

/// `c (α k + β)^(−j)` as a [`Term`], `c`, `α` and `β` free of `k`.
fn reciprocal_power(ev: &mut Eval<'_>, piece: &Expr, k: &Expr) -> Result<Option<Term>, EvalError> {
    let factors: Vec<Expr> = if piece.head() == Some(&Expr::from(Symbol::TIMES)) {
        piece.args().to_vec()
    } else {
        vec![piece.clone()]
    };
    let mut found = None;
    let mut rest = Vec::new();
    for f in factors {
        if !f.contains(k) {
            rest.push(f);
            continue;
        }
        let (Some(Symbol::POWER), [base, n]) = (f.head().and_then(Expr::as_symbol), f.args())
        else {
            return Ok(None);
        };
        let Some(j) = n.as_i64().filter(|&j| j < 0) else {
            return Ok(None);
        };
        if found.is_some() {
            return Ok(None);
        }
        found = Some((base.clone(), -j));
    }
    let Some((linear, j)) = found else {
        return Ok(None);
    };
    let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [linear.clone(), k.clone()]))?;
    if degree.as_i64() != Some(1) {
        return Ok(None);
    }
    let coefficient = |n: i64| {
        Expr::call(
            Symbol::COEFFICIENT,
            [linear.clone(), k.clone(), Expr::from(n)],
        )
    };
    let alpha = ev.eval(&coefficient(1))?;
    let beta = ev.eval(&coefficient(0))?;
    if alpha.contains(k) || beta.contains(k) || alpha == Expr::from(0) {
        return Ok(None);
    }
    let c: Expr = rest.into_iter().product();
    Ok(Some(Term {
        c: ev.eval(&(c / alpha.pow(j)))?,
        shift: ev.eval(&(beta / alpha))?,
        j,
    }))
}

/// Whether `p` is a polynomial in `k`: its denominator is free of it.
fn polynomial(ev: &mut Eval<'_>, p: &Expr, k: &Expr) -> Result<bool, EvalError> {
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [p.clone()]))?;
    let denominator = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    Ok(!denominator.contains(k))
}

/// The polynomial parts and the terms of a split.
type Split = (Vec<Expr>, Vec<Term>);

/// `f` split by `Apart`: its polynomial part and its terms; `None` for
/// another shape.
fn split(ev: &mut Eval<'_>, f: &Expr, k: &Expr) -> Result<Option<Split>, EvalError> {
    let parts = ev.eval(&Expr::call(Symbol::APART, [f.clone(), k.clone()]))?;
    if parts.head() == Some(&Expr::from(Symbol::APART)) {
        return Ok(None);
    }
    let pieces: Vec<Expr> = if parts.head() == Some(&Expr::from(Symbol::PLUS)) {
        parts.args().to_vec()
    } else {
        vec![parts]
    };
    if pieces.len() > MAX_TERMS {
        return Ok(None);
    }
    let (mut polynomials, mut terms) = (Vec::new(), Vec::new());
    for piece in pieces {
        if let Some(t) = reciprocal_power(ev, &piece, k)? {
            terms.push(t);
        } else if polynomial(ev, &piece, k)? {
            polynomials.push(piece);
        } else {
            return Ok(None);
        }
    }
    Ok((!terms.is_empty()).then_some((polynomials, terms)))
}

/// Whether `a + shift` is an integer `≤ 0`: a pole at some `k ≥ a`.
fn meets_pole(ev: &mut Eval<'_>, a: &Expr, shift: &Expr) -> Result<bool, EvalError> {
    let start = ev.eval(&(a.clone() + shift.clone()))?;
    Ok(start
        .to_rational()
        .is_some_and(|q| q.is_integer() && !q.is_positive()))
}

/// `c (−1)^(j−1)/(j − 1)! ψ^(j−1)(x)`.
fn weighted_polygamma(t: &Term, x: Expr) -> Expr {
    let order = Expr::from(t.j - 1);
    let factorial = Expr::call(Symbol::FACTORIAL, [order.clone()]);
    let sign = if t.j % 2 == 0 { -1 } else { 1 };
    Expr::from(sign) * t.c.clone() / factorial * Expr::call(Symbol::POLY_GAMMA, [x, order])
}

/// `Σ_{k=a}^{b} f`, or `None`.
pub(super) fn finite(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some((polynomials, terms)) = split(ev, f, k)? else {
        return Ok(None);
    };
    let mut total = Vec::new();
    if !polynomials.is_empty() {
        let p: Expr = polynomials.into_iter().sum();
        let Some(s) = super::antidifference(ev, e, &p, k, a, b)? else {
            return Ok(None);
        };
        total.push(s);
    }
    for t in &terms {
        if meets_pole(ev, a, &t.shift)? {
            return Ok(None);
        }
        let upper = b.clone() + Expr::from(1) + t.shift.clone();
        let lower = a.clone() + t.shift.clone();
        total.push(weighted_polygamma(t, upper) - weighted_polygamma(t, lower));
    }
    Ok(Some(ev.eval(&total.into_iter().sum::<Expr>())?))
}

/// `Σ_{k=a}^{∞} f` when it converges, or `None`.
pub(super) fn infinite(
    ev: &mut Eval<'_>,
    f: &Expr,
    k: &Expr,
    a: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some((polynomials, terms)) = split(ev, f, k)? else {
        return Ok(None);
    };
    if !polynomials.is_empty() {
        return Ok(None);
    }
    let first: Expr = terms.iter().filter(|t| t.j == 1).map(|t| t.c.clone()).sum();
    let first = ev.eval(&Expr::call(Symbol::TOGETHER, [first]))?;
    if first != Expr::from(0) {
        return Ok(None);
    }
    let mut total = Vec::new();
    for t in &terms {
        if meets_pole(ev, a, &t.shift)? {
            return Ok(None);
        }
        total.push(-weighted_polygamma(t, a.clone() + t.shift.clone()));
    }
    Ok(Some(ev.eval(&total.into_iter().sum::<Expr>())?))
}

#[cfg(test)]
mod tests;
