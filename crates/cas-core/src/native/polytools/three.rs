//! Polynomial tools III (CAS-540, CAS-543): Xcas's `rootof`, `fcoeff`,
//! `homogeneize`, `total_degree`, `factor_xn` and `is_irreducible`. The
//! variable is `x` unless said otherwise, as in giac.
//!
//! - **`RootOf[p, q]`**, `rootof`: the value `p(α)` at the root `α` of
//!   `q`, both lists of coefficients, the leading one first; it stays a
//!   call, as there is no algebraic-number type (R1.2-003). Polynomials in
//!   one variable are rewritten to their lists (giac's
//!   `rootof(poly1[1,2],poly1[1,0,-2])` for `x + 2` and `x² − 2`), and
//!   `RootOf[q]` is `RootOf[{1, 0}, q]`, `α` itself.
//! - **`FactoredFromRoots[l, v]`**, `fcoeff`: `Π (v − rᵢ)^mᵢ` for `l = {r₁,
//!   m₁, r₂, m₂, …}`, `v` by default `x`, a multiplicity negative or 0
//!   allowed: giac's `(x-1)^2*x*(x-3)^-1` for `[1,2,0,1,3,-1]`.
//! - **`Homogenize[P, v]`**, `homogeneize`: with `d` the total degree of
//!   `P` in all its symbols, each monomial `m` times `v^(d − deg m)`, `v` by
//!   default `t` and any expression, substituted as given; a rational
//!   function's numerator and denominator each on its own: giac's
//!   `x^3+x*y*z+z^3` for `x³ + x y + 1` and `z`, `(x+t)/x` for `1/x + 1`.
//! - **`TotalDegree[P, vars]`**, `total_degree`: the largest sum of the
//!   exponents of `vars` (a list or one symbol, `{x}` when left out) in a
//!   monomial of the expanded `P`, 0 for 0; a `P` that is no polynomial in
//!   `vars` leaves the call (giac's `1` for `1/x + y`).
//! - **`FactorXn[P]`**, `factor_xn`: `xⁿ (P x⁻ⁿ)`, `n` the degree of `P`
//!   in `x`, the second factor expanded in powers `x⁻ᵏ` and left so (`Times`
//!   does not distribute): giac's `x^3*(1+x^-1)` for `x³ + x²`; a `P` free
//!   of `x` as it is; a `P` that is no polynomial in `x` leaves the call
//!   (giac's `1+x^-1` for `1/x + 1`).
//! - **`IrreduciblePolynomialQ[P]`**, `is_irreducible` (CAS-543): 1 when
//!   `FactorList` (over ℚ, whatever `with_sqrt`) finds exactly one factor
//!   that is not constant, of multiplicity 1, the constant content aside
//!   (`2 x + 2` is irreducible), and for a constant, as giac; 0 otherwise.
//!   giac's `1` for `x y` and `x y + x` is not ours. A `P` that is no
//!   polynomial leaves the call (giac errs on `1/x`).
//!
//! Each call that answers takes one step, giac's name.

use super::two::{answer, is_list};
use super::{call, coefficients};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The terms of an expanded expression; none for 0.
fn terms(p: &Expr) -> Vec<Expr> {
    if *p == Expr::from(0) {
        Vec::new()
    } else if p.head() == Some(&Expr::from(Symbol::PLUS)) {
        p.args().to_vec()
    } else {
        vec![p.clone()]
    }
}

/// The degree of a monomial in the variables `var` tells apart, the
/// factors `free` tells free of them counting 0; `None` for a factor
/// neither, or a power of a variable that is no natural number.
fn monomial_degree(
    m: &Expr,
    var: &dyn Fn(&Expr) -> bool,
    free: &dyn Fn(&Expr) -> bool,
) -> Option<i64> {
    let factors = if m.head() == Some(&Expr::from(Symbol::TIMES)) {
        m.args()
    } else {
        std::slice::from_ref(m)
    };
    let mut d = 0i64;
    for f in factors {
        if var(f) {
            d += 1;
            continue;
        }
        if free(f) {
            continue;
        }
        match (f.head().and_then(Expr::as_symbol), f.args()) {
            (Some(Symbol::POWER), [b, n]) if var(b) => d += n.as_i64().filter(|n| *n >= 0)?,
            _ => return None,
        }
    }
    Some(d)
}

/// A symbol that is no numeric constant (`Pi` is one).
fn is_variable(e: &Expr) -> bool {
    e.as_symbol().is_some() && !e.is_numeric()
}

/// The symbols of `e` that are no numeric constant, heads aside.
fn symbols(e: &Expr, out: &mut Vec<Expr>) {
    if is_variable(e) {
        if !out.contains(e) {
            out.push(e.clone());
        }
    } else {
        for a in e.args() {
            symbols(a, out);
        }
    }
}

/// The degrees of the terms of `Expand[p]` in all its symbols, with the
/// terms; `None` for a `p` that is no polynomial.
fn monomials(ev: &mut Eval<'_>, p: &Expr) -> Result<Option<Vec<(Expr, i64)>>, EvalError> {
    let expanded = ev.eval(&call(Symbol::EXPAND, vec![p.clone()]))?;
    Ok(terms(&expanded)
        .into_iter()
        .map(|t| monomial_degree(&t, &is_variable, &Expr::is_numeric).map(|d| (t, d)))
        .collect())
}

/// `{c_n, …, c_0}`: the coefficients of `p` in `v`, the leading one first.
fn descending(ev: &mut Eval<'_>, p: &Expr, v: &Expr) -> Result<Option<Expr>, EvalError> {
    if is_list(p) {
        return Ok(Some(p.clone()));
    }
    Ok(coefficients(ev, p, v)?.map(|cs| call(Symbol::LIST, cs.into_iter().rev().collect())))
}

/// `RootOf[p, q]` and `RootOf[q]`.
pub(in crate::native) fn root_of(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, q) = match e.args() {
        [q] => (
            call(Symbol::LIST, vec![Expr::from(1), Expr::from(0)]),
            q.clone(),
        ),
        [p, q] => (p.clone(), q.clone()),
        _ => return Ok(None),
    };
    if is_list(&p) && is_list(&q) && e.args().len() == 2 {
        return Ok(None);
    }
    let mut found = Vec::new();
    for a in [&p, &q] {
        if !is_list(a) {
            symbols(a, &mut found);
        }
    }
    let v = match &found[..] {
        [] => Expr::symbol("x"),
        [v] => v.clone(),
        _ => return Ok(None),
    };
    let (Some(p), Some(q)) = (descending(ev, &p, &v)?, descending(ev, &q, &v)?) else {
        return Ok(None);
    };
    answer(ev, "rootof", e, &call(Symbol::ROOT_OF, vec![p, q]))
}

/// `FactoredFromRoots[l, v]`.
pub(in crate::native) fn factored_from_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, v) = match e.args() {
        [l] => (l, Expr::symbol("x")),
        [l, v] if v.as_symbol().is_some() => (l, v.clone()),
        _ => return Ok(None),
    };
    if !is_list(l) || l.args().len() % 2 != 0 {
        return Ok(None);
    }
    let factors: Vec<Expr> = l
        .args()
        .chunks(2)
        .map(|rm| {
            call(
                Symbol::POWER,
                vec![v.clone() - rm[0].clone(), rm[1].clone()],
            )
        })
        .collect();
    answer(ev, "fcoeff", e, &call(Symbol::TIMES, factors))
}

/// `p` homogeneous in `v` of the total degree of `p`; `None` for a `p`
/// that is no polynomial.
fn homogeneous(ev: &mut Eval<'_>, p: &Expr, v: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some(ms) = monomials(ev, p)? else {
        return Ok(None);
    };
    let d = ms.iter().map(|(_, k)| *k).max().unwrap_or(0);
    Ok(Some(call(
        Symbol::PLUS,
        ms.into_iter()
            .map(|(m, k)| m * call(Symbol::POWER, vec![v.clone(), Expr::from(d - k)]))
            .collect(),
    )))
}

/// `Homogenize[P, v]`.
pub(in crate::native) fn homogenize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, v) = match e.args() {
        [p] => (p, Expr::symbol("t")),
        [p, v] => (p, v.clone()),
        _ => return Ok(None),
    };
    if is_list(p) {
        return Ok(None);
    }
    let t = ev.eval(&call(Symbol::TOGETHER, vec![p.clone()]))?;
    let num = ev.eval(&call(Symbol::NUMERATOR, vec![t.clone()]))?;
    let den = ev.eval(&call(Symbol::DENOMINATOR, vec![t]))?;
    let (Some(num), Some(den)) = (homogeneous(ev, &num, &v)?, homogeneous(ev, &den, &v)?) else {
        return Ok(None);
    };
    answer(ev, "homogeneize", e, &(num * den.pow(-1)))
}

/// `TotalDegree[P, vars]`.
pub(in crate::native) fn total_degree(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, vars) = match e.args() {
        [p] => (p, vec![Expr::symbol("x")]),
        [p, l] if is_list(l) => (p, l.args().to_vec()),
        [p, v] => (p, vec![v.clone()]),
        _ => return Ok(None),
    };
    if vars.iter().any(|v| v.as_symbol().is_none()) {
        return Ok(None);
    }
    let expanded = ev.eval(&call(Symbol::EXPAND, vec![p.clone()]))?;
    let var = |f: &Expr| vars.contains(f);
    let free = |f: &Expr| !vars.iter().any(|v| f.contains(v));
    let mut d = 0;
    for t in terms(&expanded) {
        let Some(k) = monomial_degree(&t, &var, &free) else {
            return Ok(None);
        };
        d = d.max(k);
    }
    let result = Expr::from(d);
    ev.step("total_degree", e, &result);
    Ok(Some(result))
}

/// `FactorXn[P]`.
pub(in crate::native) fn factor_xn(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let Some(cs) = coefficients(ev, p, &x)? else {
        return Ok(None);
    };
    let n = i64::try_from(cs.len()).expect("small") - 1;
    if n <= 0 {
        return answer(ev, "factor_xn", e, p);
    }
    let rest = call(
        Symbol::PLUS,
        cs.into_iter()
            .zip(0i64..)
            .map(|(c, k)| c * x.pow(k - n))
            .collect(),
    );
    let result = call(Symbol::TIMES, vec![x.pow(n), rest]);
    answer(ev, "factor_xn", e, &result)
}

/// `IrreduciblePolynomialQ[P]` and `IrreduciblePolynomialQ[P, x]`.
pub(in crate::native) fn irreducible_polynomial_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let p = match e.args() {
        [p] => p,
        [p, x] if x.as_symbol().is_some() => p,
        _ => return Ok(None),
    };
    if monomials(ev, p)?.is_none() {
        return Ok(None);
    }
    let factors = ev.eval(&call(Symbol::FACTOR_LIST, vec![p.clone()]))?;
    if !is_list(&factors) || factors.args().len() % 2 != 0 {
        return Ok(None);
    }
    let mut count = 0i64;
    for fm in factors.args().chunks(2) {
        if fm[0].is_numeric() {
            continue;
        }
        let Some(m) = fm[1].as_i64().filter(|m| *m > 0) else {
            return Ok(None);
        };
        count += m;
    }
    // A constant has no factor and is irreducible, as giac's.
    let result = Expr::from(i64::from(count <= 1));
    ev.step("is_irreducible", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
