//! Gröbner bases (CAS-DOC-21 §1, CAS-288): `GroebnerBasis[{p₁, …}, {x₁,
//! …}]` and `GroebnerReduce[p, {p₁, …}, {x₁, …}]`, Xcas's `gbasis` and
//! `greduce`, with an order after the variables: `Lexicographic` (`plex`),
//! `DegreeLexicographic` (`tdeg`) or `DegreeReverseLexicographic`
//! (`revlex`), giac's default.
//!
//! - The polynomials are expressions or equations `l = r`, read as `l − r`,
//!   expanded: sums of rational multiples of products of powers of the
//!   variables, distinct symbols, `x₁ > x₂ > …`. Another symbol, a function
//!   of a variable or another number leaves the call: giac's answers with
//!   parameters are no reduced bases.
//! - `GroebnerBasis` is cas-poly's reduced basis, each element primitive over
//!   ℤ with a positive leading coefficient, in increasing order of the
//!   leading monomials: `{x − y, 2y² − 1}` for `x² + y² − 1, x − y` under
//!   `revlex`, giac's `[2*y^2-1,x-y]`; `{1}` for the whole ring, `{}` for
//!   the zero ideal.
//! - `GroebnerReduce` is the exact remainder of the division by the list,
//!   in its order: `y/2 + 1/2` for `x³ + y²` by `{2y² − 1, x − y}`, giac's.
//! - Past the bounds of a computation the call stays, with a warning.

use cas_num::Rational;
use cas_poly::{
    Cancel, Cancelled, Distributed, GroebnerError, MAX_EXPONENT, MonomialOrder,
    groebner_basis as basis_of, groebner_reduce as reduce_by,
};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The order a symbol names.
fn monomial_order(e: &Expr) -> Option<MonomialOrder> {
    match e.as_symbol()? {
        Symbol::LEXICOGRAPHIC => Some(MonomialOrder::Lex),
        Symbol::DEGREE_LEXICOGRAPHIC => Some(MonomialOrder::DegLex),
        Symbol::DEGREE_REVERSE_LEXICOGRAPHIC => Some(MonomialOrder::DegRevLex),
        _ => None,
    }
}

/// The variables of a call: a non-empty list of distinct symbols, none a
/// builtin.
fn variables(vars: &Expr) -> Option<Vec<Expr>> {
    let vs = vars.args();
    let distinct = vs
        .iter()
        .enumerate()
        .all(|(i, v)| v.as_symbol().is_some_and(|s| !s.is_builtin()) && !vs[..i].contains(v));
    (is_list(vars) && !vs.is_empty() && distinct).then(|| vs.to_vec())
}

/// `p`, or `l − r` for an equation, expanded, as a polynomial in `vars`
/// with rational coefficients; `None` when it is not one, or an exponent is
/// past the bound of a computation. Also `Solve`'s, for polynomial systems
/// (CAS-289).
pub(super) fn distributed(
    ev: &mut Eval<'_>,
    p: &Expr,
    vars: &[Expr],
    order: MonomialOrder,
) -> Result<Option<Distributed>, EvalError> {
    let p = match (p.head().and_then(Expr::as_symbol), p.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
        _ => p.clone(),
    };
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [p]))?;
    let terms = if expanded.head() == Some(&Expr::from(Symbol::PLUS)) {
        expanded.args().to_vec()
    } else {
        vec![expanded]
    };
    let mut out = Vec::with_capacity(terms.len());
    for t in &terms {
        match term(t, vars) {
            Some(term) if term.0.iter().all(|&e| e <= MAX_EXPONENT) => out.push(term),
            _ => return Ok(None),
        }
    }
    Ok(Some(Distributed::new(out, order)))
}

/// A term `c x₁^e₁ ⋯`, `c` rational: its exponents and coefficient.
fn term(t: &Expr, vars: &[Expr]) -> Option<(Vec<u64>, Rational)> {
    let factors = if t.head() == Some(&Expr::from(Symbol::TIMES)) {
        t.args().to_vec()
    } else {
        vec![t.clone()]
    };
    let mut exps = vec![0u64; vars.len()];
    let mut c = Rational::ONE;
    for f in &factors {
        if let Some(q) = f.to_rational() {
            c = &c * &q;
            continue;
        }
        let (base, e) = match (f.head().and_then(Expr::as_symbol), f.args()) {
            (Some(Symbol::POWER), [b, n]) => (b, u64::try_from(n.as_i64()?).ok()?),
            _ => (f, 1),
        };
        let k = vars.iter().position(|v| v == base)?;
        exps[k] = exps[k].checked_add(e)?;
    }
    Some((exps, c))
}

/// `g` as an expression in `vars`.
pub(super) fn expression(g: &Distributed, vars: &[Expr]) -> Expr {
    let terms = g.terms().iter().map(|(m, c)| {
        let mut factors = vec![Expr::from(c.clone())];
        for (v, &e) in vars.iter().zip(m) {
            if e > 0 {
                factors.push(v.pow(i64::try_from(e).expect("within the bound")));
            }
        }
        Expr::call(Symbol::TIMES, factors)
    });
    Expr::call(Symbol::PLUS, terms)
}

/// The polynomials of a list, in `vars`; `None` when one is not.
fn list(
    ev: &mut Eval<'_>,
    polys: &Expr,
    vars: &[Expr],
    order: MonomialOrder,
) -> Result<Option<Vec<Distributed>>, EvalError> {
    if !is_list(polys) {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(polys.args().len());
    for p in polys.args() {
        let Some(d) = distributed(ev, p, vars, order)? else {
            return Ok(None);
        };
        out.push(d);
    }
    Ok(Some(out))
}

/// `GroebnerBasis[{p₁, …}, {x₁, …}]` or `GroebnerBasis[{…}, {…}, order]`.
pub(super) fn groebner_basis(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (polys, vars, order) = match e.args() {
        [polys, vars] => (polys, vars, Some(MonomialOrder::DegRevLex)),
        [polys, vars, order] => (polys, vars, monomial_order(order)),
        _ => return Ok(None),
    };
    let (Some(vars), Some(order)) = (variables(vars), order) else {
        return Ok(None);
    };
    let Some(polys) = list(ev, polys, &vars, order)? else {
        return Ok(None);
    };
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let basis = match basis_of(&polys, order, cancel) {
        Ok(basis) => basis,
        Err(GroebnerError::Cancelled) => return Err(EvalError::Interrupted),
        Err(GroebnerError::TooLarge) => {
            ev.warn(Warning::Native {
                message: "the Gröbner basis is past the bounds of a computation".into(),
                expr: e.clone(),
            });
            return Ok(None);
        }
    };
    let elements: Vec<Expr> = basis
        .iter()
        .map(|g| expression(&g.primitive(), &vars))
        .collect();
    let result = ev.eval(&Expr::call(Symbol::LIST, elements))?;
    ev.step("GroebnerBasis", e, &result);
    Ok(Some(result))
}

/// `GroebnerReduce[p, {p₁, …}, {x₁, …}]` or with an order after the
/// variables.
pub(super) fn groebner_reduce(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, polys, vars, order) = match e.args() {
        [p, polys, vars] => (p, polys, vars, Some(MonomialOrder::DegRevLex)),
        [p, polys, vars, order] => (p, polys, vars, monomial_order(order)),
        _ => return Ok(None),
    };
    let (Some(vars), Some(order)) = (variables(vars), order) else {
        return Ok(None);
    };
    let Some(p) = distributed(ev, p, &vars, order)? else {
        return Ok(None);
    };
    let Some(divisors) = list(ev, polys, &vars, order)? else {
        return Ok(None);
    };
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let remainder =
        reduce_by(&p, &divisors, order, cancel).map_err(|Cancelled| EvalError::Interrupted)?;
    let result = ev.eval(&expression(&remainder, &vars))?;
    ev.step("GroebnerReduce", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
