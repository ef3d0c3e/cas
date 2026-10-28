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
//! - `Eliminate[{p₁, …}, {t₁, …}]`, or with one symbol `t`, Xcas's
//!   `eliminate` (CAS-292): the elimination ideal with respect to the `tᵢ`
//!   in the input's other symbols, in canonical order, as `GroebnerBasis`
//!   gives it under `revlex`: `{x³ − y²}` for `x = t², y = t³`; `{}` when
//!   no relation is left, `{1}` for inconsistent equations.
//! - Past the bounds of a computation, or past its limits ([`limits`]): work
//!   from the limit on terms and coefficients of the limit on numbers, the
//!   call stays, with a warning.

use cas_num::Rational;
use cas_poly::{
    Distributed, GroebnerError, GroebnerLimits, MAX_EXPONENT, MonomialOrder,
    groebner_basis as basis_of, groebner_reduce as reduce_by,
};

use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The work a computation may do for each term a polynomial builtin's
/// product may make ([`Limits::terms`](crate::Limits)).
const WORK_PER_TERM: u64 = 64;

/// The limits of a computation (CAS-291): [`WORK_PER_TERM`] for each term
/// the limits allow, so that it is bounded as the other polynomial builtins
/// are (CAS-DOC-8 §15), and coefficients within `Limits::number_bits`, as
/// those of their powers are (CAS-156). The session's interrupt is the
/// computation's `Cancel`. Also `Solve`'s, for polynomial systems
/// (CAS-289).
pub(super) fn limits(ctx: &Context) -> GroebnerLimits {
    let terms = u64::try_from(ctx.limits().terms).unwrap_or(u64::MAX);
    GroebnerLimits {
        work: terms.saturating_mul(WORK_PER_TERM),
        coefficient_bits: ctx.limits().number_bits,
    }
}

/// Says that a computation went past its bounds, and leaves the call.
fn too_large(ev: &mut Eval<'_>, e: &Expr) -> Native {
    ev.warn(Warning::Native {
        message: "the Gröbner basis computation is past its bounds".into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// The reduced basis of `polys` under `order`; `None` past the bounds or
/// the limits, with a warning on `e`.
fn reduced_basis(
    ev: &mut Eval<'_>,
    e: &Expr,
    polys: &[Distributed],
    order: MonomialOrder,
) -> Result<Option<Vec<Distributed>>, EvalError> {
    let ctx = ev.context();
    match basis_of(polys, order, limits(ctx), ctx.interrupt_flag()) {
        Ok(basis) => Ok(Some(basis)),
        Err(GroebnerError::Cancelled) => Err(EvalError::Interrupted),
        Err(GroebnerError::TooLarge) => too_large(ev, e).map(|_| None),
    }
}

/// The list of the elements of `basis`, each primitive, in `vars`.
fn basis_list(ev: &mut Eval<'_>, basis: &[Distributed], vars: &[Expr]) -> Result<Expr, EvalError> {
    let elements: Vec<Expr> = basis
        .iter()
        .map(|g| expression(&g.primitive(), vars))
        .collect();
    ev.eval(&Expr::call(Symbol::LIST, elements))
}

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
    let Some(basis) = reduced_basis(ev, e, &polys, order)? else {
        return Ok(None);
    };
    let result = basis_list(ev, &basis, &vars)?;
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
    let ctx = ev.context();
    let remainder = match reduce_by(&p, &divisors, order, limits(ctx), ctx.interrupt_flag()) {
        Ok(remainder) => remainder,
        Err(GroebnerError::Cancelled) => return Err(EvalError::Interrupted),
        Err(GroebnerError::TooLarge) => return too_large(ev, e),
    };
    let result = ev.eval(&expression(&remainder, &vars))?;
    ev.step("GroebnerReduce", e, &result);
    Ok(Some(result))
}

/// `Eliminate[{p₁, …}, {t₁, …}]`, `Eliminate[{…}, t]`, or of one
/// polynomial or equation: the elimination ideal with respect to the `tᵢ`,
/// in the other symbols of the input, `x₁ > x₂ > …` in canonical order.
///
/// The elements free of the `tᵢ` of a `plex` basis in which they are the
/// greatest variables generate it \[CLO15, §3.1, Theorem 2\]. The answer
/// is its reduced basis under `revlex`, as `GroebnerBasis` gives it, which
/// is canonical. giac's answers are mostly resultants, which it says "may
/// reduce further", and not normalized (`[2*x-4]`).
pub(super) fn eliminate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [polys, vars] = e.args() else {
        return Ok(None);
    };
    let one = |x: &Expr| Expr::call(Symbol::LIST, [x.clone()]);
    let vars = if vars.as_symbol().is_some() {
        one(vars)
    } else {
        vars.clone()
    };
    let polys = if is_list(polys) {
        polys.clone()
    } else {
        one(polys)
    };
    let Some(eliminated) = variables(&vars) else {
        return Ok(None);
    };
    let others: Vec<Expr> = polys
        .symbols()
        .into_iter()
        .filter(|s| !s.is_builtin())
        .map(Expr::from)
        .filter(|v| !eliminated.contains(v))
        .collect();
    let all = [eliminated.as_slice(), &others].concat();
    let Some(polys) = list(ev, &polys, &all, MonomialOrder::Lex)? else {
        return Ok(None);
    };
    let Some(basis) = reduced_basis(ev, e, &polys, MonomialOrder::Lex)? else {
        return Ok(None);
    };
    // The elements free of the `tᵢ`, without their exponents.
    let k = eliminated.len();
    let ideal: Vec<Distributed> = basis
        .iter()
        .filter(|g| {
            g.terms()
                .iter()
                .all(|(m, _)| m[..k].iter().all(|&d| d == 0))
        })
        .map(|g| {
            let terms = g.terms().iter().map(|(m, c)| (m[k..].to_vec(), c.clone()));
            Distributed::new(terms.collect(), MonomialOrder::DegRevLex)
        })
        .collect();
    let Some(reduced) = reduced_basis(ev, e, &ideal, MonomialOrder::DegRevLex)? else {
        return Ok(None);
    };
    let result = basis_list(ev, &reduced, &others)?;
    ev.step("Eliminate", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
