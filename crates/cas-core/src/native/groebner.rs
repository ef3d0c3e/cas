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
//! - `GroebnerBasis[{…}, {…}, RationalUnivariateRepresentation]`, the option
//!   in any order with the monomial order, giac's `rur` (CAS-541): for a
//!   zero-dimensional ideal in two or more variables, `{rur, l, f(t), f'(t),
//!   g₁(t), …}`, the points `xᵢ = gᵢ(l)/f'(l)` at the roots of `f`, always
//!   from `revlex`; `t` is `l` when `l` is a variable, else `x₁`. Otherwise
//!   the basis.
//! - `IdealMemberQ[p, {p₁, …}, {x₁, …}]`, with an order after the
//!   variables, or without them (the input's symbols in canonical order),
//!   giac's `in_ideal`: 1 when `p` reduces to 0 by the reduced basis, else
//!   0; for a list of `p`, the list. giac reduces by the list as given, so
//!   its 0 for `x²` in `(x + y, x − y)` is 1 here.
//! - `GroebnerVariables[{p₁, …}]`, giac's `gbasis_vars`: the symbols in
//!   order of first appearance, stably sorted by `∏ (deg_x pᵢ + 1)`; `{}` for
//!   anything but a list, and the call stays when a `pᵢ` is no polynomial in
//!   symbols.
//! - `Syzygies[{p₁, …, pₘ}]` or `Syzygies[{…}, {x₁, …}]`, giac's `syzygy`:
//!   `{S, V, A, X}` with `X` the variables (else the symbols in canonical
//!   order), `A = {A0, …, A(m−1)}`, `V = A ∪ X`, and `S` cas-poly's
//!   [`syzygies`](cas_poly::syzygies) in `V`, each primitive: the elements
//!   linear in the `Aᵢ` generate the module of syzygies. The call stays when
//!   a variable is one of the `Aᵢ`.
//! - Past the bounds of a computation, or past its limits ([`limits`]): work
//!   from the limit on terms and coefficients of the limit on numbers, the
//!   call stays, with a warning.

use cas_num::Rational;
use cas_poly::{
    Distributed, GroebnerError, GroebnerLimits, MAX_EXPONENT, MonomialOrder, Rur,
    groebner_basis as basis_of, groebner_reduce as reduce_by, rational_univariate,
    syzygies as syzygies_of,
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

/// The value of a computation; `None` past the bounds or the limits, with a
/// warning on `e`.
fn computed<T>(
    ev: &mut Eval<'_>,
    e: &Expr,
    outcome: Result<T, GroebnerError>,
) -> Result<Option<T>, EvalError> {
    match outcome {
        Ok(value) => Ok(Some(value)),
        Err(GroebnerError::Cancelled) => Err(EvalError::Interrupted),
        Err(GroebnerError::TooLarge) => too_large(ev, e).map(|_| None),
    }
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

/// `GroebnerBasis[{p₁, …}, {x₁, …}]`, or with an order and
/// `RationalUnivariateRepresentation` after the variables, in any order.
pub(super) fn groebner_basis(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [polys, vars, options @ ..] = e.args() else {
        return Ok(None);
    };
    let (mut order, mut rur) = (None, false);
    for option in options {
        if option.as_symbol() == Some(Symbol::RATIONAL_UNIVARIATE_REPRESENTATION) && !rur {
            rur = true;
        } else if let (Some(o), None) = (monomial_order(option), order) {
            order = Some(o);
        } else {
            return Ok(None);
        }
    }
    let order = order.unwrap_or(MonomialOrder::DegRevLex);
    let Some(vars) = variables(vars) else {
        return Ok(None);
    };
    if rur {
        let Some(polys) = list(ev, polys, &vars, MonomialOrder::DegRevLex)? else {
            return Ok(None);
        };
        let ctx = ev.context();
        let outcome = rational_univariate(&polys, limits(ctx), ctx.interrupt_flag());
        let Some(found) = computed(ev, e, outcome)? else {
            return Ok(None);
        };
        if let Some(r) = found {
            let result = ev.eval(&representation(&r, &vars))?;
            ev.step("GroebnerBasis", e, &result);
            return Ok(Some(result));
        }
    }
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

/// `{rur, l, f(t), f'(t), g₁(t), …}`, giac's form of `r` in `vars`: `t` is
/// `l` when `l` is a variable, else the first variable.
fn representation(r: &Rur, vars: &[Expr]) -> Expr {
    let form = r
        .form
        .iter()
        .zip(vars)
        .map(|(c, v)| Expr::from(c.clone()) * v.clone());
    let l = Expr::call(Symbol::PLUS, form);
    let mut nonzero = r.form.iter().zip(vars).filter(|(c, _)| !c.is_zero());
    let t = match (nonzero.next(), nonzero.next()) {
        (Some((c, v)), None) if *c == Rational::ONE => v.clone(),
        _ => vars[0].clone(),
    };
    let univariate = |coefficients: &[Rational]| {
        let terms = coefficients.iter().enumerate().map(|(k, c)| {
            let k = i64::try_from(k).expect("a normal set's size");
            Expr::from(c.clone()) * t.pow(k)
        });
        Expr::call(Symbol::PLUS, terms)
    };
    let mut elements = vec![
        Expr::from(Symbol::RATIONAL_UNIVARIATE_REPRESENTATION),
        l,
        univariate(&r.minimal),
        univariate(&r.derivative),
    ];
    elements.extend(r.numerators.iter().map(|g| univariate(g)));
    Expr::call(Symbol::LIST, elements)
}

/// The symbols of `e` that are no builtins, in canonical order.
fn user_symbols(e: &Expr) -> Vec<Expr> {
    e.symbols()
        .into_iter()
        .filter(|s| !s.is_builtin())
        .map(Expr::from)
        .collect()
}

/// `IdealMemberQ[p, {p₁, …}, {x₁, …}]`, with an order after the variables,
/// or `IdealMemberQ[p, {p₁, …}]` in the symbols of both; `p` may be a list.
pub(super) fn ideal_member_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, polys, vars, order) = match e.args() {
        [p, polys] => {
            let both = Expr::call(Symbol::LIST, [p.clone(), polys.clone()]);
            let vars = user_symbols(&both);
            let vars = (!vars.is_empty()).then_some(vars);
            (p, polys, vars, Some(MonomialOrder::DegRevLex))
        }
        [p, polys, vars] => (p, polys, variables(vars), Some(MonomialOrder::DegRevLex)),
        [p, polys, vars, order] => (p, polys, variables(vars), monomial_order(order)),
        _ => return Ok(None),
    };
    let (Some(vars), Some(order)) = (vars, order) else {
        return Ok(None);
    };
    let members = if is_list(p) {
        p.args().to_vec()
    } else {
        vec![p.clone()]
    };
    let mut candidates = Vec::with_capacity(members.len());
    for q in &members {
        let Some(d) = distributed(ev, q, &vars, order)? else {
            return Ok(None);
        };
        candidates.push(d);
    }
    let Some(polys) = list(ev, polys, &vars, order)? else {
        return Ok(None);
    };
    let Some(basis) = reduced_basis(ev, e, &polys, order)? else {
        return Ok(None);
    };
    let mut answers = Vec::with_capacity(candidates.len());
    for q in &candidates {
        let ctx = ev.context();
        let outcome = reduce_by(q, &basis, order, limits(ctx), ctx.interrupt_flag());
        let Some(remainder) = computed(ev, e, outcome)? else {
            return Ok(None);
        };
        answers.push(Expr::from(i64::from(remainder.is_zero())));
    }
    let result = if is_list(p) {
        Expr::call(Symbol::LIST, answers)
    } else {
        answers.swap_remove(0)
    };
    ev.step("IdealMemberQ", e, &result);
    Ok(Some(result))
}

/// `GroebnerVariables[{p₁, …}]`: the symbols in order of first appearance,
/// stably sorted by `∏ (deg_x pᵢ + 1)`, as giac orders them; `{}` for
/// anything but a list.
pub(super) fn groebner_variables(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [polys] = e.args() else {
        return Ok(None);
    };
    if !is_list(polys) {
        return Ok(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new())));
    }
    let mut vars: Vec<Expr> = Vec::new();
    for sub in polys.post_order() {
        if sub.as_symbol().is_some_and(|s| !s.is_builtin()) && !vars.contains(sub) {
            vars.push(sub.clone());
        }
    }
    let Some(polys) = list(ev, polys, &vars, MonomialOrder::DegRevLex)? else {
        return Ok(None);
    };
    let weight = |k: usize| {
        polys.iter().fold(1u128, |w, g| {
            let degree = g.terms().iter().map(|(m, _)| m[k]).max().unwrap_or(0);
            w.saturating_mul(u128::from(degree) + 1)
        })
    };
    let mut weighted: Vec<(u128, Expr)> = vars
        .into_iter()
        .enumerate()
        .map(|(k, v)| (weight(k), v))
        .collect();
    weighted.sort_by_key(|(w, _)| *w);
    let result = Expr::call(Symbol::LIST, weighted.into_iter().map(|(_, v)| v));
    ev.step("GroebnerVariables", e, &result);
    Ok(Some(result))
}

/// `Syzygies[{p₁, …, pₘ}]` or `Syzygies[{…}, {x₁, …}]`: `{S, V, A, X}`.
/// The `Aᵢ` are the symbols the reader makes of `A0`, `A1`, …, a bounded
/// set: one per generator of the largest call.
pub(super) fn syzygies(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (polys, vars) = match e.args() {
        [polys] => (polys, Some(user_symbols(polys))),
        [polys, vars] => (polys, variables(vars)),
        _ => return Ok(None),
    };
    let Some(vars) = vars else {
        return Ok(None);
    };
    let m = polys.args().len();
    if !is_list(polys) || m == 0 || vars.is_empty() {
        return Ok(None);
    }
    let a: Vec<Expr> = (0..m)
        .map(|i| Expr::from(Symbol::new(&format!("$A{i}"))))
        .collect();
    if vars.iter().any(|v| a.contains(v)) {
        return Ok(None);
    }
    let Some(polys) = list(ev, polys, &vars, MonomialOrder::DegRevLex)? else {
        return Ok(None);
    };
    let ctx = ev.context();
    let outcome = syzygies_of(&polys, limits(ctx), ctx.interrupt_flag());
    let Some(found) = computed(ev, e, outcome)? else {
        return Ok(None);
    };
    let all = [a.as_slice(), &vars].concat();
    let s = basis_list(ev, &found, &all)?;
    let lists = [all, a, vars].map(|l| Expr::call(Symbol::LIST, l));
    let result = ev.eval(&Expr::call(Symbol::LIST, [s].into_iter().chain(lists)))?;
    ev.step("Syzygies", e, &result);
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
