//! Xcas's calculus helpers: integration by parts step by step, the
//! difference of a value at two bounds, and the derivative of a function
//! object (CAS-574, CAS-577, CAS-578). Ports of giac's own code.
//!
//! - [`by_parts_dv`]: `IntegrateByPartsDV[f, v]`, giac's `ibpdv`
//!   (intg.cc `_ibpdv`). With `f = u v′` it gives `[u v, −u′ v]`: `u` is
//!   `f / v′` brought `Together` (giac's `normal`). `f` may be a pair
//!   `[F, r]` from an earlier step: `F` is added to the product term and
//!   `r` is integrated by parts. When `v′ = 0` (giac's `0` to finish) the
//!   answer is `F + ∫ r dx`. Optional arguments: the variable (giac's `x`
//!   when omitted), then bounds `a, b`, which evaluate each product term
//!   and the last integral at `b` minus at `a`.
//! - [`by_parts_u`]: `IntegrateByPartsU[f, u]`, giac's `ibpu` (misc.cc
//!   `_ibpu`). With `f = u v′` it gives `[u v, −u′ v]`, where `v = ∫ f / u
//!   dx`. A second argument of 0 or 1 finishes, as above.
//! - [`preval`]: `Preval[F, a, b]` or `Preval[F, a, b, x]`, giac's `preval`
//!   (misc.cc `_preval`, intg.cc `preval`): `F(b) − F(a)` in `x` (giac's
//!   `x` when omitted). giac always takes one-sided limits, from inside the
//!   interval. Ours substitutes, and takes such a limit only when a
//!   substituted value is indeterminate or infinite, since our `Limit` goes
//!   by series and leaves some calls as they are (an `Integrate` that stays,
//!   for one).
//! - [`function_diff`]: `FunctionDiff[g]`, giac's `function_diff` and
//!   `fonction_derivee` (derive.cc `_function_diff`): the derivative of a
//!   function object as a `Function`. `g` is applied to a private symbol.
//!   Builtins and `Function`s are called on it, `Plus`, `Times` and `Power`
//!   apply to each of their operands, and numbers stay as they are. Then it
//!   is differentiated. A derivative that keeps a `Derivative` (a function
//!   we do not know, giac's `D` left as it is) leaves the call as it is.
//!   The parameter is `x`, or `x1`, `x2`, … when the derivative already
//!   has an `x` (giac's parameter is a hidden `` ` x` ``). giac's `@@`
//!   powers of `function_diff` and its `@` composition have no form in our
//!   language.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: F, u, v, x and the bounds a, b"
)]

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::pattern::substitute;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's default variable `x`.
fn default_variable() -> Expr {
    Expr::from(Symbol::new("x"))
}

/// The parts of a by-parts call: `F` (0 when the first argument is not a
/// pair), the integrand left, the second argument, the variable and the
/// bounds.
struct Step<'e> {
    sum: Expr,
    integrand: Expr,
    second: &'e Expr,
    x: Expr,
    bounds: Option<[&'e Expr; 2]>,
}

/// The [`Step`] of `ibpdv` or `ibpu`: `[f, w]`, `[f, w, x]` or `[f, w, x,
/// a, b]`, `f` an integrand or a pair `[F, r]`, `x` a symbol.
fn step(e: &Expr) -> Option<Step<'_>> {
    let (first, second, x, bounds) = match e.args() {
        [f, w] => (f, w, default_variable(), None),
        [f, w, x] => (f, w, x.clone(), None),
        [f, w, x, a, b] => (f, w, x.clone(), Some([a, b])),
        _ => return None,
    };
    x.as_symbol()?;
    let (sum, integrand) = if first.head() == Some(&Expr::from(Symbol::LIST)) {
        let [sum, rest] = first.args() else {
            return None;
        };
        (sum.clone(), rest.clone())
    } else {
        (Expr::from(0), first.clone())
    };
    Some(Step {
        sum,
        integrand,
        second,
        x,
        bounds,
    })
}

/// `F + w` where `w` is `term`, or its difference at the bounds.
fn added(ev: &mut Eval<'_>, s: &Step<'_>, term: Expr) -> Result<Expr, EvalError> {
    let term = match s.bounds {
        Some([a, b]) => between(ev, &term, &s.x, a, b)?,
        None => term,
    };
    ev.eval(&(s.sum.clone() + term))
}

/// `[F + u v, −u′ v]` brought `Together`, the last step of both.
fn pair(ev: &mut Eval<'_>, s: &Step<'_>, u: &Expr, v: &Expr) -> Native {
    let product = ev.eval(&(u.clone() * v.clone()))?;
    let sum = added(ev, s, product)?;
    let dv = ev.eval(&Expr::call(Symbol::D, [v.clone(), s.x.clone()]))?;
    let rest = ev.eval(&Expr::call(Symbol::TOGETHER, [-(u.clone() * dv)]))?;
    Ok(Some(Expr::call(Symbol::LIST, [sum, rest])))
}

/// `F + ∫ r dx`, at the bounds when there are any.
fn finished(ev: &mut Eval<'_>, s: &Step<'_>) -> Native {
    let integral = ev.eval(&Expr::call(
        Symbol::INTEGRATE,
        [s.integrand.clone(), s.x.clone()],
    ))?;
    Ok(Some(added(ev, s, integral)?))
}

/// `IntegrateByPartsDV[f, v, …]`, giac's `ibpdv`.
pub(super) fn by_parts_dv(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(s) = step(e) else {
        return Ok(None);
    };
    let dv = ev.eval(&Expr::call(Symbol::D, [s.second.clone(), s.x.clone()]))?;
    if dv == Expr::from(0) {
        return finished(ev, &s);
    }
    let u = ev.eval(&Expr::call(Symbol::TOGETHER, [s.integrand.clone() / dv]))?;
    pair(ev, &s, s.second, &u)
}

/// `IntegrateByPartsU[f, u, …]`, giac's `ibpu`.
pub(super) fn by_parts_u(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(s) = step(e) else {
        return Ok(None);
    };
    if *s.second == Expr::from(0) || *s.second == Expr::from(1) {
        return finished(ev, &s);
    }
    let dv = ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [s.integrand.clone() / s.second.clone()],
    ))?;
    let v = ev.eval(&Expr::call(Symbol::INTEGRATE, [dv, s.x.clone()]))?;
    pair(ev, &s, &v, s.second)
}

/// Whether `e` has an indeterminate or infinite part.
fn singular(e: &Expr) -> bool {
    [
        Symbol::INDETERMINATE,
        Symbol::INFINITY,
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
    ]
    .into_iter()
    .any(|s| e.contains(&Expr::from(s)))
}

/// `f` at `x = c`; when that is singular, its limit from the side `side`
/// (1 from the right, -1 from the left).
fn value_at(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: Symbol,
    c: &Expr,
    side: i64,
) -> Result<Expr, EvalError> {
    let value = ev.eval(&substitute(f, &[(x, c.clone())]))?;
    if !singular(&value) || singular(c) {
        return Ok(value);
    }
    ev.eval(&Expr::call(
        Symbol::LIMIT,
        [f.clone(), Expr::from(x), c.clone(), Expr::from(side)],
    ))
}

/// `f(b) − f(a)` in `x`, giac's `preval`: the sides of its limits from
/// inside `[a, b]`, both sides when the order of `a` and `b` is unknown.
fn between(ev: &mut Eval<'_>, f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Result<Expr, EvalError> {
    let Some(x) = x.as_symbol() else {
        return Ok(f.clone());
    };
    let less = ev.eval(&Expr::call(Symbol::LESS, [a.clone(), b.clone()]))?;
    let (from_a, from_b) = match less.as_symbol() {
        Some(Symbol::TRUE) => (1, -1),
        Some(Symbol::FALSE) if a != b => (-1, 1),
        _ => (0, 0),
    };
    let fb = value_at(ev, f, x, b, from_b)?;
    let fa = value_at(ev, f, x, a, from_a)?;
    ev.eval(&(fb - fa))
}

/// `Preval[F, a, b]` or `Preval[F, a, b, x]`, giac's `preval`; giac's
/// `preval(F, x = a, b)` too.
pub(super) fn preval(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let equation = |s: &Expr| -> Option<(Expr, Expr)> {
        match (s.head()?.as_symbol()?, s.args()) {
            (Symbol::EQUATION | Symbol::EQUAL, [x, a]) => Some((x.clone(), a.clone())),
            _ => None,
        }
    };
    let (f, x, a, b) = match e.args() {
        [f, a, b] => match equation(a) {
            Some((x, a)) => {
                let b = equation(b).map_or_else(|| b.clone(), |(_, b)| b);
                (f, x, a, b)
            }
            None => (f, default_variable(), a.clone(), b.clone()),
        },
        [f, a, b, x] => (f, x.clone(), a.clone(), b.clone()),
        _ => return Ok(None),
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    Ok(Some(between(ev, f, &x, &a, &b)?))
}

/// `g` applied to `t`: builtins and functions called, `Plus`, `Times` and
/// `Power` through their operands, numbers kept.
fn applied(g: &Expr, t: &Expr) -> Expr {
    if g.is_number() {
        return g.clone();
    }
    match g.head().and_then(Expr::as_symbol) {
        Some(h @ (Symbol::PLUS | Symbol::TIMES | Symbol::POWER)) => {
            Expr::call(h, g.args().iter().map(|a| applied(a, t)))
        }
        _ => Expr::apply(g, std::slice::from_ref(t)),
    }
}

/// The names `x1`, `x2`, … tried for the parameter before a private one.
const MAX_NAMES: u64 = 1000;

/// Whether `e` has a `Derivative`, or a `D` that stays, in it.
fn has_derivative(e: &Expr) -> bool {
    e.contains(&Expr::from(Symbol::DERIVATIVE)) || e.contains(&Expr::from(Symbol::D))
}

/// `FunctionDiff[g]`, giac's `function_diff`.
pub(super) fn function_diff(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [g] = e.args() else {
        return Ok(None);
    };
    let t = Expr::from(Symbol::new("FunctionDiff$t"));
    if g.contains(&t) {
        return Ok(None);
    }
    let value = ev.eval(&applied(g, &t))?;
    let derivative = ev.eval(&Expr::call(Symbol::D, [value, t.clone()]))?;
    if has_derivative(&derivative) {
        return Ok(None);
    }
    let x = std::iter::once("x".to_owned())
        .chain((1..=MAX_NAMES).map(|k| format!("x{k}")))
        .map(|name| Symbol::new(&name))
        .find(|x| !derivative.contains(&Expr::from(*x)))
        .unwrap_or_else(|| Symbol::new("FunctionDiff$x"));
    let t = t.as_symbol().unwrap_or_else(|| Symbol::new("x"));
    let parameters = Expr::call(Symbol::LIST, [Expr::from(x)]);
    let body = substitute(&derivative, &[(t, Expr::from(x))]);
    let function = Expr::call(Symbol::FUNCTION, [parameters.clone(), body.clone()]);
    if maps(ev, &derivative, t)? {
        return Ok(Some(function));
    }
    // giac's `1/x` of a list is the list of inverses; ours is not, so a
    // list argument goes through `Map`.
    let x = Expr::from(x);
    let listed = Expr::call(
        Symbol::IF,
        [
            Expr::call(Symbol::LIST_Q, [x.clone()]),
            Expr::call(Symbol::MAP, [function, x]),
            body,
        ],
    );
    Ok(Some(Expr::call(Symbol::FUNCTION, [parameters, listed])))
}

/// Whether `f`, in `t`, at a list of two symbols is the list of its values
/// at each.
fn maps(ev: &mut Eval<'_>, f: &Expr, t: Symbol) -> Result<bool, EvalError> {
    let [p, q] = ["FunctionDiff$p", "FunctionDiff$q"].map(|s| Expr::from(Symbol::new(s)));
    let at = |c: Expr| substitute(f, &[(t, c)]);
    let together = ev.eval(&at(Expr::call(Symbol::LIST, [p.clone(), q.clone()])))?;
    let apart = ev.eval(&Expr::call(Symbol::LIST, [at(p), at(q)]))?;
    Ok(together == apart)
}

#[cfg(test)]
mod tests;
