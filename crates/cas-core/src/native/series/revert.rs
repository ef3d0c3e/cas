//! `Revert[f]`, Xcas's `revert` (CAS-579): the inverse of a series for
//! composition. A port of giac's series.cc `_revert` and `prevert`.
//!
//! - `Revert[f]` in `x` (giac's `x` when no variable is given),
//!   `Revert[f, x]` and `Revert[f, x, n]`. A polynomial is taken as it is,
//!   to its degree (giac's `revert(x+x^2,x,8)` is `x-x^2`); any other `f`
//!   by its series to `x^n`, `n` giac's default order 5 when it is not an
//!   integer or not given. The exponent past the known terms is `N`.
//! - The constant term is dropped, as prevert drops it. The first term
//!   left, `a x^k` with `k ≥ 1`, gives the inverse
//!   `Σ_{m=1}^{kN−1} g_m y^(m/k)` in the same variable, by Lagrange's
//!   inversion of `φ = P^(1/k) = b x (1 + h)^(1/k)`, `b = a^(1/k)`:
//!   `g_m = (1/m) b^(−m) [x^(m−1)] (1 + h)^(−m/k)`, each power by
//!   J. C. P. Miller's recurrence. These are prevert's coefficients, which
//!   giac finds by composing term by term; `k > 1` gives a Puiseux series
//!   (`revert(x^2+x^3)` starts `sqrt(x)`). No order term is written, as in
//!   giac.
//! - The cost is about `(kN)³` coefficient operations: past the
//!   evaluator's `terms` limit the call stays (C0019). A series with a pole
//!   or no term past the constant has no inverse here, as in giac.

use cas_num::{Integer, Rational};

use super::{expansion, fresh};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `series_default_order`.
const DEFAULT_ORDER: i64 = 5;

/// Whether `e` is a polynomial in `t`: `t` only through sums, products and
/// powers by nonnegative integers.
fn polynomial_in(e: &Expr, t: &Expr) -> bool {
    if e == t || !e.contains(t) {
        return true;
    }
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), args) => args.iter().all(|a| polynomial_in(a, t)),
        (Some(Symbol::POWER), [b, n]) => n.as_i64().is_some_and(|n| n >= 0) && polynomial_in(b, t),
        _ => false,
    }
}

fn q(n: i64, d: i64) -> Expr {
    Expr::from(Rational::new(Integer::from(n), Integer::from(d)))
}

pub(crate) fn revert(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let order = |n: &Expr| n.as_i64().unwrap_or(DEFAULT_ORDER);
    let (f, x, order) = match e.args() {
        [f] => (f, Expr::from(Symbol::new("x")), DEFAULT_ORDER),
        [f, x] => (f, x.clone(), DEFAULT_ORDER),
        [f, x, n] => (f, x.clone(), order(n)),
        _ => return Ok(None),
    };
    let t = fresh();
    if x.as_symbol().is_none() || order < 1 || f.contains(&t) {
        return Ok(None);
    }
    let local = ev.eval(&f.substitute(&[(x.clone(), t.clone())]))?;
    let known = if polynomial_in(&local, &t) {
        let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [local.clone(), t.clone()]))?;
        match degree.as_i64() {
            Some(d) if d >= 1 => d + 1,
            _ => return Ok(None),
        }
    } else {
        order + 1
    };
    let limit = u128::try_from(ev.context().limits().terms).unwrap_or(u128::MAX);
    let cube = |n: i64| u128::from(n.unsigned_abs()).pow(3);
    if cube(known) > limit {
        return Ok(None);
    }
    let Some(s) = expansion(ev, &local, &t, &Rational::from(Integer::from(known)))? else {
        return Ok(None);
    };
    if s.ramification != 1 || s.valuation < 0 {
        return Ok(None);
    }
    let zero = Expr::from(0);
    let a: Vec<Expr> = (0..known)
        .map(|j| s.at(j).cloned().unwrap_or_else(|| zero.clone()))
        .collect();
    let Some(k) = (1..known).find(|&j| a[usize::try_from(j).expect("small")] != zero) else {
        return Ok(None);
    };
    let terms = k * known - 1;
    if cube(terms) > limit {
        return Ok(None);
    }
    let lead = a[usize::try_from(k).expect("small")].clone();
    // h_i = a_{k+i}/a_k, 0 past the known terms; h[0] unused.
    let h: Vec<Expr> = (0..terms)
        .map(|i| {
            usize::try_from(k + i)
                .ok()
                .and_then(|j| a.get(j))
                .filter(|_| i > 0)
                .map_or_else(|| zero.clone(), |c| c.clone() / lead.clone())
        })
        .collect();
    let h = ev.eval(&Expr::call(Symbol::TOGETHER, [Expr::call(Symbol::LIST, h)]))?;
    let h = h.args().to_vec();
    let b = lead.pow(q(1, k));
    let mut out = Vec::new();
    for m in 1..=terms {
        // [x^(m−1)] (1 + h)^α, α = −m/k: p_j = (1/j) Σ ((α + 1) i − j) h_i p_{j−i}.
        let mut p = vec![Expr::from(1)];
        for j in 1..m {
            let sum: Expr = (1..=j)
                .filter(|&i| h[usize::try_from(i).expect("small")] != zero)
                .map(|i| {
                    let c = q((k - m) * i - j * k, k * j);
                    c * h[usize::try_from(i).expect("small")].clone()
                        * p[usize::try_from(j - i).expect("small")].clone()
                })
                .sum();
            p.push(ev.eval(&Expr::call(Symbol::TOGETHER, [sum]))?);
        }
        let g = p.pop().expect("one coefficient at least") / (Expr::from(m) * b.clone().pow(m));
        out.push(g * x.clone().pow(q(m, k)));
    }
    let result = ev.eval(&out.into_iter().sum::<Expr>())?;
    ev.step("Revert", e, &result);
    Ok(Some(result))
}
