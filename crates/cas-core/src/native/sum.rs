//! `Sum` and `Product` (CAS-DOC-13 §1, §2, CAS-224): definite sums and
//! products in an index, Xcas's `sum` and `product`.
//!
//! - `Sum[f, k, a, b]` and `Sum[f, Equation[k, Range[a, b]]]`, giac's
//!   `sum(f, k, a, b)` and `sum(f, k=a..b)`; `Product` likewise.
//! - Karr's convention: for `b = a − 1` the empty sum 0 and product 1; for
//!   `b < a − 1`, `−Σ_{k=b+1}^{a−1}` and `1/Π_{k=b+1}^{a−1}`, as giac's
//!   `sum(k,k,5,1)` is `-9`.
//! - Integer bounds with at most [`MAX_TERMS`] terms: the terms, built for
//!   the loop to add or multiply (`sum(k,k,1,10)` is `55`).
//! - Otherwise a sum by linearity: a term free of `k` times the number of
//!   terms, a factor free of `k` out of a product, a sum term by term; a
//!   polynomial in `k` by its antidifference ([`antidifference`]); the
//!   terms without a closed form stay under one `Sum`, and when none has
//!   one the call stays, as giac's `sum(k!,k,0,n)`. A product of a factor
//!   free of `k` is its power; another stays (CAS-227).
//!
//! Each call that computes takes one step, `"sum"` or `"product"`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of sums: f in k from a to b, its antidifference P"
)]

use cas_num::algebra::{GaussianIntegers, Integers};
use cas_num::{Integer, Rational};
use cas_poly::{Sparse, SparseRing};

use super::polynomial::{Numbers, convert};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most terms a sum or product with integer bounds is written out with.
const MAX_TERMS: i64 = 10_000;

/// The highest degree of a polynomial summed by its antidifference.
const MAX_DEGREE: u64 = 1000;

/// The highest degree of a polynomial sum written factored, past which it
/// is expanded.
const MAX_FACTORED: u64 = 24;

/// The summand, the index and the bounds of `Sum` or `Product`.
fn parts(e: &Expr) -> Option<(&Expr, &Expr, &Expr, &Expr)> {
    let (f, k, a, b) = match e.args() {
        [f, k, a, b] => (f, k, a, b),
        [f, spec] => {
            let (Some(Symbol::EQUATION | Symbol::EQUAL), [k, range]) =
                (spec.head().and_then(Expr::as_symbol), spec.args())
            else {
                return None;
            };
            let (Some(Symbol::RANGE), [a, b]) =
                (range.head().and_then(Expr::as_symbol), range.args())
            else {
                return None;
            };
            (f, k, a, b)
        }
        _ => return None,
    };
    k.as_symbol().is_some().then_some((f, k, a, b))
}

/// Whether this is a sum or a product.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Sum,
    Product,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Sum => "sum",
            Kind::Product => "product",
        }
    }
}

/// `Sum[f, k, a, b]`.
pub(super) fn sum(ev: &mut Eval<'_>, e: &Expr) -> Native {
    definite(ev, e, Kind::Sum)
}

/// `Product[f, k, a, b]`.
pub(super) fn product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    definite(ev, e, Kind::Product)
}

/// [`sum`] or [`product`].
fn definite(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((f, k, a, b)) = parts(e) else {
        return Ok(None);
    };
    let result = match (a.to_rational(), b.to_rational()) {
        (Some(a), Some(b)) if a.is_integer() && b.is_integer() => {
            let (a, b) = (a.numer().clone(), b.numer().clone());
            match finite(f, k, &a, &b, kind) {
                Some(result) => result,
                None => match symbolic(ev, e, f, k, &Expr::from(a), &Expr::from(b), kind)? {
                    Some(result) => result,
                    None => return Ok(None),
                },
            }
        }
        _ => match symbolic(ev, e, f, k, a, b, kind)? {
            Some(result) => result,
            None => return Ok(None),
        },
    };
    ev.step(kind.name(), e, &result);
    Ok(Some(result))
}

/// The terms of a sum or product with integer bounds, by Karr's convention;
/// `None` past [`MAX_TERMS`].
fn finite(f: &Expr, k: &Expr, a: &Integer, b: &Integer, kind: Kind) -> Option<Expr> {
    let count = &(b - a) + &Integer::ONE;
    if count.is_negative() {
        // Σ_a^b = −Σ_{b+1}^{a−1}, Π_a^b = 1/Π_{b+1}^{a−1}.
        let reversed = finite(f, k, &(b + &Integer::ONE), &(a - &Integer::ONE), kind)?;
        return Some(match kind {
            Kind::Sum => -reversed,
            Kind::Product => reversed.pow(-1),
        });
    }
    let count = count.to_i64().filter(|&c| c <= MAX_TERMS)?;
    let terms: Vec<Expr> = (0..count)
        .map(|i| {
            let index = Expr::from(a + &Integer::from(i));
            f.substitute(&[(k.clone(), index)])
        })
        .collect();
    Some(match kind {
        Kind::Sum => Expr::call(Symbol::PLUS, terms),
        Kind::Product => Expr::call(Symbol::TIMES, terms),
    })
}

/// A sum or product with symbolic bounds, or too many terms: by linearity
/// and antidifferences; `None` when nothing has a closed form.
fn symbolic(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
    kind: Kind,
) -> Result<Option<Expr>, EvalError> {
    let count = b.clone() - a.clone() + Expr::from(1);
    if !f.contains(k) {
        return Ok(Some(match kind {
            Kind::Sum => f.clone() * count,
            Kind::Product => f.clone().pow(count),
        }));
    }
    if kind == Kind::Product {
        return Ok(None);
    }
    let terms: Vec<Expr> = match f.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => f.args().to_vec(),
        _ => vec![f.clone()],
    };
    let mut closed = Vec::new();
    let mut open = Vec::new();
    for term in terms {
        match term_sum(ev, e, &term, k, a, b)? {
            Some(s) => closed.push(s),
            None => open.push(term),
        }
    }
    if closed.is_empty() {
        return Ok(None);
    }
    if !open.is_empty() {
        let rest = Expr::call(Symbol::PLUS, open);
        closed.push(Expr::call(
            Symbol::SUM,
            [rest, k.clone(), a.clone(), b.clone()],
        ));
    }
    Ok(Some(Expr::call(Symbol::PLUS, closed)))
}

/// The sum of one term from `a` to `b`: a factor free of `k` out of a
/// product, a polynomial by its antidifference; `None` otherwise.
fn term_sum(
    ev: &mut Eval<'_>,
    e: &Expr,
    term: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    if !term.contains(k) {
        return Ok(Some(term.clone() * (b.clone() - a.clone() + Expr::from(1))));
    }
    if term.head().and_then(Expr::as_symbol) == Some(Symbol::TIMES) {
        let (free, bound): (Vec<Expr>, Vec<Expr>) =
            term.args().iter().cloned().partition(|x| !x.contains(k));
        if !free.is_empty() {
            let rest = Expr::call(Symbol::TIMES, bound);
            return Ok(term_sum(ev, e, &rest, k, a, b)?.map(|s| {
                Expr::call(
                    Symbol::TIMES,
                    free.into_iter().chain([s]).collect::<Vec<_>>(),
                )
            }));
        }
    }
    antidifference(ev, e, term, k, a, b)
}

/// The sum from `a` to `b` of `p`, a polynomial in `k` whose other
/// generators are free of it: `P(b + 1) − P(a)` for its antidifference
/// `P(k) = Σ cₘ Σⱼ S(m, j) k^(j+1) / (j + 1)`, `k^(j)` the falling
/// factorial and `S` Stirling's numbers of the second kind, factored up to
/// [`MAX_FACTORED`]; `None` when `p` is no such polynomial.
fn antidifference(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some(conversion) = convert(ev, e, std::slice::from_ref(p), true)? else {
        return Ok(None);
    };
    let Some(kv) = conversion.variable(k) else {
        return Ok(None);
    };
    if conversion
        .generators()
        .iter()
        .enumerate()
        .any(|(j, g)| j != kv && g.contains(k))
    {
        return Ok(None);
    }
    let n = conversion.generators().len();
    let coefficients = match &conversion.fractions()[0] {
        Fraction::Integers { num, den } => {
            coefficients(&conversion, &SparseRing::new(Integers, n), num, den, kv)
        }
        Fraction::GaussianIntegers { num, den } => coefficients(
            &conversion,
            &SparseRing::new(GaussianIntegers, n),
            num,
            den,
            kv,
        ),
    };
    let Some(coefficients) = coefficients else {
        return Ok(None);
    };
    let degree = u64::try_from(coefficients.len()).unwrap_or(u64::MAX);
    if degree > MAX_DEGREE + 1 {
        return Ok(None);
    }
    let stirling = stirling(coefficients.len());
    let at = |x: &Expr| -> Expr {
        Expr::call(
            Symbol::PLUS,
            coefficients
                .iter()
                .enumerate()
                .filter(|(_, c)| **c != Expr::from(0))
                .map(|(m, c)| {
                    let antiderivative: Vec<Expr> = (0..=m)
                        .filter(|&j| !stirling[m][j].is_zero())
                        .map(|j| {
                            let s = Rational::new(stirling[m][j].clone(), Integer::from(j + 1));
                            Expr::from(s) * falling(x, j + 1)
                        })
                        .collect();
                    c.clone() * Expr::call(Symbol::PLUS, antiderivative)
                })
                .collect::<Vec<_>>(),
        )
    };
    let upper = b.clone() + Expr::from(1);
    let total = ev.eval(&(at(&upper) - at(a)))?;
    let head = if degree <= MAX_FACTORED {
        Symbol::FACTOR
    } else {
        Symbol::EXPAND
    };
    Ok(Some(ev.eval(&Expr::call(head, [total]))?))
}

/// The coefficients of `num / den` in the variable `kv`, from the constant,
/// built; `None` when the denominator has `kv`.
fn coefficients<R: Numbers>(
    p: &Polynomials,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
    kv: usize,
) -> Option<Vec<Expr>> {
    if ring.degree(den, kv).unwrap_or(0) > 0 {
        return None;
    }
    let degree = ring.degree(num, kv).unwrap_or(0);
    if degree > MAX_DEGREE {
        return None;
    }
    Some(
        (0..=degree)
            .map(|d| p.ratio_expr(&ring.coeff(num, kv, d), den))
            .collect(),
    )
}

/// The falling factorial `x (x − 1) ⋯ (x − j + 1)`, built.
fn falling(x: &Expr, j: usize) -> Expr {
    Expr::call(
        Symbol::TIMES,
        (0..j)
            .map(|i| x.clone() - Expr::from(i64::try_from(i).expect("a small degree")))
            .collect::<Vec<_>>(),
    )
}

/// Stirling's numbers of the second kind `S(m, j)` for `m, j < n`:
/// `S(m, j) = j S(m − 1, j) + S(m − 1, j − 1)`, `S(0, 0) = 1`.
fn stirling(n: usize) -> Vec<Vec<Integer>> {
    let mut s = vec![vec![Integer::ZERO; n.max(1)]; n.max(1)];
    s[0][0] = Integer::ONE;
    for m in 1..n {
        for j in 1..=m {
            s[m][j] = &(&Integer::from(j) * &s[m - 1][j]) + &s[m - 1][j - 1];
        }
    }
    s
}

#[cfg(test)]
mod tests;
