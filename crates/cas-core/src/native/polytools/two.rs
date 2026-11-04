//! Polynomial tools II (CAS-DOC-48 §2, CAS-482): Xcas's `proot`, `pcoeff`,
//! `peval` and `horner`, `ptayl` and `froot`. A polynomial is given as an
//! expression in `x` or as giac's list of coefficients, the leading one
//! first.
//!
//! - **`PolynomialRoots[p]`**: every complex root of `p` over ℚ, each as
//!   often as its multiplicity, by the square-free parts: the roots of a
//!   linear part and of a quadratic part whose discriminant is a square or
//!   minus one exactly (giac's `proot([1,0,1])` is `[-i,i]`), the others
//!   balls by cas-num's certified isolation (CAS-211); ordered by real part,
//!   then imaginary part.
//! - **`FromRoots[l]`**: the coefficients of `Π (x − rᵢ)`, the leading one
//!   first, giac's `poly1[1,-6,11,-6]` for `[1,2,3]`.
//! - **`HornerEval[p, a]`**: `p(a)` by Horner's scheme.
//! - **`TaylorShift[p, a]`**: `p(x + a)`, a list for a list (giac's
//!   `[1,6,14,13]` for `[1,0,2,1]` at 2), expanded for an expression.
//! - **`RootMultiplicities[f]`**: the roots of the numerator and the
//!   denominator of `f` over ℚ, each followed by its multiplicity, negated
//!   for a pole: giac's `froot(x^2*(x-3)^3)` is `[0,2,3,3]`. Each
//!   irreducible factor's roots are `Solve`'s, rational ones first,
//!   ascending; a factor whose roots `Solve` leaves leaves the call.

use std::cmp::Ordering;

use cas_num::algebra::Rationals;
use cas_num::ball::{ComplexBall, Float, Round, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Dense, DenseRing};

use super::super::integrate::polynomial;
use super::{call, coefficients};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The coefficients of `p` from `x⁰` up: of a list (leading first), or of
/// an expression in `x`.
fn ascending(ev: &mut Eval<'_>, p: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    if is_list(p) {
        return Ok(Some(p.args().iter().rev().cloned().collect()));
    }
    coefficients(ev, p, &Expr::symbol("x"))
}

/// `Σ cₖ xᵏ` from ascending coefficients.
fn expression(cs: &[Expr], x: &Expr) -> Expr {
    Expr::call(
        Symbol::PLUS,
        cs.iter()
            .enumerate()
            .map(|(k, c)| c.clone() * x.pow(i64::try_from(k).expect("small")))
            .collect::<Vec<_>>(),
    )
}

/// `p` over ℚ as a dense polynomial.
fn dense(ev: &mut Eval<'_>, e: &Expr, p: &Expr) -> Result<Option<Dense<Rational>>, EvalError> {
    let x = Expr::symbol("x");
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    polynomial(ev, e, &expression(&cs, &x), &x)
}

/// The exact square root of a rational, when it has one.
fn rational_sqrt(q: &Rational) -> Option<Rational> {
    if q.is_negative() {
        return None;
    }
    let (n, d) = (q.numer().isqrt(), q.denom().isqrt());
    let r = Rational::new(n, d);
    (&r * &r == *q).then_some(r)
}

/// A root, exact or a ball, with its real and imaginary parts' midpoints
/// for the order.
struct Root {
    value: Expr,
    key: (f64, f64),
}

/// A rational's nearest double, for the order.
fn approximate(q: &Rational) -> f64 {
    Float::from_rational(q, 64, Round::Nearest).0.to_f64()
}

fn exact(z: Complex) -> Root {
    let key = (approximate(z.re()), approximate(z.im()));
    Root {
        value: Expr::from(z),
        key,
    }
}

/// The roots of a square-free part over ℚ, `None` when they are not
/// isolated.
fn part_roots(factor: &Dense<Rational>, bits: u64) -> Option<Vec<Root>> {
    let cs = factor.coeffs();
    match cs {
        [c, b] => return Some(vec![exact(Complex::from(-(c / b)))]),
        [c, b, a] => {
            // (−b ± √Δ)/(2a), exact when Δ or −Δ is a square.
            let delta = &(b * b) - &(&Rational::from(Integer::from(4)) * &(a * c));
            let two_a = a + a;
            let centre = -(b / &two_a);
            if let Some(r) = rational_sqrt(&delta) {
                let r = &r / &two_a;
                return Some(vec![
                    exact(Complex::from(&centre - &r)),
                    exact(Complex::from(&centre + &r)),
                ]);
            }
            if let Some(r) = rational_sqrt(&-delta) {
                let r = (&r / &two_a).abs();
                return Some(vec![
                    exact(Complex::new(centre.clone(), -r.clone())),
                    exact(Complex::new(centre, r)),
                ]);
            }
        }
        _ => {}
    }
    let complex: Vec<Complex> = cs.iter().cloned().map(Complex::from).collect();
    let roots = polynomial_roots(&complex, bits, 64 * bits.max(64))?;
    Some(
        roots
            .into_iter()
            .map(|z: ComplexBall| {
                let key = (z.re().mid().to_f64(), z.im().mid().to_f64());
                let value = if z.is_real() {
                    Expr::from(z.re().clone())
                } else {
                    Expr::from(z)
                };
                Root { value, key }
            })
            .collect(),
    )
}

/// `PolynomialRoots[p]`.
pub(in crate::native) fn polynomial_roots_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let Some(d) = dense(ev, e, p)? else {
        return Ok(None);
    };
    if d.degree().is_none_or(|n| n == 0) {
        return Ok(None);
    }
    let bits = ev.context().working_precision();
    let qx = DenseRing::new(Rationals);
    let mut roots: Vec<Root> = Vec::new();
    for (factor, m) in qx.squarefree(&d).factors() {
        if factor.degree().is_none_or(|n| n == 0) {
            continue;
        }
        let Some(found) = part_roots(factor, bits) else {
            return Ok(None);
        };
        for r in found {
            for _ in 0..*m {
                roots.push(Root {
                    value: r.value.clone(),
                    key: r.key,
                });
            }
        }
    }
    roots.sort_by(|a, b| {
        a.key
            .0
            .partial_cmp(&b.key.0)
            .unwrap_or(Ordering::Equal)
            .then(a.key.1.partial_cmp(&b.key.1).unwrap_or(Ordering::Equal))
    });
    let result = Expr::call(
        Symbol::LIST,
        roots.into_iter().map(|r| r.value).collect::<Vec<_>>(),
    );
    ev.step("proot", e, &result);
    Ok(Some(result))
}

/// `FromRoots[l]`.
pub(in crate::native) fn from_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let x = Expr::symbol("x");
    let product = Expr::call(
        Symbol::TIMES,
        l.args()
            .iter()
            .map(|r| x.clone() - r.clone())
            .collect::<Vec<_>>(),
    );
    let Some(cs) = coefficients(ev, &product, &x)? else {
        return Ok(None);
    };
    let result = Expr::call(Symbol::LIST, cs.into_iter().rev().collect::<Vec<_>>());
    ev.step("pcoeff", e, &result);
    Ok(Some(result))
}

/// `HornerEval[p, a]`.
pub(in crate::native) fn horner_eval(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, a] = e.args() else {
        return Ok(None);
    };
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    let value = cs
        .iter()
        .rev()
        .fold(Expr::from(0), |acc, c| acc * a.clone() + c.clone());
    let result = ev.eval(&value)?;
    ev.step("horner", e, &result);
    Ok(Some(result))
}

/// `TaylorShift[p, a]`.
pub(in crate::native) fn taylor_shift(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, a] = e.args() else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    let shifted = expression(&cs, &x).substitute(&[(x.clone(), x.clone() + a.clone())]);
    let expanded = ev.eval(&call(Symbol::EXPAND, vec![shifted]))?;
    let result = if is_list(p) {
        let Some(new) = coefficients(ev, &expanded, &x)? else {
            return Ok(None);
        };
        Expr::call(Symbol::LIST, new.into_iter().rev().collect::<Vec<_>>())
    } else {
        expanded
    };
    ev.step("ptayl", e, &result);
    Ok(Some(result))
}

/// `RootMultiplicities[f]`.
pub(in crate::native) fn root_multiplicities(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f] = e.args() else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let together = ev.eval(&call(Symbol::TOGETHER, vec![f.clone()]))?;
    let numerator = ev.eval(&call(Symbol::NUMERATOR, vec![together.clone()]))?;
    let denominator = ev.eval(&call(Symbol::DENOMINATOR, vec![together]))?;
    let mut pairs: Vec<(Expr, i64, Option<Rational>)> = Vec::new();
    for (part, sign) in [(numerator, 1), (denominator, -1)] {
        let factors = ev.eval(&call(Symbol::FACTOR_LIST, vec![part]))?;
        if !is_list(&factors) {
            return Ok(None);
        }
        for item in factors.args() {
            let [g, m] = item.args() else {
                return Ok(None);
            };
            if !g.contains(&x) {
                continue;
            }
            let Some(m) = m.as_i64() else {
                return Ok(None);
            };
            let roots = ev.eval(&call(
                Symbol::SOLVE,
                vec![
                    Expr::call(Symbol::EQUATION, [g.clone(), Expr::from(0)]),
                    x.clone(),
                ],
            ))?;
            if !is_list(&roots) {
                return Ok(None);
            }
            for r in roots.args() {
                pairs.push((r.clone(), sign * m, r.to_rational()));
            }
        }
    }
    // Rational roots first, ascending.
    pairs.sort_by(|a, b| match (&a.2, &b.2) {
        (Some(p), Some(q)) => p.cmp(q),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    });
    let result = Expr::call(
        Symbol::LIST,
        pairs
            .into_iter()
            .flat_map(|(r, m, _)| [r, Expr::from(m)])
            .collect::<Vec<_>>(),
    );
    ev.step("froot", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
