//! Polynomials over a prime field (CAS-DOC-43 §5, CAS-453): `Factor`,
//! `GCD`, `PolynomialQuotient`, `PolynomialRemainder` and
//! `PolynomialQuotientRemainder` of `Modular[f, p]` in one variable, `p` a
//! prime that fits a word.
//!
//! - `Factor`: cas-poly's factoring over `𝔽ₚ` (square-free parts, then
//!   the distinct-degree factorization and the equal-degree splitting),
//!   the factors with integer coefficients, the symmetric residues, as
//!   Mathematica's `Factor[f, Modulus -> p]`: a product of modular values
//!   would expand again into one. giac's `factor((x^4+1) % 3)` is
//!   `(x^2+x-1)*(x^2-x-1) % 3`, here `(x^2+x-1)(x^2-x-1)`.
//! - `GCD`: monic, giac's `gcd((x^2-1) % 5,(x^2+3*x+2) % 5)` being `(x+1) %
//!   5`; the quotients and remainders by Euclid, modular values.
//!
//! Another modulus, another variable count or coefficients that are not
//! integers leave the call.

use cas_num::algebra::{EuclideanDomain, WordPrimeField};
use cas_num::{Integer, Natural, Primality};
use cas_poly::{Dense, DenseRing};

use super::parts;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// A polynomial over `𝔽ₚ` in the symbol `x`.
struct Univariate {
    p: u64,
    x: Option<Expr>,
    f: Dense<u64>,
}

/// The integer coefficient and the power of `x` of a term `c xᵏ`.
fn term(t: &Expr, x: &mut Option<Expr>) -> Option<(Integer, usize)> {
    if let Some(c) = t.to_integer() {
        return Some((c, 0));
    }
    let (c, rest): (Integer, &[Expr]) = match (t.head().and_then(Expr::as_symbol), t.args()) {
        (Some(Symbol::TIMES), [c, rest @ ..]) if c.to_integer().is_some() => {
            (c.to_integer()?, rest)
        }
        _ => (Integer::ONE, std::slice::from_ref(t)),
    };
    let [power] = rest else {
        return None;
    };
    let (base, k) = match (power.head().and_then(Expr::as_symbol), power.args()) {
        (Some(Symbol::POWER), [b, k]) => (b, usize::try_from(k.as_i64()?).ok()?),
        _ => (power, 1),
    };
    if base.as_symbol().is_none_or(Symbol::is_builtin) {
        return None;
    }
    match x {
        Some(x) if x != base => return None,
        _ => *x = Some(base.clone()),
    }
    Some((c, k))
}

/// `e`, a `Modular[f, p]`, as a polynomial over `𝔽ₚ`.
fn univariate(e: &Expr) -> Option<Univariate> {
    let (f, n) = parts(e)?;
    let p = u64::try_from(&n).ok()?;
    if Natural::from(p).primality() != Primality::Prime {
        return None;
    }
    let terms: &[Expr] = match f.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => f.args(),
        _ => std::slice::from_ref(f),
    };
    let mut x = None;
    let mut coeffs: Vec<u64> = Vec::new();
    for t in terms {
        let (c, k) = term(t, &mut x)?;
        if coeffs.len() <= k {
            coeffs.resize(k + 1, 0);
        }
        let (_, r) = c.div_rem_euclid(&n);
        coeffs[k] = (coeffs[k] + u64::try_from(&r).ok()?) % p;
    }
    let ring = DenseRing::new(WordPrimeField::new(p)?);
    Some(Univariate {
        p,
        x,
        f: ring.poly(coeffs),
    })
}

/// The symmetric residue of `c` modulo `p`.
fn symmetric(c: u64, p: u64) -> Integer {
    if 2 * u128::from(c) > u128::from(p) {
        Integer::from(c) - Integer::from(p)
    } else {
        Integer::from(c)
    }
}

/// `g` as an expression in `x`, its coefficients symmetric residues.
fn expression(g: &Dense<u64>, x: Option<&Expr>, p: u64) -> Expr {
    let x = x.cloned().unwrap_or_else(|| Expr::symbol("x"));
    Expr::call(
        Symbol::PLUS,
        g.coeffs()
            .iter()
            .enumerate()
            .filter(|(_, c)| **c != 0)
            .map(|(k, c)| {
                Expr::from(symmetric(*c, p)) * x.pow(i64::try_from(k).expect("a degree"))
            }),
    )
}

fn modular(g: Expr, p: u64) -> Expr {
    Expr::call(Symbol::MODULAR, [g, Expr::from(Integer::from(p))])
}

/// `Factor[Modular[f, p]]`.
pub(in crate::native) fn factor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    let Some(u) = univariate(arg) else {
        return Ok(None);
    };
    if u.f.is_zero() {
        return Ok(Some(Expr::from(0)));
    }
    let ring = DenseRing::new(WordPrimeField::new(u.p).expect("a prime"));
    let (lead, factors) = ring.factor_mod(&u.f);
    let mut product = vec![Expr::from(symmetric(lead, u.p))];
    for (g, k) in &factors {
        product.push(expression(g, u.x.as_ref(), u.p).pow(i64::try_from(*k).expect("small")));
    }
    let result = Expr::call(Symbol::TIMES, product);
    ev.step("Factor", e, &result);
    Ok(Some(result))
}

/// Two polynomials over one prime field, in one variable.
type Pair = (u64, Option<Expr>, Dense<u64>, Dense<u64>);

/// The polynomials of `Modular` arguments, one prime and one variable.
fn pair(args: &[Expr]) -> Option<Pair> {
    let [a, b, ..] = args else {
        return None;
    };
    let (a, b) = (univariate(a)?, univariate(b)?);
    if a.p != b.p || (a.x.is_some() && b.x.is_some() && a.x != b.x) {
        return None;
    }
    Some((a.p, a.x.or(b.x), a.f, b.f))
}

/// `GCD[Modular[f, p], Modular[g, p]]`, monic.
pub(in crate::native) fn gcd(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, x, f, g)) = pair(e.args()) else {
        return Ok(None);
    };
    if e.args().len() != 2 {
        return Ok(None);
    }
    let ring = DenseRing::new(WordPrimeField::new(p).expect("a prime"));
    let (h, _, _) = ring.gcd_cofactors(&f, &g);
    let result = modular(expression(&h, x.as_ref(), p), p);
    ev.step("gcd", e, &result);
    Ok(Some(result))
}

/// What a division gives.
#[derive(Clone, Copy)]
pub(in crate::native) enum Part {
    Quotient,
    Remainder,
    Both,
}

/// `PolynomialQuotient`, `PolynomialRemainder` and
/// `PolynomialQuotientRemainder` of modular polynomials.
pub(in crate::native) fn divide(ev: &mut Eval<'_>, e: &Expr, part: Part) -> Native {
    if !(2..=3).contains(&e.args().len()) {
        return Ok(None);
    }
    let Some((p, x, f, g)) = pair(e.args()) else {
        return Ok(None);
    };
    if g.is_zero() {
        return Ok(None);
    }
    let ring = DenseRing::new(WordPrimeField::new(p).expect("a prime"));
    let (q, r) = ring.div_rem(&f, &g);
    let wrap = |h: &Dense<u64>| modular(expression(h, x.as_ref(), p), p);
    let result = match part {
        Part::Quotient => wrap(&q),
        Part::Remainder => wrap(&r),
        Part::Both => Expr::call(Symbol::LIST, [wrap(&q), wrap(&r)]),
    };
    ev.step("quorem", e, &result);
    Ok(Some(result))
}
