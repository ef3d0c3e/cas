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
//! - `Solve` (CAS-459): the roots in `𝔽ₚ` of `l − r`, from its linear
//!   factors, each once, ascending by their symmetric residues, as modular
//!   values; giac's `solve((x^2-2) % 7, x)` is `[-3, 3]`. A difference
//!   that vanishes is solved by every element and leaves the call.
//!
//! - `PowerMod[a, k, p, q, x]` (CAS-460), giac's `powmod(a, k, p, q, x)`:
//!   `a^k` in `𝔽ₚ[x]/(q)` by repeated squaring, `k ≥ 0`, a polynomial with
//!   symmetric coefficients, as giac's `powmod(x+1,10,3,x^2+1,x)` is `-x`;
//!   `x` is `x` when left out.
//!
//! Another modulus, another variable count or coefficients that are not
//! integers leave the call.

use cas_num::algebra::{EuclideanDomain, Ring, WordPrimeField};
use cas_num::{Integer, Natural, Primality};
use cas_poly::{Dense, DenseRing};

use super::parts;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// A polynomial over `𝔽ₚ` in the symbol `x`.
pub(in crate::native) struct Univariate {
    pub(in crate::native) p: u64,
    pub(in crate::native) x: Option<Expr>,
    pub(in crate::native) f: Dense<u64>,
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
pub(in crate::native) fn univariate(e: &Expr) -> Option<Univariate> {
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

/// The roots of `u.f`, not 0, in `𝔽ₚ` with their multiplicities, ascending
/// by their symmetric residues; `None` when a factor of degree 2 or more
/// keeps it from splitting.
pub(in crate::native) fn split(u: &Univariate) -> Option<Vec<(Integer, usize)>> {
    let ring = DenseRing::new(WordPrimeField::new(u.p).expect("a prime"));
    let (_, factors) = ring.factor_mod(&u.f);
    let mut roots = Vec::with_capacity(factors.len());
    for (g, k) in &factors {
        match g.coeffs() {
            // A monic x + c: the root −c.
            [c, 1] => roots.push((
                symmetric((u.p - c) % u.p, u.p),
                usize::try_from(*k).expect("a multiplicity"),
            )),
            _ => return None,
        }
    }
    roots.sort();
    Some(roots)
}

/// The roots in `𝔽ₚ` of the linear factors of `u.f`, not 0, each once.
fn linear_roots(u: &Univariate) -> Vec<Integer> {
    let ring = DenseRing::new(WordPrimeField::new(u.p).expect("a prime"));
    let (_, factors) = ring.factor_mod(&u.f);
    let mut roots: Vec<Integer> = factors
        .iter()
        .filter_map(|(g, _)| match g.coeffs() {
            [c, 1] => Some(symmetric((u.p - c) % u.p, u.p)),
            _ => None,
        })
        .collect();
    roots.sort();
    roots
}

/// `Solve[eqn, x]` for an equation of modular values modulo a word prime,
/// in one unknown: `None` for another, `Some(None)` for one that stays.
pub(in crate::native) fn solve(
    ev: &mut Eval<'_>,
    e: &Expr,
    eqn: &Expr,
    x: &Expr,
) -> super::matrices::Hooked {
    if !eqn.contains(&Expr::from(Symbol::MODULAR)) {
        return Ok(None);
    }
    let difference = match (eqn.head().and_then(Expr::as_symbol), eqn.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => ev.eval(&Expr::call(
            Symbol::PLUS,
            [l.clone(), Expr::from(-1) * r.clone()],
        ))?,
        _ => eqn.clone(),
    };
    let Some(u) = univariate(&difference) else {
        return Ok(None);
    };
    if u.x.as_ref().is_some_and(|v| v != x) {
        return Ok(None);
    }
    if u.f.is_zero() {
        ev.warn(crate::trace::Warning::Native {
            message: format!("every element of the field of {} elements solves it", u.p).into(),
            expr: e.clone(),
        });
        return Ok(Some(None));
    }
    let roots: Vec<Expr> = linear_roots(&u)
        .into_iter()
        .map(|r| modular(Expr::from(r), u.p))
        .collect();
    let result = Expr::call(Symbol::LIST, roots);
    ev.step("Solve", e, &result);
    Ok(Some(Some(result)))
}

/// `PowerMod[a, k, p, q, x]` and `PowerMod[a, k, p, q]`.
pub(in crate::native) fn power_mod(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (a, k, p, q, x) = match e.args() {
        [a, k, p, q, x] => (a, k, p, q, x.clone()),
        [a, k, p, q] => (a, k, p, q, Expr::symbol("x")),
        _ => return Ok(None),
    };
    let (Some(mut k), Some(p)) = (k.to_integer(), p.to_integer()) else {
        return Ok(None);
    };
    if x.as_symbol().is_none() || k.is_negative() {
        return Ok(None);
    }
    let mut polynomial = |f: &Expr| -> Result<Option<Univariate>, EvalError> {
        let m = ev.eval(&Expr::call(
            Symbol::MODULAR,
            [f.clone(), Expr::from(p.clone())],
        ))?;
        Ok(univariate(&m).filter(|u| u.x.as_ref().is_none_or(|v| *v == x)))
    };
    let (Some(a), Some(q)) = (polynomial(a)?, polynomial(q)?) else {
        return Ok(None);
    };
    if q.f.is_zero() {
        return Ok(None);
    }
    let ring = DenseRing::new(WordPrimeField::new(a.p).expect("a prime"));
    let reduce = |f: &Dense<u64>| ring.div_rem(f, &q.f).1;
    let mut base = reduce(&a.f);
    let mut result = reduce(&ring.poly(vec![1]));
    let two = Integer::from(2);
    while !k.is_zero() {
        let (half, bit) = k.div_rem_euclid(&two);
        if !bit.is_zero() {
            result = reduce(&ring.mul(&result, &base));
        }
        k = half;
        if !k.is_zero() {
            base = reduce(&ring.mul(&base, &base));
        }
    }
    let answer = expression(&result, Some(&x), a.p);
    ev.step("powmod", e, &answer);
    Ok(Some(answer))
}
