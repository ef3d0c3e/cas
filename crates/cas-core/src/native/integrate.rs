//! `Integrate` (CAS-DOC-9, CAS-193): antiderivatives and definite integrals
//! in one symbol, Xcas's `integrate` and `int`, without a constant of
//! integration.
//!
//! - Lists thread over the integrand (CAS-194).
//! - What is free of the variable `x` is `f x`.
//! - Sums go term by term, and the terms that do not integrate stay under
//!   one `Integrate`, as giac's `x + f(x)` gives `x^2/2 + integrate(f(x), x)`.
//! - The factors of a product free of `x` come out of its antiderivative;
//!   when the rest does not integrate, the product stays whole, as giac's
//!   `integrate(2*f(x),x)`.
//! - The table at linear arguments `u = a x + b`, `a` and `b` free of `x`
//!   (CAS-DOC-9 §4, CAS-195), each entry over `a`: `u^n` is
//!   `u^(n+1) / (n+1)` for `n` free of `x`, giac's answer even for a
//!   symbolic `n`, and `u^-1` is `Log[Abs[u]]`; `Exp[u]` and `E^u` are
//!   themselves, `c^u` is `c^u / Log[c]`; `Sin`, `Cos`, `Tan` (`−ln|cos u|`)
//!   and `Log` (`u ln u − u`); `1/√(α + β x²)` with `α > 0 > β` is an
//!   arcsine, `1/(α + β x²)` with `α, β > 0` and `1/(1 + u²)` arctangents.
//! - Polynomials in `x`, through the conversion of CAS-DOC-8 §2 (CAS-194):
//!   each term `c x^k` is `c x^(k+1) / (k+1)`, over the common denominator.
//! - Rational functions of `x` with numeric coefficients (CAS-196), through
//!   cas-poly: the polynomial part, Hermite's rational part, and the
//!   logarithms with rational residues, `ln|S|` where `S` has a real root
//!   and `ln S` where it has none; the part whose residues are not rational
//!   stays under `Integrate`.
//!
//! `Integrate[f, x, a, b]`, and Mathematica's `Integrate[f, {x, a, b}]`, is
//! `F(b) − F(a)` where the antiderivative `F` is a polynomial in `x`, which
//! is continuous; other definite integrals wait for CAS-198. Each call that
//! rewrites takes one step, `"Integrate"`, whose identity `D[F, x] = f`
//! `Trace::verify` checks.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the calculus: f in x from a to b, n and m exponents"
)]

use cas_num::algebra::{GaussianIntegers, Integers, Rationals};
use cas_num::{Integer, Rational};
use cas_poly::{Cancel, Cancelled, Dense, DenseRing, RationalIntegral, Sparse, SparseRing};

use super::polynomial::{Numbers, convert};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The integrand, the variable, and the bounds of a definite integral.
type Parts<'e> = (&'e Expr, &'e Expr, Option<[&'e Expr; 2]>);

/// The parts of `Integrate[f, x]`, `Integrate[f, x, a, b]` or
/// `Integrate[f, {x, a, b}]`, the variable a symbol.
fn parts(e: &Expr) -> Option<Parts<'_>> {
    let (f, x, bounds) = match e.args() {
        [f, x] if x.head() == Some(&Expr::from(Symbol::LIST)) => match x.args() {
            [x, a, b] => (f, x, Some([a, b])),
            _ => return None,
        },
        [f, x] => (f, x, None),
        [f, x, a, b] => (f, x, Some([a, b])),
        _ => return None,
    };
    x.as_symbol().is_some().then_some((f, x, bounds))
}

pub(super) fn integrate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x, bounds)) = parts(e) else {
        return Ok(None);
    };
    let result = if f.head() == Some(&Expr::from(Symbol::LIST)) {
        // Lists thread over the integrand alone.
        let mut args = e.args().to_vec();
        Expr::call(
            Symbol::LIST,
            f.args().iter().map(|fi| {
                args[0] = fi.clone();
                Expr::call(Symbol::INTEGRATE, args.clone())
            }),
        )
    } else {
        let Some(antiderivative) = integral(ev, e, f, x)? else {
            return Ok(None);
        };
        match bounds {
            None => antiderivative,
            Some([a, b]) => {
                // Evaluated, as the conversion reads it: x^(2 + 1) is x^3.
                let antiderivative = ev.eval(&antiderivative)?;
                if !polynomial_in(ev, &antiderivative, x)? {
                    return Ok(None);
                }
                let at = |v: &Expr| antiderivative.substitute(&[(x.clone(), v.clone())]);
                at(b) - at(a)
            }
        }
    };
    ev.step("Integrate", e, &result);
    Ok(Some(result))
}

/// An antiderivative of `f` in `x`, the parts that do not integrate kept
/// under `Integrate`; `None` when no part integrates.
fn integral(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr) -> Result<Option<Expr>, EvalError> {
    if !f.contains(x) {
        return Ok(Some(f.clone() * x.clone()));
    }
    if f == x {
        return Ok(Some(x.pow(2) / 2));
    }
    match f.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => {
            let mut done = Vec::new();
            let mut rest = Vec::new();
            for t in f.args() {
                match integral(ev, e, t, x)? {
                    Some(g) => done.push(g),
                    None => rest.push(t.clone()),
                }
            }
            if done.is_empty() {
                return Ok(None);
            }
            if !rest.is_empty() {
                let rest: Expr = rest.into_iter().sum();
                done.push(Expr::call(Symbol::INTEGRATE, [rest, x.clone()]));
            }
            Ok(Some(done.into_iter().sum()))
        }
        Some(Symbol::TIMES) => {
            let (constant, varying): (Vec<Expr>, Vec<Expr>) =
                f.args().iter().cloned().partition(|a| !a.contains(x));
            if constant.is_empty() {
                return rational(ev, e, f, x);
            }
            // When the rest does not integrate, the product stays whole, as
            // giac leaves integrate(2*f(x),x).
            let varying: Expr = varying.into_iter().product();
            Ok(integral(ev, e, &varying, x)?.map(|g| constant.into_iter().product::<Expr>() * g))
        }
        _ => match table(f, x) {
            Some(g) => Ok(Some(g)),
            None => rational(ev, e, f, x),
        },
    }
}

/// The coefficient `a` of `u = a x + b`, `a` and `b` free of `x`; `None`
/// when `u` is not linear in `x`.
fn slope(u: &Expr, x: &Expr) -> Option<Expr> {
    // The coefficient of `x` in a term `c x`.
    let coefficient = |t: &Expr| -> Option<Expr> {
        if t == x {
            return Some(Expr::from(1));
        }
        if t.head() != Some(&Expr::from(Symbol::TIMES)) {
            return None;
        }
        let mut rest = Vec::new();
        let mut found = false;
        for factor in t.args() {
            if factor == x && !found {
                found = true;
            } else if factor.contains(x) {
                return None;
            } else {
                rest.push(factor.clone());
            }
        }
        found.then(|| rest.into_iter().product())
    };
    if let Some(a) = coefficient(u) {
        return Some(a);
    }
    if u.head() != Some(&Expr::from(Symbol::PLUS)) {
        return None;
    }
    let mut a = Vec::new();
    for t in u.args() {
        if t.contains(x) {
            a.push(coefficient(t)?);
        }
    }
    (!a.is_empty()).then(|| a.into_iter().sum())
}

fn call(head: Symbol, arg: Expr) -> Expr {
    Expr::call(head, [arg])
}

/// `f` from the table of antiderivatives at linear arguments `u = a x + b`
/// (CAS-DOC-9 §4), each entry over `a`.
fn table(f: &Expr, x: &Expr) -> Option<Expr> {
    let log_abs = |u: Expr| call(Symbol::LOG, call(Symbol::ABS, u));
    match (f.head()?.as_symbol()?, f.args()) {
        (Symbol::POWER, [u, n]) if !n.contains(x) => {
            let Some(a) = slope(u, x) else {
                return quadratic(u, n, x);
            };
            if n.as_i64() == Some(-1) {
                return Some(log_abs(u.clone()) / a);
            }
            let m = n.clone() + 1;
            Some(u.pow(m.clone()) / (m * a))
        }
        (Symbol::POWER, [c, u]) if !c.contains(x) => {
            let a = slope(u, x)?;
            // Log[E] stays as it is: E^u is the exponential's case.
            if c.as_symbol() == Some(Symbol::E) {
                return Some(f.clone() / a);
            }
            Some(f.clone() / (a * call(Symbol::LOG, c.clone())))
        }
        (Symbol::EXP, [u]) => Some(f.clone() / slope(u, x)?),
        (Symbol::SIN, [u]) => Some(-call(Symbol::COS, u.clone()) / slope(u, x)?),
        (Symbol::COS, [u]) => Some(call(Symbol::SIN, u.clone()) / slope(u, x)?),
        (Symbol::TAN, [u]) => {
            let a = slope(u, x)?;
            Some(-log_abs(call(Symbol::COS, u.clone())) / a)
        }
        (Symbol::LOG, [u]) => {
            let a = slope(u, x)?;
            Some((u.clone() * f.clone() - u.clone()) / a)
        }
        _ => None,
    }
}

/// `q^n` with `q = α + β x²`, `α` and `β` rationals, as an arcsine for
/// `n = −1/2` and `α > 0 > β`, `asin(x √(−β/α)) / √(−β)`, or an
/// arctangent for `n = −1` and `α, β > 0`, `atan(x √(β/α)) / √(α β)`; and
/// `(1 + u²)^−1` for a linear `u`, `atan(u) / a`.
fn quadratic(q: &Expr, n: &Expr, x: &Expr) -> Option<Expr> {
    let n = n.to_rational()?;
    let minus_one = Rational::from(Integer::from(-1));
    let minus_half = Rational::new(Integer::from(-1), Integer::from(2));
    let sqrt =
        |r: Rational| Expr::from(r).pow(Expr::from(Rational::new(Integer::ONE, Integer::from(2))));
    if n == minus_one
        && let [one, square] = q.args()
        && q.head() == Some(&Expr::from(Symbol::PLUS))
        && one.as_i64() == Some(1)
        && square.head() == Some(&Expr::from(Symbol::POWER))
        && square.args()[1].as_i64() == Some(2)
        && let Some(a) = slope(&square.args()[0], x)
    {
        return Some(call(Symbol::ARC_TAN, square.args()[0].clone()) / a);
    }
    let (alpha, beta) = even_quadratic(q, x)?;
    if n == minus_half && alpha.is_positive() && beta.is_negative() {
        let scale = sqrt(&-&beta / &alpha);
        return Some(call(Symbol::ARC_SIN, x.clone() * scale) / sqrt(-beta));
    }
    if n == minus_one && alpha.is_positive() && beta.is_positive() {
        let scale = sqrt(&beta / &alpha);
        return Some(call(Symbol::ARC_TAN, x.clone() * scale) / sqrt(&alpha * &beta));
    }
    None
}

/// `q` as `α + β x²`, `α` and `β` rationals, `β` not 0.
fn even_quadratic(q: &Expr, x: &Expr) -> Option<(Rational, Rational)> {
    if q.head() != Some(&Expr::from(Symbol::PLUS)) {
        return None;
    }
    let square = x.pow(2);
    let (mut alpha, mut beta) = (None, None);
    for t in q.args() {
        if let Some(r) = t.to_rational() {
            alpha = Some(r);
        } else if *t == square {
            beta = Some(Rational::ONE);
        } else if let [c, s] = t.args()
            && t.head() == Some(&Expr::from(Symbol::TIMES))
            && *s == square
            && let Some(r) = c.to_rational()
        {
            beta = Some(r);
        } else {
            return None;
        }
    }
    Some((alpha?, beta?))
}

/// `f` as a rational function whose generators but `x` are free of `x`,
/// with `x`'s index; `None` when it is not one.
fn in_x(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
) -> Result<Option<(Polynomials, usize)>, EvalError> {
    let Some(p) = convert(ev, e, std::slice::from_ref(f), false)? else {
        return Ok(None);
    };
    let Some(k) = p.variable(x) else {
        return Ok(None);
    };
    if p.generators()
        .iter()
        .enumerate()
        .any(|(j, g)| j != k && g.contains(x))
    {
        return Ok(None);
    }
    Ok(Some((p, k)))
}

/// The antiderivative of `f` as a rational function of `x` (CAS-DOC-9 §3):
/// over a denominator free of `x`, term by term, `c x^k` to
/// `c x^(k+1) / (k+1)` (CAS-194); otherwise, with `x` the only generator,
/// by cas-poly's integration of rational functions over ℚ (CAS-196), the
/// part whose residues are not rational kept under `Integrate`. `None` when
/// nothing integrates.
fn rational(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some((p, k)) = in_x(ev, e, f, x)? else {
        return Ok(None);
    };
    let n = p.generators().len();
    let (num, den) = match &p.fractions()[0] {
        Fraction::Integers { num, den } => {
            let ring = SparseRing::new(Integers, n);
            if let Some(g) = terms_integrated(&p, &ring, num, den, k) {
                return Ok(Some(g));
            }
            (num, den)
        }
        Fraction::GaussianIntegers { num, den } => {
            let ring = SparseRing::new(GaussianIntegers, n);
            return Ok(terms_integrated(&p, &ring, num, den, k));
        }
    };
    // Coefficients with other generators wait for CAS-DOC-9 §3.4.
    if n != 1 {
        return Ok(None);
    }
    let ring = SparseRing::new(Integers, 1);
    let qx = DenseRing::new(Rationals);
    let to_q = |s: &Sparse<Integer>| {
        let dense = ring.to_dense(s, 0).expect("one variable");
        qx.poly(
            dense
                .coeffs()
                .iter()
                .map(|c| Rational::from(c.clone()))
                .collect(),
        )
    };
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let integral = qx
        .integrate_rational_with(&to_q(num), &to_q(den), cancel)
        .map_err(|Cancelled| EvalError::Interrupted)?;
    Ok(rational_expr(&integral, x))
}

/// The expression of cas-poly's `∫ A / D`: the polynomial and rational
/// parts, the logarithms of the primitive integer multiples of their
/// arguments, of their absolute values where they have real roots, as
/// giac's `ln(abs(x-1))` and `ln(x^2+1)`, and what is left under
/// `Integrate`. `None` when nothing integrated, so that what is left does
/// not come back as itself.
fn rational_expr(r: &RationalIntegral, x: &Expr) -> Option<Expr> {
    let (gn, gd) = &r.rational;
    if r.polynomial.is_zero() && gn.is_zero() && r.logarithms.is_empty() {
        return None;
    }
    let mut terms = vec![dense_expr(&r.polynomial, x)];
    if !gn.is_zero() {
        terms.push(dense_expr(gn, x) / dense_expr(gd, x));
    }
    for log in &r.logarithms {
        let argument = dense_expr(&primitive(&log.argument), x);
        let argument = if log.real_root {
            call(Symbol::ABS, argument)
        } else {
            argument
        };
        terms.push(Expr::from(log.coefficient.clone()) * call(Symbol::LOG, argument));
    }
    if let Some((p, q)) = &r.rest {
        let rest = dense_expr(p, x) / dense_expr(q, x);
        terms.push(Expr::call(Symbol::INTEGRATE, [rest, x.clone()]));
    }
    Some(terms.into_iter().sum())
}

/// `p` over ℚ as a polynomial in `x`, built.
fn dense_expr(p: &Dense<Rational>, x: &Expr) -> Expr {
    p.coeffs()
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| {
            let k = i64::try_from(k).expect("a degree in memory");
            Expr::from(c.clone()) * x.pow(k)
        })
        .sum()
}

/// The primitive integer multiple of `p` with a positive leading
/// coefficient: `2x − 1` for `x − 1/2`.
fn primitive(p: &Dense<Rational>) -> Dense<Rational> {
    let l = p
        .coeffs()
        .iter()
        .fold(Integer::ONE, |l, c| l.lcm(c.denom()));
    let scaled: Vec<Integer> = p
        .coeffs()
        .iter()
        .map(|c| (c * &Rational::from(l.clone())).numer().clone())
        .collect();
    let mut g = scaled.iter().fold(Integer::ZERO, |g, c| g.gcd(c));
    if p.lead().is_some_and(Rational::is_negative) {
        g = -g;
    }
    DenseRing::new(Rationals).poly(
        scaled
            .into_iter()
            .map(|c| Rational::new(c, g.clone()))
            .collect(),
    )
}

/// `num / den` integrated term by term in the variable `k`, `den` free of
/// it: with `l` the lcm of the `e + 1`, `Σ (l / (e + 1)) c x^(e+1) / (l den)`,
/// in integers.
fn terms_integrated<R: Numbers>(
    p: &Polynomials,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
    k: usize,
) -> Option<Expr> {
    if ring.degree(den, k).unwrap_or(0) > 0 {
        return None;
    }
    let l = num
        .terms()
        .map(|(e, _)| Integer::from(e[k] + 1))
        .fold(Integer::ONE, |l, m| l.lcm(&m));
    let base = ring.base();
    let shifted = ring.from_terms(num.terms().map(|(mut e, c)| {
        let m = Integer::from(e[k] + 1);
        e[k] += 1;
        (e, base.mul(c, &base.elem_integer(&l.div_exact(&m))))
    }));
    let den = ring.scale(den, &base.elem_integer(&l));
    Some(p.ratio_expr(&shifted, &den))
}

/// Whether `f` is a polynomial in `x`: no other generator, no denominator,
/// contains it, nor does an unevaluated `Integrate`.
fn polynomial_in(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<bool, EvalError> {
    if f.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(false);
    }
    let Some((p, k)) = in_x(ev, f, f, x)? else {
        return Ok(false);
    };
    let n = p.generators().len();
    Ok(match &p.fractions()[0] {
        Fraction::Integers { den, .. } => {
            SparseRing::new(Integers, n).degree(den, k).unwrap_or(0) == 0
        }
        Fraction::GaussianIntegers { den, .. } => {
            SparseRing::new(GaussianIntegers, n)
                .degree(den, k)
                .unwrap_or(0)
                == 0
        }
    })
}

#[cfg(test)]
mod tests;
