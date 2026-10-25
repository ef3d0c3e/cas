//! `Integrate` (CAS-DOC-9, CAS-193): antiderivatives and definite integrals
//! in one symbol, Xcas's `integrate` and `int`, without a constant of
//! integration. Here, CAS-194:
//!
//! - lists thread over the integrand;
//! - what is free of the variable `x` is `f x`;
//! - sums go term by term, and the terms that do not integrate stay under
//!   one `Integrate`, as giac's `x + f(x)` gives `x^2/2 + integrate(f(x), x)`;
//! - the factors of a product free of `x` come out of its antiderivative;
//!   when the rest does not integrate, the product stays whole, as giac's
//!   `integrate(2*f(x),x)`;
//! - `x^n` with `n` free of `x` is `x^(n+1) / (n+1)`, giac's answer even for
//!   a symbolic `n`, and `x^-1` is `Log[Abs[x]]`;
//! - polynomials in `x`, through the conversion of CAS-DOC-8 §2: each term
//!   `c x^k` is `c x^(k+1) / (k+1)`, over the common denominator.
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

use cas_num::Integer;
use cas_poly::{Sparse, SparseRing};

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
                return polynomial(ev, e, f, x);
            }
            // When the rest does not integrate, the product stays whole, as
            // giac leaves integrate(2*f(x),x).
            let varying: Expr = varying.into_iter().product();
            Ok(integral(ev, e, &varying, x)?.map(|g| constant.into_iter().product::<Expr>() * g))
        }
        Some(Symbol::POWER) if f.args()[0] == *x && !f.args()[1].contains(x) => {
            let n = &f.args()[1];
            if n.as_i64() == Some(-1) {
                let abs = Expr::call(Symbol::ABS, [x.clone()]);
                return Ok(Some(Expr::call(Symbol::LOG, [abs])));
            }
            let m = n.clone() + 1;
            Ok(Some(x.pow(m.clone()) / m))
        }
        _ => polynomial(ev, e, f, x),
    }
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

/// The antiderivative of `f` as a polynomial in `x`, when it is one:
/// `c x^k` to `c x^(k+1) / (k+1)`, over the denominator, which is free of
/// `x`.
fn polynomial(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some((p, k)) = in_x(ev, e, f, x)? else {
        return Ok(None);
    };
    let n = p.generators().len();
    Ok(match &p.fractions()[0] {
        Fraction::Integers { num, den } => terms_integrated(
            &p,
            &SparseRing::new(cas_num::algebra::Integers, n),
            num,
            den,
            k,
        ),
        Fraction::GaussianIntegers { num, den } => terms_integrated(
            &p,
            &SparseRing::new(cas_num::algebra::GaussianIntegers, n),
            num,
            den,
            k,
        ),
    })
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
            SparseRing::new(cas_num::algebra::Integers, n)
                .degree(den, k)
                .unwrap_or(0)
                == 0
        }
        Fraction::GaussianIntegers { den, .. } => {
            SparseRing::new(cas_num::algebra::GaussianIntegers, n)
                .degree(den, k)
                .unwrap_or(0)
                == 0
        }
    })
}

#[cfg(test)]
mod tests;
