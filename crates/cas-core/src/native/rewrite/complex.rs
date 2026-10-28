//! Complex parts (CAS-DOC-36 §7, CAS-401): `ComplexExpand`, giac's
//! `evalc`, and `MultiplyConjugate` and `MultiplyComplexConjugate`, giac's
//! `mult_conjugate` and `mult_c_conjugate`.
//!
//! `ComplexExpand` takes every symbol as real: it works in a scope that
//! assumes `Element[s, Reals]` for each, and writes `z = a + i b` by parts:
//!
//! - numbers by theirs, and what the scope makes real as itself;
//! - sums by parts, products by complex arithmetic;
//! - `Exp[a + i b] = e^a cos b + i e^a sin b`, and `E^` the same;
//! - `Log[z] = ln|z| + i Arg z`;
//! - `z^n` for an integer `|n| ≤` [`MAX_POWER`] by products, through the
//!   conjugate for a negative one; an odd root of a real, real in the Real
//!   domain (O7); other powers `z^(c + i d)` by their modulus `|z|^c e^(−d
//!   Arg z)` and angle `c Arg z + d ln|z|`;
//! - `Sin`, `Cos`, `Sinh` and `Cosh` by the addition formulas, `Tan[a + i
//!   b] = (sin 2a + i sinh 2b)/(cos 2a + cosh 2b)` and `Tanh[a + i b] =
//!   (sinh 2a + i sin 2b)/(cosh 2a + cos 2b)`;
//! - `Abs` as `√(a² + b²)`, and `Re`, `Im`, `Conjugate` and `Arg` by the
//!   parts;
//! - any other part `w` as `Re[w] + i Im[w]`, or as itself when it is the
//!   whole expression.
//!
//! The answer `a + i b` is evaluated in the scope.
//!
//! `MultiplyConjugate`: the expression as `N/D`, `D` the first factor of
//! the denominator, to the power 1, that is a sum with a square root among
//! its terms, `D = r + s` with `r` the first such term: `N (s − r)/(s² −
//! r²)`, the new denominator expanded. Without one, a sum `r + s` with a
//! square root is `(s² − r²)/(s − r)`. `MultiplyComplexConjugate`: the
//! same with a factor `D = a + i b` of the denominator that is not real,
//! the symbols real, as `N (a − i b)/(a² + b²)`; without one, an expression
//! `a + i b` that is not real is `(a² + b²)/(a − i b)`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of complex numbers: z = a + i b, to the power w = c + i d"
)]

use cas_num::{Integer, Rational};

use super::{factors, head, terms};
use crate::assume::{Domain, Set, Truth, element};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The largest integer power expanded by products.
const MAX_POWER: i64 = 16;

/// `(a, b)` for `a + i b`, both real.
type Parts = (Expr, Expr);

fn zero() -> Expr {
    Expr::from(0)
}

fn i() -> Expr {
    Expr::from(Symbol::I)
}

fn half() -> Expr {
    Expr::from(Rational::new(Integer::ONE, Integer::from(2)))
}

fn call(f: Symbol, u: Expr) -> Expr {
    Expr::call(f, [u])
}

fn is_zero(e: &Expr) -> bool {
    *e == zero()
}

/// The scope of `ctx` where every symbol of `x` is real.
fn real_scope(ctx: &Context, x: &Expr) -> Context {
    let mut scope = ctx.clone();
    for s in x.symbols() {
        if !s.is_builtin() {
            scope = scope.assuming(Expr::call(
                Symbol::ELEMENT,
                [Expr::from(s), Expr::from(Symbol::REALS)],
            ));
        }
    }
    scope
}

/// `(a + i b)(c + i d)`, built.
fn times((a, b): &Parts, (c, d): &Parts) -> Parts {
    if is_zero(b) && is_zero(d) {
        return (a.clone() * c.clone(), zero());
    }
    (
        a.clone() * c.clone() - b.clone() * d.clone(),
        a.clone() * d.clone() + b.clone() * c.clone(),
    )
}

/// `e^(a + i b)`, built.
fn exp((a, b): &Parts) -> Parts {
    let scale = call(Symbol::EXP, a.clone());
    if is_zero(b) {
        return (scale, zero());
    }
    (
        scale.clone() * call(Symbol::COS, b.clone()),
        scale * call(Symbol::SIN, b.clone()),
    )
}

/// `(ln|z|, Arg z)` for `z = a + i b`, built.
fn log((a, b): &Parts) -> Parts {
    if is_zero(b) {
        return (
            call(Symbol::LOG, call(Symbol::ABS, a.clone())),
            call(Symbol::ARG, a.clone()),
        );
    }
    let norm = a.pow(2) + b.pow(2);
    (
        half() * call(Symbol::LOG, norm),
        call(Symbol::ARG, a.clone() + i() * b.clone()),
    )
}

/// `z^n` for an integer `n`, `|n| ≤` [`MAX_POWER`], `z = (a, b)`, built.
fn integer_power((a, b): &Parts, n: i64) -> Parts {
    let mut acc = (Expr::from(1), zero());
    for _ in 0..n.unsigned_abs() {
        acc = times(&acc, &(a.clone(), b.clone()));
    }
    if n >= 0 {
        return acc;
    }
    if is_zero(b) {
        return (acc.0.pow(-1), zero());
    }
    // 1/(a + i b)^m = (a − i b)^m/(a² + b²)^m: acc is (a + i b)^m.
    let norm = (a.pow(2) + b.pow(2)).pow(n);
    (acc.0 * norm.clone(), -acc.1 * norm)
}

/// `z^w`, `z = (a, b)` and `w = (c, d)`, by modulus and angle, built.
fn general_power((a, b): &Parts, (c, d): &Parts) -> Parts {
    let (size, angle, length) = if is_zero(b) {
        let size = call(Symbol::ABS, a.clone());
        (
            size.pow(c.clone()),
            call(Symbol::ARG, a.clone()),
            call(Symbol::LOG, size),
        )
    } else {
        let norm = a.pow(2) + b.pow(2);
        (
            norm.pow(c.clone() * half()),
            call(Symbol::ARG, a.clone() + i() * b.clone()),
            half() * call(Symbol::LOG, norm),
        )
    };
    let (modulus, phase) = if is_zero(d) {
        (size, c.clone() * angle)
    } else {
        (
            size * call(Symbol::EXP, -(d.clone() * angle.clone())),
            c.clone() * angle + d.clone() * length,
        )
    };
    (
        modulus.clone() * call(Symbol::COS, phase.clone()),
        modulus * call(Symbol::SIN, phase),
    )
}

/// The parts of `z = b^n`, built.
fn power(ctx: &Context, z: &Expr, b: &Expr, n: &Expr) -> Parts {
    if let Some(k) = n.as_i64().filter(|n| n.abs() <= MAX_POWER) {
        return integer_power(&parts(ctx, b), k);
    }
    let base = parts(ctx, b);
    // An odd root of a real is real in the Real domain (O7).
    let odd_root = ctx.domain() == Domain::Real
        && is_zero(&base.1)
        && n.to_rational()
            .is_some_and(|q| q.denom().unsigned_abs().bit(0));
    if odd_root {
        (z.clone(), zero())
    } else {
        general_power(&base, &parts(ctx, n))
    }
}

/// `f(a + i b)` for the circular and hyperbolic functions, built.
fn trigonometric(f: Symbol, (a, b): &Parts) -> Option<Parts> {
    let g = |h: Symbol, u: &Expr| call(h, u.clone());
    let two = |u: &Expr| Expr::from(2) * u.clone();
    Some(match f {
        Symbol::SIN => (
            g(Symbol::SIN, a) * g(Symbol::COSH, b),
            g(Symbol::COS, a) * g(Symbol::SINH, b),
        ),
        Symbol::COS => (
            g(Symbol::COS, a) * g(Symbol::COSH, b),
            -(g(Symbol::SIN, a) * g(Symbol::SINH, b)),
        ),
        Symbol::SINH => (
            g(Symbol::SINH, a) * g(Symbol::COS, b),
            g(Symbol::COSH, a) * g(Symbol::SIN, b),
        ),
        Symbol::COSH => (
            g(Symbol::COSH, a) * g(Symbol::COS, b),
            g(Symbol::SINH, a) * g(Symbol::SIN, b),
        ),
        Symbol::TAN => {
            let den = g(Symbol::COS, &two(a)) + g(Symbol::COSH, &two(b));
            (
                g(Symbol::SIN, &two(a)) / den.clone(),
                g(Symbol::SINH, &two(b)) / den,
            )
        }
        Symbol::TANH => {
            let den = g(Symbol::COSH, &two(a)) + g(Symbol::COS, &two(b));
            (
                g(Symbol::SINH, &two(a)) / den.clone(),
                g(Symbol::SIN, &two(b)) / den,
            )
        }
        _ => return None,
    })
}

/// `z`'s parts in `ctx`, where every symbol is real, built.
fn parts(ctx: &Context, z: &Expr) -> Parts {
    match z.view() {
        View::Complex(c) => return (Expr::from(c.re().clone()), Expr::from(c.im().clone())),
        View::ComplexBall(c) => return (Expr::from(c.re().clone()), Expr::from(c.im().clone())),
        View::Int(_) | View::Integer(_) | View::Rational(_) | View::Ball(_) => {
            return (z.clone(), zero());
        }
        _ => {}
    }
    if z.as_symbol() == Some(Symbol::I) {
        return (zero(), Expr::from(1));
    }
    // The parts themselves before the test of reality, as `Abs[x + i y]` is
    // real but expands.
    if let (Some(f @ (Symbol::ABS | Symbol::RE | Symbol::IM | Symbol::CONJUGATE)), [w]) =
        (head(z), z.args())
    {
        let (a, b) = parts(ctx, w);
        return match f {
            Symbol::ABS if is_zero(&b) => (call(Symbol::ABS, a), zero()),
            Symbol::ABS => ((a.pow(2) + b.pow(2)).pow(half()), zero()),
            Symbol::RE => (a, zero()),
            Symbol::IM => (b, zero()),
            _ => (a, -b),
        };
    }
    if element(ctx, z, Set::Reals) == Truth::True {
        return (z.clone(), zero());
    }
    let opaque = || (call(Symbol::RE, z.clone()), call(Symbol::IM, z.clone()));
    let Some(f) = head(z) else {
        return opaque();
    };
    match (f, z.args()) {
        (Symbol::PLUS, args) => {
            let (re, im): (Vec<Expr>, Vec<Expr>) = args.iter().map(|t| parts(ctx, t)).unzip();
            (Expr::call(Symbol::PLUS, re), Expr::call(Symbol::PLUS, im))
        }
        (Symbol::TIMES, args) => args.iter().fold((Expr::from(1), zero()), |acc, f| {
            times(&acc, &parts(ctx, f))
        }),
        (Symbol::POWER, [b, n]) if b.as_symbol() == Some(Symbol::E) => exp(&parts(ctx, n)),
        (Symbol::POWER, [b, n]) => power(ctx, z, b, n),
        (Symbol::EXP, [u]) => exp(&parts(ctx, u)),
        (Symbol::LOG, [w]) => log(&parts(ctx, w)),
        (Symbol::ARG, [_]) => (z.clone(), zero()),
        (g, [u]) => trigonometric(g, &parts(ctx, u)).unwrap_or_else(opaque),
        _ => opaque(),
    }
}

/// `ComplexExpand` of `x`, not a list or a relation, evaluated.
pub(super) fn complex_expand(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let scope = real_scope(ev.context(), x);
    let (a, b) = parts(&scope, x);
    if (a.clone(), b.clone()) == (call(Symbol::RE, x.clone()), call(Symbol::IM, x.clone())) {
        return Ok(x.clone());
    }
    ev.eval_in(&scope, &(a + i() * b)).0
}

/// `e`'s factors: those of the numerator, and the bases `D` of its factors
/// `D^(−1)`, the denominator's; other powers `D^(−k)` stay in the
/// numerator.
fn fraction(e: &Expr) -> (Vec<Expr>, Vec<Expr>) {
    let mut numerator = Vec::new();
    let mut denominator = Vec::new();
    for f in factors(e) {
        match (head(f), f.args()) {
            (Some(Symbol::POWER), [d, n]) if *n == Expr::from(-1) => denominator.push(d.clone()),
            _ => numerator.push(f.clone()),
        }
    }
    (numerator, denominator)
}

/// Whether `t` has a square root: a power whose exponent is a fraction of
/// denominator 2.
fn has_root(t: &Expr) -> bool {
    t.distinct().any(|part| {
        matches!(
            (head(part), part.args()),
            (Some(Symbol::POWER), [_, r])
                if r.to_rational().is_some_and(|q| *q.denom() == Integer::from(2))
        )
    })
}

/// `d = r + s` at its first term `r` with a square root, as `(r, s)`.
fn split_root(d: &Expr) -> Option<(Expr, Expr)> {
    let ts = terms(d);
    if ts.len() < 2 {
        return None;
    }
    let k = ts.iter().position(has_root)?;
    let rest: Vec<Expr> = ts[..k].iter().chain(&ts[k + 1..]).cloned().collect();
    Some((ts[k].clone(), Expr::call(Symbol::PLUS, rest)))
}

fn expand(e: Expr) -> Expr {
    Expr::call(Symbol::EXPAND, [e])
}

/// `N/D` times `conjugate/conjugate`, `D` the `k`-th factor of the
/// denominator, and `product` the expanded `D conjugate`, built.
fn multiplied(
    numerator: Vec<Expr>,
    denominator: &[Expr],
    k: usize,
    conjugate: Expr,
    product: &Expr,
) -> Expr {
    let mut fs = numerator;
    fs.push(conjugate);
    fs.extend(
        denominator
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != k)
            .map(|(_, d)| d.pow(-1)),
    );
    fs.push(product.pow(-1));
    Expr::call(Symbol::TIMES, fs)
}

/// `MultiplyConjugate` of `x`, not a list or a relation, built.
pub(super) fn multiply_conjugate(x: &Expr) -> Expr {
    let (numerator, denominator) = fraction(x);
    if let Some((k, (r, s))) = denominator
        .iter()
        .enumerate()
        .find_map(|(k, d)| split_root(d).map(|rs| (k, rs)))
    {
        let product = expand(s.pow(2) - r.pow(2));
        return multiplied(numerator, &denominator, k, s - r, &product);
    }
    match split_root(x) {
        Some((r, s)) if denominator.is_empty() => expand(s.pow(2) - r.pow(2)) / (s - r),
        _ => x.clone(),
    }
}

/// The parts of `d` in `scope`, where every symbol is real; `None` for a
/// real `d`.
fn complex(ev: &mut Eval<'_>, scope: &Context, d: &Expr) -> Result<Option<Parts>, EvalError> {
    let (a, b) = parts(scope, d);
    let b = ev.eval_in(scope, &b).0?;
    Ok((!is_zero(&b)).then_some((a, b)))
}

/// `MultiplyComplexConjugate` of `x`, not a list or a relation, built.
pub(super) fn multiply_complex_conjugate(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let scope = real_scope(ev.context(), x);
    let (numerator, denominator) = fraction(x);
    for (k, d) in denominator.iter().enumerate() {
        if let Some((a, b)) = complex(ev, &scope, d)? {
            let product = expand(a.pow(2) + b.pow(2));
            let conjugate = a - i() * b;
            return Ok(multiplied(numerator, &denominator, k, conjugate, &product));
        }
    }
    if denominator.is_empty()
        && let Some((a, b)) = complex(ev, &scope, x)?
    {
        return Ok(expand(a.pow(2) + b.pow(2)) / (a - i() * b));
    }
    Ok(x.clone())
}

#[cfg(test)]
mod tests;
