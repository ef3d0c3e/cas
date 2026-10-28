//! Exact values of the elementary functions (CAS-200), as giac gives them:
//!
//! - `Sin`, `Cos` and `Tan` at rational multiples of `Pi` with a
//!   denominator 1, 2, 3, 4 or 6, in radicals: `Sin[Pi/3]` is `√3/2`,
//!   `Tan[Pi/2]` is `ComplexInfinity`; at an argument plus a multiple of
//!   `Pi/2` (of `Pi` for `Tan`), the function of the argument:
//!   `Sin[x + Pi/2]` is `Cos[x]`, `Sin[x + Pi]` is `-Sin[x]`;
//! - the sign of a negative number or of a product with a negative
//!   coefficient out of the odd functions `Sin`, `Tan`, `ArcSin` and
//!   `ArcTan`, and dropped by `Cos`: `Sin[-x]` is `-Sin[x]`;
//! - `Exp[0]` is 1 and `Exp[Log[u]]` is `u`; `Log[1]` is 0, `Log[E]` is 1
//!   and `Log[E^q]` is `q` for a rational `q`;
//! - `ArcSin`, `ArcCos` and `ArcTan` at the values the others take:
//!   `ArcSin[1/2]` is `Pi/6`, `ArcCos[-1]` is `Pi`, `ArcTan[√3]` is `Pi/3`.
//!
//! Other exact arguments stay as they are, as giac's `sin(1)` and
//! `sin(pi/5)` would not in radicals of this table.
//!
//! The hyperbolic functions (CAS-DOC-15 §2, CAS-238): 0 at 0 (`Cosh` 1,
//! `ArcCosh[1]` 0 and `ArcCosh[0]` `i π/2`), `ArcTanh[±1]` `±Infinity`,
//! `Sinh`, `Cosh` and `Tanh` of `Log[u]` in `u`, of `i y` as `i Sin[y]`,
//! `Cos[y]` and `i Tan[y]`, `ArcSinh` and `ArcTanh` of `i y` as `i ArcSin[y]`
//! and `i ArcTan[y]`, and the sign of a negative argument out (`Cosh` even).
//!
//! Compositions with the inverse functions ([`composition`], CAS-DOC-14 §5,
//! CAS-231), as giac evaluates them: `Sin`, `Cos` and `Tan` of `ArcSin`,
//! `ArcCos` and `ArcTan` for every argument (`Cos[ArcSin[x]]` is
//! `Sqrt[1 − x^2]`), `Log` of `Exp[u]` or `E^u` for a real `u`, and
//! `ArcSin`, `ArcCos` and `ArcTan` of `Sin`, `Cos` and `Tan` for an
//! argument the assumptions put in the principal range: within
//! `[−157/100, 157/100]`, inside `(−π/2, π/2)`, for `ArcSin` and `ArcTan`,
//! and `[0, 314/100]` for `ArcCos`.
//!
//! Integer multiples of `π` ([`periodic`], CAS-DOC-25 §2, CAS-313), where
//! the assumptions tell that `k` is an integer: `sin(kπ + y)` is
//! `(−1)^k sin y`, `cos(kπ + y)` is `(−1)^k cos y` and `tan(kπ + y)` is
//! `tan y`, even coefficients dropped and odd ones reduced (`(−1)^(3n)` is
//! `(−1)^n`); in the exponent of `Exp` and `E^u` ([`exponential`]), `iπk`
//! for a half-integer or such a `k`: `e^(iπ)` is `−1`, `e^(iπ/2)` is `i`,
//! `e^(x + 2iπn)` is `e^x`. Other multiples stay, as `e^(iπ/3)`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of trigonometry: f(x) = sin(r π), s and c"
)]

use cas_num::{Integer, Rational};

use crate::assume::{Set, Truth, element, range};
use crate::context::Context;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn call(head: Symbol, arg: Expr) -> Expr {
    Expr::call(head, [arg])
}

/// `c √k`, built.
fn surd(c: Rational, k: i64) -> Expr {
    if k == 1 {
        return Expr::from(c);
    }
    Expr::from(c) * Expr::from(k).pow(Expr::from(q(1, 2)))
}

/// `x` as `c √k` for a rational `c` and `k` 1, 2 or 3: the values the
/// inverse functions know.
fn as_surd(x: &Expr) -> Option<(Rational, i64)> {
    let root = |p: &Expr| -> Option<i64> {
        match p.args() {
            [k, h]
                if p.head() == Some(&Expr::from(Symbol::POWER)) && h.to_rational()? == q(1, 2) =>
            {
                k.as_i64().filter(|k| matches!(k, 2 | 3))
            }
            _ => None,
        }
    };
    if let Some(c) = x.to_rational() {
        return Some((c, 1));
    }
    if let Some(k) = root(x) {
        return Some((Rational::ONE, k));
    }
    match x.args() {
        [c, p] if x.head() == Some(&Expr::from(Symbol::TIMES)) => {
            Some((c.to_rational()?, root(p)?))
        }
        _ => None,
    }
}

/// `x` as `q Pi + rest` for a rational `q`: `rest` is `None` for `q Pi`
/// alone, 0 included; `None` without such a term.
fn pi_multiple(x: &Expr) -> Option<(Rational, Option<Expr>)> {
    let term = |t: &Expr| -> Option<Rational> {
        if t.as_symbol() == Some(Symbol::PI) {
            return Some(Rational::ONE);
        }
        match t.args() {
            [c, p]
                if t.head() == Some(&Expr::from(Symbol::TIMES))
                    && p.as_symbol() == Some(Symbol::PI) =>
            {
                c.to_rational()
            }
            _ => None,
        }
    };
    if *x == Expr::from(0) {
        return Some((Rational::ZERO, None));
    }
    if let Some(c) = term(x) {
        return Some((c, None));
    }
    if x.head() != Some(&Expr::from(Symbol::PLUS)) {
        return None;
    }
    let args = x.args();
    let i = args.iter().position(|t| term(t).is_some())?;
    let mut rest: Vec<Expr> = args.to_vec();
    let c = term(&rest.remove(i)).expect("a multiple of Pi");
    let rest = match rest.len() {
        1 => rest.pop().expect("one"),
        _ => Expr::call(Symbol::PLUS, rest),
    };
    Some((c, Some(rest)))
}

/// `r` reduced into `[0, m)`.
fn modulo(r: &Rational, m: &Rational) -> Rational {
    let k = Rational::from((r / m).floor());
    r - &(&k * m)
}

/// `sin(r π)` for a denominator of `r` 1, 2, 3, 4 or 6.
fn sin_pi(r: &Rational) -> Option<Expr> {
    let r = modulo(r, &Rational::from(Integer::from(2)));
    let sixths = &r * &Rational::from(Integer::from(6));
    if sixths.is_integer() {
        // sin(m π/6), m from 0 to 11.
        let (c, k) = [
            (q(0, 1), 1),
            (q(1, 2), 1),
            (q(1, 2), 3),
            (q(1, 1), 1),
            (q(1, 2), 3),
            (q(1, 2), 1),
            (q(0, 1), 1),
            (q(-1, 2), 1),
            (q(-1, 2), 3),
            (q(-1, 1), 1),
            (q(-1, 2), 3),
            (q(-1, 2), 1),
        ][usize::try_from(sixths.numer().to_i64()?).ok()?]
        .clone();
        return Some(surd(c, k));
    }
    let quarters = &r * &Rational::from(Integer::from(4));
    if quarters.is_integer() {
        // sin(m π/4) for an odd m from 1 to 7.
        let m = quarters.numer().to_i64()?;
        let sign = if m < 4 { 1 } else { -1 };
        return Some(surd(q(sign, 2), 2));
    }
    None
}

/// `Sin`, `Cos` or `Tan` at `x` (the module's first rule).
fn trigonometric(f: Symbol, x: &Expr) -> Option<Expr> {
    let half = q(1, 2);
    let (r, rest) = pi_multiple(x)?;
    let Some(rest) = rest else {
        return match f {
            Symbol::SIN => sin_pi(&r),
            Symbol::COS => sin_pi(&(&r + &half)),
            _ => {
                let (s, c) = (sin_pi(&r)?, sin_pi(&(&r + &half))?);
                Some(if c == Expr::from(0) {
                    Expr::from(Symbol::COMPLEX_INFINITY)
                } else {
                    s / c
                })
            }
        };
    };
    // f(rest + r π) for a multiple r of 1/2: quarter turns of sin.
    let turns = |r: &Rational| -> Option<i64> {
        let halves = r * &Rational::from(Integer::from(2));
        halves
            .is_integer()
            .then(|| halves.numer().to_i64())
            .flatten()
            .map(|k| k.rem_euclid(4))
    };
    let sine = |k: i64| match k {
        0 => call(Symbol::SIN, rest.clone()),
        1 => call(Symbol::COS, rest.clone()),
        2 => -call(Symbol::SIN, rest.clone()),
        _ => -call(Symbol::COS, rest.clone()),
    };
    match f {
        Symbol::SIN => Some(sine(turns(&r)?)),
        Symbol::COS => Some(sine(turns(&(&r + &half))?)),
        _ => r.is_integer().then(|| call(Symbol::TAN, rest.clone())),
    }
}

/// `-x` when `x` is a negative number or a product with a negative
/// rational coefficient.
pub(super) fn negated(x: &Expr) -> Option<Expr> {
    if let Some(r) = x.to_rational() {
        return r.is_negative().then(|| Expr::from(-r));
    }
    let [c, rest @ ..] = x.args() else {
        return None;
    };
    if x.head() != Some(&Expr::from(Symbol::TIMES)) {
        return None;
    }
    let c = c.to_rational().filter(Rational::is_negative)?;
    let c = -c;
    let mut factors: Vec<Expr> = Vec::with_capacity(rest.len() + 1);
    if c != Rational::ONE {
        factors.push(Expr::from(c));
    }
    factors.extend(rest.iter().cloned());
    Some(match factors.len() {
        1 => factors.pop().expect("one"),
        _ => Expr::call(Symbol::TIMES, factors),
    })
}

/// `asin x / π` at the values `sin` takes at multiples of `π/6` and `π/4`.
fn asin_pi(x: &Expr) -> Option<Rational> {
    let (c, k) = as_surd(x)?;
    let (sign, c) = if c.is_negative() { (-1, -c) } else { (1, c) };
    let r = match (k, c) {
        (1, c) if c.is_zero() => q(0, 1),
        (1, c) if c == q(1, 2) => q(1, 6),
        (2, c) if c == q(1, 2) => q(1, 4),
        (3, c) if c == q(1, 2) => q(1, 3),
        (1, c) if c == q(1, 1) => q(1, 2),
        _ => return None,
    };
    Some(&r * &Rational::from(Integer::from(sign)))
}

/// `atan x / π` at the values `tan` takes at multiples of `π/6` and `π/4`.
fn atan_pi(x: &Expr) -> Option<Rational> {
    let (c, k) = as_surd(x)?;
    let (sign, c) = if c.is_negative() { (-1, -c) } else { (1, c) };
    let r = match (k, c) {
        (1, c) if c.is_zero() => q(0, 1),
        (3, c) if c == q(1, 3) => q(1, 6),
        (1, c) if c == q(1, 1) => q(1, 4),
        (3, c) if c == q(1, 1) => q(1, 3),
        _ => return None,
    };
    Some(&r * &Rational::from(Integer::from(sign)))
}

/// `r Pi`, built.
fn pi_times(r: Rational) -> Expr {
    Expr::from(r) * Expr::from(Symbol::PI)
}

/// The exact value of `f[x]` (the module's rules), `None` to leave it to
/// the rest of the native.
pub(super) fn value(f: Symbol, x: &Expr) -> Option<Expr> {
    match f {
        Symbol::SIN | Symbol::COS | Symbol::TAN => {
            if let Some(v) = trigonometric(f, x) {
                return Some(v);
            }
            let y = negated(x)?;
            Some(match f {
                Symbol::COS => call(f, y),
                _ => -call(f, y),
            })
        }
        Symbol::EXP => {
            if *x == Expr::from(0) {
                return Some(Expr::from(1));
            }
            match x.args() {
                [u] if x.head() == Some(&Expr::from(Symbol::LOG)) => Some(u.clone()),
                _ => None,
            }
        }
        Symbol::LOG => {
            if *x == Expr::from(1) {
                return Some(Expr::from(0));
            }
            if x.as_symbol() == Some(Symbol::E) {
                return Some(Expr::from(1));
            }
            match (x.head().and_then(Expr::as_symbol), x.args()) {
                (Some(Symbol::POWER), [b, r]) if b.as_symbol() == Some(Symbol::E) => {
                    r.to_rational().map(Expr::from)
                }
                (Some(Symbol::EXP), [r]) => r.to_rational().map(Expr::from),
                _ => None,
            }
        }
        Symbol::ARC_SIN => match asin_pi(x) {
            Some(r) => Some(pi_times(r)),
            None => negated(x).map(|y| -call(f, y)),
        },
        Symbol::ARC_COS => asin_pi(x).map(|r| pi_times(&q(1, 2) - &r)),
        Symbol::ARC_TAN => match atan_pi(x) {
            Some(r) => Some(pi_times(r)),
            None => negated(x).map(|y| -call(f, y)),
        },
        Symbol::SINH | Symbol::COSH | Symbol::TANH => hyperbolic(f, x),
        Symbol::ARC_SINH | Symbol::ARC_TANH => {
            if *x == Expr::from(0) {
                return Some(Expr::from(0));
            }
            if f == Symbol::ARC_TANH && x.to_rational().is_some_and(|q| q.abs() == Rational::ONE) {
                let infinity = Expr::from(Symbol::INFINITY);
                return Some(if x.to_rational()?.is_negative() {
                    -infinity
                } else {
                    infinity
                });
            }
            // asinh(i y) = i asin(y), atanh(i y) = i atan(y).
            if let Some(y) = imaginary(x) {
                let g = if f == Symbol::ARC_SINH {
                    Symbol::ARC_SIN
                } else {
                    Symbol::ARC_TAN
                };
                return Some(Expr::from(Symbol::I) * call(g, y));
            }
            negated(x).map(|y| -call(f, y))
        }
        Symbol::ARC_COSH => {
            if *x == Expr::from(1) {
                return Some(Expr::from(0));
            }
            // The principal branch, log(z + √(z+1) √(z−1)): i π/2 at 0.
            (*x == Expr::from(0)).then(|| Expr::from(Symbol::I) * pi_times(q(1, 2)))
        }
        _ => None,
    }
}

/// `y` for `x = i y`, `y` real in form: a pure imaginary number, or a
/// product with one as its coefficient.
fn imaginary(x: &Expr) -> Option<Expr> {
    let pure = |c: &Expr| match c.view() {
        View::Complex(z) if z.re().is_zero() => Some(z.im().clone()),
        _ => None,
    };
    if let Some(b) = pure(x) {
        return Some(Expr::from(b));
    }
    let (Some(Symbol::TIMES), [c, rest @ ..]) = (x.head().and_then(Expr::as_symbol), x.args())
    else {
        return None;
    };
    let b = pure(c)?;
    let mut factors: Vec<Expr> = Vec::with_capacity(rest.len() + 1);
    if b != Rational::ONE {
        factors.push(Expr::from(b));
    }
    factors.extend(rest.iter().cloned());
    Some(match factors.len() {
        1 => factors.pop().expect("one"),
        _ => Expr::call(Symbol::TIMES, factors),
    })
}

/// `Sinh`, `Cosh` or `Tanh` of `x`: at 0, of `Log[u]` (`sinh(ln u) = (u −
/// 1/u)/2`), of `i y` (`sinh(i y) = i sin y`, `cosh(i y) = cos y`), and the
/// sign of a negative `x` out (`Cosh` even, the others odd).
fn hyperbolic(f: Symbol, x: &Expr) -> Option<Expr> {
    if *x == Expr::from(0) {
        return Some(Expr::from(i64::from(f == Symbol::COSH)));
    }
    if let (Some(Symbol::LOG), [u]) = (x.head().and_then(Expr::as_symbol), x.args()) {
        let (u, inverse) = (u.clone(), u.pow(-1));
        let half = Expr::from(q(1, 2));
        return Some(match f {
            Symbol::SINH => half * (u - inverse),
            Symbol::COSH => half * (u + inverse),
            _ => (u.pow(2) - Expr::from(1)) / (u.pow(2) + Expr::from(1)),
        });
    }
    if let Some(y) = imaginary(x) {
        return Some(match f {
            Symbol::SINH => Expr::from(Symbol::I) * call(Symbol::SIN, y),
            Symbol::COSH => call(Symbol::COS, y),
            _ => Expr::from(Symbol::I) * call(Symbol::TAN, y),
        });
    }
    let y = negated(x)?;
    Some(match f {
        Symbol::COSH => call(f, y),
        _ => -call(f, y),
    })
}

/// `f[x]` for `x` an inverse function's value, or the other way (the
/// module's compositions); `None` to leave it to the rest of the native.
pub(super) fn composition(ctx: &Context, f: Symbol, x: &Expr) -> Option<Expr> {
    let real = |u: &Expr| element(ctx, u, Set::Reals) == Truth::True;
    if f == Symbol::LOG
        && let (Some(Symbol::POWER), [b, u]) = (x.head().and_then(Expr::as_symbol), x.args())
    {
        return (b.as_symbol() == Some(Symbol::E) && real(u)).then(|| u.clone());
    }
    let (Some(g), [u]) = (x.head().and_then(Expr::as_symbol), x.args()) else {
        return None;
    };
    // √(1 ± u²).
    let root = |sign: i64| (Expr::from(1) + Expr::from(sign) * u.pow(2)).pow(Expr::from(q(1, 2)));
    let within = |low: Rational, high: Rational| match range(ctx, u).interval() {
        Some((Some(l), Some(h))) => *l >= low && *h <= high,
        _ => false,
    };
    let half_pi = || within(q(-157, 100), q(157, 100));
    match (f, g) {
        (Symbol::SIN, Symbol::ARC_SIN)
        | (Symbol::COS, Symbol::ARC_COS)
        | (Symbol::TAN, Symbol::ARC_TAN)
        | (Symbol::SINH, Symbol::ARC_SINH)
        | (Symbol::COSH, Symbol::ARC_COSH)
        | (Symbol::TANH, Symbol::ARC_TANH) => Some(u.clone()),
        (Symbol::COS, Symbol::ARC_SIN) | (Symbol::SIN, Symbol::ARC_COS) => Some(root(-1)),
        (Symbol::TAN, Symbol::ARC_SIN) => Some(u.clone() / root(-1)),
        (Symbol::TAN, Symbol::ARC_COS) => Some(root(-1) / u.clone()),
        (Symbol::SIN, Symbol::ARC_TAN) => Some(u.clone() / root(1)),
        (Symbol::COS, Symbol::ARC_TAN) => Some(root(1).pow(-1)),
        (Symbol::LOG, Symbol::EXP) => real(u).then(|| u.clone()),
        (Symbol::ARC_SIN, Symbol::SIN) | (Symbol::ARC_TAN, Symbol::TAN) => {
            half_pi().then(|| u.clone())
        }
        (Symbol::ARC_COS, Symbol::COS) => within(q(0, 1), q(314, 100)).then(|| u.clone()),
        _ => None,
    }
}

/// `t/π`, or `t/(iπ)` when `imaginary`, for a term `t` with the factor
/// `π` (and a pure imaginary coefficient when `imaginary`); `None` for
/// another term.
fn over_pi(t: &Expr, imaginary: bool) -> Option<Expr> {
    let factors: &[Expr] = match (t.head().and_then(Expr::as_symbol), t.args()) {
        (Some(Symbol::TIMES), args) => args,
        _ => std::slice::from_ref(t),
    };
    let i = factors
        .iter()
        .position(|f| f.as_symbol() == Some(Symbol::PI))?;
    let mut unit = !imaginary;
    let mut rest: Vec<Expr> = Vec::with_capacity(factors.len());
    for f in factors[..i].iter().chain(&factors[i + 1..]) {
        match f.view() {
            View::Complex(z) if !unit && z.re().is_zero() => {
                unit = true;
                if *z.im() != Rational::ONE {
                    rest.push(Expr::from(z.im().clone()));
                }
            }
            View::Symbol(Symbol::I) if !unit => unit = true,
            View::Complex(_) | View::Symbol(Symbol::I) => return None,
            _ => rest.push(f.clone()),
        }
    }
    unit.then(|| match rest.len() {
        0 => Expr::from(1),
        1 => rest.pop().expect("one"),
        _ => Expr::call(Symbol::TIMES, rest),
    })
}

/// An argument as `π (q + Σ k) + Σ rest`, with `iπ` for `π` when
/// imaginary: a rational `q`, terms `k` known integers, and the other terms.
struct Turns {
    q: Rational,
    integers: Vec<Expr>,
    rest: Vec<Expr>,
}

impl Turns {
    fn of(ctx: &Context, x: &Expr, imaginary: bool) -> Turns {
        let terms: &[Expr] = match x.head().and_then(Expr::as_symbol) {
            Some(Symbol::PLUS) => x.args(),
            _ => std::slice::from_ref(x),
        };
        let mut turns = Turns {
            q: Rational::ZERO,
            integers: Vec::new(),
            rest: Vec::new(),
        };
        for t in terms {
            match over_pi(t, imaginary) {
                Some(m) => match m.to_rational() {
                    Some(r) => turns.q = &turns.q + &r,
                    None if element(ctx, &m, Set::Integers) == Truth::True => {
                        turns.integers.push(m);
                    }
                    None => turns.rest.push(t.clone()),
                },
                None => turns.rest.push(t.clone()),
            }
        }
        turns
    }

    /// `(−1)^(Σ k)`, `None` for 1: the terms with an even coefficient
    /// dropped, an odd one reduced to 1 where the rest is an integer.
    fn sign(&self, ctx: &Context) -> Option<Expr> {
        let mut odd: Vec<Expr> = Vec::new();
        for k in &self.integers {
            let (c, rest) = match k.args() {
                [c, rest @ ..] if k.head() == Some(&Expr::from(Symbol::TIMES)) => {
                    match c.to_rational().filter(Rational::is_integer) {
                        Some(c) => (c, rest),
                        None => (Rational::ONE, k.args()),
                    }
                }
                _ => (Rational::ONE, std::slice::from_ref(k)),
            };
            if c.numer()
                .unsigned_abs()
                .trailing_zeros()
                .is_none_or(|t| t > 0)
            {
                continue;
            }
            let m = match rest {
                [m] => m.clone(),
                _ => Expr::call(Symbol::TIMES, rest.to_vec()),
            };
            odd.push(if element(ctx, &m, Set::Integers) == Truth::True {
                m
            } else {
                k.clone()
            });
        }
        match odd.len() {
            0 => None,
            1 => Some(Expr::from(-1).pow(odd.pop().expect("one"))),
            _ => Some(Expr::from(-1).pow(Expr::call(Symbol::PLUS, odd))),
        }
    }
}

fn sum(mut terms: Vec<Expr>) -> Expr {
    match terms.len() {
        0 => Expr::from(0),
        1 => terms.pop().expect("one"),
        _ => Expr::call(Symbol::PLUS, terms),
    }
}

/// `Sin`, `Cos` or `Tan` of `x` with terms `kπ`, `k` known an integer
/// (CAS-DOC-25 §2): `(−1)^k` out of `Sin` and `Cos`, nothing out of `Tan`;
/// `Exp` of `x` by [`exponential`]; `None` without such a term.
pub(super) fn periodic(ctx: &Context, f: Symbol, x: &Expr) -> Option<Expr> {
    if f == Symbol::EXP {
        let (factor, rest) = exponential(ctx, x)?;
        return Some(factor * call(Symbol::EXP, rest));
    }
    if !matches!(f, Symbol::SIN | Symbol::COS | Symbol::TAN) {
        return None;
    }
    let turns = Turns::of(ctx, x, false);
    if turns.integers.is_empty() {
        return None;
    }
    let sign = turns.sign(ctx);
    let mut terms = turns.rest;
    if !turns.q.is_zero() {
        terms.push(pi_times(turns.q));
    }
    let g = call(f, sum(terms));
    Some(match sign {
        Some(s) if f != Symbol::TAN => s * g,
        _ => g,
    })
}

/// `e^u` as `c e^v` (CAS-DOC-25 §2), for the terms `iπk` of `u` with `k` a
/// half-integer or known an integer: `c = i^(2q) (−1)^(Σ k)` and `v` the
/// rest of `u`; `None` without such a term.
pub(super) fn exponential(ctx: &Context, u: &Expr) -> Option<(Expr, Expr)> {
    let turns = Turns::of(ctx, u, true);
    let twice = &turns.q * &Rational::from(Integer::from(2));
    let quarter_turns = twice.is_integer() && !turns.q.is_zero();
    if turns.integers.is_empty() && !quarter_turns {
        return None;
    }
    let sign = turns.sign(ctx);
    let mut rest = turns.rest;
    let mut factor = sign.unwrap_or_else(|| Expr::from(1));
    if quarter_turns {
        let (_, k) = twice.numer().div_rem_euclid(&Integer::from(4));
        let unit = match k.to_i64() {
            Some(0) => Expr::from(1),
            Some(1) => Expr::from(Symbol::I),
            Some(2) => Expr::from(-1),
            _ => -Expr::from(Symbol::I),
        };
        factor = unit * factor;
    } else if !turns.q.is_zero() {
        rest.push(Expr::from(Symbol::I) * pi_times(turns.q));
    }
    Some((factor, sum(rest)))
}

#[cfg(test)]
mod tests;
