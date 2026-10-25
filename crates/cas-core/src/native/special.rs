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

#![allow(
    clippy::many_single_char_names,
    reason = "the names of trigonometry: f(x) = sin(r π), s and c"
)]

use cas_num::{Integer, Rational};

use crate::expr::Expr;
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
fn negated(x: &Expr) -> Option<Expr> {
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
        _ => None,
    }
}

#[cfg(test)]
mod tests;
