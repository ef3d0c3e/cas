//! Laplace transforms (CAS-DOC-18 §3, CAS-274): `LaplaceTransform[f, t,
//! s]`, Xcas's `laplace(f, t, s)`; `laplace(f, t)` transforms into `t`
//! again and `laplace(f)` in `x` into `x`, as giac's (`laplace(1)` is
//! `1/x`).
//!
//! By the first rule that applies, recursively:
//!
//! - **Linearity**: sums term by term, factors free of `t` out; `c` free of
//!   `t` is `c/s`.
//! - **Table**: `t^n` is `Γ(n + 1)/s^(n+1)` for a number `n > −1` or a
//!   symbol (giac's `laplace(t^(1/2),t,s)` is 0); `e^(a t + b)` is `e^b/(s −
//!   a)`; `sin(b t)`, `cos(b t)`, `sinh(b t)`, `cosh(b t)` the classical
//!   ones, a phase expanded first; `HeavisideTheta[t − c]` is `e^(−c s)/s`
//!   and `DiracDelta[t − c]` `e^(−c s)` for `c ≥ 0` (`1/s` and 0 for `c <
//!   0`); `J_n(b t)` and `I_n(b t)`, `(√(s² + b²) − s)^n/(b^n √(s² + b²))`
//!   and `(s − √(s² − b²))^n/(b^n √(s² − b²))` (CAS-DOC-34 §5).
//! - **Shifts**: `e^(a t) g(t)` is `G(s − a)`; `HeavisideTheta[t − c] g(t)`
//!   is `e^(−c s) L[g(t + c)]` for `c ≥ 0`.
//! - **Powers of `t`**: `tⁿ g(t)`, `n` a positive integer, is `(−1)ⁿ
//!   G⁽ⁿ⁾(s)`.
//! - Products and powers of sines and cosines are reduced (`TrigReduce`)
//!   first.
//! - **Unknown functions**: `f(t)` stays as `LaplaceTransform[f[t], t, s]`,
//!   `F`, and `f⁽ⁿ⁾(t)` is `sⁿ F − Σ_{k<n} s^(n−1−k) f⁽ᵏ⁾(0)`.
//!
//! The answer is `Together`'s, as giac's `1/(s^2-2*s+1)`; what no rule
//! transforms stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the transform's names: f(t), F(s), e^(a t), sin(b t), t − c"
)]

use std::cmp::Ordering;

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `LaplaceTransform[f, t, s]`, `[f, t]` and `[f]`.
pub(super) fn laplace_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let x = Expr::symbol("x");
    let (f, t, s) = match e.args() {
        [f] => (f, x.clone(), x),
        [f, t] => (f, t.clone(), t.clone()),
        [f, t, s] => (f, t.clone(), s.clone()),
        _ => return Ok(None),
    };
    if t.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    // Transformed into a fresh variable, put in for s after: s may be t.
    let fresh = Expr::from(Symbol::new("LaplaceTransform$s"));
    let Some(image) = transform(ev, f, &t, &fresh)? else {
        return Ok(None);
    };
    let image = image.substitute(&[(fresh, s)]);
    // An unknown function's transform is the call itself, which stays:
    // evaluated, it would come back here.
    if image == *e {
        return Ok(None);
    }
    let result = ev.eval(&Expr::call(Symbol::TOGETHER, [image]))?;
    ev.step("LaplaceTransform", e, &result);
    Ok(Some(result))
}

/// `(a, b)` with `u = a t + b`, both free of `t`.
fn linear(ev: &mut Eval<'_>, u: &Expr, t: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    let a = ev.eval(&Expr::call(Symbol::D, [u.clone(), t.clone()]))?;
    if a.contains(t) {
        return Ok(None);
    }
    let b = ev.eval(&u.substitute(&[(t.clone(), Expr::from(0))]))?;
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - a.clone() * t.clone() - b.clone()],
    ))?;
    Ok((gap == Expr::from(0)).then_some((a, b)))
}

/// The exponent `u` of `e^u`.
fn exponent(e: &Expr) -> Option<&Expr> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => Some(u),
        (Some(Symbol::EXP), [u]) => Some(u),
        _ => None,
    }
}

/// Whether `e` is a sine or a cosine, or a power of one.
fn trigonometric(e: &Expr) -> bool {
    let head = |e: &Expr| e.head().and_then(Expr::as_symbol);
    match (head(e), e.args()) {
        (Some(Symbol::SIN | Symbol::COS), _) => true,
        (Some(Symbol::POWER), [b, _]) => matches!(head(b), Some(Symbol::SIN | Symbol::COS)),
        _ => false,
    }
}

/// `c` with `HeavisideTheta[u]` the step at `t = c`, `u = α t − β`, `α > 0`.
fn step_at(ev: &mut Eval<'_>, u: &Expr, t: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some((a, b)) = linear(ev, u, t)? else {
        return Ok(None);
    };
    if real_sign(ev.context(), &a) != Some(Ordering::Greater) {
        return Ok(None);
    }
    ev.eval(&(-(b / a))).map(Some)
}

/// The transform of `f` in `t` into `s`; `None` when no rule applies.
fn transform(ev: &mut Eval<'_>, f: &Expr, t: &Expr, s: &Expr) -> Result<Option<Expr>, EvalError> {
    if !f.contains(t) {
        return Ok(Some(f.clone() / s.clone()));
    }
    if f == t {
        return Ok(Some(s.pow(-2)));
    }
    let (head, args) = (f.head().and_then(Expr::as_symbol), f.args());
    match (head, args) {
        (Some(Symbol::PLUS), terms) => {
            let mut images = Vec::with_capacity(terms.len());
            for term in terms {
                let Some(image) = transform(ev, term, t, s)? else {
                    return Ok(None);
                };
                images.push(image);
            }
            Ok(Some(Expr::call(Symbol::PLUS, images)))
        }
        (Some(Symbol::TIMES), factors) => product(ev, factors, t, s),
        (Some(Symbol::POWER), [b, n]) if b == t && !n.contains(t) => {
            // Γ(n + 1)/s^(n+1) for n > −1, or a symbol.
            let n1 = ev.eval(&(n.clone() + Expr::from(1)))?;
            let defined = n
                .to_rational()
                .is_none_or(|q| q > cas_num::Rational::from(cas_num::Integer::from(-1)));
            if !defined {
                return Ok(None);
            }
            Ok(Some(Expr::call(Symbol::GAMMA, [n1.clone()]) * s.pow(-n1)))
        }
        _ if trigonometric(f) && head == Some(Symbol::POWER) => reduced(ev, f, t, s),
        _ => {
            if let Some(u) = exponent(f) {
                let Some((a, b)) = linear(ev, u, t)? else {
                    return Ok(None);
                };
                return Ok(Some(Expr::from(Symbol::E).pow(b) / (s.clone() - a)));
            }
            if let Some(image) = elementary(ev, f, t, s)? {
                return Ok(Some(image));
            }
            // Last, the expanded form: (t + 1)², sums under a product.
            let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone()]))?;
            if expanded == *f {
                return Ok(None);
            }
            transform(ev, &expanded, t, s)
        }
    }
}

/// The transform of a sine, cosine, hyperbolic sine or cosine, step, delta,
/// or unknown function and its derivatives.
fn elementary(ev: &mut Eval<'_>, f: &Expr, t: &Expr, s: &Expr) -> Result<Option<Expr>, EvalError> {
    let View::Apply { head, args } = f.view() else {
        return Ok(None);
    };
    let square = |b: &Expr| b.pow(2);
    match (head.as_symbol(), args) {
        (Some(g @ (Symbol::SIN | Symbol::COS | Symbol::SINH | Symbol::COSH)), [u]) => {
            let Some((b, c)) = linear(ev, u, t)? else {
                return Ok(None);
            };
            if c != Expr::from(0) {
                // A phase: sin(b t + c) = sin(b t) cos c + cos(b t) sin c.
                let expanded = if matches!(g, Symbol::SIN | Symbol::COS) {
                    Expr::call(Symbol::TRIG_EXPAND, [f.clone()])
                } else {
                    Expr::call(Symbol::TRIG_TO_EXP, [f.clone()])
                };
                let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [expanded]))?;
                if expanded == *f {
                    return Ok(None);
                }
                return transform(ev, &expanded, t, s);
            }
            let s2 = square(s);
            Ok(Some(match g {
                Symbol::SIN => b.clone() / (s2 + square(&b)),
                Symbol::COS => s.clone() / (s2 + square(&b)),
                Symbol::SINH => b.clone() / (s2 - square(&b)),
                _ => s.clone() / (s2 - square(&b)),
            }))
        }
        // The Bessel functions of the first kind (CAS-DOC-34 §5):
        // `(√(s² + b²) − s)^n/(b^n √(s² + b²))` and `(s − √(s² −
        // b²))^n/(b^n √(s² − b²))`.
        (Some(g @ (Symbol::BESSEL_J | Symbol::BESSEL_I)), [n, u]) => {
            let Some(n) = n.as_i64().filter(|&n| n >= 0) else {
                return Ok(None);
            };
            let Some((b, c)) = linear(ev, u, t)? else {
                return Ok(None);
            };
            if c != Expr::from(0) {
                return Ok(None);
            }
            let s2 = square(s);
            let (root, top) = if g == Symbol::BESSEL_J {
                let root = Expr::call(Symbol::SQRT, [s2 + square(&b)]);
                (root.clone(), root - s.clone())
            } else {
                let root = Expr::call(Symbol::SQRT, [s2 - square(&b)]);
                (root.clone(), s.clone() - root)
            };
            Ok(Some(top.pow(n) / (b.pow(n) * root)))
        }
        (Some(Symbol::HEAVISIDE_THETA), [u]) => {
            let Some(c) = step_at(ev, u, t)? else {
                return Ok(None);
            };
            Ok(match real_sign(ev.context(), &c) {
                Some(Ordering::Greater | Ordering::Equal) => {
                    Some(Expr::from(Symbol::E).pow(-(c * s.clone())) / s.clone())
                }
                Some(Ordering::Less) => Some(s.pow(-1)),
                None => None,
            })
        }
        (Some(Symbol::DIRAC_DELTA), [u]) => {
            let Some((a, b)) = linear(ev, u, t)? else {
                return Ok(None);
            };
            let c = ev.eval(&(-(b / a.clone())))?;
            let scale = Expr::call(Symbol::ABS, [a]);
            Ok(match real_sign(ev.context(), &c) {
                Some(Ordering::Greater | Ordering::Equal) => {
                    Some(Expr::from(Symbol::E).pow(-(c * s.clone())) / scale)
                }
                Some(Ordering::Less) => Some(Expr::from(0)),
                None => None,
            })
        }
        // An unknown function of t, and its derivatives.
        (Some(g), [u]) if !g.is_builtin() && u == t => Ok(Some(Expr::call(
            Symbol::LAPLACE_TRANSFORM,
            [f.clone(), t.clone(), s.clone()],
        ))),
        (None, [u]) if u == t => derivative(ev, head, t, s),
        _ => Ok(None),
    }
}

/// `L[f⁽ⁿ⁾(t)] = sⁿ F − Σ_{k<n} s^(n−1−k) f⁽ᵏ⁾(0)` for the head
/// `Derivative[n][f]`, `f` an unknown function.
fn derivative(
    ev: &mut Eval<'_>,
    head: &Expr,
    t: &Expr,
    s: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let View::Apply { head: d, args: [f] } = head.view() else {
        return Ok(None);
    };
    let View::Apply { head: d, args: [n] } = d.view() else {
        return Ok(None);
    };
    if d.as_symbol() != Some(Symbol::DERIVATIVE) || f.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    let Some(n) = n.as_i64().filter(|&n| n >= 1) else {
        return Ok(None);
    };
    let image = Expr::call(
        Symbol::LAPLACE_TRANSFORM,
        [
            Expr::apply(f, std::slice::from_ref(t)),
            t.clone(),
            s.clone(),
        ],
    );
    let mut terms = vec![s.pow(n) * image];
    for k in 0..n {
        let at_zero = if k == 0 {
            Expr::apply(f, &[Expr::from(0)])
        } else {
            let dk = Expr::call(Symbol::DERIVATIVE, [Expr::from(k)]);
            Expr::apply(&Expr::apply(&dk, std::slice::from_ref(f)), &[Expr::from(0)])
        };
        terms.push(-(s.pow(n - 1 - k) * at_zero));
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// `f` reduced (`TrigReduce`) and expanded, then transformed; `None` when
/// that changes nothing.
fn reduced(ev: &mut Eval<'_>, f: &Expr, t: &Expr, s: &Expr) -> Result<Option<Expr>, EvalError> {
    let reduced = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [Expr::call(Symbol::TRIG_REDUCE, [f.clone()])],
    ))?;
    if reduced == *f {
        return Ok(None);
    }
    transform(ev, &reduced, t, s)
}

/// The transform of a product: factors free of `t` out, then a shift by an
/// exponential or a step, a power of `t`, or products of sines and cosines
/// reduced.
fn product(
    ev: &mut Eval<'_>,
    factors: &[Expr],
    t: &Expr,
    s: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let (constant, varying): (Vec<Expr>, Vec<Expr>) =
        factors.iter().cloned().partition(|g| !g.contains(t));
    let times = |v: Vec<Expr>| match v.len() {
        1 => v.into_iter().next().expect("one"),
        _ => Expr::call(Symbol::TIMES, v),
    };
    if !constant.is_empty() {
        let rest = times(varying);
        return Ok(transform(ev, &rest, t, s)?.map(|g| times(constant) * g));
    }
    // A shift by an exponential: e^(a t + b) g(t) is e^b G(s − a).
    if let Some(i) = varying.iter().position(|g| exponent(g).is_some()) {
        let u = exponent(&varying[i]).expect("an exponential").clone();
        if let Some((a, b)) = linear(ev, &u, t)? {
            let mut rest = varying.clone();
            rest.remove(i);
            let Some(g) = transform(ev, &times(rest), t, s)? else {
                return Ok(None);
            };
            let shifted = g.substitute(&[(s.clone(), s.clone() - a)]);
            return Ok(Some(Expr::from(Symbol::E).pow(b) * shifted));
        }
    }
    // A shift by a step: θ(t − c) g(t) is e^(−c s) L[g(t + c)] for c ≥ 0.
    let step = Expr::from(Symbol::HEAVISIDE_THETA);
    if let Some(i) = varying.iter().position(|g| g.head() == Some(&step)) {
        let [u] = varying[i].args() else {
            return Ok(None);
        };
        let Some(c) = step_at(ev, u, t)? else {
            return Ok(None);
        };
        let mut rest = varying.clone();
        rest.remove(i);
        let rest = times(rest);
        return match real_sign(ev.context(), &c) {
            Some(Ordering::Greater | Ordering::Equal) => {
                let moved = ev.eval(&rest.substitute(&[(t.clone(), t.clone() + c.clone())]))?;
                let Some(g) = transform(ev, &moved, t, s)? else {
                    return Ok(None);
                };
                Ok(Some(Expr::from(Symbol::E).pow(-(c * s.clone())) * g))
            }
            // The step is 1 for every t ≥ 0.
            Some(Ordering::Less) => transform(ev, &rest, t, s),
            None => Ok(None),
        };
    }
    // A power of t: tⁿ g(t) is (−1)ⁿ G⁽ⁿ⁾(s).
    let power_of_t = |g: &Expr| -> Option<i64> {
        if g == t {
            return Some(1);
        }
        match (g.head().and_then(Expr::as_symbol), g.args()) {
            (Some(Symbol::POWER), [b, n]) if b == t => n.as_i64().filter(|&n| n >= 1),
            _ => None,
        }
    };
    if let Some(i) = varying.iter().position(|g| power_of_t(g).is_some()) {
        let n = power_of_t(&varying[i]).expect("a power of t");
        let mut rest = varying.clone();
        rest.remove(i);
        let Some(g) = transform(ev, &times(rest), t, s)? else {
            return Ok(None);
        };
        let order = Expr::call(Symbol::LIST, [s.clone(), Expr::from(n)]);
        let derivative = ev.eval(&Expr::call(Symbol::D, [g, order]))?;
        let sign = if n % 2 == 0 { 1 } else { -1 };
        return Ok(Some(Expr::from(sign) * derivative));
    }
    // Products of sines and cosines.
    if varying.iter().filter(|g| trigonometric(g)).count() >= 2 {
        return reduced(ev, &times(varying), t, s);
    }
    // A sum among the factors.
    let f = times(varying);
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone()]))?;
    if expanded == f {
        return Ok(None);
    }
    transform(ev, &expanded, t, s)
}

mod inverse;

pub(super) use inverse::inverse_laplace_transform;

#[cfg(test)]
mod tests;
