//! Z-transforms (CAS-DOC-24 §2, CAS-308): `ZTransform[f, n, z]`, Xcas's
//! `ztrans(f, n, z)`, the one-sided `Σ_{n≥0} f(n) z⁻ⁿ`, and
//! `KroneckerDelta[n]`, Xcas's `Kronecker(n)`: 1 at 0, 0 at the other
//! numbers, a symbol's staying.
//!
//! By the first rule that applies, recursively, as `LaplaceTransform`'s:
//!
//! - **Linearity**: sums term by term, factors free of `n` out; `c` free of
//!   `n` is `c z/(z − 1)`.
//! - **Table**: `aⁿ` is `z/(z − a)` (`a^(αn+β)` and `e^(αn+β)` included);
//!   `sin(ωn)` is `z sin ω/(z² − 2z cos ω + 1)` and `cos(ωn)` is `z(z − cos
//!   ω)/(z² − 2z cos ω + 1)`, a phase expanded first; `1/n!` is `e^(1/z)`;
//!   `C(n, k)` for a natural `k` is `z/(z − 1)^(k+1)`; `KroneckerDelta[n − m]`
//!   is `z⁻ᵐ` for a natural `m` and 0 for a negative integer;
//!   `HeavisideTheta[n − m]` for an integer `m ≥ 0` is `z⁻ᵐ/2 + z⁻ᵐ/(z −
//!   1)`, this CAS's step being 1/2 at 0 (giac's sequences take 1).
//! - **Products**: `aⁿ g(n)` is `G(z/a)`; `n g(n)` is `−z G′(z)`;
//!   `KroneckerDelta[n − m] g(n)` is `g(m) z⁻ᵐ`; `HeavisideTheta[n − m]
//!   g(n)` is `z⁻ᵐ (Z[g(n + m)] − g(m)/2)`; products and powers of sines and
//!   cosines reduced (`TrigReduce`).
//! - **Unknown sequences**: `u(n)` stays as `ZTransform[u[n], n, z]`, `U`,
//!   and `u(n + k)` for a natural `k` is `zᵏ (U − Σ_{j<k} u(j) z⁻ʲ)`, which
//!   makes recurrences equations in `U` (CAS-DOC-24 §4).
//!
//! The answer is `Together`'s, as giac's `z/(z^2-2*z+1)`; what no rule
//! transforms stays, where giac answers an unevaluated sum. The inverse is
//! [`inverse`]'s (CAS-309).

#![allow(
    clippy::many_single_char_names,
    reason = "the transform's names: f(n), F(z), aⁿ, sin(ωn), n − m"
)]

use std::cmp::Ordering;

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `KroneckerDelta[x]`: 1 at 0, 0 at another number, as comparisons decide.
pub(super) fn kronecker_delta(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if *x == Expr::from(0) {
        return Ok(Some(Expr::from(1)));
    }
    if !x.is_number() {
        return Ok(None);
    }
    Ok(match real_sign(ev.context(), x) {
        Some(Ordering::Greater | Ordering::Less) => Some(Expr::from(0)),
        Some(Ordering::Equal) => Some(Expr::from(1)),
        None => None,
    })
}

/// `ZTransform[f, n, z]`.
pub(super) fn z_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, n, z] = e.args() else {
        return Ok(None);
    };
    if n.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    // Transformed into a fresh variable, put in for z after: z may be n.
    let fresh = Expr::from(Symbol::new("ZTransform$z"));
    let Some(image) = transform(ev, f, n, &fresh)? else {
        return Ok(None);
    };
    let image = image.substitute(&[(fresh, z.clone())]);
    // An unknown sequence's transform is the call itself, which stays:
    // evaluated, it would come back here.
    if image == *e {
        return Ok(None);
    }
    let result = ev.eval(&Expr::call(Symbol::TOGETHER, [image]))?;
    ev.step("ZTransform", e, &result);
    Ok(Some(result))
}

/// `(a, b)` with `u = a n + b`, both free of `n`.
fn linear(ev: &mut Eval<'_>, u: &Expr, n: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    let a = ev.eval(&Expr::call(Symbol::D, [u.clone(), n.clone()]))?;
    if a.contains(n) {
        return Ok(None);
    }
    let b = ev.eval(&u.substitute(&[(n.clone(), Expr::from(0))]))?;
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - a.clone() * n.clone() - b.clone()],
    ))?;
    Ok((gap == Expr::from(0)).then_some((a, b)))
}

/// The shift `m` of `u = n − m`, an integer; `None` otherwise.
fn shift(ev: &mut Eval<'_>, u: &Expr, n: &Expr) -> Result<Option<i64>, EvalError> {
    let Some((a, b)) = linear(ev, u, n)? else {
        return Ok(None);
    };
    if a != Expr::from(1) {
        return Ok(None);
    }
    Ok(b.as_i64().and_then(i64::checked_neg))
}

/// `(base, exponent)` of `bᵘ`, `e^u` as `(E, u)`.
fn power_of(e: &Expr) -> Option<(Expr, &Expr)> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::POWER), [b, u]) => Some((b.clone(), u)),
        (Some(Symbol::EXP), [u]) => Some((Expr::from(Symbol::E), u)),
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

/// `aⁿ` with `a` free of `n`, from `b^(αn + β)` or `e^(αn + β)`: `(a,
/// b^β)`.
fn exponential(ev: &mut Eval<'_>, e: &Expr, n: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    let Some((b, u)) = power_of(e) else {
        return Ok(None);
    };
    if b.contains(n) || !u.contains(n) {
        return Ok(None);
    }
    let Some((alpha, beta)) = linear(ev, u, n)? else {
        return Ok(None);
    };
    let a = ev.eval(&b.pow(alpha))?;
    let scale = ev.eval(&b.pow(beta))?;
    Ok(Some((a, scale)))
}

/// `nᵏ` for a natural `k ≥ 1`.
fn power_of_n(g: &Expr, n: &Expr) -> Option<i64> {
    if g == n {
        return Some(1);
    }
    match (g.head().and_then(Expr::as_symbol), g.args()) {
        (Some(Symbol::POWER), [b, k]) if b == n => k.as_i64().filter(|&k| k >= 1),
        _ => None,
    }
}

/// `(−z d/dz)ᵏ G`.
fn times_n(ev: &mut Eval<'_>, g: Expr, z: &Expr, k: i64) -> Result<Expr, EvalError> {
    let mut image = g;
    for _ in 0..k {
        let derivative = ev.eval(&Expr::call(Symbol::D, [image, z.clone()]))?;
        image = ev.eval(&(-(z * derivative)))?;
    }
    Ok(image)
}

/// The transform of `f` in `n` into `z`; `None` when no rule applies.
fn transform(ev: &mut Eval<'_>, f: &Expr, n: &Expr, z: &Expr) -> Result<Option<Expr>, EvalError> {
    if !f.contains(n) {
        return Ok(Some(f.clone() * z / (z - 1)));
    }
    let (head, args) = (f.head().and_then(Expr::as_symbol), f.args());
    match (head, args) {
        (Some(Symbol::PLUS), terms) => {
            let mut images = Vec::with_capacity(terms.len());
            for term in terms {
                let Some(image) = transform(ev, term, n, z)? else {
                    return Ok(None);
                };
                images.push(image);
            }
            Ok(Some(Expr::call(Symbol::PLUS, images)))
        }
        (Some(Symbol::TIMES), factors) => product(ev, factors, n, z),
        _ if power_of_n(f, n).is_some() => {
            let k = power_of_n(f, n).expect("a power of n");
            let one = z / (z - 1);
            times_n(ev, one, z, k).map(Some)
        }
        _ if trigonometric(f) && head == Some(Symbol::POWER) => reduced(ev, f, n, z),
        _ => {
            if let Some((a, scale)) = exponential(ev, f, n)? {
                return Ok(Some(scale * z / (z - a)));
            }
            if let Some(image) = elementary(ev, f, n, z)? {
                return Ok(Some(image));
            }
            let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone()]))?;
            if expanded == *f {
                return Ok(None);
            }
            transform(ev, &expanded, n, z)
        }
    }
}

/// The transform of a sine, cosine, reciprocal factorial, binomial,
/// Kronecker delta, step, or unknown sequence and its shifts.
fn elementary(ev: &mut Eval<'_>, f: &Expr, n: &Expr, z: &Expr) -> Result<Option<Expr>, EvalError> {
    let View::Apply { head, args } = f.view() else {
        return Ok(None);
    };
    let square = |b: &Expr| b.pow(2);
    match (head.as_symbol(), args) {
        (Some(g @ (Symbol::SIN | Symbol::COS)), [u]) => {
            let Some((w, phase)) = linear(ev, u, n)? else {
                return Ok(None);
            };
            if phase != Expr::from(0) {
                let expanded = ev.eval(&Expr::call(
                    Symbol::EXPAND,
                    [Expr::call(Symbol::TRIG_EXPAND, [f.clone()])],
                ))?;
                if expanded == *f {
                    return Ok(None);
                }
                return transform(ev, &expanded, n, z);
            }
            let (sin, cos) = (
                Expr::call(Symbol::SIN, [w.clone()]),
                Expr::call(Symbol::COS, [w]),
            );
            let denominator = square(z) - 2 * z * &cos + 1;
            Ok(Some(match g {
                Symbol::SIN => z * sin / denominator,
                _ => z * (z - cos) / denominator,
            }))
        }
        // 1/n!.
        (Some(Symbol::POWER), [b, k])
            if *k == Expr::from(-1)
                && b.head() == Some(&Expr::from(Symbol::FACTORIAL))
                && b.args() == std::slice::from_ref(n) =>
        {
            Ok(Some(Expr::call(Symbol::EXP, [z.pow(-1)])))
        }
        (Some(Symbol::BINOMIAL), [m, k]) if m == n && k.as_i64().is_some_and(|k| k >= 0) => {
            Ok(Some(z / (z - 1).pow(k.clone() + 1)))
        }
        (Some(Symbol::KRONECKER_DELTA), [u]) => {
            let Some(m) = shift(ev, u, n)? else {
                return Ok(None);
            };
            Ok(Some(if m >= 0 { z.pow(-m) } else { Expr::from(0) }))
        }
        (Some(Symbol::HEAVISIDE_THETA), [u]) => {
            let Some(m) = shift(ev, u, n)? else {
                return Ok(None);
            };
            if m < 0 {
                // 1 at every n ≥ 0.
                return Ok(Some(z / (z - 1)));
            }
            let at = z.pow(-m);
            Ok(Some(&at / 2 + at / (z - 1)))
        }
        // An unknown sequence and its shifts.
        (Some(g), [u]) if !g.is_builtin() => {
            if u == n {
                return Ok(Some(Expr::call(
                    Symbol::ZTRANSFORM,
                    [f.clone(), n.clone(), z.clone()],
                )));
            }
            let Some(m) = shift(ev, u, n)? else {
                return Ok(None);
            };
            // u(n + k), k = −m.
            let Some(k) = m.checked_neg().filter(|&k| k >= 1) else {
                return Ok(None);
            };
            let image = Expr::call(
                Symbol::ZTRANSFORM,
                [
                    Expr::apply(head, std::slice::from_ref(n)),
                    n.clone(),
                    z.clone(),
                ],
            );
            let mut terms = vec![image];
            for j in 0..k {
                terms.push(-(Expr::apply(head, &[Expr::from(j)]) * z.pow(-j)));
            }
            Ok(Some(z.pow(k) * Expr::call(Symbol::PLUS, terms)))
        }
        _ => Ok(None),
    }
}

/// `f` reduced (`TrigReduce`) and expanded, then transformed; `None` when
/// that changes nothing.
fn reduced(ev: &mut Eval<'_>, f: &Expr, n: &Expr, z: &Expr) -> Result<Option<Expr>, EvalError> {
    let reduced = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [Expr::call(Symbol::TRIG_REDUCE, [f.clone()])],
    ))?;
    if reduced == *f {
        return Ok(None);
    }
    transform(ev, &reduced, n, z)
}

/// The transform of a product: factors free of `n` out, then a scaling by
/// an exponential, a multiplication by a power of `n`, a Kronecker delta or
/// a step, or products of sines and cosines reduced.
fn product(
    ev: &mut Eval<'_>,
    factors: &[Expr],
    n: &Expr,
    z: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let (constant, varying): (Vec<Expr>, Vec<Expr>) =
        factors.iter().cloned().partition(|g| !g.contains(n));
    let times = |v: Vec<Expr>| match v.len() {
        0 => Expr::from(1),
        1 => v.into_iter().next().expect("one"),
        _ => Expr::call(Symbol::TIMES, v),
    };
    if !constant.is_empty() {
        let rest = times(varying);
        return Ok(transform(ev, &rest, n, z)?.map(|g| times(constant) * g));
    }
    let without = |i: usize| {
        let mut rest = varying.clone();
        rest.remove(i);
        times(rest)
    };
    // A scaling: aⁿ g(n) is G(z/a).
    for (i, g) in varying.iter().enumerate() {
        if let Some((a, scale)) = exponential(ev, g, n)? {
            let Some(image) = transform(ev, &without(i), n, z)? else {
                return Ok(None);
            };
            return Ok(Some(scale * image.substitute(&[(z.clone(), z / &a)])));
        }
    }
    // A multiplication by nᵏ: (−z d/dz)ᵏ G.
    if let Some(i) = varying.iter().position(|g| power_of_n(g, n).is_some()) {
        let k = power_of_n(&varying[i], n).expect("a power of n");
        let Some(image) = transform(ev, &without(i), n, z)? else {
            return Ok(None);
        };
        return times_n(ev, image, z, k).map(Some);
    }
    let kronecker = Expr::from(Symbol::KRONECKER_DELTA);
    let step = Expr::from(Symbol::HEAVISIDE_THETA);
    for (i, g) in varying.iter().enumerate() {
        let [u] = g.args() else {
            continue;
        };
        if g.head() != Some(&kronecker) && g.head() != Some(&step) {
            continue;
        }
        let Some(m) = shift(ev, u, n)? else {
            return Ok(None);
        };
        let rest = without(i);
        let at =
            |ev: &mut Eval<'_>, k: i64| ev.eval(&rest.substitute(&[(n.clone(), Expr::from(k))]));
        if g.head() == Some(&kronecker) {
            // Sifting: g(m) z⁻ᵐ for m ≥ 0.
            return Ok(Some(if m >= 0 {
                at(ev, m)? * z.pow(-m)
            } else {
                Expr::from(0)
            }));
        }
        if m < 0 {
            return transform(ev, &rest, n, z);
        }
        // z⁻ᵐ (Z[g(n + m)] − g(m)/2).
        let moved = ev.eval(&rest.substitute(&[(n.clone(), n + m)]))?;
        let Some(image) = transform(ev, &moved, n, z)? else {
            return Ok(None);
        };
        return Ok(Some(z.pow(-m) * (image - at(ev, m)? / 2)));
    }
    // Products of sines and cosines.
    if varying.iter().filter(|g| trigonometric(g)).count() >= 2 {
        return reduced(ev, &times(varying), n, z);
    }
    // A sum among the factors.
    let f = times(varying);
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone()]))?;
    if expanded == f {
        return Ok(None);
    }
    transform(ev, &expanded, n, z)
}

mod inverse;

pub(super) use inverse::inverse_z_transform;

#[cfg(test)]
mod tests;
