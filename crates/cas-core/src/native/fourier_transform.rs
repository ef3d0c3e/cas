//! Fourier transforms (CAS-DOC-41): `FourierTransform[f, x, s]`, Xcas's
//! `fourier(f, x, s)`, `F(s) = ∫ f(x) e^(−i s x) dx` over the real line,
//! giac's convention.
//!
//! By the first rule that applies, recursively (§2, CAS-439):
//!
//! - **Linearity**: sums term by term, factors free of `x` out; `c` free of
//!   `x` is `2π c δ(s)`.
//! - **Deltas**: `g(x) δ(a x + b)` is `g(x₀) e^(−i s x₀)/|a|`, `x₀ = −b/a`.
//! - **Powers of `x`**: `xⁿ g(x)`, `n` a positive integer, is `iⁿ G⁽ⁿ⁾(s)`.
//! - **Oscillations**: a product with sines or cosines of linear arguments
//!   is written with exponentials (`TrigToExp`) and expanded; an
//!   exponential `e^(i b x)` shifts the rest's transform to `s − b` in the
//!   table.
//! - **Unknown functions**: `g(x)` stays as `FourierTransform[g[x], x, s]`,
//!   giac's `Fourier(f(x),x,s)`, and `g⁽ⁿ⁾(x)` is `(i s)ⁿ` times it.
//! - **The table** (§3, CAS-440), on a product of exponentials, their
//!   exponents summed into `w`, and at most one other factor, at `σ = s −
//!   b` for the modulation `e^(i b x)`:
//!   - `w` a polynomial `−A x² + B x + C`, `A` real and positive:
//!     `e^C √(π/A) e^((B − i s)²/(4A))`, the Gaussians shifted, scaled and
//!     modulated at once; `A = 0` and `B = i b`: `2π e^C δ(s − b)`;
//!   - `w = −a |p x + q| + i b x + C`, `a > 0`: `e^C (1/|p|) e^(i σ q/p)
//!     2a/(a² + (σ/p)²)`;
//!   - the step `θ(p x + q)` with `w = B x + C` linear: for `Re B < 0` and
//!     `p > 0`, `−e^C e^((B − i s) c)/(B − i s)` from `c = −q/p` on, and
//!     for `Re B > 0` and `p < 0`, `e^C e^((B − i s) c)/(B − i s)` up to
//!     `c`; for `Re B = 0`, `B = i b`, the distribution `π δ(σ) ∓ i
//!     e^(−i σ c)/σ` (`θ(x)` is giac's `pi*Dirac(s)+(-i)/s`);
//!   - `sign(p x + q)`: `sign(p) e^(i σ q/p) (−2i/σ)`; `1/(p x + q)`, the
//!     principal value: `(1/p) e^(i σ q/p) (−i π sign σ)`;
//!   - the box `UnitBox[p x + q]`, giac's `rect`: `(1/|p|) e^(i σ q/p) 2
//!     sin(σ/(2p))/(σ/p)`.
//!
//! What no rule transforms stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the transform: f in x into s, w = −A x² + B x + C"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

type Found = Result<Option<Expr>, EvalError>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn call(h: Symbol, u: Expr) -> Expr {
    Expr::call(h, [u])
}

fn i() -> Expr {
    Expr::from(Symbol::I)
}

fn pi() -> Expr {
    Expr::from(Symbol::PI)
}

fn delta(u: Expr) -> Expr {
    call(Symbol::DIRAC_DELTA, u)
}

fn exp(u: Expr) -> Expr {
    Expr::from(Symbol::E).pow(u)
}

/// `FourierTransform[f, x, s]`.
pub(super) fn fourier_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, x, s] = e.args() else {
        return Ok(None);
    };
    let user = |v: &Expr| v.as_symbol().is_some_and(|v| !v.is_builtin());
    if !user(x) || !user(s) || x == s {
        return Ok(None);
    }
    // Into a fresh variable, put in for s after.
    let fresh = Expr::from(Symbol::new("FourierTransform$s"));
    let Some(image) = transform(ev, f, x, &fresh)? else {
        return Ok(None);
    };
    let image = image.substitute(&[(fresh, s.clone())]);
    // An unknown function's transform is the call itself, which stays.
    if image == *e {
        return Ok(None);
    }
    let result = ev.eval(&image)?;
    ev.step("FourierTransform", e, &result);
    Ok(Some(result))
}

/// `(a, b)` with `u = a x + b`, both free of `x`.
fn linear(ev: &mut Eval<'_>, u: &Expr, x: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    let a = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
    if a.contains(x) {
        return Ok(None);
    }
    let b = ev.eval(&u.substitute(&[(x.clone(), Expr::from(0))]))?;
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - a.clone() * x.clone() - b.clone()],
    ))?;
    Ok((gap == Expr::from(0)).then_some((a, b)))
}

/// The exponent `u` of `e^u` or `Exp[u]`.
fn exponent(e: &Expr) -> Option<&Expr> {
    match (head(e), e.args()) {
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => Some(u),
        (Some(Symbol::EXP), [u]) => Some(u),
        _ => None,
    }
}

/// The transform of `f` in `x` into `s`; `None` when no rule applies.
fn transform(ev: &mut Eval<'_>, f: &Expr, x: &Expr, s: &Expr) -> Found {
    if !f.contains(x) {
        return Ok(Some(Expr::from(2) * pi() * f.clone() * delta(s.clone())));
    }
    match head(f) {
        Some(Symbol::PLUS) => {
            let mut images = Vec::with_capacity(f.args().len());
            for term in f.args() {
                let Some(image) = transform(ev, term, x, s)? else {
                    return Ok(None);
                };
                images.push(image);
            }
            Ok(Some(Expr::call(Symbol::PLUS, images)))
        }
        Some(Symbol::TIMES) => product(ev, f.args(), x, s),
        _ => product(ev, std::slice::from_ref(f), x, s),
    }
}

/// The positive integer `n` of a factor `x` or `xⁿ`.
fn power_of(g: &Expr, x: &Expr) -> Option<i64> {
    if g == x {
        return Some(1);
    }
    match (head(g), g.args()) {
        (Some(Symbol::POWER), [b, n]) if b == x => n.as_i64().filter(|&n| n >= 1),
        _ => None,
    }
}

/// The transform of the product of `factors`.
fn product(ev: &mut Eval<'_>, factors: &[Expr], x: &Expr, s: &Expr) -> Found {
    let (constants, varying): (Vec<Expr>, Vec<Expr>) =
        factors.iter().cloned().partition(|g| !g.contains(x));
    if !constants.is_empty() {
        let rest: Expr = varying.into_iter().product();
        let constant: Expr = constants.into_iter().product();
        return Ok(transform(ev, &rest, x, s)?.map(|image| constant * image));
    }
    let others = |k: usize| -> Expr {
        varying
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != k)
            .map(|(_, g)| g.clone())
            .product()
    };
    // A delta sifts the rest out.
    for (k, g) in varying.iter().enumerate() {
        if let (Some(Symbol::DIRAC_DELTA), [u]) = (head(g), g.args())
            && let Some((a, b)) = linear(ev, u, x)?
        {
            let x0 = ev.eval(&(-(b / a.clone())))?;
            let rest = others(k).substitute(&[(x.clone(), x0.clone())]);
            let image = rest * exp(-(i() * s.clone() * x0)) / call(Symbol::ABS, a);
            return Ok(Some(image));
        }
    }
    // xⁿ g(x): iⁿ G⁽ⁿ⁾(s).
    for (k, g) in varying.iter().enumerate() {
        if let Some(n) = power_of(g, x) {
            let Some(image) = transform(ev, &others(k), x, s)? else {
                return Ok(None);
            };
            let derivative = Expr::call(
                Symbol::D,
                [image, Expr::call(Symbol::LIST, [s.clone(), Expr::from(n)])],
            );
            return Ok(Some(i().pow(n) * ev.eval(&derivative)?));
        }
    }
    let f: Expr = varying.iter().cloned().product();
    // Oscillations as exponentials.
    if varying
        .iter()
        .any(|g| matches!(head(g), Some(Symbol::SIN | Symbol::COS)) || trigonometric_power(g))
    {
        let rewritten = ev.eval(&Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::TRIG_TO_EXP, [f.clone()])],
        ))?;
        if rewritten != f {
            return transform(ev, &rewritten, x, s);
        }
    }
    table(ev, &f, x, s)
}

/// A power of a sine or cosine.
fn trigonometric_power(g: &Expr) -> bool {
    matches!((head(g), g.args()), (Some(Symbol::POWER), [b, _]) if matches!(head(b), Some(Symbol::SIN | Symbol::COS)))
}

/// The real and imaginary parts of `b`, evaluated.
fn parts(ev: &mut Eval<'_>, b: &Expr) -> Result<(Expr, Expr), EvalError> {
    let re = ev.eval(&call(Symbol::RE, b.clone()))?;
    let im = ev.eval(&call(Symbol::IM, b.clone()))?;
    Ok((re, im))
}

/// The table on `f`: its exponentials, exponents summed, and at most one
/// other factor.
fn table(ev: &mut Eval<'_>, f: &Expr, x: &Expr, s: &Expr) -> Found {
    let factors = match head(f) {
        Some(Symbol::TIMES) => f.args().to_vec(),
        _ => vec![f.clone()],
    };
    let mut w = Vec::new();
    let mut rest = Vec::new();
    for g in factors {
        match exponent(&g) {
            Some(u) => w.push(u.clone()),
            None => rest.push(g),
        }
    }
    let w = ev.eval(&Expr::call(Symbol::PLUS, w))?;
    match rest.as_slice() {
        [] => exponential(ev, &w, x, s),
        [h] => with_factor(ev, h, &w, x, s),
        _ => Ok(None),
    }
}

/// `e^w`: a Gaussian, an oscillation, or `e^(−a|p x + q|)`, modulated.
fn exponential(ev: &mut Eval<'_>, w: &Expr, x: &Expr, s: &Expr) -> Found {
    let absolute: Vec<Expr> = w
        .distinct()
        .filter(|g| head(g) == Some(Symbol::ABS) && g.contains(x))
        .cloned()
        .collect();
    if let [absolute] = absolute.as_slice() {
        return two_sided(ev, w, absolute, x, s);
    }
    if !absolute.is_empty() {
        return Ok(None);
    }
    // w = −A x² + B x + C.
    let second = ev.eval(&Expr::call(
        Symbol::D,
        [
            w.clone(),
            Expr::call(Symbol::LIST, [x.clone(), Expr::from(2)]),
        ],
    ))?;
    if second.contains(x) {
        return Ok(None);
    }
    let a = ev.eval(&(-(second / Expr::from(2))))?;
    let first = ev.eval(&Expr::call(Symbol::D, [w.clone(), x.clone()]))?;
    let b = ev.eval(&first.substitute(&[(x.clone(), Expr::from(0))]))?;
    let c = ev.eval(&w.substitute(&[(x.clone(), Expr::from(0))]))?;
    let rebuilt = -(a.clone()) * x.pow(2) + b.clone() * x.clone() + c.clone();
    let gap = ev.eval(&Expr::call(Symbol::EXPAND, [w.clone() - rebuilt]))?;
    if gap != Expr::from(0) {
        return Ok(None);
    }
    if a == Expr::from(0) {
        // e^(i b x + C): 2π e^C δ(s − b).
        let (re, im) = parts(ev, &b)?;
        if re != Expr::from(0) {
            return Ok(None);
        }
        return Ok(Some(Expr::from(2) * pi() * exp(c) * delta(s.clone() - im)));
    }
    if real_sign(ev.context(), &a) != Some(Ordering::Greater) {
        return Ok(None);
    }
    let shifted = b - i() * s.clone();
    let half = Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
    Ok(Some(
        exp(c) * (pi() / a.clone()).pow(half) * exp(shifted.pow(2) / (Expr::from(4) * a)),
    ))
}

/// `e^w` with `w = α |p x + q| + i b x + C`, `α < 0`.
fn two_sided(ev: &mut Eval<'_>, w: &Expr, absolute: &Expr, x: &Expr, s: &Expr) -> Found {
    let [u] = absolute.args() else {
        return Ok(None);
    };
    let Some((p, q)) = linear(ev, u, x)? else {
        return Ok(None);
    };
    let t = Expr::from(Symbol::new("FourierTransform$t"));
    let in_t = ev.eval(&w.substitute(&[(absolute.clone(), t.clone())]))?;
    let alpha = ev.eval(&Expr::call(Symbol::D, [in_t.clone(), t.clone()]))?;
    if alpha.contains(x) || alpha.contains(&t) {
        return Ok(None);
    }
    if real_sign(ev.context(), &alpha) != Some(Ordering::Less) {
        return Ok(None);
    }
    let remainder = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [in_t - alpha.clone() * t.clone()],
    ))?;
    let Some((slope, c)) = linear(ev, &remainder, x)? else {
        return Ok(None);
    };
    let (re, b) = parts(ev, &slope)?;
    if re != Expr::from(0) {
        return Ok(None);
    }
    let a = -alpha;
    let sigma = s.clone() - b;
    let scaled = sigma.clone() / p.clone();
    Ok(Some(
        exp(c) / call(Symbol::ABS, p.clone())
            * exp(i() * sigma * q / p)
            * Expr::from(2)
            * a.clone()
            / (a.pow(2) + scaled.pow(2)),
    ))
}

/// `e^w h(x)` for a step, a sign, a box, `1/(p x + q)` or an unknown
/// function `h`.
fn with_factor(ev: &mut Eval<'_>, h: &Expr, w: &Expr, x: &Expr, s: &Expr) -> Found {
    // An unknown function and its derivatives, unmodulated.
    if *w == Expr::from(0)
        && let Some(image) = unknown(h, x, s)
    {
        return Ok(Some(image));
    }
    // w = B x + C.
    let Some((b, c)) = linear(ev, w, x)? else {
        return Ok(None);
    };
    let (re, im) = parts(ev, &b)?;
    let sigma = s.clone() - im;
    let (Some(g), [u]) = (head(h), h.args()) else {
        // 1/(p x + q), the principal value.
        if let (Some(Symbol::POWER), [base, n]) = (head(h), h.args())
            && *n == Expr::from(-1)
            && re == Expr::from(0)
            && let Some((p, q)) = linear(ev, base, x)?
        {
            return Ok(Some(
                exp(c) / p.clone()
                    * exp(i() * sigma.clone() * q / p)
                    * -(i() * pi() * call(Symbol::SIGN, sigma)),
            ));
        }
        return Ok(None);
    };
    let Some((p, q)) = linear(ev, u, x)? else {
        return Ok(None);
    };
    let shift = exp(i() * sigma.clone() * q.clone() / p.clone());
    match g {
        Symbol::HEAVISIDE_THETA => step(ev, &b, &c, &re, &sigma, &p, &q, s),
        Symbol::SIGN if re == Expr::from(0) => Ok(Some(
            exp(c) * call(Symbol::SIGN, p) * shift * Expr::from(-2) * i() / sigma,
        )),
        Symbol::UNIT_BOX if re == Expr::from(0) => {
            let scaled = sigma / p.clone();
            Ok(Some(
                exp(c) / call(Symbol::ABS, p)
                    * shift
                    * Expr::from(2)
                    * call(Symbol::SIN, scaled.clone() / Expr::from(2))
                    / scaled,
            ))
        }
        _ => Ok(None),
    }
}

/// `e^(B x + C) θ(p x + q)`.
#[allow(clippy::too_many_arguments, reason = "the parts of one formula")]
fn step(
    ev: &mut Eval<'_>,
    b: &Expr,
    c: &Expr,
    re: &Expr,
    sigma: &Expr,
    p: &Expr,
    q: &Expr,
    s: &Expr,
) -> Found {
    let Some(direction) = real_sign(ev.context(), p) else {
        return Ok(None);
    };
    // The step at x = cut: from it on for p > 0, up to it for p < 0.
    let cut = ev.eval(&(-(q.clone() / p.clone())))?;
    let up = direction == Ordering::Greater;
    let exponent = b.clone() - i() * s.clone();
    match (real_sign(ev.context(), re), up) {
        (Some(Ordering::Less), true) => Ok(Some(
            -exp(c.clone()) * exp(exponent.clone() * cut) / exponent,
        )),
        (Some(Ordering::Greater), false) => Ok(Some(
            exp(c.clone()) * exp(exponent.clone() * cut) / exponent,
        )),
        (Some(Ordering::Equal), _) => {
            // π δ(σ) ∓ i e^(−i σ cut)/σ.
            let jump = i() * exp(-(i() * sigma.clone() * cut)) / sigma.clone();
            let image = if up {
                pi() * delta(sigma.clone()) - jump
            } else {
                pi() * delta(sigma.clone()) + jump
            };
            Ok(Some(exp(c.clone()) * image))
        }
        _ => Ok(None),
    }
}

/// `g(x)` and `g⁽ⁿ⁾(x)` for an unknown function `g`.
fn unknown(h: &Expr, x: &Expr, s: &Expr) -> Option<Expr> {
    let View::Apply { head: g, args: [u] } = h.view() else {
        return None;
    };
    if u != x {
        return None;
    }
    if g.as_symbol().is_some_and(|g| !g.is_builtin()) {
        return Some(Expr::call(
            Symbol::FOURIER_TRANSFORM,
            [h.clone(), x.clone(), s.clone()],
        ));
    }
    // Derivative[n][g][x]: (i s)ⁿ times g's transform.
    let View::Apply {
        head: d,
        args: [function],
    } = g.view()
    else {
        return None;
    };
    let View::Apply { head: d, args: [n] } = d.view() else {
        return None;
    };
    if d.as_symbol() != Some(Symbol::DERIVATIVE)
        || function.as_symbol().is_none_or(Symbol::is_builtin)
    {
        return None;
    }
    let n = n.as_i64().filter(|&n| n >= 1)?;
    let image = Expr::call(
        Symbol::FOURIER_TRANSFORM,
        [
            Expr::apply(function, std::slice::from_ref(x)),
            x.clone(),
            s.clone(),
        ],
    );
    Some((i() * s.clone()).pow(n) * image)
}

#[cfg(test)]
mod tests;
