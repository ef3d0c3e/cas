//! Heuristics (CAS-DOC-9 §6, CAS-199), tried in order on what the table
//! and rational functions leave, each nesting the integration once more,
//! [`MAX_DEPTH`] deep at most:
//!
//! - an exponential times a sine or cosine of linear arguments, from the
//!   two-term system: with `u′ = α` and `v′ = γ`,
//!   `∫ e^u sin v = e^u (α sin v − γ cos v)/(α² + γ²)`, giac's
//!   `exp(x)*(-1/2*cos(x)+1/2*sin(x))`; at resonance, `α² + γ² = 0`, by
//!   `e^u sin v = (e^(u+iv) − e^(u−iv))/(2i)` and `e^u cos v = (e^(u+iv) +
//!   e^(u−iv))/2`, the exponential whose exponent is constant integrating
//!   to `x` times its value (CAS-DOC-25 §6, CAS-318): `∫ e^(ix) cos x =
//!   e^(2ix)/(4i) + x/2`;
//! - a polynomial `P` times an exponential, sine or cosine of a linear
//!   argument, by parts until `P` is gone, the tabular form:
//!   `∫ P e^u = e^u Σ (−1)^k P^(k)/α^(k+1)`;
//! - derivative-divides: for a part `g` of `f`, `f/g′` with a fresh `u`
//!   for `g` free of `x` is `q(u)`, and `∫ f` is `∫ q du` at `u = g`:
//!   `2x exp(x²)` is `exp(u)` in `u = x²`, `1/(x ln x)` is `1/u` in
//!   `u = ln x`, so `ln|ln x|`, where giac writes `ln(ln(x))`, complex for
//!   `0 < x < 1`;
//! - a polynomial (1 included) times a logarithm, a power of one, or an
//!   inverse tangent, sine or cosine of a linear argument, by parts:
//!   `∫ P L = Q L − ∫ Q L′` with `Q = ∫ P`;
//! - products of sines and cosines of linear arguments, linearised:
//!   `sin A sin B = (cos(A − B) − cos(A + B))/2`,
//!   `cos A cos B = (cos(A − B) + cos(A + B))/2`,
//!   `sin A cos B = (sin(A + B) + sin(A − B))/2`, so `cos² x` is
//!   `(1 + cos 2x)/2`; the sum goes term by term.

use cas_num::algebra::Rationals;
use cas_num::{Integer, Rational};
use cas_poly::{Dense, DenseRing};

use super::{call, dense_expr, integral_at, polynomial, slope};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// How deep the heuristics nest.
pub(super) const MAX_DEPTH: u8 = 3;

/// The parts of an integrand that derivative-divides tries, at most.
const MAX_CANDIDATES: usize = 12;

/// The sines and cosines a product may have to be linearised, at most,
/// powers counted: `2^6` terms before like ones add up.
const MAX_TRIGONOMETRIC: usize = 6;

type Found = Result<Option<Expr>, EvalError>;

/// An antiderivative of `f` by the module's heuristics, at `depth`.
pub(super) fn integral(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr, depth: u8) -> Found {
    if depth >= MAX_DEPTH {
        return Ok(None);
    }
    let factors: Vec<Expr> = if f.head() == Some(&Expr::from(Symbol::TIMES)) {
        f.args().to_vec()
    } else {
        vec![f.clone()]
    };
    if let Some(g) = exponential_trigonometric(ev, &factors, x)? {
        return Ok(Some(g));
    }
    if let Some(g) = tabular(ev, e, &factors, x)? {
        return Ok(Some(g));
    }
    if let Some(g) = substitution(ev, e, f, x, depth)? {
        return Ok(Some(g));
    }
    if let Some(g) = inverse_parts(ev, e, &factors, x, depth)? {
        return Ok(Some(g));
    }
    linearised(ev, e, &factors, x, depth)
}

/// What is left of `factors` without the one at `i`: their product.
fn others(factors: &[Expr], i: usize) -> Expr {
    factors
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, f)| f.clone())
        .product()
}

/// The rate `α` of an exponential `t`, `t′ = α t`: `Exp[u]` and `E^u` for
/// a linear `u = a x + b`, `α = a`, and `c^u` for a positive rational `c`,
/// `α = a ln c`.
fn exponential(t: &Expr, x: &Expr) -> Option<Expr> {
    match (t.head()?.as_symbol()?, t.args()) {
        (Symbol::EXP, [u]) => slope(u, x),
        (Symbol::POWER, [c, u]) if !c.contains(x) => {
            let a = slope(u, x)?;
            if c.as_symbol() == Some(Symbol::E) {
                return Some(a);
            }
            c.to_rational()
                .filter(Rational::is_positive)
                .map(|_| a * call(Symbol::LOG, c.clone()))
        }
        _ => None,
    }
}

/// A sine or cosine of a linear argument: whether a sine, the argument,
/// its slope.
fn trigonometric(t: &Expr, x: &Expr) -> Option<(bool, Expr, Expr)> {
    match (t.head()?.as_symbol()?, t.args()) {
        (Symbol::SIN, [u]) => Some((true, u.clone(), slope(u, x)?)),
        (Symbol::COS, [u]) => Some((false, u.clone(), slope(u, x)?)),
        _ => None,
    }
}

fn power(e: &Expr, k: usize) -> Expr {
    e.pow(i64::try_from(k).expect("a degree in memory"))
}

/// The exponent `w` of an exponential `t = e^w`: `u` of `Exp[u]` and
/// `E^u`, `u ln c` of `c^u`.
fn exponent(t: &Expr) -> Option<Expr> {
    match (t.head()?.as_symbol()?, t.args()) {
        (Symbol::EXP, [u]) => Some(u.clone()),
        (Symbol::POWER, [c, u]) if c.as_symbol() == Some(Symbol::E) => Some(u.clone()),
        (Symbol::POWER, [c, u]) => Some(u.clone() * call(Symbol::LOG, c.clone())),
        _ => None,
    }
}

/// `∫ e^u sin v` and `∫ e^u cos v` (the module's first rule).
fn exponential_trigonometric(ev: &mut Eval<'_>, factors: &[Expr], x: &Expr) -> Found {
    let [s, t] = factors else {
        return Ok(None);
    };
    let (exp, trig) = if exponential(s, x).is_some() {
        (s, t)
    } else {
        (t, s)
    };
    let (Some(alpha), Some((sine, v, gamma))) = (exponential(exp, x), trigonometric(trig, x))
    else {
        return Ok(None);
    };
    let denominator = ev.eval(&(alpha.pow(2) + gamma.pow(2)))?;
    if denominator == Expr::from(0) {
        let Some(w) = exponent(exp) else {
            return Ok(None);
        };
        return resonant(ev, &w, sine, &v, x).map(Some);
    }
    let (sin, cos) = (call(Symbol::SIN, v.clone()), call(Symbol::COS, v));
    let numerator = if sine {
        alpha * sin - gamma * cos
    } else {
        alpha * cos + gamma * sin
    };
    Ok(Some(exp.clone() * numerator / denominator))
}

/// `∫ e^w sin v` or `∫ e^w cos v` at resonance (the module's first rule):
/// `(e^(w+iv) ∓ e^(w−iv))/(2i)` or `/2`, each exponential `e^z` by its
/// slope `z′`, `x e^z` for `z′ = 0`.
fn resonant(
    ev: &mut Eval<'_>,
    w: &Expr,
    sine: bool,
    v: &Expr,
    x: &Expr,
) -> Result<Expr, EvalError> {
    let i = Expr::from(Symbol::I);
    let mut terms = Vec::with_capacity(2);
    for sign in [1, -1] {
        let z = ev.eval(&(w.clone() + Expr::from(sign) * i.clone() * v.clone()))?;
        let rate = ev.eval(&Expr::call(Symbol::D, [z.clone(), x.clone()]))?;
        let e = call(Symbol::EXP, z);
        let integral = if rate == Expr::from(0) {
            x.clone() * e
        } else {
            e / rate
        };
        // cos v = (e^(iv) + e^(−iv))/2, sin v = (e^(iv) − e^(−iv))/(2i).
        let weight = if sine {
            Expr::from(sign) / (Expr::from(2) * i.clone())
        } else {
            Expr::from(1) / Expr::from(2)
        };
        terms.push(weight * integral);
    }
    Ok(terms.into_iter().sum())
}

/// `∫ P t` for a polynomial `P` in `x` over ℚ of degree 1 at least and an
/// exponential, sine or cosine `t` of a linear argument: with `t′ = α t`,
/// `t Σ (−1)^k P^(k)/α^(k+1)`; with `v′ = a`, `∫ P sin v` is `−cos v`
/// times `Σ_{k even} (−1)^(k/2) P^(k)/a^(k+1)` plus `sin v` times
/// `Σ_{k odd} (−1)^((k−1)/2) P^(k)/a^(k+1)`, and `∫ P cos v` the same with
/// `sin v` and `cos v` for `−cos v` and `sin v`.
fn tabular(ev: &mut Eval<'_>, e: &Expr, factors: &[Expr], x: &Expr) -> Found {
    for (i, t) in factors.iter().enumerate() {
        let rate = exponential(t, x);
        let trig = trigonometric(t, x);
        if rate.is_none() && trig.is_none() {
            continue;
        }
        let Some(p) = polynomial(ev, e, &others(factors, i), x)? else {
            continue;
        };
        if p.degree().unwrap_or(0) == 0 {
            continue;
        }
        let qx = DenseRing::new(Rationals);
        let mut derivatives: Vec<Dense<Rational>> = vec![p];
        loop {
            let next = qx.derivative(derivatives.last().expect("one"));
            if next.is_zero() {
                break;
            }
            derivatives.push(next);
        }
        // (−1)^(k/2) for an even k, (−1)^((k−1)/2) for an odd one.
        let signed = |k: usize| if (k / 2) % 2 == 1 { -1 } else { 1 };
        if let Some(alpha) = rate {
            let sum: Expr = derivatives
                .iter()
                .enumerate()
                .map(|(k, d)| {
                    let sign = if k % 2 == 0 { 1 } else { -1 };
                    Expr::from(sign) * dense_expr(d, x) / power(&alpha, k + 1)
                })
                .sum();
            return Ok(Some(t.clone() * sum));
        }
        let (sine, v, a) = trig.expect("one of the two");
        let (mut even, mut odd) = (Vec::new(), Vec::new());
        for (k, d) in derivatives.iter().enumerate() {
            let term = Expr::from(signed(k)) * dense_expr(d, x) / power(&a, k + 1);
            if k % 2 == 0 {
                even.push(term);
            } else {
                odd.push(term);
            }
        }
        let (even, odd): (Expr, Expr) = (even.into_iter().sum(), odd.into_iter().sum());
        let (sin, cos) = (call(Symbol::SIN, v.clone()), call(Symbol::COS, v));
        return Ok(Some(if sine {
            -cos * even + sin * odd
        } else {
            sin * even + cos * odd
        }));
    }
    Ok(None)
}

/// The fresh variable of a substitution at `depth`, one a depth so that a
/// nested substitution has its own.
fn fresh(depth: u8) -> Expr {
    Expr::from(Symbol::new(&format!("Integrate${depth}")))
}

/// Derivative-divides (the module's third rule): the parts `g` of `f`
/// tried are applications that are not `f`, contain `x` and are not linear
/// in it, the first [`MAX_CANDIDATES`].
fn substitution(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr, depth: u8) -> Found {
    let u = fresh(depth);
    if f.contains(&u) {
        return Ok(None);
    }
    let candidates: Vec<Expr> = f
        .distinct()
        .filter(|g| g != &f && g.head().is_some() && g.contains(x) && slope(g, x).is_none())
        .take(MAX_CANDIDATES)
        .cloned()
        .collect();
    for g in candidates {
        let derivative = ev.eval(&Expr::call(Symbol::D, [g.clone(), x.clone()]))?;
        if derivative == Expr::from(0) {
            continue;
        }
        let quotient = ev.eval(&Expr::call(Symbol::TOGETHER, [f.clone() / derivative]))?;
        let in_u = quotient.substitute(&[(g.clone(), u.clone())]);
        if in_u.contains(x) {
            continue;
        }
        let in_u = ev.eval(&in_u)?;
        if let Some(antiderivative) = integral_at(ev, e, &in_u, &u, depth + 1)?
            && !antiderivative.contains(&Expr::from(Symbol::INTEGRATE))
        {
            return Ok(Some(antiderivative.substitute(&[(u, g)])));
        }
    }
    Ok(None)
}

/// `L′` for `L` a logarithm, a power of one by an integer at least 2, or
/// an inverse tangent, sine or cosine of a linear argument `u = a x + b`.
fn inverse_derivative(l: &Expr, x: &Expr) -> Option<Expr> {
    let root = |u: &Expr| {
        let half = Expr::from(Rational::new(Integer::from(-1), Integer::from(2)));
        (Expr::from(1) - u.pow(2)).pow(half)
    };
    match (l.head()?.as_symbol()?, l.args()) {
        (Symbol::LOG, [u]) => Some(slope(u, x)? / u.clone()),
        (Symbol::POWER, [g, n]) if g.head() == Some(&Expr::from(Symbol::LOG)) => {
            let n = n.as_i64().filter(|&n| n >= 2)?;
            let [u] = g.args() else {
                return None;
            };
            Some(Expr::from(n) * g.pow(n - 1) * slope(u, x)? / u.clone())
        }
        (Symbol::ARC_TAN, [u]) => Some(slope(u, x)? / (Expr::from(1) + u.pow(2))),
        (Symbol::ARC_SIN, [u]) => Some(slope(u, x)? * root(u)),
        (Symbol::ARC_COS, [u]) => Some(-slope(u, x)? * root(u)),
        _ => None,
    }
}

/// `∫ P L = Q L − ∫ Q L′` (the module's fourth rule), when the second
/// integral is found whole.
fn inverse_parts(ev: &mut Eval<'_>, e: &Expr, factors: &[Expr], x: &Expr, depth: u8) -> Found {
    for (i, l) in factors.iter().enumerate() {
        let Some(derivative) = inverse_derivative(l, x) else {
            continue;
        };
        let others = others(factors, i);
        let q = if others.contains(x) {
            let Some(p) = polynomial(ev, e, &others, x)? else {
                continue;
            };
            let qx = DenseRing::new(Rationals);
            let mut coeffs = vec![Rational::ZERO];
            coeffs.extend(
                p.coeffs()
                    .iter()
                    .enumerate()
                    .map(|(k, c)| c / &Rational::from(Integer::from(k + 1))),
            );
            dense_expr(&qx.poly(coeffs), x)
        } else {
            others * x.clone()
        };
        let rest = ev.eval(&(q.clone() * derivative))?;
        let Some(r) = integral_at(ev, e, &rest, x, depth + 1)? else {
            continue;
        };
        if r.contains(&Expr::from(Symbol::INTEGRATE)) {
            continue;
        }
        return Ok(Some(q * l.clone() - r));
    }
    Ok(None)
}

/// A term of a linearised product: a coefficient, and a sine (`Some(true)`)
/// or cosine (`Some(false)`) of an argument, or 1 (`None`).
type Term = (Rational, Option<bool>, Expr);

/// Products of sines and cosines (the module's last rule), times what
/// else `factors` has, when the sum integrates whole.
fn linearised(ev: &mut Eval<'_>, e: &Expr, factors: &[Expr], x: &Expr, depth: u8) -> Found {
    let mut trig = Vec::new();
    let mut rest = Vec::new();
    for f in factors {
        let (base, n) = match (f.head().and_then(Expr::as_symbol), f.args()) {
            (Some(Symbol::POWER), [b, n]) => {
                match n.as_i64().and_then(|n| usize::try_from(n).ok()) {
                    Some(n) if (1..=MAX_TRIGONOMETRIC).contains(&n) => (b, n),
                    _ => (f, 1),
                }
            }
            _ => (f, 1),
        };
        match trigonometric(base, x) {
            Some((sine, v, _)) => trig.extend(std::iter::repeat_n((sine, v), n)),
            None => rest.push(f.clone()),
        }
    }
    if trig.len() < 2 || trig.len() > MAX_TRIGONOMETRIC {
        return Ok(None);
    }
    let mut terms: Vec<Term> = vec![(Rational::ONE, None, Expr::from(0))];
    for (sine, b) in trig {
        let mut next = Vec::new();
        for (c, kind, a) in terms {
            let Some(s) = kind else {
                next.push((c, Some(sine), b.clone()));
                continue;
            };
            let sum = ev.eval(&(a.clone() + b.clone()))?;
            let difference = ev.eval(&(a - b.clone()))?;
            let half = &c / &Rational::from(Integer::from(2));
            let minus = -&half;
            next.extend(match (s, sine) {
                (true, true) => [(half, Some(false), difference), (minus, Some(false), sum)],
                (false, false) => [
                    (half.clone(), Some(false), difference),
                    (half, Some(false), sum),
                ],
                (true, false) => [
                    (half.clone(), Some(true), sum),
                    (half, Some(true), difference),
                ],
                (false, true) => [(half, Some(true), sum), (minus, Some(true), difference)],
            });
        }
        terms = next;
    }
    let zero = Expr::from(0);
    let rest: Expr = rest.into_iter().product();
    let sum: Expr = terms
        .into_iter()
        .filter_map(|(c, kind, a)| {
            let t = match kind {
                Some(true) if a == zero => return None,
                Some(false) if a == zero => Expr::from(1),
                None => Expr::from(1),
                Some(true) => call(Symbol::SIN, a),
                Some(false) => call(Symbol::COS, a),
            };
            Some(Expr::from(c) * rest.clone() * t)
        })
        .sum();
    let sum = ev.eval(&sum)?;
    let found = integral_at(ev, e, &sum, x, depth + 1)?;
    Ok(found.filter(|g| !g.contains(&Expr::from(Symbol::INTEGRATE))))
}
