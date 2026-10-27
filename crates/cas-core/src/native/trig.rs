//! The trigonometric rewritings (CAS-DOC-14 §7, CAS-233): `TrigExpand`,
//! giac's `texpand` and `trigexpand`; `TrigReduce`, giac's `tlin`;
//! `TrigToExp`, giac's `trig2exp`; `ExpToTrig`, giac's `exp2trig`. Each
//! rewrites every part of its argument, lists and equations among them,
//! leaves what it does not rewrite, and takes one step a call.
//!
//! - `TrigExpand`: `Sin`, `Cos` and `Tan` of a sum of at most
//!   [`MAX_SUMMANDS`] terms by the addition formulas, the first term against
//!   the rest, and of an integer multiple `n v`, `2 ≤ n ≤` [`MAX_MULTIPLE`],
//!   `v` not a sum, by Chebyshev's polynomials: `cos(n v) = T_n(cos v)`,
//!   `sin(n v) = sin v U_(n−1)(cos v)`, and `tan(n v)` the ratio of the
//!   odd and the even parts of `(1 + i tan v)^n` (giac's
//!   `texpand(cos(3*x))` is `4*cos(x)^3-3*cos(x)`).
//! - `TrigReduce`: a product of `Sin`, `Cos` and their positive integer
//!   powers, of total degree 2 to [`MAX_DEGREE`], as a sum by `sin a sin b
//!   = (cos(a − b) − cos(a + b))/2`, `cos a cos b = (cos(a − b) + cos(a +
//!   b))/2` and `sin a cos b = (sin(a + b) + sin(a − b))/2`, the angles
//!   expanded, a negative one's sign out (a sum's when its first term is
//!   negative), and the terms of one function at one angle together; the
//!   other factors stay as they are (giac's
//!   `tlin(cos(x)^3)` is `3/4*cos(x)+1/4*cos(3*x)`).
//! - `TrigToExp`: `sin u = (e^(iu) − e^(−iu))/(2i)`, `cos u = (e^(iu) +
//!   e^(−iu))/2` and `tan u = −i (e^(iu) − e^(−iu))/(e^(iu) + e^(−iu))`;
//!   `sinh`, `cosh` and `tanh` the same with `e^(±u)` (CAS-DOC-15 §2).
//! - `ExpToTrig`: `Exp[a + i b]` and `E^(a + i b)`, the terms of the
//!   exponent with a non-real coefficient its `i b`, as `e^a (cos b + i
//!   sin b)` when every term of `a` and `b` is real by the domain or the
//!   assumptions.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of trigonometry: sin a, cos b, angles u and v"
)]

use cas_num::{Integer, Rational};

use super::special::negated;
use crate::assume::{Set, Truth, element};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The terms of a sum expanded at most.
const MAX_SUMMANDS: usize = 8;

/// The largest multiple expanded.
const MAX_MULTIPLE: i64 = 64;

/// The largest degree of a product reduced.
const MAX_DEGREE: usize = 64;

/// `TrigExpand[x]`.
pub(super) fn trig_expand(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = expand(x);
    ev.step("TrigExpand", e, &result);
    Ok(Some(result))
}

/// `TrigReduce[x]`.
pub(super) fn trig_reduce(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = reduce(ev, x)?;
    ev.step("TrigReduce", e, &result);
    Ok(Some(result))
}

/// `TrigToExp[x]`.
pub(super) fn trig_to_exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = to_exp(x);
    ev.step("TrigToExp", e, &result);
    Ok(Some(result))
}

/// `ExpToTrig[x]`.
pub(super) fn exp_to_trig(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let result = to_trig(ev.context(), x);
    ev.step("ExpToTrig", e, &result);
    Ok(Some(result))
}

fn sin(u: Expr) -> Expr {
    Expr::call(Symbol::SIN, [u])
}

fn cos(u: Expr) -> Expr {
    Expr::call(Symbol::COS, [u])
}

fn tan(u: Expr) -> Expr {
    Expr::call(Symbol::TAN, [u])
}

fn i() -> Expr {
    Expr::from(Symbol::I)
}

/// `x` with `Sin`, `Cos` and `Tan` of sums and multiples expanded, built.
pub(super) fn expand(x: &Expr) -> Expr {
    x.map_bottom_up(
        |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(Symbol::SIN), [u]) if expandable(u) => sin_cos(u).0,
            (Some(Symbol::COS), [u]) if expandable(u) => sin_cos(u).1,
            (Some(Symbol::TAN), [u]) if expandable(u) => tangent(u),
            _ => part.clone(),
        },
    )
}

/// Whether `u` is a sum or a multiple that [`sin_cos`] expands.
fn expandable(u: &Expr) -> bool {
    split(u).is_some() || multiple(u).is_some()
}

/// `u` as its first term and the rest, for a sum of at most
/// [`MAX_SUMMANDS`] terms.
fn split(u: &Expr) -> Option<(Expr, Expr)> {
    let (Some(Symbol::PLUS), [a, rest @ ..]) = (u.head().and_then(Expr::as_symbol), u.args())
    else {
        return None;
    };
    if rest.is_empty() || rest.len() >= MAX_SUMMANDS {
        return None;
    }
    let b = match rest {
        [b] => b.clone(),
        _ => Expr::call(Symbol::PLUS, rest.to_vec()),
    };
    Some((a.clone(), b))
}

/// `u` as `n v`, `2 ≤ n ≤` [`MAX_MULTIPLE`], `v` not a sum.
fn multiple(u: &Expr) -> Option<(i64, Expr)> {
    let (Some(Symbol::TIMES), [n, rest @ ..]) = (u.head().and_then(Expr::as_symbol), u.args())
    else {
        return None;
    };
    let n = n.as_i64().filter(|n| (2..=MAX_MULTIPLE).contains(n))?;
    let v = match rest {
        [v] => v.clone(),
        _ => Expr::call(Symbol::TIMES, rest.to_vec()),
    };
    (v.head().and_then(Expr::as_symbol) != Some(Symbol::PLUS)).then_some((n, v))
}

/// `(sin u, cos u)` expanded, built.
fn sin_cos(u: &Expr) -> (Expr, Expr) {
    if let Some((a, b)) = split(u) {
        let (sa, ca) = sin_cos(&a);
        let (sb, cb) = sin_cos(&b);
        return (
            sa.clone() * cb.clone() + ca.clone() * sb.clone(),
            ca * cb - sa * sb,
        );
    }
    if let Some((n, v)) = multiple(u) {
        let (s, c) = (sin(v.clone()), cos(v));
        let (tn, un) = chebyshev(n);
        return (s * polynomial(&un, &c), polynomial(&tn, &c));
    }
    (sin(u.clone()), cos(u.clone()))
}

/// `tan u` expanded, built.
fn tangent(u: &Expr) -> Expr {
    if let Some((a, b)) = split(u) {
        let (ta, tb) = (tangent(&a), tangent(&b));
        return (ta.clone() + tb.clone()) / (Expr::from(1) - ta * tb);
    }
    if let Some((n, v)) = multiple(u) {
        // (1 + i t)^n: its imaginary part over its real part.
        let t = tan(v);
        let (mut odd, mut even) = (Vec::new(), Vec::new());
        let mut binomial = Integer::ONE;
        for k in 0..=n {
            if k > 0 {
                binomial = (&binomial * &Integer::from(n - k + 1)).div_exact(&Integer::from(k));
            }
            let sign = if (k / 2) % 2 == 0 { 1 } else { -1 };
            let term = Expr::from(&binomial * &Integer::from(sign)) * t.pow(k);
            if k % 2 == 1 {
                odd.push(term);
            } else {
                even.push(term);
            }
        }
        return Expr::call(Symbol::PLUS, odd) / Expr::call(Symbol::PLUS, even);
    }
    tan(u.clone())
}

/// The coefficients of Chebyshev's `T_n` and `U_(n−1)`, from the constant.
fn chebyshev(n: i64) -> (Vec<Integer>, Vec<Integer>) {
    let next = |p: &[Integer], q: &[Integer]| -> Vec<Integer> {
        // 2 x p − q.
        let mut r = vec![Integer::ZERO; p.len() + 1];
        for (k, c) in p.iter().enumerate() {
            r[k + 1] = &r[k + 1] + &(c * &Integer::from(2));
        }
        for (k, c) in q.iter().enumerate() {
            r[k] = &r[k] - c;
        }
        r
    };
    let (mut t0, mut t1) = (vec![Integer::ONE], vec![Integer::ZERO, Integer::ONE]);
    let (mut u0, mut u1) = (vec![Integer::ONE], vec![Integer::ZERO, Integer::from(2)]);
    for _ in 1..n {
        let t2 = next(&t1, &t0);
        let u2 = next(&u1, &u0);
        (t0, t1) = (t1, t2);
        (u0, u1) = (u1, u2);
    }
    // T_n, and U_(n−1): one step behind.
    (t1, u0)
}

/// `Σ cₖ x^k`, built.
fn polynomial(coefficients: &[Integer], x: &Expr) -> Expr {
    let terms: Vec<Expr> = coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(k, c)| Expr::from(c.clone()) * x.pow(i64::try_from(k).expect("a small degree")))
        .collect();
    Expr::call(Symbol::PLUS, terms)
}

/// A function of a reduced term.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    One,
    Sin,
    Cos,
}

/// `c f(θ)`.
type Term = (Rational, Kind, Expr);

/// `x` with its products of `Sin` and `Cos` as sums, evaluated.
pub(super) fn reduce(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let mut failure = None;
    let r = x.map_bottom_up(|part| {
        if failure.is_some() {
            return part.clone();
        }
        match linearized(ev, part) {
            Ok(Some(s)) => s,
            Ok(None) => part.clone(),
            Err(e) => {
                failure = Some(e);
                part.clone()
            }
        }
    });
    if let Some(e) = failure {
        return Err(e);
    }
    ev.eval(&r)
}

/// `Sin[u]`, `Cos[u]` or one to a positive integer power: its function,
/// angle and power.
fn trig_power(f: &Expr) -> Option<(Kind, Expr, usize)> {
    let kind = |g: &Expr| match (g.head().and_then(Expr::as_symbol), g.args()) {
        (Some(Symbol::SIN), [u]) => Some((Kind::Sin, u.clone())),
        (Some(Symbol::COS), [u]) => Some((Kind::Cos, u.clone())),
        _ => None,
    };
    if let Some((k, u)) = kind(f) {
        return Some((k, u, 1));
    }
    let (Some(Symbol::POWER), [g, n]) = (f.head().and_then(Expr::as_symbol), f.args()) else {
        return None;
    };
    let n = n
        .as_i64()
        .and_then(|n| usize::try_from(n).ok())
        .filter(|n| (1..=MAX_DEGREE).contains(n))?;
    kind(g).map(|(k, u)| (k, u, n))
}

/// A product, or a power, with its factors `Sin` and `Cos` as a sum; `None`
/// for one of degree below 2 or past [`MAX_DEGREE`].
fn linearized(ev: &mut Eval<'_>, part: &Expr) -> Result<Option<Expr>, EvalError> {
    let factors: &[Expr] = match part.head().and_then(Expr::as_symbol) {
        Some(Symbol::TIMES) => part.args(),
        Some(Symbol::POWER) => std::slice::from_ref(part),
        _ => return Ok(None),
    };
    let mut trig: Vec<(Kind, Expr)> = Vec::new();
    let mut rest = Vec::new();
    for f in factors {
        match trig_power(f) {
            Some((k, u, n)) => trig.extend(std::iter::repeat_n((k, u), n)),
            None => rest.push(f.clone()),
        }
    }
    if trig.len() < 2 || trig.len() > MAX_DEGREE {
        return Ok(None);
    }
    let mut terms: Vec<Term> = vec![(Rational::ONE, Kind::One, Expr::from(0))];
    for (g, phi) in trig {
        terms = times(ev, &terms, g, &phi)?;
    }
    let sum: Vec<Expr> = terms
        .into_iter()
        .map(|(c, k, theta)| {
            Expr::from(c)
                * match k {
                    Kind::One => Expr::from(1),
                    Kind::Sin => sin(theta),
                    Kind::Cos => cos(theta),
                }
        })
        .collect();
    rest.push(Expr::call(Symbol::PLUS, sum));
    Ok(Some(Expr::call(Symbol::TIMES, rest)))
}

/// `Σ c f(θ)` times `g(φ)`, by the product-to-sum formulas.
fn times(ev: &mut Eval<'_>, terms: &[Term], g: Kind, phi: &Expr) -> Result<Vec<Term>, EvalError> {
    let half = Rational::new(Integer::ONE, Integer::from(2));
    let mut out: Vec<Term> = Vec::new();
    for (c, f, theta) in terms {
        let sum = || Expr::call(Symbol::EXPAND, [theta.clone() + phi.clone()]);
        let difference = || Expr::call(Symbol::EXPAND, [theta.clone() - phi.clone()]);
        let h = &half * c;
        let new: Vec<(Rational, Kind, Expr)> = match (f, g) {
            (Kind::One, _) => vec![(c.clone(), g, phi.clone())],
            (Kind::Sin, Kind::Sin) => {
                vec![(h.clone(), Kind::Cos, difference()), (-h, Kind::Cos, sum())]
            }
            (Kind::Cos, Kind::Cos) => {
                vec![(h.clone(), Kind::Cos, difference()), (h, Kind::Cos, sum())]
            }
            (Kind::Sin, Kind::Cos) => {
                vec![(h.clone(), Kind::Sin, sum()), (h, Kind::Sin, difference())]
            }
            (Kind::Cos, Kind::Sin) => {
                vec![(h.clone(), Kind::Sin, sum()), (-h, Kind::Sin, difference())]
            }
            (_, Kind::One) => vec![(c.clone(), *f, theta.clone())],
        };
        for (c, k, angle) in new {
            let angle = ev.eval(&angle)?;
            let term = normalized(ev, c, k, angle)?;
            add(&mut out, term);
        }
    }
    out.retain(|(c, _, _)| !c.is_zero());
    Ok(out)
}

/// `c f(θ)` with `sin 0 = 0`, `cos 0 = 1`, and the sign of a negative angle
/// out: a negative number, a product with a negative coefficient, or a sum
/// whose first term is one, so that `sin(y − x)` and `sin(x − y)` meet.
fn normalized(ev: &mut Eval<'_>, c: Rational, k: Kind, angle: Expr) -> Result<Term, EvalError> {
    if angle == Expr::from(0) {
        return Ok(match k {
            Kind::Sin => (Rational::ZERO, Kind::One, angle),
            _ => (c, Kind::One, angle),
        });
    }
    let opposite = match negated(&angle) {
        Some(positive) => Some(positive),
        None => match (angle.head().and_then(Expr::as_symbol), angle.args()) {
            (Some(Symbol::PLUS), [first, ..]) if negated(first).is_some() => {
                Some(ev.eval(&Expr::call(Symbol::EXPAND, [-angle.clone()]))?)
            }
            _ => None,
        },
    };
    Ok(match opposite {
        Some(positive) if k == Kind::Sin => (-c, k, positive),
        Some(positive) => (c, k, positive),
        None => (c, k, angle),
    })
}

/// `terms` plus `c f(θ)`, with the term of the same function and angle.
fn add(terms: &mut Vec<Term>, (c, k, angle): Term) {
    match terms
        .iter_mut()
        .find(|(_, f, theta)| *f == k && (k == Kind::One || *theta == angle))
    {
        Some((d, _, _)) => *d = &*d + &c,
        None => terms.push((c, k, angle)),
    }
}

/// `x` with `Sin`, `Cos`, `Tan` and the hyperbolic `Sinh`, `Cosh`, `Tanh`
/// as exponentials, built.
pub(super) fn to_exp(x: &Expr) -> Expr {
    x.map_bottom_up(|part| {
        let (Some(f), [u]) = (part.head().and_then(Expr::as_symbol), part.args()) else {
            return part.clone();
        };
        // e^(±iu) for the circular functions, e^(±u) for the hyperbolic.
        let circular = matches!(f, Symbol::SIN | Symbol::COS | Symbol::TAN);
        let hyperbolic = matches!(f, Symbol::SINH | Symbol::COSH | Symbol::TANH);
        if !circular && !hyperbolic {
            return part.clone();
        }
        let v = if circular { i() * u.clone() } else { u.clone() };
        let plus = Expr::call(Symbol::EXP, [v.clone()]);
        let minus = Expr::call(Symbol::EXP, [-v]);
        match f {
            Symbol::SIN => (plus - minus) / (Expr::from(2) * i()),
            Symbol::COS | Symbol::COSH => (plus + minus) / Expr::from(2),
            Symbol::TAN => -i() * (plus.clone() - minus.clone()) / (plus + minus),
            Symbol::SINH => (plus - minus) / Expr::from(2),
            _ => (plus.clone() - minus.clone()) / (plus + minus),
        }
    })
}

/// `x` with its exponentials of complex exponents as `Cos` and `Sin`,
/// built.
pub(super) fn to_trig(ctx: &Context, x: &Expr) -> Expr {
    x.map_bottom_up(
        |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
            (Some(Symbol::EXP), [u]) => trig_of(ctx, u).unwrap_or_else(|| part.clone()),
            (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => {
                trig_of(ctx, u).unwrap_or_else(|| part.clone())
            }
            _ => part.clone(),
        },
    )
}

/// `e^u` as `e^a (cos b + i sin b)`, `u = a + i b`; `None` when `u` has no
/// non-real coefficient or a part not known real.
fn trig_of(ctx: &Context, u: &Expr) -> Option<Expr> {
    let terms: &[Expr] = match u.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => u.args(),
        _ => std::slice::from_ref(u),
    };
    let real = |x: &Expr| element(ctx, x, Set::Reals) == Truth::True;
    let (mut re, mut im) = (Vec::new(), Vec::new());
    for t in terms {
        // A complex coefficient and the rest, or `t` itself.
        let (c, rest): (&Expr, Vec<Expr>) = match (t.head().and_then(Expr::as_symbol), t.args()) {
            (Some(Symbol::TIMES), [c, rest @ ..]) if matches!(c.view(), View::Complex(_)) => {
                (c, rest.to_vec())
            }
            _ => (t, Vec::new()),
        };
        let rest = match rest.len() {
            0 => Expr::from(1),
            1 => rest[0].clone(),
            _ => Expr::call(Symbol::TIMES, rest),
        };
        match c.view() {
            View::Complex(z) if real(&rest) => {
                re.push(Expr::from(z.re().clone()) * rest.clone());
                im.push(Expr::from(z.im().clone()) * rest);
            }
            View::Complex(_) => return None,
            _ if real(t) => re.push(t.clone()),
            _ => return None,
        }
    }
    if im.is_empty() {
        return None;
    }
    let (a, b) = (Expr::call(Symbol::PLUS, re), Expr::call(Symbol::PLUS, im));
    Some(Expr::call(Symbol::EXP, [a]) * (cos(b.clone()) + i() * sin(b)))
}

#[cfg(test)]
mod tests;
