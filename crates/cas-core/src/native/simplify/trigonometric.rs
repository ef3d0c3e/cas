//! Trigonometric simplification (CAS-DOC-14 §8, CAS-234): candidates for an
//! expression with `Sin`, `Cos` or `Tan`, from it and from it with `Tan[u]`
//! as `Sin[u]/Cos[u]` (evaluated), its bases:
//!
//! - `TrigReduce` and `TrigExpand` of each base, as they are and made
//!   `Together`;
//! - the Pythagorean forms of each base, `Sin[u]^n` as `Sin[u]^(n mod 2)
//!   (1 − Cos[u]^2)^(n div 2)` for `|n| ≥ 2`, and the same with `Sin` and
//!   `Cos` the other way, made `Together` and `Factor`;
//! - the second base made `Together`, and `Sin[u]^p Cos[u]^q` with `p` and
//!   `q` of opposite signs with `Tan[u]` out: `sin(x)/cos(x)` is `tan(x)`.
//!
//! The search's rounds compose them: `(1+tan(x)^2)*cos(x)^2` is
//! `cos(x)^2+sin(x)^2` in the first, 1 in the second.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of trigonometry: sin u and cos u to the powers p and q"
)]

use super::evaluated;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The trigonometric candidates of `e`, evaluated.
pub(super) fn candidates(ev: &mut Eval<'_>, e: &Expr) -> Result<Vec<Expr>, EvalError> {
    if ![Symbol::SIN, Symbol::COS, Symbol::TAN]
        .into_iter()
        .any(|s| e.may_contain(s))
    {
        return Ok(Vec::new());
    }
    let mut bases = vec![e.clone()];
    let tan_free = tan_out(e);
    if tan_free != *e
        && let Some(b) = evaluated(ev, &tan_free)?
    {
        bases.push(b);
    }
    let call = |head: Symbol, x: Expr| Expr::call(head, [x]);
    let mut calls = Vec::new();
    for b in &bases {
        for head in [Symbol::TRIG_REDUCE, Symbol::TRIG_EXPAND] {
            let r = call(head, b.clone());
            calls.push(r.clone());
            calls.push(call(Symbol::TOGETHER, r));
        }
        for p in [
            pythagorean(b, Symbol::SIN, Symbol::COS),
            pythagorean(b, Symbol::COS, Symbol::SIN),
        ] {
            if p != *b {
                calls.push(call(Symbol::TOGETHER, p.clone()));
                calls.push(call(Symbol::FACTOR, p));
            }
        }
    }
    if let [_, b] = bases.as_slice() {
        calls.push(call(Symbol::TOGETHER, b.clone()));
    }
    let tan = tan_in(e);
    if tan != *e {
        calls.push(tan);
    }
    let mut out = Vec::new();
    for c in calls {
        if let Some(x) = evaluated(ev, &c)? {
            out.push(x);
        }
    }
    Ok(out)
}

/// `f[u]`'s argument, for `f` one of the two given.
fn argument(e: &Expr, f: Symbol) -> Option<&Expr> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(g), [u]) if g == f => Some(u),
        _ => None,
    }
}

/// `e` with `Tan[u]` as `Sin[u]/Cos[u]`, built.
fn tan_out(e: &Expr) -> Expr {
    e.map_bottom_up(|part| match argument(part, Symbol::TAN) {
        Some(u) => {
            Expr::call(Symbol::SIN, [u.clone()]) * Expr::call(Symbol::COS, [u.clone()]).pow(-1)
        }
        None => part.clone(),
    })
}

/// `e` with `from[u]^n`, `|n| ≥ 2`, as `from[u]^(n mod 2) (1 −
/// to[u]^2)^(n div 2)`, built.
fn pythagorean(e: &Expr, from: Symbol, to: Symbol) -> Expr {
    e.map_bottom_up(|part| {
        let (Some(Symbol::POWER), [b, n]) = (part.head().and_then(Expr::as_symbol), part.args())
        else {
            return part.clone();
        };
        let (Some(u), Some(n)) = (argument(b, from), n.as_i64().filter(|n| n.abs() >= 2)) else {
            return part.clone();
        };
        let square = Expr::from(1) - Expr::call(to, [u.clone()]).pow(2);
        b.pow(n.rem_euclid(2)) * square.pow(n.div_euclid(2))
    })
}

/// `e` with `Sin[u]^p Cos[u]^q`, `p` and `q` of opposite signs, in a
/// product as `Tan[u]^t Sin[u]^(p − t) Cos[u]^(q + t)`, `t = ±min(|p|,
/// |q|)` of the sign of `p`; built.
fn tan_in(e: &Expr) -> Expr {
    e.map_bottom_up(|part| {
        if part.head().and_then(Expr::as_symbol) != Some(Symbol::TIMES) {
            return part.clone();
        }
        // The angles with the powers of their Sin and Cos, and the rest.
        let mut angles: Vec<(Expr, i64, i64)> = Vec::new();
        let mut rest = Vec::new();
        for f in part.args() {
            let (base, n) = match (f.head().and_then(Expr::as_symbol), f.args()) {
                (Some(Symbol::POWER), [b, n]) if n.as_i64().is_some() => {
                    (b, n.as_i64().unwrap_or(1))
                }
                _ => (f, 1),
            };
            let (sin, cos) = (argument(base, Symbol::SIN), argument(base, Symbol::COS));
            let Some(u) = sin.or(cos) else {
                rest.push(f.clone());
                continue;
            };
            let (p, q) = if sin.is_some() { (n, 0) } else { (0, n) };
            match angles.iter_mut().find(|(v, _, _)| v == u) {
                Some((_, sp, cq)) => {
                    *sp += p;
                    *cq += q;
                }
                None => angles.push((u.clone(), p, q)),
            }
        }
        if !angles.iter().any(|&(_, p, q)| p * q < 0) {
            return part.clone();
        }
        for (u, p, q) in angles {
            let t = if p * q < 0 {
                p.signum() * p.abs().min(q.abs())
            } else {
                0
            };
            for (f, n) in [(Symbol::TAN, t), (Symbol::SIN, p - t), (Symbol::COS, q + t)] {
                if n != 0 {
                    rest.push(Expr::call(f, [u.clone()]).pow(n));
                }
            }
        }
        Expr::call(Symbol::TIMES, rest)
    })
}
