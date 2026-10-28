//! `TrigCollect` (CAS-DOC-36 §5, CAS-399), giac's `tcollect`: the argument
//! linearized by `TrigReduce` (CAS-DOC-14 §7) and expanded, as `TrigReduce`
//! leaves `−(1/2 − cos(2x)/2)` undistributed, then in each sum `a sin u +
//! b cos u` as `R cos(u − φ)`, with `R = √(a² + b²)` and `φ = Arg(b + i
//! a)`, where `a` and `b` are real by the domain or the assumptions and
//! free of `u`'s symbols. `Arg` gives the exact angle where it has one:
//! `sin(x)+cos(x)` is `sqrt(2)*cos(x-pi/4)`, and `3*sin(x)+4*cos(x)` is
//! `5*cos(x-atan(3/4))`. giac writes symbolic coefficients with `sign`;
//! here they stay unless the assumptions make them real.

use cas_num::{Integer, Rational};

use super::head;
use crate::assume::{Set, Truth, element};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::trig::reduce;
use crate::symbol::Symbol;

/// `t = c f[u]`, `f` `Sin` or `Cos` and `c` the other factors, as `(f, u,
/// c)`.
fn trig_term(t: &Expr) -> Option<(Symbol, Expr, Expr)> {
    let simple = |f: &Expr| match (head(f), f.args()) {
        (Some(g @ (Symbol::SIN | Symbol::COS)), [u]) => Some((g, u.clone())),
        _ => None,
    };
    if let Some((g, u)) = simple(t) {
        return Some((g, u, Expr::from(1)));
    }
    if head(t) != Some(Symbol::TIMES) {
        return None;
    }
    let mut found = None;
    let mut others = Vec::new();
    for f in t.args() {
        match simple(f) {
            Some(_) if found.is_some() => return None,
            Some(gu) => found = Some(gu),
            None => others.push(f.clone()),
        }
    }
    let (g, u) = found?;
    Some((g, u, Expr::call(Symbol::TIMES, others)))
}

/// Whether `c` is a coefficient `TrigCollect` takes for the angle `u`: real
/// and free of `u`'s symbols.
fn coefficient(ctx: &Context, c: &Expr, u: &Expr) -> bool {
    u.symbols()
        .into_iter()
        .filter(|s| !s.is_builtin())
        .all(|s| c.free_of(&Expr::from(s)))
        && element(ctx, c, Set::Reals) == Truth::True
}

/// The sum with its pairs `a sin u + b cos u` as `R cos(u − φ)`, built;
/// `None` without a pair.
fn collected(ctx: &Context, sum: &Expr) -> Option<Expr> {
    let ts = sum.args();
    let found: Vec<Option<(Symbol, Expr, Expr)>> = ts.iter().map(trig_term).collect();
    let mut used = vec![false; ts.len()];
    let mut out = Vec::new();
    for (i, sine) in found.iter().enumerate() {
        let Some((Symbol::SIN, u, a)) = sine else {
            continue;
        };
        let cosine = found.iter().enumerate().find_map(|(j, t)| match t {
            Some((Symbol::COS, v, b)) if v == u && !used[j] => Some((j, b)),
            _ => None,
        });
        let Some((j, b)) = cosine else {
            continue;
        };
        if !coefficient(ctx, a, u) || !coefficient(ctx, b, u) {
            continue;
        }
        used[i] = true;
        used[j] = true;
        let half = Rational::new(Integer::ONE, Integer::from(2));
        let amplitude = (a.pow(2) + b.pow(2)).pow(Expr::from(half));
        let phase = Expr::call(Symbol::ARG, [b.clone() + Expr::from(Symbol::I) * a.clone()]);
        out.push(amplitude * Expr::call(Symbol::COS, [u.clone() - phase]));
    }
    if out.is_empty() {
        return None;
    }
    out.extend(
        ts.iter()
            .zip(&used)
            .filter(|&(_, &u)| !u)
            .map(|(t, _)| t.clone()),
    );
    Some(Expr::call(Symbol::PLUS, out))
}

/// `x` linearized, with its sums' pairs of a sine and a cosine of one angle
/// collected, built.
pub(super) fn collect(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let reduced = reduce(ev, x)?;
    let linear = ev.eval(&Expr::call(Symbol::EXPAND, [reduced]))?;
    let ctx = ev.context();
    Ok(linear.map_bottom_up(|part| match head(part) {
        Some(Symbol::PLUS) => collected(ctx, part).unwrap_or_else(|| part.clone()),
        _ => part.clone(),
    }))
}

#[cfg(test)]
mod tests;
