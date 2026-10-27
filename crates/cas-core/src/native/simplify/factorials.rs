//! Factorials and binomials in the search (CAS-DOC-14 §9, CAS-235).
//!
//! `Binomial[n, k]` as `n!/(k! (n − k)!)`; then the factorials whose
//! arguments differ by an integer, at most [`MAX_OFFSET`], are written
//! over one factorial of their class: the smallest argument `u` shifted by
//! `s` of −1, 0 and 1, each a candidate, `(u + j)! = (u + s)! (u + s + 1) ⋯
//! (u + j)` above it and `(u + s)!/((u + j + 1) ⋯ (u + s))` below; each as
//! it evaluates and made `Together`. So `n!/(n-1)!` is `n`,
//! `binomial(n,k)/binomial(n,k-1)` is `(n-k+1)/k`, and `(n+2)/(n+2)!` is
//! `1/(n+1)!` (the shift −1).

use super::evaluated;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The largest integer difference of two arguments in one class.
const MAX_OFFSET: i64 = 64;

/// A class: its first argument, and each argument with its offset from it.
type Class = (Expr, Vec<(Expr, i64)>);

/// The factorial candidates of `e`, evaluated.
pub(super) fn candidates(ev: &mut Eval<'_>, e: &Expr) -> Result<Vec<Expr>, EvalError> {
    if !(e.may_contain(Symbol::FACTORIAL) || e.may_contain(Symbol::BINOMIAL)) {
        return Ok(Vec::new());
    }
    let factorial = |u: Expr| Expr::call(Symbol::FACTORIAL, [u]);
    let as_factorials =
        e.map_bottom_up(
            |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
                (Some(Symbol::BINOMIAL), [n, k]) => {
                    factorial(n.clone()) / (factorial(k.clone()) * factorial(n.clone() - k.clone()))
                }
                _ => part.clone(),
            },
        );
    let Some(f) = evaluated(ev, &as_factorials)? else {
        return Ok(Vec::new());
    };
    let classes = classes(ev, &f)?;
    let mut calls = Vec::new();
    if f != *e {
        calls.push(f.clone());
    }
    if !classes.is_empty() {
        for shift in [-1, 0, 1] {
            let rewritten = f.substitute(&rules(&classes, shift));
            calls.push(Expr::call(Symbol::TOGETHER, [rewritten.clone()]));
            calls.push(rewritten);
        }
    }
    let mut out = Vec::new();
    for c in calls {
        if let Some(x) = evaluated(ev, &c)? {
            out.push(x);
        }
    }
    Ok(out)
}

/// The arguments of the factorials of `f` but numbers, in classes.
fn classes(ev: &mut Eval<'_>, f: &Expr) -> Result<Vec<Class>, EvalError> {
    let mut arguments: Vec<Expr> = Vec::new();
    let mut stack = vec![f];
    while let Some(x) = stack.pop() {
        if let View::Apply { head, args } = x.view() {
            if head.as_symbol() == Some(Symbol::FACTORIAL)
                && let [u] = args
                && u.to_rational().is_none()
                && !arguments.contains(u)
            {
                arguments.push(u.clone());
            }
            stack.extend(args);
        }
    }
    let mut classes: Vec<Class> = Vec::new();
    'arguments: for u in arguments {
        for (base, members) in &mut classes {
            let difference = Expr::call(Symbol::EXPAND, [u.clone() - base.clone()]);
            if let Some(d) = evaluated(ev, &difference)?
                .and_then(|d| d.as_i64())
                .filter(|d| d.abs() <= MAX_OFFSET)
            {
                members.push((u, d));
                continue 'arguments;
            }
        }
        classes.push((u.clone(), vec![(u, 0)]));
    }
    Ok(classes)
}

/// Each factorial of a class as one of the class's smallest argument
/// shifted by `shift`, times or over the factors between.
fn rules(classes: &[Class], shift: i64) -> Vec<(Expr, Expr)> {
    let mut rules = Vec::new();
    for (base, members) in classes {
        let low = members.iter().map(|&(_, d)| d).min().unwrap_or(0) + shift;
        let at = |j: i64| base.clone() + Expr::from(j);
        let target = Expr::call(Symbol::FACTORIAL, [at(low)]);
        for (u, d) in members {
            let between = |range: std::ops::RangeInclusive<i64>| {
                Expr::call(Symbol::TIMES, range.map(at).collect::<Vec<_>>())
            };
            let value = if *d >= low {
                target.clone() * between(low + 1..=*d)
            } else {
                target.clone() / between(*d + 1..=low)
            };
            rules.push((Expr::call(Symbol::FACTORIAL, [u.clone()]), value));
        }
    }
    rules
}
