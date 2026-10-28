//! Exponentials, powers and logarithms (CAS-DOC-36 §2, CAS-396): the rules
//! of `HyperbolicToExp`, `ExpToPower`, `PowerToExp`, `ExpExpand`,
//! `LogExpand` and `LogCombine`.
//!
//! - `HyperbolicToExp` (`hyp2exp`): `Sinh`, `Cosh` and `Tanh` as `TrigToExp`
//!   writes them; `asinh u = ln(u + √(u² + 1))`, `acosh u = ln(u + √(u + 1)
//!   √(u − 1))` and `atanh u = ln((1 + u)/(1 − u))/2`. giac's `acosh` is
//!   `ln(u + √(u² − 1))`, which is not `acosh u` for `u ≤ −1`.
//! - `ExpToPower` (`exp2pow`): `e^(c ln b + r) = b^c e^r`, every term of the
//!   exponent with a factor `Log[b]` a power: `exp(2*ln(x)+3*ln(y))` is
//!   `x^2*y^3`.
//! - `PowerToExp` (`pow2exp`): `b^w = e^(w ln b)` for an exponent that is
//!   not a rational number, as giac keeps `x^2` and `sqrt(x)`; `E^w` is
//!   `Exp[w]`.
//! - `ExpExpand` (`expexpand`): `e^(Σ t) = Π e^t`, with a term's exact
//!   numeric coefficient split as `k c`: `e^(k c u) = (e^(c u))^k`, `k` an
//!   integer and `c` primitive with its first nonzero part positive. So
//!   `exp(-2*x)` is `exp(x)^-2` and `exp(-i*x)` is `exp(i*x)^-1`. A number
//!   stays one factor.
//! - `LogExpand` (`lnexpand`), formal: the logarithm of a product is the sum
//!   of the factors', `ln(b^r) = r ln b`, with `ln|b|` for a rational `r`
//!   of even numerator, `ln(p/q) = ln p − ln q` and `ln(e^u) = u`. An integer
//!   is not factored; a product with a number that is not a positive
//!   rational stays, as giac's `ln(-2*x)` does.
//! - `LogCombine` (`lncollect`), formal: in a sum, the terms `r c Log[b]`,
//!   `r` rational, grouped by `c`, are `(c/d) Log[Π b^(d r)]`, `d` the lcm
//!   of the denominators of the `r`s, the sign out when every `r` is
//!   negative. A group of one term changes only when its numerator is not
//!   ±1: `2*ln(x)` is `ln(x^2)`, `-ln(x)` stays.

use cas_num::{Complex, Integer, Rational};

use super::{factors, head, terms};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::trig::exponential_form;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

fn log(u: Expr) -> Expr {
    Expr::call(Symbol::LOG, [u])
}

fn sqrt(u: &Expr) -> Expr {
    u.pow(Expr::from(Rational::new(Integer::ONE, Integer::from(2))))
}

/// The exponent of `Exp[u]` or `E^u`.
fn exponent(e: &Expr) -> Option<&Expr> {
    match (head(e), e.args()) {
        (Some(Symbol::EXP), [u]) => Some(u),
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => Some(u),
        _ => None,
    }
}

/// The argument of `Log[b]`, the logarithm of one argument.
fn log_of(e: &Expr) -> Option<&Expr> {
    match (head(e), e.args()) {
        (Some(Symbol::LOG), [b]) => Some(b),
        _ => None,
    }
}

/// `f` of the factors, `1` for none and the factor for one, built.
fn product(mut fs: Vec<Expr>) -> Expr {
    match fs.len() {
        0 => Expr::from(1),
        1 => fs.pop().expect("one factor"),
        _ => Expr::call(Symbol::TIMES, fs),
    }
}

/// `HyperbolicToExp`'s rule.
pub(super) fn hyperbolic(_: &mut Eval<'_>, part: &Expr) -> Found {
    if let Some(form) = exponential_form(part, false) {
        return Ok(Some(form));
    }
    let (Some(f), [u]) = (head(part), part.args()) else {
        return Ok(None);
    };
    let one = || Expr::from(1);
    Ok(match f {
        Symbol::ARC_SINH => Some(log(u.clone() + sqrt(&(u.pow(2) + one())))),
        Symbol::ARC_COSH => Some(log(
            u.clone() + sqrt(&(u.clone() + one())) * sqrt(&(u.clone() - one()))
        )),
        Symbol::ARC_TANH => Some(log((one() + u.clone()) / (one() - u.clone())) / Expr::from(2)),
        _ => None,
    })
}

/// A term `c Log[b]` of an exponent as `(b, c)`: its first factor `Log[b]`
/// and the others.
fn log_factor(t: &Expr) -> Option<(Expr, Expr)> {
    if let Some(b) = log_of(t) {
        return Some((b.clone(), Expr::from(1)));
    }
    if head(t) != Some(Symbol::TIMES) {
        return None;
    }
    let fs = t.args();
    let i = fs.iter().position(|f| log_of(f).is_some())?;
    let others: Vec<Expr> = fs[..i].iter().chain(&fs[i + 1..]).cloned().collect();
    Some((fs[i].args()[0].clone(), product(others)))
}

/// `ExpToPower`'s rule.
pub(super) fn exp_to_power(_: &mut Eval<'_>, part: &Expr) -> Found {
    let Some(u) = exponent(part) else {
        return Ok(None);
    };
    let mut powers = Vec::new();
    let mut rest = Vec::new();
    for t in terms(u) {
        match log_factor(t) {
            Some((b, c)) => powers.push(b.pow(c)),
            None => rest.push(t.clone()),
        }
    }
    if powers.is_empty() {
        return Ok(None);
    }
    if !rest.is_empty() {
        let r = Expr::call(Symbol::PLUS, rest);
        powers.push(match head(part) {
            Some(Symbol::EXP) => Expr::call(Symbol::EXP, [r]),
            _ => Expr::from(Symbol::E).pow(r),
        });
    }
    Ok(Some(product(powers)))
}

/// `PowerToExp`'s rule.
pub(super) fn power_to_exp(_: &mut Eval<'_>, part: &Expr) -> Found {
    let (Some(Symbol::POWER), [b, w]) = (head(part), part.args()) else {
        return Ok(None);
    };
    if b.as_symbol() == Some(Symbol::E) {
        return Ok(Some(Expr::call(Symbol::EXP, [w.clone()])));
    }
    if w.to_rational().is_some() {
        return Ok(None);
    }
    Ok(Some(Expr::call(Symbol::EXP, [w.clone() * log(b.clone())])))
}

/// A term `k c u` of an exponent, `k c` its exact numeric coefficient, as
/// `(k, c u)`: `k` an integer, `c` with parts of coprime numerators and its
/// first nonzero part positive; `None` for `k = 1` and for a term without
/// such a coefficient.
fn multiple(term: &Expr) -> Option<(Integer, Expr)> {
    let (Some(Symbol::TIMES), [number, rest @ ..]) = (head(term), term.args()) else {
        return None;
    };
    let (re, im) = match number.view() {
        View::Complex(z) => (z.re().clone(), z.im().clone()),
        _ => (number.to_rational()?, Rational::ZERO),
    };
    let gcd = re.numer().gcd(im.numer());
    let negative = if re.is_zero() {
        im.is_negative()
    } else {
        re.is_negative()
    };
    let k = if negative { -gcd } else { gcd };
    if k == Integer::ONE {
        return None;
    }
    let scale = Rational::from(k.clone());
    let c = Complex::new(&re / &scale, &im / &scale);
    let c = if c.is_real() {
        Expr::from(c.re().clone())
    } else {
        Expr::from(c)
    };
    let mut fs = vec![c];
    fs.extend(rest.iter().cloned());
    Some((k, Expr::call(Symbol::TIMES, fs)))
}

/// `ExpExpand`'s rule.
pub(super) fn exp_expand(_: &mut Eval<'_>, part: &Expr) -> Found {
    let Some(u) = exponent(part) else {
        return Ok(None);
    };
    let exp = |v: Expr| Expr::call(Symbol::EXP, [v]);
    let ts = terms(u);
    let mut out = Vec::with_capacity(ts.len());
    let mut changed = ts.len() > 1;
    for t in ts {
        match multiple(t) {
            Some((k, v)) => {
                changed = true;
                out.push(exp(v).pow(Expr::from(k)));
            }
            None => out.push(exp(t.clone())),
        }
    }
    Ok(changed.then(|| product(out)))
}

/// The logarithm of one factor expanded, built; `None` for one that stays
/// `Log[f]`.
fn factor_log(f: &Expr) -> Option<Expr> {
    if let Some(u) = exponent(f) {
        return Some(u.clone());
    }
    if let Some(q) = f.to_rational() {
        return (!q.is_integer())
            .then(|| log(Expr::from(q.numer().clone())) - log(Expr::from(q.denom().clone())));
    }
    let (Some(Symbol::POWER), [b, r]) = (head(f), f.args()) else {
        return None;
    };
    let even = r
        .to_rational()
        .is_some_and(|q| !q.numer().unsigned_abs().bit(0));
    Some(if even {
        r.clone() * log(Expr::call(Symbol::ABS, [b.clone()]))
    } else {
        r.clone() * expanded_log(b).unwrap_or_else(|| log(b.clone()))
    })
}

/// `ln z` expanded, built; `None` where nothing expands.
fn expanded_log(z: &Expr) -> Option<Expr> {
    let fs = factors(z);
    // A number that is not a positive rational keeps the product whole.
    if fs.len() > 1
        && fs
            .iter()
            .any(|f| f.is_number() && !f.to_rational().is_some_and(|q| q.is_positive()))
    {
        return None;
    }
    let expanded: Vec<Option<Expr>> = fs.iter().map(factor_log).collect();
    if fs.len() == 1 && expanded[0].is_none() {
        return None;
    }
    let sum: Vec<Expr> = fs
        .iter()
        .zip(expanded)
        .map(|(f, e)| e.unwrap_or_else(|| log(f.clone())))
        .collect();
    Some(Expr::call(Symbol::PLUS, sum))
}

/// `LogExpand`'s rule.
pub(super) fn log_expand(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(log_of(part).and_then(expanded_log))
}

/// A log term `r c Log[b]`, `r` its rational coefficient, as `(r, c, b)`;
/// `None` for a term with no logarithm or more than one.
fn log_term(t: &Expr) -> Option<(Rational, Expr, Expr)> {
    if let Some(b) = log_of(t) {
        return Some((Rational::ONE, Expr::from(1), b.clone()));
    }
    if head(t) != Some(Symbol::TIMES) {
        return None;
    }
    let fs = t.args();
    let mut logs = fs.iter().enumerate().filter(|(_, f)| log_of(f).is_some());
    let (i, _) = logs.next()?;
    if logs.next().is_some() {
        return None;
    }
    let (r, from) = match fs[0].to_rational() {
        Some(q) if i != 0 => (q, 1),
        _ => (Rational::ONE, 0),
    };
    let others: Vec<Expr> = (from..fs.len())
        .filter(|&j| j != i)
        .map(|j| fs[j].clone())
        .collect();
    Some((r, product(others), fs[i].args()[0].clone()))
}

/// The log terms of one factor `c`: `(c, [(r, b, the term)])`.
type Group<'a> = (Expr, Vec<(Rational, Expr, &'a Expr)>);

/// The terms with their logarithms combined by group, built; `None` when no
/// group changes.
fn combined(ts: &[Expr]) -> Option<Expr> {
    let mut groups: Vec<Group<'_>> = Vec::new();
    let mut out = Vec::new();
    for t in ts {
        match log_term(t) {
            Some((r, c, b)) => match groups.iter_mut().find(|(d, _)| *d == c) {
                Some((_, members)) => members.push((r, b, t)),
                None => groups.push((c, vec![(r, b, t)])),
            },
            None => out.push(t.clone()),
        }
    }
    let mut changed = false;
    for (c, members) in groups {
        if let [(r, _, t)] = members.as_slice()
            && *r.abs().numer() == Integer::ONE
        {
            out.push((*t).clone());
            continue;
        }
        changed = true;
        let d = members
            .iter()
            .fold(Integer::ONE, |d, (r, _, _)| d.lcm(r.denom()));
        let negative = members.iter().all(|(r, _, _)| r.is_negative());
        let sign = Rational::from(if negative { -1 } else { 1 });
        let scale = Rational::from(d.clone()) * &sign;
        let argument: Vec<Expr> = members
            .iter()
            .map(|(r, b, _)| b.pow(Expr::from(r * &scale)))
            .collect();
        let coefficient = Expr::from(sign / Rational::from(d));
        out.push(coefficient * c * log(product(argument)));
    }
    changed.then(|| Expr::call(Symbol::PLUS, out))
}

/// `LogCombine`'s rule: at sums, and at a product that is one log term.
pub(super) fn log_combine(_: &mut Eval<'_>, part: &Expr) -> Found {
    Ok(match head(part) {
        Some(Symbol::PLUS) => combined(part.args()),
        Some(Symbol::TIMES) => combined(std::slice::from_ref(part)),
        _ => None,
    })
}

#[cfg(test)]
mod tests;
