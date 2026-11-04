//! Modular arithmetic (CAS-DOC-43): `Modular[a, n]`, giac's `a % n`.
//!
//! - **Values** (§2, CAS-450): `n` an integer `≥ 2`. An integer `a` is its
//!   symmetric residue in `(−n/2, n/2]`, giac's: `5 % 7` is `-2 % 7`; a
//!   rational `p/q` is `p q⁻¹`, `Indeterminate` with a warning when `q` is
//!   not invertible (giac's "Not invertible"); another expression is
//!   expanded and each rational coefficient reduced so, the terms that
//!   vanish dropped: `Modular[(x + 1)^5, 5]` is `Modular[x^5 + 1, 5]`. A
//!   modular value of the same modulus inside goes; of another, the call
//!   stays. An equation is the equation of its sides' values, and a list
//!   the list of its entries' (CAS-DOC-44 §1, §3).
//! - **Arithmetic** (§3): in `Plus` and `Times`, a modular argument takes
//!   in the others, `Modular[a, n] + b` being `Modular[a + b, n]`, so that
//!   giac's `(1 % 3)*x^2+(1 % 3)` reads back as one value; two moduli
//!   leave the call, with giac's "Mod are different". `Power` raises within
//!   the modulus, by repeated squaring for a number, and a negative power
//!   of a number is its inverse's.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of residues: a modulo n, the coefficient c, r its residue"
)]

use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The largest power of a modular expression that is not a number.
const MAX_POWER: i64 = 256;

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// `(a, n)` of a modular value `Modular[a, n]`, `n ≥ 2`.
pub(super) fn parts(e: &Expr) -> Option<(&Expr, Integer)> {
    // Every sum and product asks: the head's symbol, without building one.
    if e.head().and_then(Expr::as_symbol) != Some(Symbol::MODULAR) {
        return None;
    }
    let [a, n] = e.args() else {
        return None;
    };
    let n = n.to_integer().filter(|n| *n >= Integer::from(2))?;
    Some((a, n))
}

fn modular(a: Expr, n: &Integer) -> Expr {
    Expr::call(Symbol::MODULAR, [a, Expr::from(n.clone())])
}

/// The symmetric residue of `q` modulo `n`; `None` when `q`'s denominator
/// is not invertible.
fn residue(q: &Rational, n: &Integer) -> Option<Integer> {
    let (g, inverse, _) = q.denom().gcd_ext(n);
    if g != Integer::ONE {
        return None;
    }
    let (_, r) = (q.numer() * &inverse).div_rem_euclid(n);
    Some(if &(&r + &r) > n { &r - n } else { r })
}

/// `a` with the modular values of modulus `n` in it unwrapped; `None` when
/// one has another modulus.
fn unwrapped(a: &Expr, n: &Integer) -> Option<Expr> {
    if let Some((b, m)) = parts(a) {
        return (m == *n).then(|| b.clone());
    }
    match a.head() {
        Some(h)
            if a.args()
                .iter()
                .any(|x| x.contains(&Expr::from(Symbol::MODULAR))) =>
        {
            let args: Option<Vec<Expr>> = a.args().iter().map(|x| unwrapped(x, n)).collect();
            Some(Expr::apply(h, &args?))
        }
        _ => Some(a.clone()),
    }
}

/// A term as its rational coefficient and the rest.
fn coefficient(t: &Expr) -> (Rational, Expr) {
    if let Some(q) = t.to_rational() {
        return (q, Expr::from(1));
    }
    if t.head() == Some(&Expr::from(Symbol::TIMES))
        && let [c, rest @ ..] = t.args()
        && let Some(q) = c.to_rational()
    {
        return (q, Expr::call(Symbol::TIMES, rest.iter().cloned()));
    }
    (Rational::ONE, t.clone())
}

/// `Modular[a, n]`.
pub(super) fn modular_value(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((a, n)) = parts(e) else {
        return Ok(None);
    };
    // A list modulo n is the list of its entries modulo n, as giac's
    // `[[1,2],[3,4]] % 7` (CAS-DOC-44 §1): `Listable` would want `n` a list
    // too.
    if a.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
        let entries = a.args().iter().map(|x| modular(x.clone(), &n));
        return Ok(Some(Expr::call(Symbol::LIST, entries.collect::<Vec<_>>())));
    }
    // An equation modulo n is one of its sides modulo n (CAS-DOC-44 §3),
    // as giac's `[x+2*y=1,3*x+4*y=2] % 7`.
    if a.head().and_then(Expr::as_symbol) == Some(Symbol::EQUATION)
        && let [l, r] = a.args()
    {
        let sides = [modular(l.clone(), &n), modular(r.clone(), &n)];
        return Ok(Some(Expr::call(Symbol::EQUATION, sides)));
    }
    if let Some(q) = a.to_rational() {
        let Some(r) = residue(&q, &n) else {
            warn(ev, e, "not invertible");
            return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
        };
        let result = modular(Expr::from(r), &n);
        return Ok((result != *e).then_some(result));
    }
    let Some(a) = unwrapped(a, &n) else {
        warn(ev, e, "the moduli are different");
        return Ok(None);
    };
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [a.clone()]))?;
    let terms = match expanded.head() {
        Some(h) if *h == Expr::from(Symbol::PLUS) => expanded.args().to_vec(),
        _ => vec![expanded.clone()],
    };
    let mut kept = Vec::with_capacity(terms.len());
    for t in terms {
        let (c, rest) = coefficient(&t);
        let Some(r) = residue(&c, &n) else {
            warn(ev, e, "not invertible");
            return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
        };
        if !r.is_zero() {
            kept.push(Expr::from(r) * rest);
        }
    }
    let reduced = ev.eval(&Expr::call(Symbol::PLUS, kept))?;
    if reduced == *e.args().first().expect("two arguments") {
        return Ok(None);
    }
    Ok(Some(modular(reduced, &n)))
}

/// The common modulus of the modular arguments of `e`; `None` when there
/// is none, `Some(Err(()))` when they differ.
fn common(e: &Expr) -> Option<Result<Integer, ()>> {
    let mut found: Option<Integer> = None;
    for a in e.args() {
        if let Some((_, n)) = parts(a) {
            match &found {
                Some(m) if *m != n => return Some(Err(())),
                _ => found = Some(n),
            }
        }
    }
    found.map(Ok)
}

/// `Plus` or `Times` of `e`'s arguments with a modular one: the modular
/// value of the whole.
pub(super) fn absorbed(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let n = match common(e) {
        None => return Ok(None),
        Some(Err(())) => {
            warn(ev, e, "the moduli are different");
            return Ok(None);
        }
        Some(Ok(n)) => n,
    };
    let Some(inner) = unwrapped(e, &n) else {
        return Ok(None);
    };
    Ok(Some(modular(inner, &n)))
}

/// `Power[Modular[a, n], k]` for an integer `k`.
pub(super) fn power(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [base, k] = e.args() else {
        return Ok(None);
    };
    let Some((a, n)) = parts(base) else {
        return Ok(None);
    };
    let Some(k) = k.to_integer() else {
        return Ok(None);
    };
    if let Some(q) = a.to_rational() {
        let Some(r) = residue(&q, &n).and_then(|a| a.pow_mod(&k, &n)) else {
            warn(ev, e, "not invertible");
            return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
        };
        return Ok(Some(modular(Expr::from(r), &n)));
    }
    let Some(k) = k.to_i64().filter(|k| (0..=MAX_POWER).contains(k)) else {
        return Ok(None);
    };
    Ok(Some(modular(a.pow(k), &n)))
}

pub(super) mod matrices;
pub(super) mod polynomials;

#[cfg(test)]
mod tests;
