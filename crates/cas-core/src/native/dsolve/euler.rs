//! Euler's equations (CAS-DOC-16 §6, CAS-263): `Σ aₖ y⁽ᵏ⁾ = f` with
//! `aₖ = bₖ xᵏ g(x)`, every `bₖ` free of `x`, as `x² y″ + x y′ − y = 0`.
//!
//! The indicial polynomial `Σ bₖ r(r − 1)⋯(r − k + 1)`, its roots with their
//! multiplicities as for constant coefficients ([`super::linear::roots`]): a
//! root `r` of multiplicity `m` gives `x^r ln(x)^j`, `j < m`, and with real
//! `bₖ` a pair `a ± b i` gives `x^a ln(x)^j cos(b ln x)` and `x^a ln(x)^j
//! sin(b ln x)`. A right side by variation of parameters
//! ([`super::particular::variation`]). giac answers "Unable to solve" to
//! some, as `x² y″ + x y′ + y = 0`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equation: Σ bₖ xᵏ y⁽ᵏ⁾ = f, roots r"
)]

use super::linear::{Root, real, roots};
use super::particular::variation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The general solution of an Euler equation `Σ aₖ y⁽ᵏ⁾ = f`; `None` when
/// it is not one, or its roots or particular solution are not found.
pub(super) fn euler(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    f: &Expr,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let n = a.len() - 1;
    let together = |ev: &mut Eval<'_>, e: Expr| ev.eval(&Expr::call(Symbol::TOGETHER, [e]));
    let power = |k: usize| x.pow(i64::try_from(k).unwrap_or(i64::MAX));
    let lead = together(ev, a[n].clone() / power(n))?;
    let mut b = Vec::with_capacity(a.len());
    for (k, ak) in a.iter().enumerate() {
        let bk = together(ev, ak.clone() / (power(k) * lead.clone()))?;
        if bk.contains(x) {
            return Ok(None);
        }
        b.push(bk);
    }
    // Σ bₖ r(r − 1)⋯(r − k + 1).
    let r = Expr::from(Symbol::new("DSolve$r"));
    let falling = |k: usize| {
        Expr::call(
            Symbol::TIMES,
            (0..k).map(|i| r.clone() - Expr::from(i64::try_from(i).unwrap_or(0))),
        )
    };
    let indicial = Expr::call(
        Symbol::PLUS,
        b.iter().enumerate().map(|(k, bk)| bk.clone() * falling(k)),
    );
    let real = real(ev, &b);
    let Some((found, _)) = roots(ev, e, &indicial, &r, n, real)? else {
        return Ok(None);
    };
    let log = Expr::call(Symbol::LOG, [x.clone()]);
    let mut basis = Vec::with_capacity(n);
    for root in found {
        for j in 0..root.multiplicity() {
            let logs = log.pow(i64::try_from(j).unwrap_or(i64::MAX));
            match &root {
                Root::Single(s, _) => basis.push(x.pow(s.clone()) * logs),
                Root::Pair(re, im, _) => {
                    let scale = x.pow(re.clone()) * logs;
                    let angle = im.clone() * log.clone();
                    basis.push(scale.clone() * Expr::call(Symbol::COS, [angle.clone()]));
                    basis.push(scale * Expr::call(Symbol::SIN, [angle]));
                }
            }
        }
    }
    let mut evaluated = Vec::with_capacity(basis.len());
    for y in basis {
        evaluated.push(ev.eval(&y)?);
    }
    let mut terms: Vec<Expr> = constants
        .iter()
        .zip(&evaluated)
        .map(|(c, y)| c.clone() * y.clone())
        .collect();
    if *f != Expr::from(0) {
        let Some(y) = variation(ev, e, a, x, &evaluated, f)? else {
            return Ok(None);
        };
        terms.push(ev.eval(&Expr::call(Symbol::EXPAND, [y]))?);
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}
