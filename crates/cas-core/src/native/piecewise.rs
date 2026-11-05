//! Piecewise functions (CAS-DOC-53 §1, CAS-513): giac's
//! `piecewise(c₁, v₁, …, c_k, v_k, [d])` as `Piecewise`.
//!
//! The pairs are scanned in order: a condition `False` drops its pair, an
//! undecided one keeps it, and `True` ends the scan, its value the answer
//! when no undecided condition came before and the new default otherwise.
//! With no pair left the default is the answer, and with none the call is
//! `Indeterminate`, with a warning (giac's `"No case applies"`).
//!
//! Derivatives and antiderivatives go branch by branch, the conditions
//! kept (CAS-514): `D` through [`piecewise_map`], `Integrate[f, x]` in
//! `integrate`, which leaves the call when a branch does not integrate.
//! The antiderivative need not be continuous at the breakpoints, as
//! giac's.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// The condition-value pairs of `Piecewise[…]` and its default.
pub(crate) fn branches(args: &[Expr]) -> (Vec<(&Expr, &Expr)>, Option<&Expr>) {
    let pairs = args.chunks_exact(2).map(|p| (&p[0], &p[1])).collect();
    let default = (args.len() % 2 == 1).then(|| &args[args.len() - 1]);
    (pairs, default)
}

/// `Piecewise[c₁, v₁, …]` from pairs and a default.
pub(crate) fn piecewise(pairs: &[(Expr, Expr)], default: Option<Expr>) -> Expr {
    let mut args = Vec::with_capacity(2 * pairs.len() + 1);
    for (c, v) in pairs {
        args.push(c.clone());
        args.push(v.clone());
    }
    args.extend(default);
    Expr::call(Symbol::PIECEWISE, args)
}

/// `Piecewise[…]`.
pub(super) fn piecewise_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (pairs, mut default) = branches(e.args());
    let mut kept: Vec<(Expr, Expr)> = Vec::with_capacity(pairs.len());
    let mut changed = false;
    for (c, v) in pairs {
        match c.as_symbol() {
            Some(Symbol::FALSE) => changed = true,
            Some(Symbol::TRUE) => {
                changed = true;
                default = Some(v);
                break;
            }
            _ => kept.push((c.clone(), v.clone())),
        }
    }
    if kept.is_empty() {
        if let Some(d) = default {
            return Ok(Some(d.clone()));
        }
        ev.warn(Warning::Native {
            message: "no case applies".into(),
            expr: e.clone(),
        });
        return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
    }
    Ok(changed.then(|| piecewise(&kept, default.cloned())))
}

/// `PiecewiseMap[f, p]`, without a surface name: `f` applied to each value
/// and the default of the `Piecewise` `p`, its conditions kept (CAS-514);
/// `f[p]` for what is not a `Piecewise`. `D` maps through it.
pub(super) fn piecewise_map(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, p] = e.args() else {
        return Ok(None);
    };
    if p.head() != Some(&Expr::from(Symbol::PIECEWISE)) {
        return Ok(Some(Expr::apply(f, std::slice::from_ref(p))));
    }
    let (pairs, default) = branches(p.args());
    let pairs: Vec<(Expr, Expr)> = pairs
        .into_iter()
        .map(|(c, v)| (c.clone(), Expr::apply(f, std::slice::from_ref(v))))
        .collect();
    let default = default.map(|d| Expr::apply(f, std::slice::from_ref(d)));
    Ok(Some(piecewise(&pairs, default)))
}
