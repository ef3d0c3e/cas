//! The continuous Hilbert transform (CAS-648), giac's `hilbert(f, t)` and
//! `hilbert(f)` (variable `x`), its `hilbtrans` (signalprocessing.cc:
//! 4911-4955): `HilbertMatrix[f, t]` for a symbol `t`, and
//! `HilbertMatrix[f]` for `f` neither an integer (the Hilbert matrix) nor a
//! list (the analytic signal, [`super::hilbert`]).
//!
//! By the Fourier multiplier `−i sign(ω)`: `f` free of `t` is 0; else `F =
//! FourierTransform[f, t, ω]`, `G = InverseFourierTransform[−i sign(ω) F,
//! ω, t]`, the product expanded first so the transforms' rules see each
//! term (`H cos = sin`, `H sin = −cos`). An undefined or infinite `F` or
//! `G` is `Indeterminate`, giac's `undef`; so is a derivative of a delta at
//! `ω = 0` in `F` (`sign(ω) δ⁽ⁿ⁾(ω)` has no value: a polynomial `f`, giac's
//! `undef` for `hilbert(t^2,t)`). A transform that stays gives giac's inert
//! `Hilbert[f, t]` (`hilbert(exp(-t^2),t)`, an unknown `f(t)`).
//! Departures: where our `FourierTransform` or `InverseFourierTransform`
//! leaves a term giac closes (rational functions, whose `F` holds
//! `e^(−|ω|)`; a symbolic scale `sin(a t)`; `sinc`; the box), ours is the
//! inert call.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// Whether `f` holds `δ⁽ⁿ⁾(u)`, `n ≥ 1`, with `u` zero at `w = 0`.
fn delta_derivative_at_zero(ev: &mut Eval<'_>, f: &Expr, w: &Expr) -> Result<bool, EvalError> {
    let arguments: Vec<Expr> = f
        .distinct()
        .filter_map(|g| match g.view() {
            // Derivative[n][DiracDelta][u].
            View::Apply { head, args: [u] } => match head.view() {
                View::Apply {
                    head: d,
                    args: [function],
                } if function.as_symbol() == Some(Symbol::DIRAC_DELTA)
                    && d.head().and_then(Expr::as_symbol) == Some(Symbol::DERIVATIVE) =>
                {
                    Some(u.clone())
                }
                _ => None,
            },
            // DiracDelta[u, n].
            View::Apply { head, args: [u, _] } if head.as_symbol() == Some(Symbol::DIRAC_DELTA) => {
                Some(u.clone())
            }
            _ => None,
        })
        .collect();
    for u in arguments {
        if ev.eval(&u.substitute(&[(w.clone(), Expr::from(0))]))? == Expr::from(0) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `e` holds an undefined or infinite value.
fn undefined(e: &Expr) -> bool {
    [
        Symbol::INDETERMINATE,
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
        Symbol::INFINITY,
    ]
    .into_iter()
    .any(|s| e.contains(&Expr::from(s)))
}

/// `HilbertMatrix[f, t]`, giac's `hilbtrans`.
pub(in crate::native) fn hilbert_transform(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    t: &Expr,
) -> Native {
    if t.as_symbol().is_none_or(Symbol::is_builtin) {
        return Ok(None);
    }
    if !f.contains(t) {
        return Ok(Some(Expr::from(0)));
    }
    let inert = || Expr::call(Symbol::HILBERT, [f.clone(), t.clone()]);
    let w = Expr::from(Symbol::new("Hilbert$w"));
    let big_f = ev.eval(&Expr::call(
        Symbol::FOURIER_TRANSFORM,
        [f.clone(), t.clone(), w.clone()],
    ))?;
    let result = if undefined(&big_f) || delta_derivative_at_zero(ev, &big_f, &w)? {
        Expr::from(Symbol::INDETERMINATE)
    } else if big_f.contains(&Expr::from(Symbol::FOURIER_TRANSFORM)) {
        inert()
    } else {
        let multiplied = ev.eval(&Expr::call(
            Symbol::EXPAND,
            [-Expr::from(Symbol::I) * Expr::call(Symbol::SIGN, [w.clone()]) * big_f],
        ))?;
        let g = ev.eval(&Expr::call(
            Symbol::INVERSE_FOURIER_TRANSFORM,
            [multiplied, w.clone(), t.clone()],
        ))?;
        if undefined(&g) {
            Expr::from(Symbol::INDETERMINATE)
        } else if g.contains(&Expr::from(Symbol::FOURIER_TRANSFORM))
            || g.contains(&Expr::from(Symbol::INVERSE_FOURIER_TRANSFORM))
        {
            inert()
        } else {
            g
        }
    };
    ev.step("hilbert", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
