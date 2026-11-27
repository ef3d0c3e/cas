//! The continuous Hilbert transform (CAS-648), giac's `hilbert(f, t)` and
//! `hilbert(f)` (variable `x`), its `hilbtrans` (signalprocessing.cc:
//! 4911-4955): `HilbertMatrix[f, t]` for a symbol `t`, and
//! `HilbertMatrix[f]` for `f` neither an integer (the Hilbert matrix) nor a
//! list (the analytic signal, [`super::hilbert`]).
//!
//! A proper rational function first, by residues (CAS-654): `f = p/q`
//! (`Together`) with exact coefficients over ℚ or ℚ(i), `deg p < deg q`,
//! the poles `a` and their multiplicities `RootsWithMultiplicity`'s (each
//! factor of `FactorList` solved by `CSolve`). With `c_{a,k}` the
//! principal part of `f` at `a` (`p(a)/q′(a)` for a simple pole),
//! `H f = Σ_a s_a i Σ_k c_{a,k}/(t − a)^k`, `s_a = ±1` the sign of `Im a`:
//! giac's `t/(t^2+1)` for `1/(1+t^2)`, put together real by `ComplexExpand`
//! and `Together`. The rule passes when the roots fall short of `deg q`, a
//! root holds a float or a `RootOf`, or the sign of `Im a` is not exactly
//! `±1` (a real pole, an undecidable sign); so do symbolic coefficients.
//!
//! Otherwise by the Fourier multiplier `−i sign(ω)`: `f` free of `t` is 0; else `F =
//! FourierTransform[f, t, ω]`, `G = InverseFourierTransform[−i sign(ω) F,
//! ω, t]`, the product expanded first so the transforms' rules see each
//! term (`H cos = sin`, `H sin = −cos`). An undefined or infinite `F` or
//! `G` is `Indeterminate`, giac's `undef`; so is a derivative of a delta at
//! `ω = 0` in `F` (`sign(ω) δ⁽ⁿ⁾(ω)` has no value: a polynomial `f`, giac's
//! `undef` for `hilbert(t^2,t)`). A transform that stays gives giac's inert
//! `Hilbert[f, t]` (`hilbert(exp(-t^2),t)`, an unknown `f(t)`).
//! Departures: where our `FourierTransform` or `InverseFourierTransform`
//! leaves a term giac closes (a symbolic scale `sin(a t)`; `sinc`; the
//! box), ours is the inert call.

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

/// The exact coefficients of the polynomial `g` in `t`, leading first;
/// `None` for another `g` or a coefficient that is no exact number.
fn exact_coefficients(
    ev: &mut Eval<'_>,
    g: &Expr,
    t: &Expr,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let cs = ev.eval(&Expr::call(
        Symbol::COEFFICIENT_LIST,
        [g.clone(), t.clone()],
    ))?;
    if cs.head() != Some(&Expr::from(Symbol::LIST)) || !cs.args().iter().all(Expr::is_exact_number)
    {
        return Ok(None);
    }
    Ok(Some(cs.args().to_vec()))
}

/// The residue rule (CAS-654): `f = p/q` over ℚ or ℚ(i), `deg p < deg q`,
/// every pole `a` exact in radicals and off the real line. `H f = Σ_a s_a i
/// Σ_k c_{a,k}/(t − a)^k`, `s_a = sign(Im a)`, `c_{a,k}` the principal part
/// of `f` at `a`; `None` where the rule does not apply.
fn residue_rule(ev: &mut Eval<'_>, f: &Expr, t: &Expr) -> Native {
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [f.clone()]))?;
    let p = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let q = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    let (Some(ps), Some(qs)) = (
        exact_coefficients(ev, &p, t)?,
        exact_coefficients(ev, &q, t)?,
    ) else {
        return Ok(None);
    };
    if qs.len() < 2 || ps.len() >= qs.len() {
        return Ok(None);
    }
    let degree = i64::try_from(qs.len() - 1).expect("small");
    let lead = qs[0].clone();
    // The poles, each with its multiplicity and the sign of its imaginary part.
    let found = ev.eval(&Expr::call(
        Symbol::ROOTS_WITH_MULTIPLICITY,
        [q.clone(), t.clone()],
    ))?;
    if found.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let mut poles: Vec<(Expr, i64, i64)> = Vec::new();
    for pair in found.args() {
        let [a, m] = pair.args() else {
            return Ok(None);
        };
        let Some(m) = m.as_i64().filter(|m| *m >= 1) else {
            return Ok(None);
        };
        let inexact = a
            .distinct()
            .any(|g| matches!(g.view(), View::Ball(_) | View::ComplexBall(_)));
        if inexact || a.contains(&Expr::from(Symbol::ROOT_OF)) || a.contains(t) {
            return Ok(None);
        }
        let side = ev.eval(&Expr::call(
            Symbol::SIGN,
            [Expr::call(
                Symbol::IM,
                [Expr::call(Symbol::COMPLEX_EXPAND, [a.clone()])],
            )],
        ))?;
        let Some(s) = side.as_i64().filter(|s| s.abs() == 1) else {
            return Ok(None);
        };
        poles.push((a.clone(), m, s));
    }
    if poles.iter().map(|(_, m, _)| m).sum::<i64>() != degree {
        return Ok(None);
    }
    let i = Expr::from(Symbol::I);
    let mut terms = Vec::new();
    for (a, m, s) in &poles {
        let at = [(t.clone(), a.clone())];
        if *m == 1 {
            // A simple pole: p(a)/q'(a).
            let dq = ev.eval(&Expr::call(Symbol::D, [q.clone(), t.clone()]))?;
            let c = ev.eval(&(p.clone() / dq).substitute(&at))?;
            terms.push(Expr::from(*s) * i.clone() * c / (t.clone() - a.clone()));
            continue;
        }
        // g_a = p/(lc ∏_{b≠a} (t − b)^{m_b}); c_{a,k} = g_a^{(m−k)}(a)/(m−k)!.
        let others: Vec<Expr> = poles
            .iter()
            .filter(|(b, _, _)| b != a)
            .map(|(b, mb, _)| (t.clone() - b.clone()).pow(*mb))
            .collect();
        let mut g = p.clone() / (lead.clone() * Expr::call(Symbol::TIMES, others));
        let mut factorial = Expr::from(1);
        for n in 0..*m {
            if n > 0 {
                g = ev.eval(&Expr::call(Symbol::D, [g, t.clone()]))?;
                factorial *= Expr::from(n);
            }
            let c = ev.eval(&(g.substitute(&at) / factorial.clone()))?;
            terms.push(Expr::from(*s) * i.clone() * c / (t.clone() - a.clone()).pow(m - n));
        }
    }
    let expanded = ev.eval(&Expr::call(
        Symbol::COMPLEX_EXPAND,
        [Expr::call(Symbol::PLUS, terms)],
    ))?;
    let mut result = ev.eval(&Expr::call(Symbol::TOGETHER, [expanded]))?;
    if result.contains(&i) {
        result = ev.eval(&Expr::call(Symbol::SIMPLIFY, [result]))?;
    }
    Ok(Some(result))
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
    if let Some(g) = residue_rule(ev, f, t)? {
        ev.step("hilbert", e, &g);
        return Ok(Some(g));
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
