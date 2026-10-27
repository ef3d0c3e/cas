//! Linear equations with constant coefficients (CAS-DOC-16 §5, CAS-261),
//! `Σ aₖ y⁽ᵏ⁾ = f` with every `aₖ` free of `x`.
//!
//! The roots of the characteristic polynomial `Σ aₖ rᵏ`, with their
//! multiplicities, as `Eigenvalues` takes them (`roots_counted`, CAS-221),
//! give the basis: `x^j e^(r x)` for `j` below the multiplicity of `r`; with
//! real coefficients a root `a + b i`, `b` not 0, whose conjugate is a root
//! too gives `x^j e^(a x) cos(b x)` and `x^j e^(a x) sin(b x)` for both.
//! The constants follow the roots' order, `c_0` the first: `y″ + y = 0` is
//! `c_0 cos(x) + c_1 sin(x)`, as giac's.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equation: Σ aₖ y⁽ᵏ⁾ = f, roots r"
)]

use cas_num::Complex;

use super::super::parts::parts;
use super::super::solve::roots_counted;
use super::particular::particular;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The general solution of `Σ aₖ y⁽ᵏ⁾ = f` with constant coefficients, a
/// right side by [`particular`]; `None` when a coefficient holds `x`, or
/// the roots or a particular solution are not found.
pub(super) fn constant_coefficients(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    f: &Expr,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    if a.iter().any(|ak| ak.contains(x)) {
        return Ok(None);
    }
    let Some((basis, roots)) = basis(ev, e, a, x)? else {
        return Ok(None);
    };
    let mut terms: Vec<Expr> = constants
        .iter()
        .zip(&basis)
        .map(|(c, b)| c.clone() * b.clone())
        .collect();
    if *f != Expr::from(0) {
        let Some(y) = particular(ev, e, a, f, x, &basis, &roots)? else {
            return Ok(None);
        };
        terms.push(y);
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// The basis of `Σ aₖ y⁽ᵏ⁾ = 0` from the roots of `Σ aₖ rᵏ`, and the
/// distinct roots with their multiplicities; `None` when they are not
/// found.
#[allow(clippy::type_complexity, reason = "the basis and the roots")]
fn basis(
    ev: &mut Eval<'_>,
    e: &Expr,
    a: &[Expr],
    x: &Expr,
) -> Result<Option<(Vec<Expr>, Vec<(Expr, usize)>)>, EvalError> {
    let r = Expr::from(Symbol::new("DSolve$r"));
    let polynomial = Expr::call(
        Symbol::PLUS,
        a.iter().zip(0i64..).map(|(ak, k)| ak.clone() * r.pow(k)),
    );
    let polynomial = ev.eval(&polynomial)?;
    let Some(roots) = roots_counted(ev, e, &polynomial, &r)? else {
        return Ok(None);
    };
    if roots.len() + 1 != a.len() {
        return Ok(None);
    }
    // The distinct roots, each with its multiplicity, in their order.
    let mut distinct: Vec<(Expr, usize)> = Vec::new();
    for root in roots {
        match distinct.iter_mut().find(|(s, _)| *s == root) {
            Some((_, m)) => *m += 1,
            None => distinct.push((root, 1)),
        }
    }
    let real = a
        .iter()
        .all(|ak| parts(ev.context(), ak).is_some_and(|(_, im)| im == Expr::from(0)));
    let mut used = vec![false; distinct.len()];
    let mut basis = Vec::with_capacity(a.len() - 1);
    for i in 0..distinct.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let (root, m) = distinct[i].clone();
        let pair = if real {
            conjugate(ev, &root, &distinct, &used)?
        } else {
            None
        };
        for j in 0..m {
            let power = x.pow(i64::try_from(j).unwrap_or(i64::MAX));
            match &pair {
                Some((re, im, _)) => {
                    let scale = power * Expr::from(Symbol::E).pow(re.clone() * x.clone());
                    let angle = im.clone() * x.clone();
                    basis.push(scale.clone() * Expr::call(Symbol::COS, [angle.clone()]));
                    basis.push(scale * Expr::call(Symbol::SIN, [angle]));
                }
                None => {
                    basis.push(power * Expr::from(Symbol::E).pow(root.clone() * x.clone()));
                }
            }
        }
        if let Some((_, _, k)) = pair {
            used[k] = true;
        }
    }
    let mut evaluated = Vec::with_capacity(basis.len());
    for b in basis {
        evaluated.push(ev.eval(&b)?);
    }
    Ok(Some((evaluated, distinct)))
}

/// `(a, b, k)` for a root `a + b i`, `b` not 0, whose conjugate `a − b i`
/// is the unused root `k`.
fn conjugate(
    ev: &mut Eval<'_>,
    root: &Expr,
    distinct: &[(Expr, usize)],
    used: &[bool],
) -> Result<Option<(Expr, Expr, usize)>, EvalError> {
    let Some((re, im)) = parts(ev.context(), root) else {
        return Ok(None);
    };
    let (re, im) = (ev.eval(&re)?, ev.eval(&im)?);
    if im == Expr::from(0) {
        return Ok(None);
    }
    let other = ev.eval(&(re.clone() - Expr::from(Complex::I) * im.clone()))?;
    Ok(distinct
        .iter()
        .enumerate()
        .find(|(k, (s, _))| !used[*k] && *s == other)
        .map(|(k, _)| (re, im, k)))
}
