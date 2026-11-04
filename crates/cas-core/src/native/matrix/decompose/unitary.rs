//! Unitary decompositions (CAS-DOC-51): Schur's form (§2, CAS-496).
//!
//! - **`SchurDecomposition[A]`**, Xcas's `schur`: `{Q, T}` with `Q` unitary,
//!   `T` upper triangular and `A = Q T Qᴴ`. From the Jordan form `A = P J
//!   P⁻¹` (the chains of [`super::chains`], exact) and the QR decomposition
//!   `P = Q R`: `A = Q (R J R⁻¹) Qᴴ`, and `T = R J R⁻¹` is a product of
//!   upper triangular matrices, so that no entry has to be shown 0: its
//!   strictly lower part is 0 by construction. `T`'s diagonal holds the
//!   eigenvalues in the order of the chains, `Eigenvalues`' order. Ours is
//!   the complex form: giac's real one has 2 × 2 blocks for complex
//!   eigenvalues, and its `schur([[1,2],[0,3]])` starts from 3.

use super::{assembled, chains, list, matrix, product, rows, simplified, square_matrix};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

pub(in crate::native) fn schur_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if e.args().len() != 1 {
        return Ok(None);
    }
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let Some(chains) = chains(ev, &a, n)? else {
        return Ok(None);
    };
    let mut block = |_: &mut Eval<'_>, lambda: &Expr, k: usize| -> Result<Expr, EvalError> {
        Ok(match k {
            0 => lambda.clone(),
            1 => Expr::from(1),
            _ => Expr::from(0),
        })
    };
    let (p, j) = assembled(ev, &chains, n, &mut block)?;
    let qr = ev.eval(&Expr::call(Symbol::QR_DECOMPOSITION, [p]))?;
    let [q, r] = qr.args() else {
        return Ok(None);
    };
    if qr.head() != Some(&Expr::from(Symbol::LIST)) || rows(q).is_none() {
        return Ok(None);
    }
    let inverse = ev.eval(&Expr::call(Symbol::INVERSE, [r.clone()]))?;
    if rows(&inverse).is_none() {
        return Ok(None);
    }
    let rj = product(ev, r, &j)?;
    let Some(t) = rows(&product(ev, &rj, &inverse)?) else {
        return Ok(None);
    };
    let mut upper = Vec::with_capacity(n);
    for (i, row) in t.into_iter().enumerate() {
        let mut entries = Vec::with_capacity(n);
        for (k, x) in row.into_iter().enumerate() {
            entries.push(if k < i {
                Expr::from(0)
            } else {
                simplified(ev, x)?
            });
        }
        upper.push(entries);
    }
    let result = list(vec![q.clone(), matrix(upper)]);
    ev.step("schur", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
