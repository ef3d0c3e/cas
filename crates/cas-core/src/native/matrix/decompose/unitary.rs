//! Unitary decompositions (CAS-DOC-51): Schur's form (§2, CAS-496) and
//! the singular value decomposition (§3, CAS-497).
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
//! - **`SingularValueDecomposition[A]`**, Xcas's `svd`, for a square `A`:
//!   `{U, s, V}` with `U` and `V` unitary, `s` the singular values
//!   descending (giac's order; `SingularValues` ascends) and `A = U diag(s)
//!   Vᴴ`. The eigenvectors of `Aᴴ A` (its chains, each of length 1 as it is
//!   Hermitian), ordered by their eigenvalues `σᵢ²` descending, made
//!   orthonormal by `GramSchmidt`, are `V`'s columns: eigenvectors of
//!   distinct eigenvalues are orthogonal already, so that each stays in its
//!   eigenspace. `uᵢ = A vᵢ / σᵢ` for `σᵢ ≠ 0`, and `GramSchmidt` of these
//!   and `e₁, …, eₙ` completes `U`. giac's `svd([[3,0],[4,5]])` has `s =
//!   [3√5, √5]`.
//!
//! Matrices of balls take their eigenvectors from `Eigenvectors`'
//! certified eigenpairs (CAS-DOC-51 §4, CAS-498), through the chains'
//! fallback to it; each eigenvalue must then be shown simple, and each
//! entry of the answer is a ball holding the exact decomposition's of
//! every matrix in `A`'s balls.

use std::cmp::Ordering;

use super::forms::sign;
use super::{
    assembled, chains, identity, is_zero, list, matrix, product, rows, simplified, square_matrix,
};
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

/// The columns of the matrix whose rows are `vectors`.
fn columns(vectors: &[Vec<Expr>]) -> Expr {
    let n = vectors.first().map_or(0, Vec::len);
    matrix(
        (0..n)
            .map(|i| vectors.iter().map(|v| v[i].clone()).collect())
            .collect(),
    )
}

/// `GramSchmidt` of the vectors, as rows; `None` unless it gives `n` of
/// them.
fn orthonormal(
    ev: &mut Eval<'_>,
    vectors: Vec<Expr>,
    n: usize,
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let basis = ev.eval(&Expr::call(Symbol::GRAM_SCHMIDT, [list(vectors)]))?;
    Ok(rows(&basis).filter(|r| r.len() >= n).map(|mut r| {
        r.truncate(n);
        r
    }))
}

pub(in crate::native) fn singular_value_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if e.args().len() != 1 {
        return Ok(None);
    }
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let adjoint = Expr::call(
        Symbol::CONJUGATE,
        [Expr::call(Symbol::TRANSPOSE, [a.clone()])],
    );
    let gram = product(ev, &adjoint, &a)?;
    let Some(chains) = chains(ev, &gram, n)? else {
        return Ok(None);
    };
    let mut pairs: Vec<(Expr, Expr)> = Vec::with_capacity(n);
    for (lambda, chain) in chains {
        let [v] = chain.as_slice() else {
            return Ok(None);
        };
        pairs.push((simplified(ev, lambda)?, v.clone()));
    }
    // Descending, by insertion: each comparison decided.
    for i in 1..pairs.len() {
        let mut k = i;
        while k > 0 {
            let gap = simplified(ev, pairs[k].0.clone() - pairs[k - 1].0.clone())?;
            match sign(ev, &gap)? {
                Some(Ordering::Greater) => pairs.swap(k, k - 1),
                Some(_) => break,
                None => return Ok(None),
            }
            k -= 1;
        }
    }
    let Some(v) = orthonormal(ev, pairs.iter().map(|(_, v)| v.clone()).collect(), n)? else {
        return Ok(None);
    };
    let mut s = Vec::with_capacity(n);
    let mut u = Vec::with_capacity(2 * n);
    for ((lambda, _), vi) in pairs.iter().zip(&v) {
        let sigma = simplified(ev, Expr::call(Symbol::SQRT, [lambda.clone()]))?;
        match is_zero(ev, lambda)? {
            Some(false) => {
                let image = product(ev, &a, &list(vi.clone()))?;
                if image.head() != Some(&Expr::from(Symbol::LIST)) || image.args().len() != n {
                    return Ok(None);
                }
                let image = image.args().to_vec();
                let mut ui = Vec::with_capacity(n);
                for x in image {
                    ui.push(simplified(ev, x / sigma.clone())?);
                }
                u.push(list(ui));
            }
            Some(true) => {}
            None => return Ok(None),
        }
        s.push(sigma);
    }
    // Completed only past a zero singular value: of balls, `GramSchmidt`
    // cannot drop the unit vectors a full basis already spans.
    if u.len() < n {
        u.extend(identity(n).into_iter().map(list));
    }
    let Some(u) = orthonormal(ev, u, n)? else {
        return Ok(None);
    };
    let result = list(vec![columns(&u), list(s), columns(&v)]);
    ev.step("svd", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
