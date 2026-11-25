//! HP Prime's least squares and decompositions in giac (CAS-623, giac
//! 2.0.0.21 rpn.cc `_LSQ`, `_LQ`, `_SVD`), on this module's exact
//! decompositions:
//!
//! - **`LeastSquares[A, b]`** (`LSQ`, `lsq`): the `x` minimizing
//!   `‖A x − b‖`, `b` a vector or a matrix (a column per right-hand side),
//!   answered as a matrix of columns, `[[1],[2]]` for a vector `b`, as
//!   giac's. The normal equations `Aᴴ A x = Aᴴ b` solved by `LinSolve`,
//!   exact on exact input (giac's QR on decimals, exact otherwise). `A` must
//!   have at least as many rows as columns and independent columns: giac's
//!   smallest solution of an under-determined system leaves the call here,
//!   as does a singular `Aᴴ A` (giac "Unable to transpose").
//! - **`LQDecomposition[A]`** (`LQ`): `{L, Q, P}` with `A = L Q`, from
//!   `QRDecomposition` of `Aᵀ = Q₁ R₁`: `L = R₁ᵀ` lower triangular, `Q =
//!   Q₁ᵀ` with orthonormal rows, and `P` the identity of `A`'s column count,
//!   giac's (`v[2]`, never a permutation). giac's Householder `L` and `Q`
//!   are decimals and may differ from ours by the signs of `Q`'s rows; the
//!   products agree. Ours is exact and the reduced form: `L` is `m × m` and
//!   `Q` `m × n` where giac's are `m × n` and `n × n`. `A` with more rows
//!   than columns (giac's error too) or dependent rows (giac answers) leaves
//!   the call.
//! - **`AscendingSingularValueDecomposition[A]`** (`SVD`): giac's `SVD`,
//!   `{U, s, V}` with `s` ascending and `A = U diag(s) Vᴴ`: `svd`'s
//!   (`SingularValueDecomposition`) with `s` and the columns of `U` and `V`
//!   reversed. giac's `SVD` and `svd` differ in that order, so `SVD` is no
//!   alias; `SVL` is `svl`'s alias, one function in giac.

use super::{inner, list, matrix, qr_decomposition, rows};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::matrix::decompose::unitary::singular_value_decomposition;
use crate::native::matrix::warn;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn column(a: &[Vec<Expr>], j: usize) -> Vec<Expr> {
    a.iter().map(|r| r[j].clone()).collect()
}

fn transposed(a: &[Vec<Expr>]) -> Vec<Vec<Expr>> {
    (0..a[0].len()).map(|j| column(a, j)).collect()
}

/// `LeastSquares[A, b]`.
pub(in crate::native) fn least_squares(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    let (m, n) = (a.len(), a[0].len());
    // The right-hand sides, as columns.
    let sides: Vec<Vec<Expr>> = if let Some(b) = rows(b) {
        if b.len() != m {
            return Ok(None);
        }
        transposed(&b)
    } else if b.head() == Some(&Expr::from(Symbol::LIST))
        && b.args().len() == m
        && b.args()
            .iter()
            .all(|x| x.head() != Some(&Expr::from(Symbol::LIST)))
    {
        vec![b.args().to_vec()]
    } else {
        return Ok(None);
    };
    if m < n {
        warn(ev, "fewer equations than unknowns", e);
        return Ok(None);
    }
    let columns: Vec<Vec<Expr>> = (0..n).map(|j| column(&a, j)).collect();
    let mut normal = Vec::with_capacity(n);
    for ci in &columns {
        let mut row = Vec::with_capacity(n);
        for cj in &columns {
            let Some(x) = inner(ev, ci, cj)? else {
                return Ok(None);
            };
            row.push(x);
        }
        normal.push(row);
    }
    let normal = matrix(normal);
    let mut solutions = Vec::with_capacity(sides.len());
    for side in &sides {
        let mut c = Vec::with_capacity(n);
        for ci in &columns {
            let Some(x) = inner(ev, ci, side)? else {
                return Ok(None);
            };
            c.push(x);
        }
        let x = ev.eval(&Expr::call(Symbol::LIN_SOLVE, [normal.clone(), list(c)]))?;
        if x.head() != Some(&Expr::from(Symbol::LIST)) || x.args().len() != n {
            return Ok(None);
        }
        solutions.push(x.args().to_vec());
    }
    let result = matrix(transposed(&solutions));
    ev.step("lsq", e, &result);
    Ok(Some(result))
}

/// `LQDecomposition[A]`.
pub(in crate::native) fn lq_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    let n = a[0].len();
    let qr = Expr::call(Symbol::QR_DECOMPOSITION, [matrix(transposed(&a))]);
    let Some(qr) = qr_decomposition(ev, &qr)? else {
        return Ok(None);
    };
    let [q, r] = qr.args() else {
        return Ok(None);
    };
    let (Some(q), Some(r)) = (rows(q), rows(r)) else {
        return Ok(None);
    };
    let identity = matrix(super::identity(n));
    let result = list(vec![
        matrix(transposed(&r)),
        matrix(transposed(&q)),
        identity,
    ]);
    ev.step("LQ", e, &result);
    Ok(Some(result))
}

/// The matrix with its columns in reverse order.
fn reversed_columns(m: &Expr) -> Option<Expr> {
    let rows = rows(m)?;
    Some(matrix(
        rows.into_iter()
            .map(|mut r| {
                r.reverse();
                r
            })
            .collect(),
    ))
}

/// `AscendingSingularValueDecomposition[A]`.
pub(in crate::native) fn ascending_svd(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let svd = Expr::call(Symbol::SINGULAR_VALUE_DECOMPOSITION, e.args().to_vec());
    let Some(svd) = singular_value_decomposition(ev, &svd)? else {
        return Ok(None);
    };
    let [u, s, v] = svd.args() else {
        return Ok(None);
    };
    let (Some(u), Some(v)) = (reversed_columns(u), reversed_columns(v)) else {
        return Ok(None);
    };
    let mut s = s.args().to_vec();
    s.reverse();
    let result = list(vec![u, list(s), v]);
    ev.step("SVD", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
