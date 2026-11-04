//! The Hessenberg form and matrix constructions (CAS-DOC-47 §2, CAS-476):
//! Xcas's `hessenberg`, `companion`, `sylvester` and `hadamard`.
//!
//! - **`HessenbergDecomposition[A]`**: `{P, H}` with `A = P H P⁻¹`, `H`
//!   zero below its subdiagonal, by Gaussian similarities column by column:
//!   the pivot below the subdiagonal (the largest in absolute value when
//!   the candidates are rational, as giac's, otherwise the first that is
//!   not 0) swapped into place, rows and columns at once; then each entry
//!   below it cleared by a row operation `E`, followed by `E⁻¹` on the
//!   columns, `P` gathering the inverses. giac's
//!   `hessenberg([[1,2,3],[4,5,6],[7,8,10]])` has `P =
//!   [[1,0,0],[0,4/7,1],[0,1,0]]`.
//! - **`CompanionMatrix[p, x]`**: ones below the diagonal and `−cₖ/cₙ` in
//!   the last column, row `k`, for `p = Σ cₖ xᵏ` of degree `n ≥ 1`.
//! - **`SylvesterMatrix[p, q, x]`**: for degrees `m` and `n`, `n` rows of
//!   `p`'s coefficients from the leading one, each shifted one further,
//!   then `m` of `q`'s: giac's `[[1,0,1],[1,2,0],[0,1,2]]` for `x² + 1` and
//!   `x + 2`. Its determinant is the resultant.
//! - **`HadamardProduct[A, B]`**: entrywise, of one shape.

use std::cmp::Ordering;

use super::{identity, is_zero, list, matrix, rows, simplified, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The index of the pivot among `candidates` (row, entry): the largest
/// rational in absolute value, or the first entry not 0; `Err` when one
/// cannot be told from 0.
fn pivot(
    ev: &mut Eval<'_>,
    candidates: &[(usize, Expr)],
) -> Result<Result<Option<usize>, ()>, EvalError> {
    let rationals: Option<Vec<_>> = candidates
        .iter()
        .map(|(i, x)| x.to_rational().map(|q| (*i, q.abs())))
        .collect();
    if let Some(values) = rationals {
        return Ok(Ok(values
            .into_iter()
            .filter(|(_, q)| !q.is_zero())
            .max_by(|a, b| a.1.cmp(&b.1).then(Ordering::Greater))
            .map(|(i, _)| i)));
    }
    for (i, x) in candidates {
        match is_zero(ev, x)? {
            Some(false) => return Ok(Ok(Some(*i))),
            Some(true) => {}
            None => return Ok(Err(())),
        }
    }
    Ok(Ok(None))
}

/// `HessenbergDecomposition[A]`.
pub(in crate::native) fn hessenberg_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(mut h) = rows(a) else {
        return Ok(None);
    };
    let n = h.len();
    if h[0].len() != n {
        return Ok(None);
    }
    let mut p = identity(n);
    for k in 0..n.saturating_sub(2) {
        let candidates: Vec<(usize, Expr)> = (k + 1..n).map(|i| (i, h[i][k].clone())).collect();
        let Ok(found) = pivot(ev, &candidates)? else {
            warn(ev, "a pivot cannot be told from 0", e);
            return Ok(None);
        };
        let Some(i) = found else {
            continue;
        };
        let s = k + 1;
        if i != s {
            h.swap(i, s);
            for row in h.iter_mut().chain(p.iter_mut()) {
                row.swap(i, s);
            }
        }
        for r in s + 1..n {
            if is_zero(ev, &h[r][k])? == Some(true) {
                continue;
            }
            let f = simplified(ev, h[r][k].clone() / h[s][k].clone())?;
            // Row r −= f row s, then column s += f column r.
            let pivot_row = h[s].clone();
            for (x, y) in h[r].iter_mut().zip(&pivot_row) {
                *x = simplified(ev, x.clone() - f.clone() * y.clone())?;
            }
            for row in h.iter_mut().chain(p.iter_mut()) {
                let v = row[s].clone() + f.clone() * row[r].clone();
                row[s] = simplified(ev, v)?;
            }
        }
    }
    let result = list(vec![matrix(p), matrix(h)]);
    ev.step("hessenberg", e, &result);
    Ok(Some(result))
}

/// The coefficients of `p` in `x`, the leading one first; `None` for no
/// polynomial.
fn coefficients(ev: &mut Eval<'_>, p: &Expr, x: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    let c = ev.eval(&Expr::call(
        Symbol::COEFFICIENT_LIST,
        [p.clone(), x.clone()],
    ))?;
    Ok(
        (c.head() == Some(&Expr::from(Symbol::LIST)) && !c.args().is_empty())
            .then(|| c.args().to_vec()),
    )
}

/// `CompanionMatrix[p, x]`.
pub(in crate::native) fn companion_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, x] = e.args() else {
        return Ok(None);
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let Some(c) = coefficients(ev, p, x)? else {
        return Ok(None);
    };
    let n = c.len() - 1;
    if n == 0 {
        return Ok(None);
    }
    let lead = c[0].clone();
    let mut m = vec![vec![Expr::from(0); n]; n];
    for (k, row) in m.iter_mut().enumerate() {
        if k > 0 {
            row[k - 1] = Expr::from(1);
        }
        // c lists cₙ … c₀: cₖ is c[n − k].
        row[n - 1] = simplified(ev, -(c[n - k].clone() / lead.clone()))?;
    }
    let result = matrix(m);
    ev.step("companion", e, &result);
    Ok(Some(result))
}

/// `SylvesterMatrix[p, q, x]`.
pub(in crate::native) fn sylvester_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, q, x] = e.args() else {
        return Ok(None);
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (coefficients(ev, p, x)?, coefficients(ev, q, x)?) else {
        return Ok(None);
    };
    let (m, n) = (a.len() - 1, b.len() - 1);
    let size = m + n;
    if size == 0 {
        return Ok(None);
    }
    let shifted = |c: &[Expr], shift: usize| {
        let mut row = vec![Expr::from(0); size];
        for (j, v) in c.iter().enumerate() {
            row[shift + j] = v.clone();
        }
        row
    };
    let rows: Vec<Vec<Expr>> = (0..n)
        .map(|i| shifted(&a, i))
        .chain((0..m).map(|i| shifted(&b, i)))
        .collect();
    let result = matrix(rows);
    ev.step("sylvester", e, &result);
    Ok(Some(result))
}

/// `a` and `b` multiplied entrywise, of one shape.
fn entrywise(a: &Expr, b: &Expr) -> Option<Expr> {
    let is_list = |x: &Expr| x.head() == Some(&Expr::from(Symbol::LIST));
    match (is_list(a), is_list(b)) {
        (true, true) if a.args().len() == b.args().len() => {
            let items: Option<Vec<Expr>> = a
                .args()
                .iter()
                .zip(b.args())
                .map(|(x, y)| entrywise(x, y))
                .collect();
            Some(list(items?))
        }
        (false, false) => Some(a.clone() * b.clone()),
        _ => None,
    }
}

/// `HadamardProduct[A, B]`.
pub(in crate::native) fn hadamard_product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let is_list = |x: &Expr| x.head() == Some(&Expr::from(Symbol::LIST));
    if !is_list(a) || !is_list(b) {
        return Ok(None);
    }
    let Some(result) = entrywise(a, b) else {
        return Err(EvalError::Shape);
    };
    ev.step("hadamard", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
