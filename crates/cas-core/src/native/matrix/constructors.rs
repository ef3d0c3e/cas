//! Matrix constructors and norms (CAS-DOC-30 §4, CAS-347):
//!
//! - `DiagonalMatrix[l]` (`diag`): the diagonal matrix of a list, and the
//!   diagonal of a matrix, as giac's `diag([[1,2],[3,4]])` is `[1,4]`;
//! - `HilbertMatrix[n]` (`hilbert`): `[1/(i + j − 1)]`, `1 ≤ n` with `n²`
//!   within `Limits::terms`;
//! - `VandermondeMatrix[l]` (`vandermonde`): the rows `[1, l_i, l_i², …]`,
//!   as giac's;
//! - `Norm[v]` (`l2norm`): `√(Σ |v_i|²)`, and `FrobeniusNorm[A]`
//!   (`frobenius_norm`) over every entry, `|z|² = z z̄` with the conjugate
//!   of each entry known (CAS-DOC-12);
//! - `Normalize[v]` (`normalize`, `unitV`, CAS-564): `v/Norm[v]` for a
//!   vector, `z/Abs[z]` otherwise, as giac's `normalize(x)` is `x/abs(x)`;
//!   the call stays where the norm does or is 0;
//!
//! and giac's constructors (CAS-531), for `1 ≤ p, q` with `p q` entries
//! within `Limits::terms`, else the call stays:
//!
//! - `Matrix[p, q, a]` (`matrix`) and `MakeMatrix[a, p, q]` (`makemat`):
//!   zeros without `a`; `a[j, k]`, 0-based, for a function; a vector's
//!   entries row by row, padded with 0, extras dropped; `a` everywhere
//!   otherwise;
//! - `ListToMatrix[l, n]` (`list2mat`): `l` in rows of `n ≥ 1`, padded with
//!   0; `MatrixToList[A]` (`mat2list`): the rows of `A` one after the
//!   other, a vector itself;
//! - `BlockMatrix[n, m, l]` (`blockmatrix`): the `n m` matrices of `l` as
//!   `n` block rows of `m`, whose shapes must fit;
//! - `JordanBlock[a, n]` (`JordanBlock`): `a` on the diagonal, 1 above it;
//! - `BlockDiagonalMatrix[l]` (`BlockDiagonal`): the square matrices of `l`
//!   along the diagonal; as `diag` for a matrix or a vector.

use super::conjugate;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// Whether `e` is a list of lists of one length.
fn is_matrix(e: &Expr) -> bool {
    is_list(e)
        && !e.args().is_empty()
        && e.args()
            .iter()
            .all(|r| is_list(r) && r.args().len() == e.args()[0].args().len())
}

pub(in crate::native) fn diagonal_matrix(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if is_matrix(l) {
        let n = l.args().len().min(l.args()[0].args().len());
        return Ok(Some(list(
            (0..n).map(|i| l.args()[i].args()[i].clone()).collect(),
        )));
    }
    if !is_list(l) || l.args().is_empty() {
        return Ok(None);
    }
    let n = l.args().len();
    let rows = (0..n)
        .map(|i| {
            list(
                (0..n)
                    .map(|j| {
                        if i == j {
                            l.args()[i].clone()
                        } else {
                            Expr::from(0)
                        }
                    })
                    .collect(),
            )
        })
        .collect();
    Ok(Some(list(rows)))
}

/// A size `1 ≤ n` whose `n²` entries are within `Limits::terms`.
fn size(ev: &Eval<'_>, n: &Expr) -> Option<usize> {
    let n = usize::try_from(n.as_i64()?).ok()?;
    (n >= 1 && n.checked_mul(n)? <= ev.context().limits().terms).then_some(n)
}

pub(in crate::native) fn hilbert_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    // A list of numbers: giac's analytic signal (CAS-635).
    if n.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
        return crate::native::analytic::hilbert(ev, e, n);
    }
    let Some(n) = size(ev, n) else {
        return Ok(None);
    };
    let entry = |i: usize, j: usize| Expr::from(i64::try_from(i + j + 1).expect("small")).pow(-1);
    let rows = (0..n)
        .map(|i| list((0..n).map(|j| entry(i, j)).collect()))
        .collect();
    ev.eval(&list(rows)).map(Some)
}

pub(in crate::native) fn vandermonde_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) || l.args().is_empty() {
        return Ok(None);
    }
    let n = l.args().len();
    if size(ev, &Expr::from(i64::try_from(n).expect("small"))).is_none() {
        return Ok(None);
    }
    let rows = l
        .args()
        .iter()
        .map(|x| {
            list(
                (0..n)
                    .map(|j| x.pow(i64::try_from(j).expect("small")))
                    .collect(),
            )
        })
        .collect();
    ev.eval(&list(rows)).map(Some)
}

/// `√(Σ x x̄)` over `entries`; `None` when a conjugate is not known.
fn root_of_squares(ev: &mut Eval<'_>, entries: &[Expr]) -> Native {
    let mut terms = Vec::with_capacity(entries.len());
    for x in entries {
        let Some(c) = conjugate(ev.context(), x) else {
            return Ok(None);
        };
        terms.push(x.clone() * c);
    }
    ev.eval(&Expr::call(Symbol::SQRT, [Expr::call(Symbol::PLUS, terms)]))
        .map(Some)
}

pub(in crate::native) fn norm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v] = e.args() else {
        return Ok(None);
    };
    if !is_list(v) || v.args().iter().any(is_list) {
        return Ok(None);
    }
    root_of_squares(ev, v.args())
}

pub(in crate::native) fn normalize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v] = e.args() else {
        return Ok(None);
    };
    let size = if is_list(v) {
        match norm(ev, &Expr::call(Symbol::NORM, [v.clone()]))? {
            Some(n) => n,
            None => return Ok(None),
        }
    } else {
        ev.eval(&Expr::call(Symbol::ABS, [v.clone()]))?
    };
    if size == Expr::from(0) {
        return Ok(None);
    }
    Ok(Some(v / size))
}

pub(in crate::native) fn frobenius_norm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    if !is_matrix(a) {
        return Ok(None);
    }
    let entries: Vec<Expr> = a
        .args()
        .iter()
        .flat_map(|r| r.args().iter().cloned())
        .collect();
    root_of_squares(ev, &entries)
}

fn index(i: usize) -> Expr {
    Expr::from(i64::try_from(i).expect("a size in memory"))
}

/// An integer `1 ≤ n`.
fn positive(n: &Expr) -> Option<usize> {
    usize::try_from(n.as_i64()?).ok().filter(|&n| n >= 1)
}

/// Whether `n` entries are within `Limits::terms`.
fn within(ev: &Eval<'_>, n: Option<usize>) -> bool {
    n.is_some_and(|n| n <= ev.context().limits().terms)
}

/// The `p × q` matrix of `a`: zeros without it, `a[j, k]` for a function,
/// a vector's entries row by row padded with 0, or `a` everywhere.
fn filled(ev: &mut Eval<'_>, p: &Expr, q: &Expr, a: Option<&Expr>) -> Native {
    let (Some(p), Some(q)) = (positive(p), positive(q)) else {
        return Ok(None);
    };
    if !within(ev, p.checked_mul(q)) {
        return Ok(None);
    }
    let zero = Expr::from(0);
    let a = a.unwrap_or(&zero);
    if a.head() == Some(&Expr::from(Symbol::FUNCTION)) {
        let mut rows = Vec::with_capacity(p);
        for j in 0..p {
            let mut row = Vec::with_capacity(q);
            for k in 0..q {
                row.push(ev.eval(&Expr::apply(a, &[index(j), index(k)]))?);
            }
            rows.push(list(row));
        }
        return Ok(Some(list(rows)));
    }
    if is_list(a) && a.args().iter().any(is_list) {
        return Ok(None);
    }
    let entry = |i: usize| {
        if is_list(a) {
            a.args().get(i).cloned().unwrap_or_else(|| zero.clone())
        } else {
            a.clone()
        }
    };
    let rows = (0..p)
        .map(|j| list((0..q).map(|k| entry(j * q + k)).collect()))
        .collect();
    Ok(Some(list(rows)))
}

pub(in crate::native) fn matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [p, q] => filled(ev, p, q, None),
        [p, q, a] => filled(ev, p, q, Some(a)),
        _ => Ok(None),
    }
}

pub(in crate::native) fn make_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, p, q] = e.args() else {
        return Ok(None);
    };
    filled(ev, p, q, Some(a))
}

pub(in crate::native) fn list_to_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = positive(n) else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let items = l.args();
    let rows = items.len().div_ceil(n);
    if !within(ev, rows.checked_mul(n)) {
        return Ok(None);
    }
    let rows = (0..rows)
        .map(|i| {
            list(
                (0..n)
                    .map(|k| {
                        items
                            .get(i * n + k)
                            .cloned()
                            .unwrap_or_else(|| Expr::from(0))
                    })
                    .collect(),
            )
        })
        .collect();
    Ok(Some(list(rows)))
}

pub(in crate::native) fn matrix_to_list(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    if is_matrix(m) {
        return Ok(Some(list(
            m.args()
                .iter()
                .flat_map(|r| r.args().iter().cloned())
                .collect(),
        )));
    }
    Ok((is_list(m) && !m.args().iter().any(is_list)).then(|| m.clone()))
}

pub(in crate::native) fn block_matrix(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, m, l] = e.args() else {
        return Ok(None);
    };
    let (Some(n), Some(m)) = (positive(n), positive(m)) else {
        return Ok(None);
    };
    if !is_list(l) || Some(l.args().len()) != n.checked_mul(m) || !l.args().iter().all(is_matrix) {
        return Ok(None);
    }
    let mut rows = Vec::new();
    let mut width = None;
    for blocks in l.args().chunks(m) {
        let height = blocks[0].args().len();
        if blocks.iter().any(|b| b.args().len() != height) {
            return Ok(None);
        }
        let w: usize = blocks.iter().map(|b| b.args()[0].args().len()).sum();
        if *width.get_or_insert(w) != w {
            return Ok(None);
        }
        for r in 0..height {
            rows.push(list(
                blocks
                    .iter()
                    .flat_map(|b| b.args()[r].args().iter().cloned())
                    .collect(),
            ));
        }
    }
    Ok(Some(list(rows)))
}

pub(in crate::native) fn jordan_block(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = size(ev, n) else {
        return Ok(None);
    };
    let entry = |i: usize, j: usize| {
        if i == j {
            a.clone()
        } else {
            Expr::from(i64::from(j == i + 1))
        }
    };
    let rows = (0..n)
        .map(|i| list((0..n).map(|j| entry(i, j)).collect()))
        .collect();
    Ok(Some(list(rows)))
}

pub(in crate::native) fn block_diagonal_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) || l.args().is_empty() {
        return Ok(None);
    }
    let blocks = l.args();
    if blocks.iter().all(is_matrix) {
        if blocks
            .iter()
            .any(|b| b.args().len() != b.args()[0].args().len())
        {
            return Ok(None);
        }
        let total: usize = blocks.iter().map(|b| b.args().len()).sum();
        if !within(ev, total.checked_mul(total)) {
            return Ok(None);
        }
        let mut rows = Vec::with_capacity(total);
        let mut at = 0;
        for b in blocks {
            let k = b.args().len();
            for r in b.args() {
                let mut row = vec![Expr::from(0); total];
                row[at..at + k].clone_from_slice(r.args());
                rows.push(list(row));
            }
            at += k;
        }
        return Ok(Some(list(rows)));
    }
    if is_matrix(l) || !blocks.iter().any(is_list) {
        return diagonal_matrix(ev, e);
    }
    Ok(None)
}

#[cfg(test)]
mod tests;
