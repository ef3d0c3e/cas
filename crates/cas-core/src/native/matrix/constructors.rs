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
//!   of each entry known (CAS-DOC-12).

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

#[cfg(test)]
mod tests;
