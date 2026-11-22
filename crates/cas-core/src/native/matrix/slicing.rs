//! Shapes, slices and row operations of matrices (CAS-532, CAS-533), as
//! Xcas's. Every index is 0-based, a negative one counts from the end and
//! `a..b` is `a` to `b` included, as in [`Index::of`]:
//!
//! - `RowDimension[A]` (`rowdim`, `nrows`, `rowDim`) and
//!   `ColumnDimension[A]` (`coldim`, `ncols`, `colDim`) of a matrix; a
//!   vector stays, where giac errs;
//! - `MatrixRow[A, k]` (`row`): row `k`, and a list of rows for a range
//!   (giac's is a sequence); the entry `k` of a vector. `MatrixColumn[A, k]`
//!   (`col`) likewise for the columns of a matrix;
//! - `SubMatrix[A, r1, c1, r2, c2]` (`subMat`): rows `r1` to `r2`, columns
//!   `c1` to `c2`; to the end without `r2, c2`; `A` itself alone. Reversed
//!   bounds stay;
//! - `DeleteRows[A, k]` (`delrows`) and `DeleteColumns[A, k]` (`delcols`):
//!   an index or a range deleted, `[]` when nothing is left;
//! - `MultiplyRow[x, A, n]` (`mRow`) and `ScaleRow[A, x, n]` (`scale`,
//!   `SCALE`): row `n` times `x`; `MultiplyRowAdd[x, A, n1, n2]` (`mRowAdd`)
//!   and `ScaleRowAdd[A, x, n1, n2]` (`scaleadd`, `SCALEADD`): `x` times row
//!   `n1` added to row `n2`; `RowAdd[A, n1, n2]` (`rowAdd`): row `n1` added
//!   to row `n2`; `RowSwap[A, n1, n2]` (`rowSwap`, `swaprow`, `rowswap`,
//!   `SWAPROW`) and `ColumnSwap[A, n1, n2]` (`colSwap`, `swapcol`,
//!   `colswap`, `SWAPCOL`).
//!
//! An index out of range is [`EvalError::IndexOutOfRange`]; a negative row
//! of `mRow`, where giac errs, counts from the end. Results are built
//! unevaluated, for the loop to evaluate; each call takes one step, named
//! after the giac function.

use super::{Native, Shape, list, shape, stepped};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::index::Index;
use crate::symbol::Symbol;

/// The rows and columns of a matrix.
fn dimensions(m: &Expr) -> Option<(usize, usize)> {
    match shape(m) {
        Some(Shape::Matrix(rows, columns)) => Some((rows, columns)),
        _ => None,
    }
}

/// The one row or column `index` selects; `None` for a slice or another
/// argument.
fn one(index: &Expr, len: usize) -> Result<Option<usize>, EvalError> {
    Ok(match Index::of(index, len)? {
        Some(Index::One(i)) => Some(i),
        _ => None,
    })
}

fn column(m: &Expr, j: usize) -> Expr {
    list(m.args().iter().map(|r| r.args()[j].clone()))
}

fn size(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("a size in memory"))
}

/// `RowDimension[A]`, Xcas's `rowdim`.
pub(in crate::native) fn row_dimension(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let Some((rows, _)) = dimensions(m) else {
        return Ok(None);
    };
    stepped(ev, "rowdim", e, size(rows))
}

/// `ColumnDimension[A]`, Xcas's `coldim`.
pub(in crate::native) fn column_dimension(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let Some((_, columns)) = dimensions(m) else {
        return Ok(None);
    };
    stepped(ev, "coldim", e, size(columns))
}

/// `MatrixRow[A, k]`, Xcas's `row`.
pub(in crate::native) fn matrix_row(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, k] = e.args() else {
        return Ok(None);
    };
    let result = match shape(m) {
        Some(Shape::Matrix(rows, _)) => match Index::of(k, rows)? {
            Some(Index::One(i)) => m.args()[i].clone(),
            Some(Index::Slice(r)) => list(m.args()[r].iter().cloned()),
            None => return Ok(None),
        },
        Some(Shape::Vector(n)) => match one(k, n)? {
            Some(i) => m.args()[i].clone(),
            None => return Ok(None),
        },
        None => return Ok(None),
    };
    stepped(ev, "row", e, result)
}

/// `MatrixColumn[A, k]`, Xcas's `col`.
pub(in crate::native) fn matrix_column(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, k] = e.args() else {
        return Ok(None);
    };
    let Some((_, columns)) = dimensions(m) else {
        return Ok(None);
    };
    let result = match Index::of(k, columns)? {
        Some(Index::One(j)) => column(m, j),
        Some(Index::Slice(r)) => list(r.map(|j| column(m, j))),
        None => return Ok(None),
    };
    stepped(ev, "col", e, result)
}

/// `SubMatrix[A, r1, c1, r2, c2]`, Xcas's `subMat`.
pub(in crate::native) fn sub_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, bounds)) = e.args().split_first() else {
        return Ok(None);
    };
    let Some((rows, columns)) = dimensions(m) else {
        return Ok(None);
    };
    let (r1, c1, r2, c2) = match bounds {
        [] => (0, 0, rows - 1, columns - 1),
        [r1, c1] => {
            let (Some(r1), Some(c1)) = (one(r1, rows)?, one(c1, columns)?) else {
                return Ok(None);
            };
            (r1, c1, rows - 1, columns - 1)
        }
        [r1, c1, r2, c2] => {
            let (Some(r1), Some(c1), Some(r2), Some(c2)) = (
                one(r1, rows)?,
                one(c1, columns)?,
                one(r2, rows)?,
                one(c2, columns)?,
            ) else {
                return Ok(None);
            };
            (r1, c1, r2, c2)
        }
        _ => return Ok(None),
    };
    if r2 < r1 || c2 < c1 {
        return Ok(None);
    }
    let result = list(
        m.args()[r1..=r2]
            .iter()
            .map(|r| list(r.args()[c1..=c2].iter().cloned())),
    );
    stepped(ev, "subMat", e, result)
}

/// The positions below `len` that `index` leaves; `None` when it is not an
/// index.
fn kept(index: &Expr, len: usize) -> Result<Option<Vec<usize>>, EvalError> {
    Ok(match Index::of(index, len)? {
        Some(Index::One(i)) => Some((0..len).filter(|&k| k != i).collect()),
        Some(Index::Slice(r)) => Some((0..len).filter(|k| !r.contains(k)).collect()),
        None => None,
    })
}

/// `DeleteRows[A, k]`, Xcas's `delrows`.
pub(in crate::native) fn delete_rows(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, k] = e.args() else {
        return Ok(None);
    };
    let Some((rows, _)) = dimensions(m) else {
        return Ok(None);
    };
    let Some(kept) = kept(k, rows)? else {
        return Ok(None);
    };
    let result = list(kept.into_iter().map(|i| m.args()[i].clone()));
    stepped(ev, "delrows", e, result)
}

/// `DeleteColumns[A, k]`, Xcas's `delcols`.
pub(in crate::native) fn delete_columns(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, k] = e.args() else {
        return Ok(None);
    };
    let Some((_, columns)) = dimensions(m) else {
        return Ok(None);
    };
    let Some(kept) = kept(k, columns)? else {
        return Ok(None);
    };
    let result = if kept.is_empty() {
        list(Vec::new())
    } else {
        list(
            m.args()
                .iter()
                .map(|r| list(kept.iter().map(|&j| r.args()[j].clone()))),
        )
    };
    stepped(ev, "delcols", e, result)
}

fn times(x: &Expr, y: &Expr) -> Expr {
    Expr::call(Symbol::TIMES, [x.clone(), y.clone()])
}

/// `m` with each entry `y` of row `n`, in column `j`, replaced by `f(j, y)`.
fn with_row(m: &Expr, n: usize, f: impl Fn(usize, &Expr) -> Expr) -> Expr {
    list(m.args().iter().enumerate().map(|(i, r)| {
        if i == n {
            list(r.args().iter().enumerate().map(|(j, y)| f(j, y)))
        } else {
            r.clone()
        }
    }))
}

/// Row `n` of `m` times `x`.
fn multiply(
    ev: &mut Eval<'_>,
    e: &Expr,
    name: &'static str,
    x: &Expr,
    m: &Expr,
    n: &Expr,
) -> Native {
    let Some((rows, _)) = dimensions(m) else {
        return Ok(None);
    };
    let Some(n) = one(n, rows)? else {
        return Ok(None);
    };
    stepped(ev, name, e, with_row(m, n, |_, y| times(x, y)))
}

/// `MultiplyRow[x, A, n]`, Xcas's `mRow`.
pub(in crate::native) fn multiply_row(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, m, n] = e.args() else {
        return Ok(None);
    };
    multiply(ev, e, "mRow", x, m, n)
}

/// `ScaleRow[A, x, n]`, Xcas's `scale`.
pub(in crate::native) fn scale_row(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, x, n] = e.args() else {
        return Ok(None);
    };
    multiply(ev, e, "scale", x, m, n)
}

/// Row `n1` of `m`, times `x` when there is one, added to row `n2`.
fn add(
    ev: &mut Eval<'_>,
    e: &Expr,
    name: &'static str,
    x: Option<&Expr>,
    m: &Expr,
    [n1, n2]: [&Expr; 2],
) -> Native {
    let Some((rows, _)) = dimensions(m) else {
        return Ok(None);
    };
    let (Some(n1), Some(n2)) = (one(n1, rows)?, one(n2, rows)?) else {
        return Ok(None);
    };
    let source = &m.args()[n1];
    let result = with_row(m, n2, |j, y| {
        let s = &source.args()[j];
        let added = x.map_or_else(|| s.clone(), |x| times(x, s));
        Expr::call(Symbol::PLUS, [y.clone(), added])
    });
    stepped(ev, name, e, result)
}

/// `MultiplyRowAdd[x, A, n1, n2]`, Xcas's `mRowAdd`.
pub(in crate::native) fn multiply_row_add(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, m, n1, n2] = e.args() else {
        return Ok(None);
    };
    add(ev, e, "mRowAdd", Some(x), m, [n1, n2])
}

/// `ScaleRowAdd[A, x, n1, n2]`, Xcas's `scaleadd`.
pub(in crate::native) fn scale_row_add(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, x, n1, n2] = e.args() else {
        return Ok(None);
    };
    add(ev, e, "scaleadd", Some(x), m, [n1, n2])
}

/// `RowAdd[A, n1, n2]`, Xcas's `rowAdd`.
pub(in crate::native) fn row_add(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, n1, n2] = e.args() else {
        return Ok(None);
    };
    add(ev, e, "rowAdd", None, m, [n1, n2])
}

/// `RowSwap[A, n1, n2]`, Xcas's `rowSwap`.
pub(in crate::native) fn row_swap(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, n1, n2] = e.args() else {
        return Ok(None);
    };
    let Some((rows, _)) = dimensions(m) else {
        return Ok(None);
    };
    let (Some(n1), Some(n2)) = (one(n1, rows)?, one(n2, rows)?) else {
        return Ok(None);
    };
    let mut swapped = m.args().to_vec();
    swapped.swap(n1, n2);
    stepped(ev, "rowSwap", e, list(swapped))
}

/// `ColumnSwap[A, n1, n2]`, Xcas's `colSwap`.
pub(in crate::native) fn column_swap(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, n1, n2] = e.args() else {
        return Ok(None);
    };
    let Some((_, columns)) = dimensions(m) else {
        return Ok(None);
    };
    let (Some(n1), Some(n2)) = (one(n1, columns)?, one(n2, columns)?) else {
        return Ok(None);
    };
    let result = list(m.args().iter().map(|r| {
        let mut r = r.args().to_vec();
        r.swap(n1, n2);
        list(r)
    }));
    stepped(ev, "colSwap", e, result)
}
