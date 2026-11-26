//! giac's reshaping, bordering, triangles, quadratic forms and changes of
//! basis (CAS-639), as Xcas's; where giac errs, the call stays with a
//! warning:
//!
//! - `Redimension[A, r, c]` and `Redimension[A, [r, c]]` (`REDIM`,
//!   `redim`): each row of a matrix cut or padded with 0 to `c` entries,
//!   rows of 0 to `r` rows; a vector filled into `r × c` row by row.
//!   `Redimension[A, n]` and `Redimension[A, [n]]`: the entries row by row
//!   cut or padded with 0 to `n`, `[]` for `n ≤ 0`;
//! - `Border[A, v]` (`border`): `v`'s entries, or its rows, appended to the
//!   rows of `A`;
//! - `SemiAugment[a, b]` (`semi_augment`): the entries of `a` then of `b`,
//!   one not a list counting as one entry;
//! - `Simultaneous[A, B]` (`simult`): `A⁻¹ B` for a square `A` and a matrix
//!   `B` of as many rows;
//! - `NewMatrix[r, c]` (`newMat`): the `r × c` matrix of 0, `[]` for
//!   `r = 0`; giac takes `|r|` and `|c|`, here they stay;
//! - `UpperTriangle[A]` (`upper`), `LowerTriangle[A]` (`lower`): the entries
//!   on and above, or below, the diagonal of a matrix, the others 0; a
//!   string in ASCII capitals, or small letters;
//! - `MatrixToQuadratic[A, [x…]]` (`a2q`): `Σ aᵢⱼ xᵢ xⱼ`, expanded;
//!   `QuadraticToMatrix[q, [x…]]` (`q2a`): the symmetric matrix of a
//!   quadratic form, half its Hessian;
//! - `ChangeBasis[A, P]` (`changebase`): `P⁻¹ A P`, expanded;
//! - `ConjugateTranspose[A]` (`trn`): the conjugate of the transpose;
//! - `Basis[A]` (`basis`): the rows of `RowSpace[A]` as a set, as giac's;
//! - `FillValue[v, a]`, internal: the value giac's `Fill(v, A)` stores in
//!   `A`, a matrix of `a`'s size or a list of `a`'s length with every entry
//!   `v`, else `v`.
//!
//! Sizes beyond `Limits::terms` stay. Results are built unevaluated, for
//! the loop to evaluate; each call takes one step, named after the giac
//! function.

use super::{Native, Shape, is_list, list, matrix_argument, not_square, reduced, shape, stepped};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::linear::dot;
use crate::symbol::Symbol;
use crate::trace::Warning;

/// Says what is wrong with `e`; the call stays.
fn warned(ev: &mut Eval<'_>, e: &Expr, message: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// The entries of a vector, or of a matrix row by row.
fn flattened(m: &Expr) -> Vec<Expr> {
    match shape(m) {
        Some(Shape::Matrix(..)) => m
            .args()
            .iter()
            .flat_map(|r| r.args().iter().cloned())
            .collect(),
        _ => m.args().to_vec(),
    }
}

/// `items` cut or padded with 0 to `n`.
fn padded(items: impl IntoIterator<Item = Expr>, n: usize) -> Expr {
    list(
        items
            .into_iter()
            .chain(std::iter::repeat_with(|| Expr::from(0)))
            .take(n),
    )
}

/// The integer dimensions `dims`; `Err(true)` when one is a number that is
/// not an integer, `Err(false)` for another.
fn integers(dims: &[Expr]) -> Result<Vec<i64>, bool> {
    dims.iter()
        .map(Expr::as_i64)
        .collect::<Option<Vec<i64>>>()
        .ok_or_else(|| dims.iter().any(Expr::is_number))
}

/// `Redimension[A, r, c]`, `Redimension[A, n]`, a list of dimensions too.
pub(in crate::native) fn redimension(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (m, dims) = match e.args() {
        [m, d] if is_list(d) => (m, d.args()),
        [m, rest @ ..] if (1..=2).contains(&rest.len()) => (m, rest),
        _ => return Ok(None),
    };
    let Some(sh) = shape(m) else {
        if m.is_number() {
            return warned(ev, e, "a vector or a matrix is needed");
        }
        return Ok(None);
    };
    let dims = match integers(dims) {
        Ok(dims) => dims,
        Err(true) => return warned(ev, e, "integer dimensions are needed"),
        Err(false) => return Ok(None),
    };
    let terms = ev.context().limits().terms;
    let result = match dims[..] {
        [n] => {
            let n = usize::try_from(n).unwrap_or(0);
            if n > terms {
                return Ok(None);
            }
            padded(flattened(m), n)
        }
        [r, c] => {
            let (Ok(r @ 1..), Ok(c @ 1..)) = (usize::try_from(r), usize::try_from(c)) else {
                return warned(ev, e, "dimensions of at least 1 are needed");
            };
            if r.saturating_mul(c) > terms {
                return Ok(None);
            }
            match sh {
                Shape::Matrix(..) => list((0..r).map(|i| match m.args().get(i) {
                    Some(row) => padded(row.args().iter().cloned(), c),
                    None => padded([], c),
                })),
                Shape::Vector(_) => {
                    let mut entries = flattened(m).into_iter();
                    list((0..r).map(|_| padded(entries.by_ref().take(c), c)))
                }
            }
        }
        _ => return Ok(None),
    };
    stepped(ev, "REDIM", e, result)
}

/// `Border[A, v]`, Xcas's `border`.
pub(in crate::native) fn border(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, v] = e.args() else {
        return Ok(None);
    };
    let Some(Shape::Matrix(rows, _)) = shape(m) else {
        if is_list(m) || m.is_number() {
            return warned(ev, e, "a matrix is needed");
        }
        return Ok(None);
    };
    let nested = |x: &Expr| is_list(x) && x.args().iter().any(is_list);
    if !is_list(v) || v.args().len() != rows || v.args().iter().any(nested) {
        return warned(ev, e, "a vector or a matrix of as many rows is needed");
    }
    let result = list(m.args().iter().zip(v.args()).map(|(row, x)| {
        let added = if is_list(x) {
            x.args().to_vec()
        } else {
            vec![x.clone()]
        };
        list(row.args().iter().cloned().chain(added))
    }));
    stepped(ev, "border", e, result)
}

/// `SemiAugment[a, b]`, Xcas's `semi_augment`.
pub(in crate::native) fn semi_augment(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if !is_list(a) && !is_list(b) {
        return Ok(None);
    }
    let entries = |x: &Expr| {
        if is_list(x) {
            x.args().to_vec()
        } else {
            vec![x.clone()]
        }
    };
    let result = list(entries(a).into_iter().chain(entries(b)));
    stepped(ev, "semi_augment", e, result)
}

/// `Inverse[p]` evaluated, `None` when it stays (a warning given there).
fn inverse(ev: &mut Eval<'_>, p: &Expr) -> Result<Option<Expr>, crate::context::EvalError> {
    let inv = ev.eval(&Expr::call(Symbol::INVERSE, [p.clone()]))?;
    Ok(matches!(shape(&inv), Some(Shape::Matrix(..))).then_some(inv))
}

/// The size of `m` when it is a square matrix; `Err(())` for a matrix
/// that is not, `Ok(None)` for another.
fn square_size(m: &Expr) -> Result<Option<usize>, ()> {
    match shape(m) {
        Some(Shape::Matrix(r, c)) if r == c => Ok(Some(r)),
        Some(Shape::Matrix(..)) => Err(()),
        _ => Ok(None),
    }
}

/// `Simultaneous[A, B]`, Xcas's `simult`: `A⁻¹ B`.
pub(in crate::native) fn simultaneous(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let Ok(n) = square_size(a) else {
        not_square(ev, e);
        return Ok(None);
    };
    let Some(n) = n else {
        return Ok(None);
    };
    match shape(b) {
        Some(Shape::Matrix(rows, _)) if rows == n => {}
        Some(_) => return warned(ev, e, "a matrix of as many rows is needed"),
        None => return Ok(None),
    }
    let Some(inv) = inverse(ev, a)? else {
        return Ok(None);
    };
    let result = dot(&inv, b)?;
    stepped(ev, "simult", e, result)
}

/// `ChangeBasis[A, P]`, Xcas's `changebase`: `P⁻¹ A P`, expanded.
pub(in crate::native) fn change_basis(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, p] = e.args() else {
        return Ok(None);
    };
    let (Ok(n), Ok(k)) = (square_size(a), square_size(p)) else {
        not_square(ev, e);
        return Ok(None);
    };
    let (Some(n), Some(k)) = (n, k) else {
        return Ok(None);
    };
    if n != k {
        return warned(ev, e, "matrices of one size are needed");
    }
    let Some(inv) = inverse(ev, p)? else {
        return Ok(None);
    };
    let left = ev.eval(&dot(&inv, a)?)?;
    let result = Expr::call(Symbol::EXPAND, [dot(&left, p)?]);
    stepped(ev, "changebase", e, result)
}

/// `NewMatrix[r, c]`, Xcas's `newMat`.
pub(in crate::native) fn new_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if e.args().len() != 2 {
        return Ok(None);
    }
    let dims = match integers(e.args()) {
        Ok(dims) => dims,
        Err(true) => return warned(ev, e, "integer dimensions are needed"),
        Err(false) => return Ok(None),
    };
    let [r, c] = dims[..] else {
        return Ok(None);
    };
    if r == 0 {
        return stepped(ev, "newMat", e, list([]));
    }
    let (Ok(r @ 1..), Ok(c @ 1..)) = (usize::try_from(r), usize::try_from(c)) else {
        return warned(ev, e, "dimensions of at least 1 are needed");
    };
    if r.saturating_mul(c) > ev.context().limits().terms {
        return Ok(None);
    }
    let result = list((0..r).map(|_| padded([], c)));
    stepped(ev, "newMat", e, result)
}

/// `UpperTriangle[A]` and `LowerTriangle[A]`, Xcas's `upper` and `lower`.
fn triangle(ev: &mut Eval<'_>, e: &Expr, upper: bool) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let name = if upper { "upper" } else { "lower" };
    if let Some(s) = m.as_str() {
        let changed = if upper {
            s.to_ascii_uppercase()
        } else {
            s.to_ascii_lowercase()
        };
        return stepped(ev, name, e, Expr::string(&changed));
    }
    let Some(Shape::Matrix(..)) = shape(m) else {
        if is_list(m) || m.is_number() {
            return warned(ev, e, "a matrix or a string is needed");
        }
        return Ok(None);
    };
    let result = list(m.args().iter().enumerate().map(|(i, row)| {
        list(row.args().iter().enumerate().map(|(j, x)| {
            if (upper && j < i) || (!upper && j > i) {
                Expr::from(0)
            } else {
                x.clone()
            }
        }))
    }));
    stepped(ev, name, e, result)
}

/// `UpperTriangle[A]`, Xcas's `upper`.
pub(in crate::native) fn upper_triangle(ev: &mut Eval<'_>, e: &Expr) -> Native {
    triangle(ev, e, true)
}

/// `LowerTriangle[A]`, Xcas's `lower`.
pub(in crate::native) fn lower_triangle(ev: &mut Eval<'_>, e: &Expr) -> Native {
    triangle(ev, e, false)
}

/// `MatrixToQuadratic[A, [x…]]`, Xcas's `a2q`: `Σ aᵢⱼ xᵢ xⱼ`, expanded.
pub(in crate::native) fn matrix_to_quadratic(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, vars] = e.args() else {
        return Ok(None);
    };
    if !is_list(vars) {
        return Ok(None);
    }
    let x = vars.args();
    match shape(a) {
        Some(Shape::Matrix(r, c)) if r == c && r == x.len() => {}
        Some(_) => return warned(ev, e, "a square matrix of the variables' size is needed"),
        None => return Ok(None),
    }
    let terms: Vec<Expr> = a
        .args()
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.args().iter().enumerate().map(move |(j, aij)| {
                Expr::call(Symbol::TIMES, [aij.clone(), x[i].clone(), x[j].clone()])
            })
        })
        .collect();
    let result = Expr::call(Symbol::EXPAND, [Expr::call(Symbol::PLUS, terms)]);
    stepped(ev, "a2q", e, result)
}

/// `QuadraticToMatrix[q, [x…]]`, Xcas's `q2a`: half the Hessian of a
/// quadratic form in the variables.
pub(in crate::native) fn quadratic_to_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [q, vars] = e.args() else {
        return Ok(None);
    };
    if !is_list(vars)
        || vars.args().is_empty()
        || vars.args().iter().any(|v| v.as_symbol().is_none())
    {
        return Ok(None);
    }
    let x = vars.args();
    let n = x.len();
    if n.saturating_mul(n) > ev.context().limits().terms {
        return Ok(None);
    }
    let mut rows = Vec::with_capacity(n);
    for xi in x {
        let di = ev.eval(&Expr::call(Symbol::D, [q.clone(), xi.clone()]))?;
        let mut row = Vec::with_capacity(n);
        for xj in x {
            let dij = Expr::call(Symbol::D, [di.clone(), xj.clone()]);
            let h = ev.eval(&(dij / Expr::from(2)))?;
            if x.iter().any(|v| h.contains(v)) {
                return warned(ev, e, "a quadratic form in the variables is needed");
            }
            row.push(h);
        }
        rows.push(row);
    }
    let form: Vec<Expr> = rows
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.iter().enumerate().map(move |(j, h)| {
                Expr::call(Symbol::TIMES, [h.clone(), x[i].clone(), x[j].clone()])
            })
        })
        .collect();
    let residual = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [q.clone() - Expr::call(Symbol::PLUS, form)],
    ))?;
    if residual != Expr::from(0) {
        return warned(ev, e, "a quadratic form in the variables is needed");
    }
    stepped(ev, "q2a", e, list(rows.into_iter().map(list)))
}

/// `ConjugateTranspose[A]`, Xcas's `trn`.
pub(in crate::native) fn conjugate_transpose(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    match shape(m) {
        Some(Shape::Matrix(..)) => {}
        Some(Shape::Vector(n)) if n > 0 => {}
        _ => return Ok(None),
    }
    let result = Expr::call(
        Symbol::CONJUGATE,
        [Expr::call(Symbol::TRANSPOSE, [m.clone()])],
    );
    stepped(ev, "trn", e, result)
}

/// `Basis[A]`, Xcas's `basis`: the rows of `RowSpace[A]`, a set.
pub(in crate::native) fn basis(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(m) = matrix_argument(e) else {
        return Ok(None);
    };
    let Some(r) = reduced(ev, e, m)? else {
        return Ok(None);
    };
    let rank = r.pivots.len();
    let result = Expr::call(
        Symbol::FINITE_SET,
        r.rows.into_iter().take(rank).map(list).collect::<Vec<_>>(),
    );
    stepped(ev, "basis", e, result)
}

/// `FillValue[v, a]`, what giac's `Fill(v, A)` stores: a matrix or list of
/// `a`'s shape filled with `v`, else `v`.
pub(in crate::native) fn fill_value(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v, a] = e.args() else {
        return Ok(None);
    };
    let result = match shape(a) {
        Some(Shape::Matrix(r, c)) => {
            if r.saturating_mul(c) > ev.context().limits().terms {
                return Ok(None);
            }
            list((0..r).map(|_| list(std::iter::repeat_n(v.clone(), c))))
        }
        _ if is_list(a) => list(std::iter::repeat_n(v.clone(), a.args().len())),
        _ => v.clone(),
    };
    stepped(ev, "Fill", e, result)
}
