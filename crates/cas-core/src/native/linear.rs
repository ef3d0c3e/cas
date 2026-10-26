//! Products of lists (CAS-DOC-4 §6.2, R13, CAS-62): `*` is the product of
//! linear algebra, as in Xcas, and `.*` the element-wise one, with `./` and
//! `.^` next to it (CAS-123).
//!
//! A vector is a list whose entries are not lists, and a matrix a non-empty
//! list of vectors of one length, not zero. In a product, a scalar scales
//! every entry; two vectors of one length give their dot product; a matrix
//! times a vector, a vector times a matrix (as a row) and two matrices give
//! the matrix product. Other dimensions are [`EvalError::Shape`]. Several
//! lists fold from the left, `u*v*w` being `(u·v) w`, which is not
//! associative, as in Xcas. Entries are built unevaluated, for the loop to
//! evaluate.
//!
//! A square matrix to an integer power is the matrix power, by squaring,
//! each product evaluated; a negative power is that of the inverse
//! (CAS-DOC-12, CAS-219).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// Powers above this are left unevaluated: their entries could not be held.
const MAX_MATRIX_POWER: u64 = 1 << 16;

/// The dimensions of a vector or a matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Shape {
    Vector(usize),
    Matrix(usize, usize),
}

pub(super) fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// The shape of `e` if it is a vector or a matrix.
pub(super) fn shape(e: &Expr) -> Option<Shape> {
    if !is_list(e) {
        return None;
    }
    let rows = e.args();
    if !rows.iter().any(is_list) {
        return Some(Shape::Vector(rows.len()));
    }
    let n = rows.first()?.args().len();
    (n > 0
        && rows
            .iter()
            .all(|r| is_list(r) && r.args().len() == n && !r.args().iter().any(is_list)))
    .then_some(Shape::Matrix(rows.len(), n))
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    let items: Vec<Expr> = items.into_iter().collect();
    Expr::apply(&Symbol::LIST.into(), &items)
}

fn times(a: &Expr, b: &Expr) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])
}

/// The sum of `terms`: 0 for none, the term for one.
fn sum(terms: impl IntoIterator<Item = Expr>) -> Expr {
    let mut terms: Vec<Expr> = terms.into_iter().collect();
    match terms.len() {
        0 => Expr::from(0),
        1 => terms.pop().expect("one term"),
        _ => Expr::apply(&Symbol::PLUS.into(), &terms),
    }
}

/// The product of two vectors or matrices, unevaluated.
///
/// # Errors
///
/// [`EvalError::Shape`] when their dimensions do not fit.
#[allow(
    clippy::many_single_char_names,
    reason = "dimensions and indices, as written in mathematics"
)]
pub(super) fn dot(a: &Expr, b: &Expr) -> Result<Expr, EvalError> {
    let (Some(sa), Some(sb)) = (shape(a), shape(b)) else {
        return Err(EvalError::Shape);
    };
    let at = |m: &Expr, i: usize, j: usize| m.args()[i].args()[j].clone();
    Ok(match (sa, sb) {
        (Shape::Vector(n), Shape::Vector(k)) if n == k => {
            sum((0..n).map(|i| times(&a.args()[i], &b.args()[i])))
        }
        (Shape::Matrix(m, n), Shape::Vector(k)) if n == k => {
            list((0..m).map(|i| sum((0..n).map(|j| times(&at(a, i, j), &b.args()[j])))))
        }
        (Shape::Vector(k), Shape::Matrix(m, n)) if k == m => {
            list((0..n).map(|j| sum((0..m).map(|i| times(&a.args()[i], &at(b, i, j))))))
        }
        (Shape::Matrix(m, n), Shape::Matrix(k, p)) if n == k => {
            list((0..m).map(|i| {
                list((0..p).map(|j| sum((0..n).map(|l| times(&at(a, i, l), &at(b, l, j))))))
            }))
        }
        _ => return Err(EvalError::Shape),
    })
}

/// `c` times every entry of the list `l`, at every depth.
pub(super) fn scale(c: &Expr, l: &Expr) -> Expr {
    list(
        l.args()
            .iter()
            .map(|x| if is_list(x) { scale(c, x) } else { times(c, x) }),
    )
}

/// The product of `lists` from the left, or `None` when one of them is
/// neither a vector nor a matrix. After a dot product of vectors, the
/// scalar scales what follows.
///
/// # Errors
///
/// [`EvalError::Shape`] when dimensions do not fit.
pub(super) fn fold(lists: &[Expr]) -> Result<Option<Expr>, EvalError> {
    let Some((first, rest)) = lists.split_first() else {
        return Ok(None);
    };
    if lists.iter().any(|l| shape(l).is_none()) {
        return Ok(None);
    }
    let mut product = first.clone();
    for l in rest {
        product = if shape(&product).is_some() {
            dot(&product, l)?
        } else {
            scale(&product, l)
        };
    }
    Ok(Some(product))
}

/// The `n × n` identity matrix.
fn identity(n: usize) -> Expr {
    list((0..n).map(|i| list((0..n).map(|j| Expr::from(i64::from(i == j))))))
}

/// `m^k` for a matrix `m` and an integer `k`: the matrix power for a square
/// `m` and `|k| <= 2^16`, by squaring, each product evaluated, of the
/// inverse for a negative `k`; `None` for a vector, a larger `k`, or a
/// matrix without an inverse (after a warning).
///
/// # Errors
///
/// [`EvalError::Shape`] for a matrix that is not square, and what
/// evaluating the products gives.
pub(super) fn power(ev: &mut Eval<'_>, m: &Expr, k: &Expr) -> Native {
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    let Some(k) = k.as_i64() else {
        return Ok(None);
    };
    if rows != columns {
        return Err(EvalError::Shape);
    }
    if k.unsigned_abs() > MAX_MATRIX_POWER {
        return Ok(None);
    }
    // A negative power is a power of the inverse (CAS-219).
    let mut base = if k < 0 {
        match super::matrix::inverse_of(ev, m, m)? {
            Some(inverse) => inverse,
            None => return Ok(None),
        }
    } else {
        m.clone()
    };
    let mut k = k.unsigned_abs();
    let mut result: Option<Expr> = None;
    while k > 0 {
        if k & 1 == 1 {
            result = Some(match result {
                None => base.clone(),
                Some(r) => ev.eval(&dot(&r, &base)?)?,
            });
        }
        k >>= 1;
        if k > 0 {
            base = ev.eval(&dot(&base, &base)?)?;
        }
    }
    Ok(Some(result.unwrap_or_else(|| identity(rows))))
}

/// `ElementwiseTimes` (`.*`): `Listable`, so once threaded its arguments
/// are not lists, and it is their product.
pub(super) fn elementwise_times(_: &mut Eval<'_>, e: &Expr) -> Native {
    if e.args().iter().any(is_list) {
        return Ok(None);
    }
    Ok(Some(Expr::apply(&Symbol::TIMES.into(), e.args())))
}

/// `ElementwiseDivide` (`./`): `Listable`, so once threaded its arguments
/// are not lists, and it is their quotient.
pub(super) fn elementwise_divide(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if is_list(a) || is_list(b) {
        return Ok(None);
    }
    let reciprocal = Expr::apply(&Symbol::POWER.into(), &[b.clone(), Expr::from(-1)]);
    Ok(Some(times(a, &reciprocal)))
}

/// `ElementwisePower` (`.^`): `Listable`, so once threaded its arguments
/// are not lists, and it is their power.
pub(super) fn elementwise_power(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if is_list(a) || is_list(b) {
        return Ok(None);
    }
    Ok(Some(Expr::apply(&Symbol::POWER.into(), e.args())))
}

#[cfg(test)]
mod tests;
