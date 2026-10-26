//! Linear algebra (CAS-DOC-12): matrices as lists of lists, a vector a list
//! whose entries are not lists ([`linear::shape`]), and the products of
//! CAS-DOC-4 §6.2 in [`linear`].
//!
//! - `Transpose[m]`, Xcas's `tran`: the transpose, a vector as a column
//!   matrix.
//! - `IdentityMatrix[n]`, Xcas's `idn`: the `n × n` identity, for
//!   `1 ≤ n` with `n²` within `Limits::terms`.
//! - `Tr[m]`, Xcas's `trace`: the sum of the diagonal of a square matrix;
//!   another stays with a warning, as giac's stays.
//! - `Dot[u, v]`, Xcas's `dot`: `Σ conj(uₖ) vₖ` of two vectors of one
//!   length, and entrywise of two matrices of one shape, the conjugate of
//!   the first as giac's (`dot([i,0],[1,0])` is `-i`). An entry is its own
//!   conjugate when it is known real; another leaves the call.
//! - `Cross[u, v]`, Xcas's `cross`: the cross product of 3-vectors, a
//!   2-vector taken with a third entry 0 (`cross([1,2],[3,4])` is
//!   `[0,0,-2]`).
//!
//! - `Det[m]`, Xcas's `det`, and `Inverse[m]`, Xcas's `inv`, of a square
//!   matrix, by [`elimination`]: Bareiss's for the determinant, Gauss–Jordan
//!   on `[m | I]` for the inverse, each entry in lowest terms; negative
//!   powers are powers of the inverse (`linear::power`). A singular matrix,
//!   or one whose numeric pivots cannot be told from 0, stays with a
//!   warning, where giac errs.
//!
//! Dimensions that do not fit are [`EvalError::Shape`], as in products.
//! Results are built unevaluated, for the loop to evaluate; each call takes
//! one step, named after the giac function.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear algebra: vectors u and v, an m × n matrix"
)]

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_num::{Complex, Integer};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::linear::{Shape, is_list, shape};
use super::polynomial::{Numbers, cancelled, convert};
use crate::assume::{Set, Truth, element};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;
use elimination::{
    Entries, Numeric, bareiss, gauss_jordan, numeric, numeric_det, numeric_gauss_jordan,
    scaled_rows,
};

type Native = Result<Option<Expr>, EvalError>;

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn entry(m: &Expr, i: usize, j: usize) -> &Expr {
    &m.args()[i].args()[j]
}

/// A step and its result.
fn stepped(ev: &mut Eval<'_>, name: &'static str, e: &Expr, result: Expr) -> Native {
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// Says that a square matrix is needed.
fn not_square(ev: &mut Eval<'_>, e: &Expr) {
    ev.warn(Warning::Native {
        message: "a square matrix is needed".into(),
        expr: e.clone(),
    });
}

/// `Transpose[m]`.
pub(super) fn transpose(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let result = match shape(m) {
        Some(Shape::Matrix(rows, columns)) => {
            list((0..columns).map(|j| list((0..rows).map(|i| entry(m, i, j).clone()))))
        }
        Some(Shape::Vector(n)) if n > 0 => list(m.args().iter().map(|x| list([x.clone()]))),
        _ => return Ok(None),
    };
    stepped(ev, "tran", e, result)
}

/// `IdentityMatrix[n]`.
pub(super) fn identity_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.as_i64().and_then(|n| usize::try_from(n).ok()) else {
        return Ok(None);
    };
    if n == 0 || n.saturating_mul(n) > ev.context().limits().terms {
        return Ok(None);
    }
    let result = list((0..n).map(|i| list((0..n).map(|j| Expr::from(i64::from(i == j))))));
    stepped(ev, "idn", e, result)
}

/// `Tr[m]`.
pub(super) fn tr(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    if rows != columns {
        not_square(ev, e);
        return Ok(None);
    }
    let result = Expr::call(
        Symbol::PLUS,
        (0..rows)
            .map(|i| entry(m, i, i).clone())
            .collect::<Vec<_>>(),
    );
    stepped(ev, "trace", e, result)
}

/// The conjugate of the entry `x`: of a number, exactly or of its ball; `x`
/// itself when it is known real; `None` otherwise.
fn conjugate(ctx: &Context, x: &Expr) -> Option<Expr> {
    match x.view() {
        View::ComplexBall(z) => return Some(Expr::from(z.conj())),
        View::Ball(_) => return Some(x.clone()),
        _ => {}
    }
    if let Some(z) = x.to_complex() {
        return Some(Expr::from(&z.conj()));
    }
    (element(ctx, x, Set::Reals) == Truth::True).then(|| x.clone())
}

/// `Dot[u, v]`.
///
/// # Errors
///
/// [`EvalError::Shape`] for two vectors or matrices of other dimensions.
pub(super) fn dot(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [u, v] = e.args() else {
        return Ok(None);
    };
    let pairs: Vec<(&Expr, &Expr)> = match (shape(u), shape(v)) {
        (Some(Shape::Vector(n)), Some(Shape::Vector(k))) if n == k => {
            u.args().iter().zip(v.args()).collect()
        }
        (Some(Shape::Matrix(m, n)), Some(Shape::Matrix(k, p))) if (m, n) == (k, p) => u
            .args()
            .iter()
            .zip(v.args())
            .flat_map(|(a, b)| a.args().iter().zip(b.args()))
            .collect(),
        (Some(_), Some(_)) => return Err(EvalError::Shape),
        _ => return Ok(None),
    };
    let ctx = ev.context();
    let mut terms = Vec::with_capacity(pairs.len());
    for (a, b) in pairs {
        let Some(a) = conjugate(ctx, a) else {
            return Ok(None);
        };
        terms.push(Expr::call(Symbol::TIMES, [a, b.clone()]));
    }
    stepped(ev, "dot", e, Expr::call(Symbol::PLUS, terms))
}

/// `Cross[u, v]`.
///
/// # Errors
///
/// [`EvalError::Shape`] for vectors of other lengths.
pub(super) fn cross(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [u, v] = e.args() else {
        return Ok(None);
    };
    if !is_list(u) || !is_list(v) {
        return Ok(None);
    }
    let three = |w: &Expr| -> Option<[Expr; 3]> {
        match (shape(w)?, w.args()) {
            (Shape::Vector(3), [a, b, c]) => Some([a.clone(), b.clone(), c.clone()]),
            (Shape::Vector(2), [a, b]) => Some([a.clone(), b.clone(), Expr::from(0)]),
            _ => None,
        }
    };
    let (Some([u1, u2, u3]), Some([v1, v2, v3])) = (three(u), three(v)) else {
        return if shape(u).is_some() && shape(v).is_some() {
            Err(EvalError::Shape)
        } else {
            Ok(None)
        };
    };
    let minus = |a: Expr, b: Expr, c: Expr, d: Expr| a * b - c * d;
    let result = list([
        minus(u2.clone(), v3.clone(), u3.clone(), v2.clone()),
        minus(u3, v1.clone(), u1.clone(), v3),
        minus(u1, v2, u2, v1),
    ]);
    stepped(ev, "cross", e, result)
}

/// The square matrix of `e`'s one argument and its size; `None`, after a
/// warning when it is not square, otherwise.
fn square(ev: &mut Eval<'_>, e: &Expr) -> Option<(Expr, usize)> {
    let [m] = e.args() else {
        return None;
    };
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return None;
    };
    if rows != columns {
        not_square(ev, e);
        return None;
    }
    Some((m.clone(), rows))
}

/// The entries of `m` row by row.
fn entries(m: &Expr) -> Vec<Expr> {
    m.args()
        .iter()
        .flat_map(|r| r.args().iter().cloned())
        .collect()
}

/// The rows of `m` as vectors of entries.
fn rows_of(m: &Expr) -> Vec<Vec<Expr>> {
    m.args().iter().map(|r| r.args().to_vec()).collect()
}

/// Gives the warning `message` on `e`.
fn warn(ev: &mut Eval<'_>, message: impl Into<std::borrow::Cow<'static, str>>, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// Says that a product could pass the limit on terms.
fn too_large(ev: &mut Eval<'_>, e: &Expr) {
    let terms = ev.context().limits().terms;
    warn(ev, format!("eliminating could pass {terms} terms"), e);
}

/// The fractions of a conversion over ℤ.
fn integer_fractions(p: &Polynomials) -> Vec<(&Sparse<Integer>, &Sparse<Integer>)> {
    p.fractions()
        .iter()
        .map(|f| match f {
            Fraction::Integers { num, den } => (num, den),
            Fraction::GaussianIntegers { .. } => unreachable!("one ring"),
        })
        .collect()
}

/// The fractions of a conversion over ℤ\[i\].
fn gaussian_fractions(p: &Polynomials) -> Vec<(&Sparse<Complex>, &Sparse<Complex>)> {
    p.fractions()
        .iter()
        .map(|f| match f {
            Fraction::GaussianIntegers { num, den } => (num, den),
            Fraction::Integers { .. } => unreachable!("one ring"),
        })
        .collect()
}

/// `a / b` in lowest terms, built.
fn ratio<R: Numbers>(
    p: &Polynomials,
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Expr, EvalError> {
    if a.is_zero() {
        return Ok(Expr::from(0));
    }
    let (n, d) = cancelled(ring, a, b, cancel)?;
    Ok(p.ratio_expr(&n, &d))
}

/// `Det[m]`.
pub(super) fn det(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, n)) = square(ev, e) else {
        return Ok(None);
    };
    let result = if numeric(&m) {
        let Some(d) = numeric_det(ev, rows_of(&m))? else {
            warn(ev, "a pivot cannot be told from 0", e);
            return Ok(None);
        };
        d
    } else {
        let Some(p) = convert(ev, e, &entries(&m), true)? else {
            return Ok(None);
        };
        let k = p.generators().len();
        let found = match &p.fractions()[0] {
            Fraction::Integers { .. } => exact_det(
                ev,
                &p,
                &SparseRing::new(Integers, k),
                &integer_fractions(&p),
                n,
            )?,
            Fraction::GaussianIntegers { .. } => exact_det(
                ev,
                &p,
                &SparseRing::new(GaussianIntegers, k),
                &gaussian_fractions(&p),
                n,
            )?,
        };
        let Some(d) = found else {
            too_large(ev, e);
            return Ok(None);
        };
        d
    };
    stepped(ev, "det", e, result)
}

/// The determinant of the `n × n` matrix of `fractions`: Bareiss's on the
/// scaled rows, over the product of the scales.
fn exact_det<R: Numbers>(
    ev: &Eval<'_>,
    p: &Polynomials,
    ring: &SparseRing<R>,
    fractions: &Entries<'_, R::Elem>,
    n: usize,
) -> Result<Option<Expr>, EvalError> {
    let (rows, scales) = scaled_rows(ring, fractions, n);
    let Some(d) = bareiss(ring, rows, ev.context().limits().terms) else {
        return Ok(None);
    };
    let scale = scales.iter().fold(ring.one(), |a, s| ring.mul(&a, s));
    ratio(p, ring, &d, &scale, ev.context().interrupt_flag()).map(Some)
}

/// `Inverse[m]`.
pub(super) fn inverse(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, _)) = square(ev, e) else {
        return Ok(None);
    };
    match inverse_of(ev, e, &m)? {
        Some(result) => stepped(ev, "inv", e, result),
        None => Ok(None),
    }
}

/// The inverse of the square matrix `m`; `None` after a warning on `e` when
/// it is singular, its numeric pivots undecided, or eliminating too large.
pub(super) fn inverse_of(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some(Shape::Matrix(n, columns)) = shape(m) else {
        return Ok(None);
    };
    if n != columns {
        return Err(EvalError::Shape);
    }
    if numeric(m) {
        let mut rows: Vec<Vec<Expr>> = rows_of(m)
            .into_iter()
            .enumerate()
            .map(|(i, mut r)| {
                r.extend((0..n).map(|j| Expr::from(i64::from(i == j))));
                r
            })
            .collect();
        return Ok(match numeric_gauss_jordan(ev, &mut rows, n)? {
            Numeric::Pivots(pivots) if pivots.len() == n => {
                Some(list(rows.into_iter().map(|r| list(r.into_iter().skip(n)))))
            }
            Numeric::Pivots(_) => {
                warn(ev, "the matrix is singular", e);
                None
            }
            Numeric::Undecided => {
                warn(ev, "a pivot cannot be told from 0", e);
                None
            }
        });
    }
    let Some(p) = convert(ev, e, &entries(m), true)? else {
        return Ok(None);
    };
    let k = p.generators().len();
    let found = match &p.fractions()[0] {
        Fraction::Integers { .. } => exact_inverse(
            ev,
            &p,
            &SparseRing::new(Integers, k),
            &integer_fractions(&p),
            n,
        )?,
        Fraction::GaussianIntegers { .. } => exact_inverse(
            ev,
            &p,
            &SparseRing::new(GaussianIntegers, k),
            &gaussian_fractions(&p),
            n,
        )?,
    };
    Ok(match found {
        Inverted::Matrix(result) => Some(result),
        Inverted::Singular => {
            warn(ev, "the matrix is singular", e);
            None
        }
        Inverted::TooLarge => {
            too_large(ev, e);
            None
        }
    })
}

/// How an inversion ended.
enum Inverted {
    Matrix(Expr),
    Singular,
    TooLarge,
}

/// The inverse of the `n × n` matrix of `fractions`: Gauss–Jordan on the
/// scaled rows and `scale · I`, each entry of the right half over its row's
/// pivot.
fn exact_inverse<R: Numbers>(
    ev: &Eval<'_>,
    p: &Polynomials,
    ring: &SparseRing<R>,
    fractions: &Entries<'_, R::Elem>,
    n: usize,
) -> Result<Inverted, EvalError> {
    let (mut rows, scales) = scaled_rows(ring, fractions, n);
    for (i, (row, scale)) in rows.iter_mut().zip(&scales).enumerate() {
        row.extend((0..n).map(|j| if i == j { scale.clone() } else { ring.zero() }));
    }
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let Some(pivots) = gauss_jordan(ring, &mut rows, n, cancel, ev.context().limits().terms)?
    else {
        return Ok(Inverted::TooLarge);
    };
    if pivots.len() < n {
        return Ok(Inverted::Singular);
    }
    let mut result = Vec::with_capacity(n);
    for (j, row) in rows.iter().enumerate() {
        let entries = (0..n)
            .map(|k| ratio(p, ring, &row[n + k], &row[j], cancel))
            .collect::<Result<Vec<_>, _>>()?;
        result.push(list(entries));
    }
    Ok(Inverted::Matrix(list(result)))
}

mod elimination;
#[cfg(test)]
mod tests;
