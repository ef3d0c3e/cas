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
//! - `RowReduce`, Xcas's `rref`; `MatrixRank`, `rank`; `NullSpace`, `ker`;
//!   `ColumnSpace`, `image`; `LinSolve`, `linsolve`: from the reduced row
//!   echelon form of [`elimination`], generic in the parameters.
//! - `CharacteristicPolynomial`, Xcas's `charpoly`: `det(x I − m)`, by
//!   Bareiss's for exact entries and the Faddeev–Le Verrier recurrence for
//!   numeric ones; `Eigenvalues`, `eigenvals`: its roots as `CSolve`'s, each
//!   as often as its multiplicity; `Eigenvectors`, `eigenvects`, by
//!   [`eigen`].
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
    match determinant(ev, e, &m, n)? {
        Some(d) => stepped(ev, "det", e, d),
        None => Ok(None),
    }
}

/// The determinant of the `n × n` matrix `m`, built; `None` after a warning
/// on `e` when a numeric pivot is undecided or eliminating too large.
fn determinant(ev: &mut Eval<'_>, e: &Expr, m: &Expr, n: usize) -> Result<Option<Expr>, EvalError> {
    if numeric(m) {
        let d = numeric_det(ev, rows_of(m))?;
        if d.is_none() {
            warn(ev, "a pivot cannot be told from 0", e);
        }
        return Ok(d);
    }
    let Some(p) = convert(ev, e, &entries(m), true)? else {
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
    if found.is_none() {
        too_large(ev, e);
    }
    Ok(found)
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

/// A matrix in reduced row echelon form: its rows, built, and its pivot
/// columns, the pivot of the `k`-th in row `k`.
struct Reduced {
    rows: Vec<Vec<Expr>>,
    pivots: Vec<usize>,
}

/// The reduced row echelon form of the matrix `m`; `None` after a warning
/// on `e` when its numeric pivots are undecided or eliminating is too
/// large.
fn reduced(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Result<Option<Reduced>, EvalError> {
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    if numeric(m) {
        let mut entries = rows_of(m);
        return Ok(match numeric_gauss_jordan(ev, &mut entries, columns)? {
            Numeric::Pivots(pivots) => {
                for row in &mut entries[pivots.len()..] {
                    row.fill(Expr::from(0));
                }
                Some(Reduced {
                    rows: entries,
                    pivots,
                })
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
        Fraction::Integers { .. } => exact_reduced(
            ev,
            &p,
            &SparseRing::new(Integers, k),
            &integer_fractions(&p),
            rows,
            columns,
        )?,
        Fraction::GaussianIntegers { .. } => exact_reduced(
            ev,
            &p,
            &SparseRing::new(GaussianIntegers, k),
            &gaussian_fractions(&p),
            rows,
            columns,
        )?,
    };
    if found.is_none() {
        too_large(ev, e);
    }
    Ok(found)
}

/// The reduced row echelon form of the `rows × columns` matrix of
/// `fractions`: Gauss–Jordan on its scaled rows, each pivot row over its
/// pivot in lowest terms, the rows past the rank 0; `None` when too large.
fn exact_reduced<R: Numbers>(
    ev: &Eval<'_>,
    p: &Polynomials,
    ring: &SparseRing<R>,
    fractions: &Entries<'_, R::Elem>,
    rows: usize,
    columns: usize,
) -> Result<Option<Reduced>, EvalError> {
    let (mut a, _) = scaled_rows(ring, fractions, columns);
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let Some(pivots) = gauss_jordan(ring, &mut a, columns, cancel, ev.context().limits().terms)?
    else {
        return Ok(None);
    };
    let mut result = Vec::with_capacity(rows);
    for (i, row) in a.iter().enumerate() {
        result.push(match pivots.get(i) {
            Some(&j) => row
                .iter()
                .map(|x| ratio(p, ring, x, &row[j], cancel))
                .collect::<Result<Vec<_>, _>>()?,
            None => vec![Expr::from(0); columns],
        });
    }
    Ok(Some(Reduced {
        rows: result,
        pivots,
    }))
}

/// The one matrix argument of `e`.
fn matrix_argument(e: &Expr) -> Option<&Expr> {
    match e.args() {
        [m] if matches!(shape(m), Some(Shape::Matrix(..))) => Some(m),
        _ => None,
    }
}

/// `RowReduce[m]`, Xcas's `rref`.
pub(super) fn row_reduce(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(m) = matrix_argument(e) else {
        return Ok(None);
    };
    let Some(r) = reduced(ev, e, m)? else {
        return Ok(None);
    };
    stepped(ev, "rref", e, list(r.rows.into_iter().map(list)))
}

/// `MatrixRank[m]`, Xcas's `rank`: generic in the parameters.
pub(super) fn matrix_rank(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(m) = matrix_argument(e) else {
        return Ok(None);
    };
    let Some(r) = reduced(ev, e, m)? else {
        return Ok(None);
    };
    let rank = i64::try_from(r.pivots.len()).expect("a rank in memory");
    stepped(ev, "rank", e, Expr::from(rank))
}

/// `NullSpace[m]`, Xcas's `ker`: for each free column `f`, the vector with
/// `−1` there, the entries of column `f` in the pivot rows at their pivot
/// columns, and 0 elsewhere, as giac's `ker([[1,2,3],[2,4,6]])` is
/// `[[2,-1,0],[3,0,-1]]`.
pub(super) fn null_space(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(m) = matrix_argument(e) else {
        return Ok(None);
    };
    let Some(r) = reduced(ev, e, m)? else {
        return Ok(None);
    };
    let columns = m.args()[0].args().len();
    let basis = (0..columns).filter(|j| !r.pivots.contains(j)).map(|f| {
        let mut v = vec![Expr::from(0); columns];
        v[f] = Expr::from(-1);
        for (i, &p) in r.pivots.iter().enumerate() {
            v[p] = r.rows[i][f].clone();
        }
        list(v)
    });
    stepped(ev, "ker", e, list(basis.collect::<Vec<_>>()))
}

/// `ColumnSpace[m]`, Xcas's `image`: the nonzero rows of the reduced row
/// echelon form of the transpose.
pub(super) fn column_space(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(m) = matrix_argument(e) else {
        return Ok(None);
    };
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    let transposed = list((0..columns).map(|j| list((0..rows).map(|i| entry(m, i, j).clone()))));
    let Some(r) = reduced(ev, e, &transposed)? else {
        return Ok(None);
    };
    let rank = r.pivots.len();
    stepped(
        ev,
        "image",
        e,
        list(r.rows.into_iter().take(rank).map(list)),
    )
}

/// `LinSolve[eqns, vars]` and `LinSolve[a, b]`, Xcas's `linsolve`: the
/// solution of a linear system as one list, `[]` for none. Of equations,
/// that of `Solve`, free unknowns kept; of a matrix and a vector, the
/// vector `x` with `a x = b`, a family staying with a warning, where giac
/// names its constants `C_0`, ….
pub(super) fn lin_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [first, second] = e.args() else {
        return Ok(None);
    };
    if let (Some(Shape::Matrix(rows, columns)), Some(Shape::Vector(k))) =
        (shape(first), shape(second))
    {
        if k != rows {
            return Err(EvalError::Shape);
        }
        return matrix_lin_solve(ev, e, first, second, columns);
    }
    let equations = if is_list(first) {
        first.args().to_vec()
    } else {
        vec![first.clone()]
    };
    let unknowns = if is_list(second) {
        second.args().to_vec()
    } else {
        vec![second.clone()]
    };
    let distinct = unknowns
        .iter()
        .enumerate()
        .all(|(i, x)| x.as_symbol().is_some() && !unknowns[..i].contains(x));
    if unknowns.is_empty() || !distinct {
        return Ok(None);
    }
    let Some(solutions) = super::solve::linear_system(ev, e, &equations, &unknowns)? else {
        return Ok(None);
    };
    let result = match solutions.args() {
        [tuple] => tuple.clone(),
        _ => list([]),
    };
    stepped(ev, "linsolve", e, result)
}

/// `LinSolve[a, b]` for an `a` of `columns` columns: the reduced row echelon
/// form of `[a | b]`.
fn matrix_lin_solve(ev: &mut Eval<'_>, e: &Expr, a: &Expr, b: &Expr, columns: usize) -> Native {
    let augmented = list(
        a.args()
            .iter()
            .zip(b.args())
            .map(|(row, y)| list(row.args().iter().cloned().chain([y.clone()]))),
    );
    let Some(r) = reduced(ev, e, &augmented)? else {
        return Ok(None);
    };
    if r.pivots.contains(&columns) {
        return stepped(ev, "linsolve", e, list([]));
    }
    if r.pivots.len() < columns {
        warn(
            ev,
            "a family of solutions: LinSolve[eqns, vars] names its unknowns",
            e,
        );
        return Ok(None);
    }
    let x = (0..columns).map(|i| r.rows[i][columns].clone());
    stepped(ev, "linsolve", e, list(x.collect::<Vec<_>>()))
}

/// The variable of a characteristic polynomial computed for its roots or
/// coefficients.
fn fresh() -> Expr {
    Expr::from(Symbol::new("CharacteristicPolynomial$x"))
}

/// `det(x I − m)` in the variable [`fresh`], evaluated and expanded, and
/// its coefficients from the leading one, 1, to the constant; `None` after
/// a warning on `e`. Exact entries: the determinant of `x I − m` by
/// [`determinant`]; numeric ones: the Faddeev–Le Verrier recurrence.
fn char_poly(
    ev: &mut Eval<'_>,
    e: &Expr,
    m: &Expr,
    n: usize,
) -> Result<Option<(Expr, Vec<Expr>)>, EvalError> {
    let x = fresh();
    if numeric(m) {
        let c = faddeev(ev, &rows_of(m))?;
        let poly = ev.eval(&Expr::call(
            Symbol::PLUS,
            c.iter()
                .enumerate()
                .map(|(i, ci)| {
                    let degree = i64::try_from(n - i).expect("a degree in memory");
                    ci.clone() * x.clone().pow(degree)
                })
                .collect::<Vec<_>>(),
        ))?;
        return Ok(Some((poly, c)));
    }
    let shifted = list(m.args().iter().enumerate().map(|(i, row)| {
        list(row.args().iter().enumerate().map(|(j, a)| {
            if i == j {
                x.clone() - a.clone()
            } else {
                -a.clone()
            }
        }))
    }));
    let shifted = ev.eval(&shifted)?;
    let Some(d) = determinant(ev, e, &shifted, n)? else {
        return Ok(None);
    };
    let poly = ev.eval(&d)?;
    // Its coefficients in x, by the conversion.
    let Some(p) = convert(ev, e, std::slice::from_ref(&poly), true)? else {
        return Ok(None);
    };
    let Some(kx) = p.variable(&x) else {
        return Ok(None);
    };
    let k = p.generators().len();
    let coefficients = match &p.fractions()[0] {
        Fraction::Integers { num, den } => {
            coefficients_in(&p, &SparseRing::new(Integers, k), num, den, kx, n)
        }
        Fraction::GaussianIntegers { num, den } => {
            coefficients_in(&p, &SparseRing::new(GaussianIntegers, k), num, den, kx, n)
        }
    };
    let c = coefficients
        .into_iter()
        .map(|c| ev.eval(&c))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some((poly, c)))
}

/// The coefficients of `num / den` in the variable `kx`, from degree `n`
/// down to 0, built.
fn coefficients_in<R: Numbers>(
    p: &Polynomials,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
    kx: usize,
    n: usize,
) -> Vec<Expr> {
    (0..=n)
        .rev()
        .map(|d| {
            let c = ring.coeff(num, kx, u64::try_from(d).expect("a degree in memory"));
            p.ratio_expr(&c, den)
        })
        .collect()
}

/// The coefficients of the characteristic polynomial of the numbers `a`
/// by the Faddeev–Le Verrier recurrence: `c_n = 1`, `M₁ = I`,
/// `c_{n−k} = −tr(A Mₖ)/k`, `Mₖ₊₁ = A Mₖ + c_{n−k} I`.
fn faddeev(ev: &mut Eval<'_>, a: &[Vec<Expr>]) -> Result<Vec<Expr>, EvalError> {
    let n = a.len();
    let identity = |c: &Expr| -> Vec<Vec<Expr>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { c.clone() } else { Expr::from(0) })
                    .collect()
            })
            .collect()
    };
    let product = |ev: &mut Eval<'_>,
                   x: &[Vec<Expr>],
                   y: &[Vec<Expr>]|
     -> Result<Vec<Vec<Expr>>, EvalError> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let terms: Vec<Expr> =
                            (0..n).map(|l| x[i][l].clone() * y[l][j].clone()).collect();
                        ev.eval(&Expr::call(Symbol::PLUS, terms))
                    })
                    .collect()
            })
            .collect()
    };
    let mut c = vec![Expr::from(1)];
    let mut m = identity(&Expr::from(1));
    for k in 1..=n {
        let am = product(ev, a, &m)?;
        let trace = Expr::call(
            Symbol::PLUS,
            (0..n).map(|i| am[i][i].clone()).collect::<Vec<_>>(),
        );
        let k_expr = Expr::from(i64::try_from(k).expect("a size in memory"));
        let ck = ev.eval(&(-(trace / k_expr)))?;
        if k < n {
            m = am;
            for (i, row) in m.iter_mut().enumerate() {
                row[i] = ev.eval(&(row[i].clone() + ck.clone()))?;
            }
        }
        c.push(ck);
    }
    Ok(c)
}

/// `CharacteristicPolynomial[m, x]`, Xcas's `charpoly`: `det(x I − m)`,
/// monic in `x`; `CharacteristicPolynomial[m]`, its coefficients from the
/// leading one, as giac's `poly1[1,-5,-2]`.
pub(super) fn characteristic_polynomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (m, x) = match e.args() {
        [m] => (m, None),
        [m, x] if x.as_symbol().is_some() && !m.contains(x) => (m, Some(x)),
        _ => return Ok(None),
    };
    let Some(Shape::Matrix(n, columns)) = shape(m) else {
        return Ok(None);
    };
    if n != columns {
        not_square(ev, e);
        return Ok(None);
    }
    let Some((poly, c)) = char_poly(ev, e, m, n)? else {
        return Ok(None);
    };
    let result = match x {
        None => list(c),
        Some(x) => poly.substitute(&[(fresh(), x.clone())]),
    };
    stepped(ev, "charpoly", e, result)
}

/// `Eigenvalues[m]`, Xcas's `eigenvals`: the roots of the characteristic
/// polynomial, complex ones included, each as often as its multiplicity
/// (CAS-DOC-12 §4).
pub(super) fn eigenvalues(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, n)) = square(ev, e) else {
        return Ok(None);
    };
    let Some((poly, _)) = char_poly(ev, e, &m, n)? else {
        return Ok(None);
    };
    let Some(values) = super::solve::roots_counted(ev, e, &poly, &fresh())? else {
        return Ok(None);
    };
    stepped(ev, "eigenvals", e, list(values))
}

pub(super) mod constructors;
pub(super) mod decompose;
mod eigen;
mod elimination;
#[cfg(test)]
mod tests;

pub(super) use eigen::eigenvectors;
