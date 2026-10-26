//! Elimination on matrices (CAS-DOC-12 §2, §3), for `Det`, `Inverse` and
//! the row reductions that follow.
//!
//! - **Exact entries** are converted together (CAS-DOC-8 §2) to fractions
//!   over ℤ, or ℤ\[i\] with a complex number, in their generators, balls
//!   among them; each row is multiplied by the lcm of its denominators
//!   ([`scaled_rows`]).
//! - [`bareiss`]: the determinant by Bareiss's fraction-free elimination
//!   \[Bar68\], each step's division by the previous pivot exact, a row
//!   swap negating.
//! - [`gauss_jordan`]: for each column a nonzero pivot of the fewest terms,
//!   each other row `p r − c (pivot row)` over the gcd of its entries, as
//!   CAS-DOC-11 §5. Generic in the parameters: a pivot is a nonzero
//!   polynomial.
//! - **Numeric entries**, numbers among which a ball: elimination on their
//!   values, the pivot of largest magnitude ([`numeric_det`],
//!   [`numeric_gauss_jordan`]); a pivot whose ball holds 0 without being 0
//!   is undecided.
//!
//! Past `Limits::terms` for a product, `None`: the caller stays with a
//! warning.

use cas_num::algebra::Ring;
use cas_num::ball::Float;
use cas_poly::{Cancel, Sparse, SparseRing};

use super::super::approximation;
use super::super::polynomial::Numbers;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};

/// A row of polynomials.
pub(in crate::native) type Row<E> = Vec<Sparse<E>>;

/// The entries of a matrix row by row, as numerators and denominators.
pub(super) type Entries<'a, E> = [(&'a Sparse<E>, &'a Sparse<E>)];

/// Rows of polynomials and the scales they were multiplied by.
pub(super) type Scaled<E> = (Vec<Row<E>>, Vec<Sparse<E>>);

/// The rows of a matrix whose `fractions` are its entries row by row,
/// each row times the lcm of its denominators, and those lcms.
pub(super) fn scaled_rows<R: Numbers>(
    ring: &SparseRing<R>,
    fractions: &Entries<'_, R::Elem>,
    columns: usize,
) -> Scaled<R::Elem> {
    let mut rows = Vec::new();
    let mut scales = Vec::new();
    for row in fractions.chunks(columns) {
        let lcm = row.iter().fold(ring.one(), |l, (_, d)| ring.lcm(&l, d));
        rows.push(
            row.iter()
                .map(|(n, d)| ring.mul(n, &ring.div_exact(&lcm, d).expect("a divisor of the lcm")))
                .collect(),
        );
        scales.push(lcm);
    }
    (rows, scales)
}

/// Whether the product of `a` and `b` could pass `terms` terms.
fn large<E>(a: &Sparse<E>, b: &Sparse<E>, terms: usize) -> bool {
    a.term_count().saturating_mul(b.term_count()) > terms
}

/// The determinant of the square `a` by Bareiss's elimination; `None` when
/// a product could pass `terms` terms.
pub(super) fn bareiss<R: Numbers>(
    ring: &SparseRing<R>,
    mut a: Vec<Row<R::Elem>>,
    terms: usize,
) -> Option<Sparse<R::Elem>> {
    let n = a.len();
    let mut previous = ring.one();
    let mut negative = false;
    for k in 0..n {
        let Some(r) = (k..n)
            .filter(|&i| !a[i][k].is_zero())
            .min_by_key(|&i| a[i][k].term_count())
        else {
            return Some(ring.zero());
        };
        if r != k {
            a.swap(r, k);
            negative = !negative;
        }
        for i in k + 1..n {
            for j in k + 1..n {
                if large(&a[k][k], &a[i][j], terms) || large(&a[i][k], &a[k][j], terms) {
                    return None;
                }
                let t = ring.sub(&ring.mul(&a[k][k], &a[i][j]), &ring.mul(&a[i][k], &a[k][j]));
                a[i][j] = ring
                    .div_exact(&t, &previous)
                    .expect("Bareiss's divisions are exact");
            }
        }
        previous = a[k][k].clone();
    }
    let d = a[n - 1][n - 1].clone();
    Some(if negative { ring.neg(&d) } else { d })
}

/// Gauss–Jordan elimination of `rows` on their first `columns` entries:
/// the pivot columns in order, the pivot of the `k`-th in row `k`; `None`
/// when a product could pass `terms` terms.
///
/// # Errors
///
/// [`EvalError::Interrupted`] when `cancel` says so.
pub(in crate::native) fn gauss_jordan<R: Numbers>(
    ring: &SparseRing<R>,
    rows: &mut [Row<R::Elem>],
    columns: usize,
    cancel: &dyn Cancel,
    terms: usize,
) -> Result<Option<Vec<usize>>, EvalError> {
    let width = rows.first().map_or(0, Vec::len);
    let mut pivots: Vec<usize> = Vec::new();
    for j in 0..columns {
        let rank = pivots.len();
        let Some(best) = (rank..rows.len())
            .filter(|&i| !rows[i][j].is_zero())
            .min_by_key(|&i| (rows[i][j].term_count(), rows[i][j].degrees()))
        else {
            continue;
        };
        rows.swap(rank, best);
        for i in 0..rows.len() {
            if i == rank || rows[i][j].is_zero() {
                continue;
            }
            let (pivot, c) = (&rows[rank][j], &rows[i][j]);
            if (0..width)
                .any(|k| large(pivot, &rows[i][k], terms) || large(c, &rows[rank][k], terms))
            {
                return Ok(None);
            }
            let row: Row<R::Elem> = (0..width)
                .map(|k| ring.sub(&ring.mul(pivot, &rows[i][k]), &ring.mul(c, &rows[rank][k])))
                .collect();
            rows[i] = primitive(ring, row, cancel)?;
        }
        pivots.push(j);
    }
    Ok(Some(pivots))
}

/// `row` over the gcd of its entries.
///
/// # Errors
///
/// [`EvalError::Interrupted`] when `cancel` says so.
pub(in crate::native) fn primitive<R: Numbers>(
    ring: &SparseRing<R>,
    row: Row<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Row<R::Elem>, EvalError> {
    let mut g = ring.zero();
    for entry in &row {
        g = ring
            .gcd_with(&g, entry, cancel)
            .map_err(|_| EvalError::Interrupted)?;
        if ring.is_one(&g) {
            return Ok(row);
        }
    }
    if g.is_zero() {
        return Ok(row);
    }
    Ok(row
        .iter()
        .map(|entry| ring.div_exact(entry, &g).expect("a divisor of each entry"))
        .collect())
}

/// Whether the entries of `m` are numbers, a ball among them: the
/// numeric elimination's.
pub(super) fn numeric(m: &Expr) -> bool {
    let entries = || m.args().iter().flat_map(Expr::args);
    entries().all(Expr::is_number)
        && entries().any(|x| matches!(x.view(), View::Ball(_) | View::ComplexBall(_)))
}

/// The magnitude of a number's midpoint, for partial pivoting.
fn magnitude(ctx: &Context, x: &Expr) -> Float {
    let ball = approximation(ctx, x, 64).unwrap_or_else(|_| x.clone());
    match ball.view() {
        View::Ball(b) => b.mid().abs(),
        View::ComplexBall(z) => z.re().mid().abs().max(z.im().mid().abs()),
        _ => Float::ZERO,
    }
}

/// Whether the number `x` is exactly 0, and whether its ball holds 0.
fn zero(x: &Expr) -> (bool, bool) {
    match x.view() {
        View::Ball(b) => (
            b.is_exact() && b.mid().is_zero(),
            b.contains_rational(&cas_num::Rational::ZERO),
        ),
        View::ComplexBall(z) => (
            z.is_exact() && z.re().mid().is_zero() && z.im().mid().is_zero(),
            z.may_be_zero(),
        ),
        _ => {
            let exact = *x == Expr::from(0);
            (exact, exact)
        }
    }
}

/// How a numeric elimination ended.
pub(super) enum Numeric {
    /// The pivot columns, as [`gauss_jordan`]'s.
    Pivots(Vec<usize>),
    /// A column whose candidates hold 0 without all being 0.
    Undecided,
}

/// Gauss–Jordan elimination on the numbers of `rows`, their first `columns`
/// entries, each pivot the largest in magnitude and its row divided by it,
/// each entry evaluated.
///
/// # Errors
///
/// What evaluating the entries gives.
pub(super) fn numeric_gauss_jordan(
    ev: &mut Eval<'_>,
    rows: &mut [Vec<Expr>],
    columns: usize,
) -> Result<Numeric, EvalError> {
    let width = rows.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    for j in 0..columns {
        let rank = pivots.len();
        let ctx = ev.context();
        let Some(best) = (rank..rows.len())
            .max_by(|&a, &b| magnitude(ctx, &rows[a][j]).cmp(&magnitude(ctx, &rows[b][j])))
        else {
            break;
        };
        match zero(&rows[best][j]) {
            (true, _) => continue,
            (false, true) => return Ok(Numeric::Undecided),
            (false, false) => {}
        }
        rows.swap(rank, best);
        let pivot = rows[rank][j].clone();
        for x in &mut rows[rank] {
            *x = ev.eval(&(x.clone() / pivot.clone()))?;
        }
        for i in 0..rows.len() {
            if i == rank || zero(&rows[i][j]).0 {
                continue;
            }
            let c = rows[i][j].clone();
            for k in 0..width {
                let t = rows[i][k].clone() - c.clone() * rows[rank][k].clone();
                rows[i][k] = ev.eval(&t)?;
            }
        }
        pivots.push(j);
    }
    Ok(Numeric::Pivots(pivots))
}

/// The determinant of the square numbers `rows`, by elimination with
/// partial pivoting; `None` when a pivot is undecided.
///
/// # Errors
///
/// What evaluating the entries gives.
pub(super) fn numeric_det(
    ev: &mut Eval<'_>,
    mut rows: Vec<Vec<Expr>>,
) -> Result<Option<Expr>, EvalError> {
    let n = rows.len();
    let mut det = Expr::from(1);
    for k in 0..n {
        let ctx = ev.context();
        let best = (k..n)
            .max_by(|&a, &b| magnitude(ctx, &rows[a][k]).cmp(&magnitude(ctx, &rows[b][k])))
            .expect("a row");
        match zero(&rows[best][k]) {
            (true, _) => return Ok(Some(Expr::from(0))),
            (false, true) => return Ok(None),
            (false, false) => {}
        }
        if best != k {
            rows.swap(best, k);
            det = -det;
        }
        let pivot = rows[k][k].clone();
        det = ev.eval(&(det * pivot.clone()))?;
        for i in k + 1..n {
            let c = ev.eval(&(rows[i][k].clone() / pivot.clone()))?;
            for j in k + 1..n {
                let t = rows[i][j].clone() - c.clone() * rows[k][j].clone();
                rows[i][j] = ev.eval(&t)?;
            }
        }
    }
    Ok(Some(det))
}
