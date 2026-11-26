//! giac's minors, echelon forms and pivots (CAS-530, CAS-528).
//!
//! - `DetMinor[m]`, Xcas's `det_minor`: the determinant by cofactor
//!   (Laplace) expansion along the first row, each minor of the last rows
//!   on a set of columns computed once and expanded, for `n ≤ MAX_MINOR`.
//! - `AdjointMatrix[m]`, Xcas's `adjoint_matrix`: `[c, [B₀, …, B_{n−1}]]`
//!   of the Faddeev–Le Verrier recurrence, `c` the characteristic
//!   polynomial's coefficients from the leading 1, `B₀ = I`,
//!   `Bₖ = m Bₖ₋₁ + cₖ I`, expanded.
//! - `RowEchelon[m]`, Xcas's `ref`, and `RowReduce[m, KeepPivot]`, Xcas's
//!   `rref(m, keep_pivot)`: giac 2.0.0.21's `mrref` and `rref_reduce`
//!   (src/vecteur.cc), ported. Each row times the lcm of its denominators,
//!   over the gcd of its entries; for each column the pivot the least
//!   complex as giac's `islesscomplexthan` (src/gen.cc), from the current
//!   row down, the first among equals; each other row
//!   `(p r − c (pivot row)) / previous pivot`, fraction-free (Bareiss).
//!   `ref` eliminates below the pivots, the last row taking none, then
//!   divides each row by its first nonzero entry (among all but the last
//!   column for an `n × (n+1)` matrix, giac's `mdividebypivot`).
//!   `keep_pivot` divides nothing: with at least as many columns as rows it
//!   eliminates below, and when no column was skipped and the last pivot
//!   `d` is nonzero substitutes back, each pivot `d` (Cramer-like);
//!   otherwise it eliminates above and below from the start. Numeric or
//!   modular entries: `ref` stays, `keep_pivot` is `rref`'s. Departures
//!   kept: giac's leading coefficient of a polynomial pivot is in its own
//!   monomial order, here in ours; giac interpolates `keep_pivot` of a
//!   polynomial matrix of 5 rows and more, whose sign may differ.
//! - `Pivot[m, l, c]`, Xcas's `pivot`: each row `i ≠ l` becomes
//!   `m[l][c] row_i − m[i][c] row_l`, row `l` and rows with `m[i][c] = 0`
//!   kept as giac does (1.8 T8: `pivot([[0,0],[0,1]],0,0)`); indices from 0,
//!   negative ones from the end.
//! - `JordanMatrix[m]`, Xcas's `egvl` and `eigVl`: the `J` of
//!   `JordanDecomposition` ([`super::decompose`]), its blocks in giac 2.0's
//!   order: real part descending, then imaginary part ascending (by 64-bit
//!   balls), each eigenvalue's longest block first. A matrix already in
//!   Jordan form is its own `J`, as giac's.

use std::cmp::Ordering;

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_num::ball::{Float, Round};
use cas_num::{Complex, Rational};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::super::linear::{Shape, shape};
use super::super::polynomial::{Numbers, convert};
use super::elimination::{Entries, Row, large, numeric, primitive, scaled_rows};
use super::{
    Native, entries, faddeev, gaussian_fractions, integer_fractions, list, ratio, reduced, rows_of,
    square, stepped, too_large,
};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::native::index::Index;
use crate::poly::{Coefficient, Fraction, Polynomials};
use crate::symbol::Symbol;

/// The largest size `DetMinor` expands: `2ⁿ` minors.
const MAX_MINOR: usize = 12;

/// `DetMinor[m]`.
pub(in crate::native) fn det_minor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, n)) = square(ev, e) else {
        return Ok(None);
    };
    if n > MAX_MINOR {
        return Ok(None);
    }
    let a = rows_of(&m);
    // minors[s]: the minor of the last |s| rows on the columns in s; a
    // subset before the sets that hold it.
    let full = (1_usize << n) - 1;
    let mut minors: Vec<Expr> = Vec::with_capacity(full + 1);
    minors.push(Expr::from(1));
    for s in 1..=full {
        let row = n - s.count_ones() as usize;
        let terms: Vec<Expr> = (0..n)
            .filter(|&j| s >> j & 1 == 1)
            .enumerate()
            .map(|(t, j)| {
                let term = a[row][j].clone() * minors[s & !(1 << j)].clone();
                if t % 2 == 0 { term } else { -term }
            })
            .collect();
        let minor = ev.eval(&Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::PLUS, terms)],
        ))?;
        minors.push(minor);
    }
    let d = minors.swap_remove(full);
    stepped(ev, "det_minor", e, d)
}

/// `AdjointMatrix[m]`.
pub(in crate::native) fn adjoint_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, _)) = square(ev, e) else {
        return Ok(None);
    };
    let (c, matrices) = faddeev(ev, &rows_of(&m))?;
    let matrices = matrices.into_iter().map(|b| list(b.into_iter().map(list)));
    let result = list([list(c), list(matrices)]);
    let result = ev.eval(&Expr::call(Symbol::EXPAND, [result]))?;
    stepped(ev, "adjoint_matrix", e, result)
}

/// `Pivot[m, l, c]`.
///
/// # Errors
///
/// [`EvalError::IndexOutOfRange`] for an index past the matrix.
pub(in crate::native) fn pivot(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, l, c] = e.args() else {
        return Ok(None);
    };
    let Some(Shape::Matrix(rows, columns)) = shape(m) else {
        return Ok(None);
    };
    let (Some(Index::One(l)), Some(Index::One(c))) = (Index::of(l, rows)?, Index::of(c, columns)?)
    else {
        return Ok(None);
    };
    let a = rows_of(m);
    let p = &a[l][c];
    let result = list(a.iter().enumerate().map(|(i, row)| {
        if i == l {
            return list(row.iter().cloned());
        }
        let q = &row[c];
        // giac skips a row with 0 in the pivot's column (CAS-528, 1.8 T8).
        if *q == Expr::from(0) {
            return list(row.iter().cloned());
        }
        list((0..columns).map(|k| p.clone() * row[k].clone() - q.clone() * a[l][k].clone()))
    }));
    stepped(ev, "pivot", e, result)
}

/// `JordanMatrix[m]`.
pub(in crate::native) fn jordan_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [_] = e.args() else {
        return Ok(None);
    };
    let Some((a, n)) = super::decompose::square_matrix(e) else {
        return Ok(None);
    };
    if in_jordan_form(&rows_of(&a)) {
        return stepped(ev, "egvl", e, a);
    }
    let Some(chains) = super::decompose::chains(ev, &a, n)? else {
        return Ok(None);
    };
    let mut keyed = Vec::with_capacity(chains.len());
    for (lambda, chain) in chains {
        let key = plane(ev, &lambda)?;
        keyed.push((key, lambda, chain.len()));
    }
    // giac's order: real part descending, then imaginary part ascending;
    // the sort is stable, so each eigenvalue's chains stay longest first.
    // An eigenvalue without a numeric value keeps this CAS's order.
    if keyed.iter().all(|(key, ..)| key.is_some()) {
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-12 * x.abs().max(y.abs()).max(1.0);
        keyed.sort_by(|(u, ..), (v, ..)| {
            let (Some((ur, ui)), Some((vr, vi))) = (u, v) else {
                return Ordering::Equal;
            };
            if close(*ur, *vr) {
                if close(*ui, *vi) {
                    Ordering::Equal
                } else {
                    ui.total_cmp(vi)
                }
            } else {
                vr.total_cmp(ur)
            }
        });
    }
    let mut j = vec![vec![Expr::from(0); n]; n];
    let mut start = 0;
    for (_, lambda, size) in keyed {
        for r in 0..size {
            j[start + r][start + r] = lambda.clone();
            if r + 1 < size {
                j[start + r][start + r + 1] = Expr::from(1);
            }
        }
        start += size;
    }
    stepped(ev, "egvl", e, list(j.into_iter().map(list)))
}

/// Whether `a` is a Jordan matrix: upper bidiagonal, each entry above the
/// diagonal 0, or 1 between equal diagonal entries.
fn in_jordan_form(a: &[Vec<Expr>]) -> bool {
    let (zero, one) = (Expr::from(0), Expr::from(1));
    a.iter().enumerate().all(|(i, row)| {
        row.iter().enumerate().all(|(k, x)| {
            if k == i {
                true
            } else if k == i + 1 {
                *x == zero || (*x == one && a[i][i] == a[k][k])
            } else {
                *x == zero
            }
        })
    })
}

/// The point of the complex plane `x` is near, from its 64-bit ball;
/// `None` when it has no numeric value.
fn plane(ev: &Eval<'_>, x: &Expr) -> Result<Option<(f64, f64)>, EvalError> {
    let value = approximation(ev.context(), x, 64)?;
    let real = |q: &Rational| Float::from_rational(q, 64, Round::Nearest).0.to_f64();
    Ok(match value.view() {
        View::Ball(b) => Some((b.mid().to_f64(), 0.0)),
        View::ComplexBall(z) => Some((z.re().mid().to_f64(), z.im().mid().to_f64())),
        _ => value.to_complex().map(|z| (real(z.re()), real(z.im()))),
    })
}

/// `RowEchelon[m]`.
pub(in crate::native) fn row_echelon(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m] = e.args() else {
        return Ok(None);
    };
    if !matches!(shape(m), Some(Shape::Matrix(..)))
        || numeric(m)
        || super::super::modular::matrices::modular_matrix(m).is_some()
    {
        return Ok(None);
    }
    match echelon(ev, e, m, false)? {
        Some(result) => stepped(ev, "ref", e, result),
        None => Ok(None),
    }
}

/// `RowReduce[m, KeepPivot]`, the matrix `m`.
pub(super) fn keep_pivot(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Native {
    if numeric(m) || super::super::modular::matrices::modular_matrix(m).is_some() {
        let Some(r) = reduced(ev, e, m)? else {
            return Ok(None);
        };
        return stepped(ev, "rref", e, list(r.rows.into_iter().map(list)));
    }
    match echelon(ev, e, m, true)? {
        Some(result) => stepped(ev, "rref", e, result),
        None => Ok(None),
    }
}

/// giac's `ref` of the exact matrix `m`, or its `rref` with `keep_pivot`
/// when `kept`; `None` after a warning on `e`.
fn echelon(ev: &mut Eval<'_>, e: &Expr, m: &Expr, kept: bool) -> Result<Option<Expr>, EvalError> {
    let Some(Shape::Matrix(_, columns)) = shape(m) else {
        return Ok(None);
    };
    let Some(p) = convert(ev, e, &entries(m), true)? else {
        return Ok(None);
    };
    let k = p.generators().len();
    let found = match &p.fractions()[0] {
        Fraction::Integers { .. } => exact(
            ev,
            &p,
            &SparseRing::new(Integers, k),
            &integer_fractions(&p),
            columns,
            kept,
        )?,
        Fraction::GaussianIntegers { .. } => exact(
            ev,
            &p,
            &SparseRing::new(GaussianIntegers, k),
            &gaussian_fractions(&p),
            columns,
            kept,
        )?,
    };
    let Some(rows) = found else {
        too_large(ev, e);
        return Ok(None);
    };
    Ok(Some(list(rows.into_iter().map(list))))
}

/// The rows of `ref`, or of `rref` with `keep_pivot` when `kept`, of the
/// matrix of `fractions`, built; `None` when too large.
fn exact<R: Numbers>(
    ev: &Eval<'_>,
    p: &Polynomials,
    ring: &SparseRing<R>,
    fractions: &Entries<'_, R::Elem>,
    columns: usize,
    kept: bool,
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let terms = ev.context().limits().terms;
    let (scaled, _) = scaled_rows(ring, fractions, columns);
    let mut rows = scaled
        .into_iter()
        .map(|row| primitive(ring, row, cancel))
        .collect::<Result<Vec<_>, _>>()?;
    let n = rows.len();
    let found = if kept {
        kept_pivots(ring, rows, terms)
    } else {
        reduce(ring, &mut rows, false, terms).map(|_| rows)
    };
    let Some(rows) = found else {
        return Ok(None);
    };
    // giac's mdividebypivot: `ref` divides each row by its first nonzero
    // entry, the last column left out of an n × (n+1) matrix.
    let last = match (kept, n + 1 == columns) {
        (true, _) => 0,
        (false, true) => columns - 1,
        (false, false) => columns,
    };
    let one = ring.one();
    rows.iter()
        .map(|row| {
            let pivot = row[..last].iter().find(|x| !x.is_zero()).unwrap_or(&one);
            row.iter()
                .map(|x| ratio(p, ring, x, pivot, cancel))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// `rows` reduced as giac's `rref(…, keep_pivot)`: below the pivots and
/// substituted back when Cramer-like, else above and below; `None` when
/// too large.
fn kept_pivots<R: Numbers>(
    ring: &SparseRing<R>,
    mut rows: Vec<Row<R::Elem>>,
    terms: usize,
) -> Option<Vec<Row<R::Elem>>> {
    let (n, m) = (rows.len(), rows[0].len());
    if m < n {
        reduce(ring, &mut rows, true, terms)?;
        return Some(rows);
    }
    let saved = rows.clone();
    let skipped = reduce(ring, &mut rows, false, terms)?;
    let d = rows[n - 1][n - 1].clone();
    if skipped || d.is_zero() {
        rows = saved;
        reduce(ring, &mut rows, true, terms)?;
        return Some(rows);
    }
    // Row l: d times its right part, less its entries right of the pivot
    // times the right parts of the rows below, over its pivot.
    for l in (0..n - 1).rev() {
        let copy = rows[l].clone();
        let mut row = copy.clone();
        for x in &mut row[n..] {
            *x = ring.mul(&d, x);
        }
        for j in (l + 1..n).rev() {
            for k in n..m {
                row[k] = ring.sub(&row[k], &ring.mul(&copy[j], &rows[j][k]));
            }
        }
        for x in &mut row[l + 1..n] {
            *x = ring.zero();
        }
        row[l] = d.clone();
        for x in &mut row[n..] {
            *x = ring.div_exact(x, &copy[l])?;
        }
        rows[l] = row;
    }
    Some(rows)
}

/// giac's `rref_reduce` (fraction-free): for each column the least complex
/// nonzero entry ([`less_complex`]) from the current row down, swapped up;
/// each row below it, every other row when `full`, becomes
/// `(p r − c (pivot row)) / previous pivot`. Without `full` the last row
/// takes no pivot. Whether a column had no pivot; `None` when a product
/// could pass `terms` terms.
fn reduce<R: Numbers>(
    ring: &SparseRing<R>,
    rows: &mut [Row<R::Elem>],
    full: bool,
    terms: usize,
) -> Option<bool> {
    let (lmax, cmax) = (rows.len(), rows.first().map_or(0, Vec::len));
    let (mut l, mut c) = (0, 0);
    let mut previous = ring.one();
    let mut skipped = false;
    while l < lmax && c < cmax {
        if !full && l + 1 == lmax {
            break;
        }
        let mut best = l;
        for i in l + 1..lmax {
            let t = &rows[i][c];
            if !t.is_zero() && (rows[best][c].is_zero() || less_complex(t, &rows[best][c])) {
                best = i;
            }
        }
        if rows[best][c].is_zero() {
            skipped = true;
            c += 1;
            continue;
        }
        rows.swap(l, best);
        let pivot = rows[l][c].clone();
        for i in 0..lmax {
            if i == l || (!full && i < l) {
                continue;
            }
            let q = rows[i][c].clone();
            let mut row = Vec::with_capacity(cmax);
            for (x, y) in rows[i].iter().zip(&rows[l]) {
                if large(&pivot, x, terms) || large(&q, y, terms) {
                    return None;
                }
                let t = ring.sub(&ring.mul(&pivot, x), &ring.mul(&q, y));
                row.push(ring.div_exact(&t, &previous)?);
            }
            rows[i] = row;
        }
        previous = pivot;
        l += 1;
        c += 1;
    }
    Some(skipped)
}

/// The kind of a polynomial in giac's order of types: a real number, a
/// complex one, a polynomial.
fn kind<E: Coefficient>(a: &Sparse<E>) -> u8 {
    if a.term_count() != 1 || a.degrees().iter().any(|&d| d > 0) {
        return 2;
    }
    let real = a.coeffs()[0]
        .to_expr()
        .to_complex()
        .is_none_or(|z| z.im().is_zero());
    u8::from(!real)
}

/// giac's `islesscomplexthan` (src/gen.cc), its pivot order: 0 last; a
/// real number, then a complex one, then a polynomial; polynomials by
/// fewer terms, then by their leading coefficients ([`less_number`]).
fn less_complex<E: Coefficient + PartialEq>(a: &Sparse<E>, b: &Sparse<E>) -> bool {
    if a.is_zero() {
        return false;
    }
    if b.is_zero() {
        return true;
    }
    let (ka, kb) = (kind(a), kind(b));
    if ka != kb {
        return ka < kb;
    }
    if a == b {
        return false;
    }
    if a.term_count() != b.term_count() {
        return a.term_count() < b.term_count();
    }
    less_number(&a.coeffs()[0], &b.coeffs()[0])
}

/// giac's `islesscomplexthan` of two numbers: reals by absolute value, the
/// positive of two opposites first; a real before a complex; complexes by
/// absolute value, real part, then imaginary part, the larger first.
fn less_number<E: Coefficient>(x: &E, y: &E) -> bool {
    let (Some(x), Some(y)) = (x.to_expr().to_complex(), y.to_expr().to_complex()) else {
        return false;
    };
    if x.is_zero() {
        return false;
    }
    if y.is_zero() {
        return true;
    }
    let real = |z: &Complex| z.im().is_zero();
    match (real(&x), real(&y)) {
        (true, false) => true,
        (false, true) => false,
        (true, true) => {
            let (a, b) = (x.re(), y.re());
            if a == b {
                false
            } else if *a == -b {
                a.is_positive()
            } else {
                a.abs() < b.abs()
            }
        }
        (false, false) => {
            if x.norm() != y.norm() {
                x.norm() > y.norm()
            } else if x.re() != y.re() {
                x.re() > y.re()
            } else {
                x.im() > y.im()
            }
        }
    }
}
