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
//! - **`LDLDecomposition[A]`** (CAS-477): `{p, L, D}` for a symmetric `A =
//!   L D Lᵀ`, `L` unit lower triangular, `D` diagonal, by elimination
//!   without pivoting, `p` the identity permutation from 0, giac's form; a
//!   zero pivot with a nonzero entry below it leaves the call with a
//!   warning (giac takes a 2 × 2 block there).
//! - **`PositiveDefiniteQ[A]`**: a symmetric `A` whose pivots `d_j` are all
//!   above 0; a pivot shown 0 or below makes it `False`, as a matrix that is
//!   not symmetric does; a pivot not decided leaves the call. giac's
//!   `isposdef([[1,2],[2,1]])` is `true`, wrongly: its eigenvalues are 3
//!   and −1.
//! - **`Inertia[A]`**: `{positive, negative, zero}` eigenvalues of a real
//!   symmetric `A`, by Descartes' rule of signs on its characteristic
//!   polynomial, exact as every root is real: the zeros the trailing zero
//!   coefficients, the positive roots the sign changes of the rest, the
//!   negative ones those of `p(−x)`.

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

/// Whether `a` is square and symmetric, as far as `is_zero` decides.
fn symmetric(ev: &mut Eval<'_>, a: &[Vec<Expr>]) -> Result<Option<bool>, EvalError> {
    let n = a.len();
    if a[0].len() != n {
        return Ok(None);
    }
    for (i, row) in a.iter().enumerate() {
        for (j, x) in row.iter().enumerate().take(i) {
            match is_zero(ev, &(x.clone() - a[j][i].clone()))? {
                Some(true) => {}
                Some(false) => return Ok(Some(false)),
                None => return Ok(None),
            }
        }
    }
    Ok(Some(true))
}

/// How the elimination of a symmetric matrix ended.
enum Ldl {
    /// `L` and the pivots.
    Done(Vec<Vec<Expr>>, Vec<Expr>),
    /// A pivot shown 0 with a nonzero entry below it, at this index.
    ZeroPivot,
    /// A value not decided.
    Undecided,
}

/// `A = L D Lᵀ` for a symmetric `a`, without pivoting.
fn ldl(ev: &mut Eval<'_>, a: &[Vec<Expr>]) -> Result<Ldl, EvalError> {
    let n = a.len();
    let mut l = identity(n);
    let mut d: Vec<Expr> = Vec::with_capacity(n);
    for j in 0..n {
        let terms: Vec<Expr> = (0..j).map(|k| l[j][k].pow(2) * d[k].clone()).collect();
        let pivot = simplified(ev, a[j][j].clone() - Expr::call(Symbol::PLUS, terms))?;
        let Some(zero) = is_zero(ev, &pivot)? else {
            return Ok(Ldl::Undecided);
        };
        for i in j + 1..n {
            let terms: Vec<Expr> = (0..j)
                .map(|k| l[i][k].clone() * l[j][k].clone() * d[k].clone())
                .collect();
            let rest = simplified(ev, a[i][j].clone() - Expr::call(Symbol::PLUS, terms))?;
            if zero {
                match is_zero(ev, &rest)? {
                    Some(true) => l[i][j] = Expr::from(0),
                    Some(false) => return Ok(Ldl::ZeroPivot),
                    None => return Ok(Ldl::Undecided),
                }
            } else {
                l[i][j] = simplified(ev, rest / pivot.clone())?;
            }
        }
        d.push(pivot);
    }
    Ok(Ldl::Done(l, d))
}

/// `LDLDecomposition[A]`.
pub(in crate::native) fn ldl_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    if symmetric(ev, &a)? != Some(true) {
        warn(ev, "the matrix is not shown symmetric", e);
        return Ok(None);
    }
    let n = a.len();
    let (l, d) = match ldl(ev, &a)? {
        Ldl::Done(l, d) => (l, d),
        Ldl::ZeroPivot => {
            warn(ev, "a zero pivot with entries below it", e);
            return Ok(None);
        }
        Ldl::Undecided => {
            warn(ev, "a pivot cannot be told from 0", e);
            return Ok(None);
        }
    };
    let mut diagonal = vec![vec![Expr::from(0); n]; n];
    for (i, x) in d.into_iter().enumerate() {
        diagonal[i][i] = x;
    }
    let p = list(
        (0..n)
            .map(|i| Expr::from(i64::try_from(i).expect("small")))
            .collect(),
    );
    let result = list(vec![p, matrix(l), matrix(diagonal)]);
    ev.step("ldl", e, &result);
    Ok(Some(result))
}

/// The sign of `x` when decided: exactly, or by a ball.
fn sign(ev: &mut Eval<'_>, x: &Expr) -> Result<Option<Ordering>, EvalError> {
    if let Some(q) = x.to_rational() {
        return Ok(Some(q.cmp(&cas_num::Rational::ZERO)));
    }
    if is_zero(ev, x)? == Some(true) {
        return Ok(Some(Ordering::Equal));
    }
    let value = crate::native::approximation(ev.context(), x, 128)?;
    Ok(match value.view() {
        crate::expr::View::Ball(b) => b.sign(),
        _ => None,
    })
}

/// `PositiveDefiniteQ[A]`.
pub(in crate::native) fn positive_definite_q(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    let answer = |t: bool| {
        Ok(Some(Expr::from(if t {
            Symbol::TRUE
        } else {
            Symbol::FALSE
        })))
    };
    match symmetric(ev, &a)? {
        Some(true) => {}
        Some(false) => return answer(false),
        None => return Ok(None),
    }
    let d = match ldl(ev, &a)? {
        Ldl::Done(_, d) => d,
        // A zero pivot: not definite.
        Ldl::ZeroPivot => return answer(false),
        Ldl::Undecided => return Ok(None),
    };
    for x in &d {
        match sign(ev, x)? {
            Some(Ordering::Greater) => {}
            Some(_) => return answer(false),
            None => return Ok(None),
        }
    }
    answer(true)
}

/// The sign changes of nonzero rational coefficients.
fn changes(c: &[cas_num::Rational]) -> i64 {
    let signs: Vec<bool> = c
        .iter()
        .filter(|q| !q.is_zero())
        .map(cas_num::Rational::is_negative)
        .collect();
    let count = signs.windows(2).filter(|w| w[0] != w[1]).count();
    i64::try_from(count).expect("small")
}

/// `Inertia[A]`.
pub(in crate::native) fn inertia(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(m) = rows(a) else {
        return Ok(None);
    };
    if symmetric(ev, &m)? != Some(true) {
        return Ok(None);
    }
    let x = Expr::symbol("Inertia$x");
    let p = ev.eval(&Expr::call(
        Symbol::CHARACTERISTIC_POLYNOMIAL,
        [a.clone(), x.clone()],
    ))?;
    let Some(c) = coefficients(ev, &p, &x)? else {
        return Ok(None);
    };
    let Some(mut c) = c.iter().map(Expr::to_rational).collect::<Option<Vec<_>>>() else {
        return Ok(None);
    };
    // Leading first: the trailing zeros are the zero roots.
    let mut zeros = 0i64;
    while c.len() > 1 && c.last().is_some_and(cas_num::Rational::is_zero) {
        c.pop();
        zeros += 1;
    }
    let positive = changes(&c);
    // p(−x): the coefficient of xᵏ times (−1)ᵏ, k counted from the end.
    let degree = c.len() - 1;
    let mirrored: Vec<cas_num::Rational> = c
        .iter()
        .enumerate()
        .map(|(i, q)| {
            if (degree - i) % 2 == 1 {
                -q.clone()
            } else {
                q.clone()
            }
        })
        .collect();
    let negative = changes(&mirrored);
    let result = list(vec![
        Expr::from(positive),
        Expr::from(negative),
        Expr::from(zeros),
    ]);
    ev.step("inertia", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
