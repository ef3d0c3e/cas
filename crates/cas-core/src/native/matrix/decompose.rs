//! Matrix decompositions and functions (CAS-DOC-30 §2, §3, CAS-344 to
//! CAS-346), at the level of expressions, each entry put together
//! (`Together`), each zero test exact or by a ball told from 0:
//!
//! - **`LUDecomposition[A]`** (`lu`): `{p, L, U}` for a square `A`, by
//!   Gaussian elimination with a row exchange where a pivot is 0; `p` lists
//!   the rows from 0, giac's `[1,0]`, so that `A[p] = L U`.
//! - **`QRDecomposition[A]`** (`qr`): `{Q, R}` by Gram–Schmidt on the
//!   columns, with the conjugate dot product and the norms `√(vᴴv)`; a
//!   column that depends on the ones before leaves the call.
//! - **`CholeskyDecomposition[A]`** (`cholesky`): `L` with `A = L Lᵀ` for a
//!   symmetric `A` whose pivots `a_jj − Σ l_jk²` are shown positive.
//! - **`JordanDecomposition[A]`** (`jordan`): `{P, J}` with `A = P J P⁻¹`,
//!   from the exact eigenvalues (`Eigenvalues`) and, for each `λ`, chains of
//!   generalized eigenvectors from the kernels of `(A − λI)^k`, longest
//!   first; `P` holds the chains, eigenvector first, and `J` a block per
//!   chain.
//! - **`MinimalPolynomial[A]`** (`pmin`): `Π (x − λ)^(k_λ)`, `k_λ` the
//!   longest chain of `λ`: giac's coefficient list, or the polynomial in
//!   the variable given.
//! - **`MatrixExp[A]`** (`expm`) and **`MatrixPower[A, n]`** (`matpow`):
//!   `P f(J) P⁻¹`, `f` of a block the upper triangular Toeplitz matrix of
//!   `f(λ), f′(λ), …, f^(m−1)(λ)/(m − 1)!`: `e^λ` times `1, 1, 1/2, …`, and
//!   `λ^n` times `1, n/λ, C(n, 2)/λ², …`. A block at `λ = 0` of size `m`
//!   is 0 for `n ≥ m`, generic in `n`, which giac writes with `0^n`. An
//!   integer `n` is `Power`'s. giac's `exp` of a square matrix is its
//!   exponential, where this CAS's `Exp` of a list is entrywise: the
//!   exponential is `expm`.
//!
//! Matrices larger than [`MAX_SIZE`] leave the call, as do eigenvalues
//! without a closed form.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear algebra: A = L U = Q R, P J P⁻¹, entries i, j, k"
)]

use cas_num::Rational;

use super::{conjugate, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest size decomposed.
pub(super) const MAX_SIZE: usize = 12;

type Matrix = Vec<Vec<Expr>>;

/// Jordan chains: each eigenvalue with its vectors, eigenvector first.
type Chains = Vec<(Expr, Vec<Expr>)>;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn matrix(rows: Matrix) -> Expr {
    list(rows.into_iter().map(list).collect())
}

/// The rows of a matrix argument; `None` for another argument.
fn rows(e: &Expr) -> Option<Matrix> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) || e.args().is_empty() {
        return None;
    }
    let width = e.args()[0].args().len();
    let rows: Matrix = e
        .args()
        .iter()
        .map(|r| {
            (r.head() == Some(&Expr::from(Symbol::LIST)) && r.args().len() == width)
                .then(|| r.args().to_vec())
        })
        .collect::<Option<_>>()?;
    (width > 0 && rows.len() <= MAX_SIZE && width <= MAX_SIZE).then_some(rows)
}

/// `x` put together.
fn simplified(ev: &mut Eval<'_>, x: Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::TOGETHER, [x]))
}

/// Whether `x` is 0: exactly, or `None` when a ball cannot tell.
fn is_zero(ev: &mut Eval<'_>, x: &Expr) -> Result<Option<bool>, EvalError> {
    if *x == Expr::from(0) {
        return Ok(Some(true));
    }
    if x.is_number() {
        return Ok(Some(match x.view() {
            View::Ball(b) => {
                if b.contains_rational(&Rational::ZERO) {
                    return Ok(None);
                }
                false
            }
            View::ComplexBall(z) => {
                if z.may_be_zero() {
                    return Ok(None);
                }
                false
            }
            _ => false,
        }));
    }
    let together = simplified(ev, x.clone())?;
    if together == Expr::from(0) {
        return Ok(Some(true));
    }
    // A numeric expression that Together leaves: its ball.
    let value = approximation(ev.context(), &together, 128)?;
    Ok(match value.view() {
        View::Ball(b) if b.contains_rational(&Rational::ZERO) => None,
        View::ComplexBall(z) if z.may_be_zero() => None,
        _ => Some(false),
    })
}

fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| Expr::from(i64::from(i == j))).collect())
        .collect()
}

pub(in crate::native) fn lu_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(mut u) = rows(a) else {
        return Ok(None);
    };
    let n = u.len();
    if u[0].len() != n {
        return Ok(None);
    }
    let mut l = identity(n);
    let mut p: Vec<usize> = (0..n).collect();
    for k in 0..n {
        // A pivot: the first row from k whose entry is not 0.
        let mut pivot = None;
        for (i, row) in u.iter().enumerate().skip(k) {
            match is_zero(ev, &row[k])? {
                Some(false) => {
                    pivot = Some(i);
                    break;
                }
                Some(true) => {}
                None => {
                    warn(ev, "a pivot cannot be told from 0", e);
                    return Ok(None);
                }
            }
        }
        let Some(i) = pivot else {
            continue;
        };
        if i != k {
            u.swap(i, k);
            p.swap(i, k);
            // The multipliers so far, before column k, go with their rows;
            // the pivot's row `i` is below `k`.
            let (top, bottom) = l.split_at_mut(i);
            top[k][..k].swap_with_slice(&mut bottom[0][..k]);
        }
        for i in k + 1..n {
            let factor = simplified(ev, u[i][k].clone() / u[k][k].clone())?;
            l[i][k] = factor.clone();
            let (top, bottom) = u.split_at_mut(i);
            for (x, y) in bottom[0][k..].iter_mut().zip(&top[k][k..]) {
                *x = simplified(ev, x.clone() - factor.clone() * y.clone())?;
            }
        }
    }
    let p = list(
        p.into_iter()
            .map(|i| Expr::from(i64::try_from(i).expect("small")))
            .collect(),
    );
    let result = list(vec![p, matrix(l), matrix(u)]);
    ev.step("lu", e, &result);
    Ok(Some(result))
}

/// `vᴴ w`, the conjugate of `v`'s entries; `None` when one is not known.
fn inner(ev: &mut Eval<'_>, v: &[Expr], w: &[Expr]) -> Result<Option<Expr>, EvalError> {
    let mut terms = Vec::with_capacity(v.len());
    for (x, y) in v.iter().zip(w) {
        let Some(c) = conjugate(ev.context(), x) else {
            return Ok(None);
        };
        terms.push(c * y.clone());
    }
    simplified(ev, Expr::call(Symbol::PLUS, terms)).map(Some)
}

pub(in crate::native) fn qr_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    let (m, n) = (a.len(), a[0].len());
    if m < n {
        return Ok(None);
    }
    let column = |j: usize| -> Vec<Expr> { a.iter().map(|r| r[j].clone()).collect() };
    let mut q: Vec<Vec<Expr>> = Vec::with_capacity(n);
    // R by columns, transposed at the end.
    let mut r_columns: Matrix = Vec::with_capacity(n);
    for j in 0..n {
        let aj = column(j);
        let mut v = aj.clone();
        let mut rj = vec![Expr::from(0); n];
        for (rij, qi) in rj.iter_mut().zip(&q) {
            let Some(value) = inner(ev, qi, &aj)? else {
                return Ok(None);
            };
            for (vk, qk) in v.iter_mut().zip(qi) {
                *vk = simplified(ev, vk.clone() - value.clone() * qk.clone())?;
            }
            *rij = value;
        }
        let Some(square) = inner(ev, &v, &v)? else {
            return Ok(None);
        };
        match is_zero(ev, &square)? {
            Some(false) => {}
            Some(true) => {
                warn(ev, "a column depends on the ones before", e);
                return Ok(None);
            }
            None => {
                warn(ev, "a norm cannot be told from 0", e);
                return Ok(None);
            }
        }
        let norm = simplified(ev, Expr::call(Symbol::SQRT, [square]))?;
        let mut qj = Vec::with_capacity(m);
        for vk in v {
            qj.push(simplified(ev, vk / norm.clone())?);
        }
        rj[j] = norm;
        r_columns.push(rj);
        q.push(qj);
    }
    let q_rows: Matrix = (0..m)
        .map(|i| q.iter().map(|c| c[i].clone()).collect())
        .collect();
    let r: Matrix = (0..n)
        .map(|i| r_columns.iter().map(|c| c[i].clone()).collect())
        .collect();
    let result = list(vec![matrix(q_rows), matrix(r)]);
    ev.step("qr", e, &result);
    Ok(Some(result))
}

/// Whether `x` is shown positive: a positive number, by the assumptions,
/// or by its ball.
fn positive(ev: &mut Eval<'_>, x: &Expr) -> Result<bool, EvalError> {
    let test = ev.eval(&Expr::call(Symbol::GREATER, [x.clone(), Expr::from(0)]))?;
    if test.as_symbol() == Some(Symbol::TRUE) {
        return Ok(true);
    }
    let value = approximation(ev.context(), x, 128)?;
    Ok(matches!(value.view(), View::Ball(b) if b.sign() == Some(std::cmp::Ordering::Greater)))
}

pub(in crate::native) fn cholesky_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(a) = rows(a) else {
        return Ok(None);
    };
    let n = a.len();
    if a[0].len() != n {
        return Ok(None);
    }
    for (i, row) in a.iter().enumerate() {
        for (j, x) in row.iter().enumerate().take(i) {
            if is_zero(ev, &(x.clone() - a[j][i].clone()))? != Some(true) {
                warn(ev, "the matrix is not symmetric", e);
                return Ok(None);
            }
        }
    }
    let mut l = vec![vec![Expr::from(0); n]; n];
    for j in 0..n {
        let squares: Vec<Expr> = (0..j).map(|k| l[j][k].pow(2)).collect();
        let pivot = simplified(ev, a[j][j].clone() - Expr::call(Symbol::PLUS, squares))?;
        if !positive(ev, &pivot)? {
            warn(ev, "a pivot is not shown positive", e);
            return Ok(None);
        }
        let ljj = simplified(ev, Expr::call(Symbol::SQRT, [pivot]))?;
        for i in j + 1..n {
            let products: Vec<Expr> = (0..j).map(|k| l[i][k].clone() * l[j][k].clone()).collect();
            l[i][j] = simplified(
                ev,
                (a[i][j].clone() - Expr::call(Symbol::PLUS, products)) / ljj.clone(),
            )?;
        }
        l[j][j] = ljj;
    }
    let result = matrix(l);
    ev.step("cholesky", e, &result);
    Ok(Some(result))
}

/// The distinct eigenvalues of `a` with their multiplicities, exact;
/// `None` when one has no closed form.
fn eigenvalues(ev: &mut Eval<'_>, a: &Expr) -> Result<Option<Vec<(Expr, usize)>>, EvalError> {
    let values = ev.eval(&Expr::call(Symbol::EIGENVALUES, [a.clone()]))?;
    if values.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let mut distinct: Vec<(Expr, usize)> = Vec::new();
    for v in values.args() {
        if v.contains(&Expr::from(Symbol::new("Root")))
            || v.contains(&Expr::from(Symbol::EIGENVALUES))
        {
            return Ok(None);
        }
        let mut found = false;
        for (w, m) in &mut distinct {
            if is_zero(ev, &(v.clone() - w.clone()))? == Some(true) {
                *m += 1;
                found = true;
                break;
            }
        }
        if !found {
            distinct.push((v.clone(), 1));
        }
    }
    Ok(Some(distinct))
}

/// The product of matrices, or of a matrix and a vector: `Times`, as `Dot`
/// of two matrices is their inner product, giac's `dot`.
fn product(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::TIMES, [a.clone(), b.clone()]))
}

/// The basis of the kernel of `m`, as vectors.
fn kernel(ev: &mut Eval<'_>, m: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    let k = ev.eval(&Expr::call(Symbol::NULL_SPACE, [m.clone()]))?;
    if k.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    Ok(Some(k.args().to_vec()))
}

/// The rank of the vectors, as rows.
fn rank(ev: &mut Eval<'_>, vectors: &[Expr]) -> Result<Option<usize>, EvalError> {
    if vectors.is_empty() {
        return Ok(Some(0));
    }
    let r = ev.eval(&Expr::call(Symbol::MATRIX_RANK, [list(vectors.to_vec())]))?;
    Ok(r.as_i64().and_then(|r| usize::try_from(r).ok()))
}

/// The Jordan chains of `a`: for each one, its eigenvalue and its vectors,
/// eigenvector first; `None` where the eigenvalues or kernels fail.
fn chains(ev: &mut Eval<'_>, a: &Expr, n: usize) -> Result<Option<Chains>, EvalError> {
    let Some(values) = eigenvalues(ev, a)? else {
        return Ok(None);
    };
    let ident = matrix(identity(n));
    let mut all = Vec::new();
    for (lambda, multiplicity) in values {
        let b = ev.eval(&(a.clone() - lambda.clone() * ident.clone()))?;
        // The kernels of B^k up to the multiplicity's dimension.
        let mut kernels: Vec<Vec<Expr>> = vec![Vec::new()];
        let mut power = b.clone();
        loop {
            let Some(k) = kernel(ev, &power)? else {
                return Ok(None);
            };
            let dimension = k.len();
            kernels.push(k);
            if dimension >= multiplicity || kernels.len() > n + 1 {
                break;
            }
            power = product(ev, &power, &b)?;
        }
        if kernels.last().map(Vec::len) != Some(multiplicity) {
            return Ok(None);
        }
        // Chains from the longest: tops in K_s outside K_{s−1} and the
        // images of the chains already taken.
        let s = kernels.len() - 1;
        let mut taken: Vec<Vec<Expr>> = Vec::new();
        for size in (1..=s).rev() {
            // What K_size must be spanned beyond: K_{size−1} and the
            // vectors of the chains taken at this level.
            let mut span: Vec<Expr> = kernels[size - 1].clone();
            for chain in &taken {
                let top = chain.len() - 1;
                if top + 1 >= size {
                    span.push(chain[size - 1].clone());
                }
            }
            let Some(mut current) = rank(ev, &span)? else {
                return Ok(None);
            };
            for v in kernels[size].clone() {
                let mut with = span.clone();
                with.push(v.clone());
                let Some(r) = rank(ev, &with)? else {
                    return Ok(None);
                };
                if r > current {
                    // A new chain of length `size`, eigenvector first.
                    let mut chain = vec![v.clone()];
                    for _ in 1..size {
                        let next = product(ev, &b, chain.last().expect("one"))?;
                        chain.push(next);
                    }
                    chain.reverse();
                    span = with;
                    current = r;
                    taken.push(chain);
                }
            }
        }
        all.extend(taken.into_iter().map(|chain| (lambda.clone(), chain)));
    }
    let count: usize = all.iter().map(|(_, c)| c.len()).sum();
    Ok((count == n).then_some(all))
}

/// The square matrix argument, `n ≤ MAX_SIZE`.
fn square_matrix(e: &Expr) -> Option<(Expr, usize)> {
    let a = e.args().first()?;
    let r = rows(a)?;
    (r.len() == r[0].len()).then(|| (a.clone(), r.len()))
}

/// `P` from the chains, and the Jordan matrix of `f`'s blocks:
/// `f(λ, k)` the entry `k` places above the diagonal.
fn assembled(
    ev: &mut Eval<'_>,
    chains: &[(Expr, Vec<Expr>)],
    n: usize,
    f: &mut dyn FnMut(&mut Eval<'_>, &Expr, usize) -> Result<Expr, EvalError>,
) -> Result<(Expr, Expr), EvalError> {
    let mut columns: Vec<Expr> = Vec::with_capacity(n);
    let mut j = vec![vec![Expr::from(0); n]; n];
    let mut start = 0;
    for (lambda, chain) in chains {
        let m = chain.len();
        for r in 0..m {
            for c in r..m {
                j[start + r][start + c] = f(ev, lambda, c - r)?;
            }
        }
        columns.extend(chain.iter().cloned());
        start += m;
    }
    let p: Matrix = (0..n)
        .map(|i| columns.iter().map(|v| v.args()[i].clone()).collect())
        .collect();
    Ok((matrix(p), matrix(j)))
}

pub(in crate::native) fn jordan_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let Some(chains) = chains(ev, &a, n)? else {
        return Ok(None);
    };
    let mut block = |_: &mut Eval<'_>, lambda: &Expr, k: usize| -> Result<Expr, EvalError> {
        Ok(match k {
            0 => lambda.clone(),
            1 => Expr::from(1),
            _ => Expr::from(0),
        })
    };
    let (p, j) = assembled(ev, &chains, n, &mut block)?;
    let result = list(vec![p, j]);
    ev.step("jordan", e, &result);
    Ok(Some(result))
}

pub(in crate::native) fn minimal_polynomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let x = match e.args() {
        [_] => None,
        [_, x] if x.as_symbol().is_some() => Some(x.clone()),
        _ => return Ok(None),
    };
    let Some(chains) = chains(ev, &a, n)? else {
        return Ok(None);
    };
    let variable = x.clone().unwrap_or_else(|| Expr::symbol("x"));
    let mut longest: Vec<(Expr, usize)> = Vec::new();
    for (lambda, chain) in &chains {
        match longest.iter_mut().find(|(l, _)| l == lambda) {
            Some((_, k)) => *k = (*k).max(chain.len()),
            None => longest.push((lambda.clone(), chain.len())),
        }
    }
    let factors: Vec<Expr> = longest
        .into_iter()
        .map(|(lambda, k)| (variable.clone() - lambda).pow(i64::try_from(k).expect("small")))
        .collect();
    let p = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [Expr::call(Symbol::TIMES, factors)],
    ))?;
    let result = match x {
        Some(_) => p,
        None => ev.eval(&Expr::call(Symbol::COEFFICIENT_LIST, [p, variable]))?,
    };
    ev.step("pmin", e, &result);
    Ok(Some(result))
}

/// `P F P⁻¹`, entries put together.
fn conjugated(ev: &mut Eval<'_>, p: &Expr, f: &Expr) -> Result<Option<Expr>, EvalError> {
    let inverse = ev.eval(&Expr::call(Symbol::INVERSE, [p.clone()]))?;
    if inverse.head() == Some(&Expr::from(Symbol::INVERSE)) {
        return Ok(None);
    }
    let product = ev.eval(&Expr::call(Symbol::TIMES, [p.clone(), f.clone(), inverse]))?;
    let Some(r) = rows(&product) else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(r.len());
    for row in r {
        let mut cells = Vec::with_capacity(row.len());
        for x in row {
            cells.push(simplified(ev, x)?);
        }
        out.push(cells);
    }
    Ok(Some(matrix(out)))
}

pub(in crate::native) fn matrix_exp(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    if e.args().len() != 1 {
        return Ok(None);
    }
    let Some(chains) = chains(ev, &a, n)? else {
        return Ok(None);
    };
    // e^λ / k!.
    let mut block = |ev: &mut Eval<'_>, lambda: &Expr, k: usize| -> Result<Expr, EvalError> {
        let factorial: i64 = (1..=i64::try_from(k).expect("small")).product();
        ev.eval(&(Expr::call(Symbol::EXP, [lambda.clone()]) / Expr::from(factorial)))
    };
    let (p, f) = assembled(ev, &chains, n, &mut block)?;
    let Some(result) = conjugated(ev, &p, &f)? else {
        return Ok(None);
    };
    ev.step("expm", e, &result);
    Ok(Some(result))
}

pub(in crate::native) fn matrix_power(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, k] = e.args() else {
        return Ok(None);
    };
    if k.as_i64().is_some() {
        return ev
            .eval(&Expr::call(Symbol::POWER, [a.clone(), k.clone()]))
            .map(Some);
    }
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let Some(chains) = chains(ev, &a, n)? else {
        return Ok(None);
    };
    // C(k, j) λ^(k − j); 0 at λ = 0, for k past the block's size.
    let exponent = k.clone();
    let mut block = |ev: &mut Eval<'_>, lambda: &Expr, j: usize| -> Result<Expr, EvalError> {
        if is_zero(ev, lambda)? == Some(true) {
            return Ok(Expr::from(0));
        }
        let j = i64::try_from(j).expect("small");
        let binomial = Expr::call(Symbol::BINOMIAL, [exponent.clone(), Expr::from(j)]);
        ev.eval(&(binomial * lambda.pow(exponent.clone() - Expr::from(j))))
    };
    let (p, f) = assembled(ev, &chains, n, &mut block)?;
    let Some(result) = conjugated(ev, &p, &f)? else {
        return Ok(None);
    };
    ev.step("matpow", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
