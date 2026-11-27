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
//!   chain. Kernels over ℚ take radicals as free generators: a
//!   diagonalizable `A` whose eigenvalues have them takes its `P` from
//!   `Eigenvectors`, over ℚ\[t\]/(q(t)), giac's `jordan([[1,2],[3,4]])`.
//! - **`RationalJordanForm[A]`** (`rat_jordan`, CAS-535): `{P, R}` with
//!   `P⁻¹ A P = R` over ℚ, for a rational `A`. The characteristic
//!   polynomial (Faddeev–Le Verrier) factored over ℚ; the semisimple part
//!   `S` by Newton's iteration `S ← S − f(S) f′(S)⁻¹` from `A`, `f` the
//!   squarefree part, `N = A − S` nilpotent. For each irreducible `p` of
//!   degree `d`, nonlinear ones first by degree, then linear ones by root
//!   descending: in `V_p = ker p(S)` the chains of `N`, longest first, each
//!   top `w` the first vector of the kernel of `p(S)` and `Nʰ` (free
//!   columns in order, 1 there) outside what is spanned; its columns
//!   `N^(h−b) S^(j−1) w`, `b = 1..h`, `j = 1..d`, and its block `C(p)`
//!   along the diagonal (ones below its diagonal, `−pₖ` in its last
//!   column) with identities just above. An `A` already of that form, its
//!   blocks in any order, gives `{I, A}`, as giac's. giac's `P` differs
//!   where its own kernels do (`[[2,1],[0,3]]`: giac's `[[1,-1],[1,0]]`,
//!   here `[[1,1],[1,0]]`), and giac puts the linear factors of `x⁴ − 1`
//!   first.
//! - **`MinimalPolynomial[A]`** (`pmin`): `Π (x − λ)^(k_λ)`, `k_λ` the
//!   longest chain of `λ`: giac's coefficient list, or the polynomial in
//!   the variable given.
//! - **`MatrixExp[A]`** (`expm`) and **`MatrixPower[A, n]`** (`matpow`):
//!   `P f(J) P⁻¹`, `f` of a block the upper triangular Toeplitz matrix of
//!   `f(λ), f′(λ), …, f^(m−1)(λ)/(m − 1)!`: `Exp[λ]`, giac's `exp(λ)`,
//!   times `1, 1, 1/2, …`, and `λ^n` times `1, n/λ, C(n, 2)/λ², …`. A block at `λ = 0` of size `m`
//!   is 0 for `n ≥ m`, generic in `n`, which giac writes with `0^n`. An
//!   integer `n` is `Power`'s. giac's `exp` of a square matrix is its
//!   exponential, where this CAS's `Exp` of a list is entrywise: the
//!   exponential is `expm`.
//!
//! Matrices larger than [`MAX_SIZE`] leave the call, as do eigenvalues
//! without a closed form. A held entry is not approximated in a zero or
//! positivity test: no 2^n re-evaluation when nested (CAS-650).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear algebra: A = L U = Q R, P J P⁻¹, entries i, j, k"
)]

use std::cmp::Ordering;

use cas_num::Rational;
use cas_num::algebra::{Rationals, Ring};
use cas_poly::DenseRing;

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

/// `x` evaluated and put together; a number as it is, and a numeric
/// expression with balls as one ball: `Together` leaves a ball in its call,
/// and on a quotient of balls exceeds the recursion limit (CAS-498's
/// fuzzing: Schur and the SVD of `[[-5.]]`).
#[allow(
    clippy::needless_pass_by_value,
    reason = "the callers hand over the entries they build"
)]
fn simplified(ev: &mut Eval<'_>, x: Expr) -> Result<Expr, EvalError> {
    let x = ev.eval(&x)?;
    if x.is_number() {
        return Ok(x);
    }
    // Balls with exact radicals (a ball matrix's exact eigenvalues, the
    // fuzz target's `x³ + 27`): one ball, `Together` leaving them apart.
    if x.is_numeric() && has_ball(&x) {
        let bits = ev.context().working_precision();
        return approximation(ev.context(), &x, bits);
    }
    ev.eval(&Expr::call(Symbol::TOGETHER, [x]))
}

/// Whether a ball occurs in `x`.
fn has_ball(x: &Expr) -> bool {
    matches!(x.view(), View::Ball(_) | View::ComplexBall(_)) || x.args().iter().any(has_ball)
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
    // A numeric expression that Together leaves: its ball. Another (a
    // held call) is not evaluated again with fresh limits (CAS-650): not
    // 0, as its approximation, left as it is, says.
    if !together.is_numeric() {
        return Ok(Some(false));
    }
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
    // A held call: not evaluated again with fresh limits (CAS-650).
    if !x.is_numeric() {
        return Ok(false);
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

/// `Times[a, b, …]` applied as given: a product of matrices, which `Times`
/// keeps in its order once its operands are lists, and `Dot` does not
/// multiply (`Dot` of two matrices is their inner product, giac's `dot`).
fn times(factors: &[Expr]) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), factors)
}

/// The product of matrices, or of a matrix and a vector.
fn product(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&times(&[a.clone(), b.clone()]))
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
/// eigenvector first; `None` where the eigenvalues fail, or the kernels
/// and [`diagonal`] both.
pub(super) fn chains(ev: &mut Eval<'_>, a: &Expr, n: usize) -> Result<Option<Chains>, EvalError> {
    let Some(values) = eigenvalues(ev, a)? else {
        return Ok(None);
    };
    if let Some(all) = generalized(ev, a, n, &values)? {
        return Ok(Some(all));
    }
    diagonal(ev, a, n, &values)
}

/// The chains of a diagonalizable `a`, one vector each, from the columns of
/// `Eigenvectors`, which follow `Eigenvalues`: for eigenvalues with
/// radicals, whose kernels [`generalized`] takes over ℚ with the radicals
/// as free generators, where `Eigenvectors` takes them over ℚ\[t\]/(q(t))
/// (CAS-222). Its warnings are dropped: the call stays without them.
fn diagonal(
    ev: &mut Eval<'_>,
    a: &Expr,
    n: usize,
    values: &[(Expr, usize)],
) -> Result<Option<Chains>, EvalError> {
    let ctx = ev.context();
    let vectors = ev
        .eval_in(ctx, &Expr::call(Symbol::EIGENVECTORS, [a.clone()]))
        .0?;
    let Some(columns) = rows(&vectors).filter(|r| r.len() == n && r[0].len() == n) else {
        return Ok(None);
    };
    let mut all = Vec::with_capacity(n);
    for (lambda, multiplicity) in values {
        for _ in 0..*multiplicity {
            let j = all.len();
            if j == n {
                return Ok(None);
            }
            let v = list(columns.iter().map(|row| row[j].clone()).collect());
            all.push((lambda.clone(), vec![v]));
        }
    }
    Ok((all.len() == n).then_some(all))
}

/// The chains of `a` from the kernels of `(A − λI)^k` over the field of its
/// entries, for each of the distinct `values`; `None` where a kernel or a
/// rank fails, or the kernels fall short of a multiplicity.
fn generalized(
    ev: &mut Eval<'_>,
    a: &Expr,
    n: usize,
    values: &[(Expr, usize)],
) -> Result<Option<Chains>, EvalError> {
    let ident = matrix(identity(n));
    let mut all = Vec::new();
    for (lambda, multiplicity) in values {
        let multiplicity = *multiplicity;
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
                    // A new chain of length `size`, eigenvector first, its
                    // top's first nonzero entry positive where it is a
                    // number: giac's P of jordan([[1,1],[0,1]]) is I, where
                    // the kernels' vectors are −e₁ and −e₂.
                    let v = if leads_negative(&v) {
                        ev.eval(&-v.clone())?
                    } else {
                        v.clone()
                    };
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

/// Whether the first nonzero entry of a vector is a negative number.
fn leads_negative(v: &Expr) -> bool {
    v.args()
        .iter()
        .find(|x| **x != Expr::from(0))
        .and_then(Expr::to_rational)
        .is_some_and(|q| q.is_negative())
}

/// The square matrix argument, `n ≤ MAX_SIZE`.
pub(super) fn square_matrix(e: &Expr) -> Option<(Expr, usize)> {
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

/// A matrix of rationals, by rows.
type QMatrix = Vec<Vec<Rational>>;

fn q_identity(n: usize) -> QMatrix {
    (0..n)
        .map(|i| (0..n).map(|j| Rational::from(i64::from(i == j))).collect())
        .collect()
}

fn q_mul(a: &[Vec<Rational>], b: &[Vec<Rational>]) -> QMatrix {
    a.iter()
        .map(|row| {
            (0..b[0].len())
                .map(|j| {
                    row.iter()
                        .zip(b)
                        .fold(Rational::ZERO, |s, (x, r)| s + x * &r[j])
                })
                .collect()
        })
        .collect()
}

fn q_sub(a: &[Vec<Rational>], b: &[Vec<Rational>]) -> QMatrix {
    a.iter()
        .zip(b)
        .map(|(r, s)| r.iter().zip(s).map(|(x, y)| x - y).collect())
        .collect()
}

fn q_apply(a: &[Vec<Rational>], v: &[Rational]) -> Vec<Rational> {
    a.iter()
        .map(|row| {
            row.iter()
                .zip(v)
                .fold(Rational::ZERO, |s, (x, y)| s + x * y)
        })
        .collect()
}

/// `p(a)` for `p`'s coefficients from the constant, by Horner.
fn q_eval(p: &[Rational], a: &[Vec<Rational>]) -> QMatrix {
    let n = a.len();
    let mut r = vec![vec![Rational::ZERO; n]; n];
    for c in p.iter().rev() {
        r = q_mul(&r, a);
        for (i, row) in r.iter_mut().enumerate() {
            row[i] = &row[i] + c;
        }
    }
    r
}

/// The rows of `a` reduced (Gauss–Jordan, pivots 1) and the pivot columns.
fn q_reduced(a: &[Vec<Rational>]) -> (QMatrix, Vec<usize>) {
    let mut rows = a.to_vec();
    let columns = rows.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    for j in 0..columns {
        let rank = pivots.len();
        let Some(r) = (rank..rows.len()).find(|&i| !rows[i][j].is_zero()) else {
            continue;
        };
        rows.swap(rank, r);
        let pivot = rows[rank][j].clone();
        for x in &mut rows[rank] {
            *x = &*x / &pivot;
        }
        let pivot_row = rows[rank].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i == rank || row[j].is_zero() {
                continue;
            }
            let c = row[j].clone();
            for (x, y) in row.iter_mut().zip(&pivot_row) {
                *x = &*x - &(&c * y);
            }
        }
        pivots.push(j);
    }
    (rows, pivots)
}

/// The inverse of `a`; `None` when singular.
fn q_inverse(a: &[Vec<Rational>]) -> Option<QMatrix> {
    let n = a.len();
    let augmented: QMatrix = a
        .iter()
        .zip(q_identity(n))
        .map(|(r, e)| r.iter().cloned().chain(e).collect())
        .collect();
    let (rows, pivots) = q_reduced(&augmented);
    (pivots.len() == n && pivots[n - 1] == n - 1)
        .then(|| rows.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// A basis of the kernel of `a`, `columns` wide: for each column without a
/// pivot in order, 1 there and the negated entries of the pivot rows.
fn q_kernel(a: &[Vec<Rational>], columns: usize) -> QMatrix {
    let (rows, pivots) = q_reduced(a);
    (0..columns)
        .filter(|j| !pivots.contains(j))
        .map(|f| {
            let mut v = vec![Rational::ZERO; columns];
            v[f] = Rational::ONE;
            for (i, &p) in pivots.iter().enumerate() {
                v[p] = -&rows[i][f];
            }
            v
        })
        .collect()
}

/// Whether `v` is outside the span of `basis`, kept reduced: each vector 1
/// at its first nonzero entry, 0 there in the vectors after it. `v` joins
/// it when outside.
fn q_joins(basis: &mut Vec<Vec<Rational>>, v: &[Rational]) -> bool {
    let mut v = v.to_vec();
    for b in basis.iter() {
        let p = b.iter().position(|x| !x.is_zero()).expect("nonzero");
        if !v[p].is_zero() {
            let c = v[p].clone();
            for (x, y) in v.iter_mut().zip(b) {
                *x = &*x - &(&c * y);
            }
        }
    }
    let Some(p) = v.iter().position(|x| !x.is_zero()) else {
        return false;
    };
    let lead = v[p].clone();
    basis.push(v.iter().map(|x| x / &lead).collect());
    true
}

/// The characteristic polynomial `det(x I − a)`, coefficients from the
/// constant, by Faddeev–Le Verrier.
fn q_charpoly(a: &[Vec<Rational>]) -> Vec<Rational> {
    let n = a.len();
    let mut c = vec![Rational::ZERO; n + 1];
    c[n] = Rational::ONE;
    let mut m = vec![vec![Rational::ZERO; n]; n];
    for k in 1..=n {
        m = q_mul(a, &m);
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = &row[i] + &c[n + 1 - k];
        }
        let am = q_mul(a, &m);
        let trace = (0..n).fold(Rational::ZERO, |s, i| s + &am[i][i]);
        c[n - k] = -(trace / Rational::from(i64::try_from(k).expect("small")));
    }
    c
}

/// The companion matrix of the monic `p` (from the constant): ones below
/// the diagonal, `−pₖ` in the last column.
fn q_companion(p: &[Rational]) -> QMatrix {
    let d = p.len() - 1;
    let mut c = vec![vec![Rational::ZERO; d]; d];
    for (k, row) in c.iter_mut().enumerate() {
        if k > 0 {
            row[k - 1] = Rational::ONE;
        }
        row[d - 1] = -&p[k];
    }
    c
}

/// The rational Jordan block of `p` (monic, from the constant) with `h`
/// steps: `C(p)` along the diagonal, identities just above it.
fn q_block(p: &[Rational], h: usize) -> QMatrix {
    let d = p.len() - 1;
    let c = q_companion(p);
    let mut r = vec![vec![Rational::ZERO; h * d]; h * d];
    for b in 0..h {
        for i in 0..d {
            for j in 0..d {
                r[b * d + i][b * d + j] = c[i][j].clone();
            }
            if b + 1 < h {
                r[b * d + i][(b + 1) * d + i] = Rational::ONE;
            }
        }
    }
    r
}

/// Whether `p` (from the constant) is irreducible over ℚ.
fn q_irreducible(p: &[Rational]) -> bool {
    let qx = DenseRing::new(Rationals);
    let factors = qx.factor(&qx.poly(p.to_vec()));
    matches!(factors.factors(), [(f, 1)] if f.degree() == Some(p.len() - 1))
}

/// Whether `a` is already a rational Jordan form: its blocks along the
/// diagonal, the smallest that split it, each [`q_block`] of an
/// irreducible polynomial, in any order.
fn q_in_form(a: &[Vec<Rational>]) -> bool {
    let n = a.len();
    let mut start = 0;
    while start < n {
        let mut end = start;
        loop {
            let reach = (start..=end)
                .flat_map(|i| {
                    (end + 1..n).filter(move |&j| !a[i][j].is_zero() || !a[j][i].is_zero())
                })
                .max();
            match reach {
                Some(j) => end = j,
                None => break,
            }
        }
        let k = end + 1 - start;
        let block: QMatrix = a[start..=end]
            .iter()
            .map(|r| r[start..=end].to_vec())
            .collect();
        let fits = (1..=k).filter(|d| k % d == 0).any(|d| {
            let mut p: Vec<Rational> = (0..d).map(|i| -&block[i][d - 1]).collect();
            p.push(Rational::ONE);
            q_block(&p, k / d) == block && q_irreducible(&p)
        });
        if !fits {
            return false;
        }
        start = end + 1;
    }
    true
}

/// The semisimple part of `a`, by Newton's iteration `S ← S − f(S)
/// f′(S)⁻¹` from `a` for its squarefree `f` (from the constant); `None`
/// where an inverse fails or the iteration does not end.
fn q_semisimple(a: &[Vec<Rational>], f: &[Rational]) -> Option<QMatrix> {
    let n = a.len();
    let derivative: Vec<Rational> = f
        .iter()
        .enumerate()
        .skip(1)
        .map(|(k, c)| c * &Rational::from(i64::try_from(k).expect("small")))
        .collect();
    let mut s = a.to_vec();
    for _ in 0..=n {
        let value = q_eval(f, &s);
        if value.iter().flatten().all(Rational::is_zero) {
            break;
        }
        let step = q_mul(&value, &q_inverse(&q_eval(&derivative, &s))?);
        s = q_sub(&s, &step);
    }
    q_eval(f, &s)
        .iter()
        .flatten()
        .all(Rational::is_zero)
        .then_some(s)
}

/// The tops of the chains of `nil` in the kernel of `ps`, `p(S)` for `p`
/// of degree `d`, with their lengths, longest first: each a vector of the
/// kernel of `p(S)` and `Nʰ` outside the span of the kernel below and of
/// the chains' vectors at its level, `S` applied `d − 1` times included.
fn q_tops(
    ps: &[Vec<Rational>],
    nil: &[Vec<Rational>],
    s: &[Vec<Rational>],
    d: usize,
) -> Option<Vec<(usize, Vec<Rational>)>> {
    let n = ps.len();
    let dimension = q_kernel(ps, n).len();
    // kernels[k]: V_p ∩ ker N^k, up to V_p.
    let mut kernels: Vec<QMatrix> = vec![Vec::new()];
    let mut power = q_identity(n);
    while kernels.last().map(Vec::len) != Some(dimension) {
        if kernels.len() > n {
            return None;
        }
        power = q_mul(&power, nil);
        let stacked: QMatrix = ps.iter().chain(&power).cloned().collect();
        kernels.push(q_kernel(&stacked, n));
    }
    // Chain tops, longest first: a vector of the kernel outside the
    // span of the kernel below and the chains' vectors at its level.
    let mut tops: Vec<(usize, Vec<Rational>)> = Vec::new();
    for size in (1..kernels.len()).rev() {
        let mut span = Vec::new();
        for v in &kernels[size - 1] {
            q_joins(&mut span, v);
        }
        let level = |span: &mut Vec<Vec<Rational>>, h: usize, w: &[Rational]| {
            let mut v = w.to_vec();
            for _ in size..h {
                v = q_apply(nil, &v);
            }
            for _ in 0..d {
                q_joins(span, &v);
                v = q_apply(s, &v);
            }
        };
        for (h, w) in &tops {
            level(&mut span, *h, w);
        }
        for w in &kernels[size] {
            let mut probe = span.clone();
            if q_joins(&mut probe, w) {
                level(&mut span, size, w);
                tops.push((size, w.clone()));
            }
        }
    }
    Some(tops)
}

/// `[P, R]` with `P⁻¹ A P = R` for the rational `a`; `None` where Newton's
/// iteration or an inverse fails, which theory rules out.
fn rational_jordan(a: &[Vec<Rational>]) -> Option<(QMatrix, QMatrix)> {
    let n = a.len();
    let qx = DenseRing::new(Rationals);
    let factors = qx.factor(&qx.poly(q_charpoly(a))).into_parts().1;
    let mut factors: Vec<Vec<Rational>> = factors
        .into_iter()
        .map(|(f, _)| f.coeffs().to_vec())
        .collect();
    // Nonlinear factors first, by degree; then the linear ones by their
    // root, descending.
    factors.sort_by(|f, g| match (f.len(), g.len()) {
        (2, 2) => f[0].cmp(&g[0]),
        (2, _) => Ordering::Greater,
        (_, 2) => Ordering::Less,
        (u, v) => u.cmp(&v),
    });
    // The semisimple part S, by Newton's iteration on the squarefree f.
    let f = factors.iter().fold(qx.poly(vec![Rational::ONE]), |s, g| {
        qx.mul(&s, &qx.poly(g.clone()))
    });
    let s = q_semisimple(a, f.coeffs())?;
    let nil = q_sub(a, &s);
    let mut columns: Vec<Vec<Rational>> = Vec::with_capacity(n);
    let mut r = vec![vec![Rational::ZERO; n]; n];
    let mut start = 0;
    for p in &factors {
        let d = p.len() - 1;
        let ps = q_eval(p, &s);
        let tops = q_tops(&ps, &nil, &s, d)?;
        for (h, w) in tops {
            for b in 1..=h {
                let mut u = w.clone();
                for _ in b..h {
                    u = q_apply(&nil, &u);
                }
                for _ in 0..d {
                    let next = q_apply(&s, &u);
                    columns.push(u);
                    u = next;
                }
            }
            let block = q_block(p, h);
            for (i, row) in block.into_iter().enumerate() {
                for (j, x) in row.into_iter().enumerate() {
                    r[start + i][start + j] = x;
                }
            }
            start += h * d;
        }
    }
    if columns.len() != n {
        return None;
    }
    let p = (0..n)
        .map(|i| columns.iter().map(|v| v[i].clone()).collect())
        .collect();
    Some((p, r))
}

/// `RationalJordanForm[A]`.
pub(in crate::native) fn rational_jordan_form(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [_] = e.args() else {
        return Ok(None);
    };
    let Some((a, n)) = square_matrix(e) else {
        return Ok(None);
    };
    let Some(a) = rows(&a).and_then(|r| {
        r.iter()
            .map(|row| row.iter().map(Expr::to_rational).collect())
            .collect::<Option<QMatrix>>()
    }) else {
        return Ok(None);
    };
    let found = if q_in_form(&a) {
        Some((q_identity(n), a))
    } else {
        rational_jordan(&a)
    };
    let Some((p, r)) = found else {
        return Ok(None);
    };
    let expr = |m: QMatrix| {
        matrix(
            m.into_iter()
                .map(|row| row.into_iter().map(Expr::from).collect())
                .collect(),
        )
    };
    let result = list(vec![expr(p), expr(r)]);
    ev.step("rat_jordan", e, &result);
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
    let product = ev.eval(&times(&[p.clone(), f.clone(), inverse]))?;
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

pub(in crate::native) mod forms;
pub(in crate::native) mod iterative;
pub(in crate::native) mod least_squares;
pub(in crate::native) mod norms;
pub(in crate::native) mod unitary;

#[cfg(test)]
mod tests;
