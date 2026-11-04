//! Integer matrices (CAS-DOC-42): Xcas's `ihermite`, `ismith` and `lll`.
//!
//! - `HermiteDecomposition[A]` (§2, CAS-445): `{U, H}`, `U` unimodular and
//!   `H = U A` in Hermite normal form: row echelon, each pivot `p`
//!   positive, the entries above it in `(−p/2, p/2]` as giac's (unique as
//!   `[0, p)` would be), zero rows last. Row operations of determinant ±1
//!   on `A` and `U` at once: below a pivot, two rows `a, b` become
//!   `g = s a + t b` and 0 by the extended gcd; then the pivot made
//!   positive and the rows above reduced modulo it. giac's
//!   `ihermite([[1,2,3],[4,5,6],[7,8,10]])` has
//!   `H = [[1,-1,0],[0,3,0],[0,0,1]]`.
//! - `SmithDecomposition[A]` (§3, CAS-446): `{U, S, V}`, `U` and `V`
//!   unimodular, `S = U A V` diagonal with `d₁ | d₂ | …`, each `dᵢ ≥ 0`:
//!   at each corner, the least entry of the block moved there, its column
//!   and row cleared by gcd steps until both are, and a row whose entry the
//!   corner does not divide added to the corner's row, until it divides
//!   them all.
//! - `LatticeReduce[B]` (§4, CAS-447), giac's `lll`'s first matrix: the
//!   LLL-reduced basis of the lattice of `B`'s rows, `δ = 3/4`, in exact
//!   rationals: each row size-reduced against those before it, `|μ_kj| ≤
//!   1/2`, then swapped back while Lovász's condition `B_k ≥ (δ − μ²) B_(k−1)`
//!   fails, the Gram–Schmidt data computed again after a swap. Rows that
//!   depend on others go: the lattice's basis is then the Hermite form's
//!   rows that are not 0, reduced in turn.
//!
//! - `SmithDecomposition[A]` over ℚ\[x\] (CAS-DOC-51 §6, CAS-500), Xcas's
//!   `smith`: entries polynomials over ℚ in one variable, the same
//!   elimination with the Euclidean division of ℚ\[x\], each `dᵢ` monic;
//!   giac's `smith([[x,1],[0,x]])` has `D = [[1,0],[0,-x²]]`, ours `x²`.
//!
//! The matrices are rectangular, of integers (or for Smith's form of
//! polynomials in one variable); anything else leaves the call.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of matrix algebra: U A V = S, rows r and k, s a + t b = g"
)]

use cas_num::algebra::{EuclideanDomain, Rationals, Ring};
use cas_num::{Integer, Natural, Rational};
use cas_poly::{Dense, DenseRing};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

type Matrix = Vec<Vec<Integer>>;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `e` as a rectangular matrix of integers, one row at least.
fn integers(e: &Expr) -> Option<Matrix> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) || e.args().is_empty() {
        return None;
    }
    let rows: Matrix = e
        .args()
        .iter()
        .map(|row| {
            if row.head() != Some(&Expr::from(Symbol::LIST)) || row.args().is_empty() {
                return None;
            }
            row.args().iter().map(Expr::to_integer).collect()
        })
        .collect::<Option<_>>()?;
    let width = rows[0].len();
    rows.iter().all(|r| r.len() == width).then_some(rows)
}

fn matrix(m: &Matrix) -> Expr {
    list(
        m.iter()
            .map(|row| list(row.iter().cloned().map(Expr::from).collect()))
            .collect(),
    )
}

fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| Integer::from(i64::from(i == j))).collect())
        .collect()
}

/// Rows `r` and `k` of `m` replaced by `s r + t k` and `u r + v k`.
fn combine_rows(m: &mut Matrix, r: usize, k: usize, [s, t, u, v]: [&Integer; 4]) {
    for j in 0..m[r].len() {
        let (a, b) = (m[r][j].clone(), m[k][j].clone());
        m[r][j] = &(s * &a) + &(t * &b);
        m[k][j] = &(u * &a) + &(v * &b);
    }
}

/// The unimodular `[s, t, −b/g, a/g]` taking `(a, b)` to `(g, 0)`. When
/// `a` divides `b`, `s = ±1` and `t = 0`: `a`'s row stays but for its sign.
/// GMP's cofactors for `|a| = |b|`, `s = 0` and `t = sgn b`, swapped the
/// two, and Smith's clearing of a row and a column went round (CAS-455).
fn step(a: &Integer, b: &Integer) -> [Integer; 4] {
    if !a.is_zero() && b.div_rem_euclid(a).1.is_zero() {
        let g = a.abs();
        return [
            a.div_exact(&g),
            Integer::ZERO,
            -b.div_exact(&g),
            a.div_exact(&g),
        ];
    }
    let (g, s, t) = a.gcd_ext(b);
    [s, t, -b.div_exact(&g), a.div_exact(&g)]
}

/// The Hermite form of `h` in place, `u` along: rows from `top` on, columns
/// from `left` on.
fn hermite(h: &mut Matrix, u: &mut Matrix) {
    let (m, n) = (h.len(), h[0].len());
    let mut r = 0;
    for col in 0..n {
        if r == m {
            break;
        }
        for k in r + 1..m {
            if h[k][col].is_zero() {
                continue;
            }
            let ops = step(&h[r][col], &h[k][col]);
            let refs = [&ops[0], &ops[1], &ops[2], &ops[3]];
            combine_rows(h, r, k, refs);
            combine_rows(u, r, k, refs);
        }
        if h[r][col].is_zero() {
            continue;
        }
        if h[r][col].is_negative() {
            for row in [&mut *h, &mut *u] {
                for x in &mut row[r] {
                    *x = -x.clone();
                }
            }
        }
        let pivot = h[r][col].clone();
        for k in 0..r {
            // The symmetric remainder: q with h − q p in (−p/2, p/2].
            let (mut q, rest) = h[k][col].div_rem_euclid(&pivot);
            if (&rest + &rest) > pivot {
                q = &q + &Integer::ONE;
            }
            if q.is_zero() {
                continue;
            }
            for row in [&mut *h, &mut *u] {
                for j in 0..row[k].len() {
                    let v = &row[k][j] - &(&q * &row[r][j]);
                    row[k][j] = v;
                }
            }
        }
        r += 1;
    }
}

/// `HermiteDecomposition[A]`.
pub(super) fn hermite_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let Some(mut h) = integers(a) else {
        return Ok(None);
    };
    let mut u = identity(h.len());
    hermite(&mut h, &mut u);
    let result = list(vec![matrix(&u), matrix(&h)]);
    ev.step("HermiteDecomposition", e, &result);
    Ok(Some(result))
}

/// A Euclidean domain for Smith's form: ℤ, or ℚ\[x\] (CAS-DOC-51 §6).
trait Domain {
    type E: Clone;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn is_zero(&self, a: &Self::E) -> bool;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    /// The unimodular `[s, t, u, v]` taking `(a, b)` to `(g, 0)`.
    fn step(&self, a: &Self::E, b: &Self::E) -> [Self::E; 4];
    /// Whether `a ≠ 0` divides `b`.
    fn divides(&self, a: &Self::E, b: &Self::E) -> bool;
    /// `a`'s size in the Euclidean function: `|a|`, or the degree.
    fn size(&self, a: &Self::E) -> Natural;
    /// The unit `u` making `a u` canonical: `a`'s sign, or `1/lc(a)`.
    fn unit(&self, a: &Self::E) -> Self::E;
}

/// The integers, with GMP's cofactors as [`step`] fixes them.
struct Z;

impl Domain for Z {
    type E = Integer;
    fn zero(&self) -> Integer {
        Integer::ZERO
    }
    fn one(&self) -> Integer {
        Integer::ONE
    }
    fn is_zero(&self, a: &Integer) -> bool {
        a.is_zero()
    }
    fn add(&self, a: &Integer, b: &Integer) -> Integer {
        a + b
    }
    fn mul(&self, a: &Integer, b: &Integer) -> Integer {
        a * b
    }
    fn step(&self, a: &Integer, b: &Integer) -> [Integer; 4] {
        step(a, b)
    }
    fn divides(&self, a: &Integer, b: &Integer) -> bool {
        b.div_rem_euclid(a).1.is_zero()
    }
    fn size(&self, a: &Integer) -> Natural {
        a.unsigned_abs()
    }
    fn unit(&self, a: &Integer) -> Integer {
        if a.is_negative() {
            -Integer::ONE
        } else {
            Integer::ONE
        }
    }
}

/// ℚ\[x\]: monic gcds by the extended Euclidean algorithm.
struct Qx(DenseRing<Rationals>);

impl Domain for Qx {
    type E = Dense<Rational>;
    fn zero(&self) -> Dense<Rational> {
        self.0.zero()
    }
    fn one(&self) -> Dense<Rational> {
        self.0.one()
    }
    fn is_zero(&self, a: &Dense<Rational>) -> bool {
        a.is_zero()
    }
    fn add(&self, a: &Dense<Rational>, b: &Dense<Rational>) -> Dense<Rational> {
        self.0.add(a, b)
    }
    fn mul(&self, a: &Dense<Rational>, b: &Dense<Rational>) -> Dense<Rational> {
        self.0.mul(a, b)
    }
    fn step(&self, a: &Dense<Rational>, b: &Dense<Rational>) -> [Dense<Rational>; 4] {
        let ring = &self.0;
        let quo = |x: &Dense<Rational>, y: &Dense<Rational>| ring.div_rem(x, y).0;
        if !a.is_zero() && self.divides(a, b) {
            // `a`'s row stays but for a unit, as in ℤ.
            let (g, u) = ring.normal_form(a);
            return [u, ring.zero(), ring.neg(&quo(b, &g)), quo(a, &g)];
        }
        let (g, s, t) = EuclideanDomain::gcd_ext(ring, a, b);
        [s, t, ring.neg(&quo(b, &g)), quo(a, &g)]
    }
    fn divides(&self, a: &Dense<Rational>, b: &Dense<Rational>) -> bool {
        self.0.div_rem(b, a).1.is_zero()
    }
    fn size(&self, a: &Dense<Rational>) -> Natural {
        Natural::from(a.degree().map_or(0, |d| d as u64))
    }
    fn unit(&self, a: &Dense<Rational>) -> Dense<Rational> {
        self.0.normal_form(a).1
    }
}

/// Rows `r` and `k` of `m` replaced by `s r + t k` and `u r + v k`.
fn combine_rows_in<D: Domain>(d: &D, m: &mut [Vec<D::E>], r: usize, k: usize, ops: &[D::E; 4]) {
    let [s, t, u, v] = ops;
    for j in 0..m[r].len() {
        let (a, b) = (m[r][j].clone(), m[k][j].clone());
        m[r][j] = d.add(&d.mul(s, &a), &d.mul(t, &b));
        m[k][j] = d.add(&d.mul(u, &a), &d.mul(v, &b));
    }
}

/// Columns `c` and `k` of `m` replaced by `s c + t k` and `u c + v k`.
fn combine_columns_in<D: Domain>(d: &D, m: &mut [Vec<D::E>], c: usize, k: usize, ops: &[D::E; 4]) {
    let [s, t, u, v] = ops;
    for row in m.iter_mut() {
        let (a, b) = (row[c].clone(), row[k].clone());
        row[c] = d.add(&d.mul(s, &a), &d.mul(t, &b));
        row[k] = d.add(&d.mul(u, &a), &d.mul(v, &b));
    }
}

/// The Smith form of `s` in place, `u` and `v` along.
fn smith<D: Domain>(d: &D, s: &mut [Vec<D::E>], u: &mut [Vec<D::E>], v: &mut [Vec<D::E>]) {
    let (m, n) = (s.len(), s[0].len());
    for t in 0..m.min(n) {
        loop {
            // The least nonzero entry of the block to the corner.
            let least = (t..m)
                .flat_map(|i| (t..n).map(move |j| (i, j)))
                .filter(|&(i, j)| !d.is_zero(&s[i][j]))
                .min_by_key(|&(i, j)| d.size(&s[i][j]));
            let Some((i, j)) = least else {
                return;
            };
            s.swap(t, i);
            u.swap(t, i);
            for row in s.iter_mut().chain(v.iter_mut()) {
                row.swap(t, j);
            }
            // The column below and the row right of the corner cleared.
            let mut clear = false;
            while !clear {
                for k in t + 1..m {
                    if !d.is_zero(&s[k][t]) {
                        let ops = d.step(&s[t][t], &s[k][t]);
                        combine_rows_in(d, s, t, k, &ops);
                        combine_rows_in(d, u, t, k, &ops);
                    }
                }
                for k in t + 1..n {
                    if !d.is_zero(&s[t][k]) {
                        let ops = d.step(&s[t][t], &s[t][k]);
                        combine_columns_in(d, s, t, k, &ops);
                        combine_columns_in(d, v, t, k, &ops);
                    }
                }
                clear = (t + 1..m).all(|k| d.is_zero(&s[k][t]));
            }
            // A row whose entry the corner does not divide goes into its row.
            let corner = s[t][t].clone();
            let rest = (t + 1..m).find(|&k| (t + 1..n).any(|j| !d.divides(&corner, &s[k][j])));
            match rest {
                Some(k) => {
                    let ops = [d.one(), d.one(), d.zero(), d.one()];
                    combine_rows_in(d, s, t, k, &ops);
                    combine_rows_in(d, u, t, k, &ops);
                }
                None => break,
            }
        }
        let unit = d.unit(&s[t][t]);
        for row in [&mut *s, &mut *u] {
            for x in &mut row[t] {
                *x = d.mul(x, &unit);
            }
        }
    }
}

/// The only variable of the entries of `a`, none of them a head; `None`
/// for two or more, or none.
fn variable(a: &Expr) -> Option<Expr> {
    fn walk(e: &Expr, found: &mut Vec<Symbol>) {
        if let Some(s) = e.as_symbol() {
            if !found.contains(&s) {
                found.push(s);
            }
        }
        for x in e.args() {
            walk(x, found);
        }
    }
    let mut found = Vec::new();
    walk(a, &mut found);
    match found.as_slice() {
        [x] => Some(Expr::from(*x)),
        _ => None,
    }
}

/// `a` as a matrix of polynomials over ℚ in `x`, rectangular, from their
/// coefficient lists.
fn polynomials(
    ev: &mut Eval<'_>,
    ring: &DenseRing<Rationals>,
    a: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Vec<Dense<Rational>>>>, EvalError> {
    if a.head() != Some(&Expr::from(Symbol::LIST)) || a.args().is_empty() {
        return Ok(None);
    }
    let width = a.args()[0].args().len();
    let mut rows = Vec::with_capacity(a.args().len());
    for row in a.args() {
        if row.head() != Some(&Expr::from(Symbol::LIST)) || row.args().len() != width || width == 0
        {
            return Ok(None);
        }
        let mut entries = Vec::with_capacity(width);
        for p in row.args() {
            let c = ev.eval(&Expr::call(
                Symbol::COEFFICIENT_LIST,
                [p.clone(), x.clone()],
            ))?;
            if c.head() != Some(&Expr::from(Symbol::LIST)) {
                return Ok(None);
            }
            let Some(mut c) = c
                .args()
                .iter()
                .map(Expr::to_rational)
                .collect::<Option<Vec<_>>>()
            else {
                return Ok(None);
            };
            c.reverse();
            entries.push(ring.poly(c));
        }
        rows.push(entries);
    }
    Ok(Some(rows))
}

/// The matrix of the polynomials, in `x`.
fn polynomial_matrix(m: &[Vec<Dense<Rational>>], x: &Expr) -> Expr {
    list(
        m.iter()
            .map(|row| {
                list(
                    row.iter()
                        .map(|p| {
                            Expr::call(
                                Symbol::PLUS,
                                p.coeffs()
                                    .iter()
                                    .enumerate()
                                    .map(|(k, c)| {
                                        Expr::from(c.clone())
                                            * x.clone().pow(i64::try_from(k).expect("small"))
                                    })
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

/// `SmithDecomposition[A]`: over ℤ, or over ℚ\[x\] for entries that are
/// polynomials in one variable `x` (Xcas's `smith`), the diagonal monic.
pub(super) fn smith_decomposition(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    if let Some(mut s) = integers(a) {
        let mut u = identity(s.len());
        let mut v = identity(s[0].len());
        smith(&Z, &mut s, &mut u, &mut v);
        let result = list(vec![matrix(&u), matrix(&s), matrix(&v)]);
        ev.step("SmithDecomposition", e, &result);
        return Ok(Some(result));
    }
    let Some(x) = variable(a) else {
        return Ok(None);
    };
    let d = Qx(DenseRing::new(Rationals));
    let Some(mut s) = polynomials(ev, &d.0, a, &x)? else {
        return Ok(None);
    };
    let ident = |n: usize| -> Vec<Vec<Dense<Rational>>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { d.one() } else { d.zero() })
                    .collect()
            })
            .collect()
    };
    let mut u = ident(s.len());
    let mut v = ident(s[0].len());
    smith(&d, &mut s, &mut u, &mut v);
    let result = list(vec![
        polynomial_matrix(&u, &x),
        polynomial_matrix(&s, &x),
        polynomial_matrix(&v, &x),
    ]);
    let result = ev.eval(&result)?;
    ev.step("smith", e, &result);
    Ok(Some(result))
}

/// `⟨a, b⟩` over the rationals.
fn inner(a: &[Rational], b: &[Rational]) -> Rational {
    a.iter()
        .zip(b)
        .fold(Rational::ZERO, |sum, (x, y)| &sum + &(x * y))
}

/// The Gram–Schmidt coefficients `μ` and squared norms `B` of the rows.
fn gram_schmidt(b: &Matrix) -> (Vec<Vec<Rational>>, Vec<Rational>) {
    let n = b.len();
    let rows: Vec<Vec<Rational>> = b
        .iter()
        .map(|r| r.iter().cloned().map(Rational::from).collect())
        .collect();
    let mut star: Vec<Vec<Rational>> = Vec::with_capacity(n);
    let mut mu = vec![vec![Rational::ZERO; n]; n];
    let mut norms = Vec::with_capacity(n);
    for i in 0..n {
        let mut v = rows[i].clone();
        for j in 0..i {
            mu[i][j] = &inner(&rows[i], &star[j]) / &norms[j];
            for (x, y) in v.iter_mut().zip(&star[j]) {
                *x = &*x - &(&mu[i][j] * y);
            }
        }
        norms.push(inner(&v, &v));
        star.push(v);
    }
    (mu, norms)
}

/// The LLL reduction of independent rows.
fn lll(mut b: Matrix) -> Matrix {
    let n = b.len();
    let delta = Rational::new(Integer::from(3), Integer::from(4));
    let (mut mu, mut norms) = gram_schmidt(&b);
    let mut k = 1;
    while k < n {
        for j in (0..k).rev() {
            let q = mu[k][j].round();
            if q.is_zero() {
                continue;
            }
            for c in 0..b[k].len() {
                let v = &b[k][c] - &(&q * &b[j][c]);
                b[k][c] = v;
            }
            let qr = Rational::from(q);
            let (low, high) = mu.split_at_mut(k);
            let row = &mut high[0];
            for (x, y) in row.iter_mut().zip(&low[j]).take(j) {
                *x = &*x - &(&qr * y);
            }
            row[j] = &row[j] - &qr;
        }
        let bound = &(&delta - &(&mu[k][k - 1] * &mu[k][k - 1])) * &norms[k - 1];
        if norms[k] >= bound {
            k += 1;
        } else {
            b.swap(k, k - 1);
            (mu, norms) = gram_schmidt(&b);
            k = (k - 1).max(1);
        }
    }
    b
}

/// `LatticeReduce[B]`.
pub(super) fn lattice_reduce(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [b] = e.args() else {
        return Ok(None);
    };
    let Some(rows) = integers(b) else {
        return Ok(None);
    };
    let mut h = rows.clone();
    hermite(&mut h, &mut identity(rows.len()));
    h.retain(|r| r.iter().any(|x| !x.is_zero()));
    let basis = if h.len() == rows.len() { rows } else { h };
    let result = if basis.is_empty() {
        list(Vec::new())
    } else {
        matrix(&lll(basis))
    };
    ev.step("LatticeReduce", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
