//! `Eigenvectors[m]`, Xcas's `eigenvects` and `egv` (CAS-DOC-12 §5,
//! CAS-222): a matrix whose columns are eigenvectors, in the order of
//! `Eigenvalues`.
//!
//! The characteristic polynomial is factored over ℚ with the parameters,
//! and for each irreducible factor `q` of multiplicity `k`:
//!
//! - `q` of degree 1, its root `λ` rational in the parameters: the null
//!   space of `m − λ I` by fraction-free elimination ([`kernel`]).
//! - `q` of degree 2 or more, the entries rational numbers: the null space
//!   of `m − t I` over the field ℚ\[t\]/(q(t)), by Gauss–Jordan with
//!   inverses from the extended gcd ([`Extension`]), then `t` replaced by
//!   each root of `q`: one computation for all of them.
//!
//! Each vector has 1, or a positive integer, at its free column, the pivot
//! columns' entries from their rows, its rational coefficients cleared and
//! made primitive: giac's `egv([[1,2],[2,1]])` is `[[1,-1],[1,1]]` and
//! `eigenvects([[1,2],[3,4]])` has the column `(√33 − 3, 6)`. An
//! eigenvalue with fewer independent vectors than its multiplicity: not
//! diagonalizable, the call stays with a warning, where giac errs. A
//! numeric matrix, or a factor of degree 2 or more with parameters, stays
//! too.

use cas_num::algebra::{EuclideanDomain, Integers, Rationals, Ring};
use cas_num::{Integer, Rational};
use cas_poly::{Cancel, Dense, DenseRing, Sparse, SparseRing};

use super::super::polynomial::{Numbers, convert, factored};
use super::elimination::{Row, gauss_jordan, primitive, scaled_rows};
use super::{
    Native, char_poly, entries, fresh, integer_fractions, list, numeric, square, stepped,
    too_large, warn,
};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::Fraction;

/// An eigenvalue and a basis of its eigenvectors.
type Eigen = (Expr, Vec<Vec<Expr>>);

/// `Eigenvectors[m]`.
pub(in crate::native) fn eigenvectors(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((m, n)) = square(ev, e) else {
        return Ok(None);
    };
    if numeric(&m) {
        warn(
            ev,
            "the eigenvectors of a numeric matrix are not computed here",
            e,
        );
        return Ok(None);
    }
    let Some((poly, _)) = char_poly(ev, e, &m, n)? else {
        return Ok(None);
    };
    let x = fresh();
    let Some(values) = super::super::solve::roots_counted(ev, e, &poly, &x)? else {
        return Ok(None);
    };
    let Some(p) = convert(ev, e, std::slice::from_ref(&poly), false)? else {
        return Ok(None);
    };
    let (Some(kx), Fraction::Integers { num, den }) = (p.variable(&x), &p.fractions()[0]) else {
        warn(
            ev,
            "the eigenvectors of a matrix with complex entries are not computed here",
            e,
        );
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, p.generators().len());
    let Some(powers) = factored(ev, e, &ring, num, den)? else {
        return Ok(None);
    };
    // Each factor's roots and their vectors.
    let mut groups: Vec<Eigen> = Vec::new();
    for (f, multiplicity) in &powers.factors {
        let degree = ring.degree(f, kx).unwrap_or(0);
        if degree == 0 || *multiplicity < 0 {
            continue;
        }
        let Some(roots) = super::super::solve::complex_roots(ev, e, &p.polynomial_expr(f), &x)?
        else {
            return Ok(None);
        };
        let Some(found) = factor_vectors(ev, e, &m, n, f, kx, &ring, roots)? else {
            return Ok(None);
        };
        for (lambda, vectors) in found {
            if u64::try_from(vectors.len()).ok() != u64::try_from(*multiplicity).ok() {
                warn(ev, "the matrix is not diagonalizable", e);
                return Ok(None);
            }
            groups.push((lambda, vectors));
        }
    }
    // The columns in the order of the eigenvalues.
    let mut columns: Vec<Vec<Expr>> = Vec::with_capacity(n);
    let mut seen: Vec<&Expr> = Vec::new();
    for value in &values {
        if seen.contains(&value) {
            continue;
        }
        seen.push(value);
        let Some((_, vectors)) = groups.iter().find(|(lambda, _)| lambda == value) else {
            return Ok(None);
        };
        columns.extend(vectors.iter().cloned());
    }
    if columns.len() != n {
        return Ok(None);
    }
    let rows = (0..n).map(|i| list(columns.iter().map(|c| c[i].clone())));
    stepped(ev, "eigenvects", e, list(rows.collect::<Vec<_>>()))
}

/// The roots of the irreducible factor `f` of the characteristic polynomial
/// of the `n × n` matrix `m`, in the variable `kx` of `ring`, each with a
/// basis of its eigenvectors; `None` after a warning.
#[allow(
    clippy::too_many_arguments,
    reason = "the matrix, the factor and its ring, and the factor's roots"
)]
fn factor_vectors(
    ev: &mut Eval<'_>,
    e: &Expr,
    m: &Expr,
    n: usize,
    f: &Sparse<Integer>,
    kx: usize,
    ring: &SparseRing<Integers>,
    roots: Vec<Expr>,
) -> Result<Option<Vec<Eigen>>, EvalError> {
    if ring.degree(f, kx) == Some(1) {
        let Some(lambda) = roots.into_iter().next() else {
            return Ok(None);
        };
        let shifted = ev.eval(&shift(m, &lambda))?;
        let Some(basis) = kernel(ev, e, &shifted, n)? else {
            return Ok(None);
        };
        return Ok(Some(vec![(lambda, basis)]));
    }
    let alone = f
        .degrees()
        .iter()
        .enumerate()
        .all(|(j, &d)| j == kx || d == 0);
    let rationals: Option<Vec<Rational>> = entries(m).iter().map(Expr::to_rational).collect();
    let (true, Some(rationals)) = (alone, rationals) else {
        warn(
            ev,
            "eigenvectors for a factor of degree 2 or more with parameters are not computed here",
            e,
        );
        return Ok(None);
    };
    let dense = ring.to_dense(f, kx).expect("a polynomial in x alone");
    let qx = DenseRing::new(Rationals);
    let q = qx.poly(
        dense
            .coeffs()
            .iter()
            .map(|c| Rational::from(c.clone()))
            .collect(),
    );
    let basis = Extension { qx, q }.kernel(&rationals, n);
    let mut found = Vec::new();
    for lambda in roots {
        let mut vectors = Vec::with_capacity(basis.len());
        for v in &basis {
            let entries = v
                .iter()
                .map(|c| ev.eval(&at(c, &lambda)))
                .collect::<Result<Vec<_>, _>>()?;
            vectors.push(entries);
        }
        found.push((lambda, vectors));
    }
    Ok(Some(found))
}

/// `m − λ I`, built.
fn shift(m: &Expr, lambda: &Expr) -> Expr {
    list(m.args().iter().enumerate().map(|(i, row)| {
        list(row.args().iter().enumerate().map(|(j, a)| {
            if i == j {
                a.clone() - lambda.clone()
            } else {
                a.clone()
            }
        }))
    }))
}

/// `Σ cₖ λᵏ` for the coefficients `c` of an element of the extension.
fn at(c: &[Integer], lambda: &Expr) -> Expr {
    Expr::call(
        crate::symbol::Symbol::PLUS,
        c.iter()
            .enumerate()
            .filter(|(_, ck)| !ck.is_zero())
            .map(|(k, ck)| {
                let power = i64::try_from(k).expect("a degree in memory");
                Expr::from(ck.clone()) * lambda.clone().pow(power)
            })
            .collect::<Vec<_>>(),
    )
}

/// A basis of the null space of the `n × n` matrix `b`, exact, each vector
/// with a positive entry at its free column, primitive over the
/// polynomials of the generators; `None` after a warning when eliminating
/// is too large.
fn kernel(
    ev: &mut Eval<'_>,
    e: &Expr,
    b: &Expr,
    n: usize,
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    let Some(p) = convert(ev, e, &entries(b), true)? else {
        return Ok(None);
    };
    let ring = match &p.fractions()[0] {
        Fraction::Integers { .. } => SparseRing::new(Integers, p.generators().len()),
        Fraction::GaussianIntegers { .. } => {
            warn(
                ev,
                "the eigenvectors of a matrix with complex entries are not computed here",
                e,
            );
            return Ok(None);
        }
    };
    let (mut rows, _) = scaled_rows(&ring, &integer_fractions(&p), n);
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let Some(pivots) = gauss_jordan(&ring, &mut rows, n, cancel, ev.context().limits().terms)?
    else {
        too_large(ev, e);
        return Ok(None);
    };
    let basis = fraction_free_kernel(&ring, &rows, &pivots, n, cancel)?;
    Ok(Some(
        basis
            .into_iter()
            .map(|v| v.iter().map(|x| p.polynomial_expr(x)).collect())
            .collect(),
    ))
}

/// The null space of rows eliminated by Gauss–Jordan, their pivots
/// `pivots`, without fractions: for each free column `f`, `L` at `f` and
/// `−r[f] L / r[p]` at each pivot column `p` of a row `r`, `L` the lcm of
/// the pivots, made primitive.
fn fraction_free_kernel<R: Numbers>(
    ring: &SparseRing<R>,
    rows: &[Row<R::Elem>],
    pivots: &[usize],
    columns: usize,
    cancel: &dyn Cancel,
) -> Result<Vec<Row<R::Elem>>, EvalError> {
    let lcm = pivots
        .iter()
        .enumerate()
        .fold(ring.one(), |l, (i, &p)| ring.lcm(&l, &rows[i][p]));
    let mut basis = Vec::new();
    for f in (0..columns).filter(|j| !pivots.contains(j)) {
        let mut v: Row<R::Elem> = vec![ring.zero(); columns];
        v[f] = lcm.clone();
        for (i, &p) in pivots.iter().enumerate() {
            let scale = ring
                .div_exact(&lcm, &rows[i][p])
                .expect("a divisor of the lcm");
            v[p] = ring.neg(&ring.mul(&rows[i][f], &scale));
        }
        basis.push(primitive(ring, v, cancel)?);
    }
    Ok(basis)
}

/// The field ℚ\[t\]/(q(t)) of an irreducible `q`.
struct Extension {
    qx: DenseRing<Rationals>,
    q: Dense<Rational>,
}

impl Extension {
    fn reduce(&self, a: &Dense<Rational>) -> Dense<Rational> {
        self.qx.div_rem(a, &self.q).1
    }

    fn mul(&self, a: &Dense<Rational>, b: &Dense<Rational>) -> Dense<Rational> {
        self.reduce(&self.qx.mul(a, b))
    }

    /// The inverse of a nonzero `a`: `s` of `s a + u q = 1`.
    fn inv(&self, a: &Dense<Rational>) -> Dense<Rational> {
        let (_, s, _) = EuclideanDomain::gcd_ext(&self.qx, a, &self.q);
        self.reduce(&s)
    }

    /// A basis of the null space of `m − t I` for the `n × n` rational
    /// `m`, row by row: each vector 1 at its free column, the negated
    /// entries of its pivot rows, then over ℤ, primitive: the coefficients
    /// in `t` of each entry, from the constant.
    fn kernel(&self, m: &[Rational], n: usize) -> Vec<Vec<Vec<Integer>>> {
        let qx = &self.qx;
        let mut b: Vec<Vec<Dense<Rational>>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let c = m[i * n + j].clone();
                        if i == j {
                            qx.poly(vec![c, -Rational::ONE])
                        } else {
                            qx.poly(vec![c])
                        }
                    })
                    .collect()
            })
            .collect();
        let mut pivots = Vec::new();
        for j in 0..n {
            let rank = pivots.len();
            let Some(r) = (rank..n).find(|&i| !b[i][j].is_zero()) else {
                continue;
            };
            b.swap(rank, r);
            let inverse = self.inv(&b[rank][j]);
            for x in &mut b[rank] {
                *x = self.mul(x, &inverse);
            }
            let pivot_row = b[rank].clone();
            for (i, row) in b.iter_mut().enumerate() {
                if i == rank || row[j].is_zero() {
                    continue;
                }
                let c = row[j].clone();
                for (x, y) in row.iter_mut().zip(&pivot_row) {
                    *x = qx.sub(x, &self.mul(&c, y));
                }
            }
            pivots.push(j);
        }
        (0..n)
            .filter(|j| !pivots.contains(j))
            .map(|f| {
                let mut v = vec![qx.poly(vec![Rational::ZERO]); n];
                v[f] = qx.poly(vec![Rational::ONE]);
                for (i, &p) in pivots.iter().enumerate() {
                    v[p] = qx.neg(&b[i][f]);
                }
                integral(&v)
            })
            .collect()
    }
}

/// The entries of `v`, elements of the extension, times the lcm of their
/// coefficients' denominators and over the gcd of the numerators: their
/// coefficients from the constant.
fn integral(v: &[Dense<Rational>]) -> Vec<Vec<Integer>> {
    let coeffs = || v.iter().flat_map(|x| x.coeffs().iter());
    let l = coeffs().fold(Integer::ONE, |l, c| l.lcm(c.denom()));
    let scaled: Vec<Vec<Integer>> = v
        .iter()
        .map(|x| {
            x.coeffs()
                .iter()
                .map(|c| (c * &Rational::from(l.clone())).numer().clone())
                .collect()
        })
        .collect();
    let g = scaled.iter().flatten().fold(Integer::ZERO, |g, c| g.gcd(c));
    if g.is_zero() || g == Integer::ONE {
        return scaled;
    }
    scaled
        .into_iter()
        .map(|x| x.into_iter().map(|c| c.div_exact(&g)).collect())
        .collect()
}
