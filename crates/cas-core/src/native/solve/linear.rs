//! Systems of linear equations (CAS-DOC-11 §5, CAS-212): `Solve[{eqns},
//! {x, y, …}]`, giac's `solve([…],[…])`, by Gauss–Jordan elimination over
//! the rational functions of the parameters.
//!
//! - Each equation is `l − r = 0`, converted with the others (CAS-DOC-8 §2).
//!   The system is linear when no generator but the unknowns holds one, no
//!   denominator holds one, and no term of a numerator has two unknowns or
//!   one squared; each numerator is then a row `Σ aⱼ xⱼ = b`, its
//!   denominator, free of the unknowns, left out.
//! - The elimination is fraction-free, over ℤ or ℤ\[i\] in the parameters:
//!   for each unknown in turn a row with a nonzero coefficient, of the
//!   fewest terms, is the pivot, and each other row `r` with a coefficient
//!   `c` there becomes `p r − c (pivot row)`, over the gcd of its entries.
//!   The parameters are generic: a pivot is a nonzero polynomial, as giac's
//!   `(b+1)/(a+1)`.
//! - A row left with no unknown and a nonzero constant: no solution, `{}`.
//!   Otherwise one tuple, in the order of the unknowns: a pivot's unknown
//!   is `b/p − Σ (c/p) y` over the free unknowns `y` of its row, each
//!   fraction in lowest terms; a free unknown is itself, as giac's
//!   `[[-y+1,y]]`.

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::super::polynomial::{Dense, Numbers, cancelled, convert, within_dense};
use super::sides;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;

/// A row: the coefficients of the unknowns, then the right side.
type Row<E> = Vec<Sparse<E>>;

/// The solutions of the system `equations` in the distinct symbols
/// `unknowns`, as `{{…}}` or `{}`; `None` when it is not linear, or after
/// a warning when it is too large.
pub(super) fn system(
    ev: &mut Eval<'_>,
    e: &Expr,
    equations: &[Expr],
    unknowns: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let mut fs = Vec::with_capacity(equations.len());
    for eq in equations {
        match eq.as_symbol() {
            Some(Symbol::TRUE) => continue,
            Some(Symbol::FALSE) => return Ok(Some(none())),
            _ => {}
        }
        fs.push(match sides(eq) {
            Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
            None => eq.clone(),
        });
    }
    if fs.is_empty() {
        return Ok(Some(tuples(unknowns.to_vec())));
    }
    let Some(p) = convert(ev, e, &fs, true)? else {
        return Ok(None);
    };
    let vars: Vec<Option<usize>> = unknowns.iter().map(|x| p.variable(x)).collect();
    let inside = p
        .generators()
        .iter()
        .enumerate()
        .any(|(j, g)| !vars.contains(&Some(j)) && unknowns.iter().any(|x| g.contains(x)));
    if inside || within_dense(ev, e, &p, Dense::All).is_none() {
        return Ok(None);
    }
    let n = p.generators().len();
    let problem = Problem {
        p: &p,
        vars: &vars,
        unknowns,
    };
    match &p.fractions()[0] {
        Fraction::Integers { .. } => {
            let fractions = p.fractions().iter().map(|f| match f {
                Fraction::Integers { num, den } => (num, den),
                Fraction::GaussianIntegers { .. } => unreachable!("one ring"),
            });
            problem.solve(ev, e, &SparseRing::new(Integers, n), fractions)
        }
        Fraction::GaussianIntegers { .. } => {
            let fractions = p.fractions().iter().map(|f| match f {
                Fraction::GaussianIntegers { num, den } => (num, den),
                Fraction::Integers { .. } => unreachable!("one ring"),
            });
            problem.solve(ev, e, &SparseRing::new(GaussianIntegers, n), fractions)
        }
    }
}

/// `{}`, no solution.
fn none() -> Expr {
    Expr::call(Symbol::LIST, Vec::<Expr>::new())
}

/// `{{xs…}}`.
fn tuples(xs: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, [Expr::call(Symbol::LIST, xs)])
}

/// A system converted.
struct Problem<'a> {
    p: &'a Polynomials,
    /// Each unknown's variable, `None` for one that no equation holds.
    vars: &'a [Option<usize>],
    unknowns: &'a [Expr],
}

impl Problem<'_> {
    /// The rows of the numerators of `fractions`; `None` when one is not
    /// linear in the unknowns, or a denominator holds one.
    fn rows<'s, R: Numbers + 's>(
        &self,
        ring: &SparseRing<R>,
        fractions: impl Iterator<Item = (&'s Sparse<R::Elem>, &'s Sparse<R::Elem>)>,
    ) -> Option<Vec<Row<R::Elem>>> {
        let columns = self.vars.len();
        let mut rows = Vec::new();
        for (num, den) in fractions {
            if self
                .vars
                .iter()
                .flatten()
                .any(|&v| ring.degree(den, v).unwrap_or(0) > 0)
            {
                return None;
            }
            let mut parts: Vec<Vec<(Vec<u64>, R::Elem)>> = vec![Vec::new(); columns + 1];
            for (mut exps, c) in num.terms() {
                let mut at = columns;
                let mut degree = 0;
                for (j, var) in self.vars.iter().enumerate() {
                    if let &Some(v) = var
                        && exps[v] > 0
                    {
                        degree += exps[v];
                        exps[v] = 0;
                        at = j;
                    }
                }
                if degree > 1 {
                    return None;
                }
                parts[at].push((exps, c.clone()));
            }
            let mut row: Row<R::Elem> = parts.into_iter().map(|t| ring.from_terms(t)).collect();
            // Σ aⱼ xⱼ + c = 0: the right side is −c.
            row[columns] = ring.neg(&row[columns]);
            rows.push(row);
        }
        Some(rows)
    }

    /// The solutions of the rows of `fractions` over `ring`.
    fn solve<'s, R: Numbers + 's>(
        &self,
        ev: &mut Eval<'_>,
        e: &Expr,
        ring: &SparseRing<R>,
        fractions: impl Iterator<Item = (&'s Sparse<R::Elem>, &'s Sparse<R::Elem>)>,
    ) -> Result<Option<Expr>, EvalError> {
        let Some(mut rows) = self.rows(ring, fractions) else {
            return Ok(None);
        };
        let columns = self.vars.len();
        let cancel: &dyn Cancel = ev.context().interrupt_flag();
        let terms = ev.context().limits().terms;
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
                let large = (0..=columns).any(|k| {
                    pivot.term_count().saturating_mul(rows[i][k].term_count()) > terms
                        || c.term_count().saturating_mul(rows[rank][k].term_count()) > terms
                });
                if large {
                    ev.warn(Warning::Native {
                        message: format!("eliminating could pass {terms} terms").into(),
                        expr: e.clone(),
                    });
                    return Ok(None);
                }
                let row: Row<R::Elem> = (0..=columns)
                    .map(|k| ring.sub(&ring.mul(pivot, &rows[i][k]), &ring.mul(c, &rows[rank][k])))
                    .collect();
                rows[i] = primitive(ring, row, cancel)?;
            }
            pivots.push(j);
        }
        let rank = pivots.len();
        if rows[rank..].iter().any(|row| !row[columns].is_zero()) {
            return Ok(Some(none()));
        }
        let mut solution = Vec::with_capacity(columns);
        for (j, x) in self.unknowns.iter().enumerate() {
            let Some(i) = pivots.iter().position(|&c| c == j) else {
                solution.push(x.clone());
                continue;
            };
            let row = &rows[i];
            let pivot = &row[j];
            let mut parts = vec![self.ratio(ring, &row[columns], pivot, cancel)?];
            for (k, y) in self.unknowns.iter().enumerate() {
                if !pivots.contains(&k) && !row[k].is_zero() {
                    let c = self.ratio(ring, &ring.neg(&row[k]), pivot, cancel)?;
                    parts.push(Expr::call(Symbol::TIMES, [c, y.clone()]));
                }
            }
            solution.push(ev.eval(&Expr::call(Symbol::PLUS, parts))?);
        }
        Ok(Some(tuples(solution)))
    }

    /// `a / b` in lowest terms, built.
    fn ratio<R: Numbers>(
        &self,
        ring: &SparseRing<R>,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Expr, EvalError> {
        if a.is_zero() {
            return Ok(Expr::from(0));
        }
        let (n, d) = cancelled(ring, a, b, cancel)?;
        Ok(self.p.ratio_expr(&n, &d))
    }
}

/// `row` over the gcd of its entries.
fn primitive<R: Numbers>(
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
