//! Matrices over `ℤ/nℤ` (CAS-DOC-44 §2): a matrix whose entries are
//! modular numbers of one modulus `n`, integers and rationals taken in, as
//! `[[1,2],[3,4]] % 7` is (`Modular` maps over a list).
//!
//! - `Det` and `CharacteristicPolynomial` (CAS-457), any `n`: the integer
//!   matrix of the symmetric residues has them over `ℤ`, exactly, by the
//!   rational paths; reduced modulo `n`, they are the matrix's, as no
//!   division was made. giac's `det([[1,2],[3,4]] % 6)` is `-2 % 6`.
//! - `Inverse`, any `n`: the inverse over `ℚ` of the integer matrix has
//!   denominators that divide its determinant `d`; when `d` is invertible
//!   modulo `n`, each entry reduced is `adj · d⁻¹`, the inverse modulo `n`.
//!   Otherwise the matrix stays, with a warning, as a singular rational one
//!   does (giac errs "Not invertible").
//! - `RowReduce`, `MatrixRank`, `NullSpace`, `ColumnSpace` and `LinSolve`
//!   (CAS-458), `n` a prime (or a probable one): Gauss–Jordan in `𝔽ₚ`,
//!   the entries of the reduced form modular values, so that the rational
//!   paths build kernels, images and solutions from it, giac's
//!   `ker([[1,2,3],[2,4,6]] % 5)` being `[[2 % 5,-1,0],[-2 % 5,0,-1]]`.
//!   A composite `n` leaves these with a warning: `ℤ/nℤ` has zero
//!   divisors, and giac's answers there are not reliable ("Reduction may
//!   fail").
//! - `Eigenvalues` (CAS-459), `n` a word prime: the roots in `𝔽ₚ` of the
//!   characteristic polynomial, each as often as its multiplicity,
//!   ascending, as giac's `eigenvals([[2,1],[1,2]] % 7)` is `1 % 7,3 % 7`;
//!   one that does not split over `𝔽ₚ` leaves the call with a warning,
//!   where giac builds a splitting field.
//! - `LinSolve` of equations with a modular value in them: each `l = r`
//!   as `Modular[l − r, p]`, linear in the unknowns, the augmented matrix
//!   reduced; a pivot unknown is its modular value in the free ones, a
//!   free one itself, as giac's `linsolve([x+y=1,2*x+2*y=2] % 5,[x,y])` is
//!   `[(-1 % 5)*y+1 % 5,y]`; no solution is `{}`.

use cas_num::{Integer, Primality, Rational};

use super::{modular, parts, residue};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

/// What a hook gives: `None` for a matrix that is not modular, which the
/// other paths take; `Some(None)` for a call that stays.
pub(in crate::native) type Hooked = Result<Option<Option<Expr>>, EvalError>;

/// A matrix over `ℤ/nℤ`: its entries the symmetric residues.
pub(in crate::native) struct ModularMatrix {
    pub(in crate::native) n: Integer,
    pub(in crate::native) rows: Vec<Vec<Integer>>,
}

impl ModularMatrix {
    /// The integer matrix of the residues.
    fn integers(&self) -> Expr {
        Expr::call(
            Symbol::LIST,
            self.rows
                .iter()
                .map(|r| Expr::call(Symbol::LIST, r.iter().cloned().map(Expr::from))),
        )
    }
}

/// `m`, a rectangular list of lists, as a matrix over `ℤ/nℤ`: one modular
/// entry at least, every modular entry an integer of the one modulus, the
/// others integers or rationals whose denominators are invertible.
pub(in crate::native) fn modular_matrix(m: &Expr) -> Option<ModularMatrix> {
    let rows = m.args();
    let width = rows.first()?.args().len();
    let mut n: Option<Integer> = None;
    for x in rows.iter().flat_map(Expr::args) {
        if let Some((_, k)) = parts(x) {
            match &n {
                Some(n) if *n != k => return None,
                _ => n = Some(k),
            }
        }
    }
    let n = n?;
    let rows: Vec<Vec<Integer>> = rows
        .iter()
        .map(|r| {
            if r.args().len() != width {
                return None;
            }
            r.args()
                .iter()
                .map(|x| {
                    let q: Rational = match parts(x) {
                        Some((a, _)) => Rational::from(a.to_integer()?),
                        None => x.to_rational()?,
                    };
                    residue(&q, &n)
                })
                .collect()
        })
        .collect::<Option<_>>()?;
    Some(ModularMatrix { n, rows })
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: String) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// The integer determinant of `m`'s residues, by `Det` over `ℤ`.
fn integer_det(ev: &mut Eval<'_>, m: &ModularMatrix) -> Result<Option<Integer>, EvalError> {
    let d = ev.eval(&Expr::call(Symbol::DET, [m.integers()]))?;
    Ok(d.to_integer())
}

/// `Det[m]` of a modular matrix.
pub(in crate::native) fn det(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Hooked {
    let Some(mm) = modular_matrix(m) else {
        return Ok(None);
    };
    let Some(d) = integer_det(ev, &mm)? else {
        return Ok(Some(None));
    };
    let result = modular(Expr::from(d), &mm.n);
    ev.step("det", e, &result);
    Ok(Some(Some(result)))
}

/// `CharacteristicPolynomial[m, x]`, or its coefficients without `x`, of a
/// modular matrix.
pub(in crate::native) fn characteristic_polynomial(
    ev: &mut Eval<'_>,
    e: &Expr,
    m: &Expr,
    x: Option<&Expr>,
) -> Hooked {
    let Some(mm) = modular_matrix(m) else {
        return Ok(None);
    };
    let mut args = vec![mm.integers()];
    args.extend(x.cloned());
    let p = ev.eval(&Expr::call(Symbol::CHARACTERISTIC_POLYNOMIAL, args))?;
    if p.head().and_then(Expr::as_symbol) == Some(Symbol::CHARACTERISTIC_POLYNOMIAL) {
        return Ok(Some(None));
    }
    let result = modular(p, &mm.n);
    ev.step("charpoly", e, &result);
    Ok(Some(Some(result)))
}

/// The inverse of a modular matrix, built; `Some(None)` after a warning
/// when its determinant is not invertible modulo `n`.
pub(in crate::native) fn inverse(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Hooked {
    let Some(mm) = modular_matrix(m) else {
        return Ok(None);
    };
    let Some(d) = integer_det(ev, &mm)? else {
        return Ok(Some(None));
    };
    if d.gcd(&mm.n) != Integer::ONE {
        warn(
            ev,
            e,
            format!("the matrix is not invertible modulo {}", mm.n),
        );
        return Ok(Some(None));
    }
    let inverse = ev.eval(&Expr::call(Symbol::INVERSE, [mm.integers()]))?;
    if inverse.head().and_then(Expr::as_symbol) == Some(Symbol::INVERSE) {
        return Ok(Some(None));
    }
    // A list, entry by entry: products and powers take the shape.
    let n = &mm.n;
    let entries =
        |r: &Expr| Expr::call(Symbol::LIST, r.args().iter().map(|x| modular(x.clone(), n)));
    Ok(Some(Some(Expr::call(
        Symbol::LIST,
        inverse.args().iter().map(entries),
    ))))
}

/// Whether `n` is a prime, proven or by Baillie–PSW.
fn prime(n: &Integer) -> bool {
    n.unsigned_abs().primality() != Primality::NotPrime
}

/// The symmetric residue of `x` in `0..p`.
fn symmetric(x: &Integer, p: &Integer) -> Integer {
    if &(x + x) > p { x - p } else { x.clone() }
}

/// Gauss–Jordan modulo the prime `p` on `rows`, entries in `0..p`, the
/// pivots sought among the first `columns` columns: the pivot columns, the
/// `k`-th in row `k`, each pivot 1, the rows past the rank 0.
fn gauss_jordan(rows: &mut [Vec<Integer>], columns: usize, p: &Integer) -> Vec<usize> {
    let mut pivots = Vec::new();
    for j in 0..columns {
        let r = pivots.len();
        if r == rows.len() {
            break;
        }
        let Some(k) = (r..rows.len()).find(|&k| !rows[k][j].is_zero()) else {
            continue;
        };
        rows.swap(r, k);
        let inverse = rows[r][j].mod_inverse(p).expect("a unit of a field");
        for x in &mut rows[r] {
            *x = (&*x * &inverse).div_rem_euclid(p).1;
        }
        let pivot = rows[r].clone();
        for (k, row) in rows.iter_mut().enumerate() {
            if k == r || row[j].is_zero() {
                continue;
            }
            let c = row[j].clone();
            for (x, y) in row.iter_mut().zip(&pivot) {
                *x = (&*x - &(&c * y)).div_rem_euclid(p).1;
            }
        }
        pivots.push(j);
    }
    pivots
}

/// A reduced row echelon form: its rows, built, and its pivot columns.
pub(in crate::native) type Echelon = (Vec<Vec<Expr>>, Vec<usize>);

/// What the echelon hook gives, as [`Hooked`].
pub(in crate::native) type HookedEchelon = Result<Option<Option<Echelon>>, EvalError>;

/// The reduced row echelon form of a modular matrix over `𝔽ₚ`, its entries
/// modular values, and its pivot columns; `Some(None)` after a warning for
/// a composite modulus.
pub(in crate::native) fn reduced(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> HookedEchelon {
    let Some(mm) = modular_matrix(m) else {
        return Ok(None);
    };
    let p = &mm.n;
    if !prime(p) {
        warn(ev, e, format!("row reduction modulo the composite {p}"));
        return Ok(Some(None));
    }
    let mut rows: Vec<Vec<Integer>> = mm
        .rows
        .iter()
        .map(|r| r.iter().map(|x| x.div_rem_euclid(p).1).collect())
        .collect();
    let columns = rows[0].len();
    let pivots = gauss_jordan(&mut rows, columns, p);
    let rows = rows
        .iter()
        .map(|r| {
            r.iter()
                .map(|x| modular(Expr::from(symmetric(x, p)), p))
                .collect()
        })
        .collect();
    Ok(Some(Some((rows, pivots))))
}

/// The modulus of the modular values in `e`; `None` for none, or for two.
fn modulus(e: &Expr) -> Option<Integer> {
    let mut found: Option<Integer> = None;
    let mut two = false;
    let mut stack = vec![e];
    while let Some(x) = stack.pop() {
        if let Some((_, n)) = parts(x) {
            match &found {
                Some(m) if *m != n => two = true,
                _ => found = Some(n),
            }
            continue;
        }
        stack.extend(x.args());
    }
    found.filter(|_| !two)
}

/// The integer coefficients of the unknowns in `d`, a sum of terms `c xᵢ`
/// and a constant, and the constant; `None` for another form.
fn linear(d: &Expr, unknowns: &[Expr]) -> Option<(Vec<Integer>, Integer)> {
    let terms: &[Expr] = match d.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => d.args(),
        _ => std::slice::from_ref(d),
    };
    let mut coefficients = vec![Integer::ZERO; unknowns.len()];
    let mut constant = Integer::ZERO;
    for t in terms {
        if let Some(c) = t.to_integer() {
            constant = &constant + &c;
            continue;
        }
        let (c, x) = match (t.head().and_then(Expr::as_symbol), t.args()) {
            (Some(Symbol::TIMES), [c, x]) => (c.to_integer()?, x),
            _ => (Integer::ONE, t),
        };
        let i = unknowns.iter().position(|u| u == x)?;
        coefficients[i] = &coefficients[i] + &c;
    }
    Some((coefficients, constant))
}

/// `LinSolve[eqns, vars]` with modular values in the equations.
pub(in crate::native) fn lin_solve(
    ev: &mut Eval<'_>,
    e: &Expr,
    equations: &[Expr],
    unknowns: &[Expr],
) -> Hooked {
    let Some(p) = equations.iter().find_map(modulus) else {
        return Ok(None);
    };
    if equations.iter().any(|q| modulus(q).is_some_and(|n| n != p)) {
        return Ok(None);
    }
    if !prime(&p) {
        warn(ev, e, format!("a linear system modulo the composite {p}"));
        return Ok(Some(None));
    }
    let columns = unknowns.len();
    let mut rows = Vec::with_capacity(equations.len());
    for q in equations {
        let (l, r) = match (q.head().and_then(Expr::as_symbol), q.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => (l.clone(), r.clone()),
            _ => (q.clone(), Expr::from(0)),
        };
        let difference = Expr::call(Symbol::PLUS, [l, Expr::from(-1) * r]);
        let d = ev.eval(&modular(difference, &p))?;
        let linear = match parts(&d) {
            Some((d, _)) => linear(d, unknowns),
            None => None,
        };
        let Some((coefficients, constant)) = linear else {
            warn(ev, e, format!("not a linear system modulo {p}"));
            return Ok(Some(None));
        };
        let mut row: Vec<Integer> = coefficients
            .iter()
            .map(|c| c.div_rem_euclid(&p).1)
            .collect();
        row.push((-constant).div_rem_euclid(&p).1);
        rows.push(row);
    }
    let pivots = gauss_jordan(&mut rows, columns, &p);
    if rows[pivots.len()..].iter().any(|r| !r[columns].is_zero()) {
        ev.step("linsolve", e, &Expr::call(Symbol::LIST, Vec::<Expr>::new()));
        return Ok(Some(Some(Expr::call(Symbol::LIST, Vec::<Expr>::new()))));
    }
    let free: Vec<usize> = (0..columns).filter(|j| !pivots.contains(j)).collect();
    let solution = (0..columns).map(|j| match pivots.iter().position(|&c| c == j) {
        Some(i) => {
            let mut terms = vec![Expr::from(symmetric(&rows[i][columns], &p))];
            for &f in &free {
                let c = symmetric(&rows[i][f], &p);
                terms.push(Expr::from(-c) * unknowns[f].clone());
            }
            modular(Expr::call(Symbol::PLUS, terms), &p)
        }
        None => unknowns[j].clone(),
    });
    let result = Expr::call(Symbol::LIST, solution.collect::<Vec<_>>());
    ev.step("linsolve", e, &result);
    Ok(Some(Some(result)))
}

/// `Eigenvalues[m]` of a modular matrix.
pub(in crate::native) fn eigenvalues(ev: &mut Eval<'_>, e: &Expr, m: &Expr) -> Hooked {
    let Some(mm) = modular_matrix(m) else {
        return Ok(None);
    };
    let x = Expr::symbol("Eigenvalues$x");
    let p = ev.eval(&Expr::call(
        Symbol::CHARACTERISTIC_POLYNOMIAL,
        [mm.integers(), x.clone()],
    ))?;
    let reduced = ev.eval(&modular(p, &mm.n))?;
    let Some(u) = super::polynomials::univariate(&reduced) else {
        warn(
            ev,
            e,
            format!("eigenvalues modulo {}, not a word prime", mm.n),
        );
        return Ok(Some(None));
    };
    let Some(roots) = super::polynomials::split(&u) else {
        warn(
            ev,
            e,
            format!(
                "the characteristic polynomial does not split modulo {}",
                mm.n
            ),
        );
        return Ok(Some(None));
    };
    let values = roots
        .into_iter()
        .flat_map(|(r, k)| std::iter::repeat_n(modular(Expr::from(r), &mm.n), k));
    let result = Expr::call(Symbol::LIST, values.collect::<Vec<_>>());
    ev.step("eigenvals", e, &result);
    Ok(Some(Some(result)))
}
