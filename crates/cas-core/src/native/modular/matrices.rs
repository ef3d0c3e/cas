//! Matrices over `ℤ/nℤ` (CAS-DOC-44 §2): a matrix whose entries are
//! modular numbers of one modulus `n`, integers and rationals taken in, as
//! `[[1,2],[3,4]] % 7` is (`Modular` is `Listable`).
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

use cas_num::{Integer, Rational};

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
