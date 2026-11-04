//! Matrix norms, condition numbers and Gram–Schmidt (CAS-DOC-47 §4,
//! CAS-478), and singular values and conjugate gradients (§5, CAS-479).
//!
//! - **`MatrixNorm[A, k]`**, Xcas's `matrix_norm`: `k = 1`, the largest
//!   column sum of `|a_ij|`; `Infinity`, the largest row sum; `frobenius`,
//!   `√Σ|a_ij|²`; `2`, the largest singular value. `L1Norm` and
//!   `LInfinityNorm`, `l1norm` and `linfnorm`: those of a matrix, and of a
//!   vector the sum and the largest of `|vᵢ|`.
//! - **`ConditionNumber[A, k]`**, `cond`: `‖A‖ ‖A⁻¹‖`, the 1-norm by
//!   default, giac's `cond([[1,2],[3,4]])` being 21; `Infinity` for a
//!   singular `A`.
//! - **`GramSchmidt[B]`**, `gramschmidt`: the rows orthogonalized in order,
//!   `vₖ = bₖ − Σ ⟨uⱼ, bₖ⟩ uⱼ`, and normalized, `uₖ = vₖ/√⟨vₖ, vₖ⟩`, the
//!   inner product conjugating its first vector; a row that depends on the
//!   ones before is dropped.

use super::{inner, is_zero, list, matrix, rows, simplified};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn abs(x: &Expr) -> Expr {
    Expr::call(Symbol::ABS, [x.clone()])
}

fn max(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::MAX, items)
}

fn sum(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::PLUS, items)
}

/// Which norm.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    One,
    Two,
    Infinity,
    Frobenius,
}

fn kind(k: &Expr) -> Option<Kind> {
    if *k == Expr::from(1) {
        return Some(Kind::One);
    }
    if *k == Expr::from(2) {
        return Some(Kind::Two);
    }
    match k.as_symbol() {
        Some(Symbol::INFINITY) => Some(Kind::Infinity),
        Some(s) if s.name() == "frobenius" => Some(Kind::Frobenius),
        _ => None,
    }
}

/// The singular values of the matrix `a`, from the eigenvalues of `Aᴴ A`,
/// ascending; `None` when they are not found.
pub(super) fn singular_values(ev: &mut Eval<'_>, a: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    let adjoint = Expr::call(
        Symbol::CONJUGATE,
        [Expr::call(Symbol::TRANSPOSE, [a.clone()])],
    );
    let gram = ev.eval(&Expr::call(Symbol::TIMES, [adjoint, a.clone()]))?;
    let values = ev.eval(&Expr::call(Symbol::EIGENVALUES, [gram]))?;
    if values.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let roots: Vec<Expr> = values
        .args()
        .iter()
        .map(|v| {
            ev.eval(&Expr::call(
                Symbol::SQRT,
                [Expr::call(Symbol::RE, [v.clone()])],
            ))
        })
        .collect::<Result<_, _>>()?;
    let sorted = ev.eval(&Expr::call(Symbol::SORT, [Expr::call(Symbol::LIST, roots)]))?;
    Ok((sorted.head() == Some(&Expr::from(Symbol::LIST))).then(|| sorted.args().to_vec()))
}

/// `‖A‖_k`, unevaluated but for the 2-norm; `None` when not found.
fn norm_of(
    ev: &mut Eval<'_>,
    a: &Expr,
    m: &[Vec<Expr>],
    k: Kind,
) -> Result<Option<Expr>, EvalError> {
    let width = m[0].len();
    Ok(Some(match k {
        Kind::One => max((0..width)
            .map(|j| sum(m.iter().map(|r| abs(&r[j])).collect()))
            .collect()),
        Kind::Infinity => max(m.iter().map(|r| sum(r.iter().map(abs).collect())).collect()),
        Kind::Frobenius => Expr::call(
            Symbol::SQRT,
            [sum(m.iter().flatten().map(|x| abs(x).pow(2)).collect())],
        ),
        Kind::Two => match singular_values(ev, a)? {
            Some(values) if !values.is_empty() => values.last().expect("not empty").clone(),
            _ => return Ok(None),
        },
    }))
}

/// `MatrixNorm[A, k]`.
pub(in crate::native) fn matrix_norm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, k] = e.args() else {
        return Ok(None);
    };
    let (Some(m), Some(k)) = (rows(a), kind(k)) else {
        return Ok(None);
    };
    let Some(n) = norm_of(ev, a, &m, k)? else {
        return Ok(None);
    };
    let result = ev.eval(&n)?;
    ev.step("matrix_norm", e, &result);
    Ok(Some(result))
}

/// A vector's or a matrix's norm of `kind`, one or ∞.
fn vector_or_matrix(ev: &mut Eval<'_>, e: &Expr, k: Kind, name: &'static str) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    if a.head() != Some(&Expr::from(Symbol::LIST)) || a.args().is_empty() {
        return Ok(None);
    }
    let vector = !a
        .args()
        .iter()
        .any(|x| x.head() == Some(&Expr::from(Symbol::LIST)));
    let n = if vector {
        let parts: Vec<Expr> = a.args().iter().map(abs).collect();
        if k == Kind::One {
            sum(parts)
        } else {
            max(parts)
        }
    } else {
        let Some(m) = rows(a) else {
            return Ok(None);
        };
        let Some(n) = norm_of(ev, a, &m, k)? else {
            return Ok(None);
        };
        n
    };
    let result = ev.eval(&n)?;
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// `L1Norm[A]`.
pub(in crate::native) fn l1_norm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vector_or_matrix(ev, e, Kind::One, "l1norm")
}

/// `LInfinityNorm[A]`.
pub(in crate::native) fn l_infinity_norm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vector_or_matrix(ev, e, Kind::Infinity, "linfnorm")
}

/// `ConditionNumber[A]`, `ConditionNumber[A, k]`.
pub(in crate::native) fn condition_number(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (a, k) = match e.args() {
        [a] => (a, Kind::One),
        [a, k] => match kind(k) {
            Some(k) => (a, k),
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    let Some(m) = rows(a) else {
        return Ok(None);
    };
    if m[0].len() != m.len() {
        return Ok(None);
    }
    let result = if k == Kind::Two {
        let Some(values) = singular_values(ev, a)? else {
            return Ok(None);
        };
        let (Some(small), Some(large)) = (values.first(), values.last()) else {
            return Ok(None);
        };
        if is_zero(ev, small)? == Some(true) {
            Expr::from(Symbol::INFINITY)
        } else {
            ev.eval(&(large.clone() / small.clone()))?
        }
    } else {
        let inverse = ev.eval(&Expr::call(Symbol::INVERSE, [a.clone()]))?;
        let Some(b) = rows(&inverse) else {
            // A singular matrix: its inverse stays.
            if ev.eval(&Expr::call(Symbol::DET, [a.clone()]))? == Expr::from(0) {
                return Ok(Some(Expr::from(Symbol::INFINITY)));
            }
            return Ok(None);
        };
        let (Some(x), Some(y)) = (norm_of(ev, a, &m, k)?, norm_of(ev, &inverse, &b, k)?) else {
            return Ok(None);
        };
        ev.eval(&(x * y))?
    };
    ev.step("cond", e, &result);
    Ok(Some(result))
}

/// `GramSchmidt[B]`.
pub(in crate::native) fn gram_schmidt(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [b] = e.args() else {
        return Ok(None);
    };
    let Some(m) = rows(b) else {
        return Ok(None);
    };
    let mut basis: Vec<Vec<Expr>> = Vec::with_capacity(m.len());
    for row in &m {
        let mut v = row.clone();
        for u in &basis {
            let Some(c) = inner(ev, u, row)? else {
                return Ok(None);
            };
            for (x, y) in v.iter_mut().zip(u) {
                *x = simplified(ev, x.clone() - c.clone() * y.clone())?;
            }
        }
        let mut zero = true;
        for x in &v {
            match is_zero(ev, x)? {
                Some(true) => {}
                Some(false) => zero = false,
                None => return Ok(None),
            }
        }
        if zero {
            continue;
        }
        let Some(square) = inner(ev, &v, &v)? else {
            return Ok(None);
        };
        let length = ev.eval(&Expr::call(Symbol::SQRT, [square]))?;
        let unit: Vec<Expr> = v
            .into_iter()
            .map(|x| simplified(ev, x / length.clone()))
            .collect::<Result<_, _>>()?;
        basis.push(unit);
    }
    let result = if basis.is_empty() {
        list(Vec::new())
    } else {
        matrix(basis)
    };
    ev.step("gramschmidt", e, &result);
    Ok(Some(result))
}
