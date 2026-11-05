//! Pearson's χ² (CAS-DOC-52 §3, CAS-504): `ChiSquareTest` (`chisquaret`).
//!
//! - `ChiSquareTest[o]`: the counts against the uniform law, `eᵢ = N/k`.
//! - `ChiSquareTest[o, p]`, `p` summing to 1: against the probabilities,
//!   `eᵢ = N pᵢ`; `o` of another length than `p` is raw data, counted by
//!   its distinct values ascending (`Frequencies`).
//! - `ChiSquareTest[o₁, o₂]`, `o₂` not summing to 1, and
//!   `ChiSquareTest[m]`, a matrix of two rows or more: homogeneity, `Eᵢⱼ =
//!   Rᵢ Cⱼ / N`.
//!
//! The answer is `χ² = Σ (o − e)²/e`, exact for exact counts, as giac's;
//! the step notes the degrees of freedom, `k − 1` or `(r − 1)(k − 1)`, and
//! `1 − ChiSquareCDF(df, χ²)`. A level may follow, as giac takes one; the
//! answer does not depend on it. Counts are real numbers not below 0 and the expected ones
//! above 0; otherwise the call stays with a warning.

use std::cmp::Ordering;

use super::{Native, positive, real};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::rounding::real_sign;
use crate::native::statistics::{is_list, refuse};
use crate::symbol::Symbol;

/// The rows of `m`, a matrix of two rows or more; `None` otherwise.
fn table(m: &Expr) -> Option<Vec<Vec<Expr>>> {
    let rows = m.args();
    if !is_list(m) || rows.len() < 2 || !rows.iter().all(is_list) {
        return None;
    }
    let width = rows[0].args().len();
    (width >= 1 && rows.iter().all(|r| r.args().len() == width))
        .then(|| rows.iter().map(|r| r.args().to_vec()).collect())
}

/// Whether `p`'s entries sum to 1: exactly, or within the balls.
fn sums_to_one(ev: &mut Eval<'_>, p: &Expr) -> Result<bool, crate::context::EvalError> {
    let gap = ev.eval(&(Expr::call(Symbol::PLUS, p.args().to_vec()) - Expr::from(1)))?;
    Ok(match gap.view() {
        View::Ball(b) => b.contains_rational(&cas_num::Rational::ZERO),
        _ => gap == Expr::from(0),
    })
}

/// `Σ`, unevaluated.
fn sum(terms: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::PLUS, terms.into_iter().collect::<Vec<_>>())
}

/// `Σ (o − e)²/e` over the pairs, the expected values shown positive and
/// the counts not negative; `None` after a warning.
fn statistic(ev: &mut Eval<'_>, e: &Expr, pairs: &[(Expr, Expr)]) -> Native {
    let mut terms = Vec::with_capacity(pairs.len());
    for (o, x) in pairs {
        let x = ev.eval(x)?;
        if real_sign(ev.context(), o).is_none_or(|s| s == Ordering::Less) {
            return refuse(ev, e, "the counts are real numbers not below 0");
        }
        if !positive(ev, &x) {
            return refuse(ev, e, "an expected count is not shown positive");
        }
        terms.push((o.clone() - x.clone()).pow(2) / x);
    }
    Ok(Some(ev.eval(&sum(terms))?))
}

/// The counts against the probabilities `p`, of one length.
fn adequation(ev: &mut Eval<'_>, e: &Expr, o: &[Expr], p: &[Expr]) -> Native {
    let total = ev.eval(&sum(o.to_vec()))?;
    let pairs: Vec<(Expr, Expr)> = o
        .iter()
        .zip(p)
        .map(|(o, p)| (o.clone(), total.clone() * p.clone()))
        .collect();
    statistic(ev, e, &pairs)
}

/// The homogeneity of the rows.
fn homogeneity(ev: &mut Eval<'_>, e: &Expr, rows: &[Vec<Expr>]) -> Native {
    let width = rows[0].len();
    let row_totals: Vec<Expr> = rows
        .iter()
        .map(|r| ev.eval(&sum(r.clone())))
        .collect::<Result<_, _>>()?;
    let column_totals: Vec<Expr> = (0..width)
        .map(|j| ev.eval(&sum(rows.iter().map(|r| r[j].clone()))))
        .collect::<Result<_, _>>()?;
    let total = ev.eval(&sum(row_totals.clone()))?;
    let mut pairs = Vec::with_capacity(rows.len() * width);
    for (r, rt) in rows.iter().zip(&row_totals) {
        for (o, ct) in r.iter().zip(&column_totals) {
            pairs.push((o.clone(), rt.clone() * ct.clone() / total.clone()));
        }
    }
    statistic(ev, e, &pairs)
}

/// The counts of the distinct values of `data`, ascending; `None` when
/// `Frequencies` does not give them.
fn counts(ev: &mut Eval<'_>, data: &Expr) -> Result<Option<Vec<Expr>>, crate::context::EvalError> {
    let f = ev.eval(&Expr::call(Symbol::FREQUENCIES, [data.clone()]))?;
    if !is_list(&f)
        || !f
            .args()
            .iter()
            .all(|pair| is_list(pair) && pair.args().len() == 2)
    {
        return Ok(None);
    }
    let n = Expr::from(i64::try_from(data.args().len()).unwrap_or(i64::MAX));
    f.args()
        .iter()
        .map(|pair| ev.eval(&(pair.args()[1].clone() * n.clone())))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

/// `ChiSquareTest[…]`, giac's `chisquaret`.
pub(in crate::native) fn chi_square_test(ev: &mut Eval<'_>, e: &Expr) -> Native {
    // giac's trailing level, which the answer does not take.
    let args: Vec<Expr> = match e.args() {
        [a, b, alpha] if !is_list(alpha) && positive(ev, alpha) => vec![a.clone(), b.clone()],
        args => args.to_vec(),
    };
    let (chi, df) = match args.as_slice() {
        [m] if table(m).is_some() => {
            let rows = table(m).expect("a table");
            let df = (rows.len() - 1) * (rows[0].len() - 1);
            (homogeneity(ev, e, &rows)?, df)
        }
        [o] if is_list(o) && !o.args().is_empty() => {
            let k = o.args().len();
            let p = Expr::from(cas_num::Rational::new(
                1.into(),
                i64::try_from(k).unwrap_or(1).into(),
            ));
            (adequation(ev, e, o.args(), &vec![p; k])?, k - 1)
        }
        [o, p] if is_list(o) && is_list(p) && !p.args().is_empty() && !o.args().is_empty() => {
            if !p.args().iter().all(|x| real(ev, x)) {
                return Ok(None);
            }
            if sums_to_one(ev, p)? {
                let k = p.args().len();
                let observed = if o.args().len() == k {
                    o.args().to_vec()
                } else {
                    match counts(ev, o)? {
                        Some(found) if found.len() == k => found,
                        _ => {
                            return refuse(
                                ev,
                                e,
                                "the data has not as many distinct values as probabilities",
                            );
                        }
                    }
                };
                (adequation(ev, e, &observed, p.args())?, k - 1)
            } else if o.args().len() == p.args().len() {
                let rows = vec![o.args().to_vec(), p.args().to_vec()];
                (homogeneity(ev, e, &rows)?, o.args().len() - 1)
            } else {
                return refuse(ev, e, "two samples of different lengths");
            }
        }
        _ => return Ok(None),
    };
    let Some(chi) = chi else {
        return Ok(None);
    };
    if ev.tracing() {
        let df = Expr::from(i64::try_from(df).unwrap_or(i64::MAX));
        let p = Expr::from(1) - Expr::call(Symbol::CHI_SQUARE_CDF, [df.clone(), chi.clone()]);
        let p = ev.eval(&Expr::call(Symbol::N, [p, Expr::from(10)]))?;
        ev.step(format!("χ² test: df = {df}, p = {p}"), e, &chi);
    }
    Ok(Some(chi))
}
