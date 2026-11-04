//! Frequencies, classes and other descriptive statistics (CAS-DOC-50):
//! Xcas's `frequencies`, `cumulated_frequencies`, `classes`,
//! `covariance_correlation`, `rms` and `moving_average`.
//!
//! - **`Frequencies[l]`** (§2, CAS-491): the distinct values of a list of
//!   real numbers, ascending (`Sort`), each with its count over `n`, exact
//!   where giac writes decimals: `[[1,1/4],[2,1/2],[3,1/4]]` for `[3,1,2,2]`.
//!   **`CumulatedFrequencies[l]`**: the running sums of those.
//! - **`Classes[l, a, w]`**: the classes `[a + k w, a + (k + 1) w)` from the
//!   one of the least value to the one of the largest, empty ones kept, each
//!   `{Range[lo, hi], count}`, giac's `[[0.0 .. 3.0,2],[3.0 .. 6.0,3],[6.0 ..
//!   9.0,2]]` for `1, …, 7` from 0 by 3; `Classes[l]` from the least value
//!   by 1. A value below `a`, or values that are not real numbers, leave the
//!   call.
//! - **`CovarianceCorrelation[…]`** (§3, CAS-492): `{PopulationCovariance,
//!   Correlation}` of the same data, giac's `[5/3,15/(2*sqrt(57))]`.
//! - **`RootMeanSquare[l]`**: `√(Σ xᵢ²/n)`. **`MovingAverage[l, k]`**: the
//!   means of the `n − k + 1` windows of `k` consecutive values.

use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most classes of `Classes`.
const MAX_CLASSES: i64 = 1 << 16;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The values of a non-empty list of real numbers, sorted by value.
fn sorted(ev: &mut Eval<'_>, l: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    if !is_list(l) || l.args().is_empty() || !l.args().iter().all(Expr::is_numeric) {
        return Ok(None);
    }
    let s = ev.eval(&Expr::call(Symbol::SORT, [l.clone()]))?;
    Ok(is_list(&s).then(|| s.args().to_vec()))
}

/// The distinct values with their counts, ascending.
fn counts(values: &[Expr]) -> Vec<(Expr, i64)> {
    let mut out: Vec<(Expr, i64)> = Vec::new();
    for v in values {
        match out.last_mut() {
            Some((last, k)) if last == v => *k += 1,
            _ => out.push((v.clone(), 1)),
        }
    }
    out
}

fn frequencies_of(ev: &mut Eval<'_>, e: &Expr, cumulated: bool) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    let Some(values) = sorted(ev, l)? else {
        return Ok(None);
    };
    let n = Integer::from(i64::try_from(values.len()).expect("a length"));
    let mut running = 0i64;
    let pairs = counts(&values)
        .into_iter()
        .map(|(v, k)| {
            running += k;
            let count = if cumulated { running } else { k };
            list(vec![
                v,
                Expr::from(Rational::new(Integer::from(count), n.clone())),
            ])
        })
        .collect();
    let result = list(pairs);
    ev.step(
        if cumulated {
            "cumulated_frequencies"
        } else {
            "frequencies"
        },
        e,
        &result,
    );
    Ok(Some(result))
}

/// `Frequencies[l]`.
pub(in crate::native) fn frequencies(ev: &mut Eval<'_>, e: &Expr) -> Native {
    frequencies_of(ev, e, false)
}

/// `CumulatedFrequencies[l]`.
pub(in crate::native) fn cumulated_frequencies(ev: &mut Eval<'_>, e: &Expr) -> Native {
    frequencies_of(ev, e, true)
}

/// `Classes[l, a, w]`, `Classes[l]`.
pub(in crate::native) fn classes(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, start, width) = match e.args() {
        [l] => (l, None, Expr::from(1)),
        [l, a, w] => (l, Some(a.clone()), w.clone()),
        _ => return Ok(None),
    };
    let Some(values) = sorted(ev, l)? else {
        return Ok(None);
    };
    let start = start.unwrap_or_else(|| values[0].clone());
    let positive = ev.eval(&Expr::call(Symbol::GREATER, [width.clone(), Expr::from(0)]))?;
    if positive.as_symbol() != Some(Symbol::TRUE) {
        return Ok(None);
    }
    // The class of each value, k = ⌊(x − a)/w⌋.
    let mut indices = Vec::with_capacity(values.len());
    for x in &values {
        let k = ev.eval(&Expr::call(
            Symbol::FLOOR,
            [(x.clone() - start.clone()) / width.clone()],
        ))?;
        let Some(k) = k.as_i64().filter(|k| (0..MAX_CLASSES).contains(k)) else {
            return Ok(None);
        };
        indices.push(k);
    }
    let (first, last) = (indices[0], *indices.last().expect("not empty"));
    let mut rows = Vec::with_capacity(usize::try_from(last - first + 1).expect("small"));
    for k in first..=last {
        let count = indices.iter().filter(|&&j| j == k).count();
        let lo = ev.eval(&(start.clone() + Expr::from(k) * width.clone()))?;
        let hi = ev.eval(&(start.clone() + Expr::from(k + 1) * width.clone()))?;
        rows.push(list(vec![
            Expr::call(Symbol::RANGE, [lo, hi]),
            Expr::from(i64::try_from(count).expect("a count")),
        ]));
    }
    let result = list(rows);
    ev.step("classes", e, &result);
    Ok(Some(result))
}

/// `CovarianceCorrelation[…]`.
pub(in crate::native) fn covariance_correlation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let covariance = ev.eval(&Expr::call(
        Symbol::POPULATION_COVARIANCE,
        e.args().to_vec(),
    ))?;
    let correlation = ev.eval(&Expr::call(Symbol::CORRELATION, e.args().to_vec()))?;
    if covariance.head() == Some(&Expr::from(Symbol::POPULATION_COVARIANCE))
        || correlation.head() == Some(&Expr::from(Symbol::CORRELATION))
    {
        return Ok(None);
    }
    let result = list(vec![covariance, correlation]);
    ev.step("covariance_correlation", e, &result);
    Ok(Some(result))
}

/// `RootMeanSquare[l]`.
pub(in crate::native) fn root_mean_square(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) || l.args().is_empty() || l.args().iter().any(is_list) {
        return Ok(None);
    }
    let n = Expr::from(i64::try_from(l.args().len()).expect("a length"));
    let squares = Expr::call(
        Symbol::PLUS,
        l.args().iter().map(|x| x.pow(2)).collect::<Vec<_>>(),
    );
    let result = ev.eval(&Expr::call(Symbol::SQRT, [squares / n]))?;
    ev.step("rms", e, &result);
    Ok(Some(result))
}

/// `MovingAverage[l, k]`.
pub(in crate::native) fn moving_average(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, k] = e.args() else {
        return Ok(None);
    };
    let n = l.args().len();
    let Some(k) = k
        .as_i64()
        .and_then(|k| usize::try_from(k).ok())
        .filter(|k| (1..=n).contains(k))
    else {
        return Ok(None);
    };
    if !is_list(l) || l.args().iter().any(is_list) {
        return Ok(None);
    }
    let size = Expr::from(i64::try_from(k).expect("small"));
    let means = l
        .args()
        .windows(k)
        .map(|w| ev.eval(&(Expr::call(Symbol::PLUS, w.to_vec()) / size.clone())))
        .collect::<Result<Vec<_>, _>>()?;
    let result = list(means);
    ev.step("moving_average", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
