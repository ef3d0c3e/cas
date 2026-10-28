//! Statistics (CAS-DOC-22): Xcas's descriptive statistics of lists.
//!
//! - The data (§2, CAS-294): a non-empty list of values, with a list of as
//!   many frequencies or not. A matrix, rows of one length, is data by
//!   columns: the answer is the list of each column's, and the frequencies
//!   are the rows'. Anything else leaves the call.
//! - The moments (§3, CAS-294), exact, with `x̄ = Σ f x / Σ f`:
//!   `Mean[l]` and `Mean[l, f]` (`mean`); `PopulationVariance`
//!   (`variance`), `Σ f (x − x̄)² / Σ f` expanded, a polynomial for symbols;
//!   `PopulationStandardDeviation` (`stddev`), its square root; and
//!   `StandardDeviation` (`stddevp`), the sample deviation
//!   `√(Σ f (x − x̄)² / (Σ f − 1))`. The frequencies may be any expressions,
//!   so that `mean([x, y], [p, 1 − p])` is an expected value. Frequencies
//!   whose sum is 0, or 1 for the sample deviation, leave the call, with a
//!   warning (giac's `undef`).

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// One column of data: its values, and their frequencies, or 1 each.
struct Sample {
    values: Vec<Expr>,
    frequencies: Option<Vec<Expr>>,
}

impl Sample {
    /// `Σ fᵢ g(xᵢ)`, unevaluated.
    fn sum(&self, g: impl Fn(&Expr) -> Expr) -> Expr {
        let terms: Vec<Expr> = match &self.frequencies {
            Some(f) => self.values.iter().zip(f).map(|(x, f)| f * g(x)).collect(),
            None => self.values.iter().map(g).collect(),
        };
        Expr::call(Symbol::PLUS, terms)
    }

    /// `Σ fᵢ`, unevaluated.
    fn total(&self) -> Expr {
        match &self.frequencies {
            Some(f) => Expr::call(Symbol::PLUS, f.clone()),
            None => Expr::from(i64::try_from(self.values.len()).unwrap_or(i64::MAX)),
        }
    }
}

/// The samples of the data `l` with the frequencies `f`: one of a list, one
/// per column of a matrix (the flag); `None` when they are neither.
fn samples(l: &Expr, f: Option<&Expr>) -> Option<(Vec<Sample>, bool)> {
    let rows = l.args();
    if !is_list(l) || rows.is_empty() {
        return None;
    }
    let frequencies = match f {
        Some(f) if is_list(f) && f.args().len() == rows.len() => Some(f.args().to_vec()),
        Some(_) => return None,
        None => None,
    };
    if !rows.iter().any(is_list) {
        let sample = Sample {
            values: rows.to_vec(),
            frequencies,
        };
        return Some((vec![sample], false));
    }
    let width = rows[0].args().len();
    if width == 0 || rows.iter().any(|r| !is_list(r) || r.args().len() != width) {
        return None;
    }
    let columns = (0..width)
        .map(|j| Sample {
            values: rows.iter().map(|r| r.args()[j].clone()).collect(),
            frequencies: frequencies.clone(),
        })
        .collect();
    Some((columns, true))
}

/// A statistic of one sample; `None` leaves the call.
type Statistic = fn(&mut Eval<'_>, &Expr, &Sample) -> Native;

/// `statistic` of the data of `e[l]` or `e[l, f]`: of a list, its value; of
/// a matrix, the list of its columns'.
fn of_data(ev: &mut Eval<'_>, e: &Expr, statistic: Statistic) -> Native {
    let (l, f) = match e.args() {
        [l] => (l, None),
        [l, f] => (l, Some(f)),
        _ => return Ok(None),
    };
    let Some((samples, columns)) = samples(l, f) else {
        return Ok(None);
    };
    let mut values = Vec::with_capacity(samples.len());
    for s in &samples {
        let Some(v) = statistic(ev, e, s)? else {
            return Ok(None);
        };
        values.push(v);
    }
    Ok(Some(if columns {
        Expr::call(Symbol::LIST, values)
    } else {
        values.swap_remove(0)
    }))
}

/// `Σ f`, evaluated, and `x̄`; `None` when `Σ f` is 0, with a warning.
fn total_and_mean(
    ev: &mut Eval<'_>,
    e: &Expr,
    s: &Sample,
) -> Result<Option<(Expr, Expr)>, EvalError> {
    let total = ev.eval(&s.total())?;
    if total == Expr::from(0) {
        refuse(ev, e, "the frequencies sum to 0")?;
        return Ok(None);
    }
    let mean = ev.eval(&(s.sum(Expr::clone) / &total))?;
    Ok(Some((total, mean)))
}

/// `Σ f (x − x̄)² / (Σ f − shift)`, expanded: the population variance for
/// a shift of 0, the square of the sample deviation for 1.
fn centered(ev: &mut Eval<'_>, e: &Expr, s: &Sample, shift: i64) -> Native {
    let Some((total, mean)) = total_and_mean(ev, e, s)? else {
        return Ok(None);
    };
    let denominator = ev.eval(&(total - shift))?;
    if denominator == Expr::from(0) {
        return refuse(ev, e, "the sample deviation of one value is undefined");
    }
    let squares = s.sum(|x| (x - &mean).pow(2));
    let variance = Expr::call(Symbol::EXPAND, [squares / denominator]);
    Ok(Some(ev.eval(&variance)?))
}

fn mean_of(ev: &mut Eval<'_>, e: &Expr, s: &Sample) -> Native {
    Ok(total_and_mean(ev, e, s)?.map(|(_, mean)| mean))
}

fn population_variance_of(ev: &mut Eval<'_>, e: &Expr, s: &Sample) -> Native {
    centered(ev, e, s, 0)
}

/// The square root of the centered sum over `Σ f − shift`.
fn deviation(ev: &mut Eval<'_>, e: &Expr, s: &Sample, shift: i64) -> Native {
    let Some(variance) = centered(ev, e, s, shift)? else {
        return Ok(None);
    };
    Ok(Some(ev.eval(&Expr::call(Symbol::SQRT, [variance]))?))
}

fn population_standard_deviation_of(ev: &mut Eval<'_>, e: &Expr, s: &Sample) -> Native {
    deviation(ev, e, s, 0)
}

fn standard_deviation_of(ev: &mut Eval<'_>, e: &Expr, s: &Sample) -> Native {
    deviation(ev, e, s, 1)
}

/// `Mean[l]`, `Mean[l, f]`.
pub(super) fn mean(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, mean_of)
}

/// `PopulationVariance[l]`, `PopulationVariance[l, f]`.
pub(super) fn population_variance(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, population_variance_of)
}

/// `PopulationStandardDeviation[l]`, `PopulationStandardDeviation[l, f]`.
pub(super) fn population_standard_deviation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, population_standard_deviation_of)
}

/// `StandardDeviation[l]`, `StandardDeviation[l, f]`: the sample deviation.
pub(super) fn standard_deviation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, standard_deviation_of)
}

#[cfg(test)]
mod tests;
