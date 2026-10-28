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
//! - The order statistics (§4, CAS-295), of real values that comparisons
//!   tell apart ([`order`]), with non-negative real frequencies:
//!   `Median[l]`, `Median[l, f]` (`median`), the value at which the
//!   cumulative frequency passes half the total, or the half-sum of the
//!   two values around it when it reaches half exactly (giac answers the
//!   lower one); `Quantile[l, p]`, `Quantile[l, f, p]` (`quantile`), for
//!   `0 ≤ p ≤ 1` the smallest value whose cumulative frequency is at least
//!   `p Σ f`, `x₍⌈np⌉₎` without frequencies and the minimum for `p = 0`
//!   (giac: `undef` at 0 and 1); `Quartiles` (`quartiles`), `{min, Q(1/4),
//!   median, Q(3/4), max}`, and `FirstQuartile` (`quartile1`) and
//!   `ThirdQuartile` (`quartile3`).
//! - Two variables (§3, CAS-296), as `[x, y]`, `[x, y, f]`, or a matrix of
//!   two columns `[m]` or `[m, f]`, with `S_xy = Σ f (x − x̄)(y − ȳ)`:
//!   `PopulationCovariance` (`covariance`), `S_xy / Σ f` expanded, and
//!   `Correlation` (`correlation`), `S_xy / √(S_xx S_yy)`, which is
//!   `Indeterminate` for constant data (giac's `undef`).
//! - The regressions (§5, CAS-297), least squares: [`regression`].
//! - Probability distributions (CAS-DOC-23): [`distributions`].

pub(super) mod distributions;
mod order;
pub(super) mod regression;

use cas_num::Rational;

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

/// `statistic` of the data `l` with the frequencies `f`: of a list, its
/// value; of a matrix, the list of its columns'. `None` leaves the call.
fn of_samples(
    ev: &mut Eval<'_>,
    l: &Expr,
    f: Option<&Expr>,
    mut statistic: impl FnMut(&mut Eval<'_>, &Sample) -> Native,
) -> Native {
    let Some((samples, columns)) = samples(l, f) else {
        return Ok(None);
    };
    let mut values = Vec::with_capacity(samples.len());
    for s in &samples {
        let Some(v) = statistic(ev, s)? else {
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

/// A statistic of one sample, given the call for its warnings; `None`
/// leaves the call.
type Statistic = fn(&mut Eval<'_>, &Expr, &Sample) -> Native;

/// `statistic` of the data of `e[l]` or `e[l, f]`.
fn of_data(ev: &mut Eval<'_>, e: &Expr, statistic: Statistic) -> Native {
    let (l, f) = match e.args() {
        [l] => (l, None),
        [l, f] => (l, Some(f)),
        _ => return Ok(None),
    };
    of_samples(ev, l, f, |ev, s| statistic(ev, e, s))
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

/// `Median[l]`, `Median[l, f]`.
pub(super) fn median(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, |ev, _, s| {
        let Some(sorted) = order::Sorted::of(ev, s)? else {
            return Ok(None);
        };
        sorted.median(ev)
    })
}

/// `Quantile[l, p]`, `Quantile[l, f, p]`.
pub(super) fn quantile(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, f, p) = match e.args() {
        [l, p] => (l, None, p),
        [l, f, p] => (l, Some(f), p),
        _ => return Ok(None),
    };
    if !order::is_fraction(ev, p)? {
        return Ok(None);
    }
    quantile_of_samples(ev, l, f, p)
}

/// The quantile `p` of the data `l` with the frequencies `f`.
fn quantile_of_samples(ev: &mut Eval<'_>, l: &Expr, f: Option<&Expr>, p: &Expr) -> Native {
    of_samples(ev, l, f, |ev, s| {
        let Some(sorted) = order::Sorted::of(ev, s)? else {
            return Ok(None);
        };
        sorted.quantile(ev, p)
    })
}

/// The quantile `p` of the data of `e[l]` or `e[l, f]`.
fn quantile_of(ev: &mut Eval<'_>, e: &Expr, p: &Expr) -> Native {
    match e.args() {
        [l] => quantile_of_samples(ev, l, None, p),
        [l, f] => quantile_of_samples(ev, l, Some(f), p),
        _ => Ok(None),
    }
}

/// `Quartiles[l]`, `Quartiles[l, f]`: `{min, Q(1/4), median, Q(3/4), max}`.
pub(super) fn quartiles(ev: &mut Eval<'_>, e: &Expr) -> Native {
    of_data(ev, e, |ev, _, s| {
        let Some(sorted) = order::Sorted::of(ev, s)? else {
            return Ok(None);
        };
        let (Some(first), Some(median), Some(third)) = (
            sorted.quantile(ev, &quarters(1))?,
            sorted.median(ev)?,
            sorted.quantile(ev, &quarters(3))?,
        ) else {
            return Ok(None);
        };
        let (min, max) = sorted.ends();
        Ok(Some(Expr::call(
            Symbol::LIST,
            [min, first, median, third, max],
        )))
    })
}

/// `k/4`.
fn quarters(k: i64) -> Expr {
    Expr::from(Rational::new(k.into(), 4.into()))
}

/// `FirstQuartile[l]`, `FirstQuartile[l, f]`: `Q(1/4)`.
pub(super) fn first_quartile(ev: &mut Eval<'_>, e: &Expr) -> Native {
    quantile_of(ev, e, &quarters(1))
}

/// `ThirdQuartile[l]`, `ThirdQuartile[l, f]`: `Q(3/4)`.
pub(super) fn third_quartile(ev: &mut Eval<'_>, e: &Expr) -> Native {
    quantile_of(ev, e, &quarters(3))
}

/// The two variables of the arguments `[x, y]`, `[x, y, f]`, `[m]` or
/// `[m, f]`, `m` a matrix of two columns: samples of one length, with one
/// frequencies.
fn pairs(args: &[Expr]) -> Option<(Sample, Sample)> {
    let two = |samples: Option<(Vec<Sample>, bool)>| match samples {
        Some((mut columns, true)) if columns.len() == 2 => {
            let y = columns.pop()?;
            let x = columns.pop()?;
            Some((x, y))
        }
        _ => None,
    };
    let one = |l: &Expr, f: Option<&Expr>| match samples(l, f) {
        Some((mut samples, false)) => samples.pop(),
        _ => None,
    };
    match args {
        [m] => two(samples(m, None)),
        [a, b] => two(samples(a, Some(b))).or_else(|| Some((one(a, None)?, one(b, None)?))),
        [x, y, f] => Some((one(x, Some(f))?, one(y, Some(f))?)),
        _ => None,
    }
    .filter(|(x, y)| x.values.len() == y.values.len())
}

/// `S_xy`, `S_xx`, `S_yy` and `Σ f` of the variables of `e`; `None` when
/// they are not two, or the frequencies sum to 0, with a warning.
fn moments(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<[Expr; 4]>, EvalError> {
    let Some((x, y)) = pairs(e.args()) else {
        return Ok(None);
    };
    let Some((total, x_mean)) = total_and_mean(ev, e, &x)? else {
        return Ok(None);
    };
    let Some((_, y_mean)) = total_and_mean(ev, e, &y)? else {
        return Ok(None);
    };
    let deviations =
        |s: &Sample, mean: &Expr| -> Vec<Expr> { s.values.iter().map(|v| v - mean).collect() };
    let (dx, dy) = (deviations(&x, &x_mean), deviations(&y, &y_mean));
    let sum = |a: &[Expr], b: &[Expr]| {
        let terms = a.iter().zip(b).enumerate().map(|(i, (u, v))| {
            let product = u * v;
            match &x.frequencies {
                Some(f) => &f[i] * product,
                None => product,
            }
        });
        Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::PLUS, terms.collect::<Vec<_>>())],
        )
    };
    let (xy, xx, yy) = (sum(&dx, &dy), sum(&dx, &dx), sum(&dy, &dy));
    Ok(Some([ev.eval(&xy)?, ev.eval(&xx)?, ev.eval(&yy)?, total]))
}

/// `PopulationCovariance[x, y]` and its other forms: `S_xy / Σ f`.
pub(super) fn population_covariance(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some([xy, _, _, total]) = moments(ev, e)? else {
        return Ok(None);
    };
    let covariance = Expr::call(Symbol::EXPAND, [xy / total]);
    Ok(Some(ev.eval(&covariance)?))
}

/// `Correlation[x, y]` and its other forms: `S_xy / √(S_xx S_yy)`.
pub(super) fn correlation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some([xy, xx, yy, _]) = moments(ev, e)? else {
        return Ok(None);
    };
    let root = Expr::call(Symbol::SQRT, [xx * yy]);
    Ok(Some(ev.eval(&(xy / root))?))
}

#[cfg(test)]
mod tests;
