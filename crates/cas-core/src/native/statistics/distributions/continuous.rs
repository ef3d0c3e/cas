//! The elementary continuous distributions (CAS-DOC-23 §3, CAS-303), in
//! closed form.
//!
//! - Normal `N(μ, σ)`: `NormalPDF` (`normald`), `e^(−(x − μ)²/(2σ²)) /
//!   (σ √(2π))`; `NormalCDF` (`normald_cdf`), `(1 + erf((x − μ)/(σ √2)))/2`;
//!   `NormalICDF` (`normald_icdf`), `μ + σ √2 erfinv(2q − 1)`, `μ` at 1/2
//!   (giac: `-1.9e-19`). Without `μ` and `σ`, the standard one, giac's forms
//!   told apart by their arity: `normald_cdf(x)`, `normald_cdf(a, b)`,
//!   `normald_cdf(μ, σ, x)`, `normald_cdf(μ, σ, a, b)`.
//! - Exponential `E(λ)`: `ExponentialPDF` (`exponentiald`), `λ e^(−λx)` for
//!   `x ≥ 0`; `ExponentialCDF`, `1 − e^(−λx)`; `ExponentialICDF`,
//!   `−ln(1 − q)/λ`.
//! - Uniform `U(a, b)`: `UniformPDF` (`uniformd`), `1/(b − a)` on `[a, b]`;
//!   `UniformCDF`, `(x − a)/(b − a)` held in `[0, 1]`; `UniformICDF`,
//!   `a + q (b − a)`.
//!
//! Outside the support, decided by comparisons, the density is 0 and the
//! cumulative probability 0 or 1. At a point the comparisons do not place,
//! a symbol, the expression stands, as giac's. `cdf(a, b)` is `F(b) −
//! F(a)`. The inverse at `q` takes a number in `[0, 1]`, the ends of the
//! support at 0 and 1 (giac's `±infinity`), or a symbol, for the formula.

use std::cmp::Ordering;

use super::{Fit, compare, positive, probability, refuse};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The three distributions, before their parameters are read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Normal,
    Exponential,
    Uniform,
}

/// A distribution with its parameters.
enum Elementary {
    Normal { mu: Expr, sigma: Expr },
    Exponential { lambda: Expr },
    Uniform { a: Expr, b: Expr },
}

/// Where a point lies against the support.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Place {
    Below,
    Inside,
    Above,
    /// Not decided: a symbol.
    Unknown,
}

impl Elementary {
    /// Where `x` lies.
    fn place(&self, ev: &mut Eval<'_>, x: &Expr) -> Result<Place, EvalError> {
        Ok(match self {
            Elementary::Normal { .. } => Place::Inside,
            Elementary::Exponential { .. } => match compare(ev, x, &Expr::from(0))? {
                Some(Ordering::Less) => Place::Below,
                Some(_) => Place::Inside,
                None => Place::Unknown,
            },
            Elementary::Uniform { a, b } => match (compare(ev, x, a)?, compare(ev, x, b)?) {
                (Some(Ordering::Less), _) => Place::Below,
                (_, Some(Ordering::Greater)) => Place::Above,
                (Some(_), Some(_)) => Place::Inside,
                _ => Place::Unknown,
            },
        })
    }

    /// The density at `x`.
    fn density(&self, ev: &mut Eval<'_>, x: &Expr) -> Native {
        if matches!(self.place(ev, x)?, Place::Below | Place::Above) {
            return Ok(Some(Expr::from(0)));
        }
        let value = match self {
            Elementary::Normal { mu, sigma } => {
                let z = (x - mu) / sigma;
                let root = Expr::call(Symbol::SQRT, [2 * Expr::from(Symbol::PI)]);
                Expr::call(Symbol::EXP, [-(z.pow(2)) / 2]) / (sigma * root)
            }
            Elementary::Exponential { lambda } => lambda * Expr::call(Symbol::EXP, [-(lambda * x)]),
            Elementary::Uniform { a, b } => 1 / (b - a),
        };
        Ok(Some(ev.eval(&value)?))
    }

    /// `F(x)`, `±∞` included.
    fn cumulative(&self, ev: &mut Eval<'_>, x: &Expr) -> Native {
        match self.place(ev, x)? {
            Place::Below => return Ok(Some(Expr::from(0))),
            Place::Above => return Ok(Some(Expr::from(1))),
            Place::Inside | Place::Unknown => {}
        }
        let value = match self {
            Elementary::Normal { mu, sigma } => {
                let root = Expr::call(Symbol::SQRT, [Expr::from(2)]);
                let z = (x - mu) / (sigma * root);
                (1 + Expr::call(Symbol::ERF, [z])) / 2
            }
            Elementary::Exponential { lambda } => 1 - Expr::call(Symbol::EXP, [-(lambda * x)]),
            Elementary::Uniform { a, b } => (x - a) / (b - a),
        };
        Ok(Some(ev.eval(&value)?))
    }

    /// `F⁻¹(q)`.
    fn inverse(&self, ev: &mut Eval<'_>, q: &Expr) -> Native {
        let infinity = Expr::from(Symbol::INFINITY);
        let ends = match self {
            Elementary::Normal { .. } => (-infinity.clone(), infinity),
            Elementary::Exponential { .. } => (Expr::from(0), infinity),
            Elementary::Uniform { a, b } => (a.clone(), b.clone()),
        };
        if compare(ev, q, &Expr::from(0))? == Some(Ordering::Equal) {
            return Ok(Some(ev.eval(&ends.0)?));
        }
        if compare(ev, q, &Expr::from(1))? == Some(Ordering::Equal) {
            return Ok(Some(ev.eval(&ends.1)?));
        }
        let value = match self {
            Elementary::Normal { mu, sigma } => {
                let root = Expr::call(Symbol::SQRT, [Expr::from(2)]);
                let erfinv = Expr::call(Symbol::INVERSE_ERF, [2 * q - 1]);
                mu + sigma * root * erfinv
            }
            Elementary::Exponential { lambda } => -Expr::call(Symbol::LOG, [1 - q]) / lambda,
            Elementary::Uniform { a, b } => a + q * (b - a),
        };
        Ok(Some(ev.eval(&value)?))
    }
}

/// The distribution of `kind` with the parameters `params`; `None` when
/// they do not fit, with a warning for numbers that do not.
fn distribution(
    ev: &mut Eval<'_>,
    e: &Expr,
    kind: Kind,
    params: &[Expr],
) -> Result<Option<Elementary>, EvalError> {
    let (d, fit, why) = match (kind, params) {
        (Kind::Normal, []) => (
            Elementary::Normal {
                mu: Expr::from(0),
                sigma: Expr::from(1),
            },
            Fit::Number,
            "",
        ),
        (Kind::Normal, [mu, sigma]) => (
            Elementary::Normal {
                mu: mu.clone(),
                sigma: sigma.clone(),
            },
            positive(ev, sigma)?,
            "the standard deviation must be above 0",
        ),
        (Kind::Exponential, [lambda]) => (
            Elementary::Exponential {
                lambda: lambda.clone(),
            },
            positive(ev, lambda)?,
            "the rate must be above 0",
        ),
        (Kind::Uniform, [a, b]) => {
            let fit = match compare(ev, a, b)? {
                Some(Ordering::Less) => Fit::Number,
                Some(_) => Fit::Wrong,
                None => Fit::Symbolic,
            };
            (
                Elementary::Uniform {
                    a: a.clone(),
                    b: b.clone(),
                },
                fit,
                "the interval must have a < b",
            )
        }
        _ => return Ok(None),
    };
    if fit == Fit::Wrong {
        refuse(ev, e, why)?;
        return Ok(None);
    }
    Ok(Some(d))
}

/// How many parameters `kind` takes, given its call's arity and how many
/// points or ends follow: the normal's are optional.
fn split(kind: Kind, args: &[Expr], points: usize) -> Option<(&[Expr], &[Expr])> {
    let n = args.len().checked_sub(points)?;
    let fits = match kind {
        Kind::Normal => n == 0 || n == 2,
        Kind::Exponential => n == 1,
        Kind::Uniform => n == 2,
    };
    fits.then(|| args.split_at(n))
}

/// The density of `kind` at the last argument of `e`.
fn density(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((params, [x])) = split(kind, e.args(), 1) else {
        return Ok(None);
    };
    match distribution(ev, e, kind, params)? {
        Some(d) => d.density(ev, x),
        None => Ok(None),
    }
}

/// The cumulative probability of `kind` at the last argument of `e`, or
/// of the interval of the last two.
fn cumulative(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let args = e.args();
    // A point, or else an interval: the standard normal's `(a, b)` and
    // `(μ, σ, x)` differ in arity.
    for points in [1, 2] {
        let Some((params, ends)) = split(kind, args, points) else {
            continue;
        };
        let Some(d) = distribution(ev, e, kind, params)? else {
            return Ok(None);
        };
        return match ends {
            [x] => d.cumulative(ev, x),
            [a, b] => {
                let (Some(high), Some(low)) = (d.cumulative(ev, b)?, d.cumulative(ev, a)?) else {
                    return Ok(None);
                };
                let difference = ev.eval(&Expr::call(Symbol::EXPAND, [high - low]))?;
                // An empty interval: nothing between its ends.
                Ok(Some(match compare(ev, &difference, &Expr::from(0))? {
                    Some(Ordering::Less) => Expr::from(0),
                    _ => difference,
                }))
            }
            _ => Ok(None),
        };
    }
    Ok(None)
}

/// The inverse of `kind` at the last argument of `e`.
fn inverse(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((params, [q])) = split(kind, e.args(), 1) else {
        return Ok(None);
    };
    let Some(d) = distribution(ev, e, kind, params)? else {
        return Ok(None);
    };
    match probability(ev, q)? {
        Fit::Number | Fit::Symbolic => d.inverse(ev, q),
        Fit::Wrong => refuse(ev, e, "the probability must be in [0, 1]"),
    }
}

/// `NormalPDF[x]`, `NormalPDF[μ, σ, x]`.
pub(in crate::native) fn normal_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Normal)
}

/// `NormalCDF[x]`, `NormalCDF[a, b]`, `NormalCDF[μ, σ, x]`,
/// `NormalCDF[μ, σ, a, b]`.
pub(in crate::native) fn normal_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Normal)
}

/// `NormalICDF[q]`, `NormalICDF[μ, σ, q]`.
pub(in crate::native) fn normal_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Normal)
}

/// `ExponentialPDF[λ, x]`.
pub(in crate::native) fn exponential_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Exponential)
}

/// `ExponentialCDF[λ, x]`, `ExponentialCDF[λ, a, b]`.
pub(in crate::native) fn exponential_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Exponential)
}

/// `ExponentialICDF[λ, q]`.
pub(in crate::native) fn exponential_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Exponential)
}

/// `UniformPDF[a, b, x]`.
pub(in crate::native) fn uniform_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Uniform)
}

/// `UniformCDF[a, b, x]`, `UniformCDF[a, b, c, d]`.
pub(in crate::native) fn uniform_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Uniform)
}

/// `UniformICDF[a, b, q]`.
pub(in crate::native) fn uniform_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Uniform)
}
