//! The gamma, beta, Weibull and Cauchy distributions (CAS-DOC-46 §4,
//! CAS-470): Xcas's `gammad`, `betad`, `weibulld` and `cauchyd` with their
//! `_cdf` and `_icdf`.
//!
//! - Gamma, shape `a` and rate `b` (giac's): `bᵃ x^(a−1) e^(−bx)/Γ(a)` for
//!   `x ≥ 0`; `F(x) = P(a, bx)`, `LowerGamma[a, b x, 1]`, giac's
//!   `igamma(2,3x,1)`.
//! - Beta: `x^(a−1)(1−x)^(b−1)/B(a, b)` on `[0, 1]`; `F(x) = I_x(a, b)`,
//!   `Beta[a, b, x, 1]`.
//! - Weibull, shape `k`, scale `λ` and location `θ` (0 when left out):
//!   `(k/λ) ((x−θ)/λ)^(k−1) e^(−((x−θ)/λ)^k)` for `x ≥ θ`; `F(x) = 1 −
//!   e^(−((x−θ)/λ)^k)`; `F⁻¹(q) = θ + λ (−log(1 − q))^(1/k)`.
//! - Cauchy, location `x₀` and scale `γ` (0 and 1 when left out): `γ/(π((x
//!   − x₀)² + γ²))`; `F(x) = 1/2 + atan((x − x₀)/γ)/π`; `F⁻¹(q) = x₀ + γ
//!   tan(π(q − 1/2))`, giac's `cauchyd_icdf(1,2,0.75)` being 3.
//!
//! The gamma and beta inverses have no closed form: CAS-DOC-23's monotone
//! inversion ([`invert`]) for a number `q`, a ball at the working
//! precision. Outside the support the density is 0 and `F` is 0 or 1, as
//! comparisons decide (giac's formulas go on there). Intervals: `F(b) −
//! F(a)`, for the arities that leave no doubt (giac's: the Weibull's
//! needs its `θ`).

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::inversion::invert;
use super::{Fit, compare, positive, probability, refuse};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The four distributions, before their parameters are read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Gamma,
    Beta,
    Weibull,
    Cauchy,
}

impl Kind {
    /// Whether it takes `n` parameters.
    fn takes(self, n: usize) -> bool {
        match self {
            Kind::Gamma | Kind::Beta => n == 2,
            Kind::Weibull => n == 2 || n == 3,
            Kind::Cauchy => n == 0 || n == 2,
        }
    }
}

/// A distribution with its parameters.
enum Shape {
    Gamma { a: Expr, b: Expr },
    Beta { a: Expr, b: Expr },
    Weibull { k: Expr, lambda: Expr, theta: Expr },
    Cauchy { x0: Expr, gamma: Expr },
}

/// Where a point lies against the support.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Place {
    Below,
    Inside,
    Above,
    Unknown,
}

fn half() -> Expr {
    Expr::from(Rational::new(Integer::ONE, Integer::from(2)))
}

impl Shape {
    fn place(&self, ev: &mut Eval<'_>, x: &Expr) -> Result<Place, EvalError> {
        let low = match self {
            Shape::Gamma { .. } | Shape::Beta { .. } => Expr::from(0),
            Shape::Weibull { theta, .. } => theta.clone(),
            Shape::Cauchy { .. } => return Ok(Place::Inside),
        };
        let below = compare(ev, x, &low)?;
        let above = match self {
            Shape::Beta { .. } => compare(ev, x, &Expr::from(1))?,
            _ => Some(Ordering::Less),
        };
        Ok(match (below, above) {
            (Some(Ordering::Less), _) => Place::Below,
            (_, Some(Ordering::Greater)) => Place::Above,
            (Some(_), Some(_)) => Place::Inside,
            _ => Place::Unknown,
        })
    }

    fn density(&self, ev: &mut Eval<'_>, x: &Expr) -> Native {
        if matches!(self.place(ev, x)?, Place::Below | Place::Above) {
            return Ok(Some(Expr::from(0)));
        }
        let value = match self {
            Shape::Gamma { a, b } => {
                b.pow(a.clone()) * x.pow(a - 1) * Expr::call(Symbol::EXP, [-(b * x)])
                    / Expr::call(Symbol::GAMMA, [a.clone()])
            }
            Shape::Beta { a, b } => {
                x.pow(a - 1) * (1 - x).pow(b - 1) / Expr::call(Symbol::BETA, [a.clone(), b.clone()])
            }
            Shape::Weibull { k, lambda, theta } => {
                let z = (x - theta) / lambda;
                k / lambda * z.pow(k - 1) * Expr::call(Symbol::EXP, [-(z.pow(k.clone()))])
            }
            Shape::Cauchy { x0, gamma } => {
                gamma / (Expr::from(Symbol::PI) * ((x - x0).pow(2) + gamma.pow(2)))
            }
        };
        Ok(Some(ev.eval(&value)?))
    }

    fn cumulative(&self, ev: &mut Eval<'_>, x: &Expr) -> Native {
        match self.place(ev, x)? {
            Place::Below => return Ok(Some(Expr::from(0))),
            Place::Above => return Ok(Some(Expr::from(1))),
            Place::Inside | Place::Unknown => {}
        }
        let value = self.cumulative_formula(x);
        Ok(Some(ev.eval(&value)?))
    }

    /// `F(x)` in the support, unevaluated.
    fn cumulative_formula(&self, x: &Expr) -> Expr {
        match self {
            Shape::Gamma { a, b } => {
                Expr::call(Symbol::LOWER_GAMMA, [a.clone(), b * x, Expr::from(1)])
            }
            Shape::Beta { a, b } => Expr::call(
                Symbol::BETA,
                [a.clone(), b.clone(), x.clone(), Expr::from(1)],
            ),
            Shape::Weibull { k, lambda, theta } => {
                1 - Expr::call(Symbol::EXP, [-(((x - theta) / lambda).pow(k.clone()))])
            }
            Shape::Cauchy { x0, gamma } => {
                half() + Expr::call(Symbol::ARC_TAN, [(x - x0) / gamma]) / Expr::from(Symbol::PI)
            }
        }
    }

    fn inverse(&self, ev: &mut Eval<'_>, q: &Expr) -> Native {
        let infinity = Expr::from(Symbol::INFINITY);
        let ends = match self {
            Shape::Gamma { .. } => (Expr::from(0), infinity),
            Shape::Beta { .. } => (Expr::from(0), Expr::from(1)),
            Shape::Weibull { theta, .. } => (theta.clone(), infinity),
            Shape::Cauchy { .. } => (-infinity.clone(), infinity),
        };
        if compare(ev, q, &Expr::from(0))? == Some(Ordering::Equal) {
            return Ok(Some(ev.eval(&ends.0)?));
        }
        if compare(ev, q, &Expr::from(1))? == Some(Ordering::Equal) {
            return Ok(Some(ev.eval(&ends.1)?));
        }
        let value = match self {
            Shape::Weibull { k, lambda, theta } => {
                theta + lambda * (-Expr::call(Symbol::LOG, [1 - q])).pow(1 / k.clone())
            }
            Shape::Cauchy { x0, gamma } => {
                x0 + gamma * Expr::call(Symbol::TAN, [Expr::from(Symbol::PI) * (q - half())])
            }
            Shape::Gamma { .. } | Shape::Beta { .. } => {
                // A number: the monotone inversion, on a support from 0.
                if !q.is_number() {
                    return Ok(None);
                }
                let numbers = match self {
                    Shape::Gamma { a, b } | Shape::Beta { a, b } => a.is_number() && b.is_number(),
                    _ => false,
                };
                if !numbers {
                    return Ok(None);
                }
                let t = Expr::symbol("Distribution$t");
                let f = self.cumulative_formula(&t);
                let bits = ev.context().working_precision();
                return Ok(invert(ev.context(), &f, &t, q, Some(&Rational::ZERO), bits));
            }
        };
        Ok(Some(ev.eval(&value)?))
    }
}

/// The distribution of `kind` with `params`; `None` when they do not fit,
/// with a warning for numbers that do not.
fn distribution(
    ev: &mut Eval<'_>,
    e: &Expr,
    kind: Kind,
    params: &[Expr],
) -> Result<Option<Shape>, EvalError> {
    let both = |ev: &mut Eval<'_>, x: &Expr, y: &Expr| -> Result<Fit, EvalError> {
        Ok(match (positive(ev, x)?, positive(ev, y)?) {
            (Fit::Wrong, _) | (_, Fit::Wrong) => Fit::Wrong,
            (Fit::Number, Fit::Number) => Fit::Number,
            _ => Fit::Symbolic,
        })
    };
    let (shape, fit, why) = match (kind, params) {
        (Kind::Gamma, [a, b]) => (
            Shape::Gamma {
                a: a.clone(),
                b: b.clone(),
            },
            both(ev, a, b)?,
            "the shape and the rate must be above 0",
        ),
        (Kind::Beta, [a, b]) => (
            Shape::Beta {
                a: a.clone(),
                b: b.clone(),
            },
            both(ev, a, b)?,
            "the shapes must be above 0",
        ),
        (Kind::Weibull, [k, lambda, rest @ ..]) => (
            Shape::Weibull {
                k: k.clone(),
                lambda: lambda.clone(),
                theta: rest.first().cloned().unwrap_or_else(|| Expr::from(0)),
            },
            both(ev, k, lambda)?,
            "the shape and the scale must be above 0",
        ),
        (Kind::Cauchy, []) => (
            Shape::Cauchy {
                x0: Expr::from(0),
                gamma: Expr::from(1),
            },
            Fit::Number,
            "",
        ),
        (Kind::Cauchy, [x0, gamma]) => (
            Shape::Cauchy {
                x0: x0.clone(),
                gamma: gamma.clone(),
            },
            positive(ev, gamma)?,
            "the scale must be above 0",
        ),
        _ => return Ok(None),
    };
    if fit == Fit::Wrong {
        refuse(ev, e, why)?;
        return Ok(None);
    }
    Ok(Some(shape))
}

/// The parameters and the `points` last arguments, when the arity fits.
fn split(kind: Kind, args: &[Expr], points: usize) -> Option<(&[Expr], &[Expr])> {
    let n = args.len().checked_sub(points)?;
    kind.takes(n).then(|| args.split_at(n))
}

fn density(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((params, [x])) = split(kind, e.args(), 1) else {
        return Ok(None);
    };
    match distribution(ev, e, kind, params)? {
        Some(d) => d.density(ev, x),
        None => Ok(None),
    }
}

fn cumulative(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let args = e.args();
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

/// `GammaPDF[a, b, x]`.
pub(in crate::native) fn gamma_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Gamma)
}

/// `GammaCDF[a, b, x]`, `GammaCDF[a, b, x₁, x₂]`.
pub(in crate::native) fn gamma_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Gamma)
}

/// `GammaICDF[a, b, q]`.
pub(in crate::native) fn gamma_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Gamma)
}

/// `BetaPDF[a, b, x]`.
pub(in crate::native) fn beta_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Beta)
}

/// `BetaCDF[a, b, x]`, `BetaCDF[a, b, x₁, x₂]`.
pub(in crate::native) fn beta_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Beta)
}

/// `BetaICDF[a, b, q]`.
pub(in crate::native) fn beta_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Beta)
}

/// `WeibullPDF[k, λ, x]`, `WeibullPDF[k, λ, θ, x]`.
pub(in crate::native) fn weibull_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Weibull)
}

/// `WeibullCDF[k, λ, x]`, `WeibullCDF[k, λ, θ, x]`, `WeibullCDF[k, λ, θ,
/// x₁, x₂]`.
pub(in crate::native) fn weibull_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Weibull)
}

/// `WeibullICDF[k, λ, q]`, `WeibullICDF[k, λ, θ, q]`.
pub(in crate::native) fn weibull_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Weibull)
}

/// `CauchyPDF[x]`, `CauchyPDF[x₀, γ, x]`.
pub(in crate::native) fn cauchy_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Cauchy)
}

/// `CauchyCDF[x]`, `CauchyCDF[a, b]`, `CauchyCDF[x₀, γ, x]`, `CauchyCDF[x₀,
/// γ, a, b]`.
pub(in crate::native) fn cauchy_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Cauchy)
}

/// `CauchyICDF[q]`, `CauchyICDF[x₀, γ, q]`.
pub(in crate::native) fn cauchy_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Cauchy)
}

#[cfg(test)]
mod tests;
