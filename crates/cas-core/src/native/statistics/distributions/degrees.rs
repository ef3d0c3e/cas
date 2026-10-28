//! Student's t, χ² and Fisher's F (CAS-DOC-23 §4, CAS-304), for natural
//! degrees of freedom: exact densities through `Gamma`, exact cumulative
//! functions through the closed forms of [`incomplete`](super::incomplete),
//! and inverses as certified decimals by the monotone inversion
//! ([`invert`]).
//!
//! - `StudentTPDF[ν, t]` (`student`): `Γ((ν+1)/2) / (√(νπ) Γ(ν/2)) (1 +
//!   t²/ν)^(−(ν+1)/2)`. `StudentTCDF` (`student_cdf`): `1/2 ± I_{t²/(ν+t²)}
//!   (1/2, ν/2)/2` by the sign of `t`, `3/4` at `(1, 1)` where giac answers
//!   `0.75`. `StudentTICDF` (`student_icdf`): 0 at 1/2 (giac:
//!   `-2.05e-19`), and by symmetry below it.
//! - `ChiSquarePDF[ν, x]` (`chisquare`): `x^(ν/2−1) e^(−x/2) / (2^(ν/2)
//!   Γ(ν/2))` for `x ≥ 0`. `ChiSquareCDF`: `P(ν/2, x/2)`, `1 − 2/e` at
//!   `(4, 2)`.
//! - `FisherPDF[ν₁, ν₂, x]` (`fisher`): `(ν₁/ν₂)^(ν₁/2) x^(ν₁/2−1) (1 +
//!   ν₁x/ν₂)^(−(ν₁+ν₂)/2) / B(ν₁/2, ν₂/2)` for `x ≥ 0`. `FisherCDF`:
//!   `I_{ν₁x/(ν₁x+ν₂)}(ν₁/2, ν₂/2)`, where giac leaves `Beta(…, 1)`.
//!
//! Degrees of freedom that are numbers but not natural numbers above 0
//! leave the call: those that are not positive with a warning, the others
//! as functions this plan leaves out. The densities take a symbolic
//! point, as giac's do. The cumulative functions and inverses need a real
//! number (or `±∞`, and `q` in `[0, 1]`).

#![allow(
    clippy::many_single_char_names,
    reason = "the textbook names: ν as n, m, a point x or t, an interval [a, b]"
)]

use std::cmp::Ordering;

use cas_num::Rational;

use super::incomplete::{beta, lower_gamma};
use super::inversion::invert;
use super::{Fit, compare, probability, refuse};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The three distributions, before their degrees of freedom are read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Student,
    ChiSquare,
    Fisher,
}

impl Kind {
    fn degrees(self) -> usize {
        match self {
            Kind::Student | Kind::ChiSquare => 1,
            Kind::Fisher => 2,
        }
    }
}

/// A distribution with its degrees of freedom.
#[derive(Clone, Copy, Debug)]
enum Degrees {
    Student(u64),
    ChiSquare(u64),
    Fisher(u64, u64),
}

/// `n` as an expression.
fn int(n: u64) -> Expr {
    Expr::from(cas_num::Integer::from(n))
}

/// `k/2`.
fn half(k: u64) -> Expr {
    Expr::from(Rational::new(k.into(), 2.into()))
}

impl Degrees {
    /// The least value of the support; `None` for all of the line.
    fn low(self) -> Option<Rational> {
        match self {
            Degrees::Student(_) => None,
            _ => Some(Rational::ZERO),
        }
    }

    /// The density at `x`, unevaluated.
    fn density(self, x: &Expr) -> Expr {
        match self {
            Degrees::Student(n) => {
                let scale = Expr::call(Symbol::GAMMA, [half(n + 1)])
                    / (Expr::call(Symbol::SQRT, [int(n) * Expr::from(Symbol::PI)])
                        * Expr::call(Symbol::GAMMA, [half(n)]));
                scale * (1 + x.pow(2) / int(n)).pow(-half(n + 1))
            }
            Degrees::ChiSquare(n) => {
                let power = if n == 2 {
                    Expr::from(1)
                } else {
                    x.pow(half(n) - 1)
                };
                power * Expr::call(Symbol::EXP, [-x / 2])
                    / (Expr::from(2).pow(half(n)) * Expr::call(Symbol::GAMMA, [half(n)]))
            }
            Degrees::Fisher(m, n) => {
                let ratio = Expr::from(Rational::new(m.into(), n.into()));
                let power = if m == 2 {
                    Expr::from(1)
                } else {
                    x.pow(half(m) - 1)
                };
                ratio.pow(half(m)) * power * (1 + &ratio * x).pow(-half(m + n))
                    / Expr::call(Symbol::BETA, [half(m), half(n)])
            }
        }
    }

    /// `F(x)` for `x` in the support, `x ≥ 0` for Student's, unevaluated;
    /// `None` without a closed form.
    fn cumulative(self, x: &Expr) -> Option<Expr> {
        match self {
            Degrees::Student(n) => {
                let square = x.pow(2);
                let at = &square / (int(n) + &square);
                Some(Expr::from(Rational::new(1.into(), 2.into())) + beta(1, n, &at)? / 2)
            }
            Degrees::ChiSquare(n) => lower_gamma(n, &(x / 2)),
            Degrees::Fisher(m, n) => {
                let scaled = int(m) * x;
                beta(m, n, &(&scaled / (&scaled + int(n))))
            }
        }
    }

    /// `F(x)` for a real `x`; `None` when its place is not decided.
    fn cumulative_at(self, ev: &mut Eval<'_>, x: &Expr) -> Native {
        let infinity = Expr::from(Symbol::INFINITY);
        if *x == infinity {
            return Ok(Some(Expr::from(1)));
        }
        if *x == ev.eval(&-infinity)? {
            return Ok(Some(Expr::from(0)));
        }
        let sign = compare(ev, x, &Expr::from(0))?;
        let value = match (self, sign) {
            (_, None) => return Ok(None),
            (Degrees::Student(_), Some(Ordering::Equal)) => {
                Some(Expr::from(Rational::new(1.into(), 2.into())))
            }
            (Degrees::Student(_), Some(Ordering::Less)) => self.cumulative(&-x).map(|f| 1 - f),
            (_, Some(Ordering::Less | Ordering::Equal)) if self.low().is_some() => {
                Some(Expr::from(0))
            }
            _ => self.cumulative(x),
        };
        // Expanded: the closed forms nest their sums.
        match value {
            Some(v) => Ok(Some(ev.eval(&Expr::call(Symbol::EXPAND, [v]))?)),
            None => Ok(None),
        }
    }

    /// `F⁻¹(q)` for a number `q` in `[0, 1]`.
    fn inverse(self, ev: &mut Eval<'_>, q: &Expr) -> Native {
        let infinity = Expr::from(Symbol::INFINITY);
        if compare(ev, q, &Expr::from(0))? == Some(Ordering::Equal) {
            return Ok(Some(match self.low() {
                Some(low) => Expr::from(low),
                None => ev.eval(&-infinity)?,
            }));
        }
        if compare(ev, q, &Expr::from(1))? == Some(Ordering::Equal) {
            return Ok(Some(infinity));
        }
        let half = Expr::from(Rational::new(1.into(), 2.into()));
        // Student's below 1/2 by symmetry, on t ≥ 0.
        let (target, mirrored) = match (self, compare(ev, q, &half)?) {
            (Degrees::Student(_), Some(Ordering::Equal)) => return Ok(Some(Expr::from(0))),
            (Degrees::Student(_), Some(Ordering::Less)) => (ev.eval(&(1 - q))?, true),
            (Degrees::Student(_), None) => return Ok(None),
            _ => (q.clone(), false),
        };
        let t = Expr::symbol("Distribution$t");
        let Some(f) = self.cumulative(&t) else {
            return Ok(None);
        };
        let low = self.low().unwrap_or(Rational::ZERO);
        let bits = ev.context().working_precision();
        let Some(root) = invert(ev.context(), &f, &t, &target, Some(&low), bits) else {
            return Ok(None);
        };
        Ok(Some(if mirrored { ev.eval(&-root)? } else { root }))
    }
}

/// The distribution of `kind` with the degrees of freedom `params`; `None`
/// when they are not natural numbers above 0, with a warning for those not
/// above 0.
fn distribution(
    ev: &mut Eval<'_>,
    e: &Expr,
    kind: Kind,
    params: &[Expr],
) -> Result<Option<Degrees>, EvalError> {
    let mut degrees = Vec::with_capacity(params.len());
    for p in params {
        if p.is_number() && compare(ev, p, &Expr::from(0))? != Some(Ordering::Greater) {
            refuse(ev, e, "the degrees of freedom must be above 0")?;
            return Ok(None);
        }
        let Some(n) = p.to_integer().and_then(|n| u64::try_from(n.to_i64()?).ok()) else {
            return Ok(None);
        };
        degrees.push(n);
    }
    Ok(Some(match (kind, degrees.as_slice()) {
        (Kind::Student, [n]) => Degrees::Student(*n),
        (Kind::ChiSquare, [n]) => Degrees::ChiSquare(*n),
        (Kind::Fisher, [m, n]) => Degrees::Fisher(*m, *n),
        _ => return Ok(None),
    }))
}

/// The density of `kind` at the last argument of `e`.
fn density(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((x, params)) = e.args().split_last() else {
        return Ok(None);
    };
    if params.len() != kind.degrees() {
        return Ok(None);
    }
    let Some(d) = distribution(ev, e, kind, params)? else {
        return Ok(None);
    };
    if d.low().is_some() && compare(ev, x, &Expr::from(0))? == Some(Ordering::Less) {
        return Ok(Some(Expr::from(0)));
    }
    Ok(Some(ev.eval(&d.density(x))?))
}

/// The cumulative probability of `kind` at the last argument of `e`, or
/// of the interval of the last two.
fn cumulative(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let n = kind.degrees();
    let args = e.args();
    if args.len() != n + 1 && args.len() != n + 2 {
        return Ok(None);
    }
    let Some(d) = distribution(ev, e, kind, &args[..n])? else {
        return Ok(None);
    };
    match &args[n..] {
        [x] => d.cumulative_at(ev, x),
        [a, b] => {
            let (Some(high), Some(low)) = (d.cumulative_at(ev, b)?, d.cumulative_at(ev, a)?) else {
                return Ok(None);
            };
            let difference = ev.eval(&(high - low))?;
            Ok(Some(match compare(ev, &difference, &Expr::from(0))? {
                Some(Ordering::Less) => Expr::from(0),
                _ => difference,
            }))
        }
        _ => Ok(None),
    }
}

/// The inverse of `kind` at the last argument of `e`.
fn inverse(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((q, params)) = e.args().split_last() else {
        return Ok(None);
    };
    if params.len() != kind.degrees() {
        return Ok(None);
    }
    let Some(d) = distribution(ev, e, kind, params)? else {
        return Ok(None);
    };
    match probability(ev, q)? {
        Fit::Number => d.inverse(ev, q),
        Fit::Wrong => refuse(ev, e, "the probability must be in [0, 1]"),
        Fit::Symbolic => Ok(None),
    }
}

/// `StudentTPDF[ν, t]`.
pub(in crate::native) fn student_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Student)
}

/// `StudentTCDF[ν, t]`, `StudentTCDF[ν, a, b]`.
pub(in crate::native) fn student_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Student)
}

/// `StudentTICDF[ν, q]`.
pub(in crate::native) fn student_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Student)
}

/// `ChiSquarePDF[ν, x]`.
pub(in crate::native) fn chi_square_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::ChiSquare)
}

/// `ChiSquareCDF[ν, x]`, `ChiSquareCDF[ν, a, b]`.
pub(in crate::native) fn chi_square_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::ChiSquare)
}

/// `ChiSquareICDF[ν, q]`.
pub(in crate::native) fn chi_square_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::ChiSquare)
}

/// `FisherPDF[ν₁, ν₂, x]`.
pub(in crate::native) fn fisher_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    density(ev, e, Kind::Fisher)
}

/// `FisherCDF[ν₁, ν₂, x]`, `FisherCDF[ν₁, ν₂, a, b]`.
pub(in crate::native) fn fisher_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Fisher)
}

/// `FisherICDF[ν₁, ν₂, q]`.
pub(in crate::native) fn fisher_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Fisher)
}
