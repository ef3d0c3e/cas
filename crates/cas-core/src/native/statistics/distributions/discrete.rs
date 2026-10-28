//! Discrete distributions (CAS-DOC-23 §2, CAS-301): Xcas's binomial,
//! Poisson, geometric and negative binomial distributions.
//!
//! - `Binomial[n, p, k]` (`binomial`, the coefficient's three-argument
//!   form): `C(n, k) pᵏ (1 − p)ⁿ⁻ᵏ` on `0 … n`. `BinomialCDF`
//!   (`binomial_cdf`) and `BinomialICDF` (`binomial_icdf`).
//! - `PoissonPDF[λ, k]` (`poisson`): `e^(−λ) λᵏ / k!` on the naturals;
//!   `PoissonCDF`, `PoissonICDF`.
//! - `GeometricPDF[p, k]` (`geometric`): `p (1 − p)ᵏ⁻¹`, the trials to the
//!   first success, from 1; its cumulative is `1 − (1 − p)ᵏ`.
//!   `GeometricCDF`, `GeometricICDF`.
//! - `NegativeBinomialPDF[n, p, k]` (`negbinomial`): `C(n + k − 1, k) pⁿ
//!   (1 − p)ᵏ`, the failures before the `n`-th success, from 0;
//!   `NegativeBinomialCDF`, `NegativeBinomialICDF`.
//!
//! The mass is 0 outside the support and at a number that is not an
//! integer; at a symbol it is the expression, as giac's. The cumulative
//! `F(x)` sums the masses up to `⌊x⌋`, and `cdf(a, b)` is `P(a ≤ X ≤ b) =
//! F(⌊b⌋) − F(⌈a⌉ − 1)`, as giac's. The inverse at `q` is the least `k`
//! with `F(k) ≥ q`: the least value of the support at 0, and `∞` at 1 for
//! an unbounded support (giac's `+infinity`). Sums and searches take at
//! most [`MAX_TERMS`] masses; past it, the call stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the textbook names: B(n, p), the point k, an interval [a, b]"
)]

use std::cmp::Ordering;

use cas_num::Integer;

use super::{Fit, compare, positive, probability, refuse};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most masses a cumulative sum adds, or an inverse goes through.
const MAX_TERMS: usize = 1 << 12;

/// The four distributions, before their parameters are read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Binomial,
    Poisson,
    Geometric,
    NegativeBinomial,
}

impl Kind {
    /// How many parameters come before the point.
    fn parameters(self) -> usize {
        match self {
            Kind::Binomial | Kind::NegativeBinomial => 2,
            Kind::Poisson | Kind::Geometric => 1,
        }
    }
}

/// A discrete distribution on the integers.
enum Discrete {
    /// `B(n, p)`, on `0 … n`.
    Binomial { n: Integer, p: Expr },
    /// `P(λ)`, on the naturals.
    Poisson { lambda: Expr },
    /// The trials to the first success, from 1.
    Geometric { p: Expr },
    /// The failures before the `n`-th success, from 0.
    NegativeBinomial { n: Integer, p: Expr },
}

/// `b^e`, 1 for an exponent 0 even when `b` is 0.
fn power(b: &Expr, e: Expr) -> Expr {
    if e == Expr::from(0) {
        Expr::from(1)
    } else {
        b.pow(e)
    }
}

/// `a + k` for an integer `a`, folded when `k` is an integer.
fn plus(a: &Integer, k: &Expr) -> Expr {
    match k.to_integer() {
        Some(k) => Expr::from(a + &k),
        None => Expr::from(a.clone()) + k,
    }
}

impl Discrete {
    /// The least value of the support.
    fn low(&self) -> Integer {
        match self {
            Discrete::Geometric { .. } => Integer::from(1),
            _ => Integer::ZERO,
        }
    }

    /// The greatest value of the support; `None` when it is unbounded.
    fn high(&self) -> Option<&Integer> {
        match self {
            Discrete::Binomial { n, .. } => Some(n),
            _ => None,
        }
    }

    /// The probability of success, for those that have one.
    fn success(&self) -> Option<&Expr> {
        match self {
            Discrete::Binomial { p, .. }
            | Discrete::Geometric { p }
            | Discrete::NegativeBinomial { p, .. } => Some(p),
            Discrete::Poisson { .. } => None,
        }
    }

    /// The mass at `k`, an integer of the support or a symbol, unevaluated.
    fn mass(&self, k: &Expr) -> Expr {
        let minus = |a: &Integer, k: &Expr| match k.to_integer() {
            Some(k) => Expr::from(a - &k),
            None => Expr::from(a.clone()) - k,
        };
        match self {
            Discrete::Binomial { n, p } => {
                let coefficient = Expr::call(Symbol::BINOMIAL, [Expr::from(n.clone()), k.clone()]);
                coefficient * power(p, k.clone()) * power(&(1 - p), minus(n, k))
            }
            Discrete::Poisson { lambda } => {
                Expr::call(Symbol::EXP, [-lambda]) * power(lambda, k.clone())
                    / Expr::call(Symbol::FACTORIAL, [k.clone()])
            }
            Discrete::Geometric { p } => p * power(&(1 - p), plus(&Integer::from(-1), k)),
            Discrete::NegativeBinomial { n, p } => {
                let trials = plus(&(n - &Integer::from(1)), k);
                let coefficient = Expr::call(Symbol::BINOMIAL, [trials, k.clone()]);
                coefficient * power(p, Expr::from(n.clone())) * power(&(1 - p), k.clone())
            }
        }
    }

    /// `F(m)` for an integer `m`, evaluated; `None` past [`MAX_TERMS`].
    fn cumulative(&self, ev: &mut Eval<'_>, m: &Integer) -> Result<Option<Expr>, EvalError> {
        let low = self.low();
        if *m < low {
            return Ok(Some(Expr::from(0)));
        }
        if self.high().is_some_and(|high| m >= high) {
            return Ok(Some(Expr::from(1)));
        }
        if let Discrete::Geometric { p } = self {
            let f = 1 - power(&(1 - p), Expr::from(m.clone()));
            return Ok(Some(ev.eval(&f)?));
        }
        let count = (m - &low).to_i64().and_then(|c| usize::try_from(c).ok());
        let Some(count) = count.filter(|&c| c < MAX_TERMS) else {
            return Ok(None);
        };
        let terms: Vec<Expr> = (0..=count)
            .map(|j| self.mass(&Expr::from(&low + &Integer::from(j))))
            .collect();
        Ok(Some(ev.eval(&Expr::call(Symbol::PLUS, terms))?))
    }

    /// The least `k` with `F(k) ≥ q` for a number `q` in `(0, 1)`, by
    /// adding the masses; `None` past [`MAX_TERMS`] or when a comparison is
    /// not decided.
    fn search(&self, ev: &mut Eval<'_>, q: &Expr) -> Native {
        let mut k = self.low();
        let mut sum = Expr::from(0);
        for _ in 0..MAX_TERMS {
            sum = ev.eval(&(sum + self.mass(&Expr::from(k.clone()))))?;
            match compare(ev, &sum, q)? {
                Some(Ordering::Greater | Ordering::Equal) => return Ok(Some(Expr::from(k))),
                Some(Ordering::Less) => {}
                None => return Ok(None),
            }
            if self.high().is_some_and(|high| k >= *high) {
                return Ok(None);
            }
            k = &k + &Integer::from(1);
        }
        Ok(None)
    }

    /// The least `k` with `F(k) ≥ q` for a number `q` in `[0, 1]`.
    fn inverse(&self, ev: &mut Eval<'_>, q: &Expr) -> Native {
        let one = Expr::from(1);
        if compare(ev, q, &Expr::from(0))? == Some(Ordering::Equal) {
            return Ok(Some(Expr::from(self.low())));
        }
        // A success that is certain makes the least value of the support
        // certain too, for the geometric and the negative binomial.
        let certain = match self.success() {
            Some(p) => compare(ev, p, &one)? == Some(Ordering::Equal),
            None => false,
        };
        if compare(ev, q, &one)? == Some(Ordering::Equal) && !certain {
            // `F(k) < 1` below the greatest value of positive mass: `n`
            // for a binomial with `p > 0`, 0 for `p = 0`; none for the
            // others.
            return Ok(Some(match (self.high(), self.success()) {
                (Some(high), Some(p)) => {
                    if compare(ev, p, &Expr::from(0))? == Some(Ordering::Equal) {
                        Expr::from(0)
                    } else {
                        Expr::from(high.clone())
                    }
                }
                _ => Expr::from(Symbol::INFINITY),
            }));
        }
        if let Discrete::Geometric { p } = self {
            if certain {
                return Ok(Some(Expr::from(1)));
            }
            // The least k ≥ 1 with (1 − p)ᵏ ≤ 1 − q.
            let ratio = Expr::call(Symbol::LOG, [1 - q]) / Expr::call(Symbol::LOG, [1 - p]);
            let k = ev.eval(&Expr::call(Symbol::CEILING, [ratio]))?;
            if let Some(k) = k.to_integer() {
                return Ok(Some(Expr::from(k.max(Integer::from(1)))));
            }
        }
        self.search(ev, q)
    }
}

/// The distribution of `kind` with the parameters `params`, and whether
/// they are all numbers; `None` when they do not fit, with a warning for a
/// number that does not.
fn distribution(
    ev: &mut Eval<'_>,
    e: &Expr,
    kind: Kind,
    params: &[Expr],
) -> Result<Option<(Discrete, bool)>, EvalError> {
    let count = |n: &Expr, least: i64| n.to_integer().filter(|n| *n >= Integer::from(least));
    let (distribution, fit) = match (kind, params) {
        (Kind::Binomial | Kind::NegativeBinomial, [n, p]) => {
            // A negative binomial waits for at least one success.
            let least = i64::from(kind == Kind::NegativeBinomial);
            let Some(n) = count(n, least) else {
                if n.is_number() {
                    refuse(ev, e, "the count must be a natural number")?;
                }
                return Ok(None);
            };
            let fit = probability(ev, p)?;
            let p = p.clone();
            let distribution = if kind == Kind::Binomial {
                Discrete::Binomial { n, p }
            } else {
                Discrete::NegativeBinomial { n, p }
            };
            (distribution, fit)
        }
        (Kind::Poisson, [lambda]) => (
            Discrete::Poisson {
                lambda: lambda.clone(),
            },
            positive(ev, lambda)?,
        ),
        (Kind::Geometric, [p]) => {
            let fit = match probability(ev, p)? {
                // A success that never comes leaves no distribution.
                Fit::Number if compare(ev, p, &Expr::from(0))? != Some(Ordering::Greater) => {
                    Fit::Wrong
                }
                fit => fit,
            };
            (Discrete::Geometric { p: p.clone() }, fit)
        }
        _ => return Ok(None),
    };
    match fit {
        Fit::Wrong => {
            let why = match kind {
                Kind::Poisson => "the rate must be above 0",
                Kind::Geometric => "the probability must be in (0, 1]",
                _ => "the probability must be in [0, 1]",
            };
            refuse(ev, e, why)?;
            Ok(None)
        }
        Fit::Number => Ok(Some((distribution, true))),
        Fit::Symbolic => Ok(Some((distribution, false))),
    }
}

/// The mass of `kind` at the last argument of `e`.
fn mass(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let Some((k, params)) = e.args().split_last() else {
        return Ok(None);
    };
    if params.len() != kind.parameters() {
        return Ok(None);
    }
    let Some((d, _)) = distribution(ev, e, kind, params)? else {
        return Ok(None);
    };
    if k.is_number() {
        let in_support = k
            .to_integer()
            .filter(|k| *k >= d.low() && d.high().is_none_or(|high| k <= high));
        return Ok(Some(match in_support {
            Some(_) => ev.eval(&d.mass(k))?,
            None => Expr::from(0),
        }));
    }
    Ok(Some(ev.eval(&d.mass(k))?))
}

/// `⌊x⌋`, or `⌈x⌉` when `up`, of a real number `x`; `None` otherwise.
fn rounded(ev: &mut Eval<'_>, x: &Expr, up: bool) -> Result<Option<Integer>, EvalError> {
    let head = if up { Symbol::CEILING } else { Symbol::FLOOR };
    Ok(ev.eval(&Expr::call(head, [x.clone()]))?.to_integer())
}

/// `F(x)` for a real `x`, `±∞` included.
fn cumulative_at(ev: &mut Eval<'_>, d: &Discrete, x: &Expr) -> Native {
    let infinity = Expr::from(Symbol::INFINITY);
    if *x == infinity {
        return Ok(Some(Expr::from(1)));
    }
    if *x == ev.eval(&-infinity)? {
        return Ok(Some(Expr::from(0)));
    }
    let Some(m) = rounded(ev, x, false)? else {
        return Ok(None);
    };
    d.cumulative(ev, &m)
}

/// The cumulative probability of `kind` at the last argument of `e`, or
/// of the interval of its last two.
fn cumulative(ev: &mut Eval<'_>, e: &Expr, kind: Kind) -> Native {
    let n = kind.parameters();
    let args = e.args();
    if args.len() != n + 1 && args.len() != n + 2 {
        return Ok(None);
    }
    let Some((d, _)) = distribution(ev, e, kind, &args[..n])? else {
        return Ok(None);
    };
    match &args[n..] {
        [x] => cumulative_at(ev, &d, x),
        [a, b] => {
            let Some(below) = rounded(ev, a, true)? else {
                return Ok(None);
            };
            let below = Expr::from(&below - &Integer::from(1));
            let (Some(high), Some(low)) =
                (cumulative_at(ev, &d, b)?, cumulative_at(ev, &d, &below)?)
            else {
                return Ok(None);
            };
            let difference = ev.eval(&(high - low))?;
            // An empty interval: nothing between its ends.
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
    if params.len() != kind.parameters() {
        return Ok(None);
    }
    let Some((d, true)) = distribution(ev, e, kind, params)? else {
        return Ok(None);
    };
    match probability(ev, q)? {
        Fit::Number => d.inverse(ev, q),
        Fit::Wrong => refuse(ev, e, "the probability must be in [0, 1]"),
        Fit::Symbolic => Ok(None),
    }
}

/// `Binomial[n, p, k]`, the mass: from `Binomial`'s native.
pub(in crate::native) fn binomial_mass(ev: &mut Eval<'_>, e: &Expr) -> Native {
    mass(ev, e, Kind::Binomial)
}

/// `BinomialCDF[n, p, k]`, `BinomialCDF[n, p, a, b]`.
pub(in crate::native) fn binomial_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Binomial)
}

/// `BinomialICDF[n, p, q]`.
pub(in crate::native) fn binomial_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Binomial)
}

/// `PoissonPDF[λ, k]`.
pub(in crate::native) fn poisson_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    mass(ev, e, Kind::Poisson)
}

/// `PoissonCDF[λ, k]`, `PoissonCDF[λ, a, b]`.
pub(in crate::native) fn poisson_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Poisson)
}

/// `PoissonICDF[λ, q]`.
pub(in crate::native) fn poisson_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Poisson)
}

/// `GeometricPDF[p, k]`.
pub(in crate::native) fn geometric_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    mass(ev, e, Kind::Geometric)
}

/// `GeometricCDF[p, k]`, `GeometricCDF[p, a, b]`.
pub(in crate::native) fn geometric_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::Geometric)
}

/// `GeometricICDF[p, q]`.
pub(in crate::native) fn geometric_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::Geometric)
}

/// `NegativeBinomialPDF[n, p, k]`.
pub(in crate::native) fn negative_binomial_pdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    mass(ev, e, Kind::NegativeBinomial)
}

/// `NegativeBinomialCDF[n, p, k]`, `NegativeBinomialCDF[n, p, a, b]`.
pub(in crate::native) fn negative_binomial_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    cumulative(ev, e, Kind::NegativeBinomial)
}

/// `NegativeBinomialICDF[n, p, q]`.
pub(in crate::native) fn negative_binomial_icdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    inverse(ev, e, Kind::NegativeBinomial)
}
