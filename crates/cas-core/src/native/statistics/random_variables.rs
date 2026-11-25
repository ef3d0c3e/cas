//! Random variables (CAS-624): Xcas's `randvar` (`random_variable`) and
//! `randmultinomial`, and the objects that `sample` and `randvector` draw
//! from (giac graphtheory.cc `_randvar`, moyal.cc `_randmultinomial`).
//!
//! - **The object.** giac answers a law, `normald(5, 3/2)`, that `sample`
//!   draws from; we take the same course with no new kind of expression.
//!   A named law is the parameter-only call of its density, which stays
//!   unevaluated: `NormalPDF[μ, σ]`, `ExponentialPDF[λ]`, `UniformPDF[a,
//!   b]`, `PoissonPDF[λ]`, `GeometricPDF[p]`, `NegativeBinomialPDF[n, p]`,
//!   `StudentTPDF[ν]`, `ChiSquarePDF[ν]`, `FisherPDF[m, n]`, `GammaPDF[a,
//!   b]`, `BetaPDF[a, b]`, `WeibullPDF[k, λ]`, `CauchyPDF[x₀, γ]`, exactly
//!   these arities. `Binomial[n, p]` is the coefficient, so the binomial
//!   law is `RandomVariable[Binomial, n, p]`; a discrete law is
//!   `RandomVariable[{w…}]`, which draws an index `0 … k − 1`, or
//!   `RandomVariable[{w…}, {v…}]`, which draws a `v`, with real weights
//!   `w ≥ 0`, not all 0 (rationals or decimals). These forms print as Xcas
//!   input and read back; `RandomVariable` leaves them as they are, and
//!   `randvar(X)` of an object is `X`. [`variable`] reads them; [`draw`]
//!   draws one value.
//! - **randvar** (graphtheory.cc:612-866): a law's name alone (`normal`,
//!   `uniform`, `exponential`: the standard law) or with its parameters;
//!   or with `mean=`, `stddev=`, `variance=`, `range=a..b` (or `[a, b]`)
//!   or `a..b`, the last given winning, fitted as giac fits them, exactly
//!   (giac gives decimals): `uniformd(m − √3 s, m + √3 s)`, `gammad(m²/v,
//!   m/v)`, the binomial's `n` rounded, and Weibull's shape by Newton's
//!   method from `max(1, 1/v)` on giac's equation `Γ(1 + 2/k)/Γ(1 + 1/k)
//!   − v/m² − 1 = 0`, whose denominator giac does not square (so the law
//!   found has not the variance asked; ported as it is), in certified
//!   balls. `multinomial, P(, V)`, a list of weights with values or not,
//!   a matrix of `[v, w]` rows, or a function with `range=a..b` between
//!   integers (its weights `f(a) … f(b)` over `a … b`) give discrete laws.
//!   giac's errors leave the call with giac's message as a warning; a
//!   mean or a deviation that is no number leaves it without one.
//! - **randmultinomial** (moyal.cc:713-740): `RandomMultinomial[P]` an
//!   index of `P` (from 0), `RandomMultinomial[P, V]` an element of `V`,
//!   for probabilities `p ≥ 0` that sum to 1 (exactly, or within their
//!   balls for decimals).
//! - **Draws** go through `Eval::draw`: a named law's through its variate,
//!   rational weights exactly (masses over their common denominator, one
//!   integer drawn below their sum), decimal weights in `f64` against `u
//!   Σ w`.

#![allow(
    clippy::many_single_char_names,
    reason = "m, s, v, k, p: giac's names for the moments and parameters"
)]

use std::cmp::Ordering;

use cas_num::ball::Ball;
use cas_num::{Integer, Natural, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::numeric::approximation;
use crate::native::random::unit;
use crate::native::random::variates::{Law, variate};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `gensizeerr`.
const BAD: &str = "Bad Argument Value";
/// giac's "Invalid distribution parameter(s)".
const INVALID: &str = "Invalid distribution parameter(s)";

/// A random variable read from its object.
#[derive(Clone, Copy, Debug)]
pub(in crate::native) enum Variable<'e> {
    /// A named law and its parameters.
    Law(Law, &'e [Expr]),
    /// A discrete law: its weights, and the values drawn, or the indices.
    Discrete {
        weights: &'e [Expr],
        values: Option<&'e [Expr]>,
    },
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// The number of parameters of a law's object.
pub(in crate::native) fn arity(law: Law) -> usize {
    match law {
        Law::Exponential | Law::Poisson | Law::Geometric | Law::StudentT | Law::ChiSquare => 1,
        Law::Normal
        | Law::Uniform
        | Law::Binomial
        | Law::NegativeBinomial
        | Law::Fisher
        | Law::Gamma
        | Law::Beta
        | Law::Weibull
        | Law::Cauchy => 2,
    }
}

/// A real number: a rational or a real ball.
fn real_number(e: &Expr) -> bool {
    e.to_rational().is_some() || matches!(e.view(), View::Ball(_))
}

/// A parameter as giac's positional form takes it: a real number or a
/// symbol.
fn parameter(e: &Expr) -> bool {
    real_number(e) || e.as_symbol().is_some()
}

/// The sign of a weight: a rational's, or a finite real ball's midpoint.
fn weight_sign(w: &Expr) -> Option<Ordering> {
    if let Some(q) = w.to_rational() {
        return Some(q.sign());
    }
    match w.view() {
        View::Ball(b) if b.is_finite() => b.mid().to_f64().partial_cmp(&0.0),
        _ => None,
    }
}

/// Weights of a discrete law: real numbers `≥ 0`, not all 0.
fn weights(ws: &[Expr]) -> bool {
    let signs: Option<Vec<Ordering>> = ws.iter().map(weight_sign).collect();
    signs.is_some_and(|s| s.iter().all(|&o| o != Ordering::Less) && s.contains(&Ordering::Greater))
}

/// The random variable of an object; `None` for another expression.
pub(in crate::native) fn variable(e: &Expr) -> Option<Variable<'_>> {
    let head = e.head()?.as_symbol()?;
    let args = e.args();
    if head == Symbol::RANDOM_VARIABLE {
        return match args {
            [b, n, p]
                if b.as_symbol() == Some(Symbol::BINOMIAL) && parameter(n) && parameter(p) =>
            {
                Some(Variable::Law(Law::Binomial, &args[1..]))
            }
            [w] if is_list(w) && weights(w.args()) => Some(Variable::Discrete {
                weights: w.args(),
                values: None,
            }),
            [w, v]
                if is_list(w)
                    && is_list(v)
                    && v.args().len() == w.args().len()
                    && weights(w.args()) =>
            {
                Some(Variable::Discrete {
                    weights: w.args(),
                    values: Some(v.args()),
                })
            }
            _ => None,
        };
    }
    let law = Law::of_density(head)?;
    (law != Law::Binomial && args.len() == arity(law)).then_some(Variable::Law(law, args))
}

/// One value drawn from `x`; `None` when its law's variate stays (a
/// symbol among the parameters, an invalid one with its warning).
pub(in crate::native) fn draw(ev: &mut Eval<'_>, x: &Variable<'_>) -> Native {
    match *x {
        Variable::Law(law, parameters) => variate(ev, law, parameters),
        Variable::Discrete { weights, values } => {
            let i = pick(ev, weights);
            Ok(Some(match values {
                Some(v) => v[i].clone(),
                None => Expr::from(Integer::from(i as u64)),
            }))
        }
    }
}

/// An index drawn with probabilities proportional to `weights`, valid
/// ones ([`weights`]).
fn pick(ev: &mut Eval<'_>, weights: &[Expr]) -> usize {
    let exact: Option<Vec<Rational>> = weights.iter().map(Expr::to_rational).collect();
    if let Some(q) = exact {
        // Integer masses over the common denominator; the first index
        // whose cumulative mass passes r, drawn below their sum.
        let d = q.iter().fold(Integer::ONE, |d, q| d.lcm(q.denom()));
        let masses: Vec<Integer> = q.iter().map(|q| q.numer() * &(&d / q.denom())).collect();
        let total = masses.iter().fold(Integer::ZERO, |s, m| &s + m);
        let total = Natural::try_from(&total).expect("a positive mass");
        let r = Integer::from(ev.draw(|g| g.below(&total)));
        let mut acc = Integer::ZERO;
        for (i, m) in masses.iter().enumerate() {
            acc = &acc + m;
            if r < acc {
                return i;
            }
        }
        return masses.len() - 1;
    }
    let wp = ev.context().working_precision();
    let w: Vec<f64> = weights
        .iter()
        .map(|w| match w.to_rational() {
            Some(q) => Ball::from_rational(&q, wp).mid().to_f64(),
            None => match w.view() {
                View::Ball(b) => b.mid().to_f64(),
                _ => 0.0,
            },
        })
        .collect();
    let target = unit(ev).mid().to_f64() * w.iter().sum::<f64>();
    let mut acc = 0.0;
    for (i, wi) in w.iter().enumerate() {
        acc += wi;
        if acc > target {
            return i;
        }
    }
    w.iter().rposition(|&wi| wi > 0.0).unwrap_or(0)
}

/// The probabilities of `randmultinomial` and of `randvector`'s
/// `multinomial`: real `p ≥ 0` summing to 1, exactly for rationals, within
/// their balls at `wp` bits for decimals; `values` as many.
pub(in crate::native) fn multinomial<'e>(
    p: &'e [Expr],
    values: Option<&'e [Expr]>,
    wp: u64,
) -> Option<Variable<'e>> {
    if !weights(p) || values.is_some_and(|v| v.len() != p.len()) {
        return None;
    }
    let exact: Option<Vec<Rational>> = p.iter().map(Expr::to_rational).collect();
    let one = if let Some(q) = exact {
        q.iter().fold(Rational::from(Integer::ZERO), |s, q| &s + q) == Rational::from(Integer::ONE)
    } else {
        let mut sum = Ball::from_rational(&Rational::from(Integer::ZERO), wp);
        for x in p {
            let b = match (x.to_rational(), x.view()) {
                (Some(q), _) => Ball::from_rational(&q, wp),
                (None, View::Ball(b)) => b.clone(),
                _ => return None,
            };
            sum = sum.add(&b, wp);
        }
        sum.contains_integer(&Integer::ONE)
    };
    one.then_some(Variable::Discrete { weights: p, values })
}

/// Whether `e` is the symbol `multinomial`.
pub(in crate::native) fn is_multinomial(e: &Expr) -> bool {
    e.as_symbol().is_some_and(|s| s.name() == "multinomial")
}

/// `RandomMultinomial[P]`, `RandomMultinomial[P, V]`.
pub(in crate::native) fn random_multinomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, values) = match e.args() {
        [p] => (p, None),
        [p, v] => (p, Some(v)),
        _ => return Ok(None),
    };
    if !is_list(p) {
        return Ok(None);
    }
    if values.is_some_and(|v| !is_list(v)) {
        return refuse(ev, e, BAD);
    }
    let wp = ev.context().working_precision();
    match multinomial(p.args(), values.map(Expr::args), wp) {
        Some(x) => draw(ev, &x),
        None => refuse(ev, e, BAD),
    }
}

/// What a law's name in `randvar` stands for.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Name {
    Law(Law),
    Multinomial,
}

/// giac's names of the laws (`is_density`), as they arrive: the builtins
/// that are their aliases (C0008).
fn name(e: &Expr) -> Option<Name> {
    if is_multinomial(e) {
        return Some(Name::Multinomial);
    }
    let law = match e.as_symbol()? {
        Symbol::TOGETHER | Symbol::NORMAL_PDF | Symbol::RANDOM_NORMAL => Law::Normal,
        Symbol::EXP | Symbol::EXPONENTIAL_PDF => Law::Exponential,
        Symbol::GAMMA | Symbol::GAMMA_PDF => Law::Gamma,
        Symbol::BETA | Symbol::BETA_PDF => Law::Beta,
        s => Law::of_density(s)?,
    };
    Some(Name::Law(law))
}

/// The object of `law` with `parameters`.
fn object(law: Law, parameters: Vec<Expr>) -> Expr {
    let head = match law {
        Law::Binomial => {
            let mut args = vec![Expr::from(Symbol::BINOMIAL)];
            args.extend(parameters);
            return Expr::call(Symbol::RANDOM_VARIABLE, args);
        }
        Law::Normal => Symbol::NORMAL_PDF,
        Law::Exponential => Symbol::EXPONENTIAL_PDF,
        Law::Uniform => Symbol::UNIFORM_PDF,
        Law::Poisson => Symbol::POISSON_PDF,
        Law::Geometric => Symbol::GEOMETRIC_PDF,
        Law::NegativeBinomial => Symbol::NEGATIVE_BINOMIAL_PDF,
        Law::StudentT => Symbol::STUDENT_T_PDF,
        Law::ChiSquare => Symbol::CHI_SQUARE_PDF,
        Law::Fisher => Symbol::FISHER_PDF,
        Law::Gamma => Symbol::GAMMA_PDF,
        Law::Beta => Symbol::BETA_PDF,
        Law::Weibull => Symbol::WEIBULL_PDF,
        Law::Cauchy => Symbol::CAUCHY_PDF,
    };
    Expr::call(head, parameters)
}

/// The sign of an evaluated numeric expression, when its approximation
/// decides it; `None` for another.
fn sign(ev: &Eval<'_>, x: &Expr) -> Option<Ordering> {
    if let Some(q) = x.to_rational() {
        return Some(q.sign());
    }
    if !x.is_numeric() {
        return None;
    }
    let wp = ev.context().working_precision();
    match approximation(ev.context(), x, wp).ok()?.view() {
        View::Ball(b) => b.sign(),
        _ => None,
    }
}

/// The sides of an equation (C0016).
fn equation(e: &Expr) -> Option<(&Expr, &Expr)> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => Some((l, r)),
        _ => None,
    }
}

/// The ends of `a .. b` (C0029).
fn range(e: &Expr) -> Option<(&Expr, &Expr)> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::RANGE), [a, b]) => Some((a, b)),
        _ => None,
    }
}

/// `RandomVariable[…]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's forms and options, in its order"
)]
pub(in crate::native) fn random_variable(ev: &mut Eval<'_>, e: &Expr) -> Native {
    // Its own objects, a fixed point.
    if variable(e).is_some() {
        return Ok(None);
    }
    let Some((f, rest)) = e.args().split_first() else {
        return Ok(None);
    };
    let law = name(f);
    if rest.is_empty() {
        if variable(f).is_some() {
            return Ok(Some(f.clone()));
        }
        if is_list(f) {
            return listed(ev, e, f.args());
        }
        let one = || Expr::from(1);
        return Ok(match law {
            Some(Name::Law(Law::Normal)) => Some(object(Law::Normal, vec![Expr::from(0), one()])),
            Some(Name::Law(Law::Uniform)) => Some(object(Law::Uniform, vec![Expr::from(0), one()])),
            Some(Name::Law(Law::Exponential)) => Some(object(Law::Exponential, vec![one()])),
            _ => None,
        });
    }
    let function = f.head() == Some(&Expr::from(Symbol::FUNCTION));
    if law.is_none() && !is_list(f) && !function {
        return Ok(None);
    }
    // giac's positional form (:614-659).
    if let Some(Name::Law(law)) = law
        && rest.iter().all(parameter)
    {
        if rest.len() != arity(law) {
            return refuse(ev, e, BAD);
        }
        return Ok(Some(object(law, rest.to_vec())));
    }
    // The options (:663-703), the last given winning.
    let mut mean = None;
    let mut spread: Option<(bool, &Expr)> = None;
    let mut ends = None;
    let mut items = None;
    let mut categories = None;
    for (i, x) in rest.iter().enumerate() {
        let place = i + 1;
        if let Some((lh, rh)) = equation(x) {
            match lh.as_symbol() {
                Some(Symbol::MEAN) => mean = Some(rh),
                Some(Symbol::POPULATION_STANDARD_DEVIATION) => spread = Some((false, rh)),
                Some(Symbol::POPULATION_VARIANCE) => {
                    if sign(ev, rh) == Some(Ordering::Less) {
                        return refuse(ev, e, "Variance must be nonnegative");
                    }
                    spread = Some((true, rh));
                }
                Some(s) if s.name() == "range" => {
                    if is_list(rh) {
                        let [a, b] = rh.args() else {
                            return refuse(ev, e, BAD);
                        };
                        ends = Some((a, b));
                    } else if let Some(r) = range(rh) {
                        ends = Some(r);
                    } else {
                        return refuse(ev, e, "Invalid range specification");
                    }
                }
                _ => return refuse(ev, e, "Unknown distribution parameter"),
            }
        } else if let Some(r) = range(x) {
            ends = Some(r);
        } else if is_list(x) {
            if is_list(f) {
                if place != 1 {
                    return refuse(ev, e, BAD);
                }
                items = Some(x);
            } else if place == 1 {
                categories = Some(x);
            } else {
                items = Some(x);
            }
        } else {
            return refuse(ev, e, "Invalid distribution specification");
        }
    }
    let sdev = match spread {
        None => None,
        Some((false, s)) => Some(s.clone()),
        Some((true, v)) => Some(ev.eval(&Expr::call(Symbol::SQRT, [v.clone()]))?),
    };
    if let Some(s) = &sdev {
        match sign(ev, s) {
            Some(Ordering::Greater) => {}
            Some(_) => return refuse(ev, e, "Standard deviation must be positive"),
            None => return Ok(None),
        }
    }
    if is_list(f) {
        return discrete(ev, e, f.args(), items.map(Expr::args));
    }
    match law {
        Some(Name::Multinomial) => match categories {
            Some(p) => discrete(ev, e, p.args(), items.map(Expr::args)),
            None => refuse(ev, e, BAD),
        },
        Some(Name::Law(law)) => {
            if mean.is_some_and(|m| !m.is_numeric()) {
                return Ok(None);
            }
            fit(ev, e, law, mean, sdev.as_ref(), ends)
        }
        None => density(ev, e, f, ends, items),
    }
}

/// `randvar(l)`: a matrix of `[v, w]` rows, or the weights.
fn listed(ev: &mut Eval<'_>, e: &Expr, rows: &[Expr]) -> Native {
    let matrix = !rows.is_empty()
        && rows.iter().all(is_list)
        && rows.iter().all(|r| r.args().len() == rows[0].args().len());
    if !matrix {
        return discrete(ev, e, rows, None);
    }
    if rows[0].args().len() != 2 {
        return refuse(ev, e, "Expected a two-column matrix");
    }
    let values: Vec<Expr> = rows.iter().map(|r| r.args()[0].clone()).collect();
    let weights: Vec<Expr> = rows.iter().map(|r| r.args()[1].clone()).collect();
    discrete(ev, e, &weights, Some(&values))
}

/// The discrete law of `weights` over `values`, checked as giac checks it.
fn discrete(ev: &mut Eval<'_>, e: &Expr, ws: &[Expr], values: Option<&[Expr]>) -> Native {
    if ws.is_empty() {
        return refuse(ev, e, BAD);
    }
    if values.is_some_and(|v| v.len() != ws.len()) {
        return refuse(ev, e, "Items and weights sizes do not match");
    }
    if !weights(ws) {
        return refuse(ev, e, BAD);
    }
    let mut args = vec![list(ws.to_vec())];
    if let Some(v) = values {
        args.push(list(v.to_vec()));
    }
    Ok(Some(Expr::call(Symbol::RANDOM_VARIABLE, args)))
}

/// `randvar(f, range=a..b(, V))`: the weights `f(a) … f(b)`, over `a … b`
/// or `V`, between integers `a < b`.
fn density(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    ends: Option<(&Expr, &Expr)>,
    items: Option<&Expr>,
) -> Native {
    let Some((a, b)) = ends.and_then(|(a, b)| Some((a.to_integer()?, b.to_integer()?))) else {
        return refuse(ev, e, BAD);
    };
    if b <= a {
        return refuse(ev, e, BAD);
    }
    let most = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    let width = (&b - &a)
        .to_i64()
        .and_then(|w| u64::try_from(w).ok())
        .and_then(|w| w.checked_add(1));
    let Some(width) = width.filter(|&w| w <= most) else {
        return refuse(ev, e, "too many weights");
    };
    let mut support = Vec::new();
    let mut ws = Vec::new();
    let mut k = a;
    for _ in 0..width {
        let point = Expr::from(k.clone());
        ws.push(ev.eval(&Expr::apply(f, std::slice::from_ref(&point)))?);
        support.push(point);
        k = &k + &Integer::ONE;
    }
    let values = match items {
        Some(v) => v.args(),
        None => &support[..],
    };
    discrete(ev, e, &ws, Some(values))
}

/// A law's object from its moments or its range (:704-823).
#[allow(clippy::too_many_lines, reason = "one arm per law, as giac's")]
fn fit(
    ev: &mut Eval<'_>,
    e: &Expr,
    law: Law,
    mean: Option<&Expr>,
    sdev: Option<&Expr>,
    ends: Option<(&Expr, &Expr)>,
) -> Native {
    use Ordering::{Equal, Greater, Less};
    let n = |k: i64| Expr::from(k);
    // The mean, the deviation and the variance, all three given.
    let moments = |ev: &mut Eval<'_>| -> Result<Option<(Expr, Expr, Expr)>, EvalError> {
        match (mean, sdev) {
            (Some(m), Some(s)) if sign(ev, m) == Some(Greater) => {
                let v = ev.eval(&s.clone().pow(n(2)))?;
                Ok(Some((m.clone(), s.clone(), v)))
            }
            _ => Ok(None),
        }
    };
    let parameters: Vec<Expr> = match law {
        Law::Normal => vec![
            mean.cloned().unwrap_or_else(|| n(0)),
            sdev.cloned().unwrap_or_else(|| n(1)),
        ],
        Law::Uniform => {
            if let Some((a, b)) = ends {
                let width = ev.eval(&(b.clone() - a.clone()))?;
                if !(a.is_numeric() && b.is_numeric()) || sign(ev, &width) != Some(Greater) {
                    return refuse(ev, e, INVALID);
                }
                vec![a.clone(), b.clone()]
            } else if let (Some(m), Some(s)) = (mean, sdev) {
                let r = Expr::call(Symbol::SQRT, [n(3)]) * s.clone();
                vec![
                    ev.eval(&(m.clone() - r.clone()))?,
                    ev.eval(&(m.clone() + r))?,
                ]
            } else {
                return refuse(ev, e, BAD);
            }
        }
        Law::Poisson => {
            let lambda = match (mean, sdev) {
                (Some(m), s) => {
                    if let Some(s) = s {
                        let d = ev.eval(&(s.clone().pow(n(2)) - m.clone()))?;
                        if sign(ev, &d) != Some(Equal) {
                            return refuse(ev, e, BAD);
                        }
                    }
                    m.clone()
                }
                (None, Some(s)) => ev.eval(&s.clone().pow(n(2)))?,
                (None, None) => return refuse(ev, e, BAD),
            };
            if sign(ev, &lambda) != Some(Greater) {
                return refuse(ev, e, BAD);
            }
            vec![lambda]
        }
        Law::Exponential => {
            let lambda = match (mean, sdev) {
                (Some(m), s) => {
                    if let Some(s) = s {
                        let d = ev.eval(&(s.clone() - m.clone()))?;
                        if sign(ev, &d) != Some(Equal) {
                            return refuse(ev, e, BAD);
                        }
                    }
                    if sign(ev, m) != Some(Greater) {
                        return refuse(ev, e, BAD);
                    }
                    ev.eval(&(n(1) / m.clone()))?
                }
                (None, Some(s)) => ev.eval(&(n(1) / s.clone()))?,
                (None, None) => return refuse(ev, e, BAD),
            };
            vec![lambda]
        }
        Law::Weibull => {
            let Some((m, _, v)) = moments(ev)? else {
                return refuse(ev, e, BAD);
            };
            match weibull(ev, &m, &v)? {
                Some((k, lambda))
                    if sign(ev, &k) == Some(Greater) && sign(ev, &lambda) == Some(Greater) =>
                {
                    vec![k, lambda]
                }
                _ => return refuse(ev, e, INVALID),
            }
        }
        Law::Gamma => {
            let Some((m, _, v)) = moments(ev)? else {
                return refuse(ev, e, BAD);
            };
            vec![
                ev.eval(&(m.clone().pow(n(2)) / v.clone()))?,
                ev.eval(&(m / v))?,
            ]
        }
        Law::Beta => {
            let Some((m, _, v)) = moments(ev)? else {
                return refuse(ev, e, BAD);
            };
            let below_one = ev.eval(&(n(1) - m.clone()))?;
            if sign(ev, &below_one) != Some(Greater) {
                return refuse(ev, e, BAD);
            }
            let fac = (v.clone() + m.clone().pow(n(2)) - m.clone()) / v;
            let p1 = ev.eval(&(-m.clone() * fac.clone()))?;
            let p2 = ev.eval(&((m - n(1)) * fac))?;
            if sign(ev, &p1) != Some(Greater) || sign(ev, &p2) != Some(Greater) {
                return refuse(ev, e, INVALID);
            }
            vec![p1, p2]
        }
        Law::Geometric => {
            let p = match (mean, sdev) {
                (Some(m), s) => {
                    if sign(ev, m) != Some(Greater) {
                        return refuse(ev, e, BAD);
                    }
                    if let Some(s) = s {
                        let d = s.clone().pow(n(2)) - m.clone().pow(n(2)) - m.clone();
                        let d = ev.eval(&d)?;
                        if sign(ev, &d) != Some(Equal) {
                            return refuse(ev, e, BAD);
                        }
                    }
                    ev.eval(&(n(1) / m.clone()))?
                }
                (None, Some(s)) => {
                    let v = s.clone().pow(n(2));
                    let root = Expr::call(Symbol::SQRT, [n(4) * v.clone() + n(1)]);
                    ev.eval(&((root - n(1)) / (n(2) * v)))?
                }
                (None, None) => return refuse(ev, e, BAD),
            };
            let above_one = ev.eval(&(p.clone() - n(1)))?;
            if sign(ev, &p) == Some(Less) || sign(ev, &above_one) == Some(Greater) {
                return refuse(ev, e, INVALID);
            }
            vec![p]
        }
        Law::Fisher | Law::ChiSquare | Law::Cauchy | Law::StudentT => {
            return refuse(
                ev,
                e,
                "Specifying moments is not supported for this distribution",
            );
        }
        Law::Binomial => {
            let Some((m, _, v)) = moments(ev)? else {
                return refuse(ev, e, BAD);
            };
            let tmp = ev.eval(&(m.clone() - v))?;
            if sign(ev, &tmp) == Some(Equal) {
                return refuse(ev, e, BAD);
            }
            let size = ev.eval(&(m.clone().pow(n(2)) / tmp.clone()))?;
            let (size, p) = if size.to_integer().is_some() {
                (size, ev.eval(&(tmp / m))?)
            } else {
                let size = ev.eval(&Expr::call(Symbol::ROUND, [size]))?;
                if sign(ev, &size) == Some(Equal) {
                    return refuse(ev, e, INVALID);
                }
                let p = ev.eval(&(m / size.clone()))?;
                (size, p)
            };
            let above_one = ev.eval(&(p.clone() - n(1)))?;
            if size.to_integer().is_none()
                || sign(ev, &p) == Some(Less)
                || sign(ev, &above_one) == Some(Greater)
            {
                return refuse(ev, e, INVALID);
            }
            vec![size, p]
        }
        Law::NegativeBinomial => {
            let Some((m, _, v)) = moments(ev)? else {
                return refuse(ev, e, BAD);
            };
            let d = ev.eval(&(v.clone() - m.clone()))?;
            if sign(ev, &d) == Some(Equal) {
                return refuse(ev, e, BAD);
            }
            let size = ev.eval(&Expr::call(Symbol::ROUND, [m.clone().pow(n(2)) / d]))?;
            let p = ev.eval(&(m / v))?;
            let above_one = ev.eval(&(p.clone() - n(1)))?;
            let below_one = ev.eval(&(size.clone() - n(1)))?;
            if size.to_integer().is_none()
                || sign(ev, &below_one) == Some(Less)
                || sign(ev, &p) == Some(Less)
                || sign(ev, &above_one) == Some(Greater)
            {
                return refuse(ev, e, INVALID);
            }
            vec![size, p]
        }
    };
    Ok(Some(object(law, parameters)))
}

/// Weibull's `(k, λ)` from the mean `m` and the variance `v` as giac
/// finds them (:752-764): `Γ(1 + 2/k)/Γ(1 + 1/k) − v/m² − 1 = 0` by
/// Newton's method from `max(1, 1/v)`, certified by `FSolve`, and `λ =
/// m/Γ(1 + 1/k)`; `None` when no root is found.
fn weibull(ev: &mut Eval<'_>, m: &Expr, v: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    if sign(ev, v) != Some(Ordering::Greater) {
        return Ok(None);
    }
    let k = Expr::from(Symbol::new("RandomVariable$k"));
    let n = |k: i64| Expr::from(k);
    let gamma = |x: Expr| Expr::call(Symbol::GAMMA, [x]);
    let f = gamma(n(1) + n(2) / k.clone()) / gamma(n(1) + n(1) / k.clone())
        - v.clone() / m.clone().pow(n(2))
        - n(1);
    let inverse = ev.eval(&(n(1) / v.clone()))?;
    let over = ev.eval(&(inverse.clone() - n(1)))?;
    let start = if sign(ev, &over) == Some(Ordering::Greater) {
        inverse
    } else {
        n(1)
    };
    let solve = Expr::call(
        Symbol::FSOLVE,
        [
            Expr::call(Symbol::EQUATION, [f, n(0)]),
            k,
            start,
            Expr::from(Symbol::NEWTON_SOLVER),
        ],
    );
    let root = ev.eval(&solve)?;
    if !root.is_number() {
        return Ok(None);
    }
    let lambda = ev.eval(&(m.clone() / gamma(n(1) + n(1) / root.clone())))?;
    Ok(Some((root, lambda)))
}

#[cfg(test)]
mod tests;
