//! Definite integrals (CAS-DOC-9 §5, CAS-198): `F(b) − F(a)` for the
//! antiderivative `F` of CAS-194 to CAS-197 when the integrand is
//! continuous between the bounds, for then so is `F`: its singular points
//! are the integrand's, the poles of a rational function at the poles and
//! logarithms of `F`, `ln|cos u|` where `tan u` has its poles, `u ln u − u`
//! where `ln u` is not defined, and the arctangents of polynomials of
//! CAS-197 have no jumps.
//!
//! The interval holds the bounds: a rational as it is, a real number as
//! `Pi`, `√2` or `E` between the ends of its ball at 128 bits, another real
//! bound, as a symbol in the Real domain, between the bounds its
//! assumptions give, if any. A bound that need not be real leaves only the
//! integrands continuous on the whole plane, whatever the path.
//!
//! What is continuous on the interval, by its structure: what is free of
//! `x`; a rational function of `x` whose denominator in lowest terms has no
//! root there (Sturm's count, cas-poly); sums and products of what is;
//! `g^n` for an integer `n ≥ 0`, for an integer `n < 0` where the
//! polynomial `g` has no root, for any other rational `n` where `g` is
//! positive, or for `n > 0` not negative; `c^g` for `c > 0`; `Exp`, `Sin`,
//! `Cos`, and `ArcTan` of a real argument, of what is; `Log[g]` where the
//! polynomial `g` is positive, `Log[Abs[g]]` where it has no root;
//! `Tan[u]` for a linear `u` where no `(k + 1/2) π` lies between its values
//! at the ends.
//!
//! An improper integral (CAS-DOC-10 §7, CAS-207), with an infinite bound
//! or a singular point at a bound, is the difference of the one-sided
//! limits of `F` at the bounds, from inside, where the integrand is
//! continuous on the open interval: `∫₀^∞ e^{−x} = 1`, `∫₀¹ 1/√x = 2`,
//! `∫₀¹ ln x = −1`. A rational integrand with poles inside at rational
//! points splits there, each piece by its limits; pieces that diverge make
//! an infinity of their sign (`∫_{−1}^{1} 1/x² = Infinity`), and of both
//! signs the call stays with a warning, where giac answers the principal
//! value, `∫_{−1}^{1} 1/x = 0`. What is not shown continuous stays, with a
//! warning.

use cas_num::algebra::{EuclideanDomain, Rationals};
use cas_num::{Integer, Rational};
use cas_poly::{Dense, DenseRing};

use super::{Native, polynomial, rational_function, slope};
use crate::assume::{Sign, range, sign_of};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;
use crate::trace::Warning;

/// The bits of the balls that bound real numbers.
const BITS: u64 = 128;

/// A closed interval with rational ends, `None` at an infinite one.
#[derive(Clone, Debug)]
struct Interval {
    lo: Option<Rational>,
    hi: Option<Rational>,
}

impl Interval {
    fn point(q: Rational) -> Interval {
        Interval {
            lo: Some(q.clone()),
            hi: Some(q),
        }
    }

    /// The roots of `p` in this interval, at its finite ends but when
    /// `open`.
    fn roots(&self, p: &Dense<Rational>, open: bool) -> usize {
        let qx = DenseRing::new(Rationals);
        let all = qx.count_real_roots(p, self.lo.as_ref(), self.hi.as_ref());
        if !open {
            return all;
        }
        let mut ends = 0;
        if let Some(lo) = &self.lo {
            ends += usize::from(qx.eval(p, lo).is_zero());
        }
        if let Some(hi) = &self.hi
            && self.lo.as_ref() != Some(hi)
        {
            ends += usize::from(qx.eval(p, hi).is_zero());
        }
        all.saturating_sub(ends)
    }

    /// A point inside, away from its ends when it has two.
    fn inside(&self) -> Rational {
        let one = Rational::ONE;
        match (&self.lo, &self.hi) {
            (Some(lo), Some(hi)) => &(lo + hi) / &Rational::from(Integer::from(2)),
            (Some(lo), None) => lo + &one,
            (None, Some(hi)) => hi - &one,
            (None, None) => Rational::ZERO,
        }
    }
}

/// Where the integral runs.
enum Region {
    /// Between real bounds: `outer` holds them, and `inner`, when known,
    /// lies between them; `open` leaves the ends out.
    Real {
        outer: Interval,
        inner: Option<Interval>,
        open: bool,
    },
    /// A bound need not be real.
    Plane,
}

/// What is known of the integrand on the region.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Continuity {
    Continuous,
    /// Surely unbounded at a point, or undefined there.
    Singular,
    Unknown,
}

impl Continuity {
    /// Both.
    fn and(self, other: Continuity) -> Continuity {
        use Continuity::{Continuous, Singular, Unknown};
        match (self, other) {
            (Singular, _) | (_, Singular) => Singular,
            (Continuous, Continuous) => Continuous,
            _ => Unknown,
        }
    }
}

/// `F(b) − F(a)` for `e = Integrate[f, x, a, b]` and the antiderivative
/// `F`, or the improper integral by limits; `None` after a warning when
/// neither is shown.
pub(super) fn integral(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
    antiderivative: &Expr,
    [a, b]: [&Expr; 2],
) -> Native {
    if antiderivative.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(None);
    }
    // Evaluated first, so that no `x^0` meets `x = 0`.
    let antiderivative = ev.eval(antiderivative)?;
    let bounded = !infinite(a) && !infinite(b);
    if bounded && continuity(ev, e, f, x, &region(ev, a, b, false))? == Continuity::Continuous {
        let at = |v: &Expr| antiderivative.substitute(&[(x.clone(), v.clone())]);
        return Ok(Some(at(b) - at(a)));
    }
    // Improper: continuous on the open interval, by limits at the bounds.
    if [a, b]
        .into_iter()
        .any(|bound| infinite(bound) && end(bound).is_none())
    {
        warn(ev, "an infinite bound without a sign", e);
        return Ok(None);
    }
    let open = region(ev, a, b, true);
    let Some(ascending) = ascending(ev, a, b) else {
        warn(
            ev,
            "the integrand is not shown continuous between the bounds",
            e,
        );
        return Ok(None);
    };
    let (lo, hi) = if ascending { (a, b) } else { (b, a) };
    let points = match continuity(ev, e, f, x, &open)? {
        Continuity::Continuous => vec![lo.clone(), hi.clone()],
        Continuity::Singular => {
            let Some(poles) = interior_poles(ev, e, f, x, &open)? else {
                warn(
                    ev,
                    "an improper integral: the integrand is unbounded between the bounds",
                    e,
                );
                return Ok(None);
            };
            let mut points = vec![lo.clone()];
            points.extend(poles.into_iter().map(Expr::from));
            points.push(hi.clone());
            points
        }
        Continuity::Unknown => {
            warn(
                ev,
                "the integrand is not shown continuous between the bounds",
                e,
            );
            return Ok(None);
        }
    };
    let Some(value) = by_limits(ev, x, &antiderivative, &points)? else {
        warn(
            ev,
            "an improper integral whose parts diverge, or whose limits are not found",
            e,
        );
        return Ok(None);
    };
    Ok(Some(if ascending { value } else { -value }))
}

/// `∫ₐᵇ g(x) DiracDelta[α x + β] dx` for one delta at a linear argument
/// (CAS-DOC-18 §2, CAS-273): `g(x₀)/|α|` at `x₀ = −β/α` strictly between the
/// bounds, 0 outside them, the sign of reversed bounds; `None` at an end,
/// when undecided, or for another integrand.
pub(super) fn sifting(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    [a, b]: [&Expr; 2],
) -> Result<Option<Expr>, EvalError> {
    let delta = Expr::from(Symbol::DIRAC_DELTA);
    let factors = if f.head() == Some(&Expr::from(Symbol::TIMES)) {
        f.args().to_vec()
    } else {
        vec![f.clone()]
    };
    let (deltas, rest): (Vec<Expr>, Vec<Expr>) =
        factors.into_iter().partition(|g| g.head() == Some(&delta));
    let ([d], false) = (deltas.as_slice(), rest.iter().any(|g| g.contains(&delta))) else {
        return Ok(None);
    };
    let [u] = d.args() else {
        return Ok(None);
    };
    let alpha = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
    if alpha.contains(x) || alpha == Expr::from(0) {
        return Ok(None);
    }
    let beta = ev.eval(&u.substitute(&[(x.clone(), Expr::from(0))]))?;
    let x0 = ev.eval(&(-(beta / alpha.clone())))?;
    let Some(up) = ascending(ev, a, b) else {
        return Ok(None);
    };
    let (lo, hi) = if up { (a, b) } else { (b, a) };
    let value = match (ascending(ev, lo, &x0), ascending(ev, &x0, hi)) {
        (Some(true), Some(true)) => {
            let g = Expr::call(Symbol::TIMES, rest).substitute(&[(x.clone(), x0)]);
            g / Expr::call(Symbol::ABS, [alpha])
        }
        (Some(false), _) | (_, Some(false)) => Expr::from(0),
        _ => return Ok(None),
    };
    let signed = if up { value } else { -value };
    ev.eval(&signed).map(Some)
}

/// Whether `a < b` surely, `a > b` surely, or `None`.
fn ascending(ev: &Eval<'_>, a: &Expr, b: &Expr) -> Option<bool> {
    match (end(a), end(b)) {
        (Some(true), _) | (_, Some(false)) => return Some(false),
        (Some(false), _) | (_, Some(true)) => return Some(true),
        _ => {}
    }
    let (p, q) = (enclosure(ev, a)?, enclosure(ev, b)?);
    let below =
        |p: &Interval, q: &Interval| matches!((&p.hi, &q.lo), (Some(hi), Some(lo)) if hi < lo);
    if below(&p, &q) {
        Some(true)
    } else if below(&q, &p) {
        Some(false)
    } else {
        None
    }
}

/// `Some(true)` for `+∞`, `Some(false)` for `−∞`, `None` for a finite
/// bound.
pub(super) fn end(e: &Expr) -> Option<bool> {
    if !infinite(e) || e.contains(&Expr::from(Symbol::COMPLEX_INFINITY)) {
        return None;
    }
    let negative = e.contains(&Expr::from(-1)) || e.to_rational().is_some_and(|q| q.is_negative());
    Some(!negative)
}

/// What the parts of an integral divide into.
enum Part {
    Finite(Expr),
    Infinite(bool),
}

/// The limit of `F` at `point` from `side` (1 from the right, −1 from the
/// left), by `Limit`.
fn one_sided(
    ev: &mut Eval<'_>,
    x: &Expr,
    antiderivative: &Expr,
    point: &Expr,
    side: i64,
) -> Result<Option<Part>, EvalError> {
    let limit = Expr::call(
        Symbol::LIMIT,
        [
            antiderivative.clone(),
            x.clone(),
            point.clone(),
            Expr::from(side),
        ],
    );
    let value = ev.eval(&limit)?;
    if value.contains(&Expr::from(Symbol::LIMIT))
        || value.contains(&Expr::from(Symbol::COMPLEX_INFINITY))
    {
        return Ok(None);
    }
    Ok(Some(match end(&value) {
        Some(positive) => Part::Infinite(positive),
        None if infinite(&value) || value.contains(&Expr::from(Symbol::INDETERMINATE)) => {
            return Ok(None);
        }
        None => Part::Finite(value),
    }))
}

/// `Σ F(p_{k+1}⁻) − F(p_k⁺)` over the ascending points; an infinity when
/// the divergent pieces agree in sign; `None` when they do not, or a limit
/// is not found.
fn by_limits(
    ev: &mut Eval<'_>,
    x: &Expr,
    antiderivative: &Expr,
    points: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let mut finite = Vec::new();
    let (mut up, mut down) = (false, false);
    for pair in points.windows(2) {
        let (Some(start), Some(stop)) = (
            one_sided(ev, x, antiderivative, &pair[0], 1)?,
            one_sided(ev, x, antiderivative, &pair[1], -1)?,
        ) else {
            return Ok(None);
        };
        match (stop, start) {
            (Part::Finite(s), Part::Finite(r)) => finite.push(s - r),
            (Part::Infinite(true), Part::Finite(_) | Part::Infinite(false))
            | (Part::Finite(_), Part::Infinite(false)) => up = true,
            (Part::Infinite(false), Part::Finite(_) | Part::Infinite(true))
            | (Part::Finite(_), Part::Infinite(true)) => down = true,
            (Part::Infinite(_), Part::Infinite(_)) => return Ok(None),
        }
    }
    Ok(match (up, down) {
        (true, true) => None,
        (true, false) => Some(Expr::from(Symbol::INFINITY)),
        (false, true) => Some(-Expr::from(Symbol::INFINITY)),
        (false, false) => Some(finite.into_iter().sum()),
    })
}

/// The poles of a rational integrand strictly inside the open region, when
/// they are rational, ascending; `None` for another integrand, or a pole
/// that is not rational.
fn interior_poles(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
    region: &Region,
) -> Result<Option<Vec<Rational>>, EvalError> {
    let Region::Real {
        outer,
        inner: Some(inner),
        ..
    } = region
    else {
        return Ok(None);
    };
    let Some((num, den)) = rational_function(ev, e, f, x)? else {
        return Ok(None);
    };
    let qx = DenseRing::new(Rationals);
    let g = qx.gcd(&num, &den);
    let den = qx.div_rem(&den, &g).0;
    let factors = qx.factor(&den);
    let mut poles = Vec::new();
    for (factor, _) in factors.factors() {
        if outer.roots(factor, true) == 0 {
            continue;
        }
        // A root inside: rational, and surely inside, or not split.
        if factor.degree() != Some(1) {
            return Ok(None);
        }
        let root = -(&factor.coeffs()[0] / &factor.coeffs()[1]);
        let strictly = |i: &Interval| {
            i.lo.as_ref().is_none_or(|lo| lo < &root) && i.hi.as_ref().is_none_or(|hi| &root < hi)
        };
        if !strictly(inner) {
            return Ok(None);
        }
        poles.push(root);
    }
    poles.sort();
    Ok(Some(poles))
}

fn warn(ev: &mut Eval<'_>, message: &'static str, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// Whether `e` has an infinity in it, as `-Infinity`.
pub(super) fn infinite(e: &Expr) -> bool {
    [
        Symbol::INFINITY,
        Symbol::DIRECTED_INFINITY,
        Symbol::COMPLEX_INFINITY,
    ]
    .into_iter()
    .any(|s| e.contains(&Expr::from(s)))
}

/// The ends of the ball of a real number `e`, at [`BITS`] bits.
fn ball_ends(ev: &Eval<'_>, e: &Expr) -> Option<(Rational, Rational)> {
    let value = approximation(ev.context(), e, BITS).ok()?;
    let View::Ball(b) = value.view() else {
        return None;
    };
    let mid = b.mid().to_rational_within(4 * BITS)?;
    let rad = b.rad().to_float()?.to_rational_within(4 * BITS)?;
    Some((&mid - &rad, &mid + &rad))
}

/// An interval that holds the value of `e`: a point for a rational, the
/// ends of its ball for another real number, its assumed bounds for
/// another real value; `None` when it need not be real.
fn enclosure(ev: &Eval<'_>, e: &Expr) -> Option<Interval> {
    if let Some(q) = e.to_rational() {
        return Some(Interval::point(q));
    }
    if let Some((lo, hi)) = ball_ends(ev, e) {
        return Some(Interval {
            lo: Some(lo),
            hi: Some(hi),
        });
    }
    let r = range(ev.context(), e);
    let (lo, hi) = r.interval()?;
    Some(Interval {
        lo: lo.cloned(),
        hi: hi.cloned(),
    })
}

/// The region between `a` and `b`, without its ends when `open`: an
/// infinite bound leaves its side unbounded.
fn region(ev: &Eval<'_>, a: &Expr, b: &Expr, open: bool) -> Region {
    let bound = |e: &Expr| -> Option<(Interval, Option<bool>)> {
        match end(e) {
            Some(positive) => Some((Interval { lo: None, hi: None }, Some(positive))),
            None => enclosure(ev, e).map(|i| (i, None)),
        }
    };
    let (Some((p, pe)), Some((q, qe))) = (bound(a), bound(b)) else {
        return Region::Plane;
    };
    // Lower and upper: an infinite end is on its side.
    let (lower, upper) = match (pe, qe) {
        (Some(true), _) | (_, Some(false)) => ((q, qe), (p, pe)),
        _ => ((p, pe), (q, qe)),
    };
    let (lower, lower_end) = lower;
    let (upper, upper_end) = upper;
    let (outer, inner) = match (lower_end, upper_end) {
        (Some(false), Some(true)) => (
            Interval { lo: None, hi: None },
            Some(Interval { lo: None, hi: None }),
        ),
        (Some(false), _) => (
            Interval {
                lo: None,
                hi: upper.hi.clone(),
            },
            Some(Interval {
                lo: None,
                hi: upper.lo.clone(),
            }),
        ),
        (_, Some(true)) => (
            Interval {
                lo: lower.lo.clone(),
                hi: None,
            },
            Some(Interval {
                lo: lower.hi.clone(),
                hi: None,
            }),
        ),
        _ => {
            let min = |p: &Option<Rational>, q: &Option<Rational>| {
                p.as_ref().zip(q.as_ref()).map(|(p, q)| p.min(q).clone())
            };
            let max = |p: &Option<Rational>, q: &Option<Rational>| {
                p.as_ref().zip(q.as_ref()).map(|(p, q)| p.max(q).clone())
            };
            let outer = Interval {
                lo: min(&lower.lo, &upper.lo),
                hi: max(&lower.hi, &upper.hi),
            };
            // What lies between them when one is surely below the other.
            let between = |p: &Interval, q: &Interval| match (&p.hi, &q.lo) {
                (Some(hi), Some(lo)) if hi <= lo => Some(Interval {
                    lo: Some(hi.clone()),
                    hi: Some(lo.clone()),
                }),
                _ => None,
            };
            let inner = between(&lower, &upper).or_else(|| between(&upper, &lower));
            (outer, inner)
        }
    };
    Region::Real { outer, inner, open }
}

/// What the real roots of `p` tell: none in the region, one surely in it,
/// or neither shown.
fn roots(region: &Region, p: &Dense<Rational>) -> Continuity {
    if p.degree().unwrap_or(0) == 0 {
        return Continuity::Continuous;
    }
    let Region::Real { outer, inner, open } = region else {
        return Continuity::Unknown;
    };
    if outer.roots(p, *open) == 0 {
        Continuity::Continuous
    } else if inner.as_ref().is_some_and(|i| i.roots(p, *open) > 0) {
        Continuity::Singular
    } else {
        Continuity::Unknown
    }
}

/// Whether `p` is positive on the region, or, unless `strict`, positive
/// inside and maybe 0 at the ends.
fn positive(region: &Region, p: &Dense<Rational>, strict: bool) -> Continuity {
    let Region::Real { outer, open, .. } = region else {
        return Continuity::Unknown;
    };
    // On an open interval, positive inside is enough.
    let strict = strict && !open;
    let qx = DenseRing::new(Rationals);
    let at_ends = [&outer.lo, &outer.hi]
        .into_iter()
        .flatten()
        .filter(|end| qx.eval(p, end).is_zero())
        .count();
    let found = roots(region, p);
    let shown = match found {
        Continuity::Continuous => true,
        // Roots at the ends only.
        _ if !strict => outer.roots(p, false) == at_ends,
        _ => return found,
    };
    if shown && qx.eval(p, &outer.inside()).is_positive() {
        Continuity::Continuous
    } else {
        Continuity::Unknown
    }
}

/// What is known of `f` on the region (the module's rules).
fn continuity(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
    region: &Region,
) -> Result<Continuity, EvalError> {
    if !f.contains(x) {
        return Ok(Continuity::Continuous);
    }
    if f == x {
        return Ok(Continuity::Continuous);
    }
    let Some(head) = f.head().and_then(Expr::as_symbol) else {
        return Ok(Continuity::Unknown);
    };
    // A rational function of `x` alone, its denominator in lowest terms.
    if matches!(head, Symbol::PLUS | Symbol::TIMES | Symbol::POWER)
        && let Some((num, den)) = rational_function(ev, e, f, x)?
    {
        let qx = DenseRing::new(Rationals);
        let g = qx.gcd(&num, &den);
        return Ok(roots(region, &qx.div_rem(&den, &g).0));
    }
    let args = f.args();
    let all = |ev: &mut Eval<'_>, parts: &[Expr]| -> Result<Continuity, EvalError> {
        let mut known = Continuity::Continuous;
        for part in parts {
            known = known.and(continuity(ev, e, part, x, region)?);
            if known == Continuity::Singular {
                break;
            }
        }
        Ok(known)
    };
    Ok(match (head, args) {
        (Symbol::PLUS | Symbol::TIMES, _) | (Symbol::EXP | Symbol::SIN | Symbol::COS, [_]) => {
            all(ev, args)?
        }
        (Symbol::ARC_TAN, [g]) if matches!(region, Region::Real { .. }) && real(g) => {
            all(ev, args)?
        }
        (Symbol::POWER, [c, g]) if !c.contains(x) => {
            if sign_of(ev.context(), c) == Sign::Positive {
                all(ev, std::slice::from_ref(g))?
            } else {
                Continuity::Unknown
            }
        }
        (Symbol::POWER, [g, n]) => {
            let Some(n) = n.to_rational() else {
                return Ok(Continuity::Unknown);
            };
            let inner = all(ev, std::slice::from_ref(g))?;
            if n.is_integer() && !n.is_negative() {
                return Ok(inner);
            }
            let Some(p) = polynomial(ev, e, g, x)? else {
                return Ok(Continuity::Unknown);
            };
            inner.and(if n.is_integer() {
                roots(region, &p)
            } else {
                positive(region, &p, n.is_negative())
            })
        }
        (Symbol::LOG, [g]) => {
            if let [h] = g.args()
                && g.head() == Some(&Expr::from(Symbol::ABS))
            {
                let Some(p) = polynomial(ev, e, h, x)? else {
                    return Ok(Continuity::Unknown);
                };
                return Ok(all(ev, std::slice::from_ref(h))?.and(roots(region, &p)));
            }
            let Some(p) = polynomial(ev, e, g, x)? else {
                return Ok(Continuity::Unknown);
            };
            all(ev, args)?.and(positive(region, &p, true))
        }
        (Symbol::TAN, [u]) => tangent(ev, u, x, region),
        _ => Continuity::Unknown,
    })
}

/// Whether `g` has no complex number in it, so that it is real where `x`
/// is.
fn real(g: &Expr) -> bool {
    !g.distinct().any(|part| {
        matches!(part.view(), View::Complex(..) | View::ComplexBall(_))
            || part.as_symbol() == Some(Symbol::I)
    })
}

/// `Tan[u]` for a linear `u = s x + c` with real `s` and `c`: continuous
/// where no `(k + 1/2) π` lies between the values of `u` at the ends of the
/// interval, as `t = u/π − 1/2` passes no integer.
fn tangent(ev: &mut Eval<'_>, u: &Expr, x: &Expr, region: &Region) -> Continuity {
    let Region::Real { outer, inner, .. } = region else {
        return Continuity::Unknown;
    };
    if slope(u, x).is_none() || !real(u) {
        return Continuity::Unknown;
    }
    let Some((pi_lo, pi_hi)) = ball_ends(ev, &Expr::from(Symbol::PI)) else {
        return Continuity::Unknown;
    };
    // The values of `t` over an interval, as a rational interval holding
    // them all (`wide`) or held by them (`!wide`).
    let mut values = |interval: &Interval, wide: bool| -> Option<(Rational, Rational)> {
        let ends = [interval.lo.as_ref()?, interval.hi.as_ref()?];
        let mut ts = Vec::new();
        for end in ends {
            let at = u.substitute(&[(x.clone(), Expr::from(end.clone()))]);
            let at = ev.eval(&at).ok()?;
            let (lo, hi) = match at.to_rational() {
                Some(q) => (q.clone(), q),
                None => ball_ends(ev, &at)?,
            };
            // `v/π` over `v` in [lo, hi] and π in [pi_lo, pi_hi].
            let quotients = [&lo / &pi_lo, &lo / &pi_hi, &hi / &pi_lo, &hi / &pi_hi];
            let least = quotients.iter().min().expect("four").clone();
            let most = quotients.iter().max().expect("four").clone();
            ts.push((least, most));
        }
        let half = Rational::new(Integer::ONE, Integer::from(2));
        let (a, b) = (&ts[0], &ts[1]);
        let (lo, hi) = if wide {
            (a.0.clone().min(b.0.clone()), a.1.clone().max(b.1.clone()))
        } else if a.1 <= b.0 {
            (a.1.clone(), b.0.clone())
        } else if b.1 <= a.0 {
            (b.1.clone(), a.0.clone())
        } else {
            return None;
        };
        Some((&lo - &half, &hi - &half))
    };
    let integer_in = |(lo, hi): &(Rational, Rational)| lo.ceil() <= hi.floor();
    match values(outer, true) {
        Some(wide) if !integer_in(&wide) => Continuity::Continuous,
        Some(_) => match inner.as_ref().and_then(|i| values(i, false)) {
            Some(narrow) if integer_in(&narrow) => Continuity::Singular,
            _ => Continuity::Unknown,
        },
        None => Continuity::Unknown,
    }
}
