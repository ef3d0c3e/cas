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
//! An integrand surely unbounded at a point between the bounds makes the
//! integral improper, and one not shown continuous may be: the call stays,
//! with a warning, where giac answers `+infinity`, `undef` or the principal
//! value, `∫_{−1}^{1} 1/x = 0`. Infinite bounds wait for limits.

use cas_num::algebra::{EuclideanDomain, Integers, Rationals};
use cas_num::{Integer, Rational};
use cas_poly::{Dense, DenseRing, Ratio, Sparse, SparseRing};

use super::{Native, in_x, slope};
use crate::assume::{Sign, range, sign_of};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::poly::Fraction;
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

    /// The roots of `p` in this interval.
    fn roots(&self, p: &Dense<Rational>) -> usize {
        DenseRing::new(Rationals).count_real_roots(p, self.lo.as_ref(), self.hi.as_ref())
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
    /// lies between them.
    Real {
        outer: Interval,
        inner: Option<Interval>,
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
/// `F`, or `None` after a warning when it is not shown continuous between
/// the bounds.
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
    if infinite(a) || infinite(b) {
        warn(ev, "an infinite bound, which waits for limits", e);
        return Ok(None);
    }
    let region = region(ev, a, b);
    match continuity(ev, e, f, x, &region)? {
        Continuity::Continuous => {
            // Evaluated first, so that no `x^0` meets `x = 0`.
            let antiderivative = ev.eval(antiderivative)?;
            let at = |v: &Expr| antiderivative.substitute(&[(x.clone(), v.clone())]);
            Ok(Some(at(b) - at(a)))
        }
        Continuity::Singular => {
            warn(
                ev,
                "an improper integral: the integrand is unbounded between the bounds",
                e,
            );
            Ok(None)
        }
        Continuity::Unknown => {
            warn(
                ev,
                "the integrand is not shown continuous between the bounds",
                e,
            );
            Ok(None)
        }
    }
}

fn warn(ev: &mut Eval<'_>, message: &'static str, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// Whether `e` has an infinity in it, as `-Infinity`.
fn infinite(e: &Expr) -> bool {
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

/// The region between `a` and `b`.
fn region(ev: &Eval<'_>, a: &Expr, b: &Expr) -> Region {
    let (Some(a), Some(b)) = (enclosure(ev, a), enclosure(ev, b)) else {
        return Region::Plane;
    };
    let lower = |p: &Option<Rational>, q: &Option<Rational>| {
        p.as_ref().zip(q.as_ref()).map(|(p, q)| p.min(q).clone())
    };
    let upper = |p: &Option<Rational>, q: &Option<Rational>| {
        p.as_ref().zip(q.as_ref()).map(|(p, q)| p.max(q).clone())
    };
    let outer = Interval {
        lo: lower(&a.lo, &b.lo),
        hi: upper(&a.hi, &b.hi),
    };
    // What lies between them when one is surely below the other.
    let between = |p: &Interval, q: &Interval| match (&p.hi, &q.lo) {
        (Some(hi), Some(lo)) if hi <= lo => Some(Interval {
            lo: Some(hi.clone()),
            hi: Some(lo.clone()),
        }),
        _ => None,
    };
    let inner = between(&a, &b).or_else(|| between(&b, &a));
    Region::Real { outer, inner }
}

/// What the real roots of `p` tell: none in the region, one surely in it,
/// or neither shown.
fn roots(region: &Region, p: &Dense<Rational>) -> Continuity {
    if p.degree().unwrap_or(0) == 0 {
        return Continuity::Continuous;
    }
    let Region::Real { outer, inner } = region else {
        return Continuity::Unknown;
    };
    if outer.roots(p) == 0 {
        Continuity::Continuous
    } else if inner.as_ref().is_some_and(|i| i.roots(p) > 0) {
        Continuity::Singular
    } else {
        Continuity::Unknown
    }
}

/// Whether `p` is positive on the region, or, unless `strict`, positive
/// inside and maybe 0 at the ends.
fn positive(region: &Region, p: &Dense<Rational>, strict: bool) -> Continuity {
    let Region::Real { outer, .. } = region else {
        return Continuity::Unknown;
    };
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
        _ if !strict => outer.roots(p) == at_ends,
        _ => return found,
    };
    if shown && qx.eval(p, &outer.inside()).is_positive() {
        Continuity::Continuous
    } else {
        Continuity::Unknown
    }
}

/// `g` as a polynomial in `x` over ℚ, `None` when it is not one.
fn polynomial(
    ev: &mut Eval<'_>,
    e: &Expr,
    g: &Expr,
    x: &Expr,
) -> Result<Option<Dense<Rational>>, EvalError> {
    Ok(rational_function(ev, e, g, x)?.and_then(|(num, den)| {
        (den.degree() == Some(0)).then(|| {
            let qx = DenseRing::new(Rationals);
            let inverse = den.coeffs()[0].checked_recip().expect("not 0");
            qx.scale(&num, &inverse)
        })
    }))
}

/// `g` as a rational function of `x` alone over ℚ, `None` when it is not
/// one.
fn rational_function(
    ev: &mut Eval<'_>,
    e: &Expr,
    g: &Expr,
    x: &Expr,
) -> Result<Option<Ratio>, EvalError> {
    let Some((p, _)) = in_x(ev, e, g, x)? else {
        return Ok(None);
    };
    if p.generators().len() != 1 {
        return Ok(None);
    }
    let Fraction::Integers { num, den } = &p.fractions()[0] else {
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, 1);
    let qx = DenseRing::new(Rationals);
    let to_q = |s: &Sparse<Integer>| {
        let dense = ring.to_dense(s, 0).expect("one variable");
        qx.poly(
            dense
                .coeffs()
                .iter()
                .map(|c| Rational::from(c.clone()))
                .collect(),
        )
    };
    Ok(Some((to_q(num), to_q(den))))
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
    let Region::Real { outer, inner } = region else {
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
