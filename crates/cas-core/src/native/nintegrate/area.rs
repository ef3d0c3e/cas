//! `Area` (CAS-588, CAS-589): Xcas's `area(f, x = a..b [, n, method])`,
//! the integral of `f` from `a` to `b` by a named rule on `n` subintervals.
//!
//! - Without `n` and a method it is `Integrate[f, {x, a, b}]`, left when
//!   that stays.
//! - Each rule is a list of exact rational nodes and weights `(tᵢ, wᵢ)` on
//!   `[0, 1]` and its value `(b − a) Σ wᵢ f(a + tᵢ(b − a))`, giac's
//!   `approx_area` (plot.cc:6216): `left_rectangle`, `right_rectangle`
//!   and `middle_point` take the left ends, right ends or middles of `n`
//!   subintervals, weights `1/n`; `trapezoid` the `n + 1` ends, weights
//!   `1/n`, `1/(2n)` at `a` and `b`; `simpson` `n` panels of two
//!   subintervals, `1/6, 2/3, 1/3, 2/3, …, 2/3, 1/6` times `1/n`.
//! - `rombergt` and `rombergm` (giac's `rombergt` intg.cc:4848 and
//!   `rombergo` :4784): the trapezoid or midpoint sums `R(j, 0)` on
//!   `1, 2, …, 2ⁿ` subintervals, then Richardson's extrapolation
//!   `R(j, k) = (4ᵏ R(j, k − 1) − R(j − 1, k − 1))/(4ᵏ − 1)` on the sums,
//!   the answer `R(n, n)`, exact on polynomials of degree `2n + 1`.
//! - `gauss15` is `NIntegrate[f, {x, a, b}]`, for any `n`.
//!
//! Certification: the integrand compiles to ball operations (`analytic`)
//! and the bounds to balls; the nodes are balls computed from `a` and `b`,
//! so that the result holds the exact rule sum, not the integral. A
//! `gauss15` answer is `NIntegrate`'s certified integral.
//!
//! Departures from giac: the nodes are exact (giac steps by floats); with
//! a symbol in `f`, `a` or `b` the sum is built and evaluated exactly
//! (`a x` from 0 to 1 by the trapezoid rule is `a/2`); a negative `n`
//! leaves the call; `gauss15` never gives giac's pair of a value and an
//! error estimate. `n = 0` is `Indeterminate` for the trapezoid and
//! Simpson rules, 0 for the rectangles, one level for Romberg; a rational
//! or decimal `n` is floored. The nodes are at most [`MAX_EVALUATIONS`]
//! and the context's limit on terms, checked before any evaluation.

use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Integer, Rational};

use super::analytic::{Mode, Node};
use super::{GUARD, MAX_EVALUATIONS, Native, bound, q};
use crate::build::list;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

/// A rule of `Area`, by its option symbol.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rule {
    Left,
    Right,
    Middle,
    Trapezoid,
    Simpson,
    RombergTrapezoid,
    RombergMidpoint,
    Gauss15,
}

impl Rule {
    fn named(method: &Expr) -> Option<Rule> {
        Some(match method.as_symbol()? {
            Symbol::LEFT_RECTANGLE => Rule::Left,
            Symbol::RIGHT_RECTANGLE => Rule::Right,
            Symbol::MIDDLE_POINT => Rule::Middle,
            Symbol::TRAPEZOID => Rule::Trapezoid,
            Symbol::SIMPSON => Rule::Simpson,
            Symbol::ROMBERG_TRAPEZOID => Rule::RombergTrapezoid,
            Symbol::ROMBERG_MIDPOINT => Rule::RombergMidpoint,
            Symbol::GAUSS15 => Rule::Gauss15,
            _ => return None,
        })
    }

    /// The evaluations of `f` for `n`, at most; `None` past `u64`.
    fn nodes(self, n: u64) -> Option<u64> {
        match self {
            Rule::Left | Rule::Right | Rule::Middle => Some(n),
            Rule::Trapezoid => n.checked_add(1),
            Rule::Simpson => n.checked_mul(2)?.checked_add(1),
            // Σ_{j ≤ n} (2^j + 1) ≤ 2^(n+1) + n; n first, so no shift overflows.
            Rule::RombergTrapezoid | Rule::RombergMidpoint => {
                (n < 62).then(|| (1u64 << (n + 1)) + n)
            }
            Rule::Gauss15 => Some(0),
        }
    }

    /// The nodes and weights on `[0, 1]` of a rule that is no Romberg's,
    /// `n ≥ 1`.
    fn points(self, n: i64) -> Vec<(Rational, Rational)> {
        let w = q(1, n);
        match self {
            Rule::Left => (0..n).map(|i| (q(i, n), w.clone())).collect(),
            Rule::Right => (1..=n).map(|i| (q(i, n), w.clone())).collect(),
            Rule::Middle => (0..n).map(|i| (q(2 * i + 1, 2 * n), w.clone())).collect(),
            Rule::Trapezoid => (0..=n)
                .map(|i| {
                    let end = i == 0 || i == n;
                    (q(i, n), if end { q(1, 2 * n) } else { w.clone() })
                })
                .collect(),
            Rule::Simpson => (0..=2 * n)
                .map(|i| {
                    let weight = if i == 0 || i == 2 * n {
                        q(1, 6 * n)
                    } else if i % 2 == 1 {
                        q(2, 3 * n)
                    } else {
                        q(1, 3 * n)
                    };
                    (q(i, 2 * n), weight)
                })
                .collect(),
            Rule::RombergTrapezoid | Rule::RombergMidpoint | Rule::Gauss15 => {
                unreachable!("points of a rule that has none: {self:?}")
            }
        }
    }

    /// The rule's value for `n` (`n ≥ 1`, or `n ≥ 0` for Romberg's);
    /// `None` where an evaluation fails.
    fn apply<S: Sums>(self, sums: &mut S, n: i64) -> Result<Option<S::Value>, EvalError> {
        let base = match self {
            Rule::RombergTrapezoid => Rule::Trapezoid,
            Rule::RombergMidpoint => Rule::Middle,
            _ => return sums.sum(&self.points(n)),
        };
        let mut row: Vec<S::Value> = Vec::new();
        for j in 0..=n {
            let Some(level) = sums.sum(&base.points(1 << j))? else {
                return Ok(None);
            };
            let mut next = vec![level];
            for k in 1..=j {
                let Ok(k) = u64::try_from(k) else {
                    return Ok(None);
                };
                let four = Rational::from(Integer::from(4).pow(k));
                let less = &four - &Rational::from(Integer::from(1));
                let c = &four / &less;
                let d = -&Rational::from(Integer::from(1)) / &less;
                let value = sums.combine(&c, &next[next.len() - 1], &d, &row[next.len() - 1])?;
                next.push(value);
            }
            row = next;
        }
        Ok(row.pop())
    }
}

/// The two ways of computing a rule's sums: balls, or exact expressions.
trait Sums {
    type Value;
    /// `(b − a) Σ wᵢ f(a + tᵢ(b − a))`; `None` where a value fails.
    fn sum(&mut self, points: &[(Rational, Rational)]) -> Result<Option<Self::Value>, EvalError>;
    /// `c u + d v`.
    fn combine(
        &mut self,
        c: &Rational,
        u: &Self::Value,
        d: &Rational,
        v: &Self::Value,
    ) -> Result<Self::Value, EvalError>;
}

/// The sums in ball arithmetic, the tree in the real mode.
struct Balls<'a> {
    tree: &'a Node,
    a: Ball,
    length: Ball,
    wp: u64,
    interrupt: &'a AtomicBool,
}

impl Balls<'_> {
    fn rational(&self, q: &Rational) -> ComplexBall {
        ComplexBall::from_real(Ball::from_rational(q, self.wp))
    }
}

impl Sums for Balls<'_> {
    type Value = ComplexBall;

    fn sum(&mut self, points: &[(Rational, Rational)]) -> Result<Option<ComplexBall>, EvalError> {
        let wp = self.wp;
        let mut total = ComplexBall::ZERO;
        for (t, w) in points {
            if self.interrupt.load(AtomicOrdering::Relaxed) {
                return Err(EvalError::Interrupted);
            }
            let node = self
                .a
                .add(&Ball::from_rational(t, wp).mul(&self.length, wp), wp);
            let Some(value) = self
                .tree
                .eval(&ComplexBall::from_real(node), Mode::Real, wp)
            else {
                return Ok(None);
            };
            total = total.add(&value.mul(&self.rational(w), wp), wp);
        }
        Ok(Some(
            total.mul(&ComplexBall::from_real(self.length.clone()), wp),
        ))
    }

    fn combine(
        &mut self,
        c: &Rational,
        u: &ComplexBall,
        d: &Rational,
        v: &ComplexBall,
    ) -> Result<ComplexBall, EvalError> {
        let wp = self.wp;
        Ok(u.mul(&self.rational(c), wp)
            .add(&v.mul(&self.rational(d), wp), wp))
    }
}

/// The sums as expressions, each evaluated.
struct Exact<'e, 'c> {
    ev: &'e mut Eval<'c>,
    f: &'e Expr,
    x: &'e Expr,
    a: &'e Expr,
    length: Expr,
}

impl Sums for Exact<'_, '_> {
    type Value = Expr;

    fn sum(&mut self, points: &[(Rational, Rational)]) -> Result<Option<Expr>, EvalError> {
        let terms: Vec<Expr> = points
            .iter()
            .map(|(t, w)| {
                let node = Expr::call(
                    Symbol::PLUS,
                    [
                        self.a.clone(),
                        Expr::call(Symbol::TIMES, [Expr::from(t), self.length.clone()]),
                    ],
                );
                let value = self.f.substitute(&[(self.x.clone(), node)]);
                Expr::call(Symbol::TIMES, [Expr::from(w), value])
            })
            .collect();
        let sum = Expr::call(
            Symbol::TIMES,
            [self.length.clone(), Expr::call(Symbol::PLUS, terms)],
        );
        Ok(Some(self.ev.eval(&sum)?))
    }

    fn combine(
        &mut self,
        c: &Rational,
        u: &Expr,
        d: &Rational,
        v: &Expr,
    ) -> Result<Expr, EvalError> {
        self.ev.eval(&Expr::call(
            Symbol::PLUS,
            [
                Expr::call(Symbol::TIMES, [Expr::from(c), u.clone()]),
                Expr::call(Symbol::TIMES, [Expr::from(d), v.clone()]),
            ],
        ))
    }
}

/// The variable and the bounds of `x == a..b` (`Equation` or `Equal`).
fn range(spec: &Expr) -> Option<(&Expr, &Expr, &Expr)> {
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, range]) =
        (spec.head().and_then(Expr::as_symbol), spec.args())
    else {
        return None;
    };
    let (Some(Symbol::RANGE), [a, b]) = (range.head().and_then(Expr::as_symbol), range.args())
    else {
        return None;
    };
    x.as_symbol().is_some().then_some((x, a, b))
}

/// `n` as a count: an integer, or a rational or decimal floored; `None`
/// when negative, past `u64` or no number.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a finite float in [0, 1e18], floored"
)]
fn count(n: &Expr) -> Option<u64> {
    match n.view() {
        View::Int(i) => u64::try_from(i).ok(),
        View::Integer(i) => u64::try_from(i).ok(),
        View::Rational(r) => u64::try_from(&r.floor()).ok(),
        View::Ball(b) => {
            let v = b.mid().to_f64().floor();
            (0.0..=1e18).contains(&v).then_some(v as u64)
        }
        _ => None,
    }
}

/// Whether `e` holds a symbol other than `x` that is not numeric, heads
/// aside.
fn has_parameter(e: &Expr, x: &Expr) -> bool {
    let free = |s: &Expr| s.as_symbol().is_some() && s != x && !s.is_numeric();
    free(e) || e.distinct().any(|s| s.args().iter().any(free))
}

pub(in crate::native) fn area(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, spec, rule) = match e.args() {
        [f, spec] => (f, spec, None),
        [f, spec, n, method] => (f, spec, Some((n, method))),
        _ => return Ok(None),
    };
    let Some((x, a, b)) = range(spec) else {
        return Ok(None);
    };
    let bounds = || list([x.clone(), a.clone(), b.clone()]);
    let Some((n, method)) = rule else {
        let integral = ev.eval(&Expr::call(Symbol::INTEGRATE, [f.clone(), bounds()]))?;
        let stays = integral.head() == Some(&Expr::from(Symbol::INTEGRATE));
        return Ok((!stays).then_some(integral));
    };
    let (Some(rule), Some(n)) = (Rule::named(method), count(n)) else {
        return Ok(None);
    };
    if rule == Rule::Gauss15 {
        return super::nintegrate(ev, &Expr::call(Symbol::NINTEGRATE, [f.clone(), bounds()]));
    }
    let limit = MAX_EVALUATIONS.min(ev.context().limits().terms);
    let within = rule
        .nodes(n)
        .and_then(|k| usize::try_from(k).ok())
        .is_some_and(|k| k <= limit);
    let (true, Ok(n)) = (within, i64::try_from(n)) else {
        return Ok(None);
    };
    if n == 0 {
        match rule {
            Rule::Trapezoid | Rule::Simpson => {
                ev.warn(Warning::Indeterminate(e.clone()));
                return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
            }
            Rule::Left | Rule::Right | Rule::Middle => return Ok(Some(Expr::from(0))),
            _ => {}
        }
    }
    let wp = ev.context().working_precision() + GUARD;
    let tree = Node::compile(ev.context(), f, x, wp);
    if let (Some(tree), Some(lo), Some(hi)) = (tree, bound(ev, a, wp), bound(ev, b, wp)) {
        let mut sums = Balls {
            tree: &tree,
            length: hi.sub(&lo, wp),
            a: lo,
            wp,
            interrupt: ev.context().interrupt_flag(),
        };
        let Some(value) = rule.apply(&mut sums, n)? else {
            return Ok(None);
        };
        let result = if value.is_real() {
            Expr::from(value.re())
        } else {
            Expr::from(&value)
        };
        ev.step("Area", e, &result);
        return Ok(Some(result));
    }
    if ![f, a, b].iter().any(|part| has_parameter(part, x)) {
        return Ok(None);
    }
    let length = Expr::call(
        Symbol::PLUS,
        [
            b.clone(),
            Expr::call(Symbol::TIMES, [Expr::from(-1), a.clone()]),
        ],
    );
    let mut sums = Exact {
        ev,
        f,
        x,
        a,
        length,
    };
    rule.apply(&mut sums, n)
}

#[cfg(test)]
mod tests;
