//! Two-point boundary value problems (CAS-520): `BVPSolve`, Xcas's
//! `bvpsolve`, `y″ = f(x, y, y′)` on `[a, b]` with `y(a)` and `y(b)` given.
//!
//! `BVPSolve[f, {x == a .. b, y}, {ya, yb}, n]`, giac's
//! `bvpsolve(f, [x=a..b, y], [ya, yb], n)`:
//!
//! - `f` is written in `x`, `y` (or `y(x)`) and `y′` as `diff(y(x), x)`;
//!   it compiles as `NDSolve`'s right sides ([`rhs`](super::rhs)) for the
//!   system `y′ = p`, `p′ = f`. A third boundary value is a guess of
//!   `y′(a)`, `n` (default [`DEFAULT_POINTS`], at most [`MAX_POINTS`])
//!   the number of intervals, and `output = list` the only output option;
//!   giac's `piecewise` and `spline` leave the call.
//! - By shooting: the secant method on the slope `s = y′(a)` until
//!   `y(b; s) = yb`, each `y(b; s)` by [`taylor`](super::taylor)'s
//!   integration on balls; then two slopes about the limit at which
//!   `y(b; s) − yb` has certainly opposite signs, so that a solution's
//!   slope is between them (when `y(b; s)` is continuous, as it is where
//!   the integration succeeds).
//! - The curve is integrated from that interval of slopes, so each value
//!   holds the solution's: `{{x_k, y_k}, …}` at `n + 1` points evenly from
//!   `a` to `b`, the ends' values the given ones, as giac's list
//!   `[[1.0,17.0],[1.1,15.7554961579],…]` (its values good to about 10⁻⁵,
//!   ours to the working precision).
//! - A shooting that does not converge, or an integration that fails,
//!   leaves the call with a warning.

#![allow(
    clippy::many_single_char_names,
    reason = "the problem's names: f, x, y, a, b, the slope s and n"
)]

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::rhs::Node;
use super::{real, taylor};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The number of intervals without `n`, giac's.
pub(super) const DEFAULT_POINTS: i64 = 100;

/// The most intervals.
pub(super) const MAX_POINTS: i64 = 2000;

/// The most secant steps.
const MAX_SHOTS: usize = 60;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// What a call asks.
struct Problem {
    f: Expr,
    x: Expr,
    y: Expr,
    a: Expr,
    b: Expr,
    ya: Expr,
    yb: Expr,
    slope: Option<Expr>,
    intervals: i64,
}

fn problem(e: &Expr) -> Option<Problem> {
    let [f, spec, conditions, rest @ ..] = e.args() else {
        return None;
    };
    let (Some(Symbol::LIST), [span, y]) = (head(spec), spec.args()) else {
        return None;
    };
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, range]) = (head(span), span.args()) else {
        return None;
    };
    let (Some(Symbol::RANGE), [a, b]) = (head(range), range.args()) else {
        return None;
    };
    if x.as_symbol().is_none() || y.as_symbol().is_none() || x == y {
        return None;
    }
    let (ya, yb, slope) = match (head(conditions), conditions.args()) {
        (Some(Symbol::LIST), [ya, yb]) => (ya, yb, None),
        (Some(Symbol::LIST), [ya, yb, s]) => {
            let s = match (head(s), s.args()) {
                (Some(Symbol::LIST), [s]) => s,
                _ => s,
            };
            (ya, yb, Some(s.clone()))
        }
        _ => return None,
    };
    let mut intervals = DEFAULT_POINTS;
    for option in rest {
        if let Some(n) = option.to_rational() {
            let n = n
                .is_integer()
                .then(|| n.numer().to_i64())
                .flatten()
                .filter(|n| (1..=MAX_POINTS).contains(n))?;
            intervals = n;
            continue;
        }
        match (head(option), option.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r])
                if *l == Expr::symbol("output") && *r == Expr::symbol("list") => {}
            _ => return None,
        }
    }
    Some(Problem {
        f: f.clone(),
        x: x.clone(),
        y: y.clone(),
        a: a.clone(),
        b: b.clone(),
        ya: ya.clone(),
        yb: yb.clone(),
        slope,
        intervals,
    })
}

/// `BVPSolve[…]`.
pub(in crate::native) fn bvp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = problem(e) else {
        return Ok(None);
    };
    let (big_y, big_p) = (Expr::symbol("BVPSolve$y"), Expr::symbol("BVPSolve$p"));
    // y′ as `diff(y(x), x)` reads, `Derivative[1][y][x]`, then y(x), then y.
    let derivative = Expr::apply(
        &Expr::apply(
            &Expr::call(Symbol::DERIVATIVE, [Expr::from(1)]),
            std::slice::from_ref(&p.y),
        ),
        std::slice::from_ref(&p.x),
    );
    let f =
        p.f.substitute(&[(derivative, big_p.clone())])
            .substitute(&[(Expr::apply(&p.y, std::slice::from_ref(&p.x)), big_y.clone())])
            .substitute(&[(p.y.clone(), big_y.clone())]);
    let f = ev.eval(&f)?;
    if f.contains(&p.y) {
        return Ok(None);
    }
    let prec = ev.context().working_precision();
    let q = prec + 32;
    let variables = [p.x.clone(), big_y.clone(), big_p.clone()];
    let (Some(fy), Some(fp)) = (
        Node::compile(ev.context(), &big_p, &variables, q),
        Node::compile(ev.context(), &f, &variables, q),
    ) else {
        return Ok(None);
    };
    let system = [fy, fp];
    let ends = (
        real(ev, &p.a, q),
        real(ev, &p.b, q),
        real(ev, &p.ya, q),
        real(ev, &p.yb, q),
    );
    let (Some(a), Some(b), Some(ya), Some(yb)) = ends else {
        return Ok(None);
    };
    if b.sub(&a, q).sign() != Some(std::cmp::Ordering::Greater) {
        return Ok(None);
    }
    let guess = match &p.slope {
        Some(s) => match real(ev, s, q) {
            Some(s) => s.mid().to_rational(),
            None => return Ok(None),
        },
        None => yb.sub(&ya, q).div(&b.sub(&a, q), q).mid().to_rational(),
    };
    let shoot = Shooting {
        system: &system,
        a: &a,
        b: &b,
        ya: &ya,
        yb: &yb,
        prec,
    };
    let Some(slopes) = shoot.slopes(&guess) else {
        ev.warn(Warning::Native {
            message: "the shooting found no slope y'(a) that meets y(b)".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    let Some(values) = shoot.curve(&slopes, p.intervals) else {
        ev.warn(Warning::Native {
            message: "the integration failed: a solution that blows up, or too many steps".into(),
            expr: e.clone(),
        });
        return Ok(None);
    };
    let ya = ev.eval(&p.ya)?;
    let yb = ev.eval(&p.yb)?;
    let exact = (p.a.to_rational(), p.b.to_rational());
    let n = values.len();
    let mut rows = Vec::with_capacity(n + 1);
    for (k, value) in values.iter().enumerate() {
        let x = match &exact {
            (Some(a), Some(b)) => {
                let t = Rational::new(Integer::from(k), Integer::from(n - 1));
                Expr::from(a + &(&(b - a) * &t))
            }
            _ => Expr::from(value.0.round_to_accuracy(prec)),
        };
        let y = if k == 0 {
            ya.clone()
        } else if k + 1 == n {
            yb.clone()
        } else {
            Expr::from(value.1.round_to_accuracy(prec))
        };
        rows.push(list(vec![x, y]));
    }
    let result = list(rows);
    ev.step("BVPSolve", e, &result);
    Ok(Some(result))
}

/// A shooting problem: `y′ = p`, `p′ = f` from `a` to `b`.
struct Shooting<'a> {
    system: &'a [Node],
    a: &'a Ball,
    b: &'a Ball,
    ya: &'a Ball,
    yb: &'a Ball,
    prec: u64,
}

impl Shooting<'_> {
    /// `y(b; s) − yb` for the slope ball `s`.
    fn miss(&self, s: &Ball) -> Option<Ball> {
        let (end, _) = taylor::integrate(
            self.system,
            self.a,
            &[self.ya.clone(), s.clone()],
            self.b,
            self.prec,
        )?;
        let miss = end[0].sub(self.yb, self.prec + 32);
        miss.is_finite().then_some(miss)
    }

    fn at(&self, s: &Rational) -> Option<Ball> {
        self.miss(&Ball::from_rational(s, self.prec + 32))
    }

    /// The secant method from `guess`, then an interval of slopes about
    /// its limit across which the miss certainly changes sign.
    fn slopes(&self, guess: &Rational) -> Option<Ball> {
        let prec = self.prec + 32;
        let one = Rational::from(Integer::from(1));
        let step = guess.abs().max(one.clone()) / Rational::from(Integer::from(16));
        let (mut s0, mut s1) = (guess.clone(), guess + &step);
        let mut m0 = self.at(&s0)?.mid().to_rational();
        let exponent = i64::try_from(self.prec).unwrap_or(i64::MAX);
        let two = Rational::from(Integer::from(2));
        for _ in 0..MAX_SHOTS {
            let miss = self.at(&s1)?;
            if miss.contains_rational(&Rational::ZERO) {
                break;
            }
            let m1 = miss.mid().to_rational();
            let den = &m1 - &m0;
            if den.is_zero() {
                return None;
            }
            let next = &s1 - &(&(&m1 * &(&s1 - &s0)) / &den);
            let next = Ball::from_rational(&next, prec).mid().to_rational();
            let small =
                (&next - &s1).abs() <= &next.abs().max(one.clone()) * &two.pow(-exponent - 8);
            (s0, m0, s1) = (s1, m1, next);
            if small {
                break;
            }
        }
        let mut delta = &s1.abs().max(one) * &two.pow(8 - exponent);
        for _ in 0..6 {
            let (lo, hi) = (&s1 - &delta, &s1 + &delta);
            if let (Some(l), Some(h)) = (self.at(&lo)?.sign(), self.at(&hi)?.sign())
                && l != std::cmp::Ordering::Equal
                && h != std::cmp::Ordering::Equal
                && l != h
            {
                let mid = Ball::from_rational(&s1, prec);
                let rad = Ball::from_rational(&delta, prec);
                return Some(mid.add_error(&cas_num::ball::Mag::from_float(rad.mid())));
            }
            delta = &delta * &Rational::from(Integer::from(16));
        }
        None
    }

    /// The points `(x_k, y_k)` from `a` to `b`, `n` intervals, from the
    /// slopes `s`.
    fn curve(&self, s: &Ball, n: i64) -> Option<Vec<(Ball, Ball)>> {
        let prec = self.prec + 32;
        let width = self.b.sub(self.a, prec);
        let mut state = vec![self.ya.clone(), s.clone()];
        let mut t = self.a.clone();
        let mut points = vec![(t.clone(), state[0].clone())];
        for k in 1..=n {
            let fraction =
                Ball::from_rational(&Rational::new(Integer::from(k), Integer::from(n)), prec);
            let next = if k == n {
                self.b.clone()
            } else {
                self.a.add(&width.mul(&fraction, prec), prec)
            };
            let (y, _) = taylor::integrate(self.system, &t, &state, &next, self.prec)?;
            if !y.iter().all(Ball::is_finite) {
                return None;
            }
            state = y;
            t = next;
            points.push((t.clone(), state[0].clone()));
        }
        Some(points)
    }
}
