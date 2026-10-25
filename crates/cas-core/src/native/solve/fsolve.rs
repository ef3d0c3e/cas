//! `FSolve` (CAS-DOC-11 §7, CAS-214): numeric roots, Xcas's `fsolve`, as
//! balls at the working precision that each hold a root.
//!
//! - `FSolve[eqn, x]`: of a rational function of `x` alone with exact
//!   coefficients, the real roots of its numerator's square-free factors by
//!   `polynomial_roots` (§3.4), up to degree [`MAX_NUMERIC_DEGREE`]: the
//!   number alone when there is one, as giac's `fsolve(x^5-x-1,x)` is
//!   `1.16730397826`, else a list. Of another equation, giac's bisection
//!   after `x = tan(t)`: the sign changes of `f` on [`GRID`] points `tan(t)`,
//!   `t` evenly in `(−π/2, π/2)`, each bisected, as a list, the first
//!   [`MAX_SCANNED`] of them.
//! - `FSolve[eqn, x, x0]`: Newton's iteration from `x0`, the derivative by
//!   `D`, at the working precision and 32 bits more: the root alone.
//! - `FSolve[eqn, Equation[x, Range[a, b]]]`, giac's `fsolve(eqn,x=a..b)`:
//!   the roots in `[a, b]`, of a polynomial as above, else the sign changes
//!   on [`GRID`] points of `[a, b]`, each bisected: a list, `[]` without one.
//! - A root is certified: by `polynomial_roots`, or by values of opposite
//!   signs at the ends of an interval over which `f` is a finite ball, so
//!   that no pole is taken for a root (`tan`). Newton's limit is certified
//!   so, by a sign change across a few of its last bits; when it is not,
//!   the call stays with a warning.
//! - With parameters, `f` no number at a point, the call stays.

use std::cmp::Ordering;

use cas_num::algebra::Integers;
use cas_num::ball::{Ball, Float, Round, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::super::approximation;
use super::{MAX_NUMERIC_DEGREE, MAX_ROOT_BITS, sides};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{Bounds, Fraction, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The points a whole line or an interval is scanned at for sign changes.
const GRID: usize = 256;

/// The most steps of Newton's iteration.
const MAX_NEWTON: usize = 100;

/// The most roots a scan bisects: each takes a value of `f` a step.
const MAX_SCANNED: usize = 16;

/// The most steps of a bisection, past which its ball is as narrow as
/// they made it.
const MAX_BISECTION: u64 = 1024;

/// Where the roots are sought.
enum Where<'a> {
    Everywhere,
    From(&'a Expr),
    Between(&'a Expr, &'a Expr),
}

/// `FSolve[eqn, x]`, `FSolve[eqn, x, x0]` or `FSolve[eqn, x == a .. b]`.
pub(in crate::native) fn fsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (eqn, x, place) = match e.args() {
        [eqn, spec] => match sides(spec) {
            Some((x, r)) => match range(r) {
                Some((a, b)) => (eqn, x, Where::Between(a, b)),
                None => return Ok(None),
            },
            None => (eqn, spec, Where::Everywhere),
        },
        [eqn, x, guess] => match range(guess) {
            Some((a, b)) => (eqn, x, Where::Between(a, b)),
            None => (eqn, x, Where::From(guess)),
        },
        _ => return Ok(None),
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let f = match sides(eqn) {
        Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
        None => eqn.clone(),
    };
    if !f.contains(x) || !numeric(ev.context(), &f, x) {
        return Ok(None);
    }
    let bits = ev.context().working_precision();
    let result = match place {
        Where::Everywhere => match polynomial(ev, &f, x, bits)? {
            Some(mut roots) if roots.len() == 1 => roots.pop().expect("one root"),
            Some(roots) => list(roots),
            None => list(scanned(ev.context(), &f, x, &tangent_grid(), bits)),
        },
        Where::Between(a, b) => {
            let (Some(a), Some(b)) = (point(ev.context(), a), point(ev.context(), b)) else {
                return Ok(None);
            };
            if a > b {
                return Ok(None);
            }
            match polynomial(ev, &f, x, bits)? {
                Some(roots) => list(roots.into_iter().filter(|r| within(r, &a, &b)).collect()),
                None => list(scanned(ev.context(), &f, x, &interval_grid(&a, &b), bits)),
            }
        }
        Where::From(guess) => {
            let Some(start) = point(ev.context(), guess) else {
                return Ok(None);
            };
            let Some(root) = newton(ev, &f, x, start, bits)? else {
                ev.warn(Warning::Native {
                    message: "Newton's iteration found no root it could certify from".into(),
                    expr: guess.clone(),
                });
                return Ok(None);
            };
            root
        }
    };
    ev.step("FSolve", e, &result);
    Ok(Some(result))
}

/// Whether `f` is a number at `x = 1/3` or `x = 2/7`: no parameters.
fn numeric(ctx: &Context, f: &Expr, x: &Expr) -> bool {
    [(1, 3), (2, 7)].iter().any(|&(p, q)| {
        let at = Expr::from(Rational::new(Integer::from(p), Integer::from(q)));
        matches!(
            approximation(ctx, &f.substitute(&[(x.clone(), at)]), 64)
                .as_ref()
                .map(Expr::view),
            Ok(View::Ball(_) | View::ComplexBall(_))
        )
    })
}

/// The ends of `Range[a, b]`, the language's `a .. b`.
fn range(r: &Expr) -> Option<(&Expr, &Expr)> {
    match (r.head().and_then(Expr::as_symbol), r.args()) {
        (Some(Symbol::RANGE), [a, b]) => Some((a, b)),
        _ => None,
    }
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// A real number's value as a rational: exact, or its ball's midpoint.
fn point(ctx: &Context, e: &Expr) -> Option<Rational> {
    if let Some(q) = e.to_rational() {
        return Some(q);
    }
    match approximation(ctx, e, 128).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(b.mid().to_rational()),
        _ => None,
    }
}

/// Whether the root's ball has its midpoint in `[a, b]`.
fn within(root: &Expr, a: &Rational, b: &Rational) -> bool {
    match root.view() {
        View::Ball(ball) => {
            let m = ball.mid().to_rational();
            *a <= m && m <= *b
        }
        _ => false,
    }
}

/// The real roots of `f` when it is a rational function of `x` alone with
/// exact coefficients, ascending; `None` when it is not one, or a factor
/// is past [`MAX_NUMERIC_DEGREE`] or its roots are not validated.
fn polynomial(
    ev: &Eval<'_>,
    f: &Expr,
    x: &Expr,
    bits: u64,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let ctx = ev.context();
    let cancel: &dyn Cancel = ctx.interrupt_flag();
    let bounds = Bounds {
        terms: ctx.limits().terms,
        bits: ctx.limits().number_bits,
    };
    let Ok(p) = Polynomials::new(std::slice::from_ref(f), bounds, cancel) else {
        return Ok(None);
    };
    if p.generators() != std::slice::from_ref(x) {
        return Ok(None);
    }
    let Fraction::Integers { num, den } = &p.fractions()[0] else {
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, 1);
    let max_bits = ctx
        .limits()
        .number_bits
        .min(bits.saturating_mul(16).max(MAX_ROOT_BITS))
        .max(bits);
    let Some(mut roots) = real_roots(&ring, num, den, bits, max_bits, cancel)? else {
        return Ok(None);
    };
    roots.sort_by(|a, b| a.mid().cmp(b.mid()));
    Ok(Some(roots.into_iter().map(Expr::from).collect()))
}

/// The real roots of `num / den` in its one variable, as balls.
fn real_roots(
    ring: &SparseRing<Integers>,
    num: &Sparse<Integer>,
    den: &Sparse<Integer>,
    bits: u64,
    max_bits: u64,
    cancel: &dyn Cancel,
) -> Result<Option<Vec<Ball>>, EvalError> {
    let (_, num, _) = ring
        .gcd_cofactors_with(num, den, cancel)
        .map_err(|_| EvalError::Interrupted)?;
    let parts = ring
        .squarefree_with(&num, cancel)
        .map_err(|_| EvalError::Interrupted)?;
    let mut roots = Vec::new();
    for (part, _) in parts.factors() {
        let Some(dense) = ring.to_dense(part, 0) else {
            return Ok(None);
        };
        let degree = dense.coeffs().len().saturating_sub(1);
        if degree == 0 {
            continue;
        }
        if u64::try_from(degree).map_or(true, |d| d > MAX_NUMERIC_DEGREE) {
            return Ok(None);
        }
        let coeffs: Vec<Complex> = dense.coeffs().iter().cloned().map(Complex::from).collect();
        let Some(found) = polynomial_roots(&coeffs, bits, max_bits) else {
            return Ok(None);
        };
        roots.extend(
            found
                .into_iter()
                .filter(cas_num::ball::ComplexBall::is_real)
                .map(|z| z.re().round_to_accuracy(bits)),
        );
    }
    Ok(Some(roots))
}

/// `GRID` points `tan(t)`, `t` evenly inside `(−π/2, π/2)`, as rationals.
fn tangent_grid() -> Vec<Rational> {
    #[allow(clippy::cast_precision_loss, reason = "a small grid")]
    let n = GRID as f64;
    (0..GRID)
        .filter_map(|i| {
            #[allow(clippy::cast_precision_loss, reason = "a small grid")]
            let t = -std::f64::consts::FRAC_PI_2 + (i as f64 + 0.5) * std::f64::consts::PI / n;
            Float::from_f64(t.tan()).map(|x| x.to_rational())
        })
        .collect()
}

/// `GRID + 1` points evenly from `a` to `b`.
fn interval_grid(a: &Rational, b: &Rational) -> Vec<Rational> {
    let n = Rational::from(Integer::from(GRID));
    let width = &(b - a) / &n;
    (0..=GRID)
        .map(|i| a + &(&width * &Rational::from(Integer::from(i))))
        .collect()
}

/// `f` at `x = q`, a real ball at `bits`; `None` for no finite real value.
fn value(ctx: &Context, f: &Expr, x: &Expr, q: &Expr, bits: u64) -> Option<Ball> {
    let v = approximation(ctx, &f.substitute(&[(x.clone(), q.clone())]), bits).ok()?;
    match v.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        View::ComplexBall(z) if z.is_finite() && z.im().contains_rational(&Rational::ZERO) => {
            Some(z.re().clone())
        }
        _ => None,
    }
}

/// The sign of `f` at the rational `q`, when its ball tells it.
fn sign(ctx: &Context, f: &Expr, x: &Expr, q: &Rational, bits: u64) -> Option<Ordering> {
    value(ctx, f, x, &Expr::from(q.clone()), bits)?.sign()
}

/// The roots at the sign changes of `f` between neighbours of `grid`,
/// each bisected and certified.
fn scanned(ctx: &Context, f: &Expr, x: &Expr, grid: &[Rational], bits: u64) -> Vec<Expr> {
    let signs: Vec<Option<Ordering>> = grid.iter().map(|q| sign(ctx, f, x, q, bits)).collect();
    let mut roots = Vec::new();
    for (i, q) in grid.iter().enumerate() {
        if signs[i] == Some(Ordering::Equal) {
            roots.extend(approximation(ctx, &Expr::from(q.clone()), bits));
        }
        if roots.len() >= MAX_SCANNED {
            break;
        }
        if i + 1 < grid.len()
            && let (Some(s), Some(t)) = (signs[i], signs[i + 1])
            && s != Ordering::Equal
            && t != Ordering::Equal
            && s != t
            && let Some(root) = bisected(ctx, f, x, grid[i].clone(), grid[i + 1].clone(), s, bits)
        {
            roots.push(root);
        }
    }
    roots
}

/// The root of `f` between `lo` and `hi`, where `f` has the sign `s_lo` at
/// `lo` and the other one at `hi`, by bisection to `bits` bits: a ball, when
/// `f` over the last interval is a finite ball.
fn bisected(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    mut lo: Rational,
    mut hi: Rational,
    s_lo: Ordering,
    bits: u64,
) -> Option<Expr> {
    let two = Rational::from(Integer::from(2));
    let steps = bits
        .saturating_mul(4)
        .saturating_add(256)
        .min(MAX_BISECTION);
    for _ in 0..steps {
        if narrow(&lo, &hi, bits) {
            break;
        }
        let mid = &(&lo + &hi) / &two;
        match sign(ctx, f, x, &mid, bits + 32) {
            Some(Ordering::Equal) => return approximation(ctx, &Expr::from(mid), bits).ok(),
            Some(s) if s == s_lo => lo = mid,
            Some(_) => hi = mid,
            // Within rounding of 0: as narrow as `bits + 32` tells.
            None => break,
        }
    }
    certified(ctx, f, x, &lo, &hi, bits)
}

/// Whether `[lo, hi]` is within `bits` bits of its values, or of `2^−2bits`
/// about 0.
fn narrow(lo: &Rational, hi: &Rational, bits: u64) -> bool {
    let width = hi - lo;
    let scale = lo.abs().max(hi.abs());
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let tiny = Rational::from(Integer::from(2)).pow(-exponent.saturating_mul(2));
    let relative = &scale * &Rational::from(Integer::from(2)).pow(-exponent);
    width <= relative || width <= tiny
}

/// The ball `[lo, hi]` as a root of `f`, when `f` over it is a finite
/// ball: no pole inside.
fn certified(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    lo: &Rational,
    hi: &Rational,
    bits: u64,
) -> Option<Expr> {
    let two = Rational::from(Integer::from(2));
    let mid = &(lo + hi) / &two;
    let rad = &(hi - lo) / &two;
    let interval = Expr::call(Symbol::BALL, [Expr::from(mid), Expr::from(rad)]);
    let ball = approximation(ctx, &interval, bits.saturating_add(32)).ok()?;
    value(ctx, f, x, &ball, bits.saturating_add(32))?;
    match ball.view() {
        View::Ball(b) => Some(Expr::from(b.round_to_accuracy(bits))),
        _ => None,
    }
}

/// Newton's iteration on `f` from `start`, its limit certified by a sign
/// change around it.
fn newton(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    start: Rational,
    bits: u64,
) -> Result<Option<Expr>, EvalError> {
    let df = ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?;
    let ctx = ev.context();
    let prec = bits.saturating_add(32);
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let mut q = start;
    for _ in 0..MAX_NEWTON {
        let at = Expr::from(q.clone());
        let (Some(v), Some(d)) = (value(ctx, f, x, &at, prec), value(ctx, &df, x, &at, prec))
        else {
            return Ok(None);
        };
        if d.mid().is_zero() {
            return Ok(None);
        }
        let step = &v.mid().to_rational() / &d.mid().to_rational();
        let next = &q - &step;
        let (rounded, _) = Float::from_rational(&next, prec, Round::Nearest);
        let next = rounded.to_rational();
        let small = step.abs()
            <= &next.abs() * &Rational::from(Integer::from(2)).pow(-exponent.saturating_add(8));
        q = next;
        if small || step.is_zero() {
            return Ok(around(ctx, f, x, &q, bits));
        }
    }
    Ok(None)
}

/// A ball about `q` that holds a root of `f`: `f` of opposite signs at
/// `q ± δ`, `δ` a few of `q`'s last bits, then widened sixteenfold up to
/// four times.
fn around(ctx: &Context, f: &Expr, x: &Expr, q: &Rational, bits: u64) -> Option<Expr> {
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let two = Rational::from(Integer::from(2));
    let scale = q.abs().max(two.pow(-exponent));
    let mut delta = &scale * &two.pow(4 - exponent);
    for _ in 0..4 {
        let (lo, hi) = (q - &delta, q + &delta);
        let prec = bits.saturating_add(32);
        if let (Some(s), Some(t)) = (sign(ctx, f, x, &lo, prec), sign(ctx, f, x, &hi, prec))
            && (s != t || s == Ordering::Equal)
        {
            return certified(ctx, f, x, &lo, &hi, bits);
        }
        delta = &delta * &Rational::from(Integer::from(16));
    }
    None
}
