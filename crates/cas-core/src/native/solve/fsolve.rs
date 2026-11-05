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
//! - A last argument names giac's method (CAS-517), each answer certified
//!   as above whatever found it. The bracketing ones, `bisection_solver`,
//!   `falsepos_solver` (Illinois) and `brent_solver`, refine each sign
//!   change of the scan, on the range or the whole line, into a list, as
//!   giac's `fsolve(cos(x)=x,x=0..1,brent_solver)` is `[0.739085133215]`;
//!   with a start point, giac's error, the call stays with a warning. A
//!   false position or Brent estimate that is not certified falls back to
//!   bisection. The open ones, `newton_solver`, `secant_solver` and
//!   `steffenson_solver`, iterate from the start point, the range's
//!   midpoint (a list, `[]` when the root is outside) or 0 for the whole
//!   line; the secant from `x0` and `x0 + 1/16`, or from the two points of
//!   `[x0, x1]`, which a bracketing method takes as a range.

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
#[derive(Clone, Copy)]
enum Where<'a> {
    Everywhere,
    From(&'a Expr),
    Between(&'a Expr, &'a Expr),
    /// `[x0, x1]`, given with a method: the secant's two points, or a
    /// range.
    Pair(&'a Expr, &'a Expr),
}

/// giac's root finders, the last argument's names.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Method {
    Bisection,
    FalsePosition,
    Brent,
    Newton,
    Secant,
    Steffensen,
}

impl Method {
    const NAMES: [(&'static str, Self); 6] = [
        ("bisection_solver", Self::Bisection),
        ("falsepos_solver", Self::FalsePosition),
        ("brent_solver", Self::Brent),
        ("newton_solver", Self::Newton),
        ("secant_solver", Self::Secant),
        ("steffenson_solver", Self::Steffensen),
    ];

    fn named(e: &Expr) -> Option<Self> {
        e.as_symbol()?;
        Self::NAMES
            .iter()
            .find(|(name, _)| *e == Expr::symbol(name))
            .map(|&(_, m)| m)
    }

    fn bracketing(self) -> bool {
        matches!(self, Self::Bisection | Self::FalsePosition | Self::Brent)
    }
}

/// `FSolve[eqn, x]`, `FSolve[eqn, x, x0]` or `FSolve[eqn, x == a .. b]`,
/// each with a method last.
pub(in crate::native) fn fsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (args, method) = match e.args() {
        [rest @ .., last] if rest.len() >= 2 => match Method::named(last) {
            Some(m) => (rest, Some(m)),
            None => (e.args(), None),
        },
        args => (args, None),
    };
    let (eqn, x, place) = match args {
        [eqn, spec] => match sides(spec) {
            Some((x, r)) => match range(r) {
                Some((a, b)) => (eqn, x, Where::Between(a, b)),
                None => return Ok(None),
            },
            None => (eqn, spec, Where::Everywhere),
        },
        [eqn, x, guess] => match (range(guess), guess.args()) {
            (Some((a, b)), _) => (eqn, x, Where::Between(a, b)),
            (None, [a, b]) if method.is_some() && head(guess) == Some(Symbol::LIST) => {
                (eqn, x, Where::Pair(a, b))
            }
            _ => (eqn, x, Where::From(guess)),
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
    let result = match method {
        None => match unmethodical(ev, &f, x, place, bits)? {
            Some(r) => r,
            None => return Ok(None),
        },
        Some(m) if m.bracketing() => match bracketed(ev, &f, x, place, m, bits) {
            Some(r) => r,
            None => return Ok(None),
        },
        Some(m) => match opened(ev, &f, x, place, m, bits)? {
            Some(r) => r,
            None => return Ok(None),
        },
    };
    ev.step("FSolve", e, &result);
    Ok(Some(result))
}

/// The answer without a method: polynomial roots, a scan, or Newton's.
fn unmethodical(ev: &mut Eval<'_>, f: &Expr, x: &Expr, place: Where<'_>, bits: u64) -> Native {
    Ok(Some(match place {
        Where::Everywhere => match polynomial(ev, f, x, bits)? {
            Some(mut roots) if roots.len() == 1 => roots.pop().expect("one root"),
            Some(roots) => list(roots),
            None => list(scanned(
                ev.context(),
                f,
                x,
                &tangent_grid(),
                Method::Bisection,
                bits,
            )),
        },
        Where::Between(a, b) => {
            let Some((a, b)) = ends(ev.context(), a, b) else {
                return Ok(None);
            };
            match polynomial(ev, f, x, bits)? {
                Some(roots) => list(roots.into_iter().filter(|r| within(r, &a, &b)).collect()),
                None => list(scanned(
                    ev.context(),
                    f,
                    x,
                    &interval_grid(&a, &b),
                    Method::Bisection,
                    bits,
                )),
            }
        }
        Where::From(guess) => match opened(ev, f, x, Where::From(guess), Method::Newton, bits)? {
            Some(root) => root,
            None => return Ok(None),
        },
        Where::Pair(..) => return Ok(None),
    }))
}

/// The ends of a range as rationals, in order.
fn ends(ctx: &Context, a: &Expr, b: &Expr) -> Option<(Rational, Rational)> {
    let (a, b) = (point(ctx, a)?, point(ctx, b)?);
    (a <= b).then_some((a, b))
}

/// A bracketing method's roots: each sign change of the scan refined.
fn bracketed(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    place: Where<'_>,
    method: Method,
    bits: u64,
) -> Option<Expr> {
    let ctx = ev.context();
    let grid = match place {
        Where::Everywhere => tangent_grid(),
        Where::Between(a, b) | Where::Pair(a, b) => {
            let (a, b) = ends(ctx, a, b)?;
            interval_grid(&a, &b)
        }
        Where::From(guess) => {
            ev.warn(Warning::Native {
                message: "a bracketing method needs a range, not a start point".into(),
                expr: guess.clone(),
            });
            return None;
        }
    };
    Some(list(scanned(ctx, f, x, &grid, method, bits)))
}

/// An open method's root: a number from a start point or the whole line,
/// a list of it, or `[]`, on a range.
fn opened(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    place: Where<'_>,
    method: Method,
    bits: u64,
) -> Native {
    let ctx = ev.context();
    let half = Rational::new(Integer::from(1), Integer::from(2));
    let (q0, q1, bounds, start) = match place {
        Where::Everywhere => (None, Rational::ZERO, None, e_zero()),
        Where::From(guess) => match point(ctx, guess) {
            Some(q) => (None, q, None, guess.clone()),
            None => return Ok(None),
        },
        Where::Pair(a, b) => match (point(ctx, a), point(ctx, b)) {
            (Some(p), Some(q)) if method == Method::Secant => (Some(p), q, None, a.clone()),
            (Some(p), Some(q)) => (None, &(&p + &q) * &half, None, a.clone()),
            _ => return Ok(None),
        },
        Where::Between(a, b) => match ends(ctx, a, b) {
            Some((p, q)) if method == Method::Secant => {
                (Some(p.clone()), q.clone(), Some((p, q)), a.clone())
            }
            Some((p, q)) => (None, &(&p + &q) * &half, Some((p, q)), a.clone()),
            None => return Ok(None),
        },
    };
    let df = match method {
        Method::Newton => Some(ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?),
        _ => None,
    };
    let ctx = ev.context();
    let q0 = q0.unwrap_or_else(|| &q1 + &Rational::new(Integer::from(1), Integer::from(16)));
    let root = iterated(ctx, f, df.as_ref(), x, method, q0, q1, bits)
        .and_then(|q| around(ctx, f, x, &q, bits));
    let Some(root) = root else {
        ev.warn(Warning::Native {
            message: "the iteration found no root it could certify from".into(),
            expr: start,
        });
        return Ok(None);
    };
    Ok(Some(match bounds {
        None => root,
        Some((a, b)) if within(&root, &a, &b) => list(vec![root]),
        Some(_) => list(Vec::new()),
    }))
}

fn e_zero() -> Expr {
    Expr::from(Integer::from(0))
}

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
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

/// The sign of `f` at the rational `q`, when its ball tells it. Also the
/// monotone inversion's (CAS-302).
pub(in crate::native) fn sign(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    q: &Rational,
    bits: u64,
) -> Option<Ordering> {
    value(ctx, f, x, &Expr::from(q.clone()), bits)?.sign()
}

/// The roots at the sign changes of `f` between neighbours of `grid`,
/// each bisected and certified.
fn scanned(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    grid: &[Rational],
    method: Method,
    bits: u64,
) -> Vec<Expr> {
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
            && let Some(root) = refined(ctx, f, x, &grid[i], &grid[i + 1], s, method, bits)
        {
            roots.push(root);
        }
    }
    roots
}

/// The root of `f` between `lo` and `hi`, where `f` has the sign `s_lo` at
/// `lo` and the other one at `hi`, by bisection to `bits` bits: a ball, when
/// `f` over the last interval is a finite ball. Also the monotone
/// inversion's (CAS-302).
pub(in crate::native) fn bisected(
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
            // Within rounding of 0 at `mid`, which may come long before the
            // interval is narrow: an exact root at a midpoint, as Student's
            // t(2) inverse at 8/9 is 7/4, the third midpoint of `[1, 3]`
            // (CAS-305). The root is within the resolution around `mid`
            // when the signs there tell; otherwise the interval stands,
            // as narrow as `bits + 32` tells.
            None => {
                let delta = resolution(&mid, bits);
                let (below, above) = (&mid - &delta, &mid + &delta);
                if sign(ctx, f, x, &below, bits + 32) == Some(s_lo)
                    && sign(ctx, f, x, &above, bits + 32) == Some(s_lo.reverse())
                {
                    (lo, hi) = (below, above);
                }
                break;
            }
        }
    }
    certified(ctx, f, x, &lo, &hi, bits)
}

/// Half the width [`narrow`] allows around `x`: `2^−(bits+1)` of `|x|`, or
/// of `2^−bits` about 0. Also the monotone inversion's (CAS-306).
pub(in crate::native) fn resolution(x: &Rational, bits: u64) -> Rational {
    let two = Rational::from(Integer::from(2));
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let floor = two.pow(-exponent);
    let scale = x.abs().max(floor);
    &scale * &two.pow(-exponent.saturating_add(1))
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

/// The root of `f` in `[lo, hi]`, where its sign is `s_lo` at `lo` and
/// the other at `hi`: the method's estimate when it is certified within
/// the bracket, else by bisection.
#[allow(clippy::too_many_arguments, reason = "a bracket and its method")]
fn refined(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    lo: &Rational,
    hi: &Rational,
    s_lo: Ordering,
    method: Method,
    bits: u64,
) -> Option<Expr> {
    let estimate = match method {
        Method::FalsePosition => false_position(ctx, f, x, lo, hi, bits),
        Method::Brent => brent(ctx, f, x, lo, hi, bits),
        _ => None,
    };
    if let Some(q) = estimate
        && let Some(root) = around(ctx, f, x, &q, bits)
        && within(&root, lo, hi)
    {
        return Some(root);
    }
    bisected(ctx, f, x, lo.clone(), hi.clone(), s_lo, bits)
}

/// `q` rounded to `prec` bits, so that iterates stay short.
fn rounded(q: &Rational, prec: u64) -> Rational {
    Float::from_rational(q, prec, Round::Nearest)
        .0
        .to_rational()
}

/// The midpoint of `f`'s ball at `q`, at `prec` bits.
fn at(ctx: &Context, f: &Expr, x: &Expr, q: &Rational, prec: u64) -> Option<Rational> {
    value(ctx, f, x, &Expr::from(q.clone()), prec).map(|b| b.mid().to_rational())
}

/// Whether `step` is below `bits` bits of `q`, or of 1 about 0.
fn settled(step: &Rational, q: &Rational, bits: u64) -> bool {
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX).saturating_add(8);
    let scale = q.abs().max(Rational::from(Integer::from(1)));
    step.is_zero() || step.abs() <= &scale * &Rational::from(Integer::from(2)).pow(-exponent)
}

/// Newton's, the secant's (from `q0` and `q1`) or Steffensen's iteration
/// on `f` from `q1`: its limit, uncertified.
#[allow(clippy::too_many_arguments, reason = "the iteration's state")]
fn iterated(
    ctx: &Context,
    f: &Expr,
    df: Option<&Expr>,
    x: &Expr,
    method: Method,
    mut q0: Rational,
    mut q: Rational,
    bits: u64,
) -> Option<Rational> {
    let prec = bits.saturating_add(32);
    let mut f0 = match method {
        Method::Secant => at(ctx, f, x, &q0, prec)?,
        _ => Rational::ZERO,
    };
    for _ in 0..MAX_NEWTON {
        let v = at(ctx, f, x, &q, prec)?;
        if v.is_zero() {
            return Some(q);
        }
        let step = match method {
            Method::Newton => {
                let d = at(ctx, df?, x, &q, prec)?;
                if d.is_zero() {
                    return None;
                }
                &v / &d
            }
            Method::Secant => {
                let den = &v - &f0;
                if den.is_zero() {
                    return None;
                }
                let step = &(&v * &(&q - &q0)) / &den;
                (q0, f0) = (q.clone(), v);
                step
            }
            Method::Steffensen => {
                let w = at(ctx, f, x, &(&q + &v), prec)?;
                let den = &w - &v;
                if den.is_zero() {
                    return None;
                }
                &(&v * &v) / &den
            }
            _ => return None,
        };
        let next = rounded(&(&q - &step), prec);
        let done = settled(&step, &next, bits);
        q = next;
        if done {
            return Some(q);
        }
    }
    None
}

/// The Illinois variant of false position on `[lo, hi]`: its estimate.
pub(super) fn false_position(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    lo: &Rational,
    hi: &Rational,
    bits: u64,
) -> Option<Rational> {
    let prec = bits.saturating_add(32);
    let two = Rational::from(Integer::from(2));
    let (mut a, mut b) = (lo.clone(), hi.clone());
    let (mut fa, mut fb) = (at(ctx, f, x, &a, prec)?, at(ctx, f, x, &b, prec)?);
    let mut side = Ordering::Equal;
    let mut last: Option<Rational> = None;
    for _ in 0..MAX_BISECTION {
        let den = &fb - &fa;
        if den.is_zero() {
            return None;
        }
        let c = rounded(&(&(&(&a * &fb) - &(&b * &fa)) / &den), prec);
        if let Some(l) = &last
            && settled(&(&c - l), &c, bits)
        {
            return Some(c);
        }
        let fc = at(ctx, f, x, &c, prec)?;
        if fc.is_zero() || narrow(&a, &b, bits) {
            return Some(c);
        }
        if (fc > Rational::ZERO) == (fb > Rational::ZERO) {
            (b, fb) = (c.clone(), fc);
            if side == Ordering::Less {
                fa = &fa / &two;
            }
            side = Ordering::Less;
        } else {
            (a, fa) = (c.clone(), fc);
            if side == Ordering::Greater {
                fb = &fb / &two;
            }
            side = Ordering::Greater;
        }
        last = Some(c);
    }
    last
}

/// Brent's method on `[lo, hi]` (inverse quadratic interpolation, secant
/// and bisection steps): its estimate.
pub(super) fn brent(
    ctx: &Context,
    f: &Expr,
    x: &Expr,
    lo: &Rational,
    hi: &Rational,
    bits: u64,
) -> Option<Rational> {
    let prec = bits.saturating_add(32);
    let zero = Rational::ZERO;
    let one = Rational::from(Integer::from(1));
    let two = Rational::from(Integer::from(2));
    let three = Rational::from(Integer::from(3));
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let eps = two.pow(-exponent.saturating_add(2));
    let (mut a, mut b) = (lo.clone(), hi.clone());
    let (mut fa, mut fb) = (at(ctx, f, x, &a, prec)?, at(ctx, f, x, &b, prec)?);
    let (mut c, mut fc) = (b.clone(), fb.clone());
    let mut d = &b - &a;
    let mut e = d.clone();
    for _ in 0..MAX_BISECTION {
        if (fb > zero) == (fc > zero) {
            (c, fc) = (a.clone(), fa.clone());
            d = &b - &a;
            e = d.clone();
        }
        if fc.abs() < fb.abs() {
            (a, b, c) = (b.clone(), c.clone(), b.clone());
            (fa, fb, fc) = (fb.clone(), fc.clone(), fb.clone());
        }
        let tol = &(&two * &eps) * &b.abs().max(eps.clone());
        let xm = &(&c - &b) / &two;
        if xm.abs() <= tol || fb.is_zero() {
            return Some(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = &fb / &fa;
            let (mut p, mut q) = if a == c {
                (&(&two * &xm) * &s, &one - &s)
            } else {
                let q = &fa / &fc;
                let r = &fb / &fc;
                let p =
                    &s * &(&(&(&(&two * &xm) * &q) * &(&q - &r)) - &(&(&b - &a) * &(&r - &one)));
                (p, &(&(&q - &one) * &(&r - &one)) * &(&s - &one))
            };
            if p > zero {
                q = -q;
            }
            p = p.abs();
            let min1 = &(&(&three * &xm) * &q) - &(&tol * &q).abs();
            let min2 = (&e * &q).abs();
            if &two * &p < min1.min(min2) {
                e = d;
                d = &p / &q;
            } else {
                d = xm.clone();
                e = d.clone();
            }
        } else {
            d = xm.clone();
            e = d.clone();
        }
        (a, fa) = (b.clone(), fb.clone());
        b = if d.abs() > tol {
            rounded(&(&b + &d), prec)
        } else if xm > zero {
            &b + &tol
        } else {
            &b - &tol
        };
        fb = at(ctx, f, x, &b, prec)?;
    }
    Some(b)
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
