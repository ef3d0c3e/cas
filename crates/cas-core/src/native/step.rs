//! The Heaviside step and the Dirac delta (CAS-DOC-18 §2, CAS-273), Xcas's
//! `Heaviside` and `Dirac`.
//!
//! - `HeavisideTheta[x]`: 0 for `x < 0`, 1 for `x > 0`, `1/2` at 0 as
//!   giac's, the sign from the assumptions, a ball or the approximation, as
//!   `Sign` takes it; otherwise it stays.
//! - `DiracDelta[x]`: 0 for `x ≠ 0`; `ComplexInfinity` at 0, with the
//!   warning of `1/0` (giac's `infinity`); otherwise it stays.
//!
//! Their calculus is elsewhere: `D` of the step is the delta (`calculus`),
//! and `Integrate` has both at linear arguments and sifts the delta out of
//! definite integrals (`integrate`).
//!
//! Two rewritings in a variable `x` (giac's `x` when omitted), ports of
//! giac's signalprocessing.cc:
//!
//! - [`linear_step`]: `LinearStep[f, x]`, giac's `linstep` (CAS-570), `f`
//!   as a sum of steps `θ(x − c)` times coefficients. Each `Piecewise` of
//!   comparisons becomes steps, value `j` weighted by `θ(d_j) Π_{i<j}
//!   θ(−d_i)` (`flatten_piecewise`); `Sign[u]` is `2 θ(u) − 1`; a function
//!   of `Abs` of linear arguments with numeric zeros is split at the sorted
//!   zeros, each `Abs`, step and sign resolved on each interval and the
//!   pieces weighted by boxes, again while nested ones remain
//!   (`abs2Heaviside`, the arguments of `Abs` factored first); a step of a
//!   quadratic with two real roots is two steps (`simplify_quadratic_
//!   Heaviside`); steps of `a x + b` become `1/2 + sign(a) (θ(x + b/a) −
//!   1/2)`, the deltas `δ(x + b/a)/|a|`; the result is expanded, a product
//!   of steps `θ(x − c₁) θ(x − c₂)` is `θ(x − max(c₁, c₂))`, a positive power
//!   of a step the step (`lin_Heaviside`), and the terms are collected by
//!   step, each coefficient brought `Together` (`simpcollect`). Bounds for
//!   the fuzzer, giac has none: at most [`MAX_BREAKS`] zeros of `Abs`
//!   arguments, [`MAX_PIECES`] intervals over all rounds. Departure: for a
//!   step of a quadratic of negative leading coefficient giac writes the
//!   box with the opposite sign, ours the box itself.
//! - [`simplify_dirac`]: `SimplifyDirac[f, x]`, giac's `simplifyDirac`
//!   (CAS-572). The derivative `δ⁽ⁿ⁾(u)` is `DiracDelta[u, n]` (giac's
//!   `Dirac(u, n)`; `Derivative[n][DiracDelta][u]` is read as it). A delta
//!   of a rational argument whose denominator has no real zero splits into
//!   `Σ δ(x − r)/|u′(r)|` over the real zeros of the numerator, `0` with
//!   none, unchanged at a multiple zero (`expand_Dirac`, numerators of
//!   degree at most [`MAX_SPLIT_DEGREE`]). Then each term `g(x) δ⁽ⁿ⁾(a x +
//!   b)` of the expanded expression, `a`, `b` real, sifts at `z = −b/a`:
//!   `δ⁽ⁿ⁾(a x + b) = δ⁽ⁿ⁾(x − z)/(|a| aⁿ)` and `g(x) δ⁽ⁿ⁾(x − z) = Σ_k
//!   (−1)ᵏ C(n, k) g⁽ᵏ⁾(z) δ⁽ⁿ⁻ᵏ⁾(x − z)`; for `g = c/(x − z)ᵏ` giac's
//!   rule gives `c/((−1)ᵏ C(n + k, k) k!) δ⁽ⁿ⁺ᵏ⁾(x − z)`, another pole leaves
//!   the term (`simplify_Dirac`, `simp_dirac_rat`). The coefficients of
//!   each delta and the rest go through `Simplify`. Orders above
//!   [`MAX_DIRAC_ORDER`] stay. Departures: giac's coefficients from a split
//!   delta print as `series(c, x, z, 1)`, ours as `c`; its assumptions on
//!   `x` are not read, `x` is real.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: f in x, u = a x + b, its zero z, δ⁽ⁿ⁾ and r(x) its coefficient"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::piecewise::branches;
use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `HeavisideTheta[x]`.
pub(super) fn heaviside(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(real_sign(ev.context(), x).map(|s| match s {
        Ordering::Greater => Expr::from(1),
        Ordering::Less => Expr::from(0),
        Ordering::Equal => Expr::from(Rational::new(Integer::ONE, Integer::from(2))),
    }))
}

/// `UnitBox[x]`, the box: 1 for `|x| < 1/2`, 0 for `|x| > 1/2`, `1/2` at
/// `±1/2`, the mean of its sides as `HeavisideTheta[0]` is (CAS-DOC-41
/// §1); otherwise it stays.
pub(super) fn unit_box(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let half = Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
    let below = ev.eval(&(x.clone() + half.clone()))?;
    let above = ev.eval(&(x.clone() - half.clone()))?;
    let ctx = ev.context();
    Ok(match (real_sign(ctx, &below), real_sign(ctx, &above)) {
        (Some(Ordering::Less), _) | (_, Some(Ordering::Greater)) => Some(Expr::from(0)),
        (Some(Ordering::Greater), Some(Ordering::Less)) => Some(Expr::from(1)),
        (Some(Ordering::Equal), _) | (_, Some(Ordering::Equal)) => Some(half),
        _ => None,
    })
}

/// `DiracDelta[x]`.
pub(super) fn dirac(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    match real_sign(ev.context(), x) {
        Some(Ordering::Equal) => {
            ev.warn(Warning::Infinite(e.clone()));
            Ok(Some(Symbol::COMPLEX_INFINITY.into()))
        }
        Some(_) => Ok(Some(Expr::from(0))),
        None => Ok(None),
    }
}

/// At most this many zeros of `Abs` arguments in one split (fuzz bound).
const MAX_BREAKS: usize = 16;
/// At most this many intervals over all rounds of the `Abs` split.
const MAX_PIECES: usize = 64;
/// Delta derivatives of higher order stay.
const MAX_DIRAC_ORDER: i64 = 32;
/// A delta of a numerator of higher degree is not split at its zeros.
const MAX_SPLIT_DEGREE: i64 = 8;
/// Rounds of the loops that rewrite until nothing changes.
const MAX_ROUNDS: usize = 64;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn theta(u: Expr) -> Expr {
    Expr::call(Symbol::HEAVISIDE_THETA, [u])
}

/// `δ⁽ᵏ⁾(u)`: `DiracDelta[u]`, or `DiracDelta[u, k]` for `k > 0`.
fn delta(u: Expr, k: i64) -> Expr {
    if k == 0 {
        Expr::call(Symbol::DIRAC_DELTA, [u])
    } else {
        Expr::call(Symbol::DIRAC_DELTA, [u, Expr::from(k)])
    }
}

/// The operands of a sum, or the expression itself.
fn terms(e: &Expr) -> Vec<Expr> {
    if head(e) == Some(Symbol::PLUS) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    }
}

/// The factors of a product, or the expression itself.
fn factors(e: &Expr) -> Vec<Expr> {
    if head(e) == Some(Symbol::TIMES) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    }
}

/// The expression and the variable of `F[f, x]`, `x` when omitted.
pub(super) fn with_variable(args: &[Expr]) -> Option<(&Expr, Expr)> {
    match args {
        [f] => Some((f, Expr::symbol("x"))),
        [f, x] if matches!(x.view(), View::Symbol(s) if !s.is_builtin()) => Some((f, x.clone())),
        _ => None,
    }
}

/// `g` with `f` applied top-down as [`Expr::replace`] does, the first error
/// kept.
fn try_replace(
    ev: &mut Eval<'_>,
    g: &Expr,
    mut f: impl FnMut(&mut Eval<'_>, &Expr) -> Result<Option<Expr>, EvalError>,
) -> Result<Expr, EvalError> {
    let mut error = None;
    let r = g.replace(|e| {
        if error.is_some() {
            return None;
        }
        f(&mut *ev, e).unwrap_or_else(|err| {
            error = Some(err);
            None
        })
    });
    error.map_or(Ok(r), Err)
}

/// `(a, b)` with `u = a x + b`, both free of `x`.
pub(super) fn linear(
    ev: &mut Eval<'_>,
    u: &Expr,
    x: &Expr,
) -> Result<Option<(Expr, Expr)>, EvalError> {
    let a = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
    if a.contains(x) {
        return Ok(None);
    }
    let b = ev.eval(&u.substitute(&[(x.clone(), Expr::from(0))]))?;
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - a.clone() * x.clone() - b.clone()],
    ))?;
    Ok((gap == Expr::from(0)).then_some((a, b)))
}

/// `(a, b)` with `u = a x + b`, `a ≠ 0`.
fn proper_linear(ev: &mut Eval<'_>, u: &Expr, x: &Expr) -> Result<Option<(Expr, Expr)>, EvalError> {
    if !u.contains(x) {
        return Ok(None);
    }
    Ok(linear(ev, u, x)?.filter(|(a, _)| *a != Expr::from(0)))
}

/// The sign of `a − b`, when decided.
fn order(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<Ordering>, EvalError> {
    let d = ev.eval(&(a.clone() - b.clone()))?;
    Ok(real_sign(ev.context(), &d))
}

/// Whether `a ≥ b` is decided true.
fn at_least(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    Ok(matches!(
        order(ev, a, b)?,
        Some(Ordering::Greater | Ordering::Equal)
    ))
}

/// `LinearStep[f, x]`, giac's `linstep` (signalprocessing.cc `_linstep`).
pub(super) fn linear_step(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    if head(f) == Some(Symbol::LIST) {
        let items = f
            .args()
            .iter()
            .map(|g| Expr::call(Symbol::LINEAR_STEP, [g.clone(), x.clone()]));
        return ev.eval(&Expr::call(Symbol::LIST, items)).map(Some);
    }
    let g = steps(ev, f, &x)?;
    collect_steps(ev, &g).map(Some)
}

/// `f` as a sum of steps `θ(x − c)` times coefficients, before they are
/// collected: giac's `linH(simplify_quadratic_Heaviside(abs2Heaviside(
/// sign2Heaviside(flatten_piecewise(f)))))`, shared with `linabs`.
pub(super) fn steps(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = ev.eval(&flatten_piecewise(f))?;
    let g = ev.eval(&Expr::call(Symbol::SIGN_TO_HEAVISIDE, [g]))?;
    let g = abs_to_heaviside(ev, &g, x)?;
    let g = quadratic_steps(ev, &g, x)?;
    let g = linear_arguments(ev, &g, x)?;
    steps_of_products(ev, &g, x)
}

/// giac's `flatten_piecewise`: each `Piecewise` of comparisons as a sum of
/// steps.
fn flatten_piecewise(f: &Expr) -> Expr {
    f.replace(|p| {
        if head(p) == Some(Symbol::PIECEWISE) {
            flatten_one(p)
        } else {
            None
        }
    })
}

fn flatten_one(p: &Expr) -> Option<Expr> {
    let (pairs, default) = branches(p.args());
    let mut last = Expr::from(1);
    let mut sum = Vec::with_capacity(pairs.len() + 1);
    for (c, v) in pairs {
        let d = match (head(c), c.args()) {
            (Some(Symbol::GREATER | Symbol::GREATER_EQUAL), [u, w]) => u.clone() - w.clone(),
            (Some(Symbol::LESS | Symbol::LESS_EQUAL), [u, w]) => w.clone() - u.clone(),
            _ => return None,
        };
        sum.push(last.clone() * theta(d.clone()) * flatten_piecewise(v));
        last *= theta(-d);
    }
    if let Some(d) = default {
        sum.push(last * flatten_piecewise(d));
    }
    Some(sum.into_iter().sum())
}

/// giac's `abs2Heaviside`: `g` split at the zeros of its `Abs` arguments.
fn abs_to_heaviside(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let mut g = g.clone();
    let mut pieces = 1;
    loop {
        let factored = factor_abs(ev, &g, x)?;
        let mut zeros = Vec::new();
        if !abs_zeros(ev, &factored, x, &mut zeros)? || zeros.is_empty() {
            return Ok(g);
        }
        pieces *= zeros.len() + 1;
        if pieces > MAX_PIECES {
            return Ok(g);
        }
        let Some(zeros) = sorted_numbers(ev, zeros)? else {
            return Ok(g);
        };
        let mut sum = Vec::with_capacity(zeros.len() + 1);
        let mut lo: Option<Expr> = None;
        for hi in zeros.into_iter().map(Some).chain([None]) {
            let piece = on_interval(ev, &factored, x, lo.as_ref(), hi.as_ref())?;
            let window = match (&lo, &hi) {
                (Some(a), Some(b)) => theta(x.clone() - a.clone()) - theta(x.clone() - b.clone()),
                (Some(a), None) => theta(x.clone() - a.clone()),
                (None, Some(b)) => Expr::from(1) - theta(x.clone() - b.clone()),
                (None, None) => Expr::from(1),
            };
            sum.push(piece * window);
            lo = hi;
        }
        g = ev.eval(&sum.into_iter().sum())?;
    }
}

/// giac's `factor_abs_arg`: the argument of each `Abs` in `x` factored,
/// `Abs` of a product the product of the `Abs` of its factors.
fn factor_abs(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = try_replace(ev, g, |ev, p| {
        let (Some(Symbol::ABS), [u]) = (head(p), p.args()) else {
            return Ok(None);
        };
        if !u.contains(x) {
            return Ok(None);
        }
        let inner = factor_abs(ev, u, x)?;
        let f = ev.eval(&Expr::call(Symbol::FACTOR, [inner]))?;
        Ok(Some(split_abs(&f)))
    })?;
    ev.eval(&g)
}

fn split_abs(f: &Expr) -> Expr {
    match (head(f), f.args()) {
        (Some(Symbol::TIMES), fs) => fs.iter().map(split_abs).product(),
        (Some(Symbol::POWER), [b, n]) if n.to_integer().is_some() => split_abs(b).pow(n.clone()),
        _ => Expr::call(Symbol::ABS, [f.clone()]),
    }
}

/// giac's `find_func_linarg` for `Abs`, nested: the zeros of the linear
/// arguments in `x`, each once; `false` past [`MAX_BREAKS`].
fn abs_zeros(
    ev: &mut Eval<'_>,
    g: &Expr,
    x: &Expr,
    out: &mut Vec<Expr>,
) -> Result<bool, EvalError> {
    if !g.contains(x) {
        return Ok(true);
    }
    if let (Some(Symbol::ABS), [u]) = (head(g), g.args()) {
        if let Some((a, b)) = proper_linear(ev, u, x)? {
            let z = ev.eval(&(-b / a))?;
            if !out.contains(&z) {
                out.push(z);
            }
            return Ok(out.len() <= MAX_BREAKS);
        }
        return abs_zeros(ev, u, x, out);
    }
    for a in g.args() {
        if !abs_zeros(ev, a, x, out)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The numbers `zs` in increasing order, `None` when one is not a number
/// or two are not decided apart.
fn sorted_numbers(ev: &mut Eval<'_>, zs: Vec<Expr>) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut sorted: Vec<Expr> = Vec::with_capacity(zs.len());
    for z in zs {
        if !z.is_numeric() {
            return Ok(None);
        }
        let mut at = sorted.len();
        for (i, s) in sorted.iter().enumerate() {
            match order(ev, &z, s)? {
                Some(Ordering::Less) => {
                    at = i;
                    break;
                }
                Some(Ordering::Greater) => {}
                _ => return Ok(None),
            }
        }
        sorted.insert(at, z);
    }
    Ok(Some(sorted))
}

/// giac's `eval_for_x_in_ab`: `g` for `x` between `lo` and `hi` (unbounded
/// when `None`), each `Abs`, step and sign of a linear argument whose zero
/// lies outside resolved, until none is left.
fn on_interval(
    ev: &mut Eval<'_>,
    g: &Expr,
    x: &Expr,
    lo: Option<&Expr>,
    hi: Option<&Expr>,
) -> Result<Expr, EvalError> {
    let mut g = g.clone();
    for _ in 0..MAX_ROUNDS {
        let mut rules = Vec::new();
        for part in g.distinct() {
            let (Some(h @ (Symbol::ABS | Symbol::HEAVISIDE_THETA | Symbol::SIGN)), [u]) =
                (head(part), part.args())
            else {
                continue;
            };
            let Some((c, d)) = proper_linear(ev, u, x)? else {
                continue;
            };
            let s = match real_sign(ev.context(), &c) {
                Some(Ordering::Greater) => 1,
                Some(Ordering::Less) => -1,
                _ => continue,
            };
            let z = ev.eval(&(-d / c))?;
            let right = match lo {
                Some(a) => at_least(ev, a, &z)?,
                None => false,
            };
            let left = !right
                && match hi {
                    Some(b) => at_least(ev, &z, b)?,
                    None => false,
                };
            // On the right of the zero the function is its value for a
            // positive argument times the sign of `c`, on the left the
            // other.
            let side = match (right, left) {
                (true, _) => s,
                (_, true) => -s,
                _ => continue,
            };
            let value = match h {
                Symbol::HEAVISIDE_THETA => Expr::from(i64::from(side > 0)),
                Symbol::ABS => Expr::from(side) * u.clone(),
                _ => Expr::from(side),
            };
            rules.push((part.clone(), value));
        }
        if rules.is_empty() {
            break;
        }
        g = ev.eval(&g.substitute(&rules))?;
    }
    Ok(g)
}

/// `(a, b, c)` with `u = a x² + b x + c`, all free of `x`.
fn quadratic(
    ev: &mut Eval<'_>,
    u: &Expr,
    x: &Expr,
) -> Result<Option<(Expr, Expr, Expr)>, EvalError> {
    let d = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
    let Some((a2, b)) = linear(ev, &d, x)? else {
        return Ok(None);
    };
    let a = ev.eval(&(a2 / 2))?;
    let c = ev.eval(&u.substitute(&[(x.clone(), Expr::from(0))]))?;
    let gap = ev.eval(&Expr::call(
        Symbol::EXPAND,
        [u.clone() - a.clone() * x.pow(2) - b.clone() * x.clone() - c.clone()],
    ))?;
    Ok((gap == Expr::from(0)).then_some((a, b, c)))
}

/// giac's `simplify_quadratic_Heaviside`: a step of a quadratic in `x` with
/// two real roots `r₁ < r₂` is `θ(r₁ − x) + θ(x − r₂)` or, opening
/// downwards, the box `θ(x − r₁) − θ(x − r₂)`; with none it is constant.
fn quadratic_steps(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = try_replace(ev, g, |ev, p| {
        let (Some(Symbol::HEAVISIDE_THETA), [u]) = (head(p), p.args()) else {
            return Ok(None);
        };
        if !u.contains(x) {
            return Ok(None);
        }
        let Some((a, b, c)) = quadratic(ev, u, x)? else {
            return Ok(None);
        };
        if a == Expr::from(0) {
            return Ok(None);
        }
        let d = ev.eval(&(b.clone().pow(2) - 4 * a.clone() * c))?;
        let up = match real_sign(ev.context(), &a) {
            Some(Ordering::Greater) => true,
            Some(Ordering::Less) => false,
            _ => return Ok(None),
        };
        match real_sign(ev.context(), &d) {
            Some(Ordering::Less | Ordering::Equal) => Ok(Some(theta(a))),
            Some(Ordering::Greater) => {
                let root = |s: i64| {
                    let half = Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
                    (-b.clone() + s * d.clone().pow(half)) / (2 * a.clone())
                };
                // With `a < 0` the root with `−√d` is the larger.
                let (r1, r2) = if up {
                    (root(-1), root(1))
                } else {
                    (root(1), root(-1))
                };
                Ok(Some(if up {
                    theta(r1 - x.clone()) + theta(x.clone() - r2)
                } else {
                    theta(x.clone() - r1) - theta(x.clone() - r2)
                }))
            }
            None => Ok(None),
        }
    })?;
    ev.eval(&g)
}

/// giac's `simplify_linarg`: steps and signs of `a x + b` in terms of
/// `x + b/a`, deltas `δ(x + b/a)/|a|`, `δ⁽ⁿ⁾(−u)` as `(−1)ⁿ δ⁽ⁿ⁾(u)`.
fn linear_arguments(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = try_replace(ev, g, |ev, p| {
        let h = head(p);
        match (h, p.args()) {
            (Some(Symbol::HEAVISIDE_THETA | Symbol::SIGN | Symbol::DIRAC_DELTA), [u]) => {
                let Some((a, b)) = proper_linear(ev, u, x)? else {
                    return Ok(None);
                };
                let v = ev.eval(&(x.clone() + b / a.clone()))?;
                let sign = Expr::call(Symbol::SIGN, [a.clone()]);
                let half = Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
                Ok(Some(match h {
                    Some(Symbol::HEAVISIDE_THETA) => half.clone() + sign * (theta(v) - half),
                    Some(Symbol::SIGN) => sign * Expr::call(Symbol::SIGN, [v]),
                    _ => delta(v, 0) / Expr::call(Symbol::ABS, [a]),
                }))
            }
            (Some(Symbol::DIRAC_DELTA), [u, n]) => {
                let Some(n) = n.as_i64() else {
                    return Ok(None);
                };
                match (head(u), u.args()) {
                    (Some(Symbol::TIMES), [m, ..]) if *m == Expr::from(-1) => {
                        let w = ev.eval(&-u.clone())?;
                        Ok(Some(
                            Expr::from(if n % 2 == 0 { 1 } else { -1 }) * delta(w, n),
                        ))
                    }
                    _ => Ok(None),
                }
            }
            _ => Ok(None),
        }
    })?;
    ev.eval(&g)
}

/// giac's `lin_Heaviside` on the expanded `g`: in each term, a positive
/// power of a step is the step and two steps `θ(x − c₁) θ(x − c₂)` are
/// `θ(x − max(c₁, c₂))`, until no product of steps is left.
fn steps_of_products(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = ev.eval(&Expr::call(Symbol::EXPAND, [g.clone()]))?;
    let mut sum = Vec::new();
    for t in terms(&g) {
        let mut t = t;
        for _ in 0..MAX_ROUNDS {
            let mut rest = Vec::new();
            let mut shifts: Vec<Expr> = Vec::new();
            let mut again = false;
            for f in factors(&t) {
                if let (Some(Symbol::POWER), [b, p]) = (head(&f), f.args())
                    && head(b) == Some(Symbol::HEAVISIDE_THETA)
                    && real_sign(ev.context(), p) == Some(Ordering::Greater)
                {
                    rest.push(b.clone());
                    again = true;
                    continue;
                }
                if shifts.len() < 2
                    && let (Some(Symbol::HEAVISIDE_THETA), [u]) = (head(&f), f.args())
                    && let Some((a, b)) = proper_linear(ev, u, x)?
                    && a == Expr::from(1)
                {
                    shifts.push(-b);
                    continue;
                }
                rest.push(f);
            }
            match &shifts[..] {
                [c1, c2] => {
                    let max = Expr::call(Symbol::MAX, [c1.clone(), c2.clone()]);
                    rest.push(theta(x.clone() - max));
                    again = true;
                }
                [c] => rest.push(theta(x.clone() - c.clone())),
                _ => {}
            }
            t = ev.eval(&rest.into_iter().product())?;
            if !again {
                break;
            }
        }
        sum.push(t);
    }
    ev.eval(&sum.into_iter().sum())
}

/// giac's `simpcollect` for steps: the terms with one step collected by
/// step, each coefficient and the rest brought `Together`.
fn collect_steps(ev: &mut Eval<'_>, g: &Expr) -> Result<Expr, EvalError> {
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let mut rest = Vec::new();
    for t in terms(g) {
        let fs = factors(&t);
        let steps: Vec<usize> = (0..fs.len())
            .filter(|&i| head(&fs[i]) == Some(Symbol::HEAVISIDE_THETA))
            .collect();
        let [i] = steps[..] else {
            rest.push(t);
            continue;
        };
        let coefficient: Expr = fs
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, f)| f.clone())
            .product();
        match groups.iter_mut().find(|(s, _)| *s == fs[i]) {
            Some((_, cs)) => cs.push(coefficient),
            None => groups.push((fs[i].clone(), vec![coefficient])),
        }
    }
    let together = |cs: Vec<Expr>| Expr::call(Symbol::TOGETHER, [cs.into_iter().sum::<Expr>()]);
    let mut sum: Vec<Expr> = groups.into_iter().map(|(s, cs)| s * together(cs)).collect();
    sum.push(together(rest));
    ev.eval(&sum.into_iter().sum())
}

/// `SimplifyDirac[f, x]`, giac's `simplifyDirac` (signalprocessing.cc
/// `_simplifyDirac`).
pub(super) fn simplify_dirac(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let g = f.replace(delta_derivative);
    if !g.contains(&Expr::from(Symbol::DIRAC_DELTA)) {
        return Ok(Some(g));
    }
    let g = split_deltas(ev, &g, &x)?;
    sift(ev, &g, &x).map(Some)
}

/// `Derivative[n][DiracDelta][u]` as `DiracDelta[u, n]`.
fn delta_derivative(p: &Expr) -> Option<Expr> {
    let [u] = p.args() else {
        return None;
    };
    let h = p.head()?;
    let [f] = h.args() else {
        return None;
    };
    if f.as_symbol() != Some(Symbol::DIRAC_DELTA) {
        return None;
    }
    let d = h.head()?;
    let (Some(Symbol::DERIVATIVE), [n]) = (head(d), d.args()) else {
        return None;
    };
    Some(Expr::call(Symbol::DIRAC_DELTA, [u.clone(), n.clone()]))
}

/// Whether `e` is a rational function of `x`.
fn rational_in(e: &Expr, x: &Expr) -> bool {
    if e == x || !e.contains(x) {
        return true;
    }
    match (head(e), e.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), args) => args.iter().all(|a| rational_in(a, x)),
        (Some(Symbol::POWER), [b, n]) => n.to_integer().is_some() && rational_in(b, x),
        _ => false,
    }
}

/// The real zeros of the polynomial `p` in `x` by `Solve`, `None` when it
/// does not give them or the degree is out of `1..=`[`MAX_SPLIT_DEGREE`].
fn zeros(ev: &mut Eval<'_>, p: &Expr, x: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [p.clone(), x.clone()]))?;
    if !degree
        .as_i64()
        .is_some_and(|d| (1..=MAX_SPLIT_DEGREE).contains(&d))
    {
        return Ok(None);
    }
    let equation = Expr::call(Symbol::EQUATION, [p.clone(), Expr::from(0)]);
    let found = ev.eval(&Expr::call(Symbol::SOLVE, [equation, x.clone()]))?;
    let fine = head(&found) == Some(Symbol::LIST)
        && found.args().iter().all(|z| {
            !z.contains(x)
                && !matches!(
                    head(z),
                    Some(Symbol::LIST | Symbol::EQUATION | Symbol::RULE)
                )
        });
    Ok(fine.then(|| found.args().to_vec()))
}

/// giac's `expand_Dirac`: `δ(u)` for `u` rational in `x` with no real pole
/// as `Σ δ(x − r)/|u′(r)|` over the real zeros `r` of `u`.
fn split_deltas(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = try_replace(ev, g, |ev, p| {
        let (Some(Symbol::DIRAC_DELTA), [u]) = (head(p), p.args()) else {
            return Ok(None);
        };
        if !u.contains(x) || !rational_in(u, x) {
            return Ok(None);
        }
        let r = ev.eval(&Expr::call(Symbol::TOGETHER, [u.clone()]))?;
        let den = ev.eval(&Expr::call(Symbol::DENOMINATOR, [r.clone()]))?;
        if den.contains(x) && zeros(ev, &den, x)?.is_none_or(|z| !z.is_empty()) {
            return Ok(None);
        }
        let num = ev.eval(&Expr::call(
            Symbol::EXPAND,
            [Expr::call(Symbol::NUMERATOR, [r])],
        ))?;
        let Some(zs) = zeros(ev, &num, x)? else {
            return Ok(None);
        };
        let du = ev.eval(&Expr::call(Symbol::D, [u.clone(), x.clone()]))?;
        let mut sum = Vec::with_capacity(zs.len());
        for z in zs {
            let at = du.substitute(&[(x.clone(), z.clone())]);
            let slope = ev.eval(&Expr::call(Symbol::TOGETHER, [at]))?;
            if slope == Expr::from(0) {
                return Ok(None);
            }
            sum.push(delta(x.clone() - z, 0) / Expr::call(Symbol::ABS, [slope]));
        }
        Ok(Some(sum.into_iter().sum()))
    })?;
    ev.eval(&g)
}

/// `(z, n, r)` for a term `r(x) δ⁽ⁿ⁾(x − z)` (giac's `has_Dirac`): the
/// first delta factor of a real linear argument `a x + b`, `z = −b/a`, the
/// other factors over `|a| aⁿ`; `None` when there is none or another delta.
fn delta_term(
    ev: &mut Eval<'_>,
    t: &Expr,
    x: &Expr,
) -> Result<Option<(Expr, i64, Expr)>, EvalError> {
    let fs = factors(t);
    for (i, f) in fs.iter().enumerate() {
        let (u, n) = match (head(f), f.args()) {
            (Some(Symbol::DIRAC_DELTA), [u]) => (u, Some(0)),
            (Some(Symbol::DIRAC_DELTA), [u, n]) => (u, n.as_i64()),
            _ => continue,
        };
        let Some((a, b)) = proper_linear(ev, u, x)? else {
            continue;
        };
        let real = |ev: &mut Eval<'_>, c: &Expr| -> Result<bool, EvalError> {
            Ok(ev.eval(&Expr::call(Symbol::IM, [c.clone()]))? == Expr::from(0))
        };
        if !real(ev, &a)? || !real(ev, &b)? {
            continue;
        }
        let Some(n) = n.filter(|n| (0..=MAX_DIRAC_ORDER).contains(n)) else {
            return Ok(None);
        };
        let rest: Expr = fs
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, f)| f.clone())
            .product();
        if rest.contains(&Expr::from(Symbol::DIRAC_DELTA)) {
            return Ok(None);
        }
        let z = ev.eval(&(-b / a.clone()))?;
        let scale = Expr::call(Symbol::ABS, [a.clone()]) * a.pow(n);
        let r = ev.eval(&(rest / scale))?;
        return Ok(Some((z, n, r)));
    }
    Ok(None)
}

/// `C(n, k)` for small `n`.
fn binomial(n: i64, k: i64) -> i64 {
    (0..k).fold(1, |c, i| c * (n - i) / (i + 1))
}

/// `r(x) δ⁽ⁿ⁾(x − z)` as `Σ_k (−1)ᵏ C(n, k) r⁽ᵏ⁾(z) δ⁽ⁿ⁻ᵏ⁾(x − z)`, the
/// pairs `(n − k, coefficient)`; for `r = c/(x − z)ᵏ` giac's
/// `simp_dirac_rat`; `None` at another pole.
fn sift_one(
    ev: &mut Eval<'_>,
    r: &Expr,
    x: &Expr,
    z: &Expr,
    n: i64,
) -> Result<Option<Vec<(i64, Expr)>>, EvalError> {
    let at_z = |e: &Expr| e.substitute(&[(x.clone(), z.clone())]);
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [r.clone()]))?;
    let den = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together.clone()]))?;
    if ev.eval(&at_z(&den))? == Expr::from(0) {
        return pole(ev, &together, x, z, n);
    }
    let mut out = Vec::new();
    let mut dr = together;
    for k in 0..=n {
        let value = ev.eval(&at_z(&dr))?;
        let sign = if k % 2 == 0 { 1 } else { -1 };
        out.push((n - k, Expr::from(sign * binomial(n, k)) * value));
        if k < n {
            dr = ev.eval(&Expr::call(Symbol::D, [dr, x.clone()]))?;
        }
    }
    Ok(Some(out))
}

/// `c/(x − z)ᵏ δ⁽ⁿ⁾(x − z)` as `c/((−1)ᵏ C(n + k, k) k!) δ⁽ⁿ⁺ᵏ⁾(x − z)`
/// (giac's `simp_dirac_rat`), `k` at most 8.
fn pole(
    ev: &mut Eval<'_>,
    r: &Expr,
    x: &Expr,
    z: &Expr,
    n: i64,
) -> Result<Option<Vec<(i64, Expr)>>, EvalError> {
    let mut factorial = 1;
    for k in 1..=8 {
        factorial *= k;
        if n + k > MAX_DIRAC_ORDER {
            return Ok(None);
        }
        let shifted = r.clone() * (x.clone() - z.clone()).pow(k);
        let c = ev.eval(&Expr::call(Symbol::TOGETHER, [shifted]))?;
        if !c.contains(x) {
            let sign = if k % 2 == 0 { 1 } else { -1 };
            let scale = Expr::from(sign * binomial(n + k, k) * factorial);
            return Ok(Some(vec![(n + k, ev.eval(&(c / scale))?)]));
        }
    }
    Ok(None)
}

/// `c` added to the coefficients of `key`.
fn add(groups: &mut Vec<((Expr, i64), Vec<Expr>)>, key: (Expr, i64), c: Expr) {
    match groups.iter_mut().find(|(k, _)| *k == key) {
        Some((_, cs)) => cs.push(c),
        None => groups.push((key, vec![c])),
    }
}

/// giac's `simplify_Dirac`: the expanded `g`'s terms with a delta sifted,
/// collected by delta, each coefficient and the rest simplified.
fn sift(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let g = ev.eval(&Expr::call(Symbol::EXPAND, [g.clone()]))?;
    let mut groups: Vec<((Expr, i64), Vec<Expr>)> = Vec::new();
    let mut rest = Vec::new();
    let mut found: Vec<((Expr, i64), Vec<Expr>)> = Vec::new();
    for t in terms(&g) {
        match delta_term(ev, &t, x)? {
            Some((z, n, r)) => add(&mut found, (z, n), r),
            None => rest.push(t),
        }
    }
    for ((z, n), rs) in found {
        let r: Expr = rs.into_iter().sum();
        match sift_one(ev, &r, x, &z, n)? {
            Some(cs) => {
                for (k, c) in cs {
                    add(&mut groups, (z.clone(), k), c);
                }
            }
            None => add(&mut groups, (z, n), r),
        }
    }
    let simplify = |cs: Vec<Expr>| Expr::call(Symbol::SIMPLIFY, [cs.into_iter().sum::<Expr>()]);
    let mut sum: Vec<Expr> = groups
        .into_iter()
        .map(|((z, k), cs)| simplify(cs) * delta(x.clone() - z, k))
        .collect();
    sum.push(simplify(rest));
    ev.eval(&sum.into_iter().sum())
}

#[cfg(test)]
mod tests;
