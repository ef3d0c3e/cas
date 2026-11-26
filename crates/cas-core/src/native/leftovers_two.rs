//! giac's period, linearisations, phase shift, outer normal and envelope of
//! the gap ledger (CAS-640, CAS-641), under giac's names, ports of its
//! source (lin.cc, subst.cc, signalprocessing.cc, sym2poly.cc, plot.cc):
//!
//! - `Period[e, x]` (`period`, CAS-640, signalprocessing.cc `_period`): a
//!   period of `e` in the real `x` (`x` when omitted), `Infinity` when
//!   none is found, 0 when `e` is free of `x`. `e` brought `Together`, its
//!   trigonometric functions written with `exp(i u)` (`tan u` as `(e^{2iu}
//!   − 1)/(i (e^{2iu} + 1))`, as giac's `trig2exp`) and linearised
//!   (`ExpLinearize`); each `exp(a x + b)` with `a = i c`, `c` real, has
//!   the period `|2π/c|`, a `floor(a x + b)` `1/|a|` (the period then
//!   checked by `Simplify` alone, `x` outside allowed), any other function of one argument the periods of its
//!   argument; the period is their lcm when their ratios are rational, and
//!   none when one is not, when `x` occurs outside the functions, or under
//!   `Abs` or a non-integer power (giac's `period(abs(sin(x)))` is
//!   `+infinity`);
//! - `AbsLinearize[f, x]` (`linabs`, CAS-640, `_linabs`): `f` as steps
//!   (`LinearStep`'s rewriting before it collects), each step
//!   `θ(x − a) g(x)` rewritten with `|x − a|` (`reduce_Heaviside`: of the
//!   splits of `g` into factors and terms, the one giving the smallest
//!   expression, else `(|x − a| + x − a)/(2 (x − a))`), then brought
//!   `Together` while that shrinks it and its `Abs` terms collected. Bounds
//!   for the fuzzer, giac has 10: at most [`MAX_SPLIT`] factors and terms
//!   are split;
//! - `ExpLinearize[e]` (`lin`, CAS-640, lin.cc `_lin`): products and
//!   integer powers of exponentials multiplied out into a sum `Σ cᵢ
//!   exp(aᵢ)`, equal exponents merged, `sinh` and `cosh` as exponentials,
//!   `n^u` for an integer `n` as `exp(u ln n)`, inside function arguments
//!   too. Bounds for the fuzzer, giac has none: powers above
//!   [`MAX_POWER`] and sums of more than [`MAX_TERMS`] terms stay
//!   products;
//! - `ShiftPhase[e]` (`shift_phase`, CAS-640, subst.cc `shift_phase`):
//!   `sin u` as `−cos(u + π/2)`, `cos u` as `sin(u + π/2)`, `tan u` as
//!   `−1/tan(u + π/2)`, innermost first, the shifted arguments `Together`.
//!   Evaluation brings the shifts back as it does giac's printed answers;
//! - `NonRecursiveTogether[e]` (`non_recursive_normal`, CAS-640,
//!   sym2poly.cc): `Together` without reaching inside function arguments;
//! - `Envelope[F, t]`, `Envelope[F, [x, y, t]]` (`envelope`, CAS-641,
//!   plot.cc `_enveloppe`): the envelope of the curves `F = 0` (or `F`,
//!   an equation) of parameter `t`, `x` and `y` when not named: the
//!   solution `[x(t), y(t)]` of `F = 0, ∂F/∂t = 0` (a list of them when
//!   there are several). Departure: giac answers the plot of that curve.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: e in x, steps θ(x − a), lin's coefficients c and exponents a"
)]

use super::polynomial::together_outer;
use super::step::{linear, steps, with_variable};
use crate::build::i;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `Σ cᵢ exp(aᵢ)` as the pairs `(cᵢ, aᵢ)`, giac's `lin` vector.
type Lin = Vec<(Expr, Expr)>;

/// Integer powers of sums of exponentials above this stay powers.
const MAX_POWER: i64 = 32;
/// Linearised sums of more terms stay products.
const MAX_TERMS: usize = 512;
/// Factors and terms of a step's coefficient split by `linabs`.
const MAX_SPLIT: usize = 4;
/// Rounds of `linabs`'s normalisation.
const MAX_ROUNDS: usize = 16;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn together(e: Expr) -> Expr {
    Expr::call(Symbol::TOGETHER, [e])
}

fn is_zero(e: &Expr) -> bool {
    *e == Expr::from(0)
}

/// `e` with `x` replaced by `v`, evaluated.
fn at(ev: &mut Eval<'_>, e: &Expr, x: &Expr, v: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&e.substitute(&[(x.clone(), v.clone())]))
}

/// The exponent of `exp`: of `Exp[a]` or `E^a`.
fn exponent(e: &Expr) -> Option<&Expr> {
    match (head(e), e.args()) {
        (Some(Symbol::EXP), [a]) => Some(a),
        (Some(Symbol::POWER), [b, a]) if b.as_symbol() == Some(Symbol::E) => Some(a),
        _ => None,
    }
}

/// `Σ cᵢ exp(aᵢ)`, unevaluated; giac's `unlin`.
fn unlin(v: &Lin) -> Expr {
    v.iter()
        .filter(|(c, _)| !is_zero(c))
        .map(|(c, a)| c.clone() * Expr::call(Symbol::EXP, [a.clone()]))
        .sum()
}

/// `v` with the terms of equal exponents merged, the first kept in place;
/// giac's `compress`.
fn compress(ev: &mut Eval<'_>, v: Lin) -> Result<Lin, EvalError> {
    let mut out: Lin = Vec::with_capacity(v.len());
    for (c, a) in v {
        match out.iter().position(|(_, b)| *b == a) {
            Some(k) => {
                let sum = out[k].0.clone() + c;
                out[k].0 = ev.eval(&together(sum))?;
            }
            None => out.push((c, a)),
        }
    }
    Ok(out)
}

/// The product of two linearised sums; giac's `convolution`. Past
/// [`MAX_TERMS`] the product stays one term.
fn convolution(ev: &mut Eval<'_>, v: &Lin, w: &Lin) -> Result<Lin, EvalError> {
    if v.len() * w.len() > MAX_TERMS {
        let product = ev.eval(&(unlin(v) * unlin(w)))?;
        return Ok(vec![(product, Expr::from(0))]);
    }
    let mut res = Vec::with_capacity(v.len() * w.len());
    for (c, a) in v {
        for (d, b) in w {
            let coefficient = ev.eval(&(c.clone() * d.clone()))?;
            let sum = ev.eval(&together(a.clone() + b.clone()))?;
            res.push((coefficient, sum));
        }
    }
    compress(ev, res)
}

/// `e` linearised, giac's `lin(e, v)`.
fn lin(ev: &mut Eval<'_>, e: &Expr) -> Result<Lin, EvalError> {
    let View::Apply { head: h, args } = e.view() else {
        return Ok(vec![(e.clone(), Expr::from(0))]);
    };
    if let Some(a) = exponent(e) {
        let a = lin_expr(ev, a)?;
        return Ok(vec![(Expr::from(1), a)]);
    }
    match (h.as_symbol(), args) {
        (Some(Symbol::PLUS), _) => {
            let mut v = Vec::new();
            for a in args {
                v.extend(lin(ev, a)?);
            }
            compress(ev, v)
        }
        (Some(Symbol::TIMES), [first, rest @ ..]) => {
            let mut w = lin(ev, first)?;
            for a in rest {
                let v = lin(ev, a)?;
                w = convolution(ev, &w, &v)?;
            }
            Ok(w)
        }
        (Some(Symbol::POWER), [b, k]) => {
            let w = lin(ev, b)?;
            if let [(c, a)] = &w[..] {
                if is_zero(a) && c.to_integer().is_some() {
                    let log = Expr::call(Symbol::LOG, [c.clone()]);
                    return Ok(vec![(Expr::from(1), ev.eval(&(log * k.clone()))?)]);
                }
                let c = ev.eval(&c.clone().pow(k.clone()))?;
                return Ok(vec![(c, ev.eval(&(a.clone() * k.clone()))?)]);
            }
            if let Some(n @ 0..=MAX_POWER) = k.as_i64() {
                let mut r = vec![(Expr::from(1), Expr::from(0))];
                for _ in 0..n {
                    r = convolution(ev, &r, &w)?;
                }
                return Ok(r);
            }
            let c = ev.eval(&unlin(&w).pow(k.clone()))?;
            Ok(vec![(c, Expr::from(0))])
        }
        (Some(s @ (Symbol::SINH | Symbol::COSH)), [u]) => {
            let f = lin_expr(ev, u)?;
            let half = Expr::from(1) / Expr::from(2);
            let other = if s == Symbol::COSH {
                half.clone()
            } else {
                -half.clone()
            };
            Ok(vec![
                (ev.eval(&half)?, f.clone()),
                (ev.eval(&other)?, ev.eval(&-f)?),
            ])
        }
        _ => {
            let mut new = Vec::with_capacity(args.len());
            for a in args {
                new.push(lin_expr(ev, a)?);
            }
            let call = ev.eval(&Expr::apply(h, &new))?;
            Ok(vec![(call, Expr::from(0))])
        }
    }
}

/// `e` linearised and evaluated, giac's `_lin`.
fn lin_expr(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let v = lin(ev, e)?;
    ev.eval(&unlin(&v))
}

/// `ExpLinearize[e]`, giac's `lin`.
pub(super) fn exp_linearize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    lin_expr(ev, arg).map(Some)
}

/// `e` with each `Sin`, `Cos` and `Tan` shifted by `π/2`, innermost
/// first; giac's `shift_phase`.
fn shifted(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let View::Apply { head: h, args } = e.view() else {
        return Ok(e.clone());
    };
    let mut new = Vec::with_capacity(args.len());
    for a in args {
        new.push(shifted(ev, a)?);
    }
    let s = h.as_symbol();
    if let ([u], Some(Symbol::SIN | Symbol::COS | Symbol::TAN)) = (&new[..], s) {
        let half_pi = Expr::from(Symbol::PI) / Expr::from(2);
        let v = ev.eval(&together(u.clone() + half_pi))?;
        return Ok(match s {
            Some(Symbol::SIN) => -Expr::call(Symbol::COS, [v]),
            Some(Symbol::COS) => Expr::call(Symbol::SIN, [v]),
            _ => Expr::from(-1) / Expr::call(Symbol::TAN, [v]),
        });
    }
    Ok(Expr::apply(h, &new))
}

/// `ShiftPhase[e]`, giac's `shift_phase`.
pub(super) fn shift_phase(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    shifted(ev, arg).map(Some)
}

/// `NonRecursiveTogether[e]`, giac's `non_recursive_normal`.
pub(super) fn non_recursive_together(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if head(arg) == Some(Symbol::LIST) {
        let items = arg
            .args()
            .iter()
            .map(|a| Expr::call(Symbol::NON_RECURSIVE_TOGETHER, [a.clone()]));
        return ev.eval(&Expr::call(Symbol::LIST, items)).map(Some);
    }
    together_outer(ev, e, arg)
}

/// `e` with its trigonometric functions as exponentials, innermost first;
/// giac's `trig2exp` (subst.cc `sin2exp`, `cos2exp`, `tan2exp`).
fn trig_exp(e: &Expr) -> Expr {
    let View::Apply { head: h, args } = e.view() else {
        return e.clone();
    };
    let new: Vec<Expr> = args.iter().map(trig_exp).collect();
    let [u] = &new[..] else {
        return Expr::apply(h, &new);
    };
    let iu = i() * u.clone();
    let plus = || Expr::call(Symbol::EXP, [iu.clone()]);
    let minus = || Expr::from(1) / plus();
    let tan = || {
        let q = Expr::call(Symbol::EXP, [Expr::from(2) * iu.clone()]);
        (q.clone() - Expr::from(1)) / (i() * (q + Expr::from(1)))
    };
    match h.as_symbol() {
        Some(Symbol::SIN) => (plus() - minus()) / (Expr::from(2) * i()),
        Some(Symbol::COS) => (plus() + minus()) / Expr::from(2),
        Some(Symbol::TAN) => tan(),
        Some(Symbol::COT) => Expr::from(1) / tan(),
        Some(Symbol::SEC) => Expr::from(2) / (plus() + minus()),
        Some(Symbol::CSC) => Expr::from(2) * i() / (plus() - minus()),
        _ => Expr::apply(h, &new),
    }
}

/// Whether `e` holds a non-real number or `I`.
fn imaginary(e: &Expr) -> bool {
    match e.view() {
        View::Complex(_) | View::ComplexBall(_) => true,
        View::Symbol(s) => s == Symbol::I,
        View::Apply { head: h, args } => imaginary(h) || args.iter().any(imaginary),
        _ => false,
    }
}

/// The functions of `x` in `f` outside sums, products and integer powers,
/// added to `found` once each; whether `x` occurs outside them (giac's
/// `rlvarx` and `lvar`).
fn kernels(f: &Expr, x: &Expr, found: &mut Vec<Expr>) -> bool {
    if f == x {
        return true;
    }
    if !f.contains(x) {
        return false;
    }
    match (head(f), f.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), args) => {
            // Every operand is walked: no short circuit.
            let mut outside = false;
            for a in args {
                outside |= kernels(a, x, found);
            }
            outside
        }
        (Some(Symbol::POWER), [b, k])
            if k.to_integer().is_some() && b.as_symbol() != Some(Symbol::E) =>
        {
            kernels(b, x, found)
        }
        _ => {
            if !found.contains(f) {
                found.push(f.clone());
            }
            false
        }
    }
}

/// Whether `g` is periodic in `x`, `period` updated to the lcm of the
/// periods found; giac's `is_periodic_wrt`.
fn periodic(
    ev: &mut Eval<'_>,
    g: &Expr,
    x: &Expr,
    period: &mut Option<Expr>,
) -> Result<bool, EvalError> {
    let g = ev.eval(&together(g.clone()))?;
    if !g.contains(x) {
        return Ok(true);
    }
    let f = ev.eval(&trig_exp(&g))?;
    let f = lin_expr(ev, &f)?;
    let mut found = Vec::new();
    let outside = kernels(&f, x, &mut found);
    if found.is_empty() {
        return Ok(false);
    }
    let mut check = false;
    for k in &found {
        let mut unit = None;
        let argument = exponent(k).or_else(|| match (head(k), k.args()) {
            (Some(Symbol::FLOOR), [u]) => Some(u),
            _ => None,
        });
        if let Some(u) = argument {
            match linear(ev, u, x)? {
                Some((a, _)) if !is_zero(&a) && exponent(k).is_some() => {
                    let c = ev.eval(&together(a * -i()))?;
                    if imaginary(&c) {
                        return Ok(false);
                    }
                    let two_pi = Expr::from(2) * Expr::from(Symbol::PI);
                    unit = Some(ev.eval(&Expr::call(Symbol::ABS, [two_pi / c]))?);
                }
                Some((a, _)) if !is_zero(&a) => {
                    let one = Expr::from(1);
                    unit = Some(ev.eval(&(one / Expr::call(Symbol::ABS, [a])))?);
                    check = true;
                }
                _ => {
                    if !periodic(ev, u, x, period)? {
                        return Ok(false);
                    }
                }
            }
        } else {
            match (head(k), k.args()) {
                (Some(Symbol::ABS | Symbol::POWER), _) => return Ok(false),
                (_, [u]) => {
                    if !periodic(ev, u, x, period)? {
                        return Ok(false);
                    }
                }
                _ => return Ok(false),
            }
        }
        if let Some(u) = unit {
            *period = match period.take() {
                None => Some(u),
                Some(t) => {
                    let r = ev.eval(&together(t.clone() / u))?;
                    let Some(r) = r.to_rational() else {
                        return Ok(false);
                    };
                    Some(ev.eval(&(t * Expr::from(r.denom().clone())))?)
                }
            };
        }
    }
    if check && let Some(t) = period.clone() {
        let shifted = g.substitute(&[(x.clone(), x.clone() + t)]);
        let gap = ev.eval(&Expr::call(Symbol::SIMPLIFY, [g.clone() - shifted]))?;
        // giac's check decides alone: `floor(x) - x` has the period 1.
        return Ok(is_zero(&gap));
    }
    Ok(!outside)
}

/// `Period[e, x]`, giac's `period`.
pub(super) fn period(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let mut t = None;
    if !periodic(ev, f, &x, &mut t)? {
        return Ok(Some(Symbol::INFINITY.into()));
    }
    Ok(Some(t.unwrap_or_else(|| Expr::from(0))))
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

/// The terms of the expanded `g` with one factor of head `f`, grouped by
/// its argument as `(argument, coefficient)`, and the sum of the others;
/// giac's `collect_with_func`.
fn collect(ev: &mut Eval<'_>, g: &Expr, f: Symbol) -> Result<(Lin, Expr), EvalError> {
    let g = ev.eval(&Expr::call(Symbol::EXPAND, [g.clone()]))?;
    let mut groups: Lin = Vec::new();
    let mut rest = Vec::new();
    for t in terms(&g) {
        let fs = factors(&t);
        let Some(k) = fs
            .iter()
            .position(|p| head(p) == Some(f) && p.args().len() == 1)
        else {
            rest.push(t);
            continue;
        };
        let arg = fs[k].args()[0].clone();
        let coefficient: Expr = fs
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != k)
            .map(|(_, p)| p.clone())
            .product();
        match groups.iter_mut().find(|(a, _)| *a == arg) {
            Some((_, c)) => *c = c.clone() + coefficient,
            None => groups.push((arg, coefficient)),
        }
    }
    Ok((groups, rest.into_iter().sum()))
}

/// One split of a step's coefficient `cst r1 (s1 + s2)` at `a`: the parts
/// `(h1 + p1, h2 + p2)` with and without `|x − a|`; `None` when it does not
/// split (giac's `reduce_Heaviside` loop body).
fn split(
    ev: &mut Eval<'_>,
    x: &Expr,
    a: &Expr,
    r1: &Expr,
    s1: &Expr,
    s2: &Expr,
) -> Result<Option<(Expr, Expr)>, EvalError> {
    let mirror = ev.eval(&(Expr::from(2) * a.clone() - x.clone()))?;
    let folded = Expr::call(Symbol::ABS, [x.clone() - a.clone()]) + a.clone();
    let r2 = at(ev, s1, x, &mirror)? + s2.clone();
    let r2 = ev.eval(&Expr::call(Symbol::EXPAND, [r2]))?;
    let p1 = r1.clone() * at(ev, s1, x, &folded)?;
    let p2 = -(r1.clone() * at(ev, s1, x, &mirror)?);
    let (h1, h2) = if is_zero(&r2) {
        (Expr::from(0), Expr::from(0))
    } else {
        let f2 = at(ev, &r2, x, &mirror)?;
        let d = ev.eval(&(r2.clone() - f2.clone()))?;
        if is_zero(&ev.eval(&Expr::call(Symbol::SIMPLIFY, [d.clone()]))?) {
            return Ok(None);
        }
        let mu = ev.eval(&together(r1.clone() * r2.clone() / d))?;
        let h1 = at(ev, &r2, x, &folded)? * mu.clone();
        (h1, together(-(f2 * mu)))
    };
    Ok(Some((ev.eval(&(h1 + p1))?, ev.eval(&(h2 + p2))?)))
}

/// `g`'s steps `θ(x − a)` rewritten with `|x − a|`; giac's
/// `reduce_Heaviside`.
fn reduce_steps(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let (groups, mut rest) = collect(ev, g, Symbol::HEAVISIDE_THETA)?;
    let mut abs_part = Expr::from(0);
    for (arg, coefficient) in groups {
        let lin = linear(ev, &arg, x)?;
        if let Some((_, b)) = lin.filter(|(s, _)| *s == Expr::from(1)) {
            let a = ev.eval(&-b)?;
            let f = ev.eval(&together(coefficient.clone()))?;
            let (with_x, without): (Vec<Expr>, Vec<Expr>) =
                factors(&f).into_iter().partition(|p| p.contains(x));
            let cst: Expr = without.into_iter().product();
            let mut best: Option<(u32, Expr, Expr)> = None;
            if (1..=MAX_SPLIT).contains(&with_x.len()) {
                for sel in 0..1_usize << with_x.len() {
                    let (picked, others): (Vec<_>, Vec<_>) = with_x
                        .iter()
                        .enumerate()
                        .partition(|(k, _)| sel >> k & 1 == 1);
                    if others.len() > 1 {
                        continue;
                    }
                    let r1: Expr = picked.into_iter().map(|(_, p)| p.clone()).product();
                    let r2 = others
                        .first()
                        .map_or_else(|| Expr::from(1), |(_, p)| (*p).clone());
                    let ts = match terms(&r2) {
                        ts if (2..=MAX_SPLIT).contains(&ts.len()) => ts,
                        _ => vec![r2],
                    };
                    for kt in 0..1_usize << ts.len() {
                        let (s1, s2): (Vec<_>, Vec<_>) =
                            ts.iter().enumerate().partition(|(k, _)| kt >> k & 1 == 1);
                        let s1: Expr = s1.into_iter().map(|(_, t)| t.clone()).sum();
                        let s2: Expr = s2.into_iter().map(|(_, t)| t.clone()).sum();
                        let Some((with_abs, plain)) = split(ev, x, &a, &r1, &s1, &s2)? else {
                            continue;
                        };
                        let with_abs = ev.eval(&(cst.clone() * with_abs))?;
                        let plain = ev.eval(&(cst.clone() * plain))?;
                        let all =
                            abs_part.clone() + with_abs.clone() + rest.clone() + plain.clone();
                        let size = ev.eval(&together(all))?.size();
                        if best.as_ref().is_none_or(|(s, _, _)| size < *s) {
                            best = Some((size, with_abs, plain));
                        }
                    }
                }
            }
            if let Some((_, with_abs, plain)) = best {
                abs_part += with_abs;
                rest += plain;
                continue;
            }
        }
        let abs = Expr::call(Symbol::ABS, [arg.clone()]);
        abs_part += coefficient * (abs + arg.clone()) / (Expr::from(2) * arg);
    }
    ev.eval(&(abs_part + together(rest)))
}

/// `AbsLinearize[f, x]`, giac's `linabs`.
pub(super) fn abs_linearize(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let g = steps(ev, f, &x)?;
    let first = reduce_steps(ev, &g, &x)?;
    let start = first.size();
    let (mut res, mut size) = (first.clone(), start);
    let mut kept = first;
    for _ in 0..MAX_ROUNDS {
        let smaller = ev.eval(&together(res.clone()))?;
        let s = smaller.size();
        if s == size {
            if s < start {
                kept = res;
            }
            break;
        }
        res = smaller;
        size = s;
    }
    let (groups, rest) = collect(ev, &kept, Symbol::ABS)?;
    let mut sum = rest;
    for (arg, c) in groups {
        sum += together(c) * Expr::call(Symbol::ABS, [arg]);
    }
    ev.eval(&sum).map(Some)
}

/// The variables `x`, `y` and the parameter `t` of `Envelope[F, …]`.
fn envelope_variables(v: &Expr) -> Option<[Expr; 3]> {
    let user = |e: &Expr| matches!(e.view(), View::Symbol(s) if !s.is_builtin());
    let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
    match (head(v), v.args()) {
        (Some(Symbol::LIST), [t]) if user(t) => Some([x, y, t.clone()]),
        (Some(Symbol::LIST), [x, y, t]) if [x, y, t].into_iter().all(user) => {
            Some([x.clone(), y.clone(), t.clone()])
        }
        _ if user(v) => Some([x, y, v.clone()]),
        _ => None,
    }
}

/// `Envelope[F, t]`, giac's `envelope`.
pub(super) fn envelope(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, v] = e.args() else {
        return Ok(None);
    };
    let Some([x, y, t]) = envelope_variables(v) else {
        return Ok(None);
    };
    let f = match (head(f), f.args()) {
        (Some(Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
        _ => f.clone(),
    };
    let derivative = Expr::call(Symbol::D, [f.clone(), t]);
    let zero = |g: Expr| Expr::call(Symbol::EQUAL, [g, Expr::from(0)]);
    let system = Expr::call(Symbol::LIST, [zero(f), zero(derivative)]);
    let unknowns = Expr::call(Symbol::LIST, [x, y]);
    let solutions = ev.eval(&Expr::call(Symbol::SOLVE, [system, unknowns]))?;
    if head(&solutions) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let curves: Vec<Expr> = solutions
        .args()
        .iter()
        .filter(|s| head(s) == Some(Symbol::LIST) && s.args().len() == 2)
        .cloned()
        .collect();
    Ok(Some(match &curves[..] {
        [one] => one.clone(),
        _ => Expr::call(Symbol::LIST, curves),
    }))
}

#[cfg(test)]
mod tests;
